#include "search_answer.hpp"

#include <algorithm>
#include <cstring>
#include <memory>
#include <ostream>
#include <streambuf>
#include <variant>

#include "annotation/binary_matrix/base/decode_budget.hpp"
#include "cli/query.hpp"


namespace mtg {
namespace cli {

using mtg::annot::matrix::alloc_bytes;
using mtg::annot::matrix::buffer_bytes;

namespace {

// jsoncpp 1.9.3 with libstdc++ and jemalloc, as /resolve's answer account prices it
// (cli/traverse.cpp): a value placed in an array or an object holds a node of its
// container's map (96 B), an object member's key a copy of its own (its bytes and a
// terminator), an array's or object's map a block of its own (48 B), a string a copy of its
// bytes with a 4-byte length and a terminator
constexpr uint64_t kJsonValueBytes = 96;
constexpr uint64_t kJsonContainerBytes = 48;
// libstdc++'s std::string holds 15 bytes inline; a std::string is 32 B
constexpr uint64_t kInlineString = 15;
constexpr uint64_t kStringBytes = 32;
// zlib's deflate state at the server's window and memory level (2 x 128 KiB and its stream)
constexpr uint64_t kCompressorBytes = 272 << 10;
// the response's head beside its body in the transport's buffer (response_head_bytes: 256 B
// and the header fields), with room to spare
constexpr uint64_t kResponseHeadBytes = 4 << 10;

// an object member |key| (beyond the value's own heap)
uint64_t member(const char *key) {
    return kJsonValueBytes + alloc_bytes(std::strlen(key) + 1);
}
// a string value of |length| bytes (beyond its node)
uint64_t json_string(uint64_t length) {
    return alloc_bytes(length + 5);
}
// the buffer of a std::string of |length| bytes (none while it is inline)
uint64_t string_heap(uint64_t length) {
    return length > kInlineString ? alloc_bytes(length + 1) : 0;
}
// a std::string built by appending up to |length| bytes (its last growth an old and a new
// buffer, then its copy): what it holds at most
uint64_t built_string(uint64_t length) {
    return 3 * alloc_bytes(length + 1);
}

uint64_t decimal_length(uint64_t x) {
    uint64_t digits = 1;
    while (x >= 10) {
        x /= 10;
        ++digits;
    }
    return digits;
}

// The label's object (get_label_as_json): {"sample": its first ';' part[, "properties":
// {key: value, ...}]}, each further part a member of at most its bytes and a string value of
// at most its bytes; and, while it is built, a copy of its properties
SearchRecordPrice label_price(const std::string &label) {
    SearchRecordPrice price;
    price.tree = kJsonContainerBytes + member("sample");
    const size_t first = label.find(';');
    if (first == std::string::npos) {
        price.tree += json_string(label.size());
        return price;
    }
    price.tree += json_string(first) + member("properties") + kJsonContainerBytes;
    uint64_t properties = 0;
    uint64_t parts = 0;
    for (size_t begin = first + 1; begin <= label.size(); ) {
        size_t end = label.find(';', begin);
        if (end == std::string::npos)
            end = label.size();
        const uint64_t length = end - begin;
        properties += kJsonValueBytes + alloc_bytes(length + 1) + json_string(length);
        ++parts;
        begin = end + 1;
    }
    price.tree += properties;
    // the split parts and the properties' object, copied into the label's
    price.transient = buffer_bytes(parts + 1, kStringBytes) + label.size() + parts + 1
                      + properties + kJsonContainerBytes;
    return price;
}

// the text of encode_presence_mask(|mask|) (x<len>o<len>..., the trailing zeros included)
uint64_t presence_mask_length(const graph::PresenceRuns &mask) {
    uint64_t length = 0;
    size_t pos = 0;
    mask.call_runs([&](size_t begin, size_t end) {
        if (begin > pos)
            length += 1 + decimal_length(begin - pos);
        length += 1 + decimal_length(end - begin);
        pos = end;
    });
    if (mask.size > pos)
        length += 1 + decimal_length(mask.size - pos);
    return length;
}

// the text of count_run_strings(|counts|) joined by ':' (or, |verbose|, of every position's
// count), the number of its strings and the bytes of those not held inline
struct RunsText {
    uint64_t bytes = 0;
    uint64_t strings = 0;
    uint64_t long_bytes = 0;
};
RunsText count_runs_text(const graph::CountRuns &counts, bool verbose) {
    RunsText text;
    if (verbose) {
        // fmt::join(to_dense(), ":"): every position's count, a zero one digit
        uint64_t covered = 0;
        counts.call_runs([&](size_t first, size_t end, size_t value) {
            text.bytes += (end - first) * decimal_length(value);
            covered += end - first;
        });
        text.bytes += counts.size - covered;
        text.bytes += counts.size ? counts.size - 1 : 0;
        return text;
    }
    counts.call_runs([&](size_t first, size_t end, size_t value) {
        // "<i>=<v>" or "<i>-<j>=<v>"
        uint64_t length = decimal_length(first) + 1 + decimal_length(value);
        if (end != first + 1)
            length += 1 + decimal_length(end - 1);
        text.bytes += length + (text.strings ? 1 : 0);
        ++text.strings;
        if (length > kInlineString)
            text.long_bytes += alloc_bytes(length + 1);
    });
    return text;
}

// the text of every position's coordinates (|verbose|: fmt::join(coords, ",") per position,
// the positions joined by ':') and the bytes of its strings not held inline
RunsText dense_coords_text(const graph::PositionCoords &coords) {
    RunsText text;
    coords.call_positions([&](size_t, const graph::PositionCoords::Coords &list) {
        uint64_t length = list.size() - 1;
        for (uint64_t c : list) {
            length += decimal_length(c);
        }
        text.bytes += length;
        if (length > kInlineString)
            text.long_bytes += alloc_bytes(length + 1);
    });
    text.bytes += coords.size ? coords.size - 1 : 0;
    text.strings = coords.size;
    return text;
}

// the heap of a coordinate list copied out of a PositionCoords (to_dense)
uint64_t coords_copy_bytes(const graph::PositionCoords &coords) {
    uint64_t bytes = 0;
    coords.call_positions([&](size_t, const graph::PositionCoords::Coords &list) {
        bytes += graph::PositionCoords::list_bytes(list);
    });
    return bytes;
}

// A streambuf that counts the bytes written to it and keeps none
class CountingBuffer : public std::streambuf {
  public:
    uint64_t count() const { return count_; }

  protected:
    int_type overflow(int_type c) override {
        if (!traits_type::eq_int_type(c, traits_type::eof()))
            ++count_;
        return traits_type::not_eof(c);
    }
    std::streamsize xsputn(const char *, std::streamsize n) override {
        count_ += n;
        return n;
    }

  private:
    uint64_t count_ = 0;
};

// A streambuf that appends the bytes written to it to a string
class AppendingBuffer : public std::streambuf {
  public:
    explicit AppendingBuffer(std::string *out) : out_(out) {}

  protected:
    int_type overflow(int_type c) override {
        if (!traits_type::eq_int_type(c, traits_type::eof()))
            out_->push_back(traits_type::to_char_type(c));
        return traits_type::not_eof(c);
    }
    std::streamsize xsputn(const char *s, std::streamsize n) override {
        out_->append(s, n);
        return n;
    }

  private:
    std::string *out_;
};

} // namespace


uint64_t presence_mask_text_length(const graph::PresenceRuns &mask, bool verbose) {
    return verbose ? mask.size : presence_mask_length(mask);
}

uint64_t count_runs_text_length(const graph::CountRuns &counts, bool verbose) {
    return count_runs_text(counts, verbose).bytes;
}

uint64_t coords_text_length(const graph::PositionCoords &coords, bool verbose) {
    return verbose ? dense_coords_text(coords).bytes : collapsed_coords_text(coords).bytes;
}

uint64_t search_record_rows_bytes(const SeqSearchResult &result) {
    const QuerySequence &sequence = result.get_sequence();
    // its place in the pair's results (a vector of records at its last growth: an old and a
    // new buffer), its name and sequence
    uint64_t bytes = 3 * sizeof(SeqSearchResult) + string_heap(sequence.name.size())
                     + string_heap(sequence.sequence.size());
    if (result.get_alignment())
        bytes += string_heap(result.get_alignment()->cigar.size());
    std::visit([&](const auto &rows) {
        using Row = typename std::decay_t<decltype(rows)>::value_type;
        bytes += buffer_bytes(rows.capacity(), sizeof(Row));
        for (const auto &row : rows) {
            if constexpr (std::is_same_v<Row, std::string>) {
                bytes += string_heap(row.size());
            } else {
                bytes += string_heap(std::get<0>(row).size());
            }
        }
    }, result.get_result());
    return bytes;
}

SearchRecordPrice search_record_price(const SeqSearchResult &result, bool verbose_output,
                                      bool tagged, const std::string &graph,
                                      const std::string &graph_path,
                                      const std::string &annotation_path,
                                      const std::string &index_fp) {
    const QuerySequence &sequence = result.get_sequence();
    SearchRecordPrice price;
    // the record's object in the pair's answer, its description, results and cuts
    price.tree = kJsonValueBytes + kJsonContainerBytes
            + member(SeqSearchResult::SEQ_DESCRIPTION_JSON_FIELD)
            + json_string(sequence.name.size())
            + member("results") + kJsonContainerBytes
            + member("cuts") + kJsonContainerBytes;
    if (const auto &alignment = result.get_alignment()) {
        price.tree += member(SeqSearchResult::SEQUENCE_JSON_FIELD)
                + json_string(sequence.sequence.size())
                + member(SeqSearchResult::SCORE_JSON_FIELD)
                + member(SeqSearchResult::MAX_SCORE_JSON_FIELD)
                + member(SeqSearchResult::CIGAR_JSON_FIELD)
                + json_string(alignment->cigar.size())
                + member(SeqSearchResult::ORIENTATION_JSON_FIELD);
    }
    // the cut (cut_to_json), its pair tags and index_fp (null: no heap), and on a multi-graph
    // server its node in the merged answer (the pair's own left until the merge ends)
    price.tree += kJsonValueBytes + kJsonContainerBytes
            + member("returned") + member("complete") + member("nth_score")
            + member("tied") + member("tied_returned")
            + member("graph") + member("graph_path") + member("annotation_path")
            + member("index_fp");
    if (result.get_cut().abundances)
        price.tree += member("nth_median") + member("tied_same_median");
    if (!index_fp.empty())
        price.tree += json_string(index_fp.size());
    // the tags of a row or of the cut, and its node in the merged answer
    const uint64_t tags = tagged
            ? member("graph") + json_string(graph.size())
              + member("graph_path") + json_string(graph_path.size())
              + member("annotation_path") + json_string(annotation_path.size())
              + kJsonValueBytes
            : 0;
    // the cut's tags (its members are priced above)
    if (tagged) {
        price.tree += json_string(graph.size()) + json_string(graph_path.size())
                + json_string(annotation_path.size()) + kJsonValueBytes;
    }

    auto row = [&](const std::string &label, uint64_t extra, uint64_t transient) {
        const SearchRecordPrice l = label_price(label);
        price.tree += kJsonValueBytes + l.tree + extra + tags;
        price.transient = std::max(price.transient, l.transient + transient);
    };
    using LabelVec = SeqSearchResult::LabelVec;
    using LabelCountVec = SeqSearchResult::LabelCountVec;
    using LabelSigVec = SeqSearchResult::LabelSigVec;
    using LabelCountAbundancesVec = SeqSearchResult::LabelCountAbundancesVec;
    using LabelCountCoordsVec = SeqSearchResult::LabelCountCoordsVec;
    const auto &rows = result.get_result();
    if (const auto *v = std::get_if<LabelVec>(&rows)) {
        for (const auto &label : *v) {
            row(label, 0, 0);
        }
    } else if (const auto *v = std::get_if<LabelCountVec>(&rows)) {
        for (const auto &[label, count] : *v) {
            row(label, member(SeqSearchResult::KMER_COUNT_FIELD), 0);
        }
    } else if (const auto *v = std::get_if<LabelSigVec>(&rows)) {
        for (const auto &[label, count, mask] : *v) {
            // the mask's string (or, verbose, the mask itself and its string of every bit)
            const uint64_t length = verbose_output ? mask.size : presence_mask_length(mask);
            const uint64_t transient = built_string(length)
                    + (verbose_output
                       ? alloc_bytes(graph::PresenceRuns::dense_bytes(mask.size)
                                     + sizeof(uint64_t))
                       : 0);
            row(label, member(SeqSearchResult::SIGNATURE_FIELD) + kJsonContainerBytes
                        + member("presence_mask") + json_string(length) + member("score")
                        + member(SeqSearchResult::KMER_COUNT_FIELD),
                transient);
        }
    } else if (const auto *v = std::get_if<LabelCountAbundancesVec>(&rows)) {
        for (const auto &[label, count, abundances] : *v) {
            const RunsText text = count_runs_text(abundances, verbose_output);
            // the run strings (or, verbose, every position's count) and their join
            const uint64_t transient = (verbose_output
                       ? buffer_bytes(abundances.size, sizeof(size_t))
                       : 3 * buffer_bytes(text.strings, kStringBytes) + text.long_bytes)
                    + 4 * alloc_bytes(text.bytes + 1);
            row(label, member(SeqSearchResult::KMER_COUNT_FIELD)
                        + member(SeqSearchResult::KMER_ABUNDANCE_FIELD)
                        + json_string(text.bytes),
                transient);
        }
    } else {
        for (const auto &[label, count, coords, num_kmers_in_target]
                : std::get<LabelCountCoordsVec>(rows)) {
            uint64_t length;
            uint64_t transient;
            if (verbose_output) {
                // to_dense(): a list per position with copies of the coordinates; a string per
                // position (reserved) and their join
                const RunsText text = dense_coords_text(coords);
                length = text.bytes;
                transient = buffer_bytes(coords.size, sizeof(graph::PositionCoords::Coords))
                        + coords_copy_bytes(coords)
                        + buffer_bytes(coords.size, kStringBytes) + text.long_bytes;
            } else {
                // collapse_coord_ranges: a string per range (a vector at its last growth) and
                // the ranges open (at most every coordinate of a position, twice)
                const CoordRangesText text = collapsed_coords_text(coords);
                length = text.bytes;
                uint64_t open = 0;
                coords.call_positions([&](size_t, const graph::PositionCoords::Coords &list) {
                    open = std::max<uint64_t>(open, list.size());
                });
                transient = 3 * buffer_bytes(text.ranges, kStringBytes) + text.long_ranges_bytes
                        + 2 * buffer_bytes(open, 3 * sizeof(uint64_t));
            }
            transient += 4 * alloc_bytes(length + 1);
            row(label, member(SeqSearchResult::KMER_COUNT_FIELD)
                        + (num_kmers_in_target ? member(SeqSearchResult::KMERS_IN_TARGET_FIELD) : 0)
                        + member(SeqSearchResult::KMER_COORDINATE_FIELD) + json_string(length),
                transient);
        }
    }
    return price;
}

uint64_t json_text_length(const Json::Value &root, const Json::StreamWriterBuilder &builder) {
    CountingBuffer counter;
    std::ostream out(&counter);
    const std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
    writer->write(root, &out);
    return counter.count();
}

std::string json_text_of_length(const Json::Value &root,
                                const Json::StreamWriterBuilder &builder, uint64_t length) {
    std::string text;
    text.reserve(length);
    AppendingBuffer appender(&text);
    std::ostream out(&appender);
    const std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
    writer->write(root, &out);
    return text;
}

std::string search_answer_text(const Json::Value &answer, graph::SearchAccount &account,
                               bool cached) {
    // the writer of Json::writeString(Json::StreamWriterBuilder(), answer)
    const Json::StreamWriterBuilder builder;
    const uint64_t length = json_text_length(answer, builder);
    // the text's buffer beside the tree
    const uint64_t text = alloc_bytes(length + 1);
    account.charge(text);
    // once the tree is freed: the text beside the transport's copy (its head included) or
    // gzip's compressed copy (at most about the text) and its state, and the cache's copy
    const uint64_t after_tree = text + alloc_bytes(length + kResponseHeadBytes)
                                + kCompressorBytes + (cached ? text : 0);
    if (after_tree > account.held())
        account.charge(after_tree - account.held());
    return json_text_of_length(answer, builder, length);
}

Json::Value search_answer_too_large_json(uint64_t max_mb) {
    Json::Value body(Json::objectValue);
    body["code"] = "search_answer_too_large";
    body["error"] = graph::SearchAnswerTooLarge::message(max_mb);
    body["max_memory_mb"] = Json::UInt64(max_mb);
    return body;
}

Json::Value search_capabilities_json(uint64_t max_mb) {
    Json::Value c(Json::objectValue);
    c["max_memory_mb"] = max_mb ? Json::Value(Json::UInt64(max_mb)) : Json::Value();
    c["rule"] = kSearchBudgetRule;
    Json::Value refusal(Json::objectValue);
    refusal["status"] = 400;
    refusal["code"] = "search_answer_too_large";
    Json::Value members(Json::arrayValue);
    for (const char *name : { "code", "error", "max_memory_mb" }) {
        members.append(name);
    }
    refusal["members"] = std::move(members);
    c["too_large"] = std::move(refusal);
    return c;
}

} // namespace cli
} // namespace mtg

#include "annotation/coord_to_header.hpp"

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <limits>
#include <memory>
#include <stdexcept>
#include <tuple>

#include <ips4o.hpp>
#include <sdsl/bits.hpp>
#include <tsl/hopscotch_set.h>

#include "common/logger.hpp"
#include "common/serialization.hpp"
#include "common/utils/file_utils.hpp"
#include "common/utils/string_utils.hpp"
#include "common/threads/threading.hpp"
#include "graph/representation/base/sequence_graph.hpp"

namespace mtg {
namespace annot {

using Tuple = CoordToHeader::Tuple;
using mtg::common::logger;

// The reverse index of the headers (find_header): one entry of 16 bytes per header, sorted by
// (hash, column, seq_id), so that the entries of one hash list its headers in column order
struct CoordToHeader::HeaderIndex {
    struct Entry {
        uint64_t hash;
        uint32_t column;
        uint32_t seq_id;
    };
    static_assert(sizeof(Entry) == 16, "an entry is 16 bytes");

    // not a std::vector: its value-initialisation would write every entry once more on one
    // thread before the threads fill them
    std::unique_ptr<Entry[]> entries;
    size_t size = 0;
    size_t num_distinct = 0;
    // the hash the entries were sorted by
    HeaderHash hash = nullptr;

    const Entry* begin() const { return entries.get(); }
    const Entry* end() const { return entries.get() + size; }
    // the first entry with |hash|, or end()
    const Entry* lower_bound(uint64_t hash) const {
        return std::lower_bound(begin(), end(), hash,
                                [](const Entry &e, uint64_t h) { return e.hash < h; });
    }
};

namespace {

// The number of distinct headers among entries of one hash. Nearly always one: the run is a
// header repeated in several columns (all equal to the first), or a collision of 64-bit hashes,
// whose few headers are then sorted to count them.
template <class Entry, class HeaderOf>
size_t count_distinct_in_run(const Entry *first, const Entry *last, const HeaderOf &header_of) {
    const std::string_view head = header_of(*first);
    if (std::all_of(first + 1, last, [&](const Entry &e) { return header_of(e) == head; }))
        return 1;
    std::vector<std::string_view> names;
    names.reserve(last - first);
    for (const Entry *e = first; e != last; ++e) {
        names.push_back(header_of(*e));
    }
    std::sort(names.begin(), names.end());
    return std::unique(names.begin(), names.end()) - names.begin();
}

} // namespace

CoordToHeader::CoordToHeader() {}
CoordToHeader::~CoordToHeader() {}

// The index is never shared or carried over: its entries name the headers of the object that
// built it, so a copy builds its own on first use (with the same hash)
CoordToHeader::CoordToHeader(const CoordToHeader &other)
      : headers_(other.headers_), coord_offsets_(other.coord_offsets_),
        column_labels_(other.column_labels_),
        column_order_verified_(other.column_order_verified_),
        check_record_ends_(other.check_record_ends_),
        header_hash_(other.header_hash_) {}

void CoordToHeader::set_column_labels(std::vector<std::string> &&labels) {
    if (labels.size() != num_columns()) {
        throw std::invalid_argument(fmt::format("{} column labels for a record mapping of {} "
                                                "columns", labels.size(), num_columns()));
    }
    column_labels_ = std::move(labels);
    column_order_verified_ = false;
}

std::string CoordToHeader::check_columns(const std::vector<std::string> &labels) {
    if (labels.size() != num_columns()) {
        return fmt::format("it has {} columns, the annotation {} labels", num_columns(),
                           labels.size());
    }
    if (column_labels_) {
        for (size_t i = 0; i < labels.size(); ++i) {
            if ((*column_labels_)[i] != labels[i]) {
                return fmt::format("its column {} is '{}', the annotation's column {} '{}' (the "
                                   "files were passed in another order or under other names)",
                                   i, (*column_labels_)[i], i, labels[i]);
            }
        }
        column_labels_.reset();
        column_order_verified_ = true;
    }
    return "";
}

uint64_t CoordToHeader::default_header_hash(std::string_view header) {
    // std::hash of a string is 64 bits wide on the platforms built (libc++'s CityHash64,
    // libstdc++'s 64-bit Murmur); a collision costs a string comparison, nothing else
    static_assert(sizeof(size_t) == sizeof(uint64_t), "a 64-bit hash");
    return std::hash<std::string_view>()(header);
}

void CoordToHeader::set_header_hash(HeaderHash hash) {
    assert(hash);
    std::lock_guard<std::mutex> lock(header_index_mutex_);
    header_hash_ = hash;
    header_index_.reset();
}

void CoordToHeader::check_header_index_limits(uint64_t num_columns, uint64_t max_sequences) {
    if (num_columns > kMaxHeaderIndexColumns) {
        throw std::length_error(fmt::format(
                "Cannot index the sequence headers: {} columns exceed the header index's limit "
                "of {} columns (a column is stored in 32 bits)",
                num_columns, kMaxHeaderIndexColumns));
    }
    if (max_sequences > kMaxHeaderIndexRecords) {
        throw std::length_error(fmt::format(
                "Cannot index the sequence headers: a column with {} sequences exceeds the "
                "header index's limit of {} sequences per column (a sequence is stored in 32 bits)",
                max_sequences, kMaxHeaderIndexRecords));
    }
}

size_t CoordToHeader::build_header_index(size_t num_threads) const {
    return header_index(num_threads)->num_distinct;
}

std::shared_ptr<const CoordToHeader::HeaderIndex>
CoordToHeader::header_index(size_t num_threads) const {
    std::lock_guard<std::mutex> lock(header_index_mutex_);
    if (header_index_)
        return header_index_;

    if (!num_threads)
        num_threads = get_num_threads();
    using Entry = HeaderIndex::Entry;
    // offsets[c]: the entry of column c's first header
    std::vector<uint64_t> offsets(headers_.size() + 1, 0);
    uint64_t max_sequences = 0;
    for (size_t col = 0; col < headers_.size(); ++col) {
        offsets[col + 1] = offsets[col] + headers_[col].size();
        max_sequences = std::max<uint64_t>(max_sequences, headers_[col].size());
    }
    check_header_index_limits(headers_.size(), max_sequences);

    auto index = std::make_shared<HeaderIndex>();
    index->hash = header_hash_;
    index->size = offsets.back();
    index->entries.reset(new Entry[index->size]);
    Entry *entries = index->entries.get();
    const HeaderHash hash = header_hash_;

    // hashed in blocks of the flat order, not per column: an index of one column of 541M
    // headers (UniParc) is hashed on every thread as well
    constexpr uint64_t kBlock = uint64_t(1) << 16;
    const uint64_t num_blocks = (index->size + kBlock - 1) / kBlock;
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (uint64_t b = 0; b < num_blocks; ++b) {
        const uint64_t first = b * kBlock;
        const uint64_t last = std::min<uint64_t>(first + kBlock, index->size);
        // the column holding entry |first|: the last whose first entry is at or before it
        size_t col = std::upper_bound(offsets.begin(), offsets.end(), first)
                        - offsets.begin() - 1;
        for (uint64_t i = first; i < last; ++i) {
            while (i >= offsets[col + 1]) {
                ++col;
            }
            const uint64_t seq_id = i - offsets[col];
            entries[i] = { hash(headers_[col][seq_id]), static_cast<uint32_t>(col),
                           static_cast<uint32_t>(seq_id) };
        }
    }

    ips4o::parallel::sort(entries, entries + index->size, [](const Entry &a, const Entry &b) {
        return std::tie(a.hash, a.column, a.seq_id) < std::tie(b.hash, b.column, b.seq_id);
    }, num_threads);

    // the distinct headers: every run of one hash counted where it starts
    const auto header_of = [&](const Entry &e) -> std::string_view {
        return headers_[e.column][e.seq_id];
    };
    size_t num_distinct = 0;
    #pragma omp parallel for num_threads(num_threads) schedule(static) reduction(+:num_distinct)
    for (uint64_t i = 0; i < index->size; ++i) {
        if (i && entries[i - 1].hash == entries[i].hash)
            continue;
        uint64_t j = i + 1;
        while (j < index->size && entries[j].hash == entries[i].hash) {
            ++j;
        }
        num_distinct += j == i + 1 ? 1 : count_distinct_in_run(entries + i, entries + j,
                                                                 header_of);
    }
    index->num_distinct = num_distinct;

    logger->trace("Built header index with {} sequences ({} distinct) in {} columns, {} bytes",
                  index->size, num_distinct, headers_.size(), index->size * sizeof(Entry));
    header_index_ = std::move(index);
    header_index_builds_++;
    return header_index_;
}

std::optional<std::pair<CoordToHeader::Column, size_t>>
CoordToHeader::find_header(std::string_view header) const {
    const std::shared_ptr<const HeaderIndex> index = header_index();
    const uint64_t hash = index->hash(header);
    // the entries of this hash in column order: the first whose header is |header| is its
    // first occurrence
    for (auto it = index->lower_bound(hash); it != index->end() && it->hash == hash; ++it) {
        if (headers_[it->column][it->seq_id] == header)
            return std::make_pair(Column(it->column), size_t(it->seq_id));
    }
    return std::nullopt;
}

size_t CoordToHeader::num_header_index_builds() const {
    std::lock_guard<std::mutex> lock(header_index_mutex_);
    return header_index_builds_;
}

size_t CoordToHeader::header_index_bytes() const {
    std::lock_guard<std::mutex> lock(header_index_mutex_);
    return header_index_ ? header_index_->size * sizeof(HeaderIndex::Entry) : 0;
}

uint64_t CoordToHeader::find_repeated_headers(
        const std::vector<const CoordToHeader*> &indexes,
        const std::function<void(std::string_view, uint64_t,
                                 const std::vector<HeaderOccurrence>&)> &callback,
        size_t max_listed,
        size_t num_threads) {
    using Entry = HeaderIndex::Entry;
    // each object's array, held for the merge (set_header_hash() cannot free it meanwhile;
    // load() is not thread-safe with readers)
    std::vector<std::shared_ptr<const HeaderIndex>> arrays(indexes.size());
    HeaderHash hash = nullptr;
    for (size_t i = 0; i < indexes.size(); ++i) {
        if (!indexes[i])
            continue;
        arrays[i] = indexes[i]->header_index(num_threads);
        if (hash && arrays[i]->hash != hash) {
            throw std::invalid_argument("Cannot merge header indexes sorted by different "
                                        "hashes");
        }
        hash = arrays[i]->hash;
    }

    // The entries of one hash: a span of entries per array holding it, in index order (each
    // span in column order). Spans, not a list of entries, so that a header repeated in every
    // record of a degenerate .seqs costs no copy of its occurrences.
    struct Span {
        size_t index;
        const Entry *first;
        const Entry *last;
    };
    std::vector<Span> group;
    std::vector<HeaderOccurrence> listed;
    uint64_t repeated = 0;
    const auto header_of = [&](size_t index, const Entry &e) -> std::string_view {
        return indexes[index]->headers_[e.column][e.seq_id];
    };
    // the repeats among the entries of one hash: those whose headers compare equal
    const auto report = [&]() {
        uint64_t count = 0;
        for (const Span &span : group) {
            count += span.last - span.first;
        }
        if (count < 2)
            return;
        const std::string_view head = header_of(group.front().index, *group.front().first);
        const bool one_header = std::all_of(group.begin(), group.end(), [&](const Span &span) {
            return std::all_of(span.first, span.last, [&](const Entry &e) {
                return header_of(span.index, e) == head;
            });
        });
        if (one_header) {
            listed.clear();
            for (const Span &span : group) {
                for (const Entry *e = span.first; e != span.last && listed.size() < max_listed;
                        ++e) {
                    listed.push_back({ span.index, e->column, e->seq_id });
                }
            }
            ++repeated;
            callback(head, count, listed);
            return;
        }
        // a collision of hashes: the headers sorted, each one's occurrences kept in their order
        std::vector<std::pair<std::string_view, HeaderOccurrence>> all;
        for (const Span &span : group) {
            for (const Entry *e = span.first; e != span.last; ++e) {
                all.emplace_back(header_of(span.index, *e),
                                 HeaderOccurrence { span.index, e->column, e->seq_id });
            }
        }
        std::stable_sort(all.begin(), all.end(), [](const auto &a, const auto &b) {
            return a.first < b.first;
        });
        for (size_t a = 0; a < all.size(); ) {
            size_t b = a + 1;
            while (b < all.size() && all[b].first == all[a].first) {
                ++b;
            }
            if (b - a > 1) {
                listed.clear();
                for (size_t o = a; o < b && listed.size() < max_listed; ++o) {
                    listed.push_back(all[o].second);
                }
                ++repeated;
                callback(all[a].first, b - a, listed);
            }
            a = b;
        }
    };

    // the cursor of every non-empty array, in a min-heap by (its hash, its index), so that the
    // arrays holding one hash are popped in index order
    std::vector<const Entry*> cursor(indexes.size(), nullptr);
    using Head = std::pair<uint64_t, size_t>;
    std::vector<Head> heap;
    for (size_t i = 0; i < arrays.size(); ++i) {
        if (arrays[i] && arrays[i]->size) {
            cursor[i] = arrays[i]->begin();
            heap.emplace_back(cursor[i]->hash, i);
        }
    }
    const auto greater = std::greater<Head>();
    std::make_heap(heap.begin(), heap.end(), greater);

    while (!heap.empty()) {
        std::pop_heap(heap.begin(), heap.end(), greater);
        const size_t i = heap.back().second;
        heap.pop_back();
        // below the next array's hash, array i's entries are merged with no other's, and the
        // last array left is merged with none: a single array (a single-graph server) is
        // scanned without heap operations
        const bool alone = heap.empty();
        const uint64_t next = alone ? 0 : heap.front().first;
        const Entry *end = arrays[i]->end();
        const Entry *&it = cursor[i];
        while (it != end && (alone || it->hash < next)) {
            const Entry *run = it;
            while (++it != end && it->hash == run->hash) {}
            if (it - run > 1) {
                group.assign(1, Span { i, run, it });
                report();
            }
        }
        if (it == end)
            continue;
        if (it->hash > next) {
            heap.emplace_back(it->hash, i);
            std::push_heap(heap.begin(), heap.end(), greater);
            continue;
        }
        // the hash |next| is held by array i and by others: every array's entries of it,
        // gathered in index order (the heap orders equal hashes by index)
        const uint64_t h = next;
        group.clear();
        heap.emplace_back(h, i);
        std::push_heap(heap.begin(), heap.end(), greater);
        while (!heap.empty() && heap.front().first == h) {
            std::pop_heap(heap.begin(), heap.end(), greater);
            const size_t j = heap.back().second;
            heap.pop_back();
            const Entry *&jt = cursor[j];
            const Entry *run = jt;
            while (jt != arrays[j]->end() && jt->hash == h) {
                ++jt;
            }
            group.push_back(Span { j, run, jt });
            if (jt != arrays[j]->end()) {
                heap.emplace_back(jt->hash, j);
                std::push_heap(heap.begin(), heap.end(), greater);
            }
        }
        report();
    }
    return repeated;
}

namespace {

// a record index hashed and compared by the header it names (two types: the set derives from
// both)
struct HeaderOfRecordHash {
    const std::vector<std::string> *headers;
    size_t operator()(size_t i) const {
        return std::hash<std::string_view>()((*headers)[i]);
    }
};

struct HeaderOfRecordEqual {
    const std::vector<std::string> *headers;
    bool operator()(size_t a, size_t b) const {
        return (*headers)[a] == (*headers)[b];
    }
};

} // namespace

void CoordToHeader::find_duplicate_headers(const std::vector<std::string> &headers,
                                           const std::function<bool(size_t, size_t)> &callback) {
    tsl::hopscotch_set<size_t, HeaderOfRecordHash, HeaderOfRecordEqual>
            first(std::max<size_t>(headers.size(), 1),
                  HeaderOfRecordHash { &headers }, HeaderOfRecordEqual { &headers });
    for (size_t i = 0; i < headers.size(); ++i) {
        auto [it, inserted] = first.insert(i);
        if (!inserted && !callback(*it, i))
            return;
    }
}

CoordToHeader::CoordToHeader(std::vector<std::vector<std::string>> &&headers,
                             std::vector<std::vector<uint64_t>> &&num_kmers)
      : headers_(std::move(headers)), coord_offsets_(num_kmers.size()) {
    assert(headers_.size() == num_kmers.size());
    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
    for (size_t j = 0; j < num_kmers.size(); ++j) {
        assert(num_kmers[j].size() == num_sequences(j));
        if (num_kmers[j].empty())
            continue;
        auto &offsets = num_kmers[j];
        // Ensure no zero k-mer counts (each sequence must have at least one k-mer)
        assert(std::find(offsets.begin(), offsets.end(), 0) == offsets.end());
        std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
        coord_offsets_[j] = bit_vector_sd([&](const auto &callback) {
            for (uint64_t cur_coord : offsets) {
                callback(cur_coord - 1);
            }
        }, offsets.back(), offsets.size());
    }
}

void CoordToHeader::map_to_local_coords(std::vector<RowTuples> *rows) const {
    assert(rows);
    for (auto &row : *rows) {
        for (auto &[col, coords] : row) {
            const size_t n = num_sequences(col);
            for (uint64_t &coord : coords) {
                auto [seq_id, local_coord] = map_single_coord(col, coord);
                assert(n);
                if (local_coord > std::numeric_limits<uint64_t>::max() / n) {
                    throw std::runtime_error(fmt::format("Local coordinate {} is too large to "
                            "pack with a seq_id into a single 64-bit integer "
                            "({} sequences in column {})", local_coord, n, col));
                }
                coord = local_coord * n + seq_id;
            }
        }
    }
}

// A lookup in a mapping whose record ends do not increase (loaded without check_record_ends, as
// only a test does): refused, never answered with a record outside the column or a range that does
// not hold the coordinate (the review of 9726cb420: an inverted rank's underflow named record
// 2^64 - 3, and LabelRecorder's walk over a row's coordinates never moved past a range that ended
// before its coordinate)
[[noreturn]] static void refuse_corrupt_lookup(const std::string &what) {
    throw std::runtime_error("The record mapping (.seqs) is corrupt: " + what);
}

std::pair<size_t, uint64_t>
CoordToHeader::map_single_coord(Column col, uint64_t coord) const {
    if (col >= num_columns()) {
        throw std::out_of_range(fmt::format("Column {} out of range "
                "(CoordToHeader has {} columns)", col, num_columns()));
    }
    const auto &offsets = coord_offsets_[col];
    if (coord >= offsets.size()) {
        throw std::out_of_range(fmt::format("Coordinate {} for column {} out of range "
                "(CoordToHeader has {} coordinates for that column)", coord, col, offsets.size()));
    }
    // (in a mapping that loaded the record ends strictly increase: the record is one of the
    // column's, and the one before it ends before |coord|)
    const size_t header = coord ? offsets.rank1(coord - 1) : 0;
    if (header >= num_sequences(col)) {
        refuse_corrupt_lookup(fmt::format("coordinate {} of column {} ranks as record {} of {}",
                                          coord, col, header, num_sequences(col)));
    }
    if (!header)
        return { 0, coord };
    const uint64_t previous_end = offsets.select1(header);
    if (previous_end >= coord) {
        refuse_corrupt_lookup(fmt::format("coordinate {} of column {} ranks as record {}, whose "
                                          "predecessor ends at {}", coord, col, header,
                                          previous_end));
    }
    return { header, coord - previous_end - 1 };
}

CoordToHeader::SequenceRange
CoordToHeader::sequence_range(Column col, uint64_t coord) const {
    if (col >= num_columns()) {
        throw std::out_of_range(fmt::format("Column {} out of range "
                "(CoordToHeader has {} columns)", col, num_columns()));
    }
    const auto &offsets = coord_offsets_[col];
    if (coord >= offsets.size()) {
        throw std::out_of_range(fmt::format("Coordinate {} for column {} out of range "
                "(CoordToHeader has {} coordinates for that column)", coord, col, offsets.size()));
    }
    // as map_single_coord: a set bit marks the last coordinate of each sequence (and the range
    // holds |coord|, which a caller walking sorted coordinates relies on to move past it)
    const size_t header = coord ? offsets.rank1(coord - 1) : 0;
    if (header >= num_sequences(col)) {
        refuse_corrupt_lookup(fmt::format("coordinate {} of column {} ranks as record {} of {}",
                                          coord, col, header, num_sequences(col)));
    }
    const uint64_t first = header ? offsets.select1(header) + 1 : 0;
    const uint64_t last = offsets.select1(header + 1);
    if (coord < first || coord > last) {
        refuse_corrupt_lookup(fmt::format("coordinate {} of column {} ranks as record {}, whose "
                                          "coordinates are {} to {}", coord, col, header, first,
                                          last));
    }
    return { header, first, last };
}

uint64_t CoordToHeader::last_coord(Column col, size_t seq_id) const {
    if (col >= num_columns()) {
        throw std::out_of_range(fmt::format("Column {} out of range "
                "(CoordToHeader has {} columns)", col, num_columns()));
    }
    if (seq_id >= num_sequences(col)) {
        throw std::out_of_range(fmt::format("Sequence id {} out of range for column {} "
                "({} sequences)", seq_id, col, num_sequences(col)));
    }
    return coord_offsets_[col].select1(seq_id + 1);
}

uint64_t CoordToHeader::num_kmers_in_sequence(Column col, size_t seq_id) const {
    if (col >= num_columns()) {
        throw std::out_of_range(fmt::format("Column {} out of range "
                "(CoordToHeader has {} columns)", col, num_columns()));
    }
    const auto &offsets = coord_offsets_[col];
    if (seq_id >= num_sequences(col)) {
        throw std::out_of_range(fmt::format("Sequence id {} out of range for column {} "
                "({} sequences)", seq_id, col, num_sequences(col)));
    }
    // coord_offsets_ has a set bit at the partial-sum boundary of each
    // sequence, so the k-mer count for sequence s is
    //     select1(s+1) - (s == 0 ? -1 : select1(s)).
    uint64_t end = offsets.select1(seq_id + 1);
    uint64_t start = seq_id ? offsets.select1(seq_id) + 1 : 0;
    if (end < start) {
        // (start > end >= 0: seq_id > 0, and start - 1 is its predecessor's end)
        refuse_corrupt_lookup(fmt::format("record {} of column {} ends at {}, not after its "
                                          "predecessor's end at {}", seq_id, col, end, start - 1));
    }
    return end - start + 1;
}

// The bytes of |in|, a file of |size| bytes, after its read position
static uint64_t bytes_left(std::istream &in, uint64_t size) {
    const std::streamoff at = in.tellg();
    return at < 0 || static_cast<uint64_t>(at) > size ? 0 : size - static_cast<uint64_t>(at);
}

namespace {

// An sdsl int_vector as read by check_offsets_layout: its stated size in bits, its width and its
// words (little-endian, as sdsl writes them)
struct StoredIntVector {
    uint64_t bits = 0;
    uint8_t width = 1;
    std::vector<uint64_t> words;

    // the |i|-th entry (its |width| bits: an entry below bits / width)
    uint64_t get(uint64_t i) const {
        const uint64_t bit = i * width;
        const uint64_t word = bit >> 6;
        const unsigned offset = bit & 63;
        uint64_t value = words[word] >> offset;
        if (offset + width > 64)
            value |= words[word + 1] << (64 - offset);
        return width == 64 ? value : value & ((uint64_t(1) << width) - 1);
    }
    bool bit(uint64_t i) const { return (words[i >> 6] >> (i & 63)) & 1; }
};

} // namespace

// The coordinate offsets of column |column| (a bit_vector_sd: an sdsl sd_vector, then its
// inversion flag), validated before sdsl loads them, so that no load or lookup of a corrupt file
// reads outside what it stores (train-3 verification V2 and V3, GPT review 7 finding 9). The
// sd_vector is: its length n, the low bits per position wl, the low part (an int_vector of m
// entries of wl bits, the positions' low bits), the high part (a bit vector of m + 2^logm bits:
// the i-th position's high bits plus i, a set bit each, the zeros the bucket boundaries), then a
// select support over the high part's set bits and one over its unset bits (select_support_mcl:
// their count, a superblock per 4096 of them storing the position of its first, the
// mini-or-long flags, then per superblock a block: a long one storing the position of each of
// its arguments, a short one every 64th's less the superblock's). Checked:
// - each stated size against the bytes left before it is read or skipped, so that sdsl
//   allocates nothing from a corrupt size (a select support takes 48 bytes per 4096 stated bits);
// - the sizes against one another as sd_vector_builder makes them: m <= n, wl = logn - logm (of
//   n and m), the low width wl, the high part m + 2^logm bits — so it always has a zero (V3-S1:
//   a zero-less high part segfaulted at load) —, a length within what the high part encodes (a
//   stated length of 2^20 segfaulted, V2-S3), and wl at most 57 (select_0_support_sd strides
//   64 x 2^wl positions: at 58 and more it divides by zero, a trap on x86-64);
// - the high part: m set bits, and no set bit past its length (select_0_support_sd counts them);
// - each select support: over exactly the high part's set (unset) bits; its superblocks and
//   blocks of a width in 1..64 and of the number of entries their flags make them (V3-S3: one
//   width byte changed segfaulted the load of an r237 copy); and every value its select can read
//   equal to the position it stands for, against one pass over the high part (V3-S3: a corrupt
//   value read outside the high part or answered wrong records);
// - the inversion flag 0 or 1.
// Not checked here: the low part's values (check_record_ends checks the positions they make, once
// sdsl loaded them). |in| is returned to where it was
static void check_offsets_layout(std::istream &in, uint64_t file_size, uint64_t column) {
    const std::streampos start = in.tellg();
    auto fail = [&](const std::string &why) {
        throw std::runtime_error(fmt::format("the coordinate offsets of column {} {}", column, why));
    };
    auto number = [&]() {
        uint64_t n = 0;
        if (!in.read(reinterpret_cast<char *>(&n), sizeof(n)))
            fail("end inside their sizes");
        return n;
    };
    auto byte = [&]() {
        const int c = in.get();
        if (c == std::char_traits<char>::eof())
            fail("end inside their sizes");
        return static_cast<uint8_t>(c);
    };
    // an sdsl int_vector: its size in bits, its width (unless fixed: |with_width| false), then
    // its 64-bit words, read into the result when |read| (else skipped)
    auto int_vector = [&](bool with_width, bool read) {
        StoredIntVector v;
        v.bits = number();
        if (with_width)
            v.width = byte();
        const uint64_t num_words = v.bits / 64 + (v.bits % 64 != 0);
        if (const uint64_t left = bytes_left(in, file_size); num_words > left / 8)
            fail(fmt::format("state {} bytes of data, more than the {} bytes left",
                             num_words * 8, left));
        if (read) {
            v.words.resize(num_words);
            if (!in.read(reinterpret_cast<char *>(v.words.data()), num_words * 8))
                fail("end inside their data");
        } else {
            in.seekg(static_cast<std::streamoff>(num_words * 8), std::ios_base::cur);
        }
        return v;
    };

    const uint64_t length = number();
    const uint8_t low_bits_per_position = byte();
    const StoredIntVector low = int_vector(true, false);
    const StoredIntVector high = int_vector(false, true);
    if (low_bits_per_position >= 64 || !low.width || low.width > 64 || low.bits % low.width)
        fail("state an invalid width");
    const uint64_t ones = low.bits / low.width;
    if (ones > high.bits)
        fail(fmt::format("state {} positions in a high part of {} bits", ones, high.bits));
    // a position's high part indexes the zeros of the high part, one bucket each (a high part
    // without a zero encodes no length, not even 0: V3-S1)
    const uint64_t buckets = high.bits - ones;
    if ((length >> low_bits_per_position) >= buckets) {
        fail(fmt::format("state a length of {}, more than their encoding holds ({} buckets of {})",
                         length, buckets, uint64_t(1) << low_bits_per_position));
    }
    // the layout sd_vector_builder makes for m positions in a length of n
    if (ones > length)
        fail(fmt::format("state {} positions in a length of {}", ones, length));
    uint64_t logm = sdsl::bits::hi(ones) + 1;
    const uint64_t logn = sdsl::bits::hi(length) + 1;
    if (logm == logn)
        --logm;
    if (low_bits_per_position != logn - logm || low.width != low_bits_per_position
            || high.bits != ones + (uint64_t(1) << logm)) {
        fail(fmt::format("state {} low bits of width {} per position and a high part of {} bits, "
                         "where {} positions in a length of {} take {} and {} bits",
                         low_bits_per_position, low.width, high.bits, ones, length,
                         logn - logm, ones + (uint64_t(1) << logm)));
    }
    if (low_bits_per_position > 57) {
        fail(fmt::format("state {} low bits per position, more than the 57 their select support "
                         "indexes", low_bits_per_position));
    }
    // the high part's set bits, and none past its end
    uint64_t set_bits = 0;
    for (uint64_t word : high.words) {
        set_bits += sdsl::bits::cnt(word);
    }
    if (high.bits % 64 && high.words.back() >> (high.bits % 64))
        fail("state a high part with set bits past its end");
    if (set_bits != ones) {
        fail(fmt::format("state {} positions and a high part of {} set bits", ones, set_bits));
    }

    // an sdsl select_support_mcl over the |expected| set (or, |unset|, unset) bits of the high
    // part, read whole and checked against it
    auto select_support = [&](uint64_t expected, bool unset) {
        const char *what = unset ? "unset" : "set";
        const uint64_t count = number();
        if (count != expected) {
            fail(fmt::format("state a select support over {} {} bits, for the {} of their high part",
                             count, what, expected));
        }
        if (!count)
            return;
        auto invalid = [&](const std::string &part) {
            fail(fmt::format("state a select support over their {} bits with {}", what, part));
        };
        const uint64_t num_superblocks = count / 4096 + (count % 4096 != 0);
        const StoredIntVector superblocks = int_vector(true, true);
        if (!superblocks.width || superblocks.width > 64
                || superblocks.bits != num_superblocks * superblocks.width) {
            invalid(fmt::format("superblocks of {} bits of width {}, for {} superblocks",
                                superblocks.bits, superblocks.width, num_superblocks));
        }
        // a cleared flag marks a long superblock; none stored: every one is short
        const StoredIntVector flags = int_vector(false, true);
        if (flags.bits && flags.bits != num_superblocks) {
            invalid(fmt::format("{} mini-or-long flags, for {} superblocks", flags.bits,
                                num_superblocks));
        }
        auto is_long = [&](uint64_t k) { return flags.bits && !flags.bit(k); };
        std::vector<StoredIntVector> blocks;
        blocks.reserve(num_superblocks);
        for (uint64_t k = 0; k < num_superblocks; ++k) {
            blocks.push_back(int_vector(true, true));
            const StoredIntVector &block = blocks.back();
            const uint64_t entries = is_long(k) ? 4096 : 64;
            if (!block.width || block.width > 64 || block.bits != entries * block.width) {
                invalid(fmt::format("a block {} of {} bits of width {}, for {} entries", k,
                                    block.bits, block.width, entries));
            }
        }
        // every value select reads: each entry of a long block (the position of its argument),
        // each short superblock's (the position of its first argument) and each entry of its
        // block (the position of every 64th argument, less the superblock's), for the arguments
        // there are, against the high part's arguments met in one pass over its words (a long
        // superblock's own value is never read: sdsl's init_fast leaves the last one 0)
        uint64_t next = 1;  // the next argument whose position is stored (1-based)
        uint64_t before = 0;  // the arguments before the current word
        for (uint64_t w = 0; w < high.words.size() && next <= count; ++w) {
            uint64_t word = unset ? ~high.words[w] : high.words[w];
            if (w + 1 == high.words.size() && high.bits % 64)
                word &= (uint64_t(1) << (high.bits % 64)) - 1;
            const uint64_t in_word = sdsl::bits::cnt(word);
            while (next <= count && next <= before + in_word) {
                const uint64_t position
                        = w * 64 + sdsl::bits::sel(word, static_cast<uint32_t>(next - before));
                const uint64_t k = (next - 1) >> 12;
                const uint64_t offset = (next - 1) & 4095;
                if (!offset && !is_long(k) && superblocks.get(k) != position) {
                    invalid(fmt::format("superblock {} at {}, its first argument at {}", k,
                                        superblocks.get(k), position));
                }
                if (is_long(k)) {
                    if (blocks[k].get(offset) != position) {
                        invalid(fmt::format("entry {} of the long block {} at {}, its argument at "
                                            "{}", offset, k, blocks[k].get(offset), position));
                    }
                    ++next;
                } else {
                    if (blocks[k].get(offset >> 6) != position - superblocks.get(k)) {
                        invalid(fmt::format("entry {} of the short block {} at {} past the "
                                            "superblock, its argument {} past it", offset >> 6,
                                            k, blocks[k].get(offset >> 6),
                                            position - superblocks.get(k)));
                    }
                    next += 64;
                }
            }
            before += in_word;
        }
        // (every argument is met: the high part holds |count| of them, checked above)
        assert(next > count);
    };
    select_support(ones, false);
    select_support(buckets, true);
    if (byte() > 1)
        fail("state an inversion flag that is neither 0 nor 1");
    in.seekg(start);
}

// The positions the coordinate offsets of column |column| store (the record ends, or, stored
// inverted, the other coordinates), decoded as sdsl decodes them, strictly increase, the last below
// the length. The i-th is the high part's i-th set bit's position less i (its bucket), shifted by
// the low bits per position, plus the low part's i-th value: within a bucket only the low values
// order them, and check_offsets_layout reads none. sdsl's rank and select assume a set: in the
// review of 9726cb420, a column of 12 records stored inverted with its 4 low values zeroed made the
// inverted rank underflow (record 2^64 - 3) and record ranges end before the coordinates asked
// for, on which a /traverse in annotate mode looped, allocating. One pass over the loaded high
// part's words and low values (check_offsets_layout checked that the high part has one set bit per
// low value and none past its end), nothing copied. A corrupt low value that keeps the positions
// increasing still maps a coordinate to another record or states a wrong count, inside the column
static void check_record_ends(const bit_vector_sd &offsets, uint64_t column) {
    const auto &sd = offsets.data();
    const auto words = sd.high.data();
    const uint64_t num_words = (sd.high.bit_size() + 63) / 64;
    uint64_t i = 0;
    uint64_t previous = 0;
    for (uint64_t w = 0; w < num_words; ++w) {
        for (uint64_t word = words[w]; word; word &= word - 1, ++i) {
            const uint64_t bucket = w * 64 + sdsl::bits::lo(word) - i;
            const uint64_t position = (bucket << sd.wl) | sd.low[i];
            if (i && position <= previous) {
                throw std::runtime_error(fmt::format(
                        "the coordinate offsets of column {} decode their position {} as {}, not "
                        "above the previous {} (their low part is corrupt)", column, i, position,
                        previous));
            }
            previous = position;
        }
    }
    if (i && previous >= offsets.size()) {
        throw std::runtime_error(fmt::format(
                "the coordinate offsets of column {} decode their last position as {}, past their "
                "length of {}", column, previous, offsets.size()));
    }
}

bool CoordToHeader::load(const std::string &filename_base) {
    const std::string path = utils::make_suffix(filename_base, kExtension);
    std::unique_ptr<std::ifstream> in;
    try {
        // (with --mmap the stream maps the file, which throws on an empty one: a copy cut
        // before its first byte must fail its load, naming the file, not end the process)
        in = utils::open_ifstream(path);
    } catch (const std::exception &e) {
        logger->error("Cannot load CoordToHeader from '{}': {} (caught: {})", path,
                      utils::file_read_failure_detail(path), e.what());
        return false;
    }
    if (!in || !in->good()) {
        logger->error("Cannot open CoordToHeader file '{}': {}", path,
                      utils::file_read_failure_detail(path));
        return false;
    }
    {
        // the index views the headers about to be replaced
        std::lock_guard<std::mutex> lock(header_index_mutex_);
        header_index_.reset();
    }

    column_labels_.reset();
    column_order_verified_ = false;
    headers_.clear();
    coord_offsets_.clear();
    try {
        // each count is checked against the bytes left in the file before anything is
        // allocated for it (train-3 review R4, N3: 16 bytes stating 2^22 columns reached
        // 3.6 GB before they failed): a column takes at least 16 bytes (its number of headers
        // and its offsets' length), a header at least one (its length)
        const uint64_t file_size = std::filesystem::file_size(path);
        const uint64_t num_columns = load_number(*in);
        if (num_columns > bytes_left(*in, file_size) / 16) {
            throw std::runtime_error(fmt::format("{} columns stated, more than its {} bytes hold",
                                                 num_columns, file_size));
        }
        // (reserved, not built: each column is added once read whole)
        headers_.reserve(num_columns);
        coord_offsets_.reserve(num_columns);

        for (uint64_t i = 0; i < num_columns; ++i) {
            // the headers as load_string_vector reads them, their number and each length
            // bounded first by the bytes left
            const uint64_t num_headers = load_number(*in);
            const uint64_t left = bytes_left(*in, file_size);
            if (num_headers > left) {
                throw std::runtime_error(fmt::format(
                        "column {}: {} headers stated, more than its last bytes hold", i,
                        num_headers));
            }
            // added as they are read, never constructed from the stated count: a count within
            // the bound built that many empty strings first, 24-32 bytes each (V3-S4: 41 times
            // the file). Reserved up to one string per 32 bytes left, which no file amplifies
            std::vector<std::string> headers;
            headers.reserve(std::min(num_headers, left / 32));
            for (uint64_t h = 0; h < num_headers; ++h) {
                std::string header;
                if (!load_string(*in, &header, left))
                    throw std::runtime_error(fmt::format("the headers of column {}", i));
                headers.push_back(std::move(header));
            }
            // every column read whole: a file cut inside its last column would otherwise
            // reach the test below in a failed state, where peek() answers EOF as at the end
            // of a file written before the trailer
            check_offsets_layout(*in, file_size, i);
            bit_vector_sd offsets;
            if (!offsets.load(*in) || in->fail())
                throw std::runtime_error(fmt::format("the coordinate offsets of column {}", i));
            if (check_record_ends_)
                check_record_ends(offsets, i);
            // one set bit per record, the last coordinate of each, the column's last k-mer last
            if (offsets.num_set_bits() != headers.size()) {
                throw std::runtime_error(fmt::format(
                        "column {}: {} records, {} record ends in its coordinate offsets", i,
                        headers.size(), offsets.num_set_bits()));
            }
            if (!headers.empty() && !offsets[offsets.size() - 1]) {
                throw std::runtime_error(fmt::format(
                        "column {}: its coordinate offsets state {} k-mers, past its last record's "
                        "end", i, offsets.size()));
            }
            // and a column without records states none (V3-S2: a stated length loaded, and a
            // coordinate lookup in the column read record 0 of none)
            if (headers.empty() && offsets.size()) {
                throw std::runtime_error(fmt::format(
                        "column {}: no records, and its coordinate offsets state {} k-mers", i,
                        offsets.size()));
            }
            headers_.push_back(std::move(headers));
            coord_offsets_.push_back(std::move(offsets));
        }
        // the trailer of column labels, if any (a file written before it ends here)
        if (in->peek() != std::char_traits<char>::eof()) {
            if (load_number(*in) != kColumnLabelsMarker)
                throw std::runtime_error("unexpected data after the columns");
            // the number of labels compared before anything is allocated for them
            const uint64_t num_labels = load_number(*in);
            if (num_labels != num_columns) {
                throw std::runtime_error(fmt::format("{} column labels for {} columns",
                                                     num_labels, num_columns));
            }
            // added as they are read, as the headers
            std::vector<std::string> labels;
            const uint64_t left = bytes_left(*in, file_size);
            labels.reserve(std::min(num_labels, left / 32));
            for (uint64_t l = 0; l < num_labels; ++l) {
                std::string label;
                if (!load_string(*in, &label, left))
                    throw std::runtime_error("the column labels");
                labels.push_back(std::move(label));
            }
            if (in->peek() != std::char_traits<char>::eof())
                throw std::runtime_error("unexpected data after the column labels");
            column_labels_ = std::move(labels);
        }
        return true;

    } catch (const std::exception &e) {
        logger->error("Cannot load CoordToHeader from '{}': {} (caught: {})", path,
                      utils::file_read_failure_detail(path), e.what());
        return false;
    } catch (...) {
        logger->error("Cannot load CoordToHeader from '{}': {} (caught unknown exception)",
                      path, utils::file_read_failure_detail(path));
        return false;
    }
}

void CoordToHeader::serialize(const std::string &filename_base) const {
    auto fname = utils::make_suffix(filename_base, kExtension);
    std::ofstream out = utils::open_new_ofstream(fname);
    if (!out)
        throw std::ios_base::failure("Couldn't open file " + fname + " for writing");
    serialize_number(out, headers_.size());
    for (size_t i = 0; i < headers_.size(); ++i) {
        serialize_string_vector(out, headers_[i]);
        coord_offsets_[i].serialize(out);
    }
    if (column_labels_) {
        serialize_number(out, kColumnLabelsMarker);
        serialize_string_vector(out, *column_labels_);
    }
}

} // namespace annot
} // namespace mtg

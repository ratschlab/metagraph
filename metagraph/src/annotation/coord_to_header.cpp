#include "annotation/coord_to_header.hpp"

#include <algorithm>
#include <filesystem>
#include <limits>
#include <memory>
#include <stdexcept>
#include <tuple>

#include <ips4o.hpp>
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
    size_t header = coord ? offsets.rank1(coord - 1) : 0;
    uint64_t local_coord = !header ? coord : coord - offsets.select1(header) - 1;
    return { header, local_coord };
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
    // as map_single_coord: a set bit marks the last coordinate of each sequence
    const size_t header = coord ? offsets.rank1(coord - 1) : 0;
    const uint64_t first = header ? offsets.select1(header) + 1 : 0;
    return { header, first, offsets.select1(header + 1) };
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
    return end - start + 1;
}

// The bytes of |in|, a file of |size| bytes, after its read position
static uint64_t bytes_left(std::istream &in, uint64_t size) {
    const std::streamoff at = in.tellg();
    return at < 0 || static_cast<uint64_t>(at) > size ? 0 : size - static_cast<uint64_t>(at);
}

// The coordinate offsets of column |column| (a bit_vector_sd: an sdsl sd_vector, then its
// inversion flag), read for their layout before sdsl loads them: each stated size is checked
// against the bytes left in the file, and the sizes against one another as sd_vector builds
// them, so that sdsl allocates nothing from a corrupt size (a select support takes 48 bytes per
// 4096 stated bits, touched) and no rank selects past the encoded length (it segfaulted on a
// stated length of 2^20). |in| is returned to where it was (train-3 verification V2, S2, S3)
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
    // an sdsl int_vector: its size in bits, its width (unless fixed), then its 64-bit words
    auto int_vector = [&](uint8_t *width) {
        const uint64_t bits = number();
        if (width)
            *width = byte();
        const uint64_t bytes = (bits / 64 + (bits % 64 != 0)) * 8;
        if (const uint64_t left = bytes_left(in, file_size); bytes > left)
            fail(fmt::format("state {} bytes of data, more than the {} bytes left", bytes, left));
        in.seekg(static_cast<std::streamoff>(bytes), std::ios_base::cur);
        return bits;
    };
    // an sdsl select_support_mcl over the |expected| set (or unset) bits of the high part: their
    // count, the superblocks, the mini-or-long flags, then a block per 4096 of them
    auto select_support = [&](uint64_t expected, const char *what) {
        const uint64_t count = number();
        if (count != expected) {
            fail(fmt::format("state a select support over {} {} bits, for the {} of their high part",
                             count, what, expected));
        }
        if (!count)
            return;
        uint8_t width = 0;
        int_vector(&width);
        int_vector(nullptr);
        for (uint64_t i = 0; i < count / 4096 + (count % 4096 != 0); ++i)
            int_vector(&width);
    };
    const uint64_t length = number();
    const uint8_t low_bits_per_position = byte();
    uint8_t low_width = 0;
    const uint64_t low_bits = int_vector(&low_width);
    const uint64_t high_bits = int_vector(nullptr);
    if (low_bits_per_position >= 64 || !low_width || low_width > 64 || low_bits % low_width)
        fail("state an invalid width");
    const uint64_t ones = low_bits / low_width;
    if (ones > high_bits)
        fail(fmt::format("state {} positions in a high part of {} bits", ones, high_bits));
    // a position's high part indexes the zeros of the high part, one bucket each
    const uint64_t buckets = high_bits - ones;
    if (length && (length >> low_bits_per_position) >= buckets) {
        fail(fmt::format("state a length of {}, more than their encoding holds ({} buckets of {})",
                         length, buckets, uint64_t(1) << low_bits_per_position));
    }
    select_support(ones, "set");
    select_support(buckets, "unset");
    byte();
    in.seekg(start);
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
            std::vector<std::string> headers(num_headers);
            for (std::string &header : headers) {
                if (!load_string(*in, &header, left))
                    throw std::runtime_error(fmt::format("the headers of column {}", i));
            }
            // every column read whole: a file cut inside its last column would otherwise
            // reach the test below in a failed state, where peek() answers EOF as at the end
            // of a file written before the trailer
            check_offsets_layout(*in, file_size, i);
            bit_vector_sd offsets;
            if (!offsets.load(*in) || in->fail())
                throw std::runtime_error(fmt::format("the coordinate offsets of column {}", i));
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
            std::vector<std::string> labels(num_labels);
            const uint64_t left = bytes_left(*in, file_size);
            for (std::string &label : labels) {
                if (!load_string(*in, &label, left))
                    throw std::runtime_error("the column labels");
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

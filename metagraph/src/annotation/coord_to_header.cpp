#include "annotation/coord_to_header.hpp"

#include <charconv>
#include <filesystem>

#include <tsl/hopscotch_map.h>
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

// header -> (column, seq_id) of its first occurrence in column order
struct CoordToHeader::HeaderIndex {
    tsl::hopscotch_map<std::string_view, std::pair<Column, size_t>> map;
};

CoordToHeader::CoordToHeader() {}
CoordToHeader::~CoordToHeader() {}

// The index is never shared or carried over: its keys view the headers of the object that
// built it, so a copy builds its own on first use
CoordToHeader::CoordToHeader(const CoordToHeader &other)
      : headers_(other.headers_), coord_offsets_(other.coord_offsets_),
        taxids_(other.taxids_), taxonomy_seq_ids_(other.taxonomy_seq_ids_),
        taxid_ids_(other.taxid_ids_), taxid_names_(other.taxid_names_),
        taxid_columns_(other.taxid_columns_), taxonomy_view_(other.taxonomy_view_) {}

namespace {
template <class T>
bool parse_unsigned(std::string_view text, T *value) {
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), *value);
    return !text.empty() && error == std::errc() && end == text.data() + text.size()
        && (text.size() == 1 || text.front() != '0');
}
}

std::optional<CoordToHeader::Column> CoordToHeader::find_taxid(std::string_view name) const {
    uint32_t id = 0;
    if (!parse_unsigned(name, &id))
        return std::nullopt;
    auto found = taxid_columns_.find(id);
    return found == taxid_columns_.end() ? std::nullopt : std::optional<Column>(found->second);
}

std::pair<CoordToHeader::Column, uint64_t>
CoordToHeader::taxonomy_coord(Column column, uint64_t seq_id, uint64_t local_coord) const {
    const Column target = taxids_.at(column).at(seq_id);
    const uint32_t sequence = taxonomy_seq_ids_.at(column).at(seq_id);
    return { target, local_coord + (sequence ? taxonomy_view_->last_coord(target, sequence - 1) + 1 : 0) };
}

bool CoordToHeader::load_taxids(std::istream &in) {
    taxids_.clear();
    taxonomy_seq_ids_.clear();
    taxid_ids_.clear();
    taxid_names_.clear();
    taxid_columns_.clear();
    taxonomy_view_.reset();
    std::string line;
    if (!std::getline(in, line) || line != "METAGRAPH_RECORD_TAXIDS_V1")
        return false;
    std::vector<std::vector<uint32_t>> values(headers_.size());
    std::unordered_map<uint32_t, size_t> counts;
    for (size_t col = 0; col < headers_.size(); ++col) {
        values[col].reserve(headers_[col].size());
        for (size_t seq = 0; seq < headers_[col].size(); ++seq) {
            if (!std::getline(in, line))
                return false;
            const size_t a = line.find('\t'), b = line.find('\t', a + 1),
                         c = line.find('\t', b + 1);
            if (a == std::string::npos || b == std::string::npos || c == std::string::npos
                    || line.find('\t', c + 1) != std::string::npos)
                return false;
            const std::string_view row(line);
            uint64_t column = 0, sequence = 0;
            uint32_t id = 0;
            if (!parse_unsigned(row.substr(0, a), &column) || column != col
                    || !parse_unsigned(row.substr(a + 1, b - a - 1), &sequence) || sequence != seq
                    || row.substr(b + 1, c - b - 1) != headers_[col][seq]
                    || !parse_unsigned(row.substr(c + 1), &id) || !id)
                return false;
            values[col].push_back(id);
            ++counts[id];
        }
    }
    if (std::getline(in, line) || !in.eof() || in.bad())
        return false;
    std::vector<uint32_t> ids;
    ids.reserve(counts.size());
    for (const auto &[id, count] : counts)
        ids.push_back(id);
    std::sort(ids.begin(), ids.end());
    std::unordered_map<uint32_t, Column> columns;
    std::vector<std::string> names;
    std::vector<std::vector<std::string>> headers(ids.size());
    std::vector<std::vector<uint64_t>> lengths(ids.size());
    for (size_t i = 0; i < ids.size(); ++i) {
        columns.emplace(ids[i], i);
        names.push_back(std::to_string(ids[i]));
        headers[i].reserve(counts.at(ids[i]));
        lengths[i].reserve(counts.at(ids[i]));
    }
    std::vector<std::vector<uint32_t>> sequences(headers_.size());
    for (Column col = 0; col < headers_.size(); ++col) {
        sequences[col].reserve(headers_[col].size());
        for (size_t seq = 0; seq < headers_[col].size(); ++seq) {
            const Column target = columns.at(values[col][seq]);
            if (headers[target].size() > std::numeric_limits<uint32_t>::max())
                throw std::overflow_error("Too many records in one taxonomy column");
            sequences[col].push_back(headers[target].size());
            values[col][seq] = target;
            headers[target].push_back(headers_[col][seq]);
            lengths[target].push_back(num_kmers_in_sequence(col, seq));
        }
    }
    auto view = std::make_shared<CoordToHeader>(std::move(headers), std::move(lengths));
    taxids_ = std::move(values);
    taxonomy_seq_ids_ = std::move(sequences);
    taxid_ids_ = std::move(ids);
    taxid_names_ = std::move(names);
    taxid_columns_ = std::move(columns);
    taxonomy_view_ = std::move(view);
    return true;
}

size_t CoordToHeader::build_header_index() const {
    std::lock_guard<std::mutex> lock(header_index_mutex_);
    build_header_index_locked();
    return header_index_->map.size();
}

void CoordToHeader::build_header_index_locked() const {
    if (!header_index_) {
        auto index = std::make_unique<HeaderIndex>();
        size_t total = 0;
        for (const auto &column : headers_) {
            total += column.size();
        }
        index->map.reserve(total);
        for (Column col = 0; col < headers_.size(); ++col) {
            for (size_t s = 0; s < headers_[col].size(); ++s) {
                index->map.try_emplace(std::string_view(headers_[col][s]), col, s);
            }
        }
        logger->trace("Built header index with {} sequences in {} columns", total,
                      headers_.size());
        header_index_ = std::move(index);
        header_index_builds_++;
    }
}

std::optional<std::pair<CoordToHeader::Column, size_t>>
CoordToHeader::find_header(std::string_view header) const {
    std::lock_guard<std::mutex> lock(header_index_mutex_);
    build_header_index_locked();
    auto it = header_index_->map.find(header);
    if (it == header_index_->map.end())
        return std::nullopt;
    return it->second;
}

size_t CoordToHeader::num_header_index_builds() const {
    std::lock_guard<std::mutex> lock(header_index_mutex_);
    return header_index_builds_;
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

bool CoordToHeader::load(const std::string &filename_base, bool record_taxonomy) {
    taxids_.clear();
    taxonomy_seq_ids_.clear();
    taxid_ids_.clear();
    taxid_names_.clear();
    taxid_columns_.clear();
    taxonomy_view_.reset();
    const std::string path = utils::make_suffix(filename_base, kExtension);
    std::unique_ptr<std::ifstream> in = utils::open_ifstream(path);
    if (!in->good()) {
        logger->error("Cannot open CoordToHeader file '{}': {}", path,
                      utils::file_read_failure_detail(path));
        return false;
    }
    {
        // the index views the headers about to be replaced
        std::lock_guard<std::mutex> lock(header_index_mutex_);
        header_index_.reset();
    }

    try {
        uint64_t num_columns = load_number(*in);
        headers_.resize(num_columns);
        coord_offsets_.resize(num_columns);

        for (uint64_t i = 0; i < num_columns; ++i) {
            load_string_vector(*in, &headers_[i]);
            coord_offsets_[i].load(*in);
        }
        const std::string taxids_path = path.substr(0, path.size() - std::string(kExtension).size())
                + kTaxidsExtension;
        if (record_taxonomy) {
            auto taxids = utils::open_ifstream(taxids_path);
            if (!taxids->good() || !load_taxids(*taxids)) {
                logger->error("Invalid or incomplete record taxonomy sidecar '{}'", taxids_path);
                return false;
            }
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
}

} // namespace annot
} // namespace mtg

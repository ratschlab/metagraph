#include "utils.hpp"

#include <cstdio>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <filesystem>

#include <boost/multiprecision/integer.hpp>
#include <boost/functional/hash/hash.hpp>

#include "kmer.hpp"


std::size_t SmallVectorHash::operator()(const SmallVector &vector) const {
    return boost::hash_range(vector.begin(), vector.end());
}


namespace utils {

std::string remove_suffix(const std::string &str, const std::string &suffix) {
    auto actual_suffix = str.substr(
        std::max(0, static_cast<int>(str.size())
                    - static_cast<int>(suffix.size()))
    );
    return actual_suffix == suffix
            ? str.substr(0, str.size() - suffix.size())
            : str;
}

std::string join_strings(const std::vector<std::string> &strings,
                         const std::string &delimiter) {
    if (!strings.size())
        return "";

    if (strings.size() == 1)
        return strings[0];

    std::string result = strings[0];
    for (size_t i = 1; i < strings.size(); ++i) {
        result += delimiter + strings[i];
    }
    return result;
}

std::vector<std::string> split_string(const std::string &string,
                                      const std::string &delimiter) {
    if (!string.size())
        return {};

    if (!delimiter.size())
        return { string, };

    std::vector<std::string> result;

    size_t current_pos = 0;
    size_t delimiter_pos;

    while ((delimiter_pos = string.find(delimiter, current_pos))
                                             != std::string::npos) {
        if (delimiter_pos > current_pos)
            result.push_back(string.substr(current_pos, delimiter_pos - current_pos));
        current_pos = delimiter_pos + delimiter.size();
    }
    if (current_pos < string.size()) {
        result.push_back(string.substr(current_pos));
    }
    return result;
}

bool check_if_writable(const std::string &filename) {
    std::ifstream ifstream(filename, std::ios::binary);
    bool existed = ifstream.good();
    ifstream.close();

    std::ofstream ofstream(filename, std::ios::binary
                                        | std::ofstream::ios_base::app);
    bool can_write = ofstream.good();
    ofstream.close();

    if (!can_write)
        return false;

    if (!existed)
        std::remove(filename.c_str());

    return true;
}

/**
 *  Returns the input file type, given a filename
 */
std::string get_filetype(const std::string &fname) {
    size_t dotind = fname.rfind(".");
    if (dotind == std::string::npos)
        return "";

    std::string ext = fname.substr(dotind);

    if (ext == ".gz") {
        size_t nextind = fname.substr(0, dotind - 1).rfind(".");
        if (nextind == std::string::npos)
            return "";

        ext = fname.substr(nextind, dotind - nextind);
    }

    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    if (ext == ".vcf") {
        return "VCF";
    } else if ((ext == ".fq") || (ext == ".fastq")) {
        return "FASTQ";
    } else {
        return "FASTA";
    }
}

/**
 * Given a minimum number of splits,
 * generate a list of suffices from the alphabet.
 */
std::deque<std::string> generate_strings(const std::string &alphabet,
                                         size_t length) {

    std::deque<std::string> suffices = { "" };
    while (suffices[0].length() < length) {
        for (const char c : alphabet) {
            suffices.push_back(c + suffices[0]);
        }
        suffices.pop_front();
    }
    assert(suffices.size() == std::pow(alphabet.size(), length));
    return suffices;
}

uint32_t code_length(uint64_t a) {
    return a ? boost::multiprecision::msb(a) + 1 : 1;
}


ThreadPool::ThreadPool(size_t num_workers, size_t max_num_tasks)
      : max_num_tasks_(std::min(max_num_tasks, num_workers * 5)), stop_(false) {
    initialize(num_workers);
}

void ThreadPool::join() {
    size_t num_workers = workers.size();

    if (!num_workers) {
        return;
    } else {
        std::lock_guard<std::mutex> lock(queue_mutex);
        assert(!joining_);
        joining_ = true;
    }
    empty_condition.notify_all();

    for (std::thread &worker : workers) {
        worker.join();
    }
    workers.clear();

    if (!stop_)
        initialize(num_workers);
}

ThreadPool::~ThreadPool() {
    stop_ = true;
    join();
}

void ThreadPool::initialize(size_t num_workers) {
    assert(!stop_);
    assert(workers.size() == 0);
    joining_ = false;

    if (!num_workers)
        return;

    for(size_t i = 0; i < num_workers; ++i) {
        workers.emplace_back([this]() {
            while (true) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(this->queue_mutex);
                    this->empty_condition.wait(lock, [this]() {
                        return this->joining_ || !this->tasks.empty();
                    });
                    if (this->tasks.empty())
                        return;

                    task = std::move(this->tasks.front());
                    this->tasks.pop();
                    full_condition.notify_one();
                }

                task();
            }
        });
    }
}


#if _PROTEIN_GRAPH
    const std::string KmerExtractor::alphabet = "$ABCDEFGHIJKLMNOPQRSTUVWYZX";
    const KmerExtractor::TAlphabet KmerExtractor::kCharToNucleotide[128] = {
        26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
        26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
        26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
        26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
        26,  1,  2,  3,   4,  5,  6,  7,   8,  9, 10, 11,  12, 13, 14, 15,
        16, 17, 18, 19,  20, 21, 22, 23,  26, 24, 25, 26,  26, 26, 26, 26,
        26,  1,  2,  3,   4,  5,  6,  7,   8,  9, 10, 11,  12, 13, 14, 15,
        16, 17, 18, 19,  20, 21, 22, 23,  26, 24, 25, 26,  26, 26, 26, 26
    };
#elif _DNA_CASE_SENSITIVE_GRAPH
    //for case-specific DNA and RNA (U <-> T) data
    const std::string KmerExtractor::alphabet = "$ACGTNacgt";
    const KmerExtractor::TAlphabet KmerExtractor::kCharToNucleotide[128] = {
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  4, 4, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 6, 5, 7,  5, 5, 5, 8,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  9, 9, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
    };
#elif _DNA_GRAPH
    //for DNA and RNA (U <-> T) alphabets
    const std::string KmerExtractor::alphabet = "$ACGTN";
    const KmerExtractor::TAlphabet KmerExtractor::kCharToNucleotide[128] = {
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  4, 4, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  4, 4, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
    };
#else
    static_assert(false,
        "Define an alphabet: either "
        "_DNA_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
    );
#endif

static_assert(KmerExtractor::kLogSigma <= sizeof(KmerExtractor::TAlphabet) * 8,
              "Choose type for TAlphabet properly");


KmerExtractor::KmerExtractor() {
    assert(alphabet.size() <= (1llu << kLogSigma));
}

KmerExtractor::TAlphabet KmerExtractor::encode(char s) {
    assert((s >= 0 ? kCharToNucleotide[static_cast<size_t>(s)]
                   : kCharToNucleotide[0]) < alphabet.size());

    return s >= 0 ? kCharToNucleotide[static_cast<size_t>(s)]
                  : kCharToNucleotide[0];
}

char KmerExtractor::decode(TAlphabet c) {
    assert(c < alphabet.size());
    return alphabet[c];
}

std::vector<KmerExtractor::TAlphabet>
KmerExtractor::encode(const std::string &sequence) {
    std::vector<TAlphabet> seq_encoded(sequence.size());
    std::transform(sequence.begin(), sequence.end(),
                   seq_encoded.begin(), [](char c) { return encode(c); });
    return seq_encoded;
}

std::string KmerExtractor::decode(const std::vector<TAlphabet> &sequence) {
    std::string str(sequence.size(), 0);
    std::transform(sequence.begin(), sequence.end(),
                   str.begin(), [](TAlphabet x) { return decode(x); });
    return str;
}

/**
 * Break the sequence to kmers and extend the temporary kmers storage.
 * suffix does not take into account the last character.
 */
template <class KMER>
void KmerExtractor::sequence_to_kmers(const std::vector<TAlphabet> &seq,
                                      size_t k,
                                      const std::vector<TAlphabet> &suffix,
                                      Vector<KMER> *kmers) {
    assert(k);
    assert(suffix.size() < k);
    assert(kmers);

    if (seq.size() < k)
        return;

    // based on performance comparison
    // for KMer::pack_kmer and KMer::update_kmer
    if (suffix.size() > 1) {
        for (size_t i = 0; i < seq.size() - k + 1; ++i) {
            if (std::equal(suffix.begin(), suffix.end(),
                           &seq[i + k - 1] - suffix.size())) {
                kmers->emplace_back(&seq[i], k);
            }
        }
    } else {
        // initialize and add the first kmer from sequence
        auto kmer = KMER::pack_kmer(seq.data(), k);

        if (std::equal(suffix.begin(), suffix.end(),
                       &seq[k - 1] - suffix.size())) {
            kmers->emplace_back(kmer);
        }

        // add all other kmers
        for (size_t i = 1; i < seq.size() - k + 1; ++i) {
            KMER::update_kmer(k, seq[i + k - 1], seq[i + k - 2], &kmer);

            if (std::equal(suffix.begin(), suffix.end(),
                           &seq[i + k - 1] - suffix.size())) {
                kmers->emplace_back(kmer);
            }
        }
    }
}

template void KmerExtractor::sequence_to_kmers(const std::vector<TAlphabet> &seq,
                                               size_t k,
                                               const std::vector<TAlphabet> &suffix,
                                               Vector<Kmer64> *kmers);
template void KmerExtractor::sequence_to_kmers(const std::vector<TAlphabet> &seq,
                                               size_t k,
                                               const std::vector<TAlphabet> &suffix,
                                               Vector<Kmer128> *kmers);
template void KmerExtractor::sequence_to_kmers(const std::vector<TAlphabet> &seq,
                                               size_t k,
                                               const std::vector<TAlphabet> &suffix,
                                               Vector<Kmer256> *kmers);

template <class KMER>
std::string KmerExtractor::kmer_to_sequence(const KMER &kmer) {
    return kmer.to_string(alphabet);
}

template std::string KmerExtractor::kmer_to_sequence(const Kmer64 &kmer);
template std::string KmerExtractor::kmer_to_sequence(const Kmer128 &kmer);
template std::string KmerExtractor::kmer_to_sequence(const Kmer256 &kmer);

/**
 * Break the sequence to kmers and extend the temporary kmers storage.
 */
template <typename KMER>
void KmerExtractor::sequence_to_kmers(const std::string &sequence,
                                      size_t k,
                                      const std::vector<TAlphabet> &suffix,
                                      Vector<KMER> *kmers) {
    assert(k);
    // suffix does not include the last character
    assert(suffix.size() < k);

    if (sequence.size() < k)
        return;

    // encode sequence
    size_t dummy_prefix_size = suffix.size() > 0 ? k - 1 : 1;

    std::vector<TAlphabet> seq(sequence.size() + dummy_prefix_size + 1, 0);

    std::transform(sequence.begin(), sequence.end(),
                   &seq[dummy_prefix_size],
                   [](char c) { return encode(c); });

    KmerExtractor::sequence_to_kmers(seq, k, suffix, kmers);
}

template void KmerExtractor::sequence_to_kmers(const std::string &sequence,
                                               size_t k,
                                               const std::vector<TAlphabet> &suffix,
                                               Vector<Kmer64> *kmers);
template void KmerExtractor::sequence_to_kmers(const std::string &sequence,
                                               size_t k,
                                               const std::vector<TAlphabet> &suffix,
                                               Vector<Kmer128> *kmers);
template void KmerExtractor::sequence_to_kmers(const std::string &sequence,
                                               size_t k,
                                               const std::vector<TAlphabet> &suffix,
                                               Vector<Kmer256> *kmers);


const std::string KmerExtractor2Bit::alphabet = "ACGT";
const KmerExtractor2Bit::TAlphabet KmerExtractor2Bit::kCharToNucleotide[128] = {
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 1,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  3, 3, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 1,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  3, 3, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
};

KmerExtractor2Bit::KmerExtractor2Bit() {
    assert(alphabet.size() <= (1llu << kLogSigma));
}

KmerExtractor2Bit::TAlphabet KmerExtractor2Bit::encode(char s) const {
    assert((s >= 0 ? kCharToNucleotide[static_cast<size_t>(s)]
                   : kCharToNucleotide[0]) < alphabet.size());

    return s >= 0 ? kCharToNucleotide[static_cast<size_t>(s)]
                  : kCharToNucleotide[0];
}

char KmerExtractor2Bit::decode(TAlphabet c) const {
    assert(c < alphabet.size());
    return alphabet[c];
}

std::vector<KmerExtractor2Bit::TAlphabet>
KmerExtractor2Bit::encode(const std::string &sequence) const {
    std::vector<TAlphabet> seq_encoded(sequence.size());
    std::transform(sequence.begin(), sequence.end(),
                   seq_encoded.begin(), [this](char c) { return encode(c); });
    return seq_encoded;
}

std::string KmerExtractor2Bit::decode(const std::vector<TAlphabet> &sequence) const {
    std::string str(sequence.size(), 0);
    std::transform(sequence.begin(), sequence.end(),
                   str.begin(), [this](TAlphabet x) { return decode(x); });
    return str;
}

Vector<KmerExtractor2Bit::Kmer>
KmerExtractor2Bit
::sequence_to_kmers(std::vector<TAlphabet>&& seq, size_t k) const {
    assert(k);

    if (seq.size() < k)
        return {};

    Vector<Kmer> kmers;
    kmers.reserve(seq.size() - k + 1);

    // initialize and add the first kmer from sequence
    auto kmer = Kmer::pack_kmer(seq.data(), k);
    kmers.emplace_back(kmer);

    // add all other kmers
    for (size_t i = 1; i < seq.size() - k + 1; ++i) {
        Kmer::update_kmer(k, seq[i + k - 1], seq[i + k - 2], &kmer);
        kmers.emplace_back(kmer);
    }
    return kmers;
}

Vector<KmerExtractor2Bit::Kmer>
KmerExtractor2Bit::sequence_to_kmers(const std::string &sequence, size_t k) const {
    assert(k);
    return sequence_to_kmers(encode(sequence), k);
}

std::string KmerExtractor2Bit::kmer_to_sequence(const Kmer &kmer) const {
    return kmer.to_string(alphabet);
}


// indexes - positions of inserted elements in the final vector
template <typename Index, class Vector>
void insert_default_values(const std::vector<Index> &indexes, Vector *vector) {
    assert(std::is_sorted(indexes.begin(), indexes.end()));
    assert(vector);
    assert(!indexes.size() || indexes.back() < vector->size() + indexes.size());

    vector->resize(vector->size() + indexes.size());

    uint64_t i = vector->size() - 1;
    uint64_t shift = indexes.size();

    for (auto it = indexes.rbegin(); it != indexes.rend(); ++it) {
        while (i > *it) {
            assert(i - shift >= 0 && "Invalid indexes for insertion");
            (*vector)[i] = std::move((*vector)[i - shift]);
            i--;
        }
        // insert default value
        shift--;
        (*vector)[i--] = typename Vector::value_type();
    }
}

template
void insert_default_values(const std::vector<uint64_t> &, std::vector<bool> *);

template
void insert_default_values(const std::vector<uint64_t> &, sdsl::bit_vector *);

template
void insert_default_values(const std::vector<uint64_t> &, std::vector<SmallVector> *);

template <typename T>
void erase(std::vector<T> *vector, const std::vector<bool> &erase_mask) {
    assert(vector);
    assert(vector->size() == erase_mask.size());

    size_t j = 0;
    for (size_t i = 0; i < erase_mask.size(); ++i) {
        if (!erase_mask[i])
            vector->at(j++) = vector->at(i);
    }
    vector->resize(j);
}

template
void erase(std::vector<bool> *, const std::vector<bool> &);

RowsFromColumnsTransformer
::RowsFromColumnsTransformer(uint64_t num_rows,
                             const std::vector<std::string> &files)
      : num_rows_(num_rows) {
    streams_.reserve(files.size());
    for (const auto &file : files) {
        // initialize stream
        streams_.emplace_back(new VectorFileStream(file));
    }

    init_heap();
}

template <typename BitVectorPtr>
RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const std::vector<BitVectorPtr> &columns) {
    if (!columns.size()) {
        num_rows_ = 0;
        return;
    }

    num_rows_ = columns.at(0)->size();
    for (const auto &column : columns) {
        if (column->size() != num_rows_)
            throw std::runtime_error("Error: columns are not of the same size");
    }

    streams_.reserve(columns.size());
    for (const auto &col_ptr : columns) {
        streams_.emplace_back(new VectorBitStream(*col_ptr));
    }

    init_heap();
}

template
RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const std::vector<const bit_vector_small*> &);

template
RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const std::vector<std::unique_ptr<bit_vector_sd>> &);

template
RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const std::vector<std::unique_ptr<bit_vector_small>> &);

template
RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const std::vector<std::unique_ptr<bit_vector>> &);

template
RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const std::vector<std::shared_ptr<bit_vector_small>> &);


RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const bit_vector_small &columns_concatenated,
                             uint64_t column_size) {
    assert(columns_concatenated.size() % column_size == 0);

    if (!columns_concatenated.size()) {
        num_rows_ = 0;
        return;
    }

    num_rows_ = column_size;
    uint64_t num_columns = columns_concatenated.size() / column_size;
    streams_.reserve(num_columns);
    for (uint64_t i = 0; i < num_columns; ++i) {
        streams_.emplace_back(new VectorBitStream(columns_concatenated,
                                                  column_size * i,
                                                  column_size * (i + 1)));
    }

    init_heap();
}

void RowsFromColumnsTransformer::init_heap() {
    for (uint64_t column = 0; column < streams_.size(); ++column) {
        assert(streams_[column].get());

        num_set_bits_left_ += streams_[column]->values_left();

        if (streams_[column]->values_left()) {
            index_heap_.emplace(kmer_label_pair {
                streams_[column]->next_value(),
                column
            });
        }
    }
}

void RowsFromColumnsTransformer::call_next(ValueCallback callback) {
    assert(values_left() && index_heap_.size());

    auto index = std::move(index_heap_.top());

    index_heap_.pop();
    num_set_bits_left_--;

    if (streams_[index.col_id]->values_left())
        index_heap_.emplace(kmer_label_pair {
            streams_[index.col_id]->next_value(),
            index.col_id
        });

    assert(index_heap_.size() || !values_left());
    callback(index.row_id, index.col_id);
}


void call_rows(const std::function<void(const SetBitPositions &)> &callback,
               RowsFromColumnsTransformer&& transformer) {
    uint64_t cur_row = 0;
    std::vector<uint64_t> indices;

    while (transformer.values_left()) {
        transformer.call_next([&](uint64_t row, uint64_t column) {
            while (cur_row < row) {
                callback(indices);
                indices.clear();
                cur_row++;
            }
            indices.push_back(column);
        });
    }

    while (cur_row++ < transformer.rows()) {
        callback(indices);
        indices.clear();
    }
}

template <class BitVectorType>
std::vector<std::unique_ptr<bit_vector>>
transpose(const std::vector<std::unique_ptr<bit_vector>> &matrix) {
    std::vector<std::unique_ptr<bit_vector>> transposed;
    if (!matrix.size())
        return transposed;

    uint64_t num_rows = matrix.size();
    uint64_t num_columns = matrix[0]->size();
    transposed.reserve(num_columns);

    utils::call_rows(
        [&](const std::vector<uint64_t> &column_indices) {
            sdsl::bit_vector bv(num_rows, false);
            for (const auto &row_id : column_indices) {
                bv[row_id] = true;
            }
            transposed.emplace_back(new BitVectorType(std::move(bv)));
        },
        matrix
    );
    return transposed;
}

template
std::vector<std::unique_ptr<bit_vector>>
transpose<bit_vector_stat>(const std::vector<std::unique_ptr<bit_vector>> &matrix);

template
std::vector<std::unique_ptr<bit_vector>>
transpose<bit_vector_small>(const std::vector<std::unique_ptr<bit_vector>> &matrix);

template
std::vector<std::unique_ptr<bit_vector>>
transpose<bit_vector_dyn>(const std::vector<std::unique_ptr<bit_vector>> &matrix);


TempFile::TempFile(const std::string &tmp_dir)
      : tmp_file_name_((tmp_dir.size()
                          ? tmp_dir
                          : std::filesystem::temp_directory_path().string())
                                                + std::string("/tmp.XXXXXX")) {
    // create a file
    int fd = mkstemp(const_cast<char*>(tmp_file_name_.c_str()));
    if (fd == -1)
        throw std::runtime_error("Error: temp file "
                                    + tmp_file_name_ + " creation failed");
    // close the file descriptor
    close(fd);

    tmp_ostream_.reset(new std::ofstream(tmp_file_name_,
                                         std::ios::binary | std::ios::app));
    if (!tmp_ostream_->good()) {
        unlink(tmp_file_name_.c_str());
        throw std::runtime_error("Error: temp file "
                                    + tmp_file_name_ + " open failed");
    }
    state_ = APPEND;
}

TempFile::~TempFile() {
    unlink(tmp_file_name_.c_str());
}

std::ofstream& TempFile::ofstream() {
    assert(state_ == APPEND && "Can't write after reading");
    return *tmp_ostream_;
}

std::ifstream& TempFile::ifstream() {
    if (!tmp_istream_.get()) {
        state_ = READ;
        tmp_ostream_.reset();
        tmp_istream_.reset(new std::ifstream(tmp_file_name_, std::ios::binary));
    }
    return *tmp_istream_;
}


RangePartition::RangePartition(const std::vector<uint64_t> &arrangement,
                               const std::vector<size_t> &group_sizes) {
    size_t offset = 0;
    for (size_t group_size : group_sizes) {
        partition_.emplace_back(arrangement.begin() + offset,
                                arrangement.begin() + offset + group_size);
        offset += group_size;
    }

    assert(initialize_groups_and_ranks());
    initialize_groups_and_ranks();
}

RangePartition::RangePartition(std::vector<std::vector<uint64_t>>&& partition) {
    partition_.reserve(partition.size());
    for (auto &group : partition) {
        partition_.push_back({});
        partition_.back().reserve(group.size());
        for (auto value : group) {
            assert(uint64_t(value) <= std::numeric_limits<T>::max());
            partition_.back().push_back(value);
        }
        group.clear();
    }
    partition.clear();

    assert(initialize_groups_and_ranks());
    initialize_groups_and_ranks();
}

bool RangePartition::initialize_groups_and_ranks() {
    uint64_t set_size = 0;
    for (const auto &group : partition_) {
        set_size += group.size();
        if (!group.size())
            return false;
    }

    groups_.assign(set_size, -1);
    ranks_.assign(set_size, -1);

    for (size_t g = 0; g < partition_.size(); ++g) {
        const auto &group = partition_[g];

        for (size_t i = 0; i < group.size(); ++i) {
            auto value = group[i];
            if (value >= set_size
                    || groups_.at(value) != static_cast<T>(-1)
                    || ranks_.at(value) != static_cast<T>(-1))
                return false;

            groups_[value] = g;
            ranks_[value] = i;
        }
    }
    return true;
}

RangePartition::G RangePartition::group(T value) const {
    assert(value < groups_.size());
    return groups_[value];
}

RangePartition::R RangePartition::rank(T value) const {
    assert(value < ranks_.size());
    return ranks_[value];
}

RangePartition::T RangePartition::get(G group, R rank) const {
    return partition_[group][rank];
}

uint64_t RangePartition::num_groups() const {
    return partition_.size();
}

uint64_t RangePartition::size() const {
    assert(groups_.size() == ranks_.size());
    return ranks_.size();
}

bool RangePartition::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        partition_.assign(load_number(in), {});
        for (auto &group : partition_) {
            group = load_number_vector<T>(in);
        }
        return initialize_groups_and_ranks();
    } catch (...) {
        return false;
    }
}

void RangePartition::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Bad stream");

    serialize_number(out, partition_.size());
    for (const auto &group : partition_) {
        serialize_number_vector(out, group);
    }
}


// indexes are distinct and sorted
std::vector<bool> subvector(const bit_vector &col,
                            const std::vector<uint64_t> &indexes) {
    assert(indexes.size() <= col.size());

    std::vector<bool> shrinked(indexes.size(), 0);

    uint64_t max_rank = col.num_set_bits();
    if (!max_rank)
        return shrinked;

    // the case of uncompressed vector
    if (dynamic_cast<const bit_vector_stat *>(&col)) {
        for (size_t j = 0; j < indexes.size(); ++j) {
            if (col[indexes[j]])
                shrinked[j] = true;
        }
        return shrinked;
    }

    uint64_t cur_rank = 1;
    uint64_t next_pos = col.select1(1);

    for (size_t j = 0; j < indexes.size(); ++j) {
        if (indexes[j] < next_pos)
            continue;

        if (indexes[j] == next_pos) {
            shrinked[j] = true;
            continue;
        }

        // indexes[j] > next_pos
        if (col[indexes[j]]) {
            shrinked[j] = true;
            continue;
        }

        // we found a zero, update next_pos
        cur_rank = col.rank1(indexes[j]) + 1;
        if (cur_rank > max_rank)
            return shrinked;

        next_pos = col.select1(cur_rank);
    }

    return shrinked;
}

std::vector<uint64_t> sample_indexes(uint64_t universe_size,
                                     uint64_t sample_size,
                                     std::mt19937 &gen) {
    if (!universe_size)
        return {};

    sample_size = std::min(universe_size, sample_size);

    std::vector<uint64_t> indexes;
    indexes.reserve(3 * sample_size);

    if (sample_size * 10 < universe_size) {
        std::uniform_int_distribution<uint64_t> dis(0, universe_size - 1);

        while (indexes.size() < sample_size) {
            indexes.clear();
            for (size_t i = 0; i < 1.5 * sample_size; ++i) {
                indexes.push_back(dis(gen));
            }
            std::sort(indexes.begin(), indexes.end());
            indexes.erase(std::unique(indexes.begin(), indexes.end()),
                          indexes.end());
        }
    } else {
        std::bernoulli_distribution dis(2.0 * sample_size / universe_size);

        while (indexes.size() < sample_size) {
            indexes.clear();
            for (size_t i = 0; i < universe_size; ++i) {
                if (dis(gen)) {
                    indexes.push_back(i);
                }
            }
        }
    }

    std::shuffle(indexes.begin(), indexes.end(), gen);

    return std::vector<uint64_t>(indexes.begin(), indexes.begin() + sample_size);
}

} // namespace utils


template <typename T>
std::set<T> convert_to_set(const std::vector<T> &vector) {
    return std::set<T>(vector.begin(), vector.end());
}

template
std::set<std::string>
convert_to_set(const std::vector<std::string> &vector);

template
std::set<uint64_t>
convert_to_set(const std::vector<uint64_t> &vector);


std::set<std::string> convert_to_set(const std::vector<std::string> &vector) {
    return convert_to_set<std::string>(vector);
}

std::set<std::pair<std::string, size_t>>
to_set(const std::vector<std::pair<std::string, size_t>> &vector) {
    return std::set<std::pair<std::string, size_t>>(vector.begin(), vector.end());
}

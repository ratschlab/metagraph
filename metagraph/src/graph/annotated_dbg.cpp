#include "annotated_dbg.hpp"

#include <array>

#ifdef __AVX2__
#include <immintrin.h>
#endif

#include <cstdlib>

#include <Eigen/StdVector>
#include <sdsl/uint128_t.hpp>
#include <tsl/ordered_map.h>

#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "common/algorithms.hpp"
#include "common/vectors/int_vector_algorithm.hpp"

typedef std::pair<std::string, size_t> StringCountPair;

template <typename Key, typename T>
using VectorOrderedMap = tsl::ordered_map<Key, T,
                                          std::hash<Key>, std::equal_to<Key>,
                                          std::allocator<std::pair<Key, T>>,
                                          std::vector<std::pair<Key, T>>,
                                          uint64_t>;

template <typename T>
using AlignedVector = std::vector<T, Eigen::aligned_allocator<T>>;


AnnotatedSequenceGraph
::AnnotatedSequenceGraph(std::shared_ptr<SequenceGraph> graph,
                         std::unique_ptr<Annotator>&& annotation,
                         bool force_fast)
      : graph_(graph), annotator_(std::move(annotation)),
        force_fast_(force_fast) {
    assert(graph_.get());
    assert(annotator_.get());
}

AnnotatedDBG::AnnotatedDBG(std::shared_ptr<DeBruijnGraph> dbg,
                           std::unique_ptr<Annotator>&& annotation,
                           bool force_fast)
      : AnnotatedSequenceGraph(dbg, std::move(annotation), force_fast), dbg_(*dbg) {}

void AnnotatedSequenceGraph
::annotate_sequence(const std::string &sequence,
                    const std::vector<std::string> &labels) {
    assert(check_compatibility());

    std::vector<row_index> indices;
    indices.reserve(sequence.size());

    graph_->map_to_nodes(sequence, [&](node_index i) {
        if (i > 0)
            indices.push_back(graph_to_anno_index(i));
    });

    if (!indices.size())
        return;

    std::lock_guard<std::mutex> lock(mutex_);

    if (force_fast_) {
        auto row_major = dynamic_cast<annotate::RowCompressed<std::string>*>(annotator_.get());
        if (row_major) {
            row_major->add_labels_fast(indices, labels);
            return;
        }
    }

    annotator_->add_labels(indices, labels);
}

std::vector<std::string> AnnotatedDBG::get_labels(const std::string &sequence,
                                                  double presence_ratio) const {
    assert(presence_ratio >= 0.);
    assert(presence_ratio <= 1.);
    assert(check_compatibility());

    if (sequence.size() < dbg_.get_k())
        return {};

    VectorOrderedMap<row_index, size_t> index_counts;
    index_counts.reserve(sequence.size() - dbg_.get_k() + 1);

    size_t num_present_kmers = 0;
    size_t num_missing_kmers = 0;

    graph_->map_to_nodes(sequence, [&](node_index i) {
        if (i > 0) {
            index_counts[graph_to_anno_index(i)]++;
            num_present_kmers++;
        } else {
            num_missing_kmers++;
        }
    });

    size_t min_count = std::max(1.0, std::ceil(presence_ratio
                                                 * (num_present_kmers
                                                     + num_missing_kmers)));

    if (num_present_kmers < min_count)
        return {};

    return get_labels(index_counts.values_container(), min_count);
}

std::vector<std::string>
AnnotatedDBG::get_labels(const std::vector<std::pair<row_index, size_t>> &index_counts,
                         size_t min_count) const {
    assert(check_compatibility());

    auto code_counts = annotator_->count_labels(index_counts, min_count, min_count);

    std::vector<std::string> labels;
    labels.reserve(code_counts.size());

    const auto &label_encoder = annotator_->get_label_encoder();

    for (const auto &pair : code_counts) {
        assert(pair.second >= min_count);
        labels.push_back(label_encoder.decode(pair.first));
    }

    return labels;
}

std::vector<std::string>
AnnotatedSequenceGraph::get_labels(node_index index) const {
    assert(check_compatibility());
    assert(index != SequenceGraph::npos);

    return annotator_->get(graph_to_anno_index(index));
}

std::vector<StringCountPair>
AnnotatedDBG::get_top_labels(const std::string &sequence,
                             size_t num_top_labels,
                             double presence_ratio) const {
    assert(presence_ratio >= 0.);
    assert(presence_ratio <= 1.);
    assert(check_compatibility());

    if (sequence.size() < dbg_.get_k())
        return {};

    if (presence_ratio == 1.) {
        auto labels = get_labels(sequence, presence_ratio);
        if (labels.size() > num_top_labels)
            labels.erase(labels.begin() + num_top_labels, labels.end());

        std::vector<StringCountPair> label_counts;
        label_counts.reserve(labels.size());

        for (auto&& label : labels) {
            label_counts.emplace_back(std::move(label),
                                      sequence.size() - dbg_.get_k() + 1);
        }

        return label_counts;
    }

    VectorOrderedMap<row_index, size_t> index_counts;
    index_counts.reserve(sequence.size() - dbg_.get_k() + 1);

    size_t num_present_kmers = 0;
    size_t num_missing_kmers = 0;

    graph_->map_to_nodes(sequence, [&](node_index i) {
        if (i > 0) {
            index_counts[graph_to_anno_index(i)]++;
            num_present_kmers++;
        } else {
            num_missing_kmers++;
        }
    });

    uint64_t min_count = std::max(1.0, std::ceil(presence_ratio
                                                    * (num_present_kmers
                                                        + num_missing_kmers)));
    if (num_present_kmers < min_count)
        return {};

    return get_top_labels(index_counts.values_container(), num_top_labels, min_count);
}

std::vector<StringCountPair>
AnnotatedDBG::get_top_labels(const std::vector<std::pair<row_index, size_t>> &index_counts,
                             size_t num_top_labels,
                             size_t min_count) const {
    assert(check_compatibility());

    auto code_counts = annotator_->count_labels(index_counts, min_count);

    assert(std::all_of(
        code_counts.begin(), code_counts.end(),
        [&](const auto &code_count) { return code_count.second >= min_count; }
    ));

    // leave only first |num_top_labels| top labels
    if (num_top_labels <= code_counts.size()) {
        std::sort(code_counts.begin(), code_counts.end(),
                  [](const auto &first, const auto &second) {
                      return first.second > second.second
                          || (first.second == second.second && first.first < second.first);
                  });

        code_counts.erase(code_counts.begin() + num_top_labels, code_counts.end());
    }

    const auto &label_encoder = annotator_->get_label_encoder();

    // TODO: remove this step?
    std::vector<StringCountPair> label_counts(code_counts.size());
    for (size_t i = 0; i < code_counts.size(); ++i) {
        label_counts[i].first = label_encoder.decode(code_counts[i].first);
        label_counts[i].second = code_counts[i].second;
    }

    return label_counts;
}

std::vector<std::pair<std::string, sdsl::bit_vector>>
AnnotatedDBG::get_top_label_signatures(const std::string &sequence,
                                       size_t num_top_labels,
                                       double presence_ratio) const {
    assert(presence_ratio >= 0.);
    assert(presence_ratio <= 1.);
    assert(check_compatibility());

    if (sequence.size() < dbg_.get_k())
        return {};

    std::vector<std::pair<std::string, sdsl::bit_vector>> presence_vectors;

    if (presence_ratio == 1.) {
        auto labels = get_labels(sequence, presence_ratio);
        if (labels.size() > num_top_labels)
            labels.erase(labels.begin() + num_top_labels, labels.end());

        presence_vectors.reserve(labels.size());

        for (auto&& label : labels) {
            presence_vectors.emplace_back(
                std::move(label),
                sdsl::bit_vector(sequence.size() - dbg_.get_k() + 1, true)
            );
        }

        return presence_vectors;
    }

    size_t num_present_kmers = 0;

    // map each index in the query k-mer sequence to a row in the annotation matrix
    std::vector<row_index> row_indices;
    std::vector<node_index> node_indices;
    row_indices.reserve(sequence.size() - dbg_.get_k() + 1);
    node_indices.reserve(sequence.size() - dbg_.get_k() + 1);

    graph_->map_to_nodes(sequence, [&](node_index i) {
        node_indices.push_back(i);

        if (i > 0) {
            row_indices.push_back(graph_to_anno_index(i));
            ++num_present_kmers;
        }
    });

    assert(node_indices.size() == sequence.size() - dbg_.get_k() + 1);

    uint64_t min_count = std::max(1.0, std::ceil(presence_ratio * node_indices.size()));
    if (num_present_kmers < min_count)
        return {};

    // construct the k-mer presence masks for all label codes
    using VectorCount = std::pair<std::vector<uint8_t>, size_t>;
    using Storage = std::vector<std::pair<uint64_t, VectorCount>>;

    VectorOrderedMap<uint64_t, VectorCount> label_codes_to_presence;
    static_assert(std::is_same<
        typename decltype(label_codes_to_presence)::values_container_type,
        Storage
    >::value);

    auto label_codes = annotator_->get_label_codes(row_indices);
    assert(label_codes.size() == row_indices.size());
    auto it = label_codes.begin();
    for (size_t i = 0; i < node_indices.size(); ++i) {
        // skip this k-mer since it was not found in the graph
        if (!node_indices[i])
            continue;

        assert(it != label_codes.end());
        for (const auto &code : *it) {
            // for each label code associated with this k-mer, mark the
            // corresponding k-mer presence mask at index i
            auto &vector_count = label_codes_to_presence[code];

            if (vector_count.first.empty())
                vector_count.first.resize(node_indices.size(), 0);

            vector_count.second += !vector_count.first[i];
            vector_count.first[i] = 0xFF;
        }

        ++it;
    }

    assert(it == label_codes.end());

    // sort, decode, output
    auto vector = const_cast<Storage&&>(label_codes_to_presence.values_container());

    bool is_sorted = num_top_labels <= vector.size();
    if (is_sorted) {
        std::sort(vector.begin(), vector.end(),
                  [](const auto &a, const auto &b) {
                      return a.second.second > b.second.second
                          || (a.second.second == b.second.second && a.first < b.first);
                  });
    }

    presence_vectors.reserve(std::min(num_top_labels, vector.size()));

    for (const auto &[code, mask_count] : vector) {
        assert(mask_count.second == sdsl::util::cnt_one_bits(to_sdsl(mask_count.first)));
        if (mask_count.second < min_count) {
            if (is_sorted) {
                break;
            } else {
                continue;
            }
        }

        if (presence_vectors.size() == num_top_labels)
            break;

        // TODO: remove the decoding step?
        presence_vectors.emplace_back(annotator_->get_label_encoder().decode(code),
                                      to_sdsl(mask_count.first));
        assert(sdsl::util::cnt_one_bits(presence_vectors.back().second)
            == mask_count.second);
        assert(presence_vectors.back().second.size()
            == sequence.size() - dbg_.get_k() + 1);
    }

#ifndef NDEBUG
    // sanity check, make sure that the same matches are output by get_top_labels
    auto top_labels = get_top_labels(sequence, num_top_labels, presence_ratio);
    assert(top_labels.size() == presence_vectors.size());

    if (is_sorted) {
        assert(std::equal(presence_vectors.begin(), presence_vectors.end(),
                          top_labels.begin(),
                          [](const auto &a, const auto &b) {
                              return sdsl::util::cnt_one_bits(a.second) == b.second;
                          }));
        assert(std::equal(presence_vectors.begin(), presence_vectors.end(),
                          top_labels.begin(),
                          [](const auto &a, const auto &b) { return a.first == b.first; }));
    } else {
        std::unordered_map<std::string, uint64_t> check(top_labels.begin(), top_labels.end());
        for (const auto &[label, mask] : presence_vectors) {
            auto find = check.find(label);
            assert(find != check.end());
            assert(find->second == sdsl::util::cnt_one_bits(mask));
            check.erase(find);
        }
        assert(check.empty());
    }
#endif

    return presence_vectors;
}

bool AnnotatedSequenceGraph::label_exists(const std::string &label) const {
    return annotator_->label_exists(label);
}

bool AnnotatedSequenceGraph::has_label(node_index index, const std::string &label) const {
    assert(check_compatibility());
    assert(index != SequenceGraph::npos);

    return annotator_->has_label(graph_to_anno_index(index), label);
}

void AnnotatedSequenceGraph
::call_annotated_nodes(const std::string &label,
                       std::function<void(node_index)> callback) const {
    assert(check_compatibility());

    annotator_->call_objects(
        label,
        [&](row_index index) { callback(anno_to_graph_index(index)); }
    );
}

bool AnnotatedSequenceGraph::check_compatibility() const {
    return graph_->max_index() == annotator_->num_objects();
}


/**
 * Helper functions for score_kmer_presence_mask
 */
sdsl::bit_vector smooth_bit_vector(const sdsl::bit_vector &vector) {
    if (vector.size() < 3)
        return vector;

    auto presence = vector;

    // process one word at a time
    // TODO: is it worth vectorizing this?
    size_t i = 0;
    sdsl::uint128_t dword(vector.data()[0]);
    dword <<= 64;
    for (; i + 64 <= presence.size() - 2; i += 64) {
        dword = (dword >> 64) | (sdsl::uint128_t(vector.data()[(i >> 6) + 1]) << 64);
        presence.data()[i >> 6] &= uint64_t(dword >> 1) & uint64_t(dword >> 2);
    }

    assert(presence.size() - i >= 2);
    assert(presence.size() - i <= 126);

    // handle last word
    if (vector.size() - i <= 64) {
        uint64_t word = vector.get_int(i, vector.size() - i);
        presence.set_int(i,
                         word & ((word >> 1) | (uint64_t(1) << (vector.size() - i - 1)))
                              & ((word >> 2) | (uint64_t(3) << (vector.size() - i - 2))),
                         vector.size() - i);
    } else {
        dword = (dword >> 64) | (sdsl::uint128_t(vector.data()[(i >> 6) + 1]) << 64)
            | sdsl::uint128_t(3) << (vector.size() - i);
        dword = dword & (dword >> 1) & (dword >> 2);
        presence.set_int(i, uint64_t(dword));
        presence.set_int(i + 64, uint64_t(dword >> 64), vector.size() - i - 64);
    }

#ifndef NDEBUG
    for (size_t i = 0; i < vector.size() - 2; ++i) {
        if (presence[i] != (vector[i] && vector[i + 1] && vector[i + 2])) {
            std::cerr << i << " " << vector.size() << " " << vector.capacity() << "\n"
                      << vector << "\n" << presence << std::endl;
        }

        assert(presence[i] == (vector[i] && vector[i + 1] && vector[i + 2]));
    }

    assert(presence[vector.size() - 2]
            == (vector[vector.size() - 2] && vector[vector.size() - 1]));
    assert(presence[vector.size() - 1] == vector[vector.size() - 1]);
#endif

    return presence;
}

std::array<AlignedVector<size_t>, 2>
tabulate_score(const sdsl::bit_vector &presence, size_t correction = 0) {
    std::array<AlignedVector<size_t>, 2> table;
    table[0].reserve(presence.size());
    table[1].reserve(presence.size());

    if (!presence.size())
        return table;

    bool last_block = presence[0];
    size_t last_size = 1;
    size_t i = 1;

    if (presence.size() >= 64) {
        const auto word = *presence.data();
        if (!word || word == 0xFFFFFFFFFFFFFFFF) {
            last_size = 64;
            i = 64;
        }
    }

    for (; i < presence.size(); ++i) {
        if (!(i & 0x3F) && i + 64 <= presence.size()) {
            // if at a word boundary and the next word is either all zeros or
            // all ones
            const uint64_t word = presence.get_int(i);
            if (!word || word == 0xFFFFFFFFFFFFFFFF) {
                if ((!word && last_block) || (word == 0xFFFFFFFFFFFFFFFF && !last_block)) {
                    table[last_block].push_back(last_size + correction);
                    last_block = !last_block;
                    last_size = 0;
                }

                last_size += 64;
                i += 63;
                continue;
            }
        }

        if (last_block == presence[i]) {
            ++last_size;
        } else {
            table[last_block].push_back(last_size + correction);
            last_block = !last_block;
            last_size = 1;
        }
    }

    table[last_block].push_back(last_size);

    assert(std::accumulate(table[0].begin(), table[0].end(), size_t(0))
         + std::accumulate(table[1].begin(), table[1].end(), size_t(0))
         - correction * (table[0].size() + table[1].size() - 1)
        == presence.size());

    assert(std::accumulate(table[1].begin(), table[1].end(), size_t(0))
         - correction * (table[1].size() - last_block)
        == sdsl::util::cnt_one_bits(presence));

#ifndef NDEBUG
    if (correction == 1) {
        std::array<AlignedVector<size_t>, 2> check;

        size_t cnt = 1;
        for (size_t i = 0; i < presence.size(); ++i) {
            if (i < presence.size() - 1) {
                ++cnt;
                if (presence[i] != presence[i + 1]) {
                    check[presence[i]].push_back(cnt);
                    cnt = 1;
                }
            } else {
                check[presence[i]].push_back(cnt);
            }
        }

        assert(check[0] == table[0]);
        assert(check[1] == table[1]);
    }
#endif

    return table;
}

#ifdef __AVX2__
__m256d get_penalty_bigsi_avx2(__m256d counts,
                               __m256d match_score,
                               __m256d mismatch_score,
                               __m256d SNP_t) {
    __m256d neg_ones = _mm256_set1_pd(-1.0);

    __m256d min_N_snps = _mm256_div_pd(counts, SNP_t);

    // subtract -1 instead of add 1 to reuse a register
    __m256d max_N_snps = _mm256_max_pd(
        _mm256_sub_pd(_mm256_sub_pd(counts, SNP_t), neg_ones),
        min_N_snps
    );

    // separate mul and add is faster than using FMA instructions
    __m256d mean_N_snps = _mm256_add_pd(
        _mm256_mul_pd(max_N_snps, _mm256_set1_pd(0.05)),
        min_N_snps
    );

    assert(!_mm256_movemask_pd(_mm256_cmp_pd(mean_N_snps, counts, 14)));

    __m256d mean_penalty = _mm256_mul_pd(mean_N_snps, mismatch_score);

    __m256d penalties = _mm256_mul_pd(
        _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(mean_penalty, counts), match_score),
                      mean_penalty),
        neg_ones
    );

    return penalties;
}

// TODO: create a simd_utils.hpp and move this there so it can be in unit tests
// based off of:
// https://stackoverflow.com/questions/41144668/how-to-efficiently-perform-double-int64-conversions-with-sse-avx
__m256d uint64_to_double(__m256i x) {
    __m256d mask = _mm256_set1_pd(0x0010000000000000);

    // this part only works for sizes < 2^52
    if (LIKELY(!_mm256_movemask_epi8(_mm256_cmpgt_epi64(x, _mm256_set1_epi64x(0xFFFFFFFFFFFFF)))))
        return _mm256_sub_pd(_mm256_castsi256_pd(_mm256_or_si256(x, _mm256_castpd_si256(mask))), mask);

    __m256i xH = _mm256_srli_epi64(x, 32);
    xH = _mm256_or_si256(xH, _mm256_castpd_si256(_mm256_set1_pd(19342813113834066795298816.)));          //  2^84
    __m256i xL = _mm256_blend_epi16(x, _mm256_castpd_si256(mask), 0xcc);   //  2^52
    __m256d f = _mm256_sub_pd(_mm256_castsi256_pd(xH), _mm256_set1_pd(19342813118337666422669312.));     //  2^84 + 2^52
    return _mm256_add_pd(f, _mm256_castsi256_pd(xL));
}

#endif

double add_penalty_bigsi(double cur,
                         double count,
                         double match_score,
                         double mismatch_score,
                         double SNP_t) {
    double min_N_snps = count / SNP_t;
    double max_N_snps = std::max(count - SNP_t + 1, min_N_snps);
    double mean_N_snps = min_N_snps + 0.05 * max_N_snps;

    assert(count >= mean_N_snps);

    double mean_penalty = mean_N_snps * mismatch_score;
    return cur - mean_penalty + (count - mean_penalty) * match_score;
}

#ifdef __AVX2__
double haddall_pd(__m256d v) {
    // [ a, b, c, d ] -> [ a+b, 0, c+d, 0]
    __m256d pairs = _mm256_hadd_pd(v, _mm256_setzero_pd());

    // [ a+b, 0, c+d, 0] + [ c+d, 0, 0, 0] -> a+b+c+d
    return _mm256_cvtsd_f64(_mm256_add_pd(pairs, _mm256_permute4x64_pd(pairs, 0x56)));
}
#endif

int32_t AnnotatedDBG
::score_kmer_presence_mask(const sdsl::bit_vector &kmer_presence_mask,
                           int32_t match_score,
                           int32_t mismatch_score) const {
    if (!kmer_presence_mask.size())
        return 0;

    const size_t k = dbg_.get_k();
    const int32_t kmer_adjust = 3;
    const int32_t tabulate_offset = 1;

    const size_t sequence_length = kmer_presence_mask.size() + k - 1;
    const int32_t SNP_t = k + kmer_adjust;

    auto score_counter = tabulate_score(smooth_bit_vector(kmer_presence_mask),
                                        tabulate_offset);

    double score_init = std::accumulate(score_counter[1].begin(),
                                        score_counter[1].end(),
                                        0) * match_score;

    if (score_counter[0].empty())
        return score_init * sequence_length / kmer_presence_mask.size();

    if (!std::getenv("BIGSI_SCORE") && score_init == 0)
        return 0;

    const auto *it = score_counter[0].data();
    const auto *end = it + score_counter[0].size();

#ifdef __AVX2__

    __m256d match_score_packed = _mm256_set1_pd(match_score);
    __m256d mismatch_score_packed = _mm256_set1_pd(mismatch_score);
    __m256d SNP_t_packed = _mm256_set1_pd(SNP_t);
    __m256d penalties = _mm256_setzero_pd();

    for (; it + 4 <= end; it += 4) {
        __m256d penalty_add = get_penalty_bigsi_avx2(
            uint64_to_double(_mm256_load_si256(reinterpret_cast<const __m256i*>(it))),
            match_score_packed,
            mismatch_score_packed,
            SNP_t_packed
        );

        // TODO: the least-significant bits are usually off by one, is there
        //       a way to fix this?
        assert(float(add_penalty_bigsi(0, *it, match_score, mismatch_score, SNP_t)
            + add_penalty_bigsi(0, *(it + 1), match_score, mismatch_score, SNP_t)
            + add_penalty_bigsi(0, *(it + 2), match_score, mismatch_score, SNP_t)
            + add_penalty_bigsi(0, *(it + 3), match_score, mismatch_score, SNP_t))
            == float(haddall_pd(penalty_add)));

        penalties = _mm256_add_pd(penalties, penalty_add);
    }

    // reduce
    score_init += haddall_pd(penalties);

#endif

    return std::max(double(0), std::accumulate(
        it, end, score_init,
        [&](double cur, double count) {
            return add_penalty_bigsi(cur, count, match_score, mismatch_score, SNP_t);
        }
    ) * sequence_length / kmer_presence_mask.size());
}

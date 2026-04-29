#include <gtest/gtest.h>

#include <unordered_set>

#include "../graph/all/test_dbg_helpers.hpp"
#include "../test_helpers.hpp"
#include "test_annotated_dbg_helpers.hpp"

#include "annotation/coord_to_header.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/annotation_converters.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/alignment/aligner_labeled.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "kmer/alphabets.hpp"
#include "seq_io/sequence_io.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::align;
using namespace mtg::test;
using namespace mtg::kmer;

inline std::vector<std::string> get_alignment_labels(const AnnotatedDBG &anno_graph,
                                                     const Alignment &alignment,
                                                     bool check_full_coverage = true) {
    const auto &label_encoder = anno_graph.get_annotator().get_label_encoder();
    auto labels = anno_graph.get_labels(alignment.get_sequence(),
                                        check_full_coverage ? 1.0 : 0.0);
    if (check_full_coverage) {
        EXPECT_GE(labels.size(), alignment.label_columns.size());
    }

    std::unordered_set<uint64_t> enc_labels;
    for (const auto &label : labels) {
        enc_labels.emplace(label_encoder.encode(label));
    }

    std::vector<std::string> dec_labels;
    for (uint64_t label : alignment.label_columns) {
        EXPECT_TRUE(enc_labels.count(label)) << alignment;
        dec_labels.emplace_back(label_encoder.decode(label));
    }

    return dec_labels;
}

template <typename GraphAnnotationPair>
class LabeledAlignerTest : public ::testing::Test {};

typedef ::testing::Types<std::pair<DBGHashFast, annot::ColumnCompressed<>>,
                         std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
                         std::pair<DBGSSHash,   annot::ColumnCompressed<>>,
                         std::pair<DBGHashFast, annot::RowFlatAnnotator>,
                         std::pair<DBGSuccinct, annot::RowFlatAnnotator>,
                         std::pair<DBGSuccinct, annot::RowDiffColumnAnnotator>,
                         std::pair<DBGHashFast, annot::RowDiffColumnAnnotator>> FewGraphAnnotationPairTypes;

TYPED_TEST_SUITE(LabeledAlignerTest, FewGraphAnnotationPairTypes);

TYPED_TEST(LabeledAlignerTest, SimpleLinearGraph) {
    size_t k = 4;
    /*
        A    A    B    B    B    B
        GCAA-CAAT-AATG-ATGC-TGCT-GCTT
    */
    const std::vector<std::string> sequences {
        "GCAAT",
        "AATGCTT"
    };
    const std::vector<std::string> labels { "A", "B" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(k, sequences, labels);

    DBGAlignerConfig config;
    config.max_seed_length = std::numeric_limits<size_t>::max();
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
        { std::string("GCAATGCTT"), {{ { std::string("B"), std::string("AATGCTT") }, // 2S7=
                                       { std::string("A"), std::string("GCAAT") } //5=4S
                                     }} }
    }};

    for (const auto &[query, labels] : exp_alignments) {
        auto alignments = aligner.align(query);
        EXPECT_EQ(labels.size(), alignments.size()) << query;

        for (const auto &alignment : alignments) {
            bool found = false;
            for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                auto find = labels.find(label);
                ASSERT_TRUE(find != labels.end()) << label;
                if (alignment.get_sequence() == find->second) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found) << alignment;
        }
    }
}

TYPED_TEST(LabeledAlignerTest, SimpleTangleGraph) {
    size_t k = 3;
    /*  B                  AB  AB
       CGA                 GCC-CCT
          \ BC  BC  BC ABC/
           GAA-AAT-ATG-TGC
         C/               \  C   C
       GGA                 GCA-CAT
    */
    const std::vector<std::string> sequences {
        "TGCCT",
        "CGAATGCCT",
        "GGAATGCAT"
    };
    const std::vector<std::string> labels { "A", "B", "C" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(k, sequences, labels);

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
        { std::string("CGAATGCAT"), {{ { std::string("C"), std::string("GAATGCAT") }, // 1S8=
                                       { std::string("B"), std::string("CGAATGCCT") }, // 7=1X1=
                                       { std::string("A"), std::string("TGCCT") } // 4S3=1X1=
                                     }} }
    }};

    for (const auto &[query, labels] : exp_alignments) {
        auto alignments = aligner.align(query);
        EXPECT_EQ(labels.size(), alignments.size()) << query;

        for (const auto &alignment : alignments) {
            bool found = false;
            for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                auto find = labels.find(label);
                ASSERT_TRUE(find != labels.end()) << label;
                if (alignment.get_sequence() == find->second) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found) << alignment;
        }
    }
}

TYPED_TEST(LabeledAlignerTest, SimpleTangleGraphCoords) {
    // TODO: for now, not implemented for other annotators
    if constexpr(!std::is_same_v<typename TypeParam::second_type, annot::ColumnCompressed<>>
                    && !std::is_same_v<typename TypeParam::second_type, annot::RowDiffColumnAnnotator>) {
        return;
    }

    size_t k = 3;
    /*  B                  AB  AB
       CGA                 GCC-CCT
          \ BC  BC  BC ABC/
           GAA-AAT-ATG-TGC
         C/               \  C   C
       GGA                 GCA-CAT
    */
    const std::vector<std::string> sequences {
        "TGCCT",
        "CGAATGCCT",
        "GGAATGCAT"
    };
    const std::vector<std::string> labels { "A", "B", "C" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(
        k, sequences, labels, DeBruijnGraph::BASIC, true
    );

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

    std::unordered_map<std::string, std::unordered_map<std::string, std::pair<std::string, int32_t>>> exp_alignments {{
        { std::string("CGAATGCAT"), {{ { std::string("C"), std::make_pair(std::string("GAATGCAT"), 1) }, // 1S8=
                                       { std::string("B"), std::make_pair(std::string("CGAATGCCT"), 0) }, // 7=1X1=
                                       { std::string("A"), std::make_pair(std::string("TGCCT"), 0) } // 4S3=1X1=
                                     }} }
    }};

    for (const auto &[query, labels] : exp_alignments) {
        auto alignments = aligner.align(query);
        EXPECT_EQ(labels.size(), alignments.size()) << query;

        for (const auto &alignment : alignments) {
            bool found = false;
            ASSERT_EQ(alignment.label_columns.size(), alignment.label_coordinates.size());
            size_t label_index = 0;
            for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                ASSERT_GT(alignment.label_coordinates[label_index].size(), 0);
                auto find = labels.find(label);
                ASSERT_TRUE(find != labels.end()) << label;
                if (alignment.get_sequence() == find->second.first) {
                    found = true;
                    EXPECT_EQ(find->second.second,
                              alignment.label_coordinates[label_index][0]);
                    break;
                }
                ++label_index;
            }
            EXPECT_TRUE(found) << alignment;
        }
    }
}

TYPED_TEST(LabeledAlignerTest, SimpleTangleGraphCoordsMiddle) {
    // TODO: for now, not implemented for other annotators
    if constexpr(!std::is_same_v<typename TypeParam::second_type, annot::ColumnCompressed<>>
                    && !std::is_same_v<typename TypeParam::second_type, annot::RowDiffColumnAnnotator>) {
        return;
    }

    size_t k = 3;
    /*  B                  AB  AB
       CGA                 GCC-CCT
          \ BC  BC  BC ABC/
           GAA-AAT-ATG-TGC
         C/               \  C   C
       GGA                 GCA-CAT
    */
    const std::vector<std::string> sequences {
        "TGCCT",
        "CGAATGCCT",
        "GGAATGCAT"
    };
    const std::vector<std::string> labels { "A", "B", "C" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(
        k, sequences, labels, DeBruijnGraph::BASIC, true
    );

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

    std::unordered_map<std::string, std::unordered_map<std::string, std::pair<std::string, int32_t>>> exp_alignments {{
        { std::string("CGAAAGCCT"), {{ { std::string("C"), std::make_pair(std::string("GAATGCAT"), 1) }, // 1S3=1X2=1X1=
                                       { std::string("B"), std::make_pair(std::string("CGAATGCCT"), 0) }, // 4=1X4=
                                       { std::string("A"), std::make_pair(std::string("GCCT"), 1) } // 5S4=
                                     }} }
    }};

    for (const auto &[query, labels] : exp_alignments) {
        auto alignments = aligner.align(query);
        EXPECT_EQ(labels.size(), alignments.size()) << query;

        for (const auto &alignment : alignments) {
            bool found = false;
            ASSERT_EQ(alignment.label_columns.size(), alignment.label_coordinates.size());
            size_t label_index = 0;
            for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                ASSERT_GT(alignment.label_coordinates[label_index].size(), 0);
                auto find = labels.find(label);
                ASSERT_TRUE(find != labels.end()) << label;
                if (alignment.get_sequence() == find->second.first) {
                    found = true;
                    EXPECT_EQ(find->second.second,
                              alignment.label_coordinates[label_index][0]);
                    break;
                }
                ++label_index;
            }
            EXPECT_TRUE(found) << alignment;
        }
    }
}

// Fixture for the end-to-end CoordToHeader alignment tests. Builds a
// DBGSuccinct + ColumnCompressed annotation where all sequences share a
// single label (simulating --anno-filename), with distinct coordinate
// offsets so the k-mer coord axis spans the full column.
//
// Guarded for DNA-only: these tests use the DNA scoring matrix and short
// ACGT sequences that exercise corner cases in the chainer; some of those
// trip Debug-build assertions when run against the protein DBG (which has
// a different alphabet and seed-extension behavior). The functionality
// under test (CoordToHeader formatting) is alphabet-agnostic.
#if ! _PROTEIN_GRAPH
class LabeledAlignerCoordTest : public ::testing::Test {
  protected:
    static constexpr size_t k = 5;
    std::unique_ptr<AnnotatedDBG> anno_graph;
    DBGAlignerConfig config;

    void SetUp() override {
        config.max_seed_length = std::numeric_limits<size_t>::max();
        config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    }

    void build(const std::vector<std::string> &sequences,
               const std::vector<uint64_t> &coord_starts) {
        std::vector<std::string> labels(sequences.size(), "file");
        anno_graph = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            k, sequences, labels, DeBruijnGraph::BASIC, /*coordinates=*/true, coord_starts
        );
    }

    AlignmentResults align(const std::string &query) {
        LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());
        return aligner.align(query);
    }

    // Attach the CoordToHeader + k to every alignment in `paths`.
    // `cth` must outlive `paths` (typically declared in the test scope).
    void attach(AlignmentResults &paths, const annot::CoordToHeader &cth) {
        for (auto &aln : paths) {
            aln.coord_to_header = &cth;
            aln.coord_to_header_k = k;
        }
    }
};

TEST_F(LabeledAlignerCoordTest, CrossBoundary) {
    // Two sequences sharing the boundary 4-mer 'ACGT'. Literal seq1+seq2
    // duplicates the shared 'ACGT', but the graph path does not (spells
    // seq1 + seq2[k-1:] = 20 nt). The aligner soft-clips the first 8
    // query bases and takes the higher-scoring 16 nt cross-boundary run
    // starting at seq1 k-mer 4 (8S16= vs the alternative 12=12S).
    build({ "AAAAACGTACGT", "ACGTTTTTTTTT" }, { 0, 8 });

    auto alignments = align("AAAAACGTACGTACGTTTTTTTTT");
//                           SSSSSSSS================
    ASSERT_EQ(1u, alignments.size());
    const auto &aln = alignments[0];
    EXPECT_EQ("8S16=", aln.get_cigar().to_string());
    ASSERT_EQ(1u, aln.label_coordinates.size());
    ASSERT_EQ(1u, aln.label_coordinates[0].size());
    EXPECT_EQ(4u, aln.label_coordinates[0][0]);

    annot::CoordToHeader cth({ { "seq1", "seq2" } }, { { 8, 8 } });
    attach(alignments, cth);
    EXPECT_EQ("seq1:5-12;seq2:1-8", alignments[0].format_coords());
}

TEST_F(LabeledAlignerCoordTest, CrossBoundaryThreeSequences) {
    // Three sequences chained by shared boundary 4-mers (ACGT, AGCT).
    // Query spells the 24 nt cross-boundary path from seq1 k-mer 4
    // through all of seq2 into the first 4 nt of seq3.
    build({
        "GCGCTTCGACGT",   // ends 'ACGT'
        "ACGTATGCAGCT",   // starts 'ACGT', ends 'AGCT'
        "AGCTGGGCGCAT",   // starts 'AGCT'
    }, { 0, 8, 16 });

    auto alignments = align("TTCGACGTATGCAGCTGGGCGCAT");
    ASSERT_EQ(1u, alignments.size());
    const auto &aln = alignments[0];
    EXPECT_EQ("24=", aln.get_cigar().to_string());
    ASSERT_EQ(1u, aln.label_coordinates.size());
    ASSERT_EQ(1u, aln.label_coordinates[0].size());
    EXPECT_EQ(4u, aln.label_coordinates[0][0]);

    annot::CoordToHeader cth({ { "seq1", "seq2", "seq3" } }, { { 8, 8, 8 } });
    attach(alignments, cth);
    EXPECT_EQ("seq1:5-12;seq2:1-12;seq3:1-4", alignments[0].format_coords());
}

TEST_F(LabeledAlignerCoordTest, CrossBoundaryThreeSequencesWithIndels) {
    // Same graph as the ThreeSequences test. Query has 'AA' inserted and
    // one 'G' deleted vs the 24-nt ref path; aligner deterministically
    // picks CIGAR '13=2I3=1D7='. The ref-path spelling (and per-header
    // split) is identical to the non-indel case.
    build({
        "GCGCTTCGACGT",
        "ACGTATGCAGCT",
        "AGCTGGGCGCAT",
    }, { 0, 8, 16 });
    config.min_exact_match = 0.0;

    auto alignments = align("TTCGACGTATGCAAAGCT"/*G*/"GGCGCAT");
//                           =============II===   D   =======
    ASSERT_EQ(1u, alignments.size());
    const auto &aln = alignments[0];
    EXPECT_EQ("13=2I3=1D7=", aln.get_cigar().to_string());
    EXPECT_EQ("TTCGACGTATGCAGCTGGGCGCAT", aln.get_sequence());
    ASSERT_EQ(1u, aln.label_coordinates.size());
    ASSERT_EQ(1u, aln.label_coordinates[0].size());
    EXPECT_EQ(4u, aln.label_coordinates[0][0]);

    annot::CoordToHeader cth({ { "seq1", "seq2", "seq3" } }, { { 8, 8, 8 } });
    attach(alignments, cth);
    EXPECT_EQ("seq1:5-12;seq2:1-12;seq3:1-4", alignments[0].format_coords());
}

TEST_F(LabeledAlignerCoordTest, ThreeSequencesPartialCoverage) {
    // Three sequences share the 7-mer 'CGTACGT' in their middles. A query
    // matching just the shared 7-mer yields one alignment whose coord
    // tuple holds three starts — one per sequence — so format_coords
    // reports a partial range in the middle of each.
    build({
        "TTCGTACGTAAA",
        "AACGTACGTCCC",
        "GGCGTACGTGGG",
    }, { 0, 8, 16 });

    auto alignments = align("CGTACGT");
    ASSERT_EQ(1u, alignments.size());
    const auto &aln = alignments[0];
    EXPECT_EQ("7=", aln.get_cigar().to_string());
    ASSERT_EQ(1u, aln.label_coordinates.size());
    ASSERT_EQ(3u, aln.label_coordinates[0].size());
    EXPECT_EQ(2u,  aln.label_coordinates[0][0]);
    EXPECT_EQ(10u, aln.label_coordinates[0][1]);
    EXPECT_EQ(18u, aln.label_coordinates[0][2]);

    annot::CoordToHeader cth({ { "seq1", "seq2", "seq3" } }, { { 8, 8, 8 } });
    attach(alignments, cth);
    EXPECT_EQ("seq1:3-9;seq2:3-9;seq3:3-9", alignments[0].format_coords());
}

TEST_F(LabeledAlignerCoordTest, CrossBoundaryReverseComplement) {
    // rc of the CrossBoundary query. label_coordinates and sequence_
    // stay in fwd-strand space, so format_coords emits the same ranges
    // as the non-rc case.
    build({ "AAAAACGTACGT", "ACGTTTTTTTTT" }, { 0, 8 });

    auto alignments = align("AAAAAAAAACGTACGTACGTTTTT");
    ASSERT_EQ(1u, alignments.size());
    const auto &aln = alignments[0];
    EXPECT_TRUE(aln.get_orientation());
    EXPECT_EQ("8S16=", aln.get_cigar().to_string());
    ASSERT_EQ(1u, aln.label_coordinates.size());
    ASSERT_EQ(1u, aln.label_coordinates[0].size());
    EXPECT_EQ(4u, aln.label_coordinates[0][0]);

    annot::CoordToHeader cth({ { "seq1", "seq2" } }, { { 8, 8 } });
    attach(alignments, cth);
    EXPECT_EQ("seq1:5-12;seq2:1-8", alignments[0].format_coords());
}

TEST_F(LabeledAlignerCoordTest, CrossBoundaryWithIndel) {
    // Same graph as CrossBoundary. Query has extra bases inserted near
    // the seq1/seq2 boundary; aligner deterministically picks CIGAR
    // '12=4I8=1S' — ref-path spelling length = 20 nt.
    build({ "AAAAACGTACGT", "ACGTTTTTTTTT" }, { 0, 8 });
    config.min_exact_match = 0.0;

    auto alignments = align("AAAAACGTACGTCACGTTTTTTTTT");
//                           ============IIII========S
    ASSERT_EQ(1u, alignments.size());
    const auto &aln = alignments[0];
    EXPECT_EQ("12=4I8=1S", aln.get_cigar().to_string());
    EXPECT_EQ(20u, aln.get_sequence().size());
    ASSERT_EQ(1u, aln.label_coordinates.size());
    ASSERT_EQ(1u, aln.label_coordinates[0].size());
    EXPECT_EQ(0u, aln.label_coordinates[0][0]);

    annot::CoordToHeader cth({ { "seq1", "seq2" } }, { { 8, 8 } });
    attach(alignments, cth);
    EXPECT_EQ("seq1:1-12;seq2:1-8", alignments[0].format_coords());
}
#endif  // ! _PROTEIN_GRAPH

// Protein-alphabet equivalents of the basic CoordToHeader cases. Use BLOSUM62
// scoring and amino-acid sequences with letters not present in DNA so the
// chainer's seed pattern differs from the DNA fixture's corner cases.
#if _PROTEIN_GRAPH
class LabeledAlignerProteinCoordTest : public ::testing::Test {
  protected:
    static constexpr size_t k = 5;
    std::unique_ptr<AnnotatedDBG> anno_graph;
    DBGAlignerConfig config;

    void SetUp() override {
        config.max_seed_length = std::numeric_limits<size_t>::max();
        config.score_matrix = DBGAlignerConfig::score_matrix_blosum62;
    }

    void build(const std::vector<std::string> &sequences,
               const std::vector<uint64_t> &coord_starts) {
        std::vector<std::string> labels(sequences.size(), "file");
        anno_graph = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            k, sequences, labels, DeBruijnGraph::BASIC, /*coordinates=*/true, coord_starts
        );
    }

    AlignmentResults align(const std::string &query) {
        LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());
        return aligner.align(query);
    }

    void attach(AlignmentResults &paths, const annot::CoordToHeader &cth) {
        for (auto &aln : paths) {
            aln.coord_to_header = &cth;
            aln.coord_to_header_k = k;
        }
    }
};

TEST_F(LabeledAlignerProteinCoordTest, CrossBoundary) {
    // Two amino-acid sequences sharing the boundary 4-mer 'MNPQ'. Exact-match
    // query spans the seq1/seq2 boundary; format_coords splits it across
    // both per-sequence ranges.
    build({ "EEEEEEEEMNPQ", "MNPQRRRRRRRR" }, { 0, 8 });

    auto alignments = align("EEEEMNPQRRRR");
    ASSERT_EQ(1u, alignments.size());
    const auto &aln = alignments[0];
    EXPECT_EQ("12=", aln.get_cigar().to_string());
    ASSERT_EQ(1u, aln.label_coordinates.size());
    ASSERT_EQ(1u, aln.label_coordinates[0].size());
    EXPECT_EQ(4u, aln.label_coordinates[0][0]);

    annot::CoordToHeader cth({ { "seq1", "seq2" } }, { { 8, 8 } });
    attach(alignments, cth);
    EXPECT_EQ("seq1:5-12;seq2:1-4", alignments[0].format_coords());
}

TEST_F(LabeledAlignerProteinCoordTest, SharedKmerMultipleLabels) {
    // Two amino-acid sequences sharing an internal 5-mer 'MNPQR'. Query of
    // exactly that 5-mer matches both sequences; the resulting alignment's
    // coord tuple has two entries (one per occurrence) and format_coords
    // emits a partial range for each header.
    build({ "EEEMNPQRSSSS", "WWWMNPQRTTTT" }, { 0, 8 });

    auto alignments = align("MNPQR");
    ASSERT_EQ(1u, alignments.size());
    const auto &aln = alignments[0];
    EXPECT_EQ("5=", aln.get_cigar().to_string());
    ASSERT_EQ(1u, aln.label_coordinates.size());
    ASSERT_EQ(2u, aln.label_coordinates[0].size());

    annot::CoordToHeader cth({ { "seq1", "seq2" } }, { { 8, 8 } });
    attach(alignments, cth);
    EXPECT_EQ("seq1:4-8;seq2:4-8", alignments[0].format_coords());
}
#endif  // _PROTEIN_GRAPH

TYPED_TEST(LabeledAlignerTest, SimpleTangleGraphCoordsCycle) {
    // TODO: for now, not implemented for other annotators
    if constexpr(!std::is_same_v<typename TypeParam::second_type, annot::ColumnCompressed<>>
                    && !std::is_same_v<typename TypeParam::second_type, annot::RowDiffColumnAnnotator>) {
        return;
    }

    size_t k = 4;
    /*
        A    A    B    B    B    B    B
        GCAA-CAAT-AATG-ATGC-TGCG-GCGC-CGCA
    */
    const std::vector<std::string> sequences {
        "GCAAT",
        "AATGCGCA"
    };
    const std::vector<std::string> labels { "A", "B" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(
        k, sequences, labels, DeBruijnGraph::BASIC, true
    );

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

    std::unordered_map<std::string, std::unordered_map<std::string, std::pair<std::string, int32_t>>> exp_alignments {{
        { std::string("ATGCGCAATGCG"), {{ { std::string("B"), std::make_pair(std::string("ATGCGCA"), 1) },
                                          { std::string("A"), std::make_pair(std::string("GCAAT"), 0) },
                                     }} }
    }};

    for (const auto &[query, labels] : exp_alignments) {
        auto alignments = aligner.align(query);
        EXPECT_EQ(labels.size(), alignments.size()) << query;

        for (const auto &alignment : alignments) {
            bool found = false;
            ASSERT_EQ(alignment.label_columns.size(), alignment.label_coordinates.size());
            size_t label_index = 0;
            for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                ASSERT_GT(alignment.label_coordinates[label_index].size(), 0);
                auto find = labels.find(label);
                ASSERT_TRUE(find != labels.end()) << label;
                if (alignment.get_sequence() == find->second.first) {
                    found = true;
                    EXPECT_EQ(find->second.second,
                              alignment.label_coordinates[label_index][0]);
                    break;
                }
                ++label_index;
            }
            EXPECT_TRUE(found) << alignment;
        }
    }
}

TEST(LabeledAlignerTest, SimpleTangleGraphSuffixSeed) {
    size_t k = 4;
    /*  B    B                  AB   AB
       TCGA-CGAA                TGCC-GCCT
                \ BC   BC   BC /
                 GAAT-AATG-ATGC
          C    C/              \   C    C
       TGGA-GGAA                TGCA-GCAT
    */
    const std::vector<std::string> sequences {
        "TGCCT",
        "TCGAATGCCT",
        "TGGAATGCAT"
    };
    const std::vector<std::string> labels { "A", "B", "C" };

    auto anno_graph = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(k, sequences, labels);

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    config.min_seed_length = 2;
    config.left_end_bonus = 5;
    config.right_end_bonus = 5;
    LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
        { std::string("TGAAATGCAT"), {{
#if ! _PROTEIN_GRAPH
            { std::string("C"), std::string("TGGAATGCAT") }, // 2=1X7=
            { std::string("B"), std::string("TCGAATGCCT") } // 1=2X5=1X1=
#else
            { std::string("C"), std::string("AATGCAT") }, // 3S7=
            { std::string("B"), std::string("AATGCCT") } // 3S5=1X1=
#endif
        }} }
    }};

    for (const auto &[query, labels] : exp_alignments) {
        auto alignments = aligner.align(query);
        EXPECT_EQ(labels.size(), alignments.size()) << query;

        for (const auto &alignment : alignments) {
            bool found = false;
            for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                auto find = labels.find(label);
                ASSERT_TRUE(find != labels.end());
                if (alignment.get_sequence() == find->second) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found) << alignment;
        }
    }
}

#if ! _PROTEIN_GRAPH
TYPED_TEST(LabeledAlignerTest, CanonicalTangleGraph) {
    size_t k = 5;
    /*   B     B                AB    AB
       TTAGT-TAGTC             TCGAA-CGAAA
                  \  BC   ABC /
                   AGTCG-GTCGA
          C     C /           \   C     C
       TCAGT-CAGTC             TCGAT-CGATT
        AB    AB                 B     B
       TTTCG-TTCGA             GACTA-ACTAA
                  \ ABC    BC /
                   TCGAC-CGACT
          C     C /           \   C     C
       AATCG-ATCGA             GACTG-ACTGA
    */
    const std::vector<std::string> sequences {
        "GTCGAAA", // "TTTCGAC"
        "TTAGTCGAAA", // "TTTCGACTAA"
        "TCAGTCGATT" // "AATCGACTGA"
    };
    const std::vector<std::string> labels { "A", "B", "C" };

    for (DeBruijnGraph::Mode mode : { DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                           typename TypeParam::second_type>(
            k, sequences, labels, mode
        );

        DBGAlignerConfig config;
        config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
        LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

        std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
            // r.c. TTTGAACTAA
            { std::string("TTAGTTCAAA"), {{ { std::string("B"), std::string("TTAGTCGAAA") } }} } // 5=2X3=
        }};

        for (const auto &[query, labels] : exp_alignments) {
            auto alignments = aligner.align(query);
            EXPECT_EQ(labels.size(), alignments.size()) << query;

            for (const auto &alignment : alignments) {
                bool found = false;
                for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                    auto find = labels.find(label);
                    ASSERT_TRUE(find != labels.end());
                    if (alignment.get_sequence() == find->second) {
                        found = true;
                        break;
                    }
                }
                EXPECT_TRUE(found) << alignment;
            }
        }
    }
}
#endif

} // namespace

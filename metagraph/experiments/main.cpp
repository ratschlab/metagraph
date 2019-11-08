#include <filesystem>
#include <thread>
#include <chrono>
#include <tclap/CmdLine.h>
#include <sdsl/rrr_vector.hpp>
#include <ips4o.hpp>
#include <libmaus2/util/NumberSerialisation.hpp>
#include <boost/multiprecision/integer.hpp>
#include <progress_bar.hpp>

#include "method_constructors.hpp"
#include "data_generation.hpp"
#include "annotate_static.hpp"
#include "annotate_column_compressed.hpp"
#include "unix_tools.hpp"
#include "string_utils.hpp"
#include "kmc_parser.hpp"
#include "alphabets.hpp"
#include "aligner_helper.hpp"
#include "reverse_complement.hpp"
#include "annotate_row_compressed.hpp"
#include "static_annotators_def.hpp"
#include "config.hpp"
#include "serialization.hpp"
#include "wavelet_tree.hpp"

using namespace std::chrono_literals;

using TCLAP::ValueArg;
using TCLAP::MultiArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::UnlabeledMultiArg;
using TCLAP::ValuesConstraint;


//TODO
// THINGS TO REPORT:
//  - querying time
//  - RAM usage
//  - disk usage
//  - construction time

template <class BitVector>
void test_vector_points(uint64_t n, double d, const std::string &prefix) {
    DataGenerator generator;
    generator.set_seed(42);

    auto other = generator.generate_random_column(n, d)->convert_to<BitVector>();

    if (static_cast<int>(other.rank1(1)) < -1
            || static_cast<int>(other.rank0(1)) < -1)
        throw std::runtime_error("Never happends, just initializing the rank support");

    std::filesystem::path path(std::string("test.")
                                + prefix
                                + "." + std::to_string(n) + "_" + std::to_string(d)
                                + ".bv");
    std::ofstream out(path, std::ios::binary);

    other.serialize(out);
    out.close();

    const auto serialized_size = std::filesystem::file_size(path);

    auto mem_before = get_curr_RSS();

    BitVector another;

    std::this_thread::sleep_for(5s);

    std::ifstream in(path, std::ios::binary);
    another.load(in);
    in.close();

    if (static_cast<int>(another.rank1(another.size() / 2)) < -1
            || static_cast<int>(another.rank0(another.size() / 2)) < -1)
        throw std::runtime_error("Never happends, just initializing the rank support");
    if (another.num_set_bits() && static_cast<int>(another.select1(1)) < -1)
        throw std::runtime_error("Never happends, just initializing the select support");

    auto RAM = get_curr_RSS() - mem_before;

    std::filesystem::remove(path);

    uint64_t predicted_size;
    if constexpr(!std::is_base_of<bit_vector_dyn, BitVector>::value) {
        predicted_size = predict_size<BitVector>(another.size(), another.num_set_bits());
    } else {
        predicted_size = 0;
    }

    std::cout << prefix
              << "\t" << n
              << "\t" << d
              << "\t" << 1. * another.num_set_bits() / another.size()
              << "\t" << 1. * serialized_size * 8 / n
              << "\t" << RAM
              << "\t" << 1. * predicted_size / another.size()
              << std::endl;
}

std::vector<double> get_densities(uint64_t num_cols,
                                  const std::vector<double> &vector) {
    auto densities = vector;
    if (densities.size() == 1) {
        densities.assign(num_cols, densities[0]);
    } else if (densities.size() != num_cols) {
        std::cout << "ERROR: wrong number of column counts" << std::endl;
        exit(1);
    }
    return densities;
}

double test_point_time(const BinaryMatrix &matrix,
                       uint64_t num_samples = 1000) {
    DataGenerator generator;
    generator.set_seed(42);
    //std::cout << "Generating " << num_samples << " samples";
    auto row_positions = generator.generate_random_ints(
        num_samples, 0, matrix.num_rows()
    );
    auto col_positions = generator.generate_random_ints(
        num_samples, 0, matrix.num_columns()
    );

    Timer timer;
    timer.reset();
    for (uint64_t i = 0; i < num_samples; ++i) {
        matrix.get(row_positions[i], col_positions[i]);
    }

    return timer.elapsed() * 1000 / num_samples;
}

double test_row_time(const BinaryMatrix &matrix,
                     uint64_t num_samples = 1000) {
    DataGenerator generator;
    generator.set_seed(42);
    //std::cout << "Generating " << num_samples << " samples";
    auto positions = generator.generate_random_ints(
        num_samples, 0, matrix.num_rows());

    Timer timer;
    timer.reset();
    for (auto pos : positions) {
        matrix.get_row(pos);
    }

    return timer.elapsed() * 1000 / num_samples;
}

double test_column_time(const BinaryMatrix &matrix,
                        uint64_t num_samples = 1000) {
    DataGenerator generator;
    generator.set_seed(42);
    //std::cout << "Generating " << num_samples << " samples";
    auto positions = generator.generate_random_ints(
        num_samples, 0, matrix.num_columns());

    Timer timer;
    timer.reset();
    for (auto pos : positions) {
        matrix.get_column(pos);
    }

    return timer.elapsed() * 1000 / num_samples;
}

void dump_column_slice(const bit_vector &column,
                       double begin,
                       double end,
                       const std::string &outfile) {
    if (begin >= end)
        throw std::runtime_error("Empty range");

    size_t begin_ind = begin * column.size();
    size_t end_ind = end * column.size();

    std::ofstream out(outfile, std::ios::binary);
    out << column.size() << " "
        << column.rank1(end_ind - 1) - (begin_ind ? column.rank1(begin_ind - 1) : 0)
        << "\n";

    column.call_ones_in_range(begin * column.size(),
                              end * column.size(),
                              [&](auto i) { out << i << "\n"; });
}

int cleaning_pick_kmer_threshold(const uint64_t *kmer_covg, size_t arrlen,
                                 double *alpha_est_ptr, double *beta_est_ptr,
                                 double *false_pos_ptr, double *false_neg_ptr);


Config::AnnotationType parse_annotation_type(const std::string &filename) {
    if (utils::ends_with(filename, annotate::kColumnAnnotatorExtension)) {
        return Config::AnnotationType::ColumnCompressed;

    } else if (utils::ends_with(filename, annotate::kRowAnnotatorExtension)) {
        return Config::AnnotationType::RowCompressed;

    } else if (utils::ends_with(filename, annotate::kBRWTExtension)) {
        return Config::AnnotationType::BRWT;

    } else if (utils::ends_with(filename, annotate::kBinRelWT_sdslExtension)) {
        return Config::AnnotationType::BinRelWT_sdsl;

    } else if (utils::ends_with(filename, annotate::kBinRelWTExtension)) {
        return Config::AnnotationType::BinRelWT;

    } else if (utils::ends_with(filename, annotate::kRowPackedExtension)) {
        return Config::AnnotationType::RowFlat;

    } else if (utils::ends_with(filename, annotate::kRainbowfishExtension)) {
        return Config::AnnotationType::RBFish;

    } else {
        std::cerr << "Error: unknown annotation format in "
                  << filename << std::endl;
        exit(1);
    }
}

typedef annotate::MultiLabelEncoded<uint64_t, std::string> Annotator;

std::unique_ptr<Annotator> initialize_annotation(const std::string &filename,
                                                 size_t cache_size = 0) {
    std::unique_ptr<Annotator> annotation;

    switch (parse_annotation_type(filename)) {
        case Config::ColumnCompressed: {
            annotation.reset(new annotate::ColumnCompressed<>(1));
            break;
        }
        case Config::RowCompressed: {
            annotation.reset(new annotate::RowCompressed<>(1));
            break;
        }
        case Config::BRWT: {
            annotation.reset(new annotate::BRWTCompressed<>(cache_size));
            break;
        }
        case Config::BinRelWT_sdsl: {
            annotation.reset(new annotate::BinRelWT_sdslAnnotator(cache_size));
            break;
        }
        case Config::BinRelWT: {
            annotation.reset(new annotate::BinRelWTAnnotator(cache_size));
            break;
        }
        case Config::RowFlat: {
            annotation.reset(new annotate::RowFlatAnnotator(cache_size));
            break;
        }
        case Config::RBFish: {
            annotation.reset(new annotate::RainbowfishAnnotator(cache_size));
            break;
        }
    }

    if (annotation->load(filename))
        return annotation;

    std::cerr << "ERROR: can't load annotation "
              << filename << std::endl;
    exit(1);
}


int main(int argc, char *argv[]) {
    try {
        TCLAP::CmdLine cmd("Benchmarks and data generation for metagraph", ' ', "");

        std::vector<std::string> regimes {
            "vectors",
            "matrices",
            "subsets",
            "to_rrr",
            "stats",
            "query",
            "slice",
            "estimate_abundance_threshold",
            "evaluate_alignment",
            "query_annotation_rows",
            "to_dna4"
        };

        ValuesConstraint<std::string> regime_constraint(regimes);
        UnlabeledValueArg<std::string> regime_arg("regime", "Regime", true, "", &regime_constraint, cmd);
        cmd.parse(std::min(argc, 2), argv);

        std::string regime = regime_arg.getValue();
        regime_arg.reset();

        if (regime == "vectors") {
            std::vector<std::string> vector_types {
                "small",
                "smart",
                "stat",
                "sd",
                "rrr63",
                "rrr127",
                "rrr255",
                "dyn",
            };
            ValuesConstraint<std::string> vector_type_constraint(vector_types);
            UnlabeledValueArg<std::string> vector_type_arg("vector_types", "Vector type", true, "", &vector_type_constraint, cmd);

            cmd.parse(std::min(argc, 3), argv);
            std::string vector_type = vector_type_arg.getValue();
            cmd.reset();

            ValueArg<size_t> length_arg("l", "length", "length of bit vector", true, 0, "int", cmd);
            ValueArg<double> density_arg("d", "density", "Density of bit vector", true, 1., "double", cmd);

            cmd.parse(argc, argv);

            if (vector_type == "small") {
                test_vector_points<bit_vector_small>(length_arg.getValue(),
                                                     density_arg.getValue(),
                                                     "small");
            } else if (vector_type == "smart") {
                test_vector_points<bit_vector_smart>(length_arg.getValue(),
                                                     density_arg.getValue(),
                                                     "smart");
            } else if (vector_type == "stat") {
                test_vector_points<bit_vector_stat>(length_arg.getValue(),
                                                    density_arg.getValue(),
                                                    "stat");
            } else if (vector_type == "sd") {
                test_vector_points<bit_vector_sd>(length_arg.getValue(),
                                                  density_arg.getValue(),
                                                  "sd");
            } else if (vector_type == "rrr63") {
                test_vector_points<bit_vector_rrr<>>(length_arg.getValue(),
                                                     density_arg.getValue(),
                                                     "rrr");
            } else if (vector_type == "rrr127") {
                test_vector_points<bit_vector_rrr<127>>(length_arg.getValue(),
                                                     density_arg.getValue(),
                                                     "rrr127");
            } else if (vector_type == "rrr255") {
                test_vector_points<bit_vector_rrr<255>>(length_arg.getValue(),
                                                     density_arg.getValue(),
                                                     "rrr255");
            } else if (vector_type == "dyn") {
                test_vector_points<bit_vector_dyn>(length_arg.getValue(),
                                                   density_arg.getValue(),
                                                   "dyn");
            }
        } else if (regime == "matrices") {
            std::vector<std::string> matrix_regimes {
                "simulate",
            };
            ValuesConstraint<std::string> matrix_regime_constraint(matrix_regimes);
            UnlabeledValueArg<std::string> matrix_regime_arg("matrix_regime", "Matrix regime", true, "", &matrix_regime_constraint, cmd);

            cmd.parse(std::min(argc, 3), argv);
            std::string matrix_regime = matrix_regime_arg.getValue();
            cmd.reset();

            std::vector<std::string> modes {
                "norepl",
                "uniform_rows",
                "weighted_rows",
                "uniform_columns"
            };
            ValuesConstraint<std::string> mode_constraint(modes);
            UnlabeledValueArg<std::string> mode_arg("mode", "Row replication mode", true, "", &mode_constraint, cmd);
            cmd.parse(std::min(argc, 4), argv);
            std::string mode = mode_arg.getValue();
            cmd.reset();

            std::vector<std::string> compressors {
                "row",
                "flat",
                "column",
                "brwt",
                "bin_rel_wt_sdsl",
                "bin_rel_wt",
                "rbfish"
            };
            ValuesConstraint<std::string> compresser_constraint(compressors);
            ValueArg<std::string> compressor_arg("", "anno-type", "Compression method to use", true, "brwt", &compresser_constraint, cmd);
            cmd.parse(std::min(argc, 6), argv);
            std::string compressor_name = compressor_arg.getValue();
            cmd.reset();

            std::unique_ptr<ValueArg<size_t>> arity_arg;
            std::unique_ptr<ValueArg<bool>> greedy_arg;
            std::unique_ptr<ValueArg<size_t>> relax_arg;
            if (compressor_name == "brwt") {
                arity_arg.reset(new ValueArg<size_t>("", "arity", "number of children per brwt node", false, 2, "int", cmd));
                greedy_arg.reset(new ValueArg<bool>("", "greedy", "use greedy binary partitions", false, false, "bool", cmd));
                relax_arg.reset(new ValueArg<size_t>("", "relax", "max number of children in relaxed nodes (1 -- no relax)", false, 1, "int", cmd));
            }
            ValueArg<size_t> rows_arg("n", "rows", "number of rows", true, 0, "int", cmd);
            ValueArg<size_t> cols_arg("m", "columns", "number of columns", true, 0, "int", cmd);
            MultiArg<double> density_arg("d", "density", "Density of each column", true, "double", cmd);

            std::string out_name = std::string("simulate.") + mode;

            DataGenerator generator;
            generator.set_seed(42);

            //std::cout << "Generating" << std::endl;
            std::vector<std::unique_ptr<bit_vector>> generated_columns;

            if (mode == "norepl") {
                cmd.parse(argc, argv);

                generated_columns = generator.generate_random_columns(
                    rows_arg.getValue(),
                    cols_arg.getValue(),
                    get_densities(cols_arg.getValue(), density_arg.getValue())
                );
                out_name += std::string("_")
                    + std::to_string(density_arg.getValue()[0]);
            } else if (mode == "uniform_rows") {
                ValueArg<size_t> unique_arg("u", "unique_rows", "Number of unique rows", true, 1, "int", cmd);
                cmd.parse(argc, argv);

                generated_columns = generator.generate_random_rows(
                    unique_arg.getValue(),
                    cols_arg.getValue(),
                    get_densities(cols_arg.getValue(), density_arg.getValue()),
                    std::vector<uint32_t>(unique_arg.getValue(),
                                          rows_arg.getValue() / unique_arg.getValue())
                );
                assert(generated_columns.size() == cols_arg.getValue());
                assert(generated_columns[0]->size() == rows_arg.getValue());
                out_name += std::string("_")
                    + std::to_string(density_arg.getValue()[0]) + "_"
                    + std::to_string(unique_arg.getValue());

            } else if (mode == "uniform_columns") {
                ValueArg<size_t> unique_arg("u", "unique_columns", "Number of unique columns", true, 1, "int", cmd);
                cmd.parse(argc, argv);

                generated_columns = generator.generate_random_columns(
                    rows_arg.getValue(),
                    unique_arg.getValue(),
                    get_densities(unique_arg.getValue(), density_arg.getValue()),
                    std::vector<uint32_t>(unique_arg.getValue(),
                                          cols_arg.getValue() / unique_arg.getValue())
                );
                assert(generated_columns.size() == cols_arg.getValue());
                assert(generated_columns[0]->size() == rows_arg.getValue());
                out_name += std::string("_")
                    + std::to_string(density_arg.getValue()[0]) + "_"
                    + std::to_string(unique_arg.getValue());
            } else if (mode == "weighted_rows") {
                ValueArg<size_t> unique_arg("u", "unique_rows", "Number of unique rows", true, 1, "int", cmd);
                ValueArg<size_t> count_arg("c", "init_count", "Multiplicity of the most frequent row", true, 2, "int", cmd);
                ValueArg<double> decay_arg("D", "decay", "Decay rate of row frequency curve", false, 1., "double", cmd);
                cmd.parse(argc, argv);

                generated_columns = generator.generate_random_rows(
                    unique_arg.getValue(),
                    cols_arg.getValue(),
                    get_densities(cols_arg.getValue(), density_arg.getValue()),
                        generator.generate_row_counts_index_inverse(
                            unique_arg.getValue(),
                            count_arg.getValue(),
                            rows_arg.getValue()
                        )
                );
                out_name += std::string("_")
                    + std::to_string(density_arg.getValue()[0]) + "_"
                    + std::to_string(unique_arg.getValue()) + "_"
                    + std::to_string(count_arg.getValue()) + "_"
                    + std::to_string(decay_arg.getValue());
            } else {
                throw std::runtime_error("Invalid generation mode");
            }

            const auto num_cols = generated_columns.size();

            if (num_cols == 0) {
                std::cerr << "ERROR: Empty input matrix!" << std::endl;
                exit(1);
            }

            out_name += std::string("_")
                + std::to_string(generated_columns[0]->size()) + "_"
                + std::to_string(num_cols) + "_"
                + compressor_name;

            std::vector<double> observed_densities;
            for (const auto &column_ptr : generated_columns) {
                observed_densities.push_back(
                    static_cast<double>(column_ptr->num_set_bits())
                                            / column_ptr->size()
                );
            }

            //std::cout << "Compressing" << std::endl;
            std::unique_ptr<BinaryMatrix> matrix;

            if (compressor_name == "brwt") {
                matrix = generate_brwt_from_rows(
                    std::move(generated_columns),
                    arity_arg->getValue(),
                    greedy_arg->getValue(),
                    relax_arg->getValue()
                );
                if (greedy_arg->getValue()) {
                    out_name += "_greedy";
                } else {
                    out_name += "_arity_";
                    out_name += std::to_string(arity_arg->getValue());
                }
                if (relax_arg->getValue() > 1) {
                    out_name += "_relax_";
                    out_name += std::to_string(relax_arg->getValue());
                }
            } else {
                matrix = generate_from_rows(
                    std::move(generated_columns),
                    string_to_matrix_type(compressor_name)
                );
            }

            assert(matrix.get());
            assert(matrix->num_columns() == num_cols);

            std::ofstream outfile(out_name + ".annomat", std::ios::binary);
            matrix->serialize(outfile);
            const auto serialized_size = outfile.tellp();

            {
                std::ofstream outfile(out_name + ".stats");
                outfile << compressor_name << " "
                        << matrix->num_rows() << " "
                        << matrix->num_columns() << " "
                        << matrix->num_relations() << " "
                        << serialized_size
                        << std::endl;
            }

            {
                std::ofstream outfile(out_name + ".column_densities");
                for (double density : observed_densities) {
                    outfile << density << " ";
                }
                outfile << std::endl;
            }
        } else if (regime == "stats") {
            std::vector<std::string> compressors {
                "row",
                "flat",
                "column",
                "brwt",
                "bin_rel_wt_sdsl",
                "rbfish"
            };
            ValuesConstraint<std::string> compressor_constraint(compressors);
            ValueArg<std::string> compressor_arg("", "anno-type", "Input compression method type", false, "brwt", &compressor_constraint, cmd);
            UnlabeledMultiArg<std::string> files_arg("input_file", "Input file", true, "string", cmd);
            cmd.parse(argc, argv);

            for (const auto &filename : files_arg.getValue()) {
                std::ifstream in(filename, std::ios::binary);

                auto matrix = generate_from_rows(
                    {},
                    string_to_matrix_type(compressor_arg.getValue())
                );
                if (!matrix->load(in)) {
                    std::cerr << "ERROR: can't load compressor " << filename << std::endl;
                    continue;
                }

                std::cout << "Num rows:\t" << matrix->num_rows() << std::endl;
                std::cout << "Num cols:\t" << matrix->num_columns() << std::endl;
                std::cout << "Num set bits:\t" << matrix->num_relations() << std::endl;
                std::cout << "Density:\t" << static_cast<double>(matrix->num_relations())
                                                / matrix->num_rows() / matrix->num_columns() << std::endl;

                if (dynamic_cast<Rainbowfish*>(matrix.get()))
                    std::cout << "Num distinct rows:\t" << dynamic_cast<Rainbowfish*>(matrix.get())->num_distinct_rows() << std::endl;

                auto *brwt = dynamic_cast<BRWT*>(matrix.get());
                if (!brwt)
                    continue;

                std::cout << "BRWT num nodes:\t" << brwt->num_nodes() << std::endl;
                std::cout << "BRWT avg arity:\t" << brwt->avg_arity() << std::endl;
                std::cout << "BRWT shrink rate:\t" << brwt->shrinking_rate() << std::endl;
                std::cout << "BRWT index col size:\t" << brwt->total_column_size() << std::endl;
                std::cout << "BRWT num set bits:\t" << brwt->total_num_set_bits() << std::endl;
            }
        } else if (regime == "subsets") {
            std::vector<std::string> compressors {
                "row",
                "flat",
                "column",
                "brwt",
                "bin_rel_wt_sdsl",
                "rbfish"
            };
            ValuesConstraint<std::string> compressor_constraint(compressors);
            ValueArg<std::string> compressor_arg("", "anno-type", "Input compression method type", true, "brwt", &compressor_constraint, cmd);
            ValueArg<std::string> out_arg("o", "output", "Name of output file", true, "/dev/null", "std::string", cmd);
            ValueArg<size_t> rows_arg("n", "rows", "Number of rows to subsample", true, 0, "int", cmd);
            ValueArg<size_t> threads_arg("p", "num_threads", "Number of threads", false, 1, "int", cmd);
            UnlabeledMultiArg<std::string> files_arg("input_file", "Input file", true, "string", cmd);
            cmd.parse(argc, argv);
            MatrixType compressor = string_to_matrix_type(compressor_arg.getValue());
            auto files = files_arg.getValue();
            for (const auto &file : files) {
                if (compressor == MatrixType::COLUMN) {
                    annotate::ColumnCompressed<> annotator(0, 1, false);
                    annotator.merge_load({ file });
                    const auto &source_columns = annotator.data();
                    assert(annotator.num_labels() == source_columns.size());

                    auto columns = subsample_rows(
                        source_columns,
                        annotator.num_objects(),
                        rows_arg.getValue(),
                        threads_arg.getValue()
                    );

                    for (auto it = columns.begin(); it != columns.end(); ++it) {
                        if (!(*it)->num_set_bits())
                            columns.erase(it--);
                    }

                    ColMajorCompressed matrix_subsample(convert_to<bit_vector_sd>(std::move(columns)));
                    std::cout << "Reduced matrix from ("
                              << annotator.num_objects() << ", " << annotator.num_labels()
                              << ") to ("
                              << matrix_subsample.num_rows() << ", " << matrix_subsample.num_columns()
                              << ")" << std::endl;

                    std::ofstream out(out_arg.getValue(), std::ios::binary);
                    matrix_subsample.serialize(out);
                } else {
                    throw std::runtime_error("Not supported yet");
                }
            }
        } else if (regime == "to_rrr") {
            std::vector<std::string> compressors {
                "row",
                "flat",
                "column",
                "brwt",
                "bin_rel_wt_sdsl",
                "rbfish"
            };
            ValuesConstraint<std::string> compressor_constraint(compressors);
            ValueArg<std::string> compressor_arg("", "anno-type", "Input compression method type", true, "brwt", &compressor_constraint, cmd);
            UnlabeledMultiArg<std::string> files_arg("input_file", "Input file", true, "string", cmd);
            cmd.parse(argc, argv);
            MatrixType compressor = string_to_matrix_type(compressor_arg.getValue());
            auto files = files_arg.getValue();
            for (const auto &file : files) {
                if (compressor == MatrixType::ROW_FLAT) {
                    annotate::StaticBinRelAnnotator<RowConcatenated<>> annotator;
                    std::cout << "loading\n";
                    annotator.merge_load({ file });
                    std::cout << "done\n";
                    const auto &rows = annotator.data().data();
                    std::ofstream sdout(file + ".sd", std::ios::binary);
                    rows.serialize(sdout);
                    const auto serialized_size = sdout.tellp();
                    std::cout << "Base:\t" << serialized_size << std::endl;
                    std::cout << "converting\n";
                    std::ofstream rrrout(file + ".rrr", std::ios::binary);
                    sdsl::rrr_vector<> rrr(rows.copy_to<sdsl::bit_vector>());
                    std::cout << "Dummy:\t" << rrr.serialize(rrrout) << std::endl;
                } else if (compressor == MatrixType::RAINBOWFISH) {
                    annotate::StaticBinRelAnnotator<Rainbowfish> annotator;
                    std::cout << "loading\n";
                    annotator.merge_load({ file });
                    std::cout << "done\n";
                    const auto &rainbowfish = annotator.data();
                    const auto &row_codes = rainbowfish.get_row_codes();
                    const auto &row_code_delimiters = rainbowfish.get_row_code_delimiters();
                    const auto &distinct_rows = rainbowfish.get_distinct_rows();
                    std::ofstream logout(file + ".sizes.log");
                    logout << annotator.num_objects() << " " << annotator.num_labels() << "\n";
                    uint64_t num_distinct = 0;
                    for (const auto &rowset : distinct_rows) {
                        num_distinct += rowset->num_rows();
                    }
                    logout << distinct_rows.size() << " " << num_distinct << "\n";
                    logout.close();

                    //serialize
                    {
                        std::ofstream outsd(file + ".row_codes.sd", std::ios::binary);
                        row_codes.serialize(outsd);

                        std::ofstream outrrr(file + ".row_codes.rrr", std::ios::binary);
                        row_codes.copy_to<bit_vector_rrr<>>().serialize(outrrr);
                    }
                    {
                        std::ofstream outsd(file + ".row_code_delimiters.sd", std::ios::binary);
                        row_code_delimiters.serialize(outsd);

                        std::ofstream outrrr(file + ".row_code_delimiters.rrr", std::ios::binary);
                        row_code_delimiters.copy_to<bit_vector_rrr<>>().serialize(outrrr);
                    }
                    {
                        std::ofstream outsd(file + ".distinct_rows.sd", std::ios::binary);
                        for (const auto &a : distinct_rows) {
                            a->serialize(outsd);
                        }

                        std::ofstream outrrr(file + ".distinct_rows.rrr", std::ios::binary);
                        for (const auto &a : distinct_rows) {
                            dynamic_cast<const RowConcatenated<> &>(*a).data()
                                .copy_to<bit_vector_rrr<>>().serialize(outrrr);
                        }
                    }
                }
            }
        } else if (regime == "query") {
            std::vector<std::string> query_types {
                "point",
                "row",
                "column"
            };
            ValuesConstraint<std::string> query_constraint(query_types);
            UnlabeledValueArg<std::string> query_arg("query_type", "Query type", true, "", &query_constraint, cmd);
            cmd.parse(std::min(argc, 3), argv);
            std::string query_type = query_arg.getValue();
            cmd.reset();

            std::vector<std::string> compressors {
                "row",
                "flat",
                "column",
                "brwt",
                "brwt_extra",
                "bin_rel_wt_sdsl",
                "bin_rel_wt",
                "rbfish"
            };
            ValuesConstraint<std::string> compressor_constraint(compressors);
            ValueArg<std::string> compressor_arg("", "anno-type", "Input compression method type", true, "brwt", &compressor_constraint, cmd);
            UnlabeledMultiArg<std::string> files_arg("input_file", "Input file", true, "string", cmd);
            ValueArg<size_t> samples_arg("n", "num_samples", "Number of samples to take", false, 100, "int", cmd);
            cmd.parse(argc, argv);

            MatrixType compressor = string_to_matrix_type(compressor_arg.getValue());
            auto files = files_arg.getValue();
            uint64_t num_queries = samples_arg.getValue();
            for (const auto &file : files) {
                std::ofstream query_file(file + "." + query_type + ".query_stats");
                auto matrix = matrix_type_to_data(file, compressor);
                auto file_basename = utils::split_string(file, "/").back();
                auto exp_params = utils::split_string(file_basename, ".")[1];
                auto exp_type = utils::split_string(exp_params, "_").front();
                query_file << compressor_arg.getValue() << " "
                           << exp_type << " "
                           << static_cast<double>(matrix->num_relations())
                                  / (matrix->num_rows() * matrix->num_columns()) << " "
                           << matrix->num_rows() << " "
                           << matrix->num_columns() << " "
                           << query_type << " " << num_queries << " ";
                if (query_type == "point") {
                    query_file << test_point_time(*matrix, num_queries)
                               << std::endl;
                } else if (query_type == "row") {
                    query_file << test_row_time(*matrix, num_queries)
                               << std::endl;
                } else if (query_type == "column") {
                    query_file << test_column_time(*matrix, num_queries)
                               << std::endl;
                }
            }
        } else if (regime == "slice") {
            ValueArg<double> slice_begin_arg("",
                                             "begin",
                                             "Lower bound percentage for slice",
                                             false,
                                             0.0,
                                             "double",
                                             cmd);
            ValueArg<double> slice_end_arg("",
                                           "end",
                                           "Upper bound percentage for slice",
                                           false,
                                           1.0,
                                           "double",
                                           cmd);
            ValueArg<std::string> out_prefix_arg("",
                                                 "outprefix",
                                                 "Prefix out of the output files",
                                                 false,
                                                 "",
                                                 "string",
                                                 cmd);
            UnlabeledMultiArg<std::string> files_arg("input_file",
                                                     "Input file",
                                                     true,
                                                     "string",
                                                     cmd);
            cmd.parse(argc, argv);

            double begin = slice_begin_arg.getValue();
            double end = slice_end_arg.getValue();
            std::string out_prefix = out_prefix_arg.getValue();

            if (begin < 0 || end > 1.0)
                throw std::runtime_error("Begin and end out of bounds");

            auto files = files_arg.getValue();
            annotate::ColumnCompressed<> annotator;
            for (const auto &file : files) {
                std::string outbase = out_prefix.empty()
                    ? file
                    : out_prefix + "/" + utils::split_string(file, "/").back();
                annotator.load(file);
                for (size_t i = 0; i < annotator.data().size(); ++i) {
                    dump_column_slice(
                        *annotator.data()[i],
                        begin,
                        end,
                        outbase + "."
                            + std::to_string(i) + "_"
                            + std::to_string(begin) + "_"
                            + std::to_string(end) + ".dump.txt");
                }
            }
        } else if (regime == "estimate_abundance_threshold") {
            UnlabeledMultiArg<std::string> files_arg("input_file",
                                                     "Input KMC database",
                                                     true,
                                                     "string",
                                                     cmd);
            cmd.parse(argc, argv);

            auto files = files_arg.getValue();

            std::cout << "File\tMin\tMax\tAvg\tCutoff\n" << std::flush;

            for (const auto &file : files) {
                try {
                    uint64_t num_kmers = 0;
                    uint64_t sum_counts = 0;
                    uint64_t min_count = -1;
                    uint64_t max_count = 0;

                    std::vector<uint64_t> hist;

                    kmc::read_kmers(file, [&](std::string&&, uint64_t count) {
                        min_count = std::min(min_count, count);
                        max_count = std::max(max_count, count);
                        sum_counts += count;
                        num_kmers++;

                        assert(count && "All k-mers in graph must have non-zero counts");

                        while (count >= hist.size()) {
                            hist.push_back(0);
                        }
                        hist[count]++;

                    }, true);

                    hist.resize(std::max(uint64_t(hist.size()), uint64_t(10)), 0);

                    double alpha_est_ptr, beta_est_ptr, false_pos_ptr, false_neg_ptr;
                    auto cutoff = cleaning_pick_kmer_threshold(hist.data(), hist.size(),
                                                               &alpha_est_ptr, &beta_est_ptr,
                                                               &false_pos_ptr, &false_neg_ptr);

                    std::string result = file
                            + "\t" + std::to_string(min_count)
                            + "\t" + std::to_string(max_count)
                            + "\t" + std::to_string(sum_counts / num_kmers);

                    if (cutoff != -1) {
                        result += "\t" + std::to_string(cutoff) + "\n";
                    } else {
                        result += "\tnan\n";
                    }

                    std::cout << result << std::flush;

                } catch (...) {
                    std::cerr << "Error: Can't parse file " << file << std::endl;
                }
            }
        } else if (regime == "evaluate_alignment") {
            std::vector<std::string> align_regimes { "score_cigar", };
            ValuesConstraint<std::string> align_regime_constraint(align_regimes);
            UnlabeledValueArg<std::string> align_regime_arg(
                "align_regime", "Align regime", true, "", &align_regime_constraint, cmd
            );

            std::vector<std::string> alphabet_types { "dna", "protein" };
            ValuesConstraint<std::string> alphabet_constraint(alphabet_types);
            UnlabeledValueArg<std::string> alphabet_arg(
                "alphabet", "Alphabet", true, "", &alphabet_constraint, cmd
            );

            std::vector<std::string> modes { "unit", "non-unit" };
            ValuesConstraint<std::string> mode_constraint(modes);
            UnlabeledValueArg<std::string> mode_arg(
                "matrix_type", "Score matrix type", true, "", &mode_constraint, cmd
            );

            UnlabeledValueArg<std::string> file_arg(
                "input_file",
                "Input file (one CIGAR, reference sequence, and alternative sequence per line)",
                true,
                "",
                "string",
                cmd
            );

            cmd.parse(std::min(argc, 6), argv);
            std::string align_regime = align_regime_arg.getValue();
            std::string cur_alphabet = alphabet_arg.getValue();
            const char *alphabet = cur_alphabet == "dna"
                ? alphabets::kAlphabetDNA
                : alphabets::kAlphabetProtein;
            const uint8_t *alphabet_encoding = cur_alphabet == "dna"
                ? alphabets::kCharToDNA
                : alphabets::kCharToProtein;
            std::string mode = mode_arg.getValue();
            std::string file = file_arg.getValue();
            cmd.reset();

            std::unique_ptr<DBGAlignerConfig> config;

            if (mode != "unit" && alphabet_arg.getValue() == "protein") {
                config.reset(new DBGAlignerConfig(
                    DBGAlignerConfig::ScoreMatrix(DBGAlignerConfig::score_matrix_blosum62),
                    -3,
                    -1
                ));
            } else {
                ValueArg<int> match_arg("",
                                        "match",
                                        "score for a single match",
                                        false,
                                        mode == "unit" ? 1 : 2,
                                        "int",
                                        cmd);

                if (mode == "unit") {
                    cmd.parse(argc, argv);
                    config.reset(new DBGAlignerConfig(
                        DBGAlignerConfig::unit_scoring_matrix(match_arg.getValue(),
                                                              alphabet,
                                                              alphabet_encoding),
                        -1,
                        -1
                    ));
                } else {
                    ValueArg<int> gap_open_arg(
                        "", "gap_open_penalty", "gap open penalty", false, 3, "int", cmd
                    );
                    ValueArg<int> gap_ext_arg(
                        "", "gap_ext_penalty", "gap extension penalty", false, 1, "int", cmd
                    );

                    if (cur_alphabet == "dna") {
                        ValueArg<int> mm_transition_arg(
                            "", "transition_penalty", "transition penalty", false, 1, "int", cmd
                        );
                        ValueArg<int> mm_transversion_arg(
                            "", "transversion_penalty", "transversion penalty", false, 2, "int", cmd
                        );
                        cmd.parse(argc, argv);

                        config.reset(new DBGAlignerConfig(
                            DBGAlignerConfig::dna_scoring_matrix(match_arg.getValue(),
                                                                 -mm_transition_arg.getValue(),
                                                                 -mm_transversion_arg.getValue()),
                            -gap_open_arg.getValue(),
                            -gap_ext_arg.getValue()
                        ));
                    } else {
                        throw std::runtime_error("Not implemented for " + cur_alphabet);
                    }
                }
            }

            std::string cigar, strand, ref, alt;
            std::ifstream fin(file);
            while (fin >> strand >> cigar >> ref >> alt) {
                if (strand == "-")
                    reverse_complement(ref.begin(), ref.end());

                Cigar cigar_ob(cigar);
                if (!cigar_ob.is_valid(&*ref.begin(), &*ref.end(),
                                       &*alt.begin(), &*alt.end())) {
                    std::cerr << "ERROR: invalid CIGAR" << std::endl
                              << cigar << std::endl
                              << ref << std::endl
                              << alt << std::endl;

                    exit(1);
                }

                std::cout << config->score_cigar(&*ref.begin(), &*ref.end(),
                                                 &*alt.begin(), &*alt.end(),
                                                 cigar_ob) << std::endl;
            }
        } else if (regime == "query_annotation_rows") {
            UnlabeledValueArg<std::string> annotation_arg("input_file", "Annotation", true, "", "string", cmd);
            ValueArg<int> num_rows_arg("", "num_rows", "Rows to query", true, 0, "int", cmd);
            ValueArg<int> cache_size_arg("", "cache_size", "Number of rows cached", false, 0, "int", cmd);
            ValueArg<int> batch_size_arg("", "batch_size", "Number of rows in batch (task)", false, 100'000, "int", cmd);
            ValueArg<int> num_threads_arg("", "num_threads", "Number threads", false, 1, "int", cmd);
            cmd.parse(argc, argv);

            size_t num_threads = num_threads_arg.getValue() > 0
                                    ? num_threads_arg.getValue()
                                    : 1;
            const size_t batch_size = batch_size_arg.getValue();

            Timer timer;

            auto annotation = initialize_annotation(annotation_arg.getValue(),
                                                    cache_size_arg.getValue());

            std::cout << "Annotation loaded in " << timer.elapsed() << " sec" << std::endl;
            timer.reset();

            std::mt19937 gen;
            gen.seed(17);

            auto row_indexes = utils::sample_indexes(annotation->num_objects(),
                                                     num_rows_arg.getValue(),
                                                     gen);

            std::vector<std::pair<uint64_t, uint64_t>> from_full_to_query;

            for (uint64_t i = 0; i < row_indexes.size(); ++i) {
                from_full_to_query.emplace_back(row_indexes[i], i);
            }

            ips4o::parallel::sort(from_full_to_query.begin(), from_full_to_query.end(),
                [](const auto &first, const auto &second) { return first.first < second.first; },
                num_threads
            );

            std::cout << "Indexes sampled in " << timer.elapsed() << " sec" << std::endl;
            timer.reset();

            // initialize fast query annotation
            // copy annotations from the full graph to the query graph
            auto row_annotation = std::make_unique<annotate::RowCompressed<>>(
                num_rows_arg.getValue(),
                annotation->get_label_encoder().get_labels(),
                [&](annotate::RowCompressed<>::CallRow call_row) {

                    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
                    for (uint64_t batch_begin = 0;
                                        batch_begin < from_full_to_query.size();
                                                        batch_begin += batch_size) {

                        const uint64_t batch_end
                            = std::min(batch_begin + batch_size_arg.getValue(),
                                       static_cast<uint64_t>(from_full_to_query.size()));

                        std::vector<uint64_t> row_indexes;
                        row_indexes.reserve(batch_end - batch_begin);

                        for (uint64_t i = batch_begin; i < batch_end; ++i) {
                            assert(from_full_to_query[i].first < annotation->num_objects());

                            row_indexes.push_back(from_full_to_query[i].first);
                        }

                        auto rows = annotation->get_label_codes(row_indexes);

                        assert(rows.size() == batch_end - batch_begin);

                        for (uint64_t i = batch_begin; i < batch_end; ++i) {
                            call_row(from_full_to_query[i].second,
                                     std::move(rows[i - batch_begin]));
                        }
                    }
                }
            );

            std::cout << "Submatrix constructed in " << timer.elapsed() << " sec" << std::endl;
            timer.reset();

            std::cout << "Num rows in target annotation: " << row_annotation->num_objects() << std::endl;
            std::cout << "Num columns in target annotation: " << row_annotation->num_labels() << std::endl;

        } else if (regime == "to_dna4") {
            // This converts BOSS table constructed for alphabet DNA5 {A,C,G,T,N}
            // to DNA4 {A,C,G,T}. Only tables without any characters N are supported.

            // alphabet: {$,A,C,G,T}, without N
            const size_t alph_size = 5;

            UnlabeledValueArg<std::string> in_graph_arg("input_file", "Graph in", true, "", "string", cmd);
            UnlabeledValueArg<std::string> out_graph_arg("output_file", "Graph out", true, "", "string", cmd);
            cmd.parse(argc, argv);

            std::ifstream instream(in_graph_arg.getValue(), std::ios::binary);
            if (!instream.good()) {
                std::cerr << "ERROR: can't read from file "
                          << in_graph_arg.getValue() << std::endl;
                return 1;
            }

            // load F, k, and state
            auto F_ = libmaus2::util::NumberSerialisation::deserialiseNumberVector<uint64_t>(instream);
            auto k_ = load_number(instream);
            auto state = static_cast<Config::StateType>(load_number(instream));

            // old alphabet: {$,A,C,G,T,N}
            if (F_.size() != alph_size + 1) {
                std::cerr << "ERROR: conversion to DNA4 alphabet is implemented"
                             " only for graphs over the DNA5 alphabet" << std::endl;
                return 1;
            }

            F_.pop_back();

            size_t bits_per_char_W = boost::multiprecision::msb(alph_size - 1) + 2;

            wavelet_tree *W_;
            bit_vector *last_;

            // load W and last arrays
            switch (state) {
                case Config::DYN:
                    W_ = new wavelet_tree_dyn(bits_per_char_W);
                    last_ = new bit_vector_dyn();
                    break;
                case Config::STAT:
                    W_ = new wavelet_tree_stat(bits_per_char_W);
                    last_ = new bit_vector_stat();
                    break;
                case Config::FAST:
                    W_ = new wavelet_tree_fast(bits_per_char_W);
                    last_ = new bit_vector_stat();
                    break;
                case Config::SMALL:
                    W_ = new wavelet_tree_small(bits_per_char_W);
                    last_ = new bit_vector_small();
                    break;
                default:
                    std::cerr << "ERROR: unknown state" << std::endl;
                    return 1;
            }
            if (!W_->load(instream)) {
                std::cerr << "ERROR: failed to load W vector" << std::endl;
                return 1;
            }

            if (!last_->load(instream)) {
                std::cerr << "ERROR: failed to load L vector" << std::endl;
                return 1;
            }

            bool canonical_mode_;

            try {
                canonical_mode_ = load_number(instream);
            } catch (...) {
                canonical_mode_ = false;
            }

            instream.close();

            sdsl::int_vector<> W = W_->to_vector();
            delete W_;

            ProgressBar progress_bar(W.size(), "Converting W[]", std::cerr);
            for (uint64_t i = 0; i < W.size(); ++i) {
                auto value = W[i];

                if (value == alph_size || value == 2 * alph_size + 1) {
                    std::cerr << "ERROR: detected 'N' character in W[" << i << "]" << std::endl;
                    return 1;
                }

                if (value > alph_size)
                    W[i]--;

                ++progress_bar;
            }
            W.width(bits_per_char_W);

            switch (state) {
                case Config::DYN:
                    W_ = new wavelet_tree_dyn(bits_per_char_W, std::move(W));
                    break;
                case Config::STAT:
                    W_ = new wavelet_tree_stat(bits_per_char_W, std::move(W));
                    break;
                case Config::FAST:
                    W_ = new wavelet_tree_fast(bits_per_char_W, std::move(W));
                    break;
                case Config::SMALL:
                    W_ = new wavelet_tree_small(bits_per_char_W, std::move(W));
                    break;
                default:
                    std::cerr << "ERROR: unknown state" << std::endl;
                    return 1;
            }

            std::ofstream outstream(out_graph_arg.getValue(), std::ios::binary);
            if (!outstream.good()) {
                std::cerr << "ERROR: can't write to file "
                          << out_graph_arg.getValue() << std::endl;
                return 1;
            }

            // write F values, k, and state
            libmaus2::util::NumberSerialisation::serialiseNumberVector(outstream, F_);
            serialize_number(outstream, k_);
            serialize_number(outstream, state);
            // write Wavelet Tree
            W_->serialize(outstream);
            // write last array
            last_->serialize(outstream);

            serialize_number(outstream, canonical_mode_);

            std::cout << "Converted graph to DNA4 alphabet and serialized to "
                      << out_graph_arg.getValue() << std::endl;
        }

    } catch (const TCLAP::ArgException &e) {
        std::cerr << "ERROR: " << e.error()
                  << " for arg " << e.argId() << std::endl;
        return 1;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

#include <filesystem>
#include <thread>
#include <chrono>
#include <tclap/CmdLine.h>
#include <sdsl/rrr_vector.hpp>

#include "method_constructors.hpp"
#include "data_generation.hpp"
#include "annotate_static.hpp"
#include "annotate_column_compressed.hpp"
#include "unix_tools.hpp"
#include "utils.hpp"

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
            "slice"
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

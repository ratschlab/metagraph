#include "transform_annotation.hpp"

#include <progress_bar.hpp>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/threads/threading.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/binary_matrix/multi_brwt/clustering.hpp"
#include "annotation/annotation_converters.hpp"
#include "config/config.hpp"
#include "load/load_annotation.hpp"


namespace mtg {
namespace cli {

using namespace mtg::annot;

using mtg::common::logger;
using mtg::common::get_verbose;

typedef MultiLabelEncoded<std::string> Annotator;

static const Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                       Eigen::DontAlignCols, " ", "\n");


template <class AnnotatorTo, class AnnotatorFrom>
void convert(std::unique_ptr<AnnotatorFrom> annotator,
             const Config &config,
             const Timer &timer) {
    logger->trace("Converting annotation to {}...",
                  Config::annotype_to_string(config.anno_type));

    auto target_annotator = convert<AnnotatorTo>(std::move(*annotator));
    annotator.reset();
    logger->trace("Conversion done in {} sec", timer.elapsed());

    logger->trace("Serializing annotation to '{}'...", config.outfbase);
    target_annotator->serialize(config.outfbase);
}

template <class T>
binmat::LinkageMatrix cluster_columns(const std::vector<std::string> &files,
                                      Config::AnnotationType anno_type,
                                      uint64_t num_rows_subsampled) {
    std::vector<uint64_t> row_indexes;
    std::vector<std::unique_ptr<T>> subcolumn_ptrs;
    std::vector<uint64_t> column_ids;
    uint64_t num_rows = 0;

    logger->trace("Loading annotation and sampling subcolumns of size {}",
                  num_rows_subsampled);

    ThreadPool subsampling_pool(get_num_threads(), 1);

    // Load columns from disk
    auto on_column = [&](uint64_t i, const std::string &label,
                         std::unique_ptr<bit_vector> &&column) {
        subsampling_pool.enqueue([&, i, label, column{std::move(column)}]() {
            T *subvector;
            #pragma omp critical
            {
                if (row_indexes.empty()) {
                    num_rows = column->size();
                    if (std::is_same_v<T, sdsl::bit_vector>)
                        row_indexes = binmat::sample_row_indexes(num_rows,
                                                                 num_rows_subsampled);
                } else if (column->size() != num_rows) {
                    logger->error("Size of column {} is {} != {}", label,
                                  column->size(), num_rows);
                    exit(1);
                }
                subcolumn_ptrs.emplace_back(new T());
                subvector = subcolumn_ptrs.back().get();
                column_ids.push_back(i);
                logger->trace("Column {}: {}", i, label);
            }

            if constexpr(std::is_same_v<T, sdsl::bit_vector>) {
                *subvector = sdsl::bit_vector(row_indexes.size(), false);
                for (size_t j = 0; j < row_indexes.size(); ++j) {
                    if ((*column)[row_indexes[j]])
                        (*subvector)[j] = true;
                }
            } else {
                static_assert(std::is_same_v<T, binmat::SparseColumn>);

                uint64_t &size = subvector->size;
                std::vector<uint64_t> &set_bits = subvector->set_bits;

                size = column->num_set_bits() <= num_rows_subsampled
                                    ? column->size()
                                    : column->select1(num_rows_subsampled);

                set_bits.reserve(column->rank1(size));
                column->call_ones_in_range(0, size,
                    [&](uint64_t i) { set_bits.push_back(i); }
                );

                common::logger->trace("Subsampled set bits: {:.2e}/{:.2e}"
                                      ", total size: {:.2e}/{:.2e}, column: {}",
                                      (double)set_bits.size(), (double)column->num_set_bits(),
                                      (double)size, (double)column->size(), label);
            }
        });
    };
    bool success;
    if (anno_type == Config::ColumnCompressed) {
        success = ColumnCompressed<>::merge_load(files, on_column, get_num_threads());
    } else {
        success = merge_load_row_diff(files, on_column, get_num_threads());
    }
    subsampling_pool.join();

    if (!success) {
        logger->error("Could not load annotations");
        exit(1);
    }

    // arrange the columns in their original order
    std::vector<T> subcolumns(subcolumn_ptrs.size());
    for (size_t i = 0; i < column_ids.size(); ++i) {
        subcolumns.at(column_ids[i]) = std::move(*subcolumn_ptrs[i]);
    }

    return binmat::agglomerative_greedy_linkage(std::move(subcolumns), get_num_threads());
}

uint64_t get_num_columns(const std::vector<std::string> &files,
                         Config::AnnotationType anno_type) {
    size_t num_columns = 0;
    std::string extension = anno_type == Config::ColumnCompressed
            ? ColumnCompressed<>::kExtension
            : RowDiffColumnAnnotator::kExtension;
    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
    for (size_t i = 0; i < files.size(); ++i) {
        std::string file = utils::make_suffix(files[i], extension);
        std::ifstream instream(file, std::ios::binary);
        if (!instream.good()) {
            logger->error("Can't read from {}", file);
            exit(1);
        }
        if (anno_type == Config::ColumnCompressed)
            std::ignore = load_number(instream);

        annot::LabelEncoder<std::string> label_encoder;
        if (!label_encoder.load(instream)) {
            logger->error("Can't load label encoder from {}", file);
            exit(1);
        }
        #pragma omp atomic
        num_columns += label_encoder.size();
    }
    return num_columns;
}

binmat::LinkageMatrix trivial_linkage(const std::vector<std::string> &files,
                                      Config::AnnotationType anno_type) {
    logger->trace("Computing total number of columns");
    uint64_t num_columns = get_num_columns(files, anno_type);

    logger->trace("Generating trivial linkage matrix for {} columns",
                  num_columns);

    return binmat::agglomerative_linkage_trivial(num_columns);
}

binmat::LinkageMatrix compute_linkage(const std::vector<std::string> &files,
                                      Config::AnnotationType anno_type,
                                      const Config &config) {
    if (config.greedy_brwt) {
        if (config.fast) {
            return cluster_columns<binmat::SparseColumn>(files, anno_type,
                                                         config.num_rows_subsampled);
        } else {
            return cluster_columns<sdsl::bit_vector>(files, anno_type,
                                                     config.num_rows_subsampled);
        }
    } else {
        return trivial_linkage(files, anno_type);
    }
}

std::vector<std::vector<uint64_t>>
parse_linkage_matrix(const std::string &filename) {
    std::ifstream in(filename);

    std::vector<std::vector<uint64_t>> linkage;
    std::string line;
    while (std::getline(in, line)) {
        std::vector<std::string> parts = utils::split_string(line, " ");
        if (parts.empty())
            continue;

        try {
            if (parts.size() != 4)
                throw std::runtime_error("Invalid format");

            uint64_t first = std::stoi(parts.at(0));
            uint64_t second = std::stoi(parts.at(1));
            uint64_t merged = std::stoi(parts.at(3));

            if (first == second || first >= merged || second >= merged) {
                logger->error("Invalid format of the linkage matrix."
                              " Indexes of parent clusters must be larger than"
                              " indexes of the objects/clusters the include");
                exit(1);
            }

            while (linkage.size() <= merged) {
                linkage.push_back({});
            }

            linkage[merged].push_back(first);
            linkage[merged].push_back(second);

        } catch (const std::exception &e) {
            logger->error("Possibly invalid format of the linkage matrix."
                          " Each line must contain exactly 4 values:"
                          " <cluster 1> <cluster 2> <dist> <cluster 3>"
                          "\nException: {}", e.what());
            exit(1);
        }
    }

    return linkage;
}

auto convert_to_MultiBRWT(const std::vector<std::string> &files,
                          const Config &config) {
    std::string linkage_file = config.linkage_file;
    if (!linkage_file.size()) {
        logger->trace("Generating new column linkage...");
        binmat::LinkageMatrix linkage_matrix
                = compute_linkage(files, Config::ColumnCompressed, config);
        linkage_file = config.outfbase + ".linkage";
        std::ofstream out(linkage_file);
        out << linkage_matrix.format(CSVFormat) << std::endl;
        logger->trace("Generated new linkage and saved to {}", linkage_file);
    }

    auto linkage = parse_linkage_matrix(linkage_file);
    logger->trace("Linkage loaded from {}", linkage_file);

    return convert_to_BRWT<MultiBRWTAnnotator>(
                files, linkage, config.parallel_nodes,
                get_num_threads(), config.tmp_dir);
}

IntMultiBRWTAnnotator
convert_to_IntMultiBRWT(const std::vector<std::string> &files,
                        const Config &config,
                        const Timer &timer) {
    auto brwt_annotator = convert_to_MultiBRWT(files, config);

    logger->trace("Converted to Multi-BRWT in {} sec", timer.elapsed());

    std::vector<CountsVector> column_values;
    ColumnCompressed<>::load_column_values(files,
        [&](size_t j, const std::string &label, sdsl::int_vector<>&& values) {
            if (label != brwt_annotator->get_label_encoder().decode(j)) {
                logger->error("The order of labels in Multi-BRWT is different"
                              " from the order of the input columns");
                exit(1);
            }
            CountsVector values_compressed(std::move(values));
            #pragma omp critical
            {
                while (j >= column_values.size()) {
                    column_values.push_back(sdsl::int_vector<>());
                }
                column_values[j] = std::move(values_compressed);
            }
        },
        get_num_threads()
    );

    if (column_values.size() != brwt_annotator->num_labels()) {
        logger->error("The number of annotation columns does not match"
                      " the number of columns with relation values ({} != {})",
                      brwt_annotator->num_labels(), column_values.size());
        exit(1);
    }

    auto multi_brwt = brwt_annotator->release_matrix();
    return IntMultiBRWTAnnotator(
                std::make_unique<matrix::CSCMatrix<binmat::BRWT, CountsVector>>(
                        std::move(*multi_brwt), std::move(column_values)),
                brwt_annotator->get_label_encoder());
}


int transform_annotation(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    if (config->anno_type == Config::RowDiff && !files.size()) {
        // Only prepare for the row-diff transform:
        //      Compute number of labels per row (if stage 0),
        //      Generate pred/succ/anchors (if stage 1) or optimize anchors (if stage 2).
        logger->trace("Passed no columns to transform. Only preparations will be performed.");
        auto out_dir = std::filesystem::path(config->outfbase).remove_filename();
        convert_to_row_diff({}, config->infbase, config->memory_available * 1e9,
                            config->max_path_length, out_dir, config->tmp_dir,
                            static_cast<RowDiffStage>(config->row_diff_stage));
        logger->trace("Done");
        return 0;
    }

    if (!std::filesystem::exists(files.at(0))) {
        logger->error("File {} does not exist", files.at(0));
        exit(1);
    }

    const Config::AnnotationType input_anno_type
        = parse_annotation_type(files.at(0));

    if (input_anno_type != Config::ColumnCompressed
        && input_anno_type != Config::RowDiff && files.size() > 1) {
        logger->error("Conversion of multiple annotators is only "
                      "supported for ColumnCompressed and ColumnRowDiff");
        exit(1);
    }

    Timer timer;

    /********************************************************/
    /***************** dump labels to text ******************/
    /********************************************************/

    if (config->dump_text_anno) {
        auto annotation = initialize_annotation(files.at(0), *config);

        logger->trace("Loading annotation...");

        if (config->anno_type == Config::ColumnCompressed) {
            if (!dynamic_cast<ColumnCompressed<>&>(*annotation).merge_load(files)) {
                logger->error("Cannot load annotations");
                exit(1);
            }
        } else {
            // Load annotation from disk
            if (!annotation->load(files.at(0))) {
                logger->error("Cannot load annotations from file '{}'", files.at(0));
                exit(1);
            }
        }

        logger->trace("Annotation loaded in {} sec", timer.elapsed());
        logger->trace("Dumping annotators...\t");

        if (input_anno_type == Config::ColumnCompressed) {
            assert(dynamic_cast<ColumnCompressed<>*>(annotation.get()));
            dynamic_cast<ColumnCompressed<>*>(
                annotation.get()
            )->dump_columns(config->outfbase, get_num_threads());
        } else if (input_anno_type == Config::BRWT) {
            assert(dynamic_cast<MultiBRWTAnnotator*>(annotation.get()));
            dynamic_cast<MultiBRWTAnnotator*>(
                annotation.get()
            )->dump_columns(config->outfbase, get_num_threads());
        } else {
            throw std::runtime_error("Dumping columns for this type not implemented");
        }

        logger->trace("Dumping done in {} sec", timer.elapsed());

        return 0;
    }

    /********************************************************/
    /***************** rename column labels *****************/
    /********************************************************/

    if (config->rename_instructions_file.size()) {
        tsl::hopscotch_map<std::string, std::string> dict;
        std::ifstream instream(config->rename_instructions_file);
        if (!instream.is_open()) {
            logger->error("Cannot open file '{}'", config->rename_instructions_file);
            exit(1);
        }
        std::string old_name;
        std::string new_name;
        while (instream.good() && !(instream >> old_name).eof()) {
            instream >> new_name;
            if (instream.fail() || instream.eof()) {
                logger->error("Wrong format of the rules for renaming"
                              " annotation columns passed in file '{}'",
                              config->rename_instructions_file);
                exit(1);
            }
            dict[old_name] = new_name;
        }

        auto annotation = initialize_annotation(files.at(0), *config);

        logger->trace("Loading annotation...");

        // TODO: rename columns without loading the full annotation
        if (config->anno_type == Config::ColumnCompressed) {
            if (!dynamic_cast<ColumnCompressed<>&>(*annotation).merge_load(files)) {
                logger->error("Cannot load annotations");
                exit(1);
            } else {
                logger->info("Annotation #objects: {}\t#labels: {}",
                             annotation->num_objects(), annotation->num_labels());
            }
        } else {
            // Load annotation from disk
            if (!annotation->load(files.at(0))) {
                logger->error("Cannot load annotations from file '{}'", files.at(0));
                exit(1);
            }
        }

        logger->trace("Annotation loaded in {} sec", timer.elapsed());
        logger->trace("Renaming...");

        //TODO: could be made to work with streaming
        annotation->rename_labels(dict);

        annotation->serialize(config->outfbase);
        logger->trace("Renaming done in {} sec", timer.elapsed());

        return 0;
    }

    /********************************************************/
    /**************** operations on columns *****************/
    /********************************************************/

    if (config->aggregate_columns) {
        logger->trace("Loading annotation...");

        uint64_t num_columns = get_num_columns(files, Config::ColumnCompressed);
        if (!num_columns) {
            logger->warn("No input columns to aggregate");
            exit(1);
        }

        const uint64_t min_cols
            = std::max((uint64_t)std::ceil(num_columns * config->min_fraction),
                       (uint64_t)config->min_count);
        const uint64_t max_cols
            = std::min((uint64_t)std::floor(num_columns * config->max_fraction),
                       (uint64_t)config->max_count);

        // TODO: set width to log2(max_cols + 1) but make sure atomic
        //       increments can't lead to overflow
        sdsl::int_vector<> sum(0, 0, sdsl::bits::hi(num_columns) + 1);
        ProgressBar progress_bar(num_columns, "Intersect columns",
                                 std::cerr, !get_verbose());
        ThreadPool thread_pool(get_num_threads(), 1);
        std::mutex mu;
        auto on_column = [&](uint64_t, const auto &, auto&& col) {
            #pragma omp critical
            {
                if (!sum.size()) {
                    sum.resize(col->size());
                } else if (sum.size() != col->size()) {
                    logger->error("Input columns have inconsistent size ({} != {})",
                                  sum.size(), col->size());
                    exit(1);
                }
            }
            thread_pool.enqueue([&,col{std::move(col)}]() {
                col->call_ones([&](uint64_t i) {
                    atomic_fetch_and_add(sum, i, 1, mu, __ATOMIC_RELAXED);
                });
                ++progress_bar;
            });
        };

        if (!ColumnCompressed<>::merge_load(files, on_column, get_num_threads())) {
            logger->error("Couldn't load annotations");
            exit(1);
        }
        thread_pool.join();
        std::atomic_thread_fence(std::memory_order_acquire);

        logger->trace("Selecting k-mers annotated in {} <= * <= {} (out of {}) columns",
                      min_cols, max_cols, num_columns);

        sdsl::bit_vector mask(sum.size(), false);
        for (uint64_t i = 0; i < sum.size(); ++i) {
            if (sum[i] >= min_cols && sum[i] <= max_cols) {
                mask[i] = true;
            }
        }
        sum = sdsl::int_vector<>();

        const std::string &col_name = config->anno_labels.size()
                                        ? config->anno_labels[0]
                                        : "mask";
        ColumnCompressed<> aggregated_column(std::move(mask), col_name);

        const auto &outfname = utils::make_suffix(config->outfbase,
                                                  ColumnCompressed<>::kExtension);
        aggregated_column.serialize(outfname);
        logger->trace("Columns are aggregated and the resulting column '{}'"
                      " serialized to {}", col_name, outfname);
        return 0;
    }

    // TODO: move to column_op
    if (config->intersected_columns.size()) {
        ColumnCompressed<> base_columns(0);
        if (!base_columns.load(config->intersected_columns)) {
            logger->error("Cannot load base columns from file '{}'",
                          config->intersected_columns);
            exit(1);
        }

        logger->trace("Loading input annotations and computing the inner product...");

        for (const auto &file : files) {
            std::unique_ptr<MultiLabelEncoded<std::string>> annotator
                    = initialize_annotation(file, *config);
            if (!annotator->load(file)) {
                logger->error("Cannot load annotations from file '{}'", file);
                exit(1);
            }

            for (const std::string &base_label : base_columns.get_all_labels()) {
                std::vector<std::pair<uint64_t, size_t>> col;
                base_columns.get_column(base_label).call_ones([&col](uint64_t i) {
                    col.emplace_back(i, 1);
                });
                // TODO: make parallel (call sum_rows on batches)
                auto row_sum = annotator->get_matrix().sum_rows(std::move(col),
                                                                config->min_count,
                                                                config->max_count);

                std::cout << fmt::format("({}<{}>, {}<*>):", config->intersected_columns,
                                         base_label, file);
                for (auto [j, sum] : row_sum) {
                    const std::string &label = annotator->get_label_encoder().decode(j);
                    std::cout << fmt::format("\t<{}>:{}", label, sum);
                }
                std::cout << "\n";
            }
        }

        return 0;
    }

    /********************************************************/
    /****************** convert annotation ******************/
    /********************************************************/

    if (config->cluster_linkage) {
        if (input_anno_type != Config::ColumnCompressed
            && input_anno_type != Config::RowDiff) {
            logger->error(
                    "Column clustering is only supported for ColumnCompressed and "
                    "RowDiff");
            exit(1);
        }

        binmat::LinkageMatrix linkage_matrix
                = compute_linkage(files, input_anno_type, *config);

        std::ofstream out(config->outfbase);
        out << linkage_matrix.format(CSVFormat) << std::endl;

        logger->trace("Linkage matrix is written to {}", config->outfbase);
        return 0;
    }

    if (config->anno_type == input_anno_type) {
        logger->info("Skipping conversion: same input and target type: {}",
                      Config::annotype_to_string(config->anno_type));
        return 0;
    }

    logger->trace("Converting to {} annotator...",
                  Config::annotype_to_string(config->anno_type));

    if (input_anno_type == Config::RowCompressed) {

        std::unique_ptr<const Annotator> target_annotator;

        switch (config->anno_type) {
            case Config::RowFlat: {
                auto annotator = annot::convert<RowFlatAnnotator>(files.at(0));
                target_annotator = std::move(annotator);
                break;
            }
            case Config::RowSparse: {
                auto annotator = annot::convert<RowSparseAnnotator>(files.at(0));
                target_annotator = std::move(annotator);
                break;
            }
            case Config::RBFish: {
                auto annotator = annot::convert<RainbowfishAnnotator>(files.at(0));
                target_annotator = std::move(annotator);
                break;
            }
            case Config::BinRelWT_sdsl: {
                auto annotator = annot::convert<BinRelWT_sdslAnnotator>(files.at(0));
                target_annotator = std::move(annotator);
                break;
            }
            case Config::BinRelWT: {
                auto annotator = annot::convert<BinRelWTAnnotator>(files.at(0));
                target_annotator = std::move(annotator);
                break;
            }
            default:
                logger->error(
                        "Streaming conversion from RowCompressed "
                        "annotation is not implemented for the requested "
                        "target type: {}",
                        Config::annotype_to_string(config->anno_type));
                exit(1);
        }

        logger->trace("Annotation converted in {} sec", timer.elapsed());

        logger->trace("Serializing to '{}'...", config->outfbase);

        target_annotator->serialize(config->outfbase);

        logger->trace("Serialization done in {} sec", timer.elapsed());

    } else if (input_anno_type == Config::ColumnCompressed) {
        std::unique_ptr<ColumnCompressed<>> annotator;

        // The entire annotation is loaded in all cases except for transforms
        // to BRWT or RbBRWT, for which the construction is done with streaming
        // columns from disk.
        if (config->anno_type != Config::BRWT
                && config->anno_type != Config::RbBRWT
                && config->anno_type != Config::IntBRWT
                && config->anno_type != Config::RowDiff) {
            annotator = std::make_unique<ColumnCompressed<>>(0);
            logger->trace("Loading annotation from disk...");
            if (!annotator->merge_load(files)) {
                logger->error("Cannot load annotations");
                exit(1);
            }
            logger->trace("Annotation loaded in {} sec", timer.elapsed());
        }

        switch (config->anno_type) {
            case Config::ColumnCompressed: {
                assert(false);
                break;
            }
            case Config::ColumnCoord: {
                auto label_encoder = annotator->get_label_encoder();
                auto tuple_matrix = std::make_unique<matrix::TupleCSCMatrix<binmat::ColumnMajor>>(
                                            annotator->release_matrix());
                if (files.size() > 1) {
                    logger->error("Merging coordinates from multiple columns is not supported");
                    exit(1);
                }
                auto coords_fname = utils::remove_suffix(files.at(0),
                                                         ColumnCompressed<>::kExtension)
                                        + ColumnCompressed<>::kCoordExtension;
                std::ifstream in(coords_fname);
                tuple_matrix->load_tuples(in);

                ColumnCoordAnnotator column_coord(std::move(tuple_matrix), label_encoder);

                logger->trace("Annotation converted in {} sec", timer.elapsed());
                column_coord.serialize(config->outfbase);
                logger->trace("Serialized to {}", config->outfbase);
                break;
            }
            case Config::RowDiffBRWT: {
                logger->error("Convert to row_diff first, and then to row_diff_brwt");
                return 0;

            }
            case Config::RowDiffRowSparse: {
                logger->error("Convert to row_diff first, and then to row_diff_sparse");
                return 0;

            }
            case Config::RowDiff: {
                auto out_dir = std::filesystem::path(config->outfbase).remove_filename();
                convert_to_row_diff(files, config->infbase, config->memory_available * 1e9,
                                    config->max_path_length, out_dir, config->tmp_dir,
                                    static_cast<RowDiffStage>(config->row_diff_stage),
                                    config->outfbase, config->count_kmers,
                                    config->coordinates, config->num_kmers_in_seq);
                break;
            }
            case Config::RowCompressed: {
                if (config->fast) {
                    RowCompressed<> row_annotator(annotator->num_objects());
                    convert_to_row_annotator(*annotator,
                                             &row_annotator,
                                             get_num_threads());
                    annotator.reset();

                    logger->trace("Annotation converted in {} sec", timer.elapsed());
                    logger->trace("Serializing to '{}'...", config->outfbase);

                    row_annotator.serialize(config->outfbase);

                    logger->trace("Serialization done in {} sec", timer.elapsed());

                } else {
                    convert_to_row_annotator(*annotator,
                                             config->outfbase,
                                             get_num_threads());
                    logger->trace("Annotation converted and serialized in {} sec",
                                  timer.elapsed());
                }
                break;
            }
            case Config::BRWT: {
                auto brwt_annotator = convert_to_MultiBRWT(files, *config);
                logger->trace("Annotation converted in {} sec", timer.elapsed());
                brwt_annotator->serialize(config->outfbase);
                logger->trace("Serialized to {}", config->outfbase);
                break;
            }
            case Config::IntBRWT: {
                auto annotator = convert_to_IntMultiBRWT(files, *config, timer);
                logger->trace("Annotation converted in {} sec", timer.elapsed());
                annotator.serialize(config->outfbase);
                logger->trace("Serialized to {}", config->outfbase);
                break;
            }
            case Config::BinRelWT_sdsl: {
                convert<BinRelWT_sdslAnnotator>(std::move(annotator), *config, timer);
                break;
            }
            case Config::BinRelWT: {
                convert<BinRelWTAnnotator>(std::move(annotator), *config, timer);
                break;
            }
            case Config::RowFlat: {
                convert<RowFlatAnnotator>(std::move(annotator), *config, timer);
                break;
            }
            case Config::RowSparse: {
                convert<RowSparseAnnotator>(std::move(annotator), *config, timer);
                break;
            }
            case Config::RBFish: {
                convert<RainbowfishAnnotator>(std::move(annotator), *config, timer);
                break;
            }
            case Config::RbBRWT: {
                auto rb_brwt_annotator
                    = convert_to_RbBRWT<RbBRWTAnnotator>(files, config->relax_arity_brwt);
                logger->trace("Annotation converted in {} sec", timer.elapsed());
                logger->trace("Serializing to '{}'", config->outfbase);
                rb_brwt_annotator->serialize(config->outfbase);
                break;
            }
            case Config::IntRowDiffBRWT: {
                assert(config->infbase.size());
                const std::string anchors_file = config->infbase + annot::binmat::kRowDiffAnchorExt;
                if (!std::filesystem::exists(anchors_file)) {
                    logger->error("Anchor bitmap {} does not exist. Run the row_diff"
                                  " transform followed by anchor optimization.", anchors_file);
                    std::exit(1);
                }
                const std::string fork_succ_file = config->infbase + annot::binmat::kRowDiffForkSuccExt;
                if (!std::filesystem::exists(fork_succ_file)) {
                    logger->error("Fork successor bitmap {} does not exist", fork_succ_file);
                    std::exit(1);
                }
                auto int_annotation = convert_to_IntMultiBRWT(files, *config, timer);
                logger->trace("Annotation converted in {} sec", timer.elapsed());

                using CSCMatrix = matrix::CSCMatrix<binmat::BRWT, CountsVector>;

                IntRowDiffBRWTAnnotator annotation(
                        std::make_unique<matrix::IntRowDiff<CSCMatrix>>(nullptr,
                                std::move(*int_annotation.release_matrix())),
                        int_annotation.get_label_encoder());

                const_cast<matrix::IntRowDiff<CSCMatrix> &>(annotation.get_matrix())
                        .load_anchor(anchors_file);
                const_cast<matrix::IntRowDiff<CSCMatrix> &>(annotation.get_matrix())
                        .load_fork_succ(fork_succ_file);

                logger->trace("RowDiff support bitmaps loaded");

                annotation.serialize(config->outfbase);
                logger->trace("Serialized to {}", config->outfbase);
                break;
            }
        }

    } else if (input_anno_type == Config::RowDiff) {
        if (config->anno_type != Config::RowDiffBRWT
                && config->anno_type != Config::ColumnCompressed
                && config->anno_type != Config::RowDiffRowSparse) {
            logger->error(
                    "Only conversion to 'column', 'row_diff_sparse', and 'row_diff_brwt' "
                    "supported for row_diff");
            exit(1);
        }
        if (config->anno_type == Config::ColumnCompressed) {
            convert_row_diff_to_col_compressed(files, config->outfbase);
        } else {
            assert(config->infbase.size());
            const std::string anchors_file = config->infbase + annot::binmat::kRowDiffAnchorExt;
            if (!std::filesystem::exists(anchors_file)) {
                logger->error("Anchor bitmap {} does not exist. Run the row_diff"
                              " transform followed by anchor optimization.", anchors_file);
                std::exit(1);
            }
            const std::string fork_succ_file = config->infbase + annot::binmat::kRowDiffForkSuccExt;
            if (!std::filesystem::exists(fork_succ_file)) {
                logger->error("Fork successor bitmap {} does not exist", fork_succ_file);
                std::exit(1);
            }
            if (config->anno_type == Config::RowDiffBRWT) {
                if (!config->linkage_file.size()) {
                    logger->trace("Generating new column linkage...");
                    binmat::LinkageMatrix linkage_matrix
                            = compute_linkage(files, input_anno_type, *config);
                    config->linkage_file = config->outfbase + ".linkage";
                    std::ofstream out(config->linkage_file);
                    out << linkage_matrix.format(CSVFormat) << std::endl;
                    logger->trace("Generated new linkage and saved to {}",
                                  config->linkage_file);
                }
                std::vector<std::vector<uint64_t>> linkage
                        = parse_linkage_matrix(config->linkage_file);
                logger->trace("Linkage loaded from {}", config->linkage_file);

                auto brwt_annotator = convert_to_BRWT<RowDiffBRWTAnnotator>(
                        files, linkage, config->parallel_nodes,
                        get_num_threads(), config->tmp_dir);

                logger->trace("Annotation converted in {} sec", timer.elapsed());

                logger->trace("Serializing to '{}'", config->outfbase);
                const_cast<binmat::RowDiff<binmat::BRWT> &>(brwt_annotator->get_matrix())
                        .load_anchor(anchors_file);
                const_cast<binmat::RowDiff<binmat::BRWT> &>(brwt_annotator->get_matrix())
                        .load_fork_succ(fork_succ_file);
                brwt_annotator->serialize(config->outfbase);

            } else { // RowDiff<RowSparse>
                logger->trace("Loading annotation from disk...");
                std::unique_ptr<RowDiffRowSparseAnnotator> row_sparse
                        = convert_row_diff_to_RowDiffSparse(files);
                logger->trace("Annotation converted in {} sec", timer.elapsed());
                const_cast<binmat::RowDiff<binmat::RowSparse> &>(row_sparse->get_matrix())
                        .load_anchor(anchors_file);
                const_cast<binmat::RowDiff<binmat::RowSparse> &>(row_sparse->get_matrix())
                        .load_fork_succ(fork_succ_file);
                logger->trace("Serializing to '{}'", config->outfbase);
                row_sparse->serialize(config->outfbase);
            }
        }
    } else if (config->anno_type == Config::RowDiff) {
        for (const auto &file : files) {
            std::unique_ptr<MultiLabelEncoded<std::string>> annotator
                    = initialize_annotation(file, *config);
            if (!annotator->load(file)) {
                logger->error("Cannot load annotations from file '{}'", file);
                exit(1);
            }

            using std::filesystem::path;
            path out_dir = path(config->outfbase).remove_filename();
            path file_name = path(file).filename().replace_extension("");
            std::string old_extension = file_name.extension();
            file_name = file_name.replace_extension("row_diff_" + old_extension.substr(1));
            std::string out_file = out_dir/(file_name.string() + ".annodbg");

            wrap_in_row_diff(std::move(*annotator), config->linkage_file, out_file);
        }
    } else {
        logger->error(
                "Conversion to other representations is not implemented for {} "
                "annotator", Config::annotype_to_string(input_anno_type));
        exit(1);
    }

    logger->trace("Done");

    return 0;
}

int merge_annotation(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    if (config->anno_type == Config::ColumnCompressed) {
        ColumnCompressed<> annotation(0, config->num_columns_cached);
        if (!annotation.merge_load(files)) {
            logger->error("Cannot load annotations");
            exit(1);
        }
        annotation.serialize(config->outfbase);
        return 0;
    }

    std::vector<std::unique_ptr<Annotator>> annotators;
    std::vector<std::string> stream_files;

    for (const auto &filename : files) {
        auto anno_file_type = parse_annotation_type(filename);
        if (anno_file_type == Config::AnnotationType::RowCompressed) {
            stream_files.push_back(filename);
        } else {
            auto annotator = initialize_annotation(filename, *config);
            if (!annotator->load(filename)) {
                logger->error("Cannot load annotations from file '{}'", filename);
                exit(1);
            }
            annotators.push_back(std::move(annotator));
        }
    }

    if (config->anno_type == Config::RowCompressed) {
        merge<RowCompressed<>>(std::move(annotators), stream_files, config->outfbase);
    } else if (config->anno_type == Config::RowFlat) {
        merge<RowFlatAnnotator>(std::move(annotators), stream_files, config->outfbase);
    } else if (config->anno_type == Config::RBFish) {
        merge<RainbowfishAnnotator>(std::move(annotators), stream_files, config->outfbase);
    } else if (config->anno_type == Config::BinRelWT_sdsl) {
        merge<BinRelWT_sdslAnnotator>(std::move(annotators), stream_files, config->outfbase);
    } else if (config->anno_type == Config::BinRelWT) {
        merge<BinRelWTAnnotator>(std::move(annotators), stream_files, config->outfbase);
    } else if (config->anno_type == Config::BRWT) {
        merge<MultiBRWTAnnotator>(std::move(annotators), stream_files, config->outfbase);
    } else {
        logger->error("Merging of annotations to '{}' representation is not implemented",
                      config->annotype_to_string(config->anno_type));
        exit(1);
    }

    return 0;
}

int relax_multi_brwt(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(files.size() == 1);
    assert(config->outfbase.size());

    Timer timer;

    const std::string &fname = files.at(0);

    std::unique_ptr<MultiLabelEncoded<std::string>> annotator;
    Config::AnnotationType anno_type = parse_annotation_type(fname);
    switch(anno_type) {
        case Config::BRWT:
            annotator = std::make_unique<MultiBRWTAnnotator>();
            break;
        case Config::RowDiffBRWT:
            annotator = std::make_unique<RowDiffBRWTAnnotator>();
            break;
        case Config::IntBRWT:
            annotator = std::make_unique<IntMultiBRWTAnnotator>();
            break;
        case Config::IntRowDiffBRWT:
            annotator = std::make_unique<IntRowDiffBRWTAnnotator>();
            break;
        default:
            logger->error("Relaxation only supported for BRWT and RowDiffBRWT");
            exit(1);
    }

    logger->trace("Loading annotator...");

    if (!annotator->load(fname)) {
        logger->error("Cannot load annotations from file '{}'", files.at(0));
        exit(1);
    }
    logger->trace("Annotator loaded in {} sec", timer.elapsed());

    logger->trace("Relaxing BRWT tree...");

    const binmat::BRWT &matrix = anno_type == Config::BRWT
            ? dynamic_cast<MultiBRWTAnnotator &>(*annotator).get_matrix()
            : (anno_type == Config::IntBRWT
                ? dynamic_cast<IntMultiBRWTAnnotator &>(*annotator).get_matrix().get_binary_matrix()
                : (anno_type == Config::IntRowDiffBRWT
                    ? dynamic_cast<IntRowDiffBRWTAnnotator &>(*annotator).get_matrix().diffs().get_binary_matrix()
                    : dynamic_cast<RowDiffBRWTAnnotator &>(*annotator).get_matrix().diffs()));
    relax_BRWT(const_cast<binmat::BRWT *>(&matrix), config->relax_arity_brwt,
               get_num_threads());

    annotator->serialize(config->outfbase);
    logger->trace("BRWT relaxation done in {} sec", timer.elapsed());

    return 0;
}

} // namespace cli
} // namespace mtg

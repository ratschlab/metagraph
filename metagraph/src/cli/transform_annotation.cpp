#include "transform_annotation.hpp"

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

int transform_annotation(Config *config) {
    assert(config);

    const auto &files = config->fnames;

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
            if (!annotation->merge_load(files)) {
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
            if (!annotation->merge_load(files)) {
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
    /****************** convert annotation ******************/
    /********************************************************/

    if (config->cluster_linkage) {
        if (input_anno_type != Config::ColumnCompressed) {
            logger->error("Column clustering is only supported for ColumnCompressed");
            exit(1);
        }

        logger->trace("Loading annotation and sampling subcolumns of size {}",
                      config->num_rows_subsampled);

        std::vector<uint64_t> row_indexes;
        std::vector<std::unique_ptr<sdsl::bit_vector>> subcolumn_ptrs;
        std::vector<uint64_t> column_ids;
        uint64_t num_rows = 0;
        std::mutex mu;

        ThreadPool subsampling_pool(get_num_threads(), 1);

        // Load columns from disk
        bool success = ColumnCompressed<>::merge_load(files,
            [&](uint64_t i, const std::string &label, auto&& column) {
                subsampling_pool.enqueue([&,i,label,column{std::move(column)}]() {
                    sdsl::bit_vector *subvector;
                    {
                        std::lock_guard<std::mutex> lock(mu);
                        if (row_indexes.empty()) {
                            num_rows = column->size();
                            row_indexes
                                = binmat::sample_row_indexes(num_rows,
                                                             config->num_rows_subsampled);
                        } else if (column->size() != num_rows) {
                            logger->error("Size of column {} is {} != {}",
                                          label, column->size(), num_rows);
                            exit(1);
                        }
                        subcolumn_ptrs.emplace_back(new sdsl::bit_vector());
                        subvector = subcolumn_ptrs.back().get();
                        column_ids.push_back(i);
                        logger->trace("Column {}: {}", i, label);
                    }

                    *subvector = sdsl::bit_vector(row_indexes.size(), false);
                    for (size_t j = 0; j < row_indexes.size(); ++j) {
                        if ((*column)[row_indexes[j]])
                            (*subvector)[j] = true;
                    }
                });
            },
            get_num_threads()
        );

        if (!success) {
            logger->error("Cannot load annotations");
            exit(1);
        }

        subsampling_pool.join();

        // arrange the columns in their original order
        std::vector<std::unique_ptr<sdsl::bit_vector>> permuted(subcolumn_ptrs.size());
        permuted.swap(subcolumn_ptrs);
        for (size_t i = 0; i < column_ids.size(); ++i) {
            subcolumn_ptrs[column_ids[i]] = std::move(permuted[i]);
        }

        std::vector<sdsl::bit_vector> subcolumns;
        for (auto &col_ptr : subcolumn_ptrs) {
            subcolumns.push_back(std::move(*col_ptr));
        }

        binmat::LinkageMatrix linkage_matrix
                = binmat::agglomerative_greedy_linkage(std::move(subcolumns),
                                                       get_num_threads());

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
        std::unique_ptr<annot::MultiLabelEncoded<std::string>> annotation
                = initialize_annotation(files.at(0), *config);

        // The entire annotation is loaded in all cases except for transforms
        // to BRWT with a hierarchical clustering of columns specified (infbase)
        // or RbBRWT, for which the construction is done with streaming columns
        // from disk.
        if ((config->anno_type != Config::BRWT || !config->infbase.size())
            && config->anno_type != Config::RbBRWT
            && config->anno_type != Config::RowDiff) {
            logger->trace("Loading annotation from disk...");
            if (!annotation->merge_load(files)) {
                logger->error("Cannot load annotations");
                exit(1);
            }
            logger->trace("Annotation loaded in {} sec", timer.elapsed());
        }

        std::unique_ptr<ColumnCompressed<>> annotator {
            dynamic_cast<ColumnCompressed<> *>(annotation.release())
        };
        assert(annotator);

        switch (config->anno_type) {
            case Config::ColumnCompressed: {
                assert(false);
                break;
            }
            case Config::BRWTRowDiff: {
                logger->error("Convert to row_diff first, and then to brwt_row_diff");
                return 0;

            }
            case Config::RowDiff: {
                logger->trace("Loading graph...");
                graph::DBGSuccinct graph(2);
                bool result = graph.load(config->infbase);
                if (!result) {
                    logger->error("Cannot load graph from {}", config->infbase);
                    std::exit(1);
                }

                // load as many columns as we can fit in memory, and convert them
                size_t avail_mem_bytes = config->memory_available * 1e9;
                avail_mem_bytes -= std::filesystem::file_size(config->infbase);

                for (uint32_t i = 0; i < files.size();) {
                    logger->trace("Loading columns for batch-conversion...");
                    size_t cur_mem_bytes = avail_mem_bytes;
                    std::vector<std::string> file_batch;
                    for (; i < files.size(); ++i) {
                        // *2 in order to account for constructing the sparsified column
                        size_t file_size = 2 * std::filesystem::file_size(files[i]);
                        if (file_size > avail_mem_bytes) {
                            logger->warn(
                                    "Not enough memory to process {}, requires {} MB",
                                    files[i], file_size/1e6);
                            continue;
                        }
                        if (file_size > cur_mem_bytes) {
                            break;
                        }
                        cur_mem_bytes -= file_size;
                        file_batch.push_back(files[i]);
                    }

                    std::vector<std::unique_ptr<annot::ColumnCompressed<>>> anno_batch;
                    for(const auto& fname : file_batch) {
                        auto anno = std::make_unique<annot::ColumnCompressed<>>() ;
                        anno->merge_load({fname});
                        anno_batch.push_back(std::move(anno));
                    }
                    timer.reset();
                    logger->trace("Starting converting column-batch with {} columns ...",
                                  file_batch.size());
                    std::vector<std::unique_ptr<RowDiffAnnotator>> row_diffs
                            = convert_to_row_diff(graph, anno_batch, config->infbase,
                                                     config->max_path_length);
                    logger->trace("Column-batch converted in {} sec", timer.elapsed());

                    logger->trace("Serializing columns...", config->outfbase);
                    timer.reset();
                    assert(row_diffs.size() == file_batch.size());
                    for(uint32_t idx = 0; idx < file_batch.size(); ++idx) {
                        using std::filesystem::path;
                        auto fname = path(file_batch[idx])
                                             .filename()
                                             .replace_extension()
                                             .replace_extension(
                                                     RowDiffAnnotator::kExtension);
                        auto fpath = path(config->outfbase).remove_filename()/fname;
                        row_diffs[idx]->serialize(fpath);
                        logger->trace("Serialized {}", fpath);
                    }
                    logger->trace("Serialization done in {} sec", timer.elapsed());
                }

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
                auto brwt_annotator = config->infbase.size()
                    ? convert_to_BRWT<MultiBRWTAnnotator>(
                        files, config->infbase,
                        config->parallel_nodes,
                        get_num_threads(),
                        config->tmp_dir.empty()
                            ? std::filesystem::path(config->outfbase).remove_filename()
                            : config->tmp_dir)
                    : (config->greedy_brwt
                        ? convert_to_greedy_BRWT(
                            std::move(*annotator),
                            config->parallel_nodes,
                            get_num_threads(),
                            config->num_rows_subsampled)
                        : convert_to_simple_BRWT(
                            std::move(*annotator),
                            config->arity_brwt,
                            config->parallel_nodes,
                            get_num_threads()));

                annotator.reset();
                logger->trace("Annotation converted in {} sec", timer.elapsed());

                logger->trace("Serializing to '{}'", config->outfbase);

                brwt_annotator->serialize(config->outfbase);

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
        }

    } else if (input_anno_type == Config::RowDiff) {
        if (config->anno_type != Config::BRWTRowDiff
            && config->anno_type != Config::ColumnCompressed) {
            logger->error(
                    "Only conversion to column compressed and brwt supported for "
                    "row_diff");
            exit(1);
        }
        if(config->anno_type == Config::BRWTRowDiff) {
            std::unique_ptr<BRWTRowDiffAnnotator> brwt_annotator;
            if (config->infbase.empty()) { // load all columns in memory and compute linkage on the fly
                logger->trace("Loading annotation from disk...");
                auto annotator = std::make_unique<RowDiffAnnotator>();
                if (!annotator->merge_load(files))
                    std::exit(1);
                logger->trace("Annotation loaded in {} sec", timer.elapsed());
                brwt_annotator = config->greedy_brwt
                        ? convert_to_greedy_BRWT(std::move(*annotator),
                                                 config->parallel_nodes, get_num_threads(),
                                                 config->num_rows_subsampled)
                        : convert_to_simple_BRWT(std::move(*annotator), config->arity_brwt,
                                                 config->parallel_nodes, get_num_threads());
            } else {
                std::string tmp_dir = config->tmp_dir.empty()
                        ? std::filesystem::path(config->outfbase).remove_filename()
                        : config->tmp_dir;
                brwt_annotator
                        = convert_to_BRWT<BRWTRowDiffAnnotator>(files, config->infbase,
                                                                config->parallel_nodes,
                                                                get_num_threads(), tmp_dir);
            }
            logger->trace("Annotation converted in {} sec", timer.elapsed());

            logger->trace("Serializing to '{}'", config->outfbase);
            brwt_annotator->serialize(config->outfbase);
        } else {
            convert_row_diff_to_col_compressed(files, config->outfbase);
        }

    } else {
        logger->error("Conversion to other representations"
                      " is not implemented for {} annotator",
                      Config::annotype_to_string(input_anno_type));
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

    auto annotator = std::make_unique<MultiBRWTAnnotator>();

    logger->trace("Loading annotator...");

    if (!annotator->load(files.at(0))) {
        logger->error("Cannot load annotations from file '{}'", files.at(0));
        exit(1);
    }
    logger->trace("Annotator loaded in {} sec", timer.elapsed());

    logger->trace("Relaxing BRWT tree...");

    relax_BRWT<MultiBRWTAnnotator>(annotator.get(),
                                   config->relax_arity_brwt,
                                   get_num_threads());

    annotator->serialize(config->outfbase);
    logger->trace("BRWT relaxation done in {} sec", timer.elapsed());

    return 0;
}

} // namespace cli
} // namespace mtg

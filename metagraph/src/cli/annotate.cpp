#include "annotate.hpp"

#include <algorithm>
#include <filesystem>
#include <optional>
#include <string_view>

#include <omp.h>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/batch_accumulator.hpp"
#include "common/threads/threading.hpp"
#include "common/utils/string_utils.hpp"
#include "common/utils/file_utils.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "seq_io/formats.hpp"
#include "seq_io/sequence_io.hpp"
#include "seq_io/kmc_parser.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/coord_to_header.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"


namespace mtg {
namespace cli {

namespace fs = std::filesystem;

using namespace mtg::seq_io;

using mtg::common::logger;

const size_t LOGGING_INTERVAL_SEQS = 100'000;


void push_header_labels(const kseq_t *read_stream,
                        const std::string &file,
                        size_t record,
                        const std::string &comment_delim,
                        const std::string &header_delim,
                        bool unchecked,
                        std::vector<std::string> *labels) {
    // kseq resets comment.l for every record but writes comment.s only for a header that has
    // a comment, so the buffer still spells the last comment read: the length tells whether
    // this record has one
    const std::string header = comment_delim != Config::UNINITIALIZED_STR
            ? utils::join_strings({ read_stream->name.s,
                                    read_stream->comment.l ? read_stream->comment.s : "" },
                                  comment_delim, true)
            : read_stream->name.s;
    for (auto &label : utils::split_string(header, header_delim)) {
        if (!unchecked && !utils::valid_utf8(label)) {
            logger->error("Header of record {} in '{}' is not valid UTF-8 (label bytes in hex: {}):"
                          " the label names of an index must be UTF-8 for the serving routes to"
                          " list them (--allow-unchecked-headers skips this check)",
                          record, file, utils::to_hex(label));
            exit(1);
        }
        labels->push_back(std::move(label));
    }
}

// |callback|(sequence, labels, abundance, record) for every sequence of |file|, the record its
// 1-based number in the file (with |forward_and_reverse| a record's reverse complement is
// passed under the same number)
template <class Callback>
void call_annotations(const std::string &file,
                      const std::string &ref_sequence_path,
                      const graph::DeBruijnGraph &graph,
                      bool forward_and_reverse,
                      size_t min_count,
                      size_t max_count,
                      bool filename_anno,
                      bool annotate_sequence_headers,
                      bool parse_counts_from_headers,
                      const std::string &fasta_anno_comment_delim,
                      const std::string &fasta_header_delimiter,
                      bool unchecked_headers,
                      const std::vector<std::string> &anno_labels,
                      const Callback &callback) {
    size_t total_seqs = 0;

    Timer timer;

    logger->trace("Parsing '{}'", file);

    std::vector<std::string> labels = anno_labels;
    if (filename_anno) {
        labels.push_back(file);
    }
    // remember the number of base labels to remove those unique to each sequence quickly
    const size_t num_base_labels = labels.size();

    if (file_format(file) == "VCF") {
        read_vcf_file_with_annotations_critical(file, ref_sequence_path, graph.get_k(),
            [&](auto&& seq, const auto &variant_labels) {
                labels.insert(labels.end(),
                              variant_labels.begin(), variant_labels.end());

                callback(std::move(seq), labels, 1, total_seqs + 1);

                total_seqs += 1;

                if (logger->level() <= spdlog::level::level_enum::trace
                                                && total_seqs % LOGGING_INTERVAL_SEQS == 0) {
                    logger->trace(
                        "processed {} variants, last was annotated as <{}>, {} sec",
                        total_seqs, fmt::join(labels, "><"), timer.elapsed());
                }

                labels.resize(num_base_labels);
            },
            forward_and_reverse
        );
    } else if (file_format(file) == "KMC") {
        read_kmers(
            file,
            [&](std::string_view sequence) {
                callback(std::string(sequence), labels, 1, total_seqs + 1);

                total_seqs += 1;

                if (logger->level() <= spdlog::level::level_enum::trace
                                                && total_seqs % LOGGING_INTERVAL_SEQS == 0) {
                    logger->trace("processed {} sequences, annotated as <{}>, {} sec",
                                  total_seqs, fmt::join(labels, "><"), timer.elapsed());
                }
            },
            graph.get_mode() != graph::DeBruijnGraph::CANONICAL, // |call_both_from_canonical|
            min_count,
            max_count
        );
    } else if (file_format(file) == "FASTA"
                || file_format(file) == "FASTQ") {
        bool parsed_first_header = false;
        size_t calls = 0;
        read_fasta_file_critical(
            file,
            [&](kseq_t *read_stream) {
                // two calls per record with |forward_and_reverse|, the second its reverse
                // complement: the same record, its header checked on the first call
                const size_t record = forward_and_reverse ? calls / 2 + 1 : calls + 1;
                const bool reverse = forward_and_reverse && calls % 2 == 1;
                ++calls;
                // add sequence header to labels
                if (annotate_sequence_headers) {
                    push_header_labels(read_stream, file, record,
                                       fasta_anno_comment_delim, fasta_header_delimiter,
                                       unchecked_headers || reverse, &labels);
                }

                uint64_t abundance = 1;

                if (parse_counts_from_headers) {
                    // the comment buffer outlives its record (see push_header_labels): only a
                    // record with a comment carries a count
                    if (auto count = read_stream->comment.l ? utils::parse_abundance(read_stream->comment.s) : std::nullopt) {
                        abundance = *count;
                        parsed_first_header = true;
                    } else if (parsed_first_header) {
                        mtg::common::logger->error("Inconsistent k-mer count headers in file '{}'", file);
                        exit(1);
                    } else {
                        parse_counts_from_headers = false;
                        mtg::common::logger->warn("No k-mer count found in header '{}', "
                                                  "will treat all sequences as having k-mer count 1",
                                                  read_stream->name.s);
                    }
                }

                callback(read_stream->seq.s, labels, abundance, record);

                total_seqs += 1;

                if (logger->level() <= spdlog::level::level_enum::trace
                                                && total_seqs % LOGGING_INTERVAL_SEQS == 0) {
                    logger->trace("processed {} sequences, last was {}, annotated as <{}>, {} sec",
                                  total_seqs, read_stream->name.s, fmt::join(labels, "><"), timer.elapsed());
                }

                labels.resize(num_base_labels);
            },
            forward_and_reverse
        );
    } else {
        logger->error("Unknown input format for '{}'. Use FASTA, FASTQ, VCF or KMC input.", file);
        exit(1);
    }

    logger->trace("File '{}' processed in {} sec, current mem usage: {} MiB",
                  file, timer.elapsed(), get_curr_RSS() >> 20);
}

template <class Callback>
void add_kmer_counts(const std::string &file,
                     const graph::DeBruijnGraph &graph,
                     bool forward_and_reverse,
                     bool filename_anno,
                     bool annotate_sequence_headers,
                     const std::string &fasta_anno_comment_delim,
                     const std::string &fasta_header_delimiter,
                     bool unchecked_headers,
                     const std::vector<std::string> &anno_labels,
                     const Callback &callback) {
    size_t total_seqs = 0;

    Timer timer;

    std::vector<std::string> labels = anno_labels;
    if (filename_anno) {
        labels.push_back(file);
    }
    // remember the number of base labels to remove those unique to each sequence quickly
    const size_t num_base_labels = labels.size();

    const std::string &counts_fname
            = utils::remove_suffix(file, ".gz", ".fasta") + ".kmer_counts.gz";

    logger->trace("Parsing k-mer counts from '{}'", counts_fname);

    size_t calls = 0;
    read_extended_fasta_file_critical<uint32_t>(file, "kmer_counts",
        [&](size_t k, const kseq_t *read_stream, const uint32_t *kmer_counts) {
            if (k != graph.get_k()) {
                logger->error("File '{}' contains counts for k-mers of "
                              "length {} but graph is constructed with k={}",
                              file, k, graph.get_k());
                exit(1);
            }
            assert(read_stream->seq.l >= k && "sequences can't be shorter than k-mers");

            // as in call_annotations: two calls per record with |forward_and_reverse|
            const size_t record = forward_and_reverse ? calls / 2 + 1 : calls + 1;
            const bool reverse = forward_and_reverse && calls % 2 == 1;
            ++calls;
            // add sequence header to labels
            if (annotate_sequence_headers) {
                push_header_labels(read_stream, file, record,
                                   fasta_anno_comment_delim, fasta_header_delimiter,
                                   unchecked_headers || reverse, &labels);
            }

            callback(read_stream->seq.s, labels,
                     std::vector<uint64_t>(kmer_counts, kmer_counts + read_stream->seq.l - k + 1));

            total_seqs++;

            if (logger->level() <= spdlog::level::level_enum::trace
                                            && total_seqs % LOGGING_INTERVAL_SEQS == 0) {
                logger->trace("processed {} sequences, last was {}, annotated as <{}>, {} sec",
                              total_seqs, read_stream->name.s, fmt::join(labels, "><"), timer.elapsed());
            }

            labels.resize(num_base_labels);
        },
        forward_and_reverse
    );
}

// Runs the header check of push_header_labels over every FASTA/FASTQ file of |files| (other
// formats carry no headers) before any annotation is written or the output directory is
// created: with --separately each file's annotation is serialized as soon as it is annotated,
// so a label refused during a later file's annotation would leave the earlier files'
// annotations on disk, while a refusal writes nothing
void check_header_labels(const std::vector<std::string> &files,
                         const Config &config,
                         size_t num_threads) {
    logger->trace("Checking the headers of {} files", files.size());
    #pragma omp parallel for num_threads(num_threads) default(shared) schedule(dynamic, 1)
    for (size_t i = 0; i < files.size(); ++i) {
        if (file_format(files[i]) != "FASTA" && file_format(files[i]) != "FASTQ") {
            continue;
        }
        size_t record = 0;
        std::vector<std::string> labels;
        read_fasta_file_critical(files[i],
            [&](kseq_t *read_stream) {
                push_header_labels(read_stream, files[i], ++record,
                                   config.fasta_anno_comment_delim,
                                   config.fasta_header_delimiter,
                                   false, &labels);
                labels.clear();
            },
            false
        );
    }
}

// A stretch of a file's sequences, for the spot check of a record mapping against an annotation
// (annotate --index-header-coords -a): |sequence| (at most kProbeLength characters, holding a
// k-mer the graph has; empty: none found) and the coordinate of its first k-mer in the file's
// column
struct CoordinateProbe {
    std::string sequence;
    uint64_t first = 0;
};
constexpr size_t kProbeLength = 1000;
// the end of each record kept for its file's tail probe, looked for in it from its end: an
// assembly's gap of Ns (10 kbp at a GRCh38 chromosome's end) has no k-mer of the graph
constexpr size_t kTailSearch = 1 << 16;

// The probe of |sequence|, whose first k-mer is at |coordinate|: the first stretch of
// kProbeLength characters holding a k-mer the graph has, from its start or (|from_end|) from its
// end, the stretches overlapping by k - 1 characters so that each k-mer is in one; none when the
// graph has none of its k-mers (characters it cannot hold, a gap of Ns, are passed over so)
std::optional<CoordinateProbe> find_probe(const graph::DeBruijnGraph &dbg,
                                          std::string_view sequence, uint64_t coordinate,
                                          bool from_end) {
    const size_t k = dbg.get_k();
    if (sequence.size() < k)
        return std::nullopt;
    const size_t num_kmers = sequence.size() - k + 1;
    const size_t step = std::max(kProbeLength, k) - k + 1;
    for (size_t done = 0; done < num_kmers; done += step) {
        // the k-mers [begin, end) of this stretch
        const size_t begin = from_end ? num_kmers - std::min(num_kmers, done + step) : done;
        const size_t end = from_end ? num_kmers - done : std::min(num_kmers, done + step);
        const std::string_view stretch = sequence.substr(begin, end - begin + k - 1);
        bool found = false;
        dbg.map_to_nodes(stretch, [&](auto node) { found |= node != graph::DeBruijnGraph::npos; });
        if (found)
            return CoordinateProbe { std::string(stretch), coordinate + begin };
    }
    return std::nullopt;
}

// Whether column |column| of |anno_graph|'s coordinate annotation (|matrix|) holds the coordinate
// of the first and of the last k-mer of |probe| that the graph has, at their nodes (an empty
// probe: nothing to check): -1 when it does, else the coordinate it lacks
int64_t check_probe(const graph::AnnotatedDBG &anno_graph,
                    const annot::matrix::MultiIntMatrix &matrix, uint64_t column,
                    const CoordinateProbe &probe) {
    using Row = graph::AnnotatedDBG::row_index;
    std::vector<std::pair<uint64_t, Row>> found;
    uint64_t position = 0;
    // as the coordinates were annotated and as queries map (canonical nodes on a canonical graph)
    anno_graph.get_graph().map_to_nodes(probe.sequence, [&](auto node) {
        if (node != graph::DeBruijnGraph::npos)
            found.emplace_back(probe.first + position, graph::AnnotatedDBG::graph_to_anno_index(node));
        ++position;
    });
    if (found.size() > 2)
        found.erase(found.begin() + 1, found.end() - 1);
    for (const auto &[coordinate, row] : found) {
        const auto tuples = matrix.get_row_tuples(std::vector<Row> { row });
        const bool held = std::any_of(tuples.at(0).begin(), tuples.at(0).end(), [&](const auto &t) {
            return t.first == column && std::find(t.second.begin(), t.second.end(), coordinate)
                                            != t.second.end();
        });
        if (!held)
            return static_cast<int64_t>(coordinate);
    }
    return -1;
}

void annotate_data(std::shared_ptr<graph::DeBruijnGraph> graph,
                   const Config &config,
                   const std::vector<std::string> &files,
                   const std::string &annotator_filename,
                   size_t num_chunks_open = 2000) {
    auto anno_graph = initialize_annotated_dbg(graph, config, num_chunks_open);
    const size_t k = anno_graph->get_graph().get_k();

    bool forward_and_reverse = config.forward_and_reverse;
    if (anno_graph->get_graph().get_mode() == graph::DeBruijnGraph::CANONICAL) {
        logger->trace("Annotating canonical graph");
        forward_and_reverse = false;
    }

    // not too small, not too large
    const size_t BATCH_SIZE = 5'000;
    const size_t BATCH_LENGTH = 100'000;

    if (config.index_header_coords) {
        // Collect sequence headers and compute k-mer counts
        std::vector<std::vector<std::string>> headers(files.size());
        std::vector<std::vector<uint64_t>> num_kmers(files.size());
        // with -a, each file's first stretch with a k-mer the graph has (find_probe), and the
        // last such stretch in the end of its last record, whose coordinates the annotation's
        // column of the file must hold
        const bool probe = !config.infbase_annotators.empty();
        const graph::DeBruijnGraph &probed_graph = anno_graph->get_graph();
        std::vector<CoordinateProbe> heads(files.size());
        std::vector<CoordinateProbe> tails(files.size());
        #pragma omp parallel for num_threads(get_num_threads()) default(shared) schedule(dynamic)
        for (size_t i = 0; i < files.size(); ++i) {
            size_t last_record = 0;
            uint64_t coordinate = 0;
            // the end of the last record so far (at most kTailSearch characters), from
            // |tail_coordinate| on
            std::string tail;
            uint64_t tail_coordinate = 0;
            // the header of every record, those too short to have k-mers too: a record
            // without k-mers is not in the record mapping, but a repeat of its header is
            // refused all the same, so that a file is accepted or refused by its headers alone
            std::vector<std::string> record_headers;
            call_annotations(
                files[i],
                config.refpath,
                anno_graph->get_graph(),
                forward_and_reverse,
                config.min_count,
                config.max_count,
                false, // config.filename_anno
                true, // config.annotate_sequence_headers
                false, // parse_counts_from_headers
                config.fasta_anno_comment_delim,
                config.fasta_header_delimiter,
                config.allow_unchecked_headers,
                {}, // config.anno_labels
                [&](std::string sequence, auto labels, uint64_t, size_t record) {
                    // a header of `>` alone yields no label (nothing to split)
                    if (labels.empty()) {
                        logger->error("Header of record {} in '{}' is empty: a record of the"
                                      " record mapping (.seqs) needs a header name",
                                      record, files[i]);
                        exit(1);
                    }
                    // the header names its record whole (annotate_graph refuses
                    // --header-delimiter beside --index-header-coords)
                    if (labels.size() != 1) {
                        logger->error("Header of record {} in '{}' is split into {} names: a"
                                      " record of the record mapping (.seqs) is named by its"
                                      " whole header", record, files[i], labels.size());
                        exit(1);
                    }
                    // (a record's reverse complement is passed under its number again)
                    if (!config.allow_unchecked_headers && record > record_headers.size())
                        record_headers.push_back(labels[0]);
                    if (sequence.size() < k)
                        return;
                    // a record's reverse complement (--fwd-and-reverse) continues the record's
                    // coordinates, as `annotate --coordinates` assigns them; with unchecked
                    // headers, consecutive records under one header are merged into one as well
                    const bool same_record = !headers[i].empty()
                        && (record == last_record
                            || (config.allow_unchecked_headers && headers[i].back() == labels[0]));
                    if (!same_record) {
                        headers[i].push_back(labels[0]);
                        num_kmers[i].push_back(0);
                    }
                    if (probe) {
                        // looked for in each record until found: a file with no k-mer of the
                        // graph is mapped whole, and refused below
                        if (heads[i].sequence.empty()) {
                            if (auto head = find_probe(probed_graph, sequence, coordinate, false))
                                heads[i] = std::move(*head);
                        }
                        const size_t kept = std::min(sequence.size(), kTailSearch);
                        tail.assign(sequence, sequence.size() - kept, kept);
                        tail_coordinate = coordinate + (sequence.size() - kept);
                    }
                    coordinate += sequence.size() - k + 1;
                    num_kmers[i].back() += sequence.size() - k + 1;
                    last_record = record;
                }
            );
            if (probe) {
                if (auto last = find_probe(probed_graph, tail, tail_coordinate, true))
                    tails[i] = std::move(*last);
            }
            if (!config.allow_unchecked_headers) {
                // find_header() answers a header with its first record only
                assert(record_headers.size() >= headers[i].size());
                annot::CoordToHeader::find_duplicate_headers(record_headers,
                    [&](size_t first, size_t repeat) {
                        logger->error("Header of record {} in '{}' repeats the header of record {}"
                                      " (bytes in hex: {}): the record headers of a file must be"
                                      " distinct for the record mapping (.seqs), which finds a"
                                      " header that repeats record {}'s at record {} only"
                                      " (--allow-unchecked-headers skips this check)",
                                      repeat + 1, files[i], first + 1,
                                      utils::to_hex(record_headers[repeat]), first + 1, first + 1);
                        exit(1);
                        return false;
                    }
                );
            }
        }
        // each column's label, stated in the .seqs: a load compares them with the annotation's
        // labels, so that files passed in another order are refused there. Without -a, the label
        // --anno-filename gives the file the column was made from (its name as passed); with -a,
        // the annotation's own, once each file passed the spot check against its column: the
        // column holds the coordinates of the file's probes (a renamed column, or files at
        // another path, then pass the load's check). A spot check: a file of the same sequences
        // under other headers, or with other records between its probes, passes it (what the
        // records are rests on the files passed being those the annotation was made from)
        std::vector<std::string> column_labels(files.begin(), files.end());
        if (probe) {
            const std::string &anno_file = config.infbase_annotators.at(0);
            const auto &annotation = anno_graph->get_annotator();
            const auto *matrix
                    = dynamic_cast<const annot::matrix::MultiIntMatrix *>(&annotation.get_matrix());
            if (!matrix) {
                logger->error("--index-header-coords with -a needs a coordinate annotation: '{}' "
                              "has no k-mer coordinates", anno_file);
                exit(1);
            }
            const auto &labels = annotation.get_label_encoder().get_labels();
            if (labels.size() != files.size()) {
                logger->error("The annotation '{}' has {} columns, and {} files were passed: pass "
                              "every input file of the annotation in its column order (metagraph "
                              "stats --print-col-names lists the columns)", anno_file,
                              labels.size(), files.size());
                exit(1);
            }
            for (size_t i = 0; i < files.size(); ++i) {
                // a file with k-mers none of which the graph has (before GPT review 6 it passed,
                // nothing checked); one without k-mers has no records in the mapping to bind
                if (heads[i].sequence.empty() && !headers[i].empty()) {
                    logger->error("The file '{}' has no k-mer of the graph: it cannot be checked "
                                  "against column {} ('{}') of the annotation '{}'. Pass the files "
                                  "the annotation was made from, in its column order (metagraph "
                                  "stats --print-col-names lists the columns)", files[i], i,
                                  labels[i], anno_file);
                    exit(1);
                }
                if (!heads[i].sequence.empty() && tails[i].sequence.empty()) {
                    logger->warn("The last record of '{}' has no k-mer of the graph in its last {} "
                                 "characters: the file is checked by its first k-mers alone",
                                 files[i], kTailSearch);
                }
                for (const CoordinateProbe *p : { &heads[i], &tails[i] }) {
                    if (const int64_t lacking = check_probe(*anno_graph, *matrix, i, *p);
                            lacking >= 0) {
                        logger->error("The file '{}' is not column {} ('{}') of the annotation "
                                      "'{}': the column lacks the file's k-mer at coordinate {}. "
                                      "Pass the files in the annotation's column order (metagraph "
                                      "stats --print-col-names lists the columns)", files[i], i,
                                      labels[i], anno_file, lacking);
                        exit(1);
                    }
                }
            }
            logger->trace("Each of the {} files passed the spot check against its column of the "
                          "annotation '{}': its labels are stated in the record mapping",
                          files.size(), anno_file);
            column_labels = labels;
        }
        annot::CoordToHeader headers_map(std::move(headers), std::move(num_kmers));
        headers_map.set_column_labels(std::move(column_labels));
        auto fname = utils::make_suffix(annotator_filename, annot::CoordToHeader::kExtension);
        headers_map.serialize(fname);
        logger->trace("CoordToHeader mapping serialized to {}", fname);
        return;
    }

    if (config.coordinates) {
        #pragma omp parallel num_threads(get_num_threads())
        #pragma omp single
        for (const auto &file : files) {
            using Seq = std::tuple<std::string, std::vector<std::string>, uint64_t>;
            BatchAccumulator<Seq> batcher(
                [graph_p = anno_graph.get()](std::vector<Seq>&& data) {
                    auto *data_p = new std::vector<Seq>(std::move(data));
                    // The task can outlive the batcher's callback at a file boundary.
                    #pragma omp task firstprivate(data_p, graph_p)
                    {
                        std::unique_ptr<std::vector<Seq>> data_ptr(data_p);
                        graph_p->annotate_kmer_coords(*data_ptr);
                    }
                },
                BATCH_SIZE, BATCH_LENGTH, BATCH_SIZE
            );

            logger->trace("Annotating k-mer coordinates for file {}", file);

            if (file_format(file) != "FASTA"
                    && file_format(file) != "FASTQ") {
                logger->error("--coordinates requires FASTA or FASTQ input: '{}'", file);
                exit(1);
            }

            uint64_t coord = 0;
            call_annotations(
                file,
                config.refpath,
                anno_graph->get_graph(),
                forward_and_reverse,
                config.min_count,
                config.max_count,
                config.filename_anno,
                config.annotate_sequence_headers,
                /*parse_counts_from_headers*/false,
                config.fasta_anno_comment_delim,
                config.fasta_header_delimiter,
                config.allow_unchecked_headers,
                config.anno_labels,
                [&](std::string sequence, auto labels, uint64_t, size_t) {
                    if (config.num_kmers_in_seq
                            && config.num_kmers_in_seq + k - 1 != sequence.size()) {
                        logger->error("--num-kmers-in-seq requires every sequence to contain"
                                      " that many k-mers (sequence length = value + k - 1)");
                        exit(1);
                    }
                    if (sequence.size() >= k) {
                        uint64_t num_kmers = sequence.size() - k + 1;
                        batcher.push_and_pay(sequence.size(),
                                             std::move(sequence), std::move(labels), coord);
                        if (!config.annotate_sequence_headers)
                            coord += num_kmers;
                    }
                }
            );
        }

        anno_graph->get_annotator().serialize(annotator_filename);
        return;
    }

    #pragma omp parallel num_threads(get_num_threads())
    #pragma omp single
    for (const auto &file : files) {
        using Seq = std::pair<std::string, std::vector<std::string>>;
        BatchAccumulator<Seq> batcher(
            [graph_p = anno_graph.get()](std::vector<Seq>&& data) {
                auto *data_p = new std::vector<Seq>(std::move(data));
                // The task can outlive the batcher's callback at a file boundary.
                #pragma omp task firstprivate(data_p, graph_p)
                {
                    std::unique_ptr<std::vector<Seq>> data_ptr(data_p);
                    graph_p->annotate_sequences(*data_ptr);
                }
            },
            BATCH_SIZE, BATCH_LENGTH, BATCH_SIZE
        );
        call_annotations(
            file,
            config.refpath,
            anno_graph->get_graph(),
            forward_and_reverse,
            config.min_count,
            config.max_count,
            config.filename_anno,
            config.annotate_sequence_headers,
            /*parse_counts_from_headers*/false,
            config.fasta_anno_comment_delim,
            config.fasta_header_delimiter,
            config.allow_unchecked_headers,
            config.anno_labels,
            [&](std::string sequence, auto labels, uint64_t, size_t) {
                if (sequence.size() >= k) {
                    batcher.push_and_pay(sequence.size(),
                                         std::move(sequence), std::move(labels));
                }
            }
        );
    }

    if (config.count_kmers) {
        // add k-mer counts to existing binary annotations
        #pragma omp parallel num_threads(get_num_threads())
        #pragma omp single
        for (const auto &file : files) {
            logger->trace("Annotating k-mer counts for file {}", file);

            using Seq = std::tuple<std::string, std::vector<std::string>, std::vector<uint64_t>>;
            BatchAccumulator<Seq> batcher(
                [graph_p = anno_graph.get()](std::vector<Seq>&& data) {
                    auto *data_p = new std::vector<Seq>(std::move(data));
                    // The task can outlive the batcher's callback at a file boundary.
                    #pragma omp task firstprivate(data_p, graph_p)
                    {
                        std::unique_ptr<std::vector<Seq>> data_ptr(data_p);
                        for (auto &[seq, labels, kmer_counts] : *data_ptr) {
                            graph_p->add_kmer_counts(seq, labels, std::move(kmer_counts));
                        }
                    }
                },
                BATCH_SIZE, BATCH_LENGTH, BATCH_SIZE
            );

            const std::string &counts_fname
                    = utils::remove_suffix(file, ".gz", ".fasta") + ".kmer_counts.gz";

            if (fs::exists(counts_fname)) {
                add_kmer_counts(
                    file,
                    anno_graph->get_graph(),
                    forward_and_reverse,
                    config.filename_anno,
                    config.annotate_sequence_headers,
                    config.fasta_anno_comment_delim,
                    config.fasta_header_delimiter,
                    config.allow_unchecked_headers,
                    config.anno_labels,
                    [&](std::string sequence,
                                std::vector<std::string> labels,
                                std::vector<uint64_t> kmer_counts) {
                        if (sequence.size() >= k) {
                            batcher.push_and_pay(sequence.size(),
                                                 std::move(sequence), std::move(labels),
                                                 std::move(kmer_counts));
                        }
                    }
                );
            } else {
                logger->warn("No k-mer counts found at '{}', "
                             "will try reading counts from headers",
                             counts_fname);
                call_annotations(
                    file,
                    config.refpath,
                    anno_graph->get_graph(),
                    forward_and_reverse,
                    config.min_count,
                    config.max_count,
                    config.filename_anno,
                    config.annotate_sequence_headers,
                    /*parse_counts_from_headers*/true,
                    config.fasta_anno_comment_delim,
                    config.fasta_header_delimiter,
                    config.allow_unchecked_headers,
                    config.anno_labels,
                    [&](std::string sequence, auto labels, uint64_t kmer_count, size_t) {
                        if (sequence.size() >= k) {
                            batcher.push_and_pay(sequence.size(),
                                                 std::move(sequence), std::move(labels),
                                                 std::vector<uint64_t>(sequence.size() - k + 1, kmer_count));
                        }
                    }
                );
            }
        }
    }

    anno_graph->get_annotator().serialize(annotator_filename);
}


int annotate_graph(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(config->infbase_annotators.size() <= 1);

    if (config->allow_unchecked_headers) {
        logger->warn("Header names are not checked (--allow-unchecked-headers): a label or"
                     " record header that is not valid UTF-8 cannot be listed by the serving"
                     " routes, and a repeated record header is found at its first record only");
    }

    if (config->index_header_coords) {
        for (const auto &file : files) {
            if (file_format(file) != "FASTA" && file_format(file) != "FASTQ") {
                logger->error("--index-header-coords requires FASTA or FASTQ input: {}", file);
                exit(1);
            }
        }
        if (!config->filename_anno || config->annotate_sequence_headers || config->anno_labels.size()) {
            logger->error("--index-header-coords requires --anno-filename "
                          "without --anno-header or --anno-label.");
            exit(1);
        }
        if (!config->fasta_header_delimiter.empty()) {
            logger->error("--header-delimiter cannot be used with --index-header-coords: the "
                          "record mapping (.seqs) names each record by its whole header, which "
                          "is not split");
            exit(1);
        }
        logger->trace("Parsing headers and computing coordinate offsets for {} files", files.size());
        logger->info("Note: It's essential that the files are passed in the order of the columns "
                     "in the graph annotation, under the names it labels them with (check with "
                     "`metagraph stats --print-col-names ...`): the mapping states each file's "
                     "name, and a load refuses a mapping whose names are not the annotation's");
        config->separately = false;
    }

    const auto graph = load_critical_dbg(config->infbase);

    if (!config->separately) {
        omp_set_max_active_levels(2);
        annotate_data(graph, *config, files, config->outfbase);

    } else {
        // |config->separately| is true

        // annotate multiple files in parallel, each with |parallel_each| threads
        size_t num_threads = std::max<size_t>(1, get_num_threads());
        set_num_threads(std::max(1u, config->parallel_each));
        omp_set_max_active_levels(3);

        // before anything is created: with several files an earlier file's annotation, with
        // an output directory the directory itself (a single file without one writes nothing
        // before its own check)
        if (config->annotate_sequence_headers && !config->allow_unchecked_headers
                && (files.size() > 1 || !config->outfbase.empty())) {
            check_header_labels(files, *config, num_threads);
        }

        if (!config->outfbase.empty()) {
            try {
                fs::create_directory(config->outfbase);
            } catch (...) {
                logger->error("Failed to create directory {}", config->outfbase);
                throw;
            }
        }

        #pragma omp parallel for num_threads(num_threads) default(shared) schedule(dynamic, 1)
        for (size_t i = 0; i < files.size(); ++i) {
            annotate_data(graph, *config, { files[i] },
                config->outfbase.size()
                    ? config->outfbase + "/" + utils::split_string(files[i], "/").back()
                    : files[i],
                2000 / num_threads);
        }
    }

    // Coordinate-aware annotations need a CoordToHeader (.seqs) mapping so that
    // `metagraph query --query-mode coords` and `metagraph align` can report
    // per-sequence positions and the per-target sequence length needed to
    // compute the fraction of the target covered. This mapping must be built
    // once against the FINAL merged annotation with the FULL set of input
    // FASTAs (run it after any transform/merge step). It is a separate step
    // because the column order may shift during transforms and a single
    // annotate call may only see a subset of the inputs.
    if (config->coordinates && config->filename_anno && !config->index_header_coords
            && !config->annotate_sequence_headers && config->anno_labels.empty()) {
        logger->info("To enable per-sequence coordinate reporting "
                     "(`<seq>/N:...` in TSV; JSON `kmers_in_target` for query, "
                     "`nt_length` for align), build the CoordToHeader mapping "
                     "once against the final merged annotation with all input "
                     "FASTAs in its column order (-a checks each file against its "
                     "column):\n"
                     "    metagraph annotate --anno-filename --index-header-coords "
                     "-i {} -a <final_annotation> -o {} <input_fastas>",
                     !config->infbase.empty() ? config->infbase : "<graph>",
                     !config->outfbase.empty() ? config->outfbase
                                               : "<final_anno_basename>");
    }

    return 0;
}

} // namespace cli
} // namespace mtg

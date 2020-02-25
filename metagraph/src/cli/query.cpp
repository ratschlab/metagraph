#include "query.hpp"

#include <ips4o.hpp>
#include <fmt/format.h>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/template_utils.hpp"
#include "common/threads/threading.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"
#include "align.hpp"

const size_t kRowBatchSize = 100'000;
const bool kPrefilterWithBloom = true;
const char ALIGNED_SEQ_HEADER_FORMAT[] = "{}:{}:{}:{}";

using mg::common::logger;


void execute_query(const std::string &seq_name,
                   const std::string &sequence,
                   bool count_labels,
                   bool print_signature,
                   bool suppress_unlabeled,
                   size_t num_top_labels,
                   double discovery_fraction,
                   std::string anno_labels_delimiter,
                   const AnnotatedDBG &anno_graph,
                   std::ostream &output_stream) {
    std::string output;
    output.reserve(1'000);

    if (print_signature) {
        auto top_labels
            = anno_graph.get_top_label_signatures(sequence,
                                                  num_top_labels,
                                                  discovery_fraction);

        if (!top_labels.size() && suppress_unlabeled)
            return;

        output += seq_name;

        for (const auto &[label, kmer_presence_mask] : top_labels) {
            output += fmt::format("\t<{}>:{}:{}:{}",
                label,
                sdsl::util::cnt_one_bits(kmer_presence_mask),
                sdsl::util::to_string(kmer_presence_mask),
                anno_graph.score_kmer_presence_mask(kmer_presence_mask)
            );
        }

        output += '\n';

    } else if (count_labels) {
        auto top_labels = anno_graph.get_top_labels(sequence,
                                                    num_top_labels,
                                                    discovery_fraction);

        if (!top_labels.size() && suppress_unlabeled)
            return;

        output += seq_name;

        for (const auto &[label, count] : top_labels) {
            output += "\t<";
            output += label;
            output += ">:";
            output += fmt::format_int(count).c_str();
        }

        output += '\n';

    } else {
        auto labels_discovered = anno_graph.get_labels(sequence, discovery_fraction);

        if (!labels_discovered.size() && suppress_unlabeled)
            return;

        output += seq_name;
        output += '\t';
        output += utils::join_strings(labels_discovered, anno_labels_delimiter);
        output += '\n';
    }

    output_stream << output;
}

std::unique_ptr<AnnotatedDBG>
construct_query_graph(const AnnotatedDBG &anno_graph,
                      StringGenerator call_sequences,
                      double discovery_fraction,
                      size_t num_threads) {
    const auto &full_dbg = anno_graph.get_graph();
    const auto &full_annotation = anno_graph.get_annotation();

    Timer timer;

    // construct graph storing all k-mers in query
    auto graph = std::make_shared<DBGHashOrdered>(full_dbg.get_k(), false);

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&full_dbg);
    if (kPrefilterWithBloom && dbg_succ) {
        if (dbg_succ->get_bloom_filter())
            logger->trace("[Query graph construction] Started indexing k-mers pre-filtered with Bloom filter");

        call_sequences([&graph,&dbg_succ](const std::string &sequence) {
            graph->add_sequence(sequence, get_missing_kmer_skipper(
                dbg_succ->get_bloom_filter(),
                sequence
            ));
        });
    } else {
        call_sequences([&graph](const std::string &sequence) {
            graph->add_sequence(sequence);
        });
    }

    logger->trace("[Query graph construction] k-mer indexing took {} sec", timer.elapsed());
    timer.reset();

    // pull contigs from query graph
    std::vector<std::pair<std::string, std::vector<DeBruijnGraph::node_index>>> contigs;
    graph->call_sequences(
        [&](auto&&... contig_args) { contigs.emplace_back(std::move(contig_args)...); },
        full_dbg.is_canonical_mode()
    );

    logger->trace("[Query graph construction] Contig extraction took {} sec", timer.elapsed());
    timer.reset();

    // map contigs onto the full graph
    auto index_in_full_graph = std::make_shared<std::vector<uint64_t>>();

    if (full_dbg.is_canonical_mode()) {
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 10)
        for (size_t i = 0; i < contigs.size(); ++i) {
            const std::string &contig = contigs[i].first;
            auto &nodes_in_full = contigs[i].second;

            size_t j = 0;
            full_dbg.map_to_nodes(contig,
                [&](auto node_in_full) { nodes_in_full[j++] = node_in_full; }
            );
            assert(j == nodes_in_full.size());
        }

        logger->trace("[Query graph construction] Contigs mapped to graph in {} sec", timer.elapsed());
        timer.reset();

        // construct canonical graph storing all k-mers found in the full graph
        graph = std::make_shared<DBGHashOrdered>(full_dbg.get_k(), true);

        for (size_t i = 0; i < contigs.size(); ++i) {
            const std::string &contig = contigs[i].first;
            const auto &nodes_in_full = contigs[i].second;
            size_t j = 0;
            graph->add_sequence(contig, [&]() { return nodes_in_full[j++] == 0; });
        }

        logger->trace("[Query graph construction] k-mers reindexed in canonical mode in {} sec",
                      timer.elapsed());
        timer.reset();

        index_in_full_graph->assign(graph->max_index() + 1, 0);

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 10)
        for (size_t i = 0; i < contigs.size(); ++i) {
            const std::string &contig = contigs[i].first;
            const auto &nodes_in_full = contigs[i].second;

            size_t j = 0;
            graph->map_to_nodes(contig,
                [&](auto node) { (*index_in_full_graph)[node] = nodes_in_full[j++]; }
            );
            assert(j == nodes_in_full.size());
        }

        logger->trace("[Query graph construction] Mapping between graphs constructed in {} sec",
                      timer.elapsed());
        timer.reset();

    } else {
        index_in_full_graph->assign(graph->max_index() + 1, 0);

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 10)
        for (size_t i = 0; i < contigs.size(); ++i) {
            const std::string &contig = contigs[i].first;
            const auto &path = contigs[i].second;

            size_t j = 0;
            full_dbg.map_to_nodes(contig,
                [&](auto node_in_full) { (*index_in_full_graph)[path[j++]] = node_in_full; }
            );
            assert(j == path.size());
        }

        logger->trace("[Query graph construction] Contigs mapped to graph in {} sec", timer.elapsed());
        timer.reset();
    }

    contigs = decltype(contigs)();

    assert(!(*index_in_full_graph)[0]);

    if (discovery_fraction > 0) {
        sdsl::bit_vector mask(graph->max_index() + 1, false);

        call_sequences([&](const std::string &sequence) {
            if (sequence.length() < graph->get_k())
                return;

            const size_t num_kmers = sequence.length() - graph->get_k() + 1;
            const size_t max_kmers_missing = num_kmers * (1 - discovery_fraction);
            const size_t min_kmers_discovered = num_kmers - max_kmers_missing;
            size_t num_kmers_discovered = 0;
            size_t num_kmers_missing = 0;

            std::vector<DeBruijnGraph::node_index> nodes;
            nodes.reserve(num_kmers);

            graph->map_to_nodes(sequence,
                [&](auto node) {
                    if ((*index_in_full_graph)[node]) {
                        num_kmers_discovered++;
                        nodes.push_back(node);
                    } else {
                        num_kmers_missing++;
                    }
                },
                [&]() { return num_kmers_missing > max_kmers_missing
                                || num_kmers_discovered >= min_kmers_discovered; }
            );

            if (num_kmers_missing <= max_kmers_missing) {
                for (auto node : nodes) { mask[node] = true; }
            }
        });

        // correcting the mask
        call_zeros(mask, [&](auto i) { (*index_in_full_graph)[i] = 0; });

        logger->trace("[Query graph construction] Reduced k-mer dictionary in {} sec",
                      timer.elapsed());
        timer.reset();
    }

    assert(index_in_full_graph.get());

    std::vector<std::pair<uint64_t, uint64_t>> from_full_to_query;
    from_full_to_query.reserve(index_in_full_graph->size());

    for (uint64_t node = 0; node < index_in_full_graph->size(); ++node) {
        if ((*index_in_full_graph)[node]) {
            from_full_to_query.emplace_back(
                AnnotatedDBG::graph_to_anno_index((*index_in_full_graph)[node]),
                AnnotatedDBG::graph_to_anno_index(node)
            );
        }
    }

    ips4o::parallel::sort(from_full_to_query.begin(), from_full_to_query.end(),
                          utils::LessFirst(), num_threads);

    logger->trace("[Query graph construction] Prepared row indexes for query {} sec",
                  timer.elapsed());
    timer.reset();

    // initialize fast query annotation
    // TODO: use SmallVector if it doesn't slow it down too much. Increase the batch size
    std::vector<BinaryMatrix::SetBitPositions> annotation_rows(graph->max_index());

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (uint64_t batch_begin = 0;
                        batch_begin < from_full_to_query.size();
                                        batch_begin += kRowBatchSize) {

        const uint64_t batch_end
            = std::min(batch_begin + kRowBatchSize,
                       static_cast<uint64_t>(from_full_to_query.size()));

        std::vector<uint64_t> row_indexes;
        row_indexes.reserve(batch_end - batch_begin);

        for (uint64_t i = batch_begin; i < batch_end; ++i) {
            assert(from_full_to_query[i].first < full_annotation.num_objects());

            row_indexes.push_back(from_full_to_query[i].first);
        }

        auto rows = full_annotation.get_matrix().get_rows(row_indexes);

        assert(rows.size() == batch_end - batch_begin);

        for (uint64_t i = batch_begin; i < batch_end; ++i) {
            annotation_rows[from_full_to_query[i].second]
                = std::move(rows[i - batch_begin]);
        }
    }

    // copy annotations from the full graph to the query graph
    auto annotation = std::make_unique<annotate::RowCompressed<>>(
        std::move(annotation_rows),
        full_annotation.get_label_encoder().get_labels()
    );

    logger->trace("[Query graph construction] Query annotation constructed in {} sec",
                  timer.elapsed());
    timer.reset();

    auto masked_graph = std::make_shared<MaskedDeBruijnGraph>(graph,
        [=](auto i) -> bool { return (*index_in_full_graph)[i]; }
    );

    // build annotated graph from the query graph and copied annotations
    return std::make_unique<AnnotatedDBG>(masked_graph, std::move(annotation));
}


std::string get_alignment_header_and_swap_query(const std::string &name,
                                                std::string *query_seq,
                                                QueryAlignment<> *matches) {
    std::string header;

    if (matches->size()) {
        // sequence for querying -- the best alignment
        *query_seq = const_cast<std::string&&>((*matches)[0].get_sequence());
        header = fmt::format(ALIGNED_SEQ_HEADER_FORMAT,
                             name, *query_seq, (*matches)[0].get_score(),
                             (*matches)[0].get_cigar().to_string());

    } else {
        // no alignment was found
        // the original sequence `query_seq` will be queried
        header = fmt::format(ALIGNED_SEQ_HEADER_FORMAT,
                             name, *query_seq, 0,
                             fmt::format("{}S", query_seq->length()));
    }

    return header;
}

int query_graph(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(config->infbase_annotators.size() == 1);

    auto graph = load_critical_dbg(config->infbase);
    auto anno_graph = initialize_annotated_dbg(graph, *config);

    ThreadPool thread_pool(std::max(1u, get_num_threads()) - 1, 1000);

    Timer timer;

    std::unique_ptr<IDBGAligner> aligner;
    if (config->align_sequences) {
        assert(config->alignment_num_alternative_paths == 1u
                && "only the best alignment is used in query");

        aligner = build_aligner(*graph, *config);
    }

    // iterate over input files
    for (const auto &file : files) {
        logger->trace("Parsing sequences from file '{}'", file);

        Timer curr_timer;

        size_t seq_count = 0;

        const auto *graph_to_query = anno_graph.get();

        auto execute = [&](size_t id, const std::string &name, const std::string &seq) {
            execute_query(fmt::format_int(id).str() + '\t' + name,
                          seq,
                          config->count_labels,
                          config->print_signature,
                          config->suppress_unlabeled,
                          config->num_top_labels,
                          config->discovery_fraction,
                          config->anno_labels_delimiter,
                          std::ref(*graph_to_query),
                          std::ref(std::cout));
        };

        FastaParser fasta_parser(file, config->forward_and_reverse);

        // Graph constructed from a batch of queried sequences
        // Used only in fast mode
        if (config->fast) {
            auto begin = fasta_parser.begin();
            auto end = fasta_parser.end();

            const uint64_t batch_size = config->query_batch_size_in_bytes;
            FastaParser::iterator it;

            while (begin != end) {
                Timer batch_timer;

                std::vector<std::tuple<size_t, std::string, std::string>> named_alignments;
                uint64_t num_bytes_read = 0;
                auto query_graph = construct_query_graph(*anno_graph,
                    [&](auto call_sequence) {
                        num_bytes_read = 0;

                        if (!aligner) {
                            // basic query regime
                            // the query graph is constructed directly from the input sequences
                            for (it = begin; it != end && num_bytes_read <= batch_size; ++it) {
                                call_sequence(it->seq.s);
                                num_bytes_read += it->seq.l;
                            }
                            return;
                        }

                        // Query with alignment to graph

                        // Check if the sequences have already been aligned and
                        // if the results are stored in |named_alignments|.
                        if (named_alignments.size()) {
                            for (const auto &[id, name, seq] : named_alignments) {
                                call_sequence(seq);
                            }
                            return;
                        }

                        std::mutex sequence_mutex;
                        for ( ; begin != end && num_bytes_read <= batch_size; ++begin) {
                            // Align the sequence, then add the best match
                            // in the graph to the query graph.
                            thread_pool.enqueue(
                                [&](size_t id, const std::string &name, std::string &seq) {
                                    auto matches = aligner->align(seq);

                                    std::string seq_header
                                        = get_alignment_header_and_swap_query(name, &seq, &matches);

                                    std::lock_guard<std::mutex> lock(sequence_mutex);

                                    named_alignments.emplace_back(
                                        id, std::move(seq_header), std::move(seq));

                                    call_sequence(std::get<2>(named_alignments.back()));
                                },
                                seq_count++, std::string(begin->name.s), std::string(begin->seq.s)
                            );

                            num_bytes_read += begin->seq.l;
                        }
                        thread_pool.join();
                    },
                    config->count_labels ? 0 : config->discovery_fraction,
                    get_num_threads()
                );

                graph_to_query = query_graph.get();

                logger->trace("Query graph constructed for batch of {} bytes from '{}' in {} sec",
                              num_bytes_read, file, batch_timer.elapsed());

                batch_timer.reset();

                if (!aligner) {
                    for ( ; begin != it; ++begin) {
                        assert(begin != end);

                        thread_pool.enqueue(execute, seq_count++,
                                                     std::string(begin->name.s),
                                                     std::string(begin->seq.s));
                    }
                } else {
                    for (auto&& [id, name, seq] : named_alignments) {
                        thread_pool.enqueue(execute, id, std::move(name), std::move(seq));
                    }
                }

                thread_pool.join();

                logger->trace("Batch of {} bytes from '{}' queried in {} sec",
                              num_bytes_read, file, batch_timer.elapsed());
            }

        } else {
            for (const auto &kseq : fasta_parser) {
                thread_pool.enqueue(
                    [&](size_t id, const std::string &name, std::string &seq) {
                        if (!aligner) {
                            execute(id, name, seq);
                            return;
                        }

                        // query the alignment matches against the annotator
                        auto matches = aligner->align(seq);

                        std::string seq_header
                            = get_alignment_header_and_swap_query(name, &seq, &matches);

                        execute(id, seq_header, seq);
                    },
                    seq_count++, std::string(kseq.name.s), std::string(kseq.seq.s)
                );
            }

            // wait while all threads finish processing the current file
            thread_pool.join();
        }

        logger->trace("File '{}' was processed in {} sec, total time: {}", file,
                      curr_timer.elapsed(), timer.elapsed());
    }

    return 0;
}

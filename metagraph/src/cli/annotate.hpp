#ifndef __ANNOTATE_GRAPH_HPP__
#define __ANNOTATE_GRAPH_HPP__

#include <string>
#include <vector>

#include "seq_io/sequence_io.hpp"


namespace mtg {
namespace cli {

class Config;

/**
 * The labels a FASTA header yields, appended to |labels|: the name token, or name and comment
 * joined by |comment_delim| when set (Config::UNINITIALIZED_STR: never), split by
 * |header_delim|. A record without a comment is its name alone, whatever kseq's comment buffer
 * still holds from an earlier record. A label that is not valid UTF-8 refuses the build,
 * naming |file|, the 1-based |record| and the bytes in hex (never raw): the serving routes
 * cannot list such a label, since jsoncpp would replace its bytes by U+FFFD and a valid name
 * holding U+FFFD could be another label's. |unchecked| skips the check.
 */
void push_header_labels(const seq_io::kseq_t *read_stream,
                        const std::string &file,
                        size_t record,
                        const std::string &comment_delim,
                        const std::string &header_delim,
                        bool unchecked,
                        std::vector<std::string> *labels);

int annotate_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __ANNOTATE_GRAPH_HPP__

#include "dbg_aligner.hpp"

namespace mtg {
namespace graph {
namespace align {

IDBGAligner::DBGQueryAlignment IDBGAligner::align(std::string_view query,
                                                  bool is_reverse_complement) const {
    DBGQueryAlignment result(query);
    std::string empty_header;
    align_batch(
        [&](const QueryCallback &callback) {
            callback(empty_header, query, is_reverse_complement);
        },
        [&](std::string_view, DBGQueryAlignment&& alignment) {
            result = std::move(alignment);
        }
    );

    return result;
}

} // namespace align
} // namespace graph
} // namespace mtg

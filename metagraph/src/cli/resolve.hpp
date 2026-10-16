#ifndef __CLI_RESOLVE_HPP__
#define __CLI_RESOLVE_HPP__
#include <json/json.h>
namespace mtg::graph { class AnnotatedDBG; }
namespace mtg::cli {
class Config;
Json::Value process_resolve_request(const Json::Value &, const graph::AnnotatedDBG &);
int resolve_graph(Config *);
}
#endif

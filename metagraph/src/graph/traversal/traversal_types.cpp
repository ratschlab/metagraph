#include "traversal_types.hpp"


namespace mtg {
namespace graph {
namespace traversal {

const char* to_string(Regime regime) {
    switch (regime) {
        case Regime::BASIC: return "basic";
        case Regime::PRIMARY: return "primary";
        case Regime::CANONICAL: return "canonical";
    }
    return "unknown";
}

const char* to_string(LabelKind kind) {
    switch (kind) {
        case LabelKind::COLUMN: return "column";
        case LabelKind::HEADER: return "header";
    }
    return "unknown";
}

} // namespace traversal
} // namespace graph
} // namespace mtg

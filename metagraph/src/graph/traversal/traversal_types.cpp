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

const char* to_string(Support support) {
    switch (support) {
        case Support::KMER: return "kmer";
        case Support::TRACE: return "trace";
    }
    return "unknown";
}

const char* to_string(Arm arm) {
    return arm == Arm::LEFT ? "left" : "right";
}

const char* to_string(EndReason reason) {
    switch (reason) {
        case EndReason::DEAD_END: return "dead_end";
        case EndReason::LABEL_LOST: return "label_lost";
        case EndReason::LOSS_BUDGET: return "loss_budget";
        case EndReason::BRANCH: return "branch";
        case EndReason::EDGE_REUSE: return "edge_reuse";
        case EndReason::EDGE_REUSE_RC: return "edge_reuse_rc";
        case EndReason::REACHED_SEED: return "rejoined_seed";
        case EndReason::RECORD_END: return "trace_break";
        case EndReason::MAX_EXTENSION: return "max_extension_bp";
        case EndReason::MAX_STEPS: return "max_steps";
        case EndReason::MAX_LIVE_PATHS: return "max_live_paths";
        case EndReason::MAX_PATHS: return "max_paths";
        case EndReason::MAX_OUTPUT: return "max_output_bp";
        case EndReason::TIME_BUDGET: return "time_budget";
        case EndReason::BEAM_PRUNED: return "beam_pruned";
    }
    return "unknown";
}

bool is_resource_stop(EndReason reason) {
    switch (reason) {
        case EndReason::MAX_STEPS:
        case EndReason::MAX_LIVE_PATHS:
        case EndReason::MAX_PATHS:
        case EndReason::MAX_OUTPUT:
        case EndReason::TIME_BUDGET:
        case EndReason::BEAM_PRUNED:
            return true;
        default:
            return false;
    }
}

} // namespace traversal
} // namespace graph
} // namespace mtg

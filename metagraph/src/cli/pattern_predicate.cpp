#include "pattern_predicate.hpp"

#include <set>
#include <stdexcept>
#include <string_view>

#include <tsl/hopscotch_map.h>

#include "graph/alignment/pattern_search.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "pattern.hpp"


namespace mtg {
namespace cli {
namespace predicate {

using graph::pattern::Budget;

namespace {

// a node's value: false, true, or (three-valued evaluation only) not known
constexpr uint8_t kFalse = 0;
constexpr uint8_t kTrue = 1;
constexpr uint8_t kUnknown = 2;

// a label's mark during an evaluation
constexpr uint8_t kSure = 1;
constexpr uint8_t kMaybe = 2;

// Binding reads the work time before every this many names resolved, nodes and items folded
// and built (§19.9: "every 4,096 names bound"); the echo calls the answer's delivery check
// before every this many names and nodes written
constexpr uint64_t kStride = 4096;

constexpr const char kOperators[] = "any, all, none, at_least, and, or, not";

PatternRefusal invalid(const std::string &message) {
    return PatternRefusal(400, "invalid_request", message);
}

/**
 * Strict access to one JSON object of the request, as /traverse's: every field read is
 * remembered, and finish() refuses the first one nothing read ("unknown field"), after the
 * known ones were checked — the guarantee rule: a field this increment does not know is never
 * ignored. (A copy of the route's reader in pattern.cpp, kept here so that this increment does
 * not touch the route; the two can move into one header with the route's predicate fields.)
 */
class Fields {
  public:
    Fields(const Json::Value &value, std::string path) : v_(value), path_(std::move(path)) {
        if (!v_.isObject())
            throw invalid(path_ + ": expected an object");
    }

    bool has(const std::string &key) {
        seen_.insert(key);
        return v_.isMember(key);
    }
    const Json::Value& raw(const std::string &key) {
        seen_.insert(key);
        return v_[key];
    }
    std::string path(const std::string &key) const { return path_ + "." + key; }

    void finish() const {
        for (const std::string &name : v_.getMemberNames()) {
            if (!seen_.count(name))
                throw invalid(path_ + ": unknown field '" + name + "'");
        }
    }

  private:
    const Json::Value &v_;
    std::string path_;
    std::set<std::string> seen_;
};

bool is_leaf(Op op) { return op <= Op::AT_LEAST; }

std::optional<Op> op_of(const std::string &key) {
    for (Op op : { Op::ANY, Op::ALL, Op::NONE, Op::AT_LEAST, Op::AND, Op::OR, Op::NOT }) {
        if (key == to_string(op))
            return op;
    }
    return std::nullopt;
}

std::string item_path(const std::string &path, size_t i) {
    return path + "[" + std::to_string(i) + "]";
}

/**
 * The parse of one predicate (§19.3): a walk of the JSON that checks every rule, the first
 * broken one refused, and builds the nodes children first. A name is seen through a view into
 * the JSON (no copy): the views deduplicate the names of the whole predicate, and the list a
 * name was last listed in tells a name listed twice in one list. The names are copied once
 * the walk found the predicate valid and within the cap. Past the cap's |max_labels| names
 * the walk only checks each name's form and counts it: the predicate is refused for its size
 * whatever the rest lists (a later list repeating a name included), and the index holds at
 * most |max_labels| names, whatever the body's size.
 */
class Parser {
  public:
    Parser(Predicate *out, uint64_t max_labels) : out_(*out), max_labels_(max_labels) {}

    uint32_t node(const Json::Value &v, const std::string &path, unsigned depth);

    void finish(const std::string &path) {
        if (out_.listed > max_labels_) {
            throw PatternRefusal(400, "predicate_too_large",
                path + ": " + std::to_string(out_.listed) + " names in its lists (a name in "
                "two lists counts twice), above the server's cap of "
                + std::to_string(max_labels_) + " (caps.max_predicate_labels): send at most "
                + std::to_string(max_labels_) + " names; split the cohort");
        }
        out_.names.reserve(distinct_.size());
        for (std::string_view name : distinct_) {
            out_.names.emplace_back(name);
        }
    }

  private:
    Predicate &out_;
    const uint64_t max_labels_;
    tsl::hopscotch_map<std::string_view, uint32_t> index_;
    std::vector<std::string_view> distinct_;
    // per distinct name: the list it was last listed in (1-based, 0: none yet), and where
    std::vector<uint32_t> last_list_;
    std::vector<uint32_t> last_pos_;
    uint32_t lists_ = 0;

    void names(const Json::Value &list, const std::string &path, Predicate::Node *node);
};

void Parser::names(const Json::Value &list, const std::string &path, Predicate::Node *node) {
    if (!list.isArray())
        throw invalid(path + ": expected a list of column labels");
    if (list.empty())
        throw invalid(path + ": expected a list of at least one column label");

    ++lists_;
    node->begin = out_.items.size();
    for (Json::ArrayIndex i = 0; i < list.size(); ++i) {
        const Json::Value &name = list[i];
        if (!name.isString()) {
            // P7: names are strings, a taxid among them; the message names the fix
            if (name.isNumeric()) {
                const std::string example = !name.isIntegral() ? std::string("562")
                        : name.isInt64() ? std::to_string(name.asInt64())
                                         : std::to_string(name.asUInt64());
                throw invalid(item_path(path, i) + ": expected a string (a column label; "
                              "write a taxid as \"" + example + "\")");
            }
            throw invalid(item_path(path, i) + ": expected a string (a column label)");
        }
        const char *begin = nullptr;
        const char *end = nullptr;
        if (!name.getString(&begin, &end) || begin == end)
            throw invalid(item_path(path, i) + ": expected a non-empty string (a column label)");

        if (out_.listed++ >= max_labels_)
            continue;

        const std::string_view view(begin, end - begin);
        auto [it, inserted] = index_.try_emplace(view, distinct_.size());
        const uint32_t id = it->second;
        if (inserted) {
            distinct_.push_back(view);
            last_list_.push_back(0);
            last_pos_.push_back(0);
        } else if (last_list_[id] == lists_) {
            throw invalid(item_path(path, i) + ": '" + std::string(view) + "' listed twice in "
                          "one list (first at [" + std::to_string(last_pos_[id]) + "])");
        }
        last_list_[id] = lists_;
        last_pos_[id] = i;
        out_.items.push_back(id);
    }
    node->end = out_.items.size();
}

uint32_t Parser::node(const Json::Value &v, const std::string &path, unsigned depth) {
    if (!v.isObject()) {
        throw invalid(path + ": expected a predicate: an object with one operator ("
                      + kOperators + ")");
    }
    if (v.size() != 1) {
        if (v.empty())
            throw invalid(path + ": expected one operator (" + kOperators + "), got {}");
        auto it = v.begin();
        const std::string first = it.name();
        const std::string second = (++it).name();
        throw invalid(path + ": expected exactly one operator (" + kOperators + "), got "
                      + std::to_string(v.size()) + " members ('" + first + "', '" + second
                      + "'" + (v.size() > 2 ? ", ...)" : ")"));
    }
    const std::string key = v.begin().name();
    const std::optional<Op> op = op_of(key);
    if (!op) {
        throw invalid(path + ": unknown operator '" + key + "' (expected one of "
                      + kOperators + ")");
    }
    const Json::Value &arg = *v.begin();
    const std::string at = path + "." + key;

    Predicate::Node node;
    node.op = *op;
    switch (*op) {
        case Op::ANY:
        case Op::ALL:
        case Op::NONE:
            names(arg, at, &node);
            break;
        case Op::AT_LEAST: {
            if (!arg.isObject()) {
                throw invalid(at + ": expected an object {\"n\": <integer>, \"labels\": "
                              "[<column labels>]}");
            }
            Fields f(arg, at);
            if (!f.has("n"))
                throw invalid(f.path("n") + ": required: an integer from 1 to the labels listed");
            const Json::Value &n = f.raw("n");
            if (!n.isIntegral())
                throw invalid(f.path("n") + ": expected an integer from 1 to the labels listed");
            if (!f.has("labels"))
                throw invalid(f.path("labels") + ": required: a list of column labels");
            names(f.raw("labels"), f.path("labels"), &node);
            // the list's size: past the cap its names are counted, not kept
            const uint64_t listed = f.raw("labels").size();
            if ((n.isInt64() && n.asInt64() < 1) || n.asUInt64() > listed) {
                throw invalid(f.path("n") + ": expected an integer from 1 to "
                              + std::to_string(listed) + " (the labels listed), got "
                              + (n.isInt64() ? std::to_string(n.asInt64())
                                             : std::to_string(n.asUInt64())));
            }
            f.finish();
            node.n = n.asUInt64();
            break;
        }
        case Op::AND:
        case Op::OR: {
            if (depth >= kMaxNesting) {
                throw invalid(path + ": predicates nested more than "
                              + std::to_string(kMaxNesting) + " deep (and, or, not)");
            }
            if (!arg.isArray())
                throw invalid(at + ": expected a list of predicates");
            if (arg.empty())
                throw invalid(at + ": expected a list of at least one predicate");
            // the operands first (children before their parent), then their ids, contiguous
            std::vector<uint32_t> operands;
            operands.reserve(arg.size());
            for (Json::ArrayIndex i = 0; i < arg.size(); ++i) {
                operands.push_back(this->node(arg[i], item_path(at, i), depth + 1));
            }
            node.begin = out_.items.size();
            out_.items.insert(out_.items.end(), operands.begin(), operands.end());
            node.end = out_.items.size();
            break;
        }
        case Op::NOT: {
            if (depth >= kMaxNesting) {
                throw invalid(path + ": predicates nested more than "
                              + std::to_string(kMaxNesting) + " deep (and, or, not)");
            }
            const uint32_t operand = this->node(arg, at, depth + 1);
            node.begin = out_.items.size();
            out_.items.push_back(operand);
            node.end = out_.items.size();
            break;
        }
    }
    out_.nodes.push_back(node);
    return out_.nodes.size() - 1;
}

// The bytes of |name| written as a JSON string by jsoncpp, quotes included, at most: exact for
// printable ASCII, a character it escapes as its escape, a byte outside ASCII as 6 (a \u
// escape; a code point of several bytes is written in fewer)
uint64_t quoted_bytes(const std::string &name) {
    uint64_t bytes = 2;
    for (unsigned char c : name) {
        if (c == '"' || c == '\\' || c == '\b' || c == '\f' || c == '\n' || c == '\r'
                || c == '\t') {
            bytes += 2;
        } else if (c < 0x20 || c >= 0x80) {
            bytes += 6;
        } else {
            bytes += 1;
        }
    }
    return bytes;
}

// The bytes of {"<key>":<value>} in compact JSON
uint64_t member_bytes(const char *key, uint64_t value_bytes) {
    return std::char_traits<char>::length(key) + 5 + value_bytes;
}

/**
 * The work time read before every kStride units of binding (a name resolved, a node or an
 * item folded or built): false once it passed, the budget then recording TIME (sticky).
 */
class BindClock {
  public:
    explicit BindClock(Budget &budget) : budget_(budget) {}

    bool tick() {
        if (!left_) {
            if (!budget_.check_time())
                return false;
            left_ = kStride;
        }
        --left_;
        return true;
    }

  private:
    Budget &budget_;
    uint64_t left_ = 0;
};

} // namespace


const char* to_string(Op op) {
    switch (op) {
        case Op::ANY: return "any";
        case Op::ALL: return "all";
        case Op::NONE: return "none";
        case Op::AT_LEAST: return "at_least";
        case Op::AND: return "and";
        case Op::OR: return "or";
        case Op::NOT: return "not";
    }
    return "";
}

Predicate Predicate::parse(const Json::Value &json, uint64_t max_labels,
                           const std::string &path) {
    Predicate predicate;
    Parser parser(&predicate, max_labels);
    parser.node(json, path, 0);
    parser.finish(path);
    return predicate;
}


// ---------------------------------------------------------------- binding

namespace {

// the state of a parsed node folded: a constant, or the node standing for it in the normal
// form (itself, or the one operand left of an and or an or)
constexpr uint32_t kStateFalse = UINT32_MAX;
constexpr uint32_t kStateTrue = UINT32_MAX - 1;

bool is_node(uint32_t state) { return state < kStateTrue; }

// a parsed node or a name not (yet) in the normal form, or reached by it
constexpr uint32_t kNone = UINT32_MAX;
constexpr uint32_t kReached = UINT32_MAX - 1;

} // namespace

/**
 * bind() in four steps, each polled: the names resolved; the predicate folded (SPEC §19.4),
 * children first, each node's state the constant it folds to or the node that stands for it
 * in the normal form; what the root's node reaches, counted (its nodes, the names its leaves
 * list, its labels, its bytes) without copying anything; then, once admitted, the bound
 * predicate built from those states.
 */
Binding Bound::bind(const Predicate &predicate, const ColumnLookup &find_column,
                    Budget &budget, const std::function<bool(uint64_t)> &admit) {
    Binding result;
    BindClock clock(budget);
    auto stopped = [&]() {
        result.stop = "time";
        return std::move(result);
    };

    const size_t num_names = predicate.names.size();
    const size_t num_parsed = predicate.nodes.size();
    if (!num_parsed)
        throw std::invalid_argument("predicate: bind() of a predicate never parsed");

    // the names
    std::vector<std::optional<Column>> column(num_names);
    uint64_t num_known = 0;
    for (size_t i = 0; i < num_names; ++i) {
        if (!clock.tick())
            return stopped();
        column[i] = find_column(predicate.names[i]);
        num_known += column[i].has_value();
    }

    // the folding (the table of §19.4), children first
    std::vector<uint32_t> state(num_parsed);
    for (uint32_t v = 0; v < num_parsed; ++v) {
        if (!clock.tick())
            return stopped();
        const Predicate::Node &node = predicate.nodes[v];
        const uint32_t *begin = predicate.items.data() + node.begin;
        const uint32_t *end = predicate.items.data() + node.end;
        if (is_leaf(node.op)) {
            uint64_t known = 0;
            for (const uint32_t *it = begin; it != end; ++it) {
                if (!clock.tick())
                    return stopped();
                known += column[*it].has_value();
            }
            const uint64_t listed = end - begin;
            switch (node.op) {
                case Op::ANY: state[v] = known ? v : kStateFalse; break;
                case Op::ALL: state[v] = known == listed ? v : kStateFalse; break;
                case Op::NONE: state[v] = known ? v : kStateTrue; break;
                case Op::AT_LEAST: state[v] = known >= node.n ? v : kStateFalse; break;
                default: break;
            }
        } else if (node.op == Op::NOT) {
            const uint32_t s = state[*begin];
            state[v] = s == kStateFalse ? kStateTrue : s == kStateTrue ? kStateFalse : v;
        } else {
            // and: false if an operand is false, the true ones dropped; or: the dual
            const uint32_t absorbing = node.op == Op::AND ? kStateFalse : kStateTrue;
            const uint32_t neutral = node.op == Op::AND ? kStateTrue : kStateFalse;
            uint64_t left = 0;
            uint32_t last = kNone;
            bool absorbed = false;
            for (const uint32_t *it = begin; it != end && !absorbed; ++it) {
                if (!clock.tick())
                    return stopped();
                const uint32_t s = state[*it];
                if (s == absorbing) {
                    absorbed = true;
                } else if (s != neutral) {
                    ++left;
                    last = s;
                }
            }
            state[v] = absorbed ? absorbing
                     : left == 0 ? neutral
                     : left == 1 ? last
                     : v;
        }
    }

    // what the normal form holds: the nodes the root's node reaches (parents first, so
    // descending), the names their leaves list and the labels among them
    const uint32_t root_state = state[num_parsed - 1];
    std::vector<uint32_t> new_id(num_parsed, kNone);
    // per name: kReached when the normal form lists it, then its LabelId
    std::vector<uint32_t> label_of(num_names, kNone);
    uint64_t num_nodes = 0;
    uint64_t num_listed = 0;
    uint64_t num_items = 0;
    if (is_node(root_state)) {
        new_id[root_state] = kReached;
        for (uint32_t v = root_state + 1; v-- > 0; ) {
            if (!clock.tick())
                return stopped();
            if (new_id[v] != kReached)
                continue;
            ++num_nodes;
            const Predicate::Node &node = predicate.nodes[v];
            for (uint32_t i = node.begin; i < node.end; ++i) {
                if (!clock.tick())
                    return stopped();
                const uint32_t item = predicate.items[i];
                if (is_leaf(node.op)) {
                    if (column[item]) {
                        ++num_listed;
                        ++num_items;
                        label_of[item] = kReached;
                    }
                } else if (is_node(state[item])) {
                    ++num_items;
                    new_id[state[item]] = kReached;
                }
            }
        }
    }

    // the bytes of the model (Bound::bytes), admitted before anything is copied
    uint64_t bytes = 64 * num_nodes + 8 * num_listed;
    uint64_t num_labels = 0;
    uint64_t num_unknown = 0;
    for (size_t i = 0; i < num_names; ++i) {
        if (!clock.tick())
            return stopped();
        if (label_of[i] == kReached) {
            bytes += 192 + 2 * predicate.names[i].size();
            ++num_labels;
        } else if (!column[i]) {
            bytes += 64 + predicate.names[i].size();
            ++num_unknown;
        }
    }
    result.bytes = bytes;
    if (admit && !admit(bytes)) {
        result.stop = "max_memory";
        return result;
    }
    result.admitted = admit ? bytes : 0;

    // the bound predicate
    Bound b;
    b.num_names_ = num_names;
    b.num_known_ = num_known;
    b.bytes_ = bytes;
    b.unknown_.reserve(num_unknown);
    b.labels_.reserve(num_labels);
    // the text of the echo: the unknown labels' list and the normal form
    uint64_t unknown_text = 2 + (num_unknown ? num_unknown - 1 : 0);
    std::vector<uint64_t> label_text;
    label_text.reserve(num_labels);
    for (size_t i = 0; i < num_names; ++i) {
        if (!clock.tick())
            return stopped();
        if (label_of[i] == kReached) {
            label_of[i] = b.labels_.size();
            LabelRef ref;
            ref.kind = graph::traversal::LabelKind::COLUMN;
            ref.column = *column[i];
            ref.name = predicate.names[i];
            label_text.push_back(quoted_bytes(ref.name));
            b.labels_.push_back(std::move(ref));
        } else if (!column[i]) {
            b.unknown_.push_back(predicate.names[i]);
            unknown_text += quoted_bytes(b.unknown_.back());
        }
    }

    if (!is_node(root_state)) {
        b.constant_ = root_state == kStateTrue;
        b.vacuous_ = *b.constant_;
        b.monotone_ = true;
        b.leaves_begin_.assign(1, 0);
        b.text_bytes_ = unknown_text + (*b.constant_ ? 4 : 5);
        result.bound = std::move(b);
        return result;
    }

    // the nodes, children first (the parsed order restricted to the reached ones), each
    // with the text of its JSON
    b.nodes_.reserve(num_nodes);
    b.items_.reserve(num_items);
    std::vector<uint64_t> text;
    text.reserve(num_nodes);
    for (uint32_t v = 0; v <= root_state; ++v) {
        if (!clock.tick())
            return stopped();
        if (new_id[v] != kReached)
            continue;
        const Predicate::Node &parsed = predicate.nodes[v];
        const uint32_t id = b.nodes_.size();
        new_id[v] = id;
        Bound::Node node;
        node.op = parsed.op;
        node.n = parsed.n;
        node.begin = b.items_.size();
        node.parent = Bound::kNoParent;
        uint64_t value_text = 2;   // the brackets of a list
        for (uint32_t i = parsed.begin; i < parsed.end; ++i) {
            if (!clock.tick())
                return stopped();
            const uint32_t item = predicate.items[i];
            if (is_leaf(parsed.op)) {
                if (column[item]) {
                    value_text += label_text[label_of[item]];
                    b.items_.push_back(label_of[item]);
                }
            } else if (is_node(state[item])) {
                const uint32_t operand = new_id[state[item]];
                b.nodes_[operand].parent = id;
                value_text += text[operand];
                b.items_.push_back(operand);
            }
        }
        node.end = b.items_.size();
        value_text += node.end - node.begin - 1;   // the commas
        switch (node.op) {
            case Op::AT_LEAST:
                // {"labels":[...],"n":<n>}
                value_text += 16 + std::to_string(node.n).size();
                break;
            case Op::NOT:
                value_text = text[b.items_[node.begin]];
                break;
            default:
                break;
        }
        text.push_back(member_bytes(to_string(node.op), value_text));
        b.monotone_ = b.monotone_ && node.op != Op::NONE && node.op != Op::NOT;
        b.nodes_.push_back(node);
    }
    b.text_bytes_ = unknown_text + text.back();

    // each label's leaves
    b.leaves_begin_.assign(num_labels + 1, 0);
    for (const Bound::Node &node : b.nodes_) {
        if (!is_leaf(node.op))
            continue;
        for (uint32_t i = node.begin; i < node.end; ++i) {
            if (!clock.tick())
                return stopped();
            ++b.leaves_begin_[b.items_[i] + 1];
        }
    }
    for (size_t l = 0; l < num_labels; ++l) {
        b.leaves_begin_[l + 1] += b.leaves_begin_[l];
    }
    b.leaves_.resize(num_listed);
    std::vector<uint32_t> fill(b.leaves_begin_.begin(), b.leaves_begin_.end() - 1);
    for (uint32_t id = 0; id < b.nodes_.size(); ++id) {
        const Bound::Node &node = b.nodes_[id];
        if (!is_leaf(node.op))
            continue;
        for (uint32_t i = node.begin; i < node.end; ++i) {
            if (!clock.tick())
                return stopped();
            b.leaves_[fill[b.items_[i]]++] = id;
        }
    }

    // the values on a context carrying none of the labels, children first
    b.base_value_.resize(num_nodes);
    b.base_true_.assign(num_nodes, 0);
    b.base_false_.assign(num_nodes, 0);
    for (uint32_t id = 0; id < num_nodes; ++id) {
        if (!clock.tick())
            return stopped();
        const Bound::Node &node = b.nodes_[id];
        if (is_leaf(node.op)) {
            b.base_value_[id] = b.leaf_value(node, 0, 0);
            continue;
        }
        for (uint32_t i = node.begin; i < node.end; ++i) {
            (b.base_value_[b.items_[i]] == kTrue ? b.base_true_ : b.base_false_)[id]++;
        }
        b.base_value_[id] = b.combinator_value(node, b.base_true_[id], b.base_false_[id]);
    }
    b.vacuous_ = b.base_value_.back() == kTrue;

    // the scratch, equal to the base
    b.value_ = b.base_value_;
    b.count_a_ = b.base_true_;
    b.count_b_ = b.base_false_;
    b.dirty_.assign(num_nodes, 0);
    b.dirty_list_.reserve(num_nodes);
    b.mark_.assign(num_labels, 0);
    b.marked_.reserve(num_labels);

    result.bound = std::move(b);
    return result;
}

Binding Bound::bind(const Predicate &predicate,
                    const graph::traversal::LabelOracle &oracle,
                    Budget &budget,
                    const std::function<bool(uint64_t)> &admit) {
    return bind(predicate,
                [&oracle](const std::string &name) { return oracle.find_column(name); },
                budget, admit);
}


// ---------------------------------------------------------------- evaluation

uint8_t Bound::leaf_value(const Node &node, uint32_t sure, uint32_t maybe) const {
    const uint32_t listed = node.end - node.begin;
    switch (node.op) {
        case Op::ANY:
            return sure ? kTrue : !maybe ? kFalse : kUnknown;
        case Op::ALL:
            return sure == listed ? kTrue : sure + maybe < listed ? kFalse : kUnknown;
        case Op::NONE:
            return sure ? kFalse : !maybe ? kTrue : kUnknown;
        case Op::AT_LEAST:
            return sure >= node.n ? kTrue : sure + maybe < node.n ? kFalse : kUnknown;
        default:
            throw std::logic_error("predicate: a combinator evaluated as a leaf");
    }
}

uint8_t Bound::combinator_value(const Node &node, uint32_t num_true,
                                uint32_t num_false) const {
    const uint32_t operands = node.end - node.begin;
    switch (node.op) {
        case Op::AND:
            return num_false ? kFalse : num_true == operands ? kTrue : kUnknown;
        case Op::OR:
            return num_true ? kTrue : num_false == operands ? kFalse : kUnknown;
        case Op::NOT:
            return num_true ? kFalse : num_false ? kTrue : kUnknown;
        default:
            throw std::logic_error("predicate: a leaf evaluated as a combinator");
    }
}

uint64_t Bound::mark(const LabelId *sure, size_t num_sure, const LabelId *maybe,
                     size_t num_maybe) const {
    uint64_t units = 1;
    auto add = [&](LabelId label, uint8_t kind) {
        if (label >= labels_.size()) {
            unmark();
            throw std::invalid_argument("predicate: label id " + std::to_string(label)
                                        + " outside the permitted set of "
                                        + std::to_string(labels_.size()));
        }
        if (mark_[label])
            return;
        mark_[label] = kind;
        marked_.push_back(label);
        units += leaves_begin_[label + 1] - leaves_begin_[label];
    };
    for (size_t i = 0; i < num_sure; ++i) {
        add(sure[i], kSure);
    }
    for (size_t i = 0; i < num_maybe; ++i) {
        add(maybe[i], kMaybe);
    }
    return units;
}

void Bound::unmark() const {
    for (LabelId label : marked_) {
        mark_[label] = 0;
    }
    marked_.clear();
}

// |node| takes |value|; its parent's counts follow, and so on up while a value changes
void Bound::set_value(uint32_t node, uint8_t value) const {
    while (value_[node] != value) {
        const uint8_t old = value_[node];
        value_[node] = value;
        if (!dirty_[node]) {
            dirty_[node] = 1;
            dirty_list_.push_back(node);
        }
        const uint32_t parent = nodes_[node].parent;
        if (parent == kNoParent)
            return;
        if (!dirty_[parent]) {
            dirty_[parent] = 1;
            dirty_list_.push_back(parent);
        }
        if (old == kTrue) {
            --count_a_[parent];
        } else if (old == kFalse) {
            --count_b_[parent];
        }
        if (value == kTrue) {
            ++count_a_[parent];
        } else if (value == kFalse) {
            ++count_b_[parent];
        }
        value = combinator_value(nodes_[parent], count_a_[parent], count_b_[parent]);
        node = parent;
    }
}

uint8_t Bound::evaluate() const {
    // the leaves' counts of the marked labels
    for (LabelId label : marked_) {
        const bool sure = mark_[label] == kSure;
        for (uint32_t i = leaves_begin_[label]; i < leaves_begin_[label + 1]; ++i) {
            const uint32_t leaf = leaves_[i];
            if (!dirty_[leaf]) {
                dirty_[leaf] = 1;
                dirty_list_.push_back(leaf);
            }
            ++(sure ? count_a_ : count_b_)[leaf];
        }
    }
    // their values, carried up to the root (set_value appends the ancestors it changes)
    const size_t touched = dirty_list_.size();
    for (size_t i = 0; i < touched; ++i) {
        const uint32_t leaf = dirty_list_[i];
        set_value(leaf, leaf_value(nodes_[leaf], count_a_[leaf], count_b_[leaf]));
    }
    const uint8_t result = value_.back();

    for (uint32_t node : dirty_list_) {
        value_[node] = base_value_[node];
        count_a_[node] = base_true_[node];
        count_b_[node] = base_false_[node];
        dirty_[node] = 0;
    }
    dirty_list_.clear();
    unmark();
    return result;
}

uint64_t Bound::units(const LabelId *present, size_t n) const {
    const uint64_t units = mark(present, n, nullptr, 0);
    unmark();
    return units;
}

bool Bound::eval(const LabelId *present, size_t n, uint64_t *units) const {
    const uint64_t u = mark(present, n, nullptr, 0);
    if (units)
        *units = u;
    if (constant_) {
        unmark();
        return *constant_;
    }
    const uint8_t value = evaluate();
    if (value == kUnknown)
        throw std::logic_error("predicate: no value without a label in doubt");
    return value == kTrue;
}

std::optional<bool> Bound::eval3(const LabelId *sure, size_t num_sure,
                                 const LabelId *maybe, size_t num_maybe,
                                 uint64_t *units) const {
    const uint64_t u = mark(sure, num_sure, maybe, num_maybe);
    if (units)
        *units = u;
    if (constant_) {
        unmark();
        return *constant_;
    }
    const uint8_t value = evaluate();
    if (value == kUnknown)
        return std::nullopt;
    return value == kTrue;
}


// ---------------------------------------------------------------- the answer

void Bound::write_json(uint32_t id, Json::Value *out, uint64_t *written,
                       const std::function<void()> &check) const {
    auto count = [&]() {
        if (check && *written % kStride == 0)
            check();
        ++*written;
    };
    count();
    const Node &node = nodes_[id];
    const char *key = to_string(node.op);
    if (is_leaf(node.op)) {
        Json::Value list(Json::arrayValue);
        for (uint32_t i = node.begin; i < node.end; ++i) {
            count();
            list.append(labels_[items_[i]].name);
        }
        if (node.op == Op::AT_LEAST) {
            Json::Value inner(Json::objectValue);
            inner["n"] = Json::Value(static_cast<Json::UInt>(node.n));
            inner["labels"] = std::move(list);
            (*out)[key] = std::move(inner);
        } else {
            (*out)[key] = std::move(list);
        }
    } else if (node.op == Op::NOT) {
        Json::Value operand(Json::objectValue);
        write_json(items_[node.begin], &operand, written, check);
        (*out)[key] = std::move(operand);
    } else {
        Json::Value operands(Json::arrayValue);
        for (uint32_t i = node.begin; i < node.end; ++i) {
            Json::Value operand(Json::objectValue);
            write_json(items_[i], &operand, written, check);
            operands.append(std::move(operand));
        }
        (*out)[key] = std::move(operands);
    }
}

Json::Value Bound::normal_form_json(const std::function<void()> &check) const {
    if (constant_)
        return Json::Value(*constant_);
    Json::Value out(Json::objectValue);
    uint64_t written = 0;
    write_json(nodes_.size() - 1, &out, &written, check);
    return out;
}

Json::Value Bound::unknown_labels_json(const std::function<void()> &check) const {
    Json::Value out(Json::arrayValue);
    for (size_t i = 0; i < unknown_.size(); ++i) {
        if (check && i % kStride == 0)
            check();
        out.append(unknown_[i]);
    }
    return out;
}

} // namespace predicate
} // namespace cli
} // namespace mtg

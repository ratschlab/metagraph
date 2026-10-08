#ifndef __METAGRAPH_CLI_PATTERN_PREDICATE_HPP__
#define __METAGRAPH_CLI_PATTERN_PREDICATE_HPP__

/**
 * The predicate language of POST /pattern (increment 5b, SPEC-pattern-search.md §19.3,
 * §19.4, §19.9; the owner's decisions P1, P3, P5, P7, P18-P20 of 2026-10-08): a logical
 * condition on the annotation columns of a graph context or of a supported path,
 *
 *   {"any": [n1, ...]}  {"all": [...]}  {"none": [...]}
 *   {"at_least": {"n": m, "labels": [...]}}  {"and": [p, ...]}  {"or": [p, ...]}  {"not": p}
 *
 * parsed from the request (Predicate::parse), bound to one index's columns (bind: unknown
 * names folded away into the normal form) and evaluated on the labels a context carries
 * (Bound::eval) or may carry (Bound::eval3, three-valued, for the supported-path search's
 * pruning). Labels are column names and nothing more (P6: no taxonomy; a cohort is a list).
 * Nothing here reads annotation: the selection pass and the supported-path search read the
 * rows and hand the labels they found to the evaluator.
 */

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include <json/json.h>

#include "graph/traversal/traversal_types.hpp"


namespace mtg {

namespace graph {
namespace pattern {
class Budget;
}
namespace traversal {
class LabelOracle;
}
}

namespace cli {
namespace predicate {

using graph::traversal::Column;
using graph::traversal::LabelId;
using graph::traversal::LabelRef;

// The operators of the language, in the order the capabilities list them
enum class Op : uint8_t { ANY, ALL, NONE, AT_LEAST, AND, OR, NOT };
const char* to_string(Op op);

// The deepest nesting of and, or and not (§19.3): a predicate with a combinator inside 64
// others is refused
constexpr unsigned kMaxNesting = 64;

/**
 * A predicate as the request wrote it, checked (§19.3) and not yet bound to an index. Its
 * nodes are stored children first, the root last; a node's items are the indices of the names
 * it lists (any, all, none, at_least, in the request's order) or of its operands (and, or,
 * not). It holds less than the request's JSON it was parsed from.
 */
struct Predicate {
    struct Node {
        Op op = Op::ANY;
        uint32_t n = 0;         // at_least: the threshold, 1 <= n <= its names
        uint32_t begin = 0;     // the node's items: items[begin, end)
        uint32_t end = 0;
    };

    // the distinct names of the request's lists, in the order of their first appearance
    std::vector<std::string> names;
    std::vector<Node> nodes;
    std::vector<uint32_t> items;
    // the names of all lists together, a name in two lists counted twice: what
    // caps.max_predicate_labels bounds
    uint64_t listed = 0;

    uint32_t root() const { return nodes.size() - 1; }

    /**
     * The predicate |json| of the request (|path| names it in the messages). Throws
     * PatternRefusal 400 "invalid_request" for the first rule of §19.3 it breaks, the message
     * naming the path (request.predicate.and[1].none[0]: expected a string ...), and once its
     * form is valid, 400 "predicate_too_large" when it lists more than |max_labels| names
     * (the message naming the count and the cap). The names are checked through views into
     * |json| and copied only after the cap: a body of millions of names costs no copy, and
     * the walk's index of them (about 64 bytes a distinct name, 4 a listed one) holds at most
     * |max_labels| of them: past the cap a name's form is checked and the name counted, no
     * more, so a name listed twice in one list after the cap is refused for the size, not
     * for the repetition. Not polled, as the rest of the route's parse of the request: its
     * time is linear in the body.
     */
    static Predicate parse(const Json::Value &json, uint64_t max_labels,
                           const std::string &path = "request.predicate");
};

// the column of an annotation label, or nullopt when no column has that name
using ColumnLookup = std::function<std::optional<Column>(const std::string &name)>;

struct Binding;

/**
 * A predicate bound to one index (§19.4): the names that are not columns of its annotation
 * are unknown, absent from every context, and folded away (any of unknowns false, none true,
 * all with an unknown false, at_least with fewer than n known false, constants through and,
 * or and not); what is left is the normal form, the predicate the answer echoes and the
 * selection evaluates. A normal form that is a constant is evaluated without reading any
 * annotation (P18).
 *
 * labels() is the permitted set of the selection's reads: the columns the normal form lists,
 * in the order of their first appearance in the request; a LabelId is an index into it (as a
 * LabelQuery built on it numbers them). The evaluators take such ids.
 *
 * Evaluation is incremental: each node's value on a context carrying none of the labels is
 * kept, and a call visits only the leaves its labels are listed in and their ancestors (at
 * most kMaxNesting + 1 levels each), then restores them. The scratch it uses lives here and
 * was allocated (and modelled in bytes()) when the predicate was bound: an evaluation
 * allocates nothing. Like LabelOracle, one Bound is used by one thread at a time.
 */
class Bound {
  public:
    /**
     * Binds |predicate| to the columns |find_column| knows (§19.4): resolves its distinct
     * names (the work time read through |budget| before every 4,096, and every 4,096 nodes and
     * items folded and built: a stop is "time", nothing returned), folds the unknown ones
     * away, and asks |admit| (when given; e.g. the request's memory account) for the bound
     * predicate's bytes before anything is copied ("max_memory" when refused). |budget| may
     * throw Aborted (Budget::set_abort). Nothing is charged as work: the names are bounded by
     * max_predicate_labels and each is one hash lookup. Its working arrays (about 20 bytes a
     * name and 8 a parsed node, freed on return) are not admitted: they hold less than
     * |predicate|, which the route holds before its account exists. (A static member, never
     * a free bind(): an unqualified call with a std::function argument would find std::bind.)
     */
    static Binding bind(const Predicate &predicate,
                        const ColumnLookup &find_column,
                        graph::pattern::Budget &budget,
                        const std::function<bool(uint64_t)> &admit = nullptr);

    // The same with |oracle|'s columns (LabelOracle::find_column; record headers are not
    // looked up: a header that is no column is unknown, P20)
    static Binding bind(const Predicate &predicate,
                        const graph::traversal::LabelOracle &oracle,
                        graph::pattern::Budget &budget,
                        const std::function<bool(uint64_t)> &admit = nullptr);

    // the distinct names of the request's lists, and of them the index's columns (also those
    // the folding dropped): the answer's predicate.names and predicate.known
    uint64_t num_names() const { return num_names_; }
    uint64_t num_known() const { return num_known_; }
    // the names that are not columns of the index (record headers among them, P20), in the
    // order of their first appearance in the request: predicate.unknown_labels
    const std::vector<std::string>& unknown_labels() const { return unknown_; }

    // the normal form is the constant false or true (no node, no label)
    std::optional<bool> constant() const { return constant_; }
    // the normal form holds on a context carrying none of its labels (§19.4, P5): none(A),
    // not(any(A)), or(none(A), any(B)), the constant true
    bool vacuous() const { return vacuous_; }
    // the normal form has neither none nor not: it cannot turn true when a label is taken
    // away (the supported-path search may prune a branch on which it is already false, §20.9)
    bool monotone() const { return monotone_; }

    // the permitted set: the columns of the normal form (LabelId = index). Empty for a
    // constant: then build no LabelQuery (it refuses an empty set) and evaluate nothing
    // (constant() is the answer; eval() throws for every id)
    const std::vector<LabelRef>& labels() const { return labels_; }
    // the nodes of the normal form (0 for a constant)
    uint64_t num_nodes() const { return nodes_.size(); }

    /**
     * The work units of deciding a context whose present labels are |present|[0, n) (§19.9):
     * 1, and 1 per leaf of the normal form each distinct present label is listed in. A label
     * listed twice in |present| counts once. Throws std::invalid_argument for an id outside
     * labels(). Known before the evaluation, which it bounds: the evaluation visits at most
     * kMaxNesting + 1 nodes per unit.
     */
    uint64_t units(const LabelId *present, size_t n) const;

    /**
     * Whether the normal form holds on a context whose labels among labels() are
     * |present|[0, n) (any order; a label listed twice counts once). |units|, when given,
     * receives units(present, n). Throws std::invalid_argument for an id outside labels().
     */
    bool eval(const LabelId *present, size_t n, uint64_t *units = nullptr) const;

    /**
     * Kleene's three-valued evaluation: the labels |sure| are present, the labels |maybe| may
     * be (an id in both is sure), every other label of labels() is absent. True or false when
     * the normal form takes that value whatever the maybes turn out to be, nullopt when it
     * cannot tell. Sound, not complete: a definite answer is never contradicted by a
     * completion of the maybes, but or(any(A), none(A)) with A maybe answers nullopt; a leaf
     * alone is decided exactly. |units| as eval's over the union of the two lists.
     */
    std::optional<bool> eval3(const LabelId *sure, size_t num_sure,
                              const LabelId *maybe, size_t num_maybe,
                              uint64_t *units = nullptr) const;

    /**
     * The answer's predicate.normal_form, in the request's syntax ({"any": [...]}, ...,
     * {"at_least": {"n": m, "labels": [...]}}), false or true for a constant; and its
     * unknown_labels. |check| (the answer's delivery check) is called before every 4,096
     * names or nodes written; it may throw.
     */
    Json::Value normal_form_json(const std::function<void()> &check = nullptr) const;
    Json::Value unknown_labels_json(const std::function<void()> &check = nullptr) const;
    // the bytes of the compact JSON text of the two (no indentation, no spaces), at most:
    // exact when every name is printable ASCII; a character that jsoncpp escapes counted as
    // its escape, a byte outside ASCII as 6 (its \u escape). What the answer's volume and
    // account are charged before the echo is built
    uint64_t text_bytes() const { return text_bytes_; }

    /**
     * The memory model of the bound predicate (§19.9), what bind() admits before it copies
     * anything: per label of labels() 192 + 2 x its length (its LabelRef here and in the
     * selection's LabelQuery), per unknown name 64 + its length, per node of the normal form
     * 64, and per name a leaf lists 8 (its entry in the leaf's list and in the label's list of
     * leaves).
     */
    uint64_t bytes() const { return bytes_; }

  private:
    // built by bind() only
    Bound() = default;

    struct Node {
        Op op = Op::ANY;
        uint32_t n = 0;          // at_least: the threshold
        uint32_t begin = 0;      // items_[begin, end): the LabelIds a leaf lists, as the
        uint32_t end = 0;        // request did, or a combinator's operands
        uint32_t parent = 0;     // kNoParent for the root
    };
    static constexpr uint32_t kNoParent = UINT32_MAX;

    uint64_t num_names_ = 0;
    uint64_t num_known_ = 0;
    std::vector<std::string> unknown_;
    std::optional<bool> constant_;
    bool vacuous_ = false;
    bool monotone_ = true;
    std::vector<LabelRef> labels_;

    // the normal form, children before their parent, the root last (empty for a constant)
    std::vector<Node> nodes_;
    std::vector<uint32_t> items_;
    // each label's leaves: leaves_[leaves_begin_[l], leaves_begin_[l + 1])
    std::vector<uint32_t> leaves_begin_;
    std::vector<uint32_t> leaves_;
    // per node, on a context carrying none of the labels: its value (kFalse, kTrue) and, for
    // a combinator, how many of its operands are true and false
    std::vector<uint8_t> base_value_;
    std::vector<uint32_t> base_true_;
    std::vector<uint32_t> base_false_;

    uint64_t text_bytes_ = 0;
    uint64_t bytes_ = 0;

    // The scratch of an evaluation, equal to the base between calls: per node its value
    // (kFalse, kTrue, kUnknown) and two counts (a leaf: its sure and maybe labels; a
    // combinator: its true and false operands), whether the call changed it, the nodes it
    // changed; per label its mark (sure, maybe) and the labels marked
    mutable std::vector<uint8_t> value_;
    mutable std::vector<uint32_t> count_a_;
    mutable std::vector<uint32_t> count_b_;
    mutable std::vector<uint8_t> dirty_;
    mutable std::vector<uint32_t> dirty_list_;
    mutable std::vector<uint8_t> mark_;
    mutable std::vector<LabelId> marked_;

    // marks the labels of the two lists (sure first), returns the units
    uint64_t mark(const LabelId *sure, size_t num_sure, const LabelId *maybe,
                  size_t num_maybe) const;
    // evaluates over the marked labels, restores the scratch and returns the root's value
    uint8_t evaluate() const;
    void unmark() const;
    void set_value(uint32_t node, uint8_t value) const;
    uint8_t leaf_value(const Node &node, uint32_t sure, uint32_t maybe) const;
    uint8_t combinator_value(const Node &node, uint32_t num_true, uint32_t num_false) const;
    void write_json(uint32_t node, Json::Value *out, uint64_t *written,
                    const std::function<void()> &check) const;
};

// The outcome of Bound::bind()
struct Binding {
    // the bound predicate, or empty when bind() stopped (|stop|)
    std::optional<Bound> bound;
    // why |bound| is empty: "time" (the work time passed before it was built: |budget|
    // recorded TIME) or "max_memory" (|admit| refused |bytes|); nullptr when bound
    const char *stop = nullptr;
    // the model's bytes of the bound predicate (Bound::bytes), known once every name was
    // resolved and the predicate folded; 0 before
    uint64_t bytes = 0;
    // what |admit| accepted: |bytes| once admitted (also when the time stopped the build
    // afterwards), else 0, and 0 when no |admit| was given (nothing was charged). The caller
    // releases it with the predicate
    uint64_t admitted = 0;
};

} // namespace predicate
} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_CLI_PATTERN_PREDICATE_HPP__

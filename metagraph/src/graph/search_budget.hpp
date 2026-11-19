#ifndef __SEARCH_BUDGET_HPP__
#define __SEARCH_BUDGET_HPP__

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>


namespace mtg {
namespace graph {

/**
 * A /search answer priced past the server's bound (--search-max-answer-memory-mb;
 * SPEC-multigraph-envelope.md section 3, the /search answer budget): nothing partial is sent
 * (the server: 400 search_answer_too_large). The message is the same whichever account passed
 * the bound; |account| ("values" or "answer") is the one that threw. The 400 names the
 * request's account, SearchBudget::account(), a function of the request (lane T7-P1), never
 * of the order the work charged in.
 */
class SearchAnswerTooLarge : public std::runtime_error {
  public:
    SearchAnswerTooLarge(uint64_t max_bytes, const char *account)
          : std::runtime_error(message(max_bytes >> 20)),
            max_bytes_(max_bytes), account_(account) {}

    static std::string message(uint64_t max_mb) {
        return "search: the answer would hold more than --search-max-answer-memory-mb "
               + std::to_string(max_mb) + " MiB: nothing partial is sent; ask for fewer "
               "top_labels, a larger discovery_fraction, fewer or shorter FASTA records, or "
               "no per-position values (with_signature, query_counts, query_coords)";
    }

    uint64_t max_bytes() const { return max_bytes_; }
    uint64_t max_mb() const { return max_bytes_ >> 20; }
    const char* account() const { return account_; }

  private:
    uint64_t max_bytes_;
    const char *account_;
};

/**
 * A unit of a /search request (one FASTA record on one pair) refused alone, with the request's
 * field `partial_refusal` (SPEC-multigraph-envelope.md section 3, the /search answer budget):
 * its values account passed the bound. Not a SearchAnswerTooLarge: the request goes on, and the
 * search that threw it (cli::execute_query, cli::query_sequence) answers the unit's refused stub
 * (no row, its cut marked). |account|: "values".
 */
class SearchUnitRefused : public std::runtime_error {
  public:
    explicit SearchUnitRefused(const char *account)
          : std::runtime_error(std::string("search: a record's ") + account
                               + " account on a pair passed the bound: the unit is refused alone"),
            account_(account) {}

    const char* account() const { return account_; }

  private:
    const char *account_;
};

// |a| + |b|, UINT64_MAX where the sum would pass the range
inline uint64_t saturating_add(uint64_t a, uint64_t b) {
    return a > std::numeric_limits<uint64_t>::max() - b ? std::numeric_limits<uint64_t>::max()
                                                        : a + b;
}

// |bytes| in MiB, rounded up
inline uint64_t mib_rounded_up(uint64_t bytes) {
    return (bytes >> 20) + ((bytes & ((uint64_t(1) << 20) - 1)) != 0);
}

class SearchBudget;

/**
 * One account of a /search request (SearchBudget): the request's answer (its JSON tree and
 * text; SearchBudget::answer()), or the values of one unit, a FASTA record on a pair (its rows
 * and their per-position values; lane T7-P1). Charged BEFORE what it prices is allocated, from
 * any thread; what it holds (held()) only grows, so that whether a charge passes the bound never
 * depends on the order the threads charge in. An account of one task alone (a unit's) may lower
 * its scratch where that task frees some (lower_scratch): its later charges use those bytes
 * first, so that held() is the most the task held at once, which depends on that one task's
 * order only. A model of the heap (the deployed build's sizes), never a measurement.
 *
 * A charge that would take it past the bound:
 * - a unit's values account, with `partial_refusal`: marks the unit refused (sticky: every
 *   further charge of it throws again) and throws SearchUnitRefused; the request goes on;
 * - otherwise (a unit's values without the field, or the answer account): marks the request
 *   (SearchBudget::note_passed) and throws SearchAnswerTooLarge; the request is refused 400
 *   (SearchBudget::passed), and the accounts the request's state closes (SearchBudget::
 *   check_open) throw at their next charge, so that the work still running stops.
 */
class SearchAccount {
  public:
    // An account of the request: |where|, what it is of, for the log; |threads|, how many
    // tasks of its work run at once, for charge_scratch
    SearchAccount(SearchBudget *budget, const char *name, std::string where = "",
                  size_t threads = 1)
          : budget_(budget), name_(name), where_(std::move(where)),
            threads_(std::max<size_t>(threads, 1)) {}
    // The values account of a unit: FASTA record |record| (its index in the body) on |pair|
    // ("pair (<name>, <graph>, <annotation>)", or "the graph" on a single-graph server), the
    // pair's place in the answer |pair_index| (the order of its `cuts` entries), one task (its
    // scratch charged once)
    SearchAccount(SearchBudget *budget, size_t record, const std::string &pair,
                  size_t pair_index = 0)
          : budget_(budget), name_("values"),
            where_("record " + std::to_string(record) + " on " + pair),
            threads_(1), record_(record), pair_index_(pair_index) {}
    SearchAccount(const SearchAccount &) = delete;
    SearchAccount& operator=(const SearchAccount &) = delete;

    inline void charge(uint64_t bytes);
    // Charges the scratch one task of the account's work holds while it runs and frees before
    // it ends (a record's store of its rows): the account holds |threads| times the most any
    // task's scratch took, raised (charge) before the scratch is allocated, which covers the
    // tasks running at once without charging a freed scratch again for every task; the total
    // it adds is the same in every order the tasks run in (lane T6-W2)
    inline void charge_scratch(uint64_t bytes);
    // An account of one task (threads 1): its scratch lowered to |bytes| where it holds more,
    // the difference kept as freed, which the task's later charges use before they raise held()
    // (and a scratch that grows again is charged again). Lane T7-G4E, the Opus review of
    // cf77380a0, MA-1: a column-label search's table frees its old buckets when it doubles,
    // and the unit's values and rows charged after it use them, as the heap does
    inline void lower_scratch(uint64_t bytes);

    uint64_t held() const { return held_.load(std::memory_order_relaxed); }
    // the most one task's scratch took (charge_scratch), lowered where a single task freed some
    // (lower_scratch), and how many run at once
    uint64_t scratch() const { return scratch_.load(std::memory_order_relaxed); }
    // what its one task freed and its charges have not used again (lower_scratch)
    uint64_t freed() const { return freed_.load(std::memory_order_relaxed); }
    size_t threads() const { return threads_; }
    const char* name() const { return name_; }
    const std::string& where() const { return where_; }
    SearchBudget* budget() const { return budget_; }
    // a unit's account (its record's index), and whether the unit was refused
    bool unit() const { return record_ != kNoRecord; }
    size_t record() const { return record_; }
    size_t pair_index() const { return pair_index_; }
    bool refused() const { return refused_.load(std::memory_order_relaxed); }

    static constexpr size_t kNoRecord = std::numeric_limits<size_t>::max();

  private:
    SearchBudget *budget_;
    const char *name_;
    std::string where_;
    size_t threads_;
    size_t record_ = kNoRecord;
    size_t pair_index_ = 0;
    std::atomic<uint64_t> held_ { 0 };
    std::atomic<uint64_t> scratch_ { 0 };
    std::atomic<uint64_t> freed_ { 0 };
    std::atomic<bool> refused_ { false };
};

/**
 * The accounts of one /search request, each bounded by |max_bytes| (the server's
 * --search-max-answer-memory-mb, in bytes): the answer's (answer()), and a values account per
 * unit, a FASTA record on a pair (SearchAccount(&budget, record, pair), made by the unit's
 * task). |partial_refusal|: the request's field (SPEC-multigraph-envelope.md section 3): a unit
 * whose values, or whose own object, would pass the bound is refused alone and reported in its
 * cut entry; without it any such unit refuses the whole request.
 *
 * The request's state (lane T7-P1), a function of the request, never of the order its work
 * charges in (the design note's section 4):
 * - OPEN;
 * - VALUES: refused, without the field, by a unit's values; every charge throws (the work
 *   stops early, as before the field existed);
 * - ANSWER: refused, with the field, by the request's answer account; every charge throws;
 * - ANSWER_PENDING: refused, without the field, on the answer side (a unit's object alone, or
 *   the request's answer account): no object is built or charged any more, and the units not
 *   yet finished run their searches, charging their own values accounts, so that whether some
 *   unit's values would pass is learned in every order (then VALUES).
 * The 400 names `account()`: "values" in VALUES, else "answer".
 */
class SearchBudget {
  public:
    explicit SearchBudget(uint64_t max_bytes, bool partial_refusal = false)
          : max_bytes_(max_bytes), partial_refusal_(partial_refusal), answer_(this, "answer") {}
    SearchBudget(const SearchBudget &) = delete;
    SearchBudget& operator=(const SearchBudget &) = delete;

    uint64_t max_bytes() const { return max_bytes_; }
    bool partial_refusal() const { return partial_refusal_; }
    // the request is refused (400 search_answer_too_large), whatever else its work raised
    bool passed() const { return state_.load(std::memory_order_acquire) != OPEN; }
    // the account the 400 names, once the request's work has ended: "values" or "answer"
    const char* account() const {
        return state_.load(std::memory_order_acquire) == VALUES ? "values" : "answer";
    }
    SearchAccount& answer() { return answer_; }

    // throws SearchAnswerTooLarge if the request's state closes |account| (every account once
    // the request is refused, but the units' values accounts while the account is pending)
    void check_open(const SearchAccount &account) const {
        const int state = state_.load(std::memory_order_acquire);
        if (state == OPEN || (state == ANSWER_PENDING && &account != &answer_))
            return;
        throw SearchAnswerTooLarge(max_bytes_, account.name());
    }

    // marks the request refused: |account|, holding |held| bytes, would have passed the bound
    // with a charge of |charge| bytes (the first one of each kind is kept for the log, read
    // once the request's work has ended)
    void note_passed(const SearchAccount &account, uint64_t held, uint64_t charge) {
        std::string where = std::string("the ") + account.name() + " account";
        if (!account.where().empty())
            where += " of " + account.where();
        refuse(&account != &answer_, std::move(where), held, charge);
    }
    // without the field: the object of the unit |where| (one record on one pair) alone, priced
    // |price| bytes, would pass the bound (the answer side; the request is refused)
    void note_object_passed(const std::string &where, uint64_t price) {
        refuse(false, "the answer account of " + where, 0, price);
    }

    // a unit refused alone (with the field), for the log: record |record| on its pair
    // (|where|), its |account| ("values" or "answer"), what it held before the charge refused
    // and that charge (for "answer": 0 and the unit's object's price), and the pair's place in
    // the answer (|pair|: the order of the record's `cuts` entries)
    struct UnitRefusal {
        size_t record;
        std::string where;
        const char *account;
        uint64_t held;
        uint64_t charge;
        size_t pair = 0;
    };
    void note_unit_refused(UnitRefusal refusal) {
        std::lock_guard<std::mutex> lock(mutex_);
        unit_refusals_.push_back(std::move(refusal));
    }
    // the units refused alone, in the answer's order: by record, then by the pair's place in
    // the answer (as the record's `cuts` entries; read once the request's work has ended)
    std::vector<UnitRefusal> unit_refusals() const {
        std::lock_guard<std::mutex> lock(mutex_);
        std::vector<UnitRefusal> refusals = unit_refusals_;
        std::sort(refusals.begin(), refusals.end(), [](const auto &a, const auto &b) {
            return std::make_tuple(a.record, a.pair, std::string_view(a.where),
                                   std::string_view(a.account))
                    < std::make_tuple(b.record, b.pair, std::string_view(b.where),
                                      std::string_view(b.account));
        });
        return refusals;
    }

    // the account the 400 names (account()): its first refusal, for the log
    const std::string& first_passed() const { return first(account()).where; }
    // what that account held before the charge refused, and that charge
    uint64_t first_held() const { return first(account()).held; }
    uint64_t first_charge() const { return first(account()).charge; }

    // the most a unit's values account held (noted as each unit's work ends), for the log
    void note_values(uint64_t bytes) {
        uint64_t max = values_max_.load(std::memory_order_relaxed);
        while (bytes > max
                && !values_max_.compare_exchange_weak(max, bytes, std::memory_order_relaxed)) {}
    }
    uint64_t values_max() const { return values_max_.load(std::memory_order_relaxed); }
    // the units whose object (or refused stub) was built, for the log
    void note_unit_built() { units_.fetch_add(1, std::memory_order_relaxed); }
    size_t units() const { return units_.load(std::memory_order_relaxed); }

  private:
    enum State : int { OPEN, VALUES, ANSWER, ANSWER_PENDING };
    struct First {
        bool noted = false;
        std::string where;
        uint64_t held = 0;
        uint64_t charge = 0;
    };

    void refuse(bool values, std::string where, uint64_t held, uint64_t charge) {
        std::lock_guard<std::mutex> lock(mutex_);
        First &first = values ? first_values_ : first_answer_;
        if (!first.noted)
            first = { true, std::move(where), held, charge };
        const int state = state_.load(std::memory_order_relaxed);
        if (values) {
            if (state != ANSWER)
                state_.store(VALUES, std::memory_order_release);
        } else if (state == OPEN) {
            state_.store(partial_refusal_ ? ANSWER : ANSWER_PENDING, std::memory_order_release);
        }
    }
    const First& first(const char *account) const {
        return account[0] == 'v' ? first_values_ : first_answer_;
    }

    const uint64_t max_bytes_;
    const bool partial_refusal_;
    std::atomic<int> state_ { OPEN };
    mutable std::mutex mutex_;
    First first_values_;
    First first_answer_;
    std::vector<UnitRefusal> unit_refusals_;
    std::atomic<uint64_t> values_max_ { 0 };
    std::atomic<size_t> units_ { 0 };
    SearchAccount answer_;
};

inline void SearchAccount::charge(uint64_t bytes) {
    budget_->check_open(*this);
    if (refused())
        throw SearchUnitRefused(name_);
    // what the account's one task freed first (lower_scratch; zero for an account of more
    // tasks): held already, it is used again before held() grows
    if (const uint64_t freed = freed_.load(std::memory_order_relaxed)) {
        const uint64_t used = std::min(freed, bytes);
        freed_.store(freed - used, std::memory_order_relaxed);
        bytes -= used;
        if (!bytes)
            return;
    }
    const uint64_t max = budget_->max_bytes();
    uint64_t held = held_.load(std::memory_order_relaxed);
    do {
        if (bytes > max || held > max - bytes) {
            if (unit() && budget_->partial_refusal()) {
                // the unit alone: the request goes on
                if (!refused_.exchange(true, std::memory_order_relaxed))
                    budget_->note_unit_refused({ record_, where_, name_, held, bytes,
                                                 pair_index_ });
                throw SearchUnitRefused(name_);
            }
            budget_->note_passed(*this, held, bytes);
            throw SearchAnswerTooLarge(max, name_);
        }
    } while (!held_.compare_exchange_weak(held, held + bytes, std::memory_order_relaxed));
}

// |account|.charge(|bytes|) where an account is given (null: none, the CLI's)
inline void charge_search(SearchAccount *account, uint64_t bytes) {
    if (account && bytes)
        account->charge(bytes);
}

inline void SearchAccount::charge_scratch(uint64_t bytes) {
    uint64_t most = scratch_.load(std::memory_order_relaxed);
    while (bytes > most) {
        if (scratch_.compare_exchange_weak(most, bytes, std::memory_order_relaxed)) {
            // the raise for each task that can run at once (a product past the range is past
            // every bound)
            const uint64_t raise = bytes - most;
            charge(raise > std::numeric_limits<uint64_t>::max() / threads_
                       ? std::numeric_limits<uint64_t>::max()
                       : raise * threads_);
            return;
        }
    }
}

inline void SearchAccount::lower_scratch(uint64_t bytes) {
    assert(threads_ == 1);
    const uint64_t most = scratch_.load(std::memory_order_relaxed);
    if (bytes < most) {
        scratch_.store(bytes, std::memory_order_relaxed);
        freed_.fetch_add(most - bytes, std::memory_order_relaxed);
    }
}

// |account|.charge_scratch(|bytes|) where an account is given (null: none, the CLI's)
inline void charge_search_scratch(SearchAccount *account, uint64_t bytes) {
    if (account && bytes)
        account->charge_scratch(bytes);
}

} // namespace graph
} // namespace mtg

#endif // __SEARCH_BUDGET_HPP__

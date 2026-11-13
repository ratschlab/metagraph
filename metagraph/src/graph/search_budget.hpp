#ifndef __SEARCH_BUDGET_HPP__
#define __SEARCH_BUDGET_HPP__

#include <atomic>
#include <cstdint>
#include <stdexcept>
#include <string>


namespace mtg {
namespace graph {

/**
 * A /search answer priced past the server's bound (--search-max-answer-memory-mb;
 * SPEC-multigraph-envelope.md section 3, the /search answer budget): nothing partial is sent
 * (the server: 400 search_answer_too_large). The message is the same whichever account passed
 * the bound, so that the answer never depends on the order the work charged in; |account|
 * ("values" or "answer") is for the server's log only.
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

class SearchBudget;

/**
 * One account of a /search request (SearchBudget): the values a pair's rows hold (one per
 * pair) or the request's answer (its JSON tree and text). Charged BEFORE what it prices is
 * allocated, from any thread; it only grows (nothing is released while the request runs), so
 * that whether a charge passes the bound never depends on the order the threads charge in. A
 * charge that would take it past the bound throws SearchAnswerTooLarge and marks the request
 * (SearchBudget::passed); from then on every charge of the request throws, so that the work
 * still running stops at its next charge. A model of the heap (the deployed build's sizes),
 * never a measurement.
 */
class SearchAccount {
  public:
    // |where|: what the account is of, for the log (a values account's pair)
    SearchAccount(SearchBudget *budget, const char *name, std::string where = "")
          : budget_(budget), name_(name), where_(std::move(where)) {}
    SearchAccount(const SearchAccount &) = delete;
    SearchAccount& operator=(const SearchAccount &) = delete;

    inline void charge(uint64_t bytes);

    uint64_t held() const { return held_.load(std::memory_order_relaxed); }
    const char* name() const { return name_; }
    const std::string& where() const { return where_; }
    SearchBudget* budget() const { return budget_; }

  private:
    SearchBudget *budget_;
    const char *name_;
    std::string where_;
    std::atomic<uint64_t> held_ { 0 };
};

/**
 * The accounts of one /search request, each bounded by |max_bytes| (the server's
 * --search-max-answer-memory-mb, in bytes): the answer's (answer()), and a values account per
 * pair the pair's work makes (SearchAccount(&budget, "values")). |passed| once any of them
 * would have passed the bound; the server then answers the request 400
 * search_answer_too_large, whatever else its work raised.
 */
class SearchBudget {
  public:
    explicit SearchBudget(uint64_t max_bytes) : max_bytes_(max_bytes), answer_(this, "answer") {}
    SearchBudget(const SearchBudget &) = delete;
    SearchBudget& operator=(const SearchBudget &) = delete;

    uint64_t max_bytes() const { return max_bytes_; }
    bool passed() const { return passed_.load(std::memory_order_acquire); }
    SearchAccount& answer() { return answer_; }

    // marks the request: |account| would have passed the bound (the first one is kept for
    // the log, read once the request's work has ended)
    void note_passed(const SearchAccount &account) {
        bool expected = false;
        if (passed_.compare_exchange_strong(expected, true, std::memory_order_acq_rel)) {
            first_passed_ = std::string("the ") + account.name() + " account";
            if (!account.where().empty())
                first_passed_ += " of " + account.where();
        }
    }
    const std::string& first_passed() const { return first_passed_; }

    // the most a pair's values account held (noted as each pair's work ends), for the log
    void note_values(uint64_t bytes) {
        uint64_t max = values_max_.load(std::memory_order_relaxed);
        while (bytes > max
                && !values_max_.compare_exchange_weak(max, bytes, std::memory_order_relaxed)) {}
    }
    uint64_t values_max() const { return values_max_.load(std::memory_order_relaxed); }

  private:
    const uint64_t max_bytes_;
    std::atomic<bool> passed_ { false };
    std::string first_passed_;
    std::atomic<uint64_t> values_max_ { 0 };
    SearchAccount answer_;
};

inline void SearchAccount::charge(uint64_t bytes) {
    if (budget_->passed())
        throw SearchAnswerTooLarge(budget_->max_bytes(), name_);
    const uint64_t max = budget_->max_bytes();
    uint64_t held = held_.load(std::memory_order_relaxed);
    do {
        if (bytes > max || held > max - bytes) {
            budget_->note_passed(*this);
            throw SearchAnswerTooLarge(max, name_);
        }
    } while (!held_.compare_exchange_weak(held, held + bytes, std::memory_order_relaxed));
}

// |account|.charge(|bytes|) where an account is given (null: none, the CLI's)
inline void charge_search(SearchAccount *account, uint64_t bytes) {
    if (account && bytes)
        account->charge(bytes);
}

} // namespace graph
} // namespace mtg

#endif // __SEARCH_BUDGET_HPP__

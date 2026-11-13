#include "boss_pattern.hpp"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <stdexcept>
#include <string>


namespace mtg {
namespace graph {
namespace boss {

typedef BOSS::TAlphabet TAlphabet;
typedef BOSS::edge_index edge_index;

SymbolSet all_symbols(const BOSS &boss) {
    return (SymbolSet(1) << boss.alph_size) - 2;
}

#if ! _PROTEIN_GRAPH
namespace {

const std::pair<char, const char*> kIupac[] = {
    { 'A', "A" }, { 'C', "C" }, { 'G', "G" }, { 'T', "T" },
    { 'R', "AG" }, { 'Y', "CT" }, { 'S', "CG" }, { 'W', "AT" }, { 'K', "GT" }, { 'M', "AC" },
    { 'B', "CGT" }, { 'D', "AGT" }, { 'H', "ACT" }, { 'V', "ACG" }, { 'N', "ACGT" },
};

} // namespace
#endif

Pattern parse_pattern(const BOSS &boss, std::string_view sequence) {
    Pattern pattern;
    for (char c : sequence) {
        SymbolSet set = 0;
#if _PROTEIN_GRAPH
        c = std::toupper(static_cast<unsigned char>(c));
        if (c == 'X') {
            set = all_symbols(boss);
        } else if (size_t code = boss.alphabet.find(c); code && code != std::string::npos) {
            set = SymbolSet(1) << code;
        }
#else
#if _DNA_CASE_SENSITIVE_GRAPH
        // the alphabet holds both cases: a code stands for bases of its own case
        const bool lower = std::islower(static_cast<unsigned char>(c));
#else
        const bool lower = false;
#endif
        c = std::toupper(static_cast<unsigned char>(c));
        for (const auto &[code, bases] : kIupac) {
            for (const char *b = bases; code == c && *b; ++b) {
                set |= SymbolSet(1) << boss.encode(lower ? std::tolower(*b) : *b);
            }
        }
#endif
        if (!set)
            throw std::invalid_argument(std::string("invalid pattern symbol '") + c + "'");

        pattern.push_back(set);
    }
    return pattern;
}

#if ! _PROTEIN_GRAPH
Pattern reverse_complement_pattern(Pattern pattern) {
    std::reverse(pattern.begin(), pattern.end());
    for (SymbolSet &set : pattern) {
        SymbolSet complement = 0;
        for (TAlphabet s = 1; set >> s; ++s) {
            if ((set >> s) & 1)
                complement |= SymbolSet(1) << kmer::KmerExtractorBOSS::complement(s);
        }
        set = complement;
    }
    return pattern;
}
#endif

bool SearchBudget::step() {
    if (stop_)
        return false;

    if (++steps_ > max_steps_) {
        stop_ = "max_steps";
    } else if (!(steps_ % 1024) && std::chrono::steady_clock::now() > deadline_) {
        stop_ = "time_budget";
    }
    return !stop_;
}

uint64_t count_edges(const BOSS &boss, edge_index rl, edge_index ru, SymbolSet symbols) {
    if (symbols == all_symbols(boss)) {
        return ru - rl + 1 - (boss.rank_W(ru, BOSS::kSentinelCode)
                                - boss.rank_W(rl - 1, BOSS::kSentinelCode));
    }

    uint64_t count = 0;
    for (TAlphabet s = 1; s < boss.alph_size; ++s) {
        if ((symbols >> s) & 1) {
            count += boss.rank_W(ru, s) - boss.rank_W(rl - 1, s)
                + boss.rank_W(ru, s + boss.alph_size) - boss.rank_W(rl - 1, s + boss.alph_size);
        }
    }
    return count;
}

bool call_edges(const BOSS &boss, edge_index rl, edge_index ru, SymbolSet symbols,
                const std::function<void(edge_index)> &callback,
                SearchBudget *budget) {
    for (TAlphabet s = 1; s < boss.alph_size; ++s) {
        if (!((symbols >> s) & 1))
            continue;

        for (edge_index i = rl; i <= ru; ++i) {
            i = boss.succ_W(i, s, s + boss.alph_size).first;
            if (i > ru)
                break;

            if (!budget->step())
                return false;

            callback(i);
        }
    }
    return true;
}

bool call_pattern_ranges(const BOSS &boss,
                         const Pattern &pattern,
                         const std::function<void(edge_index, edge_index,
                                                  size_t, SymbolSet)> &callback,
                         SearchBudget *budget) {
    assert(pattern.size() && pattern.size() <= boss.get_k() + 1);

    // The nodes of [rl, ru] end with |length| symbols: the pattern's first ones, or
    // the whole pattern and a part of its right flank.
    // A full node is never extended: its only edge to a successor can be marked.
    std::function<bool(edge_index, edge_index, size_t)> search
            = [&](edge_index rl, edge_index ru, size_t length) {
        if (length + 1 == pattern.size()) {
            // the edge symbol is the pattern's last one
            callback(rl, ru, boss.get_k() - length, pattern.back());
        } else if (length >= pattern.size()) {
            callback(rl, ru, boss.get_k() - length, all_symbols(boss));
        }

        if (length == boss.get_k())
            return true;

        SymbolSet next = length < pattern.size() ? pattern[length] : all_symbols(boss);
        for (TAlphabet s = 1; s < boss.alph_size; ++s) {
            if (!((next >> s) & 1))
                continue;

            if (!budget->step())
                return false;

            edge_index l = rl;
            edge_index u = ru;
            if (boss.tighten_range(&l, &u, s) && !search(l, u, length + 1))
                return false;
        }
        return true;
    };

    return search(1, boss.num_edges(), 0) && !budget->stop();
}

} // namespace boss
} // namespace graph
} // namespace mtg

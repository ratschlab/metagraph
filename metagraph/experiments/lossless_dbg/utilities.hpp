//
//  utils.hpp
//  lossless_dbg
//
//  Created by Jan Studený on 11/03/2019.
//  Copyright © 2019 Jan Studený. All rights reserved.
//

#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include "boss_construct.hpp"
#include "boss.hpp"
#include "kmer_extractor.hpp"
#include "sequence_io.hpp"
#include "unix_tools.hpp"
#include "logger.hpp"

using namespace std;
using ll = long long;

#define x first
#define y second
#define all(x) begin(x), end(x)

namespace fs = std::filesystem;

#define local_file(filename) (fs::path(__FILE__).parent_path() / (filename))

// template<typename POD>
// std::ostream& serialize(std::ostream& os, std::vector<POD> const& v);
//
// template<typename POD>
// std::istream& deserialize(std::istream& is, std::vector<POD>& v);

#ifndef RELEASE_MODE
#include <cassert>
#define alt_assert(condition) \
    { \
        if (!(condition)) { \
            std::cerr << "Assertion failed at " << __FILE__ << ":" << __LINE__; \
            std::cerr << " inside " << __FUNCTION__ << std::endl; \
            std::cerr << "Condition: " << #condition; \
            abort(); \
        } \
    }
#else
#define alt_assert(condition) (condition)
#endif

struct d_t {
    d_t() {}
    d_t(d_t &&) = default;
    // d_t(d_t &&d) : content(std::move(d.content)){}

    ~d_t() {
        if (!content.str().empty())
            //using mg::common::logger;
            mg::common::logger->debug(content.str());
        else {
            std::runtime_error("Something wrong");
        }
    }

    template <typename T>
    d_t &&operator,(const T &first) {
        content << ' ' << x;
        return std::move(*this);
    }
    stringstream content;
};

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#define PRINT_VAR(args...) \
    { d_t(), "|", __FILENAME__,__LINE__,__FUNCTION__, "|", #args, ":", args; }

template <typename POD>
inline std::istream &deserialize(std::istream &is, vector<POD> &v) {
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
                  "Can only deserialize POD types with this function");

    decltype(v.size()) size;
    is.read(reinterpret_cast<char *>(&size), sizeof(size));
    v.resize(size);
    is.read(reinterpret_cast<char *>(v.data()), v.size() * sizeof(POD));
    return is;
}

template <typename POD>
inline std::ostream &serialize(std::ostream &os, const vector<POD> &v) {
    // this only works on built in data types (PODs)
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
                  "Can only serialize POD types with this function");

    auto size = v.size();
    os.write(reinterpret_cast<char const *>(&size), sizeof(size));
    os.write(reinterpret_cast<char const *>(v.data()), v.size() * sizeof(POD));
    return os;
}

inline string decode_from_input(const vector<string> &input,
                                const pair<string, int64_t> &path_id,
                                int64_t kmer_length) {
    int64_t current_relative_id = 0;
    for (auto &sequence : input) {
        if (sequence.substr(0, kmer_length) == path_id.first) {
            if (current_relative_id == path_id.second) {
                return sequence;
            } else {
                current_relative_id++;
            }
        }
    }
    return "";
}

inline void save_string(const string &to_save, const string &filename) {
    ofstream myfile;
    myfile.open(filename);
    if (not myfile) {
        throw std::runtime_error("can't open "s + filename + "for writing");
    }
    myfile << to_save;
    myfile.close();
}

inline void transform_to_fasta(const string &filename, const vector<string> &reads) {
    ofstream myfile;
    myfile.open(filename);
    if (not myfile) {
        throw std::runtime_error("can't open "s + filename + "for reading");
    }
    for (auto &read : reads) {
        myfile << ">" << endl;
        myfile << read << endl;
    }
    myfile.close();
}


inline void write_reads_to_fasta(const vector<string> &reads, const string &filename) {
    transform_to_fasta(filename, reads);
}

inline void reduce_maps(std::map<int, int> &output, std::map<int, int> &input) {
    for (auto &X : input) {
        output[X.first] += X.second;
    }
}


// inline them
template <typename T>
inline void append_vectors(std::vector<T> &output, std::vector<T> &input) {
    output.insert(output.end(), begin(input), end(input));
}
// inline them
#pragma omp declare reduction(append : \
    std::vector<string> : \
    append_vectors(omp_out, omp_in))

#pragma omp declare reduction(append : \
    std::vector<char> : \
    append_vectors(omp_out, omp_in))

#pragma omp declare reduction(append : \
    std::vector<int64_t> : \
    append_vectors(omp_out, omp_in))

#pragma omp declare reduction(append : \
    std::vector<uint64_t> : \
    append_vectors(omp_out, omp_in))

#pragma omp declare reduction(append : \
    std::vector<bool> : \
    append_vectors(omp_out, omp_in))

const static int64_t delimiter_encoded = 6;
inline int8_t encode(char c) {
    if (c == '#')
        return delimiter_encoded; // alphabet_decoder.alph_size;
    if (c == '$')
        return 0;
    return KmerExtractorBOSS::encode(c);
}
inline char decode(int8_t c) {
    if (c == delimiter_encoded)
        return '#';
    if (c == 0)
        return '$';
    return KmerExtractorBOSS::decode(c);
}

inline string &
clamp_alphabet(string &text, const string &alphabet = "$ACGTN", char replacement = 'N') {
    for (auto &c : text) {
        if (alphabet.find(c) == string::npos) {
            c = replacement;
        }
    }
    return text;
}

// TODO: change occurences of read_reads to streaming mode (with callback function)
inline vector<string> read_reads_from_fasta(const string &filename) {
    vector<string> result;
    read_fasta_file_critical(filename, [&](kseq_t *read) {
        string read_seq = read->seq.s;
        clamp_alphabet(read_seq);
        result.push_back(read_seq);
    });
    return result;
}

template <typename... Args>
inline void doPrint(Args &&... args) {
    mg::common::logger->debug("{}", make_tuple(args...));
}

/*
    - Multiple inheritance forces ambiguity of member names.
    - SFINAE is used to make aliases to member names.
    - Expression SFINAE is used in just one generic has_member that can accept
      any alias we pass it.
*/

// Variadic to force ambiguity of class members.  C++11 and up.
template <typename... Args>
struct ambiguate : public Args... {};

// Non-variadic version of the line above.
// template <typename A, typename B> struct ambiguate : public A, public B {};

template <typename A, typename = void>
struct got_type : std::false_type {};

template <typename A>
struct got_type<A> : std::true_type {
    typedef A type;
};

template <typename T, T>
struct sig_check : std::true_type {};

template <typename Alias, typename AmbiguitySeed>
struct has_member {
    template <typename C>
    static char((&f(decltype(&C::value))))[1];
    template <typename C>
    static char((&f(...)))[2];

    // Make sure the member name is consistently spelled the same.
    static_assert((sizeof(f<AmbiguitySeed>(0)) == 1),
                  "Member name specified in AmbiguitySeed is different from member name "
                  "specified in Alias, or wrong Alias/AmbiguitySeed has been specified.");

    static bool const value = sizeof(f<Alias>(0)) == 2;
};

inline int64_t get_used_memory() {
    struct rusage my_rusage {};
    getrusage(RUSAGE_SELF, &my_rusage);
    return my_rusage.ru_maxrss;
}

// Check for function with given name, any signature.
// Check for any member with given name, whether var, func, class, union, enum.
#define CREATE_MEMBER_CHECK(member) \
\
    template <typename T, typename = std::true_type> \
    struct Alias_##member; \
\
    template <typename T> \
    struct Alias_##member<T, std::integral_constant<bool, got_type<decltype(&T::member)>::value>> { \
        static const decltype(&T::member) value; \
    }; \
\
    struct AmbiguitySeed_##member { \
        char member; \
    }; \
\
    template <typename T> \
    struct has_member_##member { \
        static const bool value \
                = has_member<Alias_##member<ambiguate<T, AmbiguitySeed_##member>>, \
                             Alias_##member<AmbiguitySeed_##member>>::value; \
    }


// TODO: move to tests
#undef protected
template <typename Test, typename Reference>
class IdentityComparator : public Reference {
  public:
    template <typename... Args>
    IdentityComparator(const Args &... args) : Reference(args...), t(args...) {}

    virtual ~IdentityComparator() {}

    CREATE_MEMBER_CHECK(offset);
    template <typename... Args>
    int64_t offset(Args... args) const {
        if constexpr (has_member_offset<Reference>::value) {
            auto target = Reference::offset(args...);
            auto value = t.offset(args...);
            if (target != value) {
                doPrint(target, value, args...);
            }
            assert(target == value);
            return target;
        } else {
            throw std::runtime_error("Can't call function offset");
        }
    }


    CREATE_MEMBER_CHECK(size);
    template <typename... Args>
    int64_t size(Args... args) const {
        if constexpr (has_member_size<Reference>::value) {
            auto target = Reference::size(args...);
            auto value = t.size(args...);
            if (target != value) {
                doPrint(target, value, args...);
            }
            assert(target == value);
            return target;
        } else {
            throw std::runtime_error("Can't call function size");
        }
    }


    CREATE_MEMBER_CHECK(select);
    template <typename... Args>
    int64_t select(Args... args) const {
        if constexpr (has_member_select<Reference>::value) {
            auto target = Reference::select(args...);
            auto value = t.select(args...);
            if (target != value) {
                doPrint(target, value, args...);
            }
            assert(target == value);
            return target;
        } else {
            throw std::runtime_error("Can't call function select");
        }
    }


    CREATE_MEMBER_CHECK(select_unchecked);
    template <typename... Args>
    int64_t select_unchecked(Args... args) const {
        if constexpr (has_member_select_unchecked<Reference>::value) {
            auto target = Reference::select_unchecked(args...);
            auto value = t.select_unchecked(args...);
            if (target != value) {
                doPrint(target, value, args...);
            }
            assert(target == value);
            return target;
        } else {
            throw std::runtime_error("Can't call function select_unchecked");
        }
    }


    CREATE_MEMBER_CHECK(rank);
    template <typename... Args>
    int64_t rank(Args... args) const {
        if constexpr (has_member_rank<Reference>::value) {
            auto target = Reference::rank(args...);
            auto value = t.rank(args...);
            if (target != value) {
                doPrint(target, value, args...);
            }
            assert(target == value);
            return target;
        } else {
            throw std::runtime_error("Can't call function rank");
        }
    }


    CREATE_MEMBER_CHECK(get);
    template <typename... Args>
    int64_t get(Args... args) const {
        if constexpr (has_member_get<Reference>::value) {
            auto target = Reference::get(args...);
            auto value = t.get(args...);
            if (target != value) {
                doPrint(target, value, args...);
            }
            assert(target == value);
            return target;
        } else {
            throw std::runtime_error("Can't call function get");
        }
    }


    CREATE_MEMBER_CHECK(traversed_base);
    template <typename... Args>
    int64_t traversed_base(Args... args) const {
        if constexpr (has_member_traversed_base<Reference>::value) {
            auto target = Reference::traversed_base(args...);
            auto value = t.traversed_base(args...);
            if (target != value) {
                doPrint(target, value, args...);
            }
            assert(target == value);
            return target;
        } else {
            throw std::runtime_error("Can't call function traversed_base");
        }
    }


    CREATE_MEMBER_CHECK(branch_offset);
    template <typename... Args>
    int64_t branch_offset(Args... args) const {
        if constexpr (has_member_branch_offset<Reference>::value) {
            auto target = Reference::branch_offset(args...);
            auto value = t.branch_offset(args...);
            if (target != value) {
                doPrint(target, value, args...);
            }
            assert(target == value);
            return target;
        } else {
            throw std::runtime_error("Can't call function branch_offset");
        }
    }


    CREATE_MEMBER_CHECK(branch_size);
    template <typename... Args>
    int64_t branch_size(Args... args) const {
        if constexpr (has_member_branch_size<Reference>::value) {
            auto target = Reference::branch_size(args...);
            auto value = t.branch_size(args...);
            if (target != value) {
                doPrint(target, value, args...);
            }
            assert(target == value);
            return target;
        } else {
            throw std::runtime_error("Can't call function branch_size");
        }
    }


    CREATE_MEMBER_CHECK(branch_offset_and_increment);
    template <typename... Args>
    int64_t branch_offset_and_increment(Args... args) {
        if constexpr (has_member_branch_offset_and_increment<Reference>::value) {
            auto target = Reference::branch_offset_and_increment(args...);
            auto value = t.branch_offset_and_increment(args...);
            if (target != value) {
                doPrint(target, value, args...);
            }
            assert(target == value);
            return target;
        } else {
            throw std::runtime_error("Can't call function branch_offset_and_increment");
        }
    }


    CREATE_MEMBER_CHECK(insert);
    template <typename... Args>
    int64_t insert(Args... args) {
        if constexpr (has_member_insert<Reference>::value) {
            Reference::insert(args...);
            auto result = t.insert(args...);
            return result;
        } else {
            throw std::runtime_error("Can't call function insert");
        }
    }

    CREATE_MEMBER_CHECK(new_relative_position);
    template <typename... Args>
    int64_t new_relative_position(Args... args) const {
        if constexpr (has_member_new_relative_position<Reference>::value) {
            auto target = Reference::new_relative_position(args...);
            auto value = t.new_relative_position(args...);
            if (target != value) {
                doPrint(target, value, args...);
            }
            assert(target == value);
            return target;
        } else {
            throw std::runtime_error("Can't call function new_relative_position");
        }
    }

    CREATE_MEMBER_CHECK(exit);
    template <typename... Args>
    int64_t exit(Args... args) {
        // not thread safe, needs additional lock here
        if constexpr (has_member_exit<Reference>::value) {
            auto target = Reference::exit(args...);
            auto value = t.exit(args...);
            if (target != value) {
                Reference::print_content(args...);
                t.print_content(args...);
                doPrint(target, value, args...);
            }
            assert(target == value);
            return target;
        } else {
            throw std::runtime_error("Can't call function exit");
        }
    }


    CREATE_MEMBER_CHECK(enter);
    template <typename... Args>
    void enter(Args... args) {
        // not thread safe, needs additional lock here
        if constexpr (has_member_enter<Reference>::value) {
            Reference::enter(args...);
            t.enter(args...);
        } else {
            throw std::runtime_error("Can't call function enter");
        }
    }

    CREATE_MEMBER_CHECK(print_content);
    template <typename... Args>
    void print_content(Args... args) const {
        // not thread safe, needs additional lock here
        if constexpr (has_member_print_content<Reference>::value) {
            Reference::print_content(args...);
            t.print_content(args...);
        } else {
            throw std::runtime_error("Can't call function print_content");
        }
    }


    Test t;
};

class VerboseTimer {
  public:
    VerboseTimer(string procedure_name) : procedure_name(procedure_name) {
        mg::common::logger->info("Started {}", procedure_name);
    }
    double finished() {
        double elapsed = timer.elapsed();
        mg::common::logger->info("Finished {} in {}  sec.", procedure_name, elapsed);
        return elapsed;
    }
    string procedure_name;
    Timer timer;
};

inline BOSS *dbg_succ_graph_constructor(const vector<string> &reads, size_t kmer_length) {
    auto graph_constructor
            = BOSSConstructor(kmer_length - 1); // because BOSS has smaller kmers

    graph_constructor.add_sequences(reads);

    return new BOSS(&graph_constructor);
}

#endif /* __UTILS_HPP__ */

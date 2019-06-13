//
//  utils.hpp
//  lossless_dbg
//
//  Created by Jan Studený on 11/03/2019.
//  Copyright © 2019 Jan Studený. All rights reserved.
//

#ifndef utils_h
#define utils_h

#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include "boss_construct.hpp"
#include "boss.hpp"
#include "kmer_extractor.hpp"
#include "sequence_io.hpp"

using namespace std;
using ll = long long;

#define x first
#define y second
#define all(x) begin(x),end(x)

namespace fs = std::filesystem;

#define local_file(filename) (fs::path(__FILE__).parent_path() / (filename))

//template<typename POD>
//std::ostream& serialize(std::ostream& os, std::vector<POD> const& v);
//
//template<typename POD>
//std::istream& deserialize(std::istream& is, std::vector<POD>& v);

struct d_t {
    template<typename T>
    d_t operator,(const T &first) {
        std::cerr << ' ' <<  x;
        return *this;
    }
};



#define PRINT_VAR(args ...) { d_t(), "|", __LINE__, "|", #args, ":", args, "\n"; }

template<typename POD>
inline std::istream &deserialize(std::istream &is, vector<POD> &v) {
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
                  "Can only deserialize POD types with this function");

    decltype(v.size()) size;
    is.read(reinterpret_cast<char*>(&size), sizeof(size));
    v.resize(size);
    is.read(reinterpret_cast<char*>(v.data()), v.size() * sizeof(POD));
    return is;
}

template<typename POD>
inline std::ostream &serialize(std::ostream &os, const vector<POD> &v) {
    // this only works on built in data types (PODs)
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
                  "Can only serialize POD types with this function");

    auto size = v.size();
    os.write(reinterpret_cast<char const*>(&size), sizeof(size));
    os.write(reinterpret_cast<char const*>(v.data()), v.size() * sizeof(POD));
    return os;
}

inline string decode_from_input(const vector<string>& input,const pair<string,int64_t>& path_id,int64_t kmer_length) {
    int64_t current_relative_id = 0;
    for(auto& sequence : input) {
        if (sequence.substr(0,kmer_length) == path_id.first) {
            if (current_relative_id == path_id.second) {
                return sequence;
            }
            else {
                current_relative_id++;
            }
        }
    }
    return "";
}

inline void save_string(const string &to_save, const string &filename) {
    ofstream myfile;
    myfile.open (filename);
    myfile << to_save;
    myfile.close();
}

inline void transform_to_fasta(const string &filename, const vector <string> &reads) {
    ofstream myfile;
    myfile.open (filename);
    for(auto& read : reads) {
        myfile << ">" << endl;
        myfile << read << endl;
    }
    myfile.close();
}


inline void write_reads_to_fasta(const vector <string> &reads, const string &filename) {
    transform_to_fasta(filename,reads);
}

inline void reduce_maps(std::map<int, int> &output, std::map<int, int> &input) {
    for (auto& X : input) {
        output[X.first] += X.second;
    }
}

inline int8_t encode(char c) {
    if (c == '#') return 6;//alphabet_decoder.alph_size;
    if (c == '$') return 0;
    return KmerExtractor::encode(c);
}
inline char decode(int8_t c) {
    if (c == 6) return '#';
    if (c == 0) return '$';
    return KmerExtractor::decode(c);
}

inline string& clamp_alphabet(string& text,const string& alphabet="$ACGTN",char replacement='N') {
    for(auto& c : text) {
        if (alphabet.find(c) == string::npos) {
            c = replacement;
        }
    }
    return text;
}

inline vector <string> read_reads_from_fasta(const string &filename) {
    vector<string> result;
    read_fasta_file_critical(
            filename,
            [&](kseq_t* read) {
                string read_seq = read->seq.s;
                clamp_alphabet(read_seq);
                result.push_back(read_seq);
            });
    return result;
}

template <typename... Args>
inline void doPrint(std::ostream& out, Args&&... args)
{
    ((out << ',' << std::forward<Args>(args)), ...);
}

/*
    - Multiple inheritance forces ambiguity of member names.
    - SFINAE is used to make aliases to member names.
    - Expression SFINAE is used in just one generic has_member that can accept
      any alias we pass it.
*/

//Variadic to force ambiguity of class members.  C++11 and up.
template <typename... Args> struct ambiguate : public Args... {};

//Non-variadic version of the line above.
//template <typename A, typename B> struct ambiguate : public A, public B {};

template<typename A, typename = void>
struct got_type : std::false_type {};

template<typename A>
struct got_type<A> : std::true_type {
    typedef A type;
};

template<typename T, T>
struct sig_check : std::true_type {};

template<typename Alias, typename AmbiguitySeed>
struct has_member {
    template<typename C> static char ((&f(decltype(&C::value))))[1];
    template<typename C> static char ((&f(...)))[2];

    //Make sure the member name is consistently spelled the same.
    static_assert(
            (sizeof(f<AmbiguitySeed>(0)) == 1)
            , "Member name specified in AmbiguitySeed is different from member name specified in Alias, or wrong Alias/AmbiguitySeed has been specified."
    );

    static bool const value = sizeof(f<Alias>(0)) == 2;
};

//Check for function with given name, any signature.
//Check for any member with given name, whether var, func, class, union, enum.
#define CREATE_MEMBER_CHECK(member)                                         \
                                                                            \
template<typename T, typename = std::true_type>                             \
struct Alias_##member;                                                      \
                                                                            \
template<typename T>                                                        \
struct Alias_##member <                                                     \
    T, std::integral_constant<bool, got_type<decltype(&T::member)>::value>  \
> { static const decltype(&T::member) value; };                             \
                                                                            \
struct AmbiguitySeed_##member { char member; };                             \
                                                                            \
template<typename T>                                                        \
struct has_member_##member {                                                \
    static const bool value                                                 \
        = has_member<                                                       \
            Alias_##member<ambiguate<T, AmbiguitySeed_##member>>            \
            , Alias_##member<AmbiguitySeed_##member>                        \
        >::value                                                            \
    ;                                                                       \
}





#undef protected
template <typename Test,typename Reference>
class IdentityComparator : public Reference {
public:
    template<typename ...Args>
    IdentityComparator(Args... args) : Reference(args...), t(args...) {}

    virtual ~IdentityComparator() {}

    CREATE_MEMBER_CHECK(offset);
    template<typename ...Args>
    int64_t offset(Args... args) const {
        if constexpr (has_member_offset<Reference>::value) {
            auto target = Reference::offset(args...);
            auto value = t.offset(args...);
            if (target != value) {
                doPrint(cout,target,value,args...);
            }
            assert(target==value);
            return target;
        }
        else {
            throw "Can't call function offset";
        }
    }


    CREATE_MEMBER_CHECK(size);
    template<typename ...Args>
    int64_t size(Args... args) const {
        if constexpr (has_member_size<Reference>::value) {
            auto target = Reference::size(args...);
            auto value = t.size(args...);
            if (target != value) {
                doPrint(cout,target,value,args...);
            }
            assert(target==value);
            return target;
        }
        else {
            throw "Can't call function size";
        }
    }


    CREATE_MEMBER_CHECK(select);
    template<typename ...Args>
    int64_t select(Args... args) const {
        if constexpr (has_member_select<Reference>::value) {
            auto target = Reference::select(args...);
            auto value = t.select(args...);
            if (target != value) {
                doPrint(cout,target,value,args...);
            }
            assert(target==value);
            return target;
        }
        else {
            throw "Can't call function select";
        }
    }


    CREATE_MEMBER_CHECK(select_unchecked);
    template<typename ...Args>
    int64_t select_unchecked(Args... args) const {
        if constexpr (has_member_select_unchecked<Reference>::value) {
            auto target = Reference::select_unchecked(args...);
            auto value = t.select_unchecked(args...);
            if (target != value) {
                doPrint(cout,target,value,args...);
            }
            assert(target==value);
            return target;
        }
        else {
            throw "Can't call function select_unchecked";
        }
    }


    CREATE_MEMBER_CHECK(rank);
    template<typename ...Args>
    int64_t rank(Args... args) const {
        if constexpr (has_member_rank<Reference>::value) {
            auto target = Reference::rank(args...);
            auto value = t.rank(args...);
            if (target != value) {
                doPrint(cout,target,value,args...);
            }
            assert(target==value);
            return target;
        }
        else {
            throw "Can't call function rank";
        }
    }


    CREATE_MEMBER_CHECK(get);
    template<typename ...Args>
    int64_t get(Args... args) const {
        if constexpr (has_member_get<Reference>::value) {
            auto target = Reference::get(args...);
            auto value = t.get(args...);
            if (target != value) {
                doPrint(cout,target,value,args...);
            }
            assert(target==value);
            return target;
        }
        else {
            throw "Can't call function get";
        }
    }


    CREATE_MEMBER_CHECK(traversed_base);
    template<typename ...Args>
    int64_t traversed_base(Args... args) const {
        if constexpr (has_member_traversed_base<Reference>::value) {
            auto target = Reference::traversed_base(args...);
            auto value = t.traversed_base(args...);
            if (target != value) {
                doPrint(cout,target,value,args...);
            }
            assert(target==value);
            return target;
        }
        else {
            throw "Can't call function traversed_base";
        }
    }


    CREATE_MEMBER_CHECK(branch_offset);
    template<typename ...Args>
    int64_t branch_offset(Args... args) const {
        if constexpr (has_member_branch_offset<Reference>::value) {
            auto target = Reference::branch_offset(args...);
            auto value = t.branch_offset(args...);
            if (target != value) {
                doPrint(cout,target,value,args...);
            }
            assert(target==value);
            return target;
        }
        else {
            throw "Can't call function branch_offset";
        }
    }


    CREATE_MEMBER_CHECK(branch_size);
    template<typename ...Args>
    int64_t branch_size(Args... args) const {
        if constexpr (has_member_branch_size<Reference>::value) {
            auto target = Reference::branch_size(args...);
            auto value = t.branch_size(args...);
            if (target != value) {
                doPrint(cout,target,value,args...);
            }
            assert(target==value);
            return target;
        }
        else {
            throw "Can't call function branch_size";
        }
    }


    CREATE_MEMBER_CHECK(branch_offset_and_increment);
    template<typename ...Args>
    int64_t branch_offset_and_increment(Args... args) {
        if constexpr (has_member_branch_offset_and_increment<Reference>::value) {
            auto target = Reference::branch_offset_and_increment(args...);
            auto value = t.branch_offset_and_increment(args...);
            if (target != value) {
                doPrint(cout,target,value,args...);
            }
            assert(target==value);
            return target;
        }
        else {
            throw "Can't call function branch_offset_and_increment";
        }
    }




    CREATE_MEMBER_CHECK(insert);
    template<typename ...Args>
    void insert(Args... args) {
        if constexpr (has_member_insert<Reference>::value) {
            Reference::insert(args...);
            t.insert(args...);
        }
        else {
            throw "Can't call function insert";
        }
    }

    CREATE_MEMBER_CHECK(new_relative_position);
    template<typename ...Args>
    int64_t new_relative_position(Args... args) const {
        if constexpr (has_member_new_relative_position<Reference>::value) {
            auto target = Reference::new_relative_position(args...);
            auto value = t.new_relative_position(args...);
            if (target != value) {
                doPrint(cout,target,value,args...);
            }
            assert(target==value);
            return target;
        }
        else {
            throw "Can't call function new_relative_position";
        }
    }

    CREATE_MEMBER_CHECK(exit);
    template<typename ...Args>
    int64_t exit(Args... args) {
        if constexpr (has_member_exit<Reference>::value) {
            auto target = Reference::exit(args...);
            auto value = t.exit(args...);
            if (target != value) {
                doPrint(cout,target,value,args...);
            }
            assert(target==value);
            return target;
        }
        else {
            throw "Can't call function exit";
        }
    }




    CREATE_MEMBER_CHECK(enter);
    template<typename ...Args>
    void enter(Args... args) {
        if constexpr (has_member_enter<Reference>::value) {
            Reference::enter(args...);
            t.enter(args...);
        }
        else {
            throw "Can't call function enter";
        }
    }


    Test t;
};

inline BOSS* dbg_succ_graph_constructor(const vector<string> &reads,
                                        size_t kmer_length) {

    auto graph_constructor = BOSSConstructor(kmer_length - 1);// because BOSS has smaller kmers

    graph_constructor.add_sequences(reads);

    return new BOSS(&graph_constructor);
}

#endif /* utils_h */

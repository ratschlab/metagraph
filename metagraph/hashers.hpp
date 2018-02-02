#ifndef __ANNOBLOOM__
#define __ANNOBLOOM__

#include <iostream>
#include <cassert>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <unordered_set>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <MurmurHash3.h>


namespace annotate {

class ExactFilter {
  public:
    ExactFilter(size_t num_hash_functions = 0, size_t size = 0) {
        num_hash_functions++;
        size++;
    }

    template <typename Object>
    bool find(Object *a, Object *b) const {
        return set.find(hash_(a, b)) != set.end();
    }

    template <typename Object>
    bool insert(Object *a, Object *b) {
        return !set.insert(hash_(a, b)).second;
    }

  private:
    template <typename Object>
    static uint64_t hash_(Object *a, Object *b) {
        __uint128_t hash = 0;
        MurmurHash3_x64_128(
            reinterpret_cast<const char*>(a),
            reinterpret_cast<const char*>(b) - reinterpret_cast<const char*>(a),
            0,
            &hash
        );
        uint64_t first_hash;
        std::copy(reinterpret_cast<const uint64_t*>(&hash),
                  reinterpret_cast<const uint64_t*>(&hash) + 1, &first_hash);
        return first_hash;
    }
    std::unordered_set<uint64_t> set;
};


class BloomFilter {
  public:
    typedef std::vector<uint64_t> big_int;

    BloomFilter(size_t num_hash_functions, size_t _n_bits = 0) {
        seeds.resize(num_hash_functions);
        std::iota(seeds.begin(), seeds.end(), 0);

        n_bits = _n_bits;
        if (n_bits > 0)
            bits.resize((n_bits >> 6) + 1);
    }

    size_t size() {
        return n_bits;
    }

    void resize(size_t new_size) {
        n_bits = new_size;
        bits.resize((n_bits >> 6) + 1);
    }

  public:
    template <typename T>
    bool find(T *a, T *b) const {
        return hash_helper(a, b);
    }

    template <typename T>
    bool find(const std::vector<T> &a) const {
        return find(&(*(a.begin())), &(*(a.end())));
    }
    template <typename T>
    bool find(const T &a) const {
        return find(&a, &a + sizeof(a));
    }

    /*
    template <typename T>
    bool annotate(T *a, T *b) {
        bool annot = 1;
        for (auto it = a; it != b; ++it) {
            annot &= find(*it);
        }
        return annot;
    }
    */

    template <typename T>
    bool insert(T *a, T *b) {
        return hash_helper_insert(a, b);
    }
    /*
    template <typename T>
    bool insert(const std::vector<T> &a) {
        return insert(&(*(a.begin())), &(*(a.end())));
    }
    */
    /*
    template <typename T>
    bool insert(const T &a) {
        return insert(&a, &a + sizeof(a));
    }
    */

    void serialize(std::ostream &out) const {
        /*
        //shift
        out.write(reinterpret_cast<const char*>(&shift), sizeof(shift));

        //seeds
        char *buffer;
        size_t size = seeds.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));
        buffer = (char*)malloc(seeds.size() * sizeof(seeds[0]));
        memcpy(buffer, reinterpret_cast<const char*>(&(*seeds.begin())), seeds.size() * sizeof(seeds[0]));
        //std::copy(seeds.begin(), seeds.end(), reinterpret_cast<size_t*>(buffer));
        out.write(buffer, seeds.size() * sizeof(seeds[0]));
        free(buffer);

        //bits
        size = bits.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));
        buffer = (char*)malloc(bits.size() * sizeof(bits[0]));
        memcpy(buffer, reinterpret_cast<const char*>(&(*bits.begin())), bits.size() * sizeof(bits[0]));
        //std::copy(bits.begin(), bits.end(), reinterpret_cast<size_t*>(buffer));
        out.write(buffer, bits.size() * sizeof(bits[0]));
        free(buffer);
        */
        out << n_bits << "\n";
        /*
        out << seeds.size() << "\n";
        for (auto it = seeds.begin(); it != seeds.end(); ++it) {
            out << *it << " ";
        }
        out << bits.size() << "\n";
        for (auto it = bits.begin(); it != bits.end(); ++it) {
            out << *it << " ";
        }
        */
        boost::archive::binary_oarchive oar(out);
        oar & seeds;
        oar & bits;
    }

    void deserialize(std::istream &in) {
        //size_t size;
        /*
        //shift
        in.read(reinterpret_cast<char*>(&shift), sizeof(shift));

        //seeds
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        seeds.resize(size);
        char *buffer;
        buffer = (char*)malloc(seeds.size() * sizeof(seeds[0]));
        in.read(buffer, seeds.size() * sizeof(seeds[0]));
        std::copy(buffer, buffer + seeds.size() * sizeof(seeds[0]), reinterpret_cast<char*>(&seeds[0]));
        free(buffer);

        //bits
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        bits.resize(size);
        buffer = (char*)malloc(bits.size() * sizeof(bits[0]));
        std::copy(buffer, buffer + bits.size() * sizeof(bits[0]), reinterpret_cast<char*>(&bits[0]));
        free(buffer);
        */
        in >> n_bits;
        /*
        in >> size;
        seeds.resize(size);
        for (auto it = seeds.begin(); it != seeds.end(); ++it) {
            in >> *it;
        }
        in >> size;
        bits.resize(size);
        for (auto it = bits.begin(); it != bits.end(); ++it) {
            in >> *it;
        }
        */
        boost::archive::binary_iarchive iar(in);
        iar & seeds;
        iar & bits;
    }

    bool operator==(const BloomFilter &a) const {
        if (n_bits != a.n_bits) {
            std::cerr << "Different number of bits\n";
            std::cerr << n_bits << " " << a.n_bits << "\n";
            return false;
        }
        if (seeds.size() != a.seeds.size()) {
            std::cerr << "Different number of seeds\n";
            std::cerr << seeds.size() << " " << a.seeds.size() << "\n";
            return false;
        }
        auto jt = a.seeds.begin();
        for (auto it = seeds.begin(); it != seeds.end(); ++it, ++jt) {
            if (*it != *jt) {
                std::cerr << "Seeds not equal\n";
                return false;
            }
        }
        if (bits.size() != a.bits.size()) {
            std::cerr << "Different number of blocks\n";
            std::cerr << bits.size() << " " << a.bits.size() << "\n";
            return false;
        }
        jt = a.bits.begin();
        for (auto it = bits.begin(); it != bits.end(); ++it, ++jt) {
            if (*it != *jt) {
                std::cerr << "Different bits\n";
                std::cerr << *it << " " << *jt << "\n";
                return false;
            }
        }
        return true;
    }

    bool operator!=(const BloomFilter &a) const {
        return !(*this == a);
    }

    double occupancy() const {
        size_t count = 0;
        for (auto it = bits.begin(); it != bits.end(); ++it) {
            count += __builtin_popcountl(*it);
        }
        return static_cast<double>(count) / (bits.size() * 64);
    }

  private:
    template <typename T>
    inline bool hash_helper(T *a, T *b) const {
        const char *begin = reinterpret_cast<const char*>(a);
        size_t len = static_cast<size_t>(b - a) * sizeof(T);

        __uint128_t hash = 0;
        uint64_t select;

        for (auto it = seeds.begin(); it != seeds.end(); ++it) {
            MurmurHash3_x64_128(begin, len, *it, &hash);
            select = hash % n_bits;
            auto jt = bits.begin() + (select >> 6);
            if (jt >= bits.end()) {
                std::cerr << "Out of bounds\n";
                exit(1);
            }
            if ((*jt | (1llu << (select % 64))) != *jt) {
                return false;
            }
        }
        return true;
    }

    template <typename T>
    inline bool hash_helper_insert(T *a, T *b) {
        assert(b >= a);

        const char *begin = reinterpret_cast<const char*>(a);
        size_t len = static_cast<size_t>(b - a) * sizeof(T);

        __uint128_t hash = 0;
        uint64_t select;
        bool might_contain = true;
        uint64_t last;
        for (auto it = seeds.begin(); it != seeds.end(); ++it) {
            MurmurHash3_x64_128(begin, len, *it, &hash);
            select = hash % n_bits;

            auto jt = bits.begin() + (select >> 6);
            assert(jt < bits.end());

            last = *jt;
            *jt |= 1llu << (select % 64);
            if (might_contain && *jt != last) {
                might_contain = false;
            }
        }
        return might_contain;
    }

    big_int bits;
    uint64_t n_bits = 0;
    std::vector<uint64_t> seeds;
};


template <class Filter = ExactFilter>
class HashAnnotation {
  public:
    HashAnnotation(size_t num_hash_functions = 0)
          : num_hash_functions_(num_hash_functions) {}

    // template <typename T>
    // HashAnnotation(T *a, T *b) {
    //     color_bits.resize(b - a - 1);
    //     for (auto it = a; it != b - 1; ++it) {
    //         color_bits[it - a] = *it;
    //     }
    // }
    // template <typename T>
    // HashAnnotation(const T *a, const T *b) {
    //     //init from vector of sizes, last is size of cont bit
    //     assert(b - a > 1);
    //     color_bits.resize(b - a - 1);
    //     for (auto it = a; it != b - 1; ++it) {
    //         color_bits[it - a] = Filter(*it);
    //     }
    // }

    void resize(size_t size) {
        assert(size > color_bits.size());
        color_bits.resize(size, Filter(num_hash_functions_));
    }

    Filter& operator[](size_t i) {
        return color_bits[i];
    }

    const Filter& operator[](size_t i) const {
        return color_bits[i];
    }

    /*
    template <class T>
    HashAnnotation(const std::vector<T> &a)
        : HashAnnotation(&(*(a.begin())), &(*(a.end()))) { }
    */

    /*
    template <typename T, typename S>
    void insert(const T &a, const S &begin, const S &end, bool&& change = false) {
        //insert a into Bloom filters with indices in *begin, ...
        for (auto it = begin; it != end; ++it) {
            if (*it >= color_bits.size()) {
                std::cerr << "Index " << *it << " >= " << color_bits.size() << "\n";
                exit(1);
            }
            color_bits[*it].insert(a);
        }
        if (change) {
            cont_bit.insert(a);
        }
    }
    */

    template <typename T, typename S>
    std::vector<size_t> insert(T *a, T *b, S begin, S end) {
        std::vector<size_t> annot((color_bits.size() >> 6) + 1);
        for (auto it = begin; it != end; ++it) {
            if (*it >= color_bits.size()) {
                std::cerr << "ERROR: Index " << *it << " >= " << color_bits.size() << "\n";
                exit(1);
            }
            if (color_bits[*it].insert(a, b)) {
                //might contain
                annot[*it >> 6] |= (1llu << (*it % 64));
            }
        }
        return annot;
    }

    template <typename T>
    std::vector<size_t> insert(T *a, T *b, size_t ind) {
        return insert(a, b, &ind, &ind + 1);
    }

    template <typename T, typename S>
    std::vector<size_t> find(T *a, T *b, S begin, S end) const {
        std::vector<size_t> annot((color_bits.size() >> 6) + 1);
        for (auto it = begin; it != end; ++it) {
            if (*it >= color_bits.size()) {
                std::cerr << "ERROR: Index " << *it << " >= " << color_bits.size() << "\n";
                exit(1);
            }
            if (color_bits[*it].find(a, b)) {
                //might contain
                annot[*it >> 6] |= (1llu << (*it % 64));
            }
        }
        return annot;
    }

    template <typename T>
    std::vector<size_t> find(T *a, T *b) const {
        std::vector<size_t> annot((color_bits.size() >> 6) + 1);
        for (size_t i = 0; i < color_bits.size(); ++i) {
            if (color_bits[i].find(a, b)) {
                annot[i >> 6] |= (1llu << (i % 64));
            }
        }
        return annot;
    }

    template <typename T>
    std::vector<size_t> find(T *a, T *b, size_t ind) const {
        return find(a, b, &ind, &ind + 1);
    }

    // struct Hash {
    //     __uint128_t hash;
    // };

    // std::vector<Hash> hash(const std::string &sequence, size_t k) const {
    //     assert(sequence.length() >= k);
    //     //TODO
    // }

    // template <typename T>
    // std::vector<size_t> annotate(const Hash &hashed_kmer) const {
    //     //TODO
    // }

    /*
    template <typename T, typename S>
    std::vector<size_t> insert(const std::vector<T> &a, const S &begin,
                               const S &end, bool&& change = false, bool query = false) {
        return insert(&(*(a.begin())), &(*(a.end())), begin, end, change, query);
    }
    template <typename T, typename S>
    std::vector<size_t> insert(T a, S begin, S end, bool change = false, bool query = false) {
        return insert(&a, &a + sizeof(a), begin, end, change, query);
    }
    */

    size_t size() const {
        return color_bits.size();
    }

    void serialize(std::ostream &out) const {
        size_t size = color_bits.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));
        //out << color_bits.size() << "\n";
        for (auto it = color_bits.begin(); it != color_bits.end(); ++it) {
            it->serialize(out);
        }
    }

    void deserialize(std::istream &in) {
        size_t size;
        //in >> size;
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        color_bits.resize(size);
        for (auto it = color_bits.begin(); it != color_bits.end(); ++it) {
            it->deserialize(in);
        }
    }

    bool operator==(const HashAnnotation<Filter> &a) {
        if (color_bits.size() != a.color_bits.size()) {
            std::cerr << "Different number of filters\n";
            return false;
        }
        auto jt = a.color_bits.begin();
        for (auto it = color_bits.begin(); it != color_bits.end(); ++it, ++jt) {
            if (*it != *jt) {
                std::cerr << "Color filter " << (it - color_bits.begin()) << " different\n";
                return false;
            }
        }
        return true;
    }

    bool operator!=(const HashAnnotation<Filter> &a) const {
        return !(*this == a);
    }

    void append_bit(size_t num_bits = 0) {
        color_bits.push_back(Filter(num_hash_functions_, num_bits));
    }

    std::vector<double> occupancy() const {
        std::vector<double> occ(color_bits.size());
        for (size_t i = 0; i < color_bits.size(); ++i) {
            occ[i] = color_bits[i].occupancy();
        }
        return occ;
    }

    static std::vector<uint64_t> merge_or(const std::vector<uint64_t> &a,
                                          const std::vector<uint64_t> &b) {
        assert(a.size() == b.size() && "ORing different sizes");

        std::vector<uint64_t> merged(a.size());
        for (size_t i = 0; i < merged.size(); ++i) {
            merged[i] = a[i] | b[i];
        }
        return merged;
    }

    static std::vector<uint64_t> merge_and(const std::vector<uint64_t> &a,
                                           const std::vector<uint64_t> &b) {
        assert(a.size() == b.size() && "ANDing different sizes");

        std::vector<uint64_t> merged(a.size());
        for (size_t i = 0; i < merged.size(); ++i) {
            merged[i] = a[i] & b[i];
        }
        return merged;
    }

    static uint64_t bigint_popcount(const std::vector<uint64_t> &a) {
        uint64_t popcount = 0;
        for (auto value : a) {
            popcount += __builtin_popcountl(value);
        }
        return popcount;
    }

  private:
    std::vector<Filter> color_bits;
    size_t num_hash_functions_;
};

} // namespace annotate

#endif // __ANNOBLOOM__

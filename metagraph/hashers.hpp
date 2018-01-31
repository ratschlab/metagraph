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
        ExactFilter(int size = 0) { size++; }
        template <typename Object>
        bool find(Object *a, Object *b) {
            return set.find(hash_(a, b)) != set.end();
        }
        template <typename Object>
        bool insert(Object *a, Object *b) {
            return !set.insert(hash_(a, b)).second;
        }
    private:
        template <typename Object>
        uint64_t hash_(Object *a, Object *b) {
            __uint128_t hash = 0;
            MurmurHash3_x64_128(
                    reinterpret_cast<char*>(a), 
                    reinterpret_cast<char*>(b) - reinterpret_cast<char*>(a), 
                    0, &hash);
            uint64_t first_hash;
            std::copy(reinterpret_cast<uint64_t*>(&hash), reinterpret_cast<uint64_t*>(&hash) + 1, &first_hash);
            return first_hash;
        }
        std::unordered_set<uint64_t> set;
};

template <int Size>
class BloomFilter {
    //Constructors and deconstructor
    public:
        typedef std::vector<uint64_t> big_int;
        BloomFilter() {
            seeds.resize(Size);
            std::iota(seeds.begin(), seeds.end(), 0);
        }
        BloomFilter(size_t n_bits) : BloomFilter() {
            if (n_bits == 0) {
                std::cerr << "Need at least one bit\n";
                exit(1);
            }
            this->n_bits = n_bits;
            bits.resize((n_bits >> 6) + 1);
        }
        ~BloomFilter() {
        }

        size_t size() {
            return n_bits;
        }

        void resize(size_t size) {
            this->n_bits = size;
            bits.resize((n_bits >> 6) + 1);
        }

    //Hasher
    private:
        template <typename T>
        inline bool hash_helper_(T *a, T *b, bool&& query = true) {
            const char *begin = reinterpret_cast<const char*>(a);
            int len = reinterpret_cast<const char*>(b) - begin;
            __uint128_t hash = 0;
            uint64_t select;
            big_int::iterator jt;
            bool might_contain = true;
            uint64_t last;
            for (auto it = seeds.begin(); it != seeds.end(); ++it) {
                MurmurHash3_x64_128(begin, len, *it, &hash);
                select = hash % n_bits;
                jt = bits.begin() + (select >> 6);
                if (jt >= bits.end()) {
                    std::cerr << "Out of bounds\n";
                    exit(1);
                }
                if (query) {
                    //if ((*jt | (1llu << (select % 64))) != *jt) {
                    if ((*jt | (1llu << (select % 64))) != *jt) {
                        return false;
                    }
                } else {
                    last = *jt;
                    *jt |= (1llu << (select % 64));
                    if (might_contain && *jt != last) {
                        might_contain = false;
                    }
                }
            }
            return might_contain;
        }
        
    public:
        template <typename T>
        bool find(T *a, T *b) {
            return hash_helper_(a, b, true);
        }

        template <typename T>
        bool find(const std::vector<T> &a) {
            return find(&(*(a.begin())), &(*(a.end())));
        }
        template <typename T>
        bool find(const T &a) {
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
            return hash_helper_(a, b, false);
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


        uint8_t* data() {
            return bits.data();
        }

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

        bool operator==(const BloomFilter<Size> &a) {
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

        bool operator!=(const BloomFilter<Size> &a) {
            return !(*this == a);
        }

        double occupancy() {
            size_t count = 0;
            for (auto it = bits.begin(); it != bits.end(); ++it) {
                count += __builtin_popcountl(*it);
            }
            return (double)count / (double)(bits.size() * 64);
        }

    private:
        big_int bits;
        uint64_t n_bits = 0;
        std::vector<uint64_t> seeds;
};

template <class Filter = ExactFilter>
class HashAnnotation {
    public:

        HashAnnotation() { }

        template <typename T>
        HashAnnotation(T *a, T *b) {
            color_bits.resize(b - a - 1);
            for (auto it = a; it != b - 1; ++it) {
                color_bits[it - a] = *it;
            }
            cont_bit = *b;
        }
        template <typename T>
        HashAnnotation(const T *a, const T *b) { 
            //init from vector of sizes, last is size of cont bit
            assert(b - a > 1);
            color_bits.resize(b - a - 1);
            for (auto it = a; it != b - 1; ++it) {
                color_bits[it - a] = Filter(*it);
            }
            cont_bit = Filter(*(b - 1));
        }

        void resize(size_t size) {
            assert(size > color_bits.size());
            color_bits.resize(size);
        }

        Filter& operator[](size_t i) {
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
        std::vector<size_t> insert(T *a, T *b, S begin, S end, bool change = false, bool query = false) {
            std::vector<size_t> annot((color_bits.size() >> 6) + 1);
            for (auto it = begin; it != end; ++it) {
                if (*it >= color_bits.size()) {
                    std::cerr << "Index " << *it << " >= " << color_bits.size() << "\n";
                    exit(1);
                }
                if (query) {
                    if (color_bits[*it].find(a, b)) {
                        //might contain
                        annot[*it >> 6] |= (1llu << (*it % 64));
                    }
                } else {
                    if (color_bits[*it].insert(a, b)) {
                        //might contain
                        annot[*it >> 6] |= (1llu << (*it % 64));
                    }
                }
            }
            if (change) {
                cont_bit.insert(a, b);
            }
            return annot;
        }

        template <typename T>
        std::vector<size_t> insert(T *a, T *b, const size_t &ind, bool change = false, bool query = false) {
            return insert(a, b, &ind, &ind + 1, change, query);
        }

        template <typename T>
        std::vector<size_t> set_junction(T *a, T *b) {
            size_t *false_ind = NULL;
            return insert(a, b, false_ind, false_ind, true, false);
        }

        template <typename T>
        std::vector<size_t> find(T *a, T *b) {
            std::vector<size_t> annot((color_bits.size() >> 6) + 1);
            for (size_t i = 0; i < color_bits.size(); ++i) {
                if (color_bits[i].find(a, b)) {
                    annot[i >> 6] |= (1llu << (i % 64));
                }
            }
            return annot;
        }

        template <typename T, typename S>
        std::vector<size_t> find(T *a, T *b, S begin, S end) {
            return insert(a, b, begin, end, false, true);
        }

        template <typename T>
        std::vector<size_t> find(T *a, T *b, const size_t &ind) {
            return insert(a, b, ind, false, true);
        }


        /*
        template <typename T, typename S>
        std::vector<size_t> insert(const std::vector<T> &a, const S &begin, const S &end, bool&& change = false, bool query = false) {
            return insert(&(*(a.begin())), &(*(a.end())), begin, end, change, query);
        }
        template <typename T, typename S>
        std::vector<size_t> insert(T a, S begin, S end, bool change = false, bool query = false) {
            return insert(&a, &a + sizeof(a), begin, end, change, query);
        }
        */

        template <typename T>
        void annotate(const T &a, std::vector<uint64_t> &annot, bool op_or_and = false) {
            for (auto it = color_bits.begin(); it != color_bits.end(); ++it) {
                if (op_or_and) {
                    //bitwise AND
                    if (!(it->find(a))) {
                        annot[(it - color_bits.begin()) >> 6] &= (-1llu ^ (1llu << ((it - color_bits.begin()) % 64)));
                    }
                } else {
                    //bitwise OR
                    if (it->find(a)) {
                        annot[(it - color_bits.begin()) >> 6] |= (1llu << ((it - color_bits.begin()) % 64));
                    }
                }
            }
        }

        size_t size() {
            return color_bits.size();
        }

        void serialize(std::ostream &out) const {
            size_t size = color_bits.size();
            out.write(reinterpret_cast<const char*>(&size), sizeof(size));
            //out << color_bits.size() << "\n";
            for (auto it = color_bits.begin(); it != color_bits.end(); ++it) {
                it->serialize(out);
            }
            //cont_bit.serialize(out);
        }

        void deserialize(std::istream &in) {
            size_t size;
            //in >> size;
            in.read(reinterpret_cast<char*>(&size), sizeof(size));
            color_bits.resize(size);
            for (auto it = color_bits.begin(); it != color_bits.end(); ++it) {
                it->deserialize(in);
            }
            //cont_bit.deserialize(in);
        }

        bool operator==(const HashAnnotation<Filter> &a) {
            if (color_bits.size() != a.color_bits.size()) {
                std::cerr << "Different number of filters\n";
                return false;
            }
            if (cont_bit != a.cont_bit) {
                std::cerr << "Continuity filter different\n";
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

        bool operator!=(const HashAnnotation<Filter> &a) {
            return !(*this == a);
        }

        void append_bit(size_t size = 0) {
            color_bits.push_back(Filter(size));
        }

        void set_cont(size_t size) {
            cont_bit = Filter(size);
        }

        std::vector<double> occupancy() {
            std::vector<double> occ(color_bits.size());
            for (size_t i = 0; i < color_bits.size(); ++i) {
                occ[i] = color_bits[i].occupancy();
            }
            return occ;
        }

        template <typename T>
        bool is_junction(T *a, T *b) {
            return cont_bit.find(a, b);
        }
    private:
        std::vector<Filter> color_bits;
        Filter cont_bit;

    public:
static bool merge_or (std::vector<uint64_t> &a, const std::vector<uint64_t> &b) {
    bool changed = false;
    if (a.size() != b.size()) {
        std::cerr << "ORing different sizes: " << a.size() << " " << b.size() << "\n";
    }
    auto jt = b.begin();
    uint64_t temp;
    for (auto it = a.begin(); it != a.end(); ++it, ++jt) {
        temp = *it | *jt;
        if (!changed && temp != *it)
            changed = true;
        *it = temp;
    }
    return changed;
}

static std::vector<uint64_t> merge_and (std::vector<uint64_t> &a, const std::vector<uint64_t> &b) {
    if (a.size() != b.size()) {
        std::cerr << "ANDing different sizes: " << a.size() << " " << b.size() << "\n";
    }
    std::vector<uint64_t> merged(a.size());
    auto jt = b.begin();
    auto kt = merged.begin();
    for (auto it = a.begin(); it != a.end(); ++it, ++jt, ++kt) {
        *kt = *it & *jt;
    }
    return merged;
}


static uint64_t bigint_popcount(std::vector<uint64_t> &a) {
    uint64_t popcount = 0;
    for (auto it = a.begin(); it != a.end(); ++it) {
        popcount += __builtin_popcountl(*it);
    }
    return popcount;
}

};

} // namespace annotate
#endif // __ANNOBLOOM__

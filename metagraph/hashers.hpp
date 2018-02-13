#ifndef __ANNOBLOOM__
#define __ANNOBLOOM__

#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>


namespace annotate {

//MERGING

typedef std::vector<uint64_t> MultiHash;

std::vector<uint64_t> merge_or(const std::vector<uint64_t> &a,
                               const std::vector<uint64_t> &b);

std::vector<uint64_t> merge_and(const std::vector<uint64_t> &a,
                                const std::vector<uint64_t> &b);

uint64_t popcount(const std::vector<uint64_t> &a);

bool test_bit(const std::vector<uint64_t> &a, size_t i);

void set_bit(std::vector<uint64_t> &a, size_t i);

bool equal(const std::vector<uint64_t> &a, const std::vector<uint64_t> &b);


class HashIterator {
  public:
    virtual ~HashIterator() {}

    virtual HashIterator& operator++() = 0;
    virtual bool is_end() const = 0;

    virtual const MultiHash& operator*() const = 0;
    virtual const MultiHash* operator->() const = 0;
};

class CyclicMultiHash {
  public:
    CyclicMultiHash(const char *data, size_t k, size_t num_hash);

    CyclicMultiHash(const std::string &sequence, size_t num_hash)
          : CyclicMultiHash(&sequence[0], sequence.length(), num_hash) {}

    CyclicMultiHash(CyclicMultiHash &&other) = default;
    CyclicMultiHash& operator=(CyclicMultiHash &&other) = default;

    ~CyclicMultiHash();

    void update(char next);
    void reverse_update(char prev);

    const MultiHash& get_hash() const { return hashes_; }

  private:
    MultiHash hashes_;
    size_t k_;

    std::string cache_;
    size_t begin_;
    //using void* to prevent including cyclichasher.h here
    std::vector<void*> chashers_;
};

class CyclicHashIterator : public HashIterator {
  public:
    CyclicHashIterator(const char *begin, const char *end,
                       size_t k, size_t num_hash);
    CyclicHashIterator(const std::string &sequence, size_t k, size_t num_hash);

    CyclicHashIterator(CyclicHashIterator &&other) = default;
    CyclicHashIterator& operator=(CyclicHashIterator &&other) = default;

    CyclicHashIterator& operator++();
    bool is_end() const { return next_ >= end_; }

    const MultiHash& operator*() const { return hasher_.get_hash(); }
    const MultiHash* operator->() const { return &hasher_.get_hash(); }

  private:
    CyclicMultiHash hasher_;
    const char *next_;
    const char *end_;
};


class ExactFilter {
  public:
    ExactFilter(size_t num_hash_functions = 0, size_t size = 0) {
        num_hash_functions++;
        size++;
    }

    template <typename T>
    bool find(const T *begin, const T *end) const {
        return set_.find(std::string(reinterpret_cast<const char*>(begin),
                                     reinterpret_cast<const char*>(end))) != set_.end();
    }

    bool find(const MultiHash &hash) const {
        std::cerr << "This function must not be used due to collisions\n";
        std::cerr << "This function must not be used due to collisions\n";
        std::cerr << &hash << "\n";
        exit(1);
        return false;
    }

    template <typename T>
    bool insert(const T *begin, const T *end) {
        return !set_.emplace(reinterpret_cast<const char*>(begin),
                             reinterpret_cast<const char*>(end)).second;
    }

    bool insert(const MultiHash &hash) {
        std::cerr << "This function must not be used due to collisions\n";
        std::cerr << &hash << "\n";
        exit(1);
        return false;
    }

  private:
    std::unordered_set<std::string> set_;
};


class BloomFilter {
  public:
    BloomFilter(size_t num_hash_functions, size_t n_bits = 0)
          : n_bits_(n_bits) {
        seeds.resize(num_hash_functions);
        std::iota(seeds.begin(), seeds.end(), 0);

        if (n_bits_ > 0)
            bits.resize((n_bits_ >> 6) + 1);
    }

    size_t size() const {
        return n_bits_;
    }

    void resize(size_t new_size) {
        n_bits_ = new_size;
        bits.resize((n_bits_ >> 6) + 1);
    }

  public:
    template <typename T>
    bool find(T *a, T *b) const {
        assert(b >= a);
        return find(compute_hash(a, b));
    }

    bool find(const MultiHash &multihash) const {
        for (auto hash : multihash) {
            if (!test_bit(bits, hash % n_bits_)) {
                return false;
            }
        }
        return true;
        //return equal(merge_or(bits, annotate(hash, n_bits_)), bits);
    }

    bool insert(const MultiHash &multihash) {
        bool might_contain = true;
        for (auto hash : multihash) {
            if (might_contain && !test_bit(bits, hash % n_bits_)) {
                might_contain = false;
            }
            set_bit(bits, hash % n_bits_);
        }
        return might_contain;
        /*
        auto merged = merge_or(bits, annotate(hash, n_bits_));
        bool might_contain = equal(merged, bits);
        std::copy(merged.begin(), merged.end(), bits.begin());
        return might_contain;
        */
    }

    template <typename T>
    bool insert(T *a, T *b) {
        assert(b >= a);
        return insert(compute_hash(a, b));
    }

    void serialize(std::ostream &out) const {
        out << n_bits_ << "\n";
        boost::archive::binary_oarchive oar(out);
        oar & seeds;
        oar & bits;
    }

    void deserialize(std::istream &in) {
        in >> n_bits_;
        boost::archive::binary_iarchive iar(in);
        iar & seeds;
        iar & bits;
    }

    bool operator==(const BloomFilter &a) const {
        if (n_bits_ != a.n_bits_) {
            std::cerr << "Different number of bits\n";
            std::cerr << n_bits_ << " " << a.n_bits_ << "\n";
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
            count += __builtin_popcountll(*it);
        }
        return static_cast<double>(count) / (bits.size() * 64);
    }

  private:
    template <typename T>
    MultiHash compute_hash(const T *begin, const T *end) const {
        CyclicMultiHash hasher(
            reinterpret_cast<const char*>(begin),
            (end - begin) * sizeof(T),
            seeds.size()
        );
        return hasher.get_hash();
    }

    std::vector<uint64_t> bits;
    uint64_t n_bits_ = 0;
    std::vector<uint64_t> seeds;
};


template <class Filter = ExactFilter>
class HashAnnotation {
  public:
    HashAnnotation(size_t num_hash_functions = 0)
          : num_hash_functions_(num_hash_functions) {}

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

    template <typename T, typename S>
    std::vector<uint64_t> insert(T *a, T *b, S begin, S end) {
        std::vector<uint64_t> annot((color_bits.size() >> 6) + 1);
        for (auto it = begin; it != end; ++it) {
            if (*it >= color_bits.size()) {
                std::cerr << "ERROR: Index " << *it << " >= " << color_bits.size() << "\n";
                exit(1);
            }
            if (color_bits[*it].insert(a, b)) {
                //might contain
                set_bit(annot, *it);
                //annot[*it >> 6] |= (1llu << (*it % 64));
            }
        }
        return annot;
    }

    template <typename T>
    std::vector<uint64_t> insert(T *a, T *b, size_t ind) {
        return insert(a, b, &ind, &ind + 1);
    }

    template <typename S>
    std::vector<uint64_t> insert(const MultiHash &hash, S begin, S end) {
        std::vector<uint64_t> annot((color_bits.size() >> 6) + 1);
        for (auto it = begin; it != end; ++it) {
            if (*it >= color_bits.size()) {
                std::cerr << "ERROR: Index " << *it << " >= " << color_bits.size() << "\n";
                exit(1);
            }
            if (color_bits[*it].insert(hash)) {
                set_bit(annot, *it);
                //annot[*it >> 6] |= (1llu << (*it % 64));
            }
        }
        return annot;
    }

    std::vector<uint64_t> insert(const MultiHash &hash, size_t ind) {
        return insert(hash, &ind, &ind + 1);
    }

    template <typename T, typename S>
    std::vector<uint64_t> find(T *a, T *b, S begin, S end) const {
        std::vector<uint64_t> annot((color_bits.size() >> 6) + 1, 0);
        for (auto it = begin; it != end; ++it) {
            if (*it >= color_bits.size()) {
                std::cerr << "ERROR: Index " << *it << " >= " << color_bits.size() << "\n";
                exit(1);
            }
            if (color_bits[*it].find(a, b)) {
                //might contain
                set_bit(annot, *it);
                //annot[*it >> 6] |= (1llu << (*it % 64));
            }
        }
        return annot;
    }

    template <typename T>
    std::vector<uint64_t> find(T *a, T *b) const {
        std::vector<uint64_t> annot((color_bits.size() >> 6) + 1);
        for (size_t i = 0; i < color_bits.size(); ++i) {
            if (color_bits[i].find(a, b)) {
                set_bit(annot, i);
                //annot[i >> 6] |= (1llu << (i % 64));
            }
        }
        return annot;
    }

    template <typename T>
    std::vector<uint64_t> find(T *a, T *b, size_t ind) const {
        return find(a, b, &ind, &ind + 1);
    }

    template <typename S>
    std::vector<uint64_t> find(const MultiHash &hash, S begin, S end) const {
        std::vector<uint64_t> annot((color_bits.size() >> 6) + 1);
        for (auto it = begin; it != end; ++it) {
            if (*it >= color_bits.size()) {
                std::cerr << "ERROR: Index " << *it << " >= " << color_bits.size() << "\n";
                exit(1);
            }
            if (color_bits[*it].find(hash)) {
                set_bit(annot, *it);
                //annot[*it >> 6] |= (1llu << (*it % 64));
            }
        }
        return annot;
    }

    std::vector<uint64_t> find(const MultiHash &hash, size_t ind) const {
        return find(hash, &ind, &ind + 1);
    }

    std::vector<uint64_t> find(const MultiHash &hash) const {
        std::vector<uint64_t> annot((color_bits.size() >> 6) + 1);
        for (size_t i = 0; i < color_bits.size(); ++i) {
            if (color_bits[i].find(hash)) {
                set_bit(annot, i);
            }
        }
        return annot;
    }


    size_t size() const {
        return color_bits.size();
    }

    size_t num_hash_functions() const {
        return num_hash_functions_;
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

    bool operator==(const HashAnnotation<Filter> &a) const {
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

  private:
    std::vector<Filter> color_bits;
    size_t num_hash_functions_;
};

} // namespace annotate

#endif // __ANNOBLOOM__

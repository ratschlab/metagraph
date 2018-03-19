#ifndef __ANNOBLOOM__
#define __ANNOBLOOM__

#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/archive/impl/basic_binary_oprimitive.ipp>
#include <boost/archive/impl/basic_binary_iprimitive.ipp>

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


template <typename MultiHash>
class MultiHasher {
  public:
    virtual ~MultiHasher() {}
    virtual bool reinitialize(const char *data, size_t k, size_t num_hash) = 0;
    virtual const MultiHash& get_hash() const = 0;
};

class DummyHasher : MultiHasher<std::string> {
  public:
    DummyHasher(const char *data, size_t k, size_t num_hash) {
        reinitialize(data, k, num_hash);
    }

    bool reinitialize(const char *data, size_t k, size_t num_hash) {
        hash_.assign(data, k);
        num_hash++;
        return true;
    }

    const std::string& get_hash() const { return hash_; }

  private:
    std::string hash_;
};

class CyclicMultiHash : public MultiHasher<MultiHash> {
  public:
    CyclicMultiHash(const char *data, size_t k, size_t num_hash);

    CyclicMultiHash(const std::string &sequence, size_t num_hash)
          : CyclicMultiHash(&sequence[0], sequence.length(), num_hash) {}

    CyclicMultiHash(CyclicMultiHash &&other) = default;
    CyclicMultiHash& operator=(CyclicMultiHash &&other) = default;

    ~CyclicMultiHash();

    bool reinitialize(const char *data, size_t k, size_t num_hash);

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
    bool is_end() const { return next_ > end_; }

    const MultiHash& operator*() const { return hasher_.get_hash(); }
    const MultiHash* operator->() const { return &hasher_.get_hash(); }

  private:
    CyclicMultiHash hasher_;
    const char *next_;
    const char *end_;
};


class ExactFilter {
  public:
    explicit ExactFilter(size_t size = 0) { size++; }

    bool find(const std::string &string) const {
        return set_.find(string) != set_.end();
    }

    bool insert(const std::string &string) {
        return !set_.emplace(string).second;
    }

  private:
    std::unordered_set<std::string> set_;
};


class BloomFilter {
  public:
    explicit BloomFilter(size_t n_bits = 0);

    size_t size() const { return n_bits_; }
    void resize(size_t new_size);

    bool find(const MultiHash &multihash) const;
    bool insert(const MultiHash &multihash);

    void serialize(std::ostream &out) const;
    void deserialize(std::istream &in);

    bool operator==(const BloomFilter &a) const;
    bool operator!=(const BloomFilter &a) const { return !operator==(a); }

    double occupancy() const;

  private:
    std::vector<uint64_t> bits;
    uint64_t n_bits_ = 0;
};


template <class Filter, class Hash, class Hasher>
class HashAnnotation {
  public:
    explicit HashAnnotation(size_t num_hash_functions = 0)
          : num_hash_functions_(num_hash_functions) {}

    ~HashAnnotation() {
        if (hasher_)
            delete hasher_;
    }

    void resize(size_t size) {
        assert(size > color_bits.size());
        color_bits.resize(size, Filter());
    }

    size_t size() const { return color_bits.size(); }

    size_t get_size(size_t i) const { return color_bits.at(i).size(); }

    size_t num_hash_functions() const { return num_hash_functions_; }

    Filter& operator[](size_t i) { return color_bits[i]; }
    const Filter& operator[](size_t i) const { return color_bits[i]; }

    template <typename T, typename S>
    std::vector<uint64_t> insert(T *a, T *b, S begin, S end) {
        std::vector<uint64_t> annot((color_bits.size() >> 6) + 1);
        for (auto it = begin; it != end; ++it) {
            if (*it >= color_bits.size()) {
                std::cerr << "ERROR: Index " << *it << " >= " << color_bits.size() << "\n";
                exit(1);
            }
            if (color_bits[*it].insert(compute_hash(a, b))) {
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
    std::vector<uint64_t> insert(const Hash &hash, S begin, S end) {
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

    std::vector<uint64_t> insert(const Hash &hash, size_t ind) {
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
            if (color_bits[*it].find(compute_hash(a, b))) {
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
            if (color_bits[i].find(compute_hash(a, b))) {
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
    std::vector<uint64_t> find(const Hash &hash, S begin, S end) const {
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

    std::vector<uint64_t> find(const Hash &hash, size_t ind) const {
        return find(hash, &ind, &ind + 1);
    }

    std::vector<uint64_t> find(const Hash &hash) const {
        std::vector<uint64_t> annot((color_bits.size() >> 6) + 1);
        for (size_t i = 0; i < color_bits.size(); ++i) {
            if (color_bits[i].find(hash)) {
                set_bit(annot, i);
            }
        }
        return annot;
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

    bool operator==(const HashAnnotation<Filter, Hash, Hasher> &a) const {
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

    bool operator!=(const HashAnnotation<Filter, Hash, Hasher> &a) const {
        return !(*this == a);
    }

    void append_bit(size_t filter_size = 0) {
        color_bits.push_back(Filter(filter_size));
    }

    std::vector<double> occupancy() const {
        std::vector<double> occ(color_bits.size());
        for (size_t i = 0; i < color_bits.size(); ++i) {
            occ[i] = color_bits[i].occupancy();
        }
        return occ;
    }

    Hash compute_hash(const std::string &sequence) const {
        return compute_hash(sequence.data(), sequence.data() + sequence.size());
    }

    template <typename T>
    Hash compute_hash(const T *begin, const T *end) const {
        Hasher *&hasher = const_cast<Hasher*&>(hasher_);
        if (hasher && !hasher->reinitialize(reinterpret_cast<const char*>(begin),
                                            (end - begin) * sizeof(T),
                                            num_hash_functions_)) {
            delete hasher;
            hasher = NULL;
        }
        if (!hasher) {
            hasher = new Hasher(reinterpret_cast<const char*>(begin),
                                (end - begin) * sizeof(T),
                                num_hash_functions_);
        }
        return hasher->get_hash();
    }

  private:
    std::vector<Filter> color_bits;
    size_t num_hash_functions_;
    Hasher *hasher_ = NULL;
};


typedef HashAnnotation<BloomFilter, MultiHash, CyclicMultiHash> BloomHashAnnotation;
//typedef HashAnnotation<ExactFilter, std::string, DummyHasher> ExactHashAnnotation;

class ExactHashAnnotation {
    public:
    ExactHashAnnotation() : num_columns_(0) { }

    template <typename T>
    std::vector<uint64_t> find(const T *begin, const T *end, long long i = -1) const {
        return find(compute_hash(begin, end), i);
    }

    template <typename T>
    std::vector<uint64_t> insert(const T *begin, const T *end, size_t i) {
        return insert(compute_hash(begin, end), i);
    }

    std::vector<uint64_t> find(const std::string &kmer, long long i = -1) const {
        std::vector<uint64_t> annot((num_columns_ + 63) >> 6, 0);
        auto it = kmer_map_.find(kmer);
        if (it == kmer_map_.end())
            return annot;
        if (i == -1) {
            for (auto &index : it->second) {
                set_bit(annot, index);
            }
        } else {
            if (it->second.find(i) != it->second.end())
                set_bit(annot, i);
        }
        return annot;
    }

    std::vector<uint64_t> insert(const std::string &kmer, size_t i) {
        if (i < -1llu)
            num_columns_ = std::max(num_columns_, i + 1);
        std::vector<uint64_t> annot((num_columns_ + 63) >> 6, 0);
        auto &indices = kmer_map_[kmer];
        if (i < -1llu)
            indices.insert(i);
        for (auto &index : indices) {
            set_bit(annot, index);
        }
        return annot;
    }

    size_t size() const {
        return num_columns_;
    }

    void resize(size_t size) {
        num_columns_ = size;
    }

    void append_bit() {
        num_columns_++;
    }

    template <typename T>
    std::string compute_hash(const T *begin, const T *end) const {
        return std::string(reinterpret_cast<const char*>(begin),
                           reinterpret_cast<const char*>(end));
    }

    void serialize(std::ostream &out) const {
        boost::archive::binary_oarchive oarch(out);
        oarch & kmer_map_;
        oarch & num_columns_;
    }
    void serialize(const std::string &filename) const {
        std::ofstream fout(filename);
        serialize(fout);
        fout.close();
    }

    void load(std::istream &in) {
        boost::archive::binary_iarchive iarch(in);
        iarch & kmer_map_;
        iarch & num_columns_;
    }

    void load(const std::string &filename) {
        std::ifstream fin(filename);
        load(fin);
        fin.close();
    }

    friend class PreciseAnnotator;
    private:
        std::unordered_map<std::string, std::set<size_t>> kmer_map_;
        size_t num_columns_;
};


} // namespace annotate

#endif // __ANNOBLOOM__

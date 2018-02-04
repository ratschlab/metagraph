#ifndef __ANNOBLOOM__
#define __ANNOBLOOM__

#include <iostream>
#include <cassert>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_set>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <MurmurHash3.h>


namespace annotate {

//MERGING

std::vector<uint64_t> merge_or(const std::vector<uint64_t> &a,
                                      const std::vector<uint64_t> &b);
std::vector<uint64_t> merge_and(const std::vector<uint64_t> &a,
                                       const std::vector<uint64_t> &b);

uint64_t popcount(const std::vector<uint64_t> &a);

//multihash container
struct MultiHash {
    MultiHash(size_t num_hash = 0) : hashes(num_hash) { }

    template <typename T>
    MultiHash(const T *data, const size_t num_hash = 0) : hashes(num_hash) {
        std::copy(data, data + num_hash, hashes.begin());
    }

    //vector to store hashes (one per seed/hash function)
    std::vector<size_t> hashes;
};

//HASH WRAPPER STRUCTS
struct {
    template <typename T, typename S>
    void operator()(const T *data, const size_t len, const uint32_t seed, S *hash) {
        //WARNING: make sure that the size of space allocated to hash is at least 8
        MurmurHash3_x64_128(
                reinterpret_cast<const char*>(data),
                len * sizeof(T),
                seed,
                reinterpret_cast<uint64_t*>(hash)
        );
    }
} Murmur3Hasher;

//same interface as ntHash, for turning any hash function into an iterative multihasher
class HashIterator {
    public:
      HashIterator& operator++() {
          for (size_t i = 0; i < hashes_.size(); ++i) {
              //use index as seed
              Murmur3Hasher(seq_cur, k_, i, &hashes_[i]);
          }
          seq_cur++;
          return *this;
      }

      HashIterator(const char *seq_begin, const char *seq_end, const size_t num_hash, const size_t k)
          : seq_begin(seq_begin),
            seq_cur(seq_begin),
            seq_end(seq_end),
            hashes_(num_hash),
            k_(k) {
      }

      HashIterator(const std::string &sequence, const size_t num_hash, const size_t k)
          : seq_begin(sequence.c_str()),
            seq_cur(seq_begin),
            seq_end(sequence.c_str() + sequence.length()),
            hashes_(num_hash),
            k_(k) {
          end_ = std::make_unique<HashIterator>(seq_end - k, seq_end, num_hash, k);
          operator++();
      }

      bool operator==(const HashIterator &that) {
          return
              seq_begin == that.seq_begin
              && seq_cur == that.seq_cur
              && seq_end == that.seq_end
              && k_ == that.k_
              && hashes_.size() == that.hashes_.size();
      }

      bool operator!=(const HashIterator &that) {
          return
              seq_begin != that.seq_begin
              || seq_cur != that.seq_cur
              || seq_end != that.seq_end
              || k_ != that.k_
              || hashes_.size() != that.hashes_.size();
      }

      const HashIterator& end() const {
          return *end_;
      }

      size_t size() const {
          return hashes_.size();
      }

      const uint64_t* operator*() const {
          return hashes_.data();
      }

      size_t pos() const {
          return seq_cur - seq_begin;
      }

    private:
      const char *seq_begin, *seq_cur, *seq_end;
      std::vector<uint64_t> hashes_;
      size_t k_;
      std::unique_ptr<HashIterator> end_;
};


template <class HashIterator>
std::vector<MultiHash> hash(HashIterator &hash_it);

std::vector<MultiHash> hash_murmur(
        const std::string &sequence,
        const size_t num_hash,
        const size_t k);

std::vector<size_t> annotate(MultiHash &multihash, const size_t max_size);

class ExactFilter {
  public:
    ExactFilter(size_t num_hash_functions = 0, size_t size = 0) {
        num_hash_functions++;
        size++;
    }

    template <typename Object>
    bool find(Object *a, Object *b) const {
        return set_.find(hash_(a, b)) != set_.end();
    }

    template <typename Object>
    bool insert(Object *a, Object *b) {
        return !set_.insert(hash_(a, b)).second;
    }

  private:
    std::unordered_set<uint64_t> set_;

    template <typename Object>
    static uint64_t hash_(Object *a, Object *b) {
        __uint128_t hash = 0;
        Murmur3Hasher(a, b - a, 0, &hash);
        uint64_t first_hash;
        std::copy(reinterpret_cast<const uint64_t*>(&hash),
                  reinterpret_cast<const uint64_t*>(&hash) + 1, &first_hash);
        return first_hash;
    }
};


class BloomFilter {
  public:
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

    template <typename T>
    bool insert(T *a, T *b) {
        return hash_helper_insert(a, b);
    }

    void serialize(std::ostream &out) const {
        out << n_bits << "\n";
        boost::archive::binary_oarchive oar(out);
        oar & seeds;
        oar & bits;
    }

    void deserialize(std::istream &in) {
        in >> n_bits;
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
        __uint128_t hash = 0;
        uint64_t select;

        for (auto it = seeds.begin(); it != seeds.end(); ++it) {
            Murmur3Hasher(a, b - a, *it, &hash);
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

        __uint128_t hash = 0;
        uint64_t select;
        bool might_contain = true;
        uint64_t last;
        for (auto it = seeds.begin(); it != seeds.end(); ++it) {
            Murmur3Hasher(a, b - a, *it, &hash);
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

    std::vector<uint64_t> bits;
    uint64_t n_bits = 0;
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
        std::vector<size_t> annot((color_bits.size() >> 6) + 1, 0);
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



  private:
    std::vector<Filter> color_bits;
    size_t num_hash_functions_;
};

} // namespace annotate

#endif // __ANNOBLOOM__

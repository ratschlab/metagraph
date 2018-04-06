#ifndef __ANNOTATE_COLOR_COMPRESSED_HPP__
#define __ANNOTATE_COLOR_COMPRESSED_HPP__

#include <cache.hpp>
#include <lru_cache_policy.hpp>

#include "annotate.hpp"
#include "dbg_succinct.hpp"


namespace annotate {

class ColorCompressed : public AnnotationCategory<std::set<std::string>> {
  public:
    typedef std::set<std::string> SetStr;

    ColorCompressed(const DBG_succ &graph, size_t cache_size = 1);

    ColorCompressed(uint64_t graph_size, size_t cache_size = 1);

    // Merge constructor
    //ColorCompressed(const DBG_succ &graph,
    //                const std::vector<ColorCompressed> &categories,
    //                const std::vector<std::vector<size_t>> &merge_plan);

    ~ColorCompressed() { release(); }

    SetStr get(Index i) const;

    std::vector<uint64_t> get_row(Index i) const;

    bool has_label(Index i, const SetStr &label) const;
    bool has_label(Index i, const std::string &label) const;

    void set_label(Index i, const SetStr &label);
    void add_label(Index i, const std::string &label);
    void add_labels(const std::string &sequence, const SetStr &labels);

    bool load(const std::string &filename);
    void serialize(const std::string &filename) const;

    std::vector<std::string> get_label_names() const { return id_to_label_; }

  private:
    const DBG_succ *graph_;
    void release();
    void flush();
    void flush(uint32_t curr_id, sdsl::bit_vector *annotation_curr);
    sdsl::bit_vector* inflate_column(const uint32_t id) const;

    uint64_t graph_size_;
    std::unordered_map<std::string, uint32_t> label_to_id_;
    std::vector<std::string> id_to_label_;
    std::vector<sdsl::sd_vector<>*> bitmatrix_;

    caches::fixed_sized_cache<uint32_t,
                              sdsl::bit_vector*,
                              caches::LRUCachePolicy<uint32_t>> cached_colors_;
};

} // namespace annotate

#endif // __ANNOTATE_COLOR_COMPRESSED_HPP__

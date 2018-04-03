#ifndef __ANNOTATE_COLOR_COMPRESSED_HPP__
#define __ANNOTATE_COLOR_COMPRESSED_HPP__

#include "annotate.hpp"
#include "dbg_succinct.hpp"


namespace annotate {

class ColorCompressed : public AnnotationCategory<std::set<std::string>> {
  public:
    typedef std::set<std::string> SetStr;

    ColorCompressed(const DBG_succ &graph)
          : graph_(&graph), graph_size_(graph.num_edges() + 1), annotation_curr_(NULL) {}

    ColorCompressed(size_t graph_size)
          : graph_(NULL), graph_size_(graph_size), annotation_curr_(NULL) {}

    // Merge constructor
    //ColorCompressed(const DBG_succ &graph,
    //                const std::vector<ColorCompressed> &categories,
    //                const std::vector<std::vector<size_t>> &merge_plan);

    ~ColorCompressed() { release(); }

    SetStr get(Index i) const;

    std::vector<uint64_t> get_row(Index i) const;

    void serialize_uncompressed_rows(const std::string &filename) const;

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
    sdsl::bit_vector* inflate_column(const uint32_t id) const;

    uint64_t graph_size_;
    std::unordered_map<std::string, uint32_t> label_to_id_;
    std::vector<std::string> id_to_label_;
    std::vector<sdsl::sd_vector<>*> bitmatrix_;

    uint32_t curr_id_;
    sdsl::bit_vector *annotation_curr_;
};

} // namespace annotate

#endif // __ANNOTATE_COLOR_COMPRESSED_HPP__

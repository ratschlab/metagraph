#ifndef __ANNOTATE_ROW_COMPRESSED_HPP__
#define __ANNOTATE_ROW_COMPRESSED_HPP__

#include "annotate.hpp"
#include "dbg_succinct.hpp"


namespace annotate {

class RowCompressed : public AnnotationCategory<std::set<std::string>> {
  public:
    typedef std::set<std::string> SetStr;

    RowCompressed(const DBG_succ &graph);

    RowCompressed(uint64_t graph_size);

    // Merge constructor
    //ColorCompressed(const DBG_succ &graph,
    //                const std::vector<ColorCompressed> &categories,
    //                const std::vector<std::vector<size_t>> &merge_plan);

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
    typedef std::list<uint32_t> IndexContainer;

    const DBG_succ *graph_;

    std::vector<IndexContainer> index_to_ids_;

    std::unordered_map<std::string, uint32_t> label_to_id_;
    std::vector<std::string> id_to_label_;

    std::pair<IndexContainer::const_iterator, bool> find(Index i, uint32_t id) const;

};

} // namespace annotate

#endif // __ANNOTATE_ROW_COMPRESSED_HPP__

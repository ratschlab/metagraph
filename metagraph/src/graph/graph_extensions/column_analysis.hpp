#ifndef __COLUMN_ANALYSIS_HPP__
#define __COLUMN_ANALYSIS_HPP__

#include "sequence_graph.hpp"
#include "annotate_column_compressed.hpp"
#include "utils/string_utils.hpp"
#include "annotation/annotate.hpp"


template <typename Label = std::string>
class ColumnAnalysis : public SequenceGraph::GraphExtension {
  public:

    ColumnAnalysis() {}

    bool load(const std::string &filename_base) {
        auto filename = utils::remove_suffix(filename_base,
            annotate::kColumnAnnotatorExtension) + annotate::kColumnAnnotatorExtension;
        auto label_encoder = annotate::ColumnCompressed<Label>::load_label_encoder(filename);
        std::ignore = label_encoder;
        return true;
    };

    bool load(const std::vector<std::string> &anno_files) {
        for (const auto &filename_base : anno_files) {
            if(!load(filename_base))
                return false;
        }

        return true;
    };

    void serialize(const std::string &filename_base) const {
        std::ignore = filename_base;
    };

    bool is_compatible(const SequenceGraph &graph, bool verbose = true) const {
        std::ignore = graph;
        std::ignore = verbose;
        return true;
    };

  private:

};

#endif // __COLUMN_ANALYSIS_HPP__

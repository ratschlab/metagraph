#ifndef __COLUMN_ANALYSIS_HPP__
#define __COLUMN_ANALYSIS_HPP__

#include "sequence_graph.hpp"
#include "annotate_column_compressed.hpp"
#include "utils/string_utils.hpp"
#include "annotation/annotate.hpp"
#include "graph/annotated_dbg.hpp"


template <typename Index = uint64_t, typename Label = std::string>
class ColumnAnalysis : public AnnotatedDBG::AnnotatedGraphExtension {
  public:

    ColumnAnalysis() {}

    bool load(const std::string &filename_base) {
        auto filename = utils::remove_suffix(filename_base,
            annotate::kColumnAnnotatorExtension) + annotate::kColumnAnnotatorExtension;

        auto label_encoder = annotate::ColumnCompressed<Label>::load_label_encoder(filename);
        label_encoders_.push_back(std::move(label_encoder));

        return true;
    };

    bool load(const std::vector<std::string> &anno_file_basenames) {
        for (const auto &filename_base : anno_file_basenames) {
            if(!load(filename_base))
                return false;
        }

        return true;
    };

    void serialize(const std::string &filename_base) const {
        std::ignore = filename_base;
    };

    bool is_compatible(const AnnotatedDBG &anno_graph, bool verbose = true) const {
        std::ignore = verbose;

        for (auto &&label_encoder : label_encoders_) {
            for (const Label &label : label_encoder->get_labels()) {
                if (!anno_graph.label_exists(label))
                    return false;
            }
        }

        return true;
    };

  private:
    std::vector<std::unique_ptr<annotate::LabelEncoder<Label>>> label_encoders_;

};

#endif // __COLUMN_ANALYSIS_HPP__

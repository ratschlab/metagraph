#ifndef __COLUMN_ANALYSIS_HPP__
#define __COLUMN_ANALYSIS_HPP__

#include "sequence_graph.hpp"
#include "annotate_column_compressed.hpp"
#include "utils/string_utils.hpp"
#include "annotation/annotate.hpp"
#include "graph/annotated_dbg.hpp"
#include "serialization.hpp"

#include <fstream>


template <typename Index = uint64_t, typename Label = std::string>
class ColumnAnalysis : public AnnotatedDBG::AnnotatedGraphExtension {
  public:

    ColumnAnalysis() {}

    std::unordered_map<Label, double> densities() {
        std::unordered_map<Label, double> result;

        for (const auto &element : column_locations_) {
            const Label &label = element.first;

            std::unique_ptr<bit_vector> column = load_column_by_label(label);

            result.insert({label, column->num_set_bits()
                    / static_cast<double>(get_base_object()->get_graph().num_nodes())});
        }
        return result;
    }

    bool load(const std::string &filename_base) {
        auto filename = utils::remove_suffix(filename_base,
            annotate::kColumnAnnotatorExtension) + annotate::kColumnAnnotatorExtension;

        std::ifstream instream(filename, std::ios::binary);
            if (!instream.good())
                throw std::ifstream::failure("can't open stream");

        const auto num_rows = load_number(instream);

        auto label_encoder_load = std::make_unique<annotate::LabelEncoder<Label>>();
        if (!label_encoder_load->load(instream))
            throw std::ifstream::failure("can't load label encoder");

        if (!label_encoder_load->size()) {
            std::cerr << "No labels in " << filename << "\n" << std::flush;
            throw std::ifstream::failure("label encoder is empty: " + filename);
        }

        for (size_t c = 0; c < label_encoder_load->size(); ++c) {

            auto pos = instream.tellg();

            std::unique_ptr<bit_vector> new_column =
                annotate::ColumnCompressed<Label>::load_column_from_stream(instream);

            if (new_column->size() != num_rows)
                throw std::ifstream::failure("inconsistent column size");

            column_locations_.emplace(label_encoder_load->decode(c),
                std::make_pair(filename, pos));
        }

        label_encoders_.push_back(std::move(label_encoder_load));

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

    std::unique_ptr<bit_vector> load_column_by_label(const Label &label) {
        const auto &[filename, offset] = column_locations_[label];

        std::ifstream instream(filename, std::ios::binary);
        load_number(instream);
        instream.seekg(offset, instream.beg);

        std::unique_ptr<bit_vector> column =
            annotate::ColumnCompressed<Label>::load_column_from_stream(instream);

        return column;
    }

    std::vector<std::unique_ptr<annotate::LabelEncoder<Label>>> label_encoders_;

    // Label -> (columncompressed filename, offset to column)
    std::map<Label, std::pair<std::string, std::streampos>> column_locations_;

};

#endif // __COLUMN_ANALYSIS_HPP__

#ifndef __ANNOTATION_CONVERTERS_HPP__
#define __ANNOTATION_CONVERTERS_HPP__

#include <memory>
#include <vector>
#include <functional>

#include "annotate.hpp"
#include "annotate_row_compressed.hpp"

#include "utils.hpp"


namespace annotate {

typedef LabelEncoder<std::string> LEncoder;


template <typename Label>
class RowCompressed;

template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert(RowCompressed<Label>&& annotation);

template <class FromAnnotation, class ToAnnotation, typename Label, bool sparse = false>
typename std::unique_ptr<ToAnnotation>
convert(const std::string &filename);


//TODO: rather than remove eigenspmat, make it identifiable by file extension?
template <class FromAnnotation, class ToAnnotation, typename Label, bool sparse = false>
uint64_t merge(const std::vector<std::string> &filenames, const std::string &outfile);

template <class FromAnnotation, class ToAnnotation, typename Label, bool sparse = false>
uint64_t merge(const std::vector<const FromAnnotation*> &annotators, const std::string &outfile);


template <typename Label>
class ColumnCompressed;

template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert(ColumnCompressed<Label>&& annotation);

template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert_to_simple_BRWT(ColumnCompressed<Label>&& annotation,
                       size_t grouping_arity = 2,
                       size_t num_threads = 1);

template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert_to_greedy_BRWT(ColumnCompressed<Label>&& annotation,
                       size_t num_threads = 1);

template <class StaticAnnotation>
void relax_BRWT(StaticAnnotation *annotation,
                size_t relax_max_arity,
                size_t num_threads = 1);


uint64_t merge_rows(
    const std::vector<const LEncoder*> &label_encoders,
    const std::function<const std::vector<uint64_t>(const uint64_t, const uint64_t)> &get_row,
    const uint64_t &num_rows,
    const std::function<uint64_t(LEncoder&, const std::function<void (const BinaryMatrix::RowCallback&)>)> &write_rows
);

//TODO: delete other merge functions?
template <class ToAnnotation, typename Label, bool sparse = false>
uint64_t
merge(const std::vector<const MultiLabelEncoded<uint64_t, Label>*> &annotators, const std::vector<std::string> &filenames, const std::string &outfile) {

    uint64_t num_rows;
    uint64_t num_relations;

    if (annotators.size()>0) {
        num_rows = annotators.at(0)->num_objects();
    } else if (filenames.size()>0) {
        RowCompressed<Label>::stream_counts(filenames.at(0), &num_rows, &num_relations, sparse);
    }
    assert(num_rows);

    std::vector<const LEncoder*> label_encoders;
    for (auto annotator : annotators) {
        label_encoders.push_back(&annotator->label_encoder_);
    }

    std::vector<std::unique_ptr<const LEncoder> > loaded_label_encoders;
    std::vector<std::unique_ptr<class RowCompressed<Label>::StreamRows> > streams;
    for (auto filename : filenames) {
        if (utils::ends_with(filename, RowCompressed<Label>::kExtension)) {

            auto label_encoder = RowCompressed<Label>::load_label_encoder(filename);
            label_encoders.push_back(label_encoder.get());
            loaded_label_encoders.push_back(std::move(label_encoder));

            auto annotator = std::make_unique<class RowCompressed<Label>::StreamRows>(filename, sparse);
            streams.push_back(std::move(annotator));

        } else {
            throw std::runtime_error("streaming only supported for rowcompressed annotator");
        }
    }

    merge_rows(
        label_encoders,
        [&](const uint64_t row_idx, const uint64_t annotator_idx) -> const std::vector<uint64_t> {
            //TODO: instead some vector of some sort of union-type elements?
            if (annotator_idx < annotators.size()) {
                return annotators.at(annotator_idx)->get_label_indexes(row_idx);
            } else {
                return *streams[annotator_idx-annotators.size()]->next_row();
            }
        },
        num_rows,
        [&](LEncoder &merged_label_enc, const std::function<void (const BinaryMatrix::RowCallback&)> &callback) {
            //TODO: add write_rows to rowflat which just writes rowcompressed then converts and serializes as below, then can use template param here and delete if statement below
            return RowCompressed<Label>::write_rows(outfile,
                                               merged_label_enc,
                                               callback,
                                               sparse);
        }
    );

    if (!std::is_same<RowCompressed<Label>, ToAnnotation>::value) {
        auto out_annotator = convert<RowCompressed<Label>, ToAnnotation, Label>(outfile);
        out_annotator->serialize(outfile);
    }

    return num_rows;
}

} // namespace annotate

#endif // __ANNOTATION_CONVERTERS_HPP__

#ifndef __TEST_ANNOTATION_HPP__
#define __TEST_ANNOTATION_HPP__

#include "gtest/gtest.h"

#include "../test_helpers.hpp"

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/annotation_converters.hpp"


namespace mtg {
namespace test {

template <typename... Args>
class RowCompressedDynamic : public annot::RowCompressed<Args...> {
  public:
    template <typename... CArgs>
    RowCompressedDynamic(CArgs&&... args)
          : annot::RowCompressed<Args...>(std::forward<CArgs>(args)...) {}
};

template <typename... Args>
class RowCompressedSparse : public annot::RowCompressed<Args...> {
  public:
    RowCompressedSparse(uint64_t num_rows = 0)
          : annot::RowCompressed<Args...>(num_rows, true) {}
};

template <typename Annotator>
class AnnotatorTest : public ::testing::Test {
  public:
    std::unique_ptr<Annotator> annotation;
    const std::string test_data_dir = "../tests/data";
    const std::string test_dump_basename = test_data_dir + "/bit_vector_dump_test";

    virtual void set(annot::ColumnCompressed<>&& column_annotator) {
        if constexpr(std::is_same_v<Annotator, annot::MultiBRWTAnnotator>) {
            annotation = annot::convert_to_simple_BRWT(std::move(column_annotator));

        } else if constexpr(std::is_same_v<Annotator, annot::RowCompressed<>>) {
            convert_to_row_annotator(column_annotator, test_dump_basename);
            annotation.reset(new annot::RowCompressed<>(0));
            annotation->load(test_dump_basename);

        } else if constexpr(std::is_same_v<Annotator, RowCompressedDynamic<>>) {
            annotation.reset(new RowCompressedDynamic<>(column_annotator.num_objects()));
            for (RowCompressedDynamic<>::Index i = 0; i < column_annotator.num_objects(); ++i) {
                annotation->add_labels({ i }, std::move(column_annotator.get_labels(i)));
            }

        } else if constexpr(std::is_same_v<Annotator, RowCompressedSparse<>>) {
            annotation.reset(new RowCompressedSparse<>(column_annotator.num_objects()));
            for (RowCompressedSparse<>::Index i = 0; i < column_annotator.num_objects(); ++i) {
                annotation->add_labels({ i }, std::move(column_annotator.get_labels(i)));
            }

        } else if constexpr(std::is_same_v<Annotator, annot::ColumnCompressed<>>) {
            // TODO: introduce move constructor for ColumnCompressed
            //annotation.reset(new annot::ColumnCompressed<>(std::move(column_annotator)));
            annotation.reset(new annot::ColumnCompressed<>(column_annotator.num_objects()));
            for (annot::ColumnCompressed<>::Index i = 0; i < column_annotator.num_objects(); ++i) {
                annotation->add_labels({ i }, std::move(column_annotator.get_labels(i)));
            }
        } else {
            annotation = annot::convert<Annotator>(std::move(column_annotator));
        }
    }

    virtual void SetUp() { set(annot::ColumnCompressed<>(0)); }
};

template <typename Annotator>
class AnnotatorStaticTest : public AnnotatorTest<Annotator> { };

template <typename Annotator>
class AnnotatorPresetTest : public AnnotatorTest<Annotator> {
  public:
    virtual void SetUp() override {
        annot::ColumnCompressed<> column_annotator(5);
        column_annotator.add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        column_annotator.add_labels({ 2 }, {"Label1", "Label2"});
        column_annotator.add_labels({ 3 }, {"Label1", "Label2", "Label8"});
        column_annotator.add_labels({ 4 }, {"Label2"});
        this->set(std::move(column_annotator));
    }
};

template <typename Annotator>
class AnnotatorPreset2Test : public AnnotatorTest<Annotator> {
  public:
    virtual void SetUp() override {
        annot::ColumnCompressed<> column_annotator(5);
        column_annotator.add_labels({ 0 }, { "Label0", "Label2", "Label8" });
        column_annotator.add_labels({ 2 }, { "Label1", "Label2" });
        column_annotator.add_labels({ 4 }, { "Label8" });
        this->set(std::move(column_annotator));
    }
};

template <typename Annotator>
class AnnotatorPreset3Test : public AnnotatorTest<Annotator> {
  public:
    virtual void SetUp() override {
        annot::ColumnCompressed<> column_annotator(5);
        column_annotator.add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        column_annotator.add_labels({ 2 }, {"Label1", "Label2"});
        column_annotator.add_labels({ 3 }, {"Label1", "Label2", "Label8"});
        column_annotator.add_labels({ 4 }, {"Label2", "Label8"});
        this->set(std::move(column_annotator));
    }
};

template <typename Annotator>
class AnnotatorDynamicTest : public AnnotatorPreset2Test<Annotator> { };

template <typename Annotator>
class AnnotatorDynamicNoSparseTest : public AnnotatorPreset2Test<Annotator> { };


typedef ::testing::Types<annot::BinRelWTAnnotator,
                         annot::RbBRWTAnnotator,
                         annot::MultiBRWTAnnotator,
                         annot::RainbowfishAnnotator,
                         annot::RowFlatAnnotator,
                         annot::RowSparseAnnotator,
                         annot::UniqueRowAnnotator,
                         annot::ColumnCompressed<>,
                         annot::RowCompressed<>,
                         RowCompressedDynamic<>,
                         RowCompressedSparse<>> AnnotatorTypes;

typedef ::testing::Types<annot::BinRelWTAnnotator,
                         annot::RbBRWTAnnotator,
                         annot::RainbowfishAnnotator,
                         annot::RowFlatAnnotator,
                         annot::RowSparseAnnotator,
                         annot::UniqueRowAnnotator,
                         annot::MultiBRWTAnnotator> AnnotatorStaticTypes;

typedef ::testing::Types<annot::ColumnCompressed<>,
                         annot::RowCompressed<>,
                         RowCompressedDynamic<>,
                         RowCompressedSparse<>> AnnotatorDynamicTypes;

typedef ::testing::Types<annot::ColumnCompressed<>,
                         annot::RowCompressed<>,
                         RowCompressedDynamic<>> AnnotatorDynamicNoSparseTypes;


TYPED_TEST_SUITE(AnnotatorTest, AnnotatorTypes);
TYPED_TEST_SUITE(AnnotatorPresetTest, AnnotatorTypes);
TYPED_TEST_SUITE(AnnotatorPreset2Test, AnnotatorTypes);
TYPED_TEST_SUITE(AnnotatorPreset3Test, AnnotatorTypes);
TYPED_TEST_SUITE(AnnotatorStaticTest, AnnotatorStaticTypes);
TYPED_TEST_SUITE(AnnotatorDynamicTest, AnnotatorDynamicTypes);
TYPED_TEST_SUITE(AnnotatorDynamicNoSparseTest, AnnotatorDynamicNoSparseTypes);

} // namespace test
} // namespace mtg

#endif // __TEST_ANNOTATION_HPP__

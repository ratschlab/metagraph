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
class RowCompressedParallel : public annotate::RowCompressed<Args...> {
  public:
    template <typename... CArgs>
    RowCompressedParallel(CArgs&&... args)
          : annotate::RowCompressed<Args...>(std::forward<CArgs>(args)...) {}
};

template <typename... Args>
class RowCompressedDynamic : public annotate::RowCompressed<Args...> {
  public:
    template <typename... CArgs>
    RowCompressedDynamic(CArgs&&... args)
          : annotate::RowCompressed<Args...>(std::forward<CArgs>(args)...) {}
};

template <typename... Args>
class RowCompressedSparse : public annotate::RowCompressed<Args...> {
  public:
    RowCompressedSparse(uint64_t num_rows = 0)
          : annotate::RowCompressed<Args...>(num_rows, true) {}
};

template <typename Annotator>
class AnnotatorTest : public ::testing::Test {
  public:
    std::unique_ptr<Annotator> annotation;

    virtual void set(annotate::ColumnCompressed<>&& column_annotator) {
        if constexpr(std::is_same_v<Annotator, annotate::MultiBRWTAnnotator>) {
            annotation = annotate::convert_to_simple_BRWT<annotate::MultiBRWTAnnotator>(
                std::move(column_annotator)
            );

        } else if constexpr(std::is_same_v<Annotator, annotate::RowCompressed<>>) {
            annotation.reset(new annotate::RowCompressed<>(column_annotator.num_objects()));
            convert_to_row_annotator(column_annotator, annotation.get());

        } else if constexpr(std::is_same_v<Annotator, RowCompressedParallel<>>) {
            annotation.reset(new RowCompressedParallel<>(column_annotator.num_objects()));
            convert_to_row_annotator(column_annotator, annotation.get(), 10);

        } else if constexpr(std::is_same_v<Annotator, RowCompressedDynamic<>>) {
            annotation.reset(new RowCompressedDynamic<>(column_annotator.num_objects()));
            for (RowCompressedDynamic<>::Index i = 0; i < column_annotator.num_objects(); ++i) {
                annotation->add_labels({ i }, std::move(column_annotator.get(i)));
            }

        } else if constexpr(std::is_same_v<Annotator, RowCompressedSparse<>>) {
            annotation.reset(new RowCompressedSparse<>(column_annotator.num_objects()));
            for (RowCompressedSparse<>::Index i = 0; i < column_annotator.num_objects(); ++i) {
                annotation->add_labels({ i }, std::move(column_annotator.get(i)));
            }

        } else if constexpr(std::is_same_v<Annotator, annotate::ColumnCompressed<>>) {
            // TODO: introduce move constructor for ColumnCompressed
            //annotation.reset(new annotate::ColumnCompressed<>(std::move(column_annotator)));
            annotation.reset(new annotate::ColumnCompressed<>(column_annotator.num_objects()));
            for (annotate::ColumnCompressed<>::Index i = 0; i < column_annotator.num_objects(); ++i) {
                annotation->add_labels({ i }, std::move(column_annotator.get(i)));
            }
        } else {
            annotation = annotate::convert<Annotator>(std::move(column_annotator));
        }
    }

    virtual void SetUp() { set(annotate::ColumnCompressed<>(0)); }
};

template <typename Annotator>
class AnnotatorStaticTest : public AnnotatorTest<Annotator> { };

template <typename Annotator>
class AnnotatorStaticLargeTest : public AnnotatorStaticTest<Annotator> { };

template <typename Annotator>
class AnnotatorPresetTest : public AnnotatorTest<Annotator> {
  public:
    virtual void SetUp() override {
        annotate::ColumnCompressed<> column_annotator(5);
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
        annotate::ColumnCompressed<> column_annotator(5);
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
        annotate::ColumnCompressed<> column_annotator(5);
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


typedef ::testing::Types<annotate::BinRelWTAnnotator,
                         annotate::BinRelWT_sdslAnnotator,
                         annotate::MultiBRWTAnnotator,
                         annotate::RainbowfishAnnotator,
                         annotate::RowFlatAnnotator,
                         annotate::UniqueRowAnnotator,
                         annotate::ColumnCompressed<>,
                         annotate::RowCompressed<>,
                         RowCompressedParallel<>,
                         RowCompressedDynamic<>,
                         RowCompressedSparse<>> AnnotatorTypes;
typedef ::testing::Types<annotate::BinRelWTAnnotator,
                         annotate::BinRelWT_sdslAnnotator,
                         annotate::RainbowfishAnnotator,
                         annotate::RowFlatAnnotator,
                         annotate::UniqueRowAnnotator,
                         annotate::MultiBRWTAnnotator> AnnotatorStaticTypes;
typedef ::testing::Types<annotate::MultiBRWTAnnotator> AnnotatorStaticLargeTypes;
typedef ::testing::Types<annotate::ColumnCompressed<>,
                         annotate::RowCompressed<>,
                         RowCompressedParallel<>,
                         RowCompressedDynamic<>,
                         RowCompressedSparse<>> AnnotatorDynamicTypes;
typedef ::testing::Types<annotate::ColumnCompressed<>,
                         annotate::RowCompressed<>,
                         RowCompressedParallel<>,
                         RowCompressedDynamic<>> AnnotatorDynamicNoSparseTypes;


TYPED_TEST_SUITE(AnnotatorTest, AnnotatorTypes);
TYPED_TEST_SUITE(AnnotatorPresetTest, AnnotatorTypes);
TYPED_TEST_SUITE(AnnotatorPreset2Test, AnnotatorTypes);
TYPED_TEST_SUITE(AnnotatorPreset3Test, AnnotatorTypes);
TYPED_TEST_SUITE(AnnotatorStaticTest, AnnotatorStaticTypes);
TYPED_TEST_SUITE(AnnotatorStaticLargeTest, AnnotatorStaticLargeTypes);
TYPED_TEST_SUITE(AnnotatorDynamicTest, AnnotatorDynamicTypes);
TYPED_TEST_SUITE(AnnotatorDynamicNoSparseTest, AnnotatorDynamicNoSparseTypes);

} // namespace test
} // namespace mtg

#endif // __TEST_ANNOTATION_HPP__

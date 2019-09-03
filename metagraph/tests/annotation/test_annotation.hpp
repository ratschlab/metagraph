#ifndef __TEST_ANNOTATION_HPP__
#define __TEST_ANNOTATION_HPP__

#include "gtest/gtest.h"

#include "../test_helpers.hpp"

#include "annotate_column_compressed.hpp"
#include "annotate_row_compressed.hpp"
#include "static_annotators_def.hpp"
#include "annotation_converters.hpp"


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
        if constexpr(std::is_same_v<Annotator, annotate::BRWTCompressed<>>) {
            annotation = annotate::convert_to_simple_BRWT<annotate::BRWTCompressed<>>(
                std::move(column_annotator)
            );

        } else if constexpr(std::is_same_v<Annotator, annotate::RowCompressed<>>) {
            annotation.reset(new annotate::RowCompressed<>(0));
            column_annotator.convert_to_row_annotator(annotation.get());

        } else if constexpr(std::is_same_v<Annotator, RowCompressedParallel<>>) {
            annotation.reset(new RowCompressedParallel<>(0));
            column_annotator.convert_to_row_annotator(annotation.get(), 10);

        } else if constexpr(std::is_same_v<Annotator, RowCompressedDynamic<>>) {
            annotation.reset(new RowCompressedDynamic<>(column_annotator.num_objects()));
            for (RowCompressedDynamic<>::Index i = 0; i < column_annotator.num_objects(); ++i) {
                annotation->add_labels(i, std::move(column_annotator.get_labels(i)));
            }

        } else if constexpr(std::is_same_v<Annotator, RowCompressedSparse<>>) {
            annotation.reset(new RowCompressedSparse<>(column_annotator.num_objects()));
            for (RowCompressedSparse<>::Index i = 0; i < column_annotator.num_objects(); ++i) {
                annotation->add_labels(i, std::move(column_annotator.get_labels(i)));
            }

        } else if constexpr(std::is_same_v<Annotator, annotate::ColumnCompressed<>>) {
            // TODO: introduce move constructor for ColumnCompressed
            //annotation.reset(new annotate::ColumnCompressed<>(std::move(column_annotator)));
            annotation.reset(new annotate::ColumnCompressed<>(column_annotator.num_objects()));
            for (annotate::ColumnCompressed<>::Index i = 0; i < column_annotator.num_objects(); ++i) {
                annotation->add_labels(i, std::move(column_annotator.get_labels(i)));
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
        column_annotator.add_labels(0, {"Label0", "Label2", "Label8"});
        column_annotator.add_labels(2, {"Label1", "Label2"});
        column_annotator.add_labels(3, {"Label1", "Label2", "Label8"});
        column_annotator.add_labels(4, {"Label2"});
        this->set(std::move(column_annotator));
    }
};

template <typename Annotator>
class AnnotatorPreset2Test : public AnnotatorTest<Annotator> {
  public:
    virtual void SetUp() override {
        annotate::ColumnCompressed<> column_annotator(5);
        column_annotator.add_labels(0, { "Label0", "Label2", "Label8" });
        column_annotator.add_labels(2, { "Label1", "Label2" });
        column_annotator.add_labels(4, { "Label8" });
        this->set(std::move(column_annotator));
    }
};

template <typename Annotator>
class AnnotatorPreset3Test : public AnnotatorTest<Annotator> {
  public:
    virtual void SetUp() override {
        annotate::ColumnCompressed<> column_annotator(5);
        column_annotator.add_labels(0, {"Label0", "Label2", "Label8"});
        column_annotator.add_labels(2, {"Label1", "Label2"});
        column_annotator.add_labels(3, {"Label1", "Label2", "Label8"});
        column_annotator.add_labels(4, {"Label2", "Label8"});
        this->set(std::move(column_annotator));
    }
};

template <typename Annotator>
class AnnotatorPresetDumpTest : public AnnotatorPreset2Test<Annotator> { };

template <typename Annotator>
class AnnotatorDynamicTest : public AnnotatorPreset2Test<Annotator> { };

template <typename Annotator>
class AnnotatorDynamicNoSparseTest : public AnnotatorPreset2Test<Annotator> { };


typedef ::testing::Types<annotate::BinRelWTAnnotator,
                         annotate::BinRelWT_sdslAnnotator,
                         annotate::BRWTCompressed<>,
                         annotate::RainbowfishAnnotator,
                         annotate::RowFlatAnnotator,
                         annotate::ColumnCompressed<>,
                         annotate::RowCompressed<>,
                         RowCompressedParallel<>,
                         RowCompressedDynamic<>,
                         RowCompressedSparse<>> AnnotatorTypes;
typedef ::testing::Types<annotate::BinRelWTAnnotator,
                         annotate::BinRelWT_sdslAnnotator,
                         annotate::RainbowfishAnnotator,
                         annotate::RowFlatAnnotator,
                         annotate::BRWTCompressed<>> AnnotatorStaticTypes;
typedef ::testing::Types<annotate::BRWTCompressed<>> AnnotatorStaticLargeTypes;
typedef ::testing::Types<annotate::ColumnCompressed<>,
                         annotate::RowCompressed<>,
                         RowCompressedParallel<>,
                         RowCompressedDynamic<>,
                         RowCompressedSparse<>> AnnotatorDynamicTypes;
typedef ::testing::Types<annotate::ColumnCompressed<>,
                         annotate::RowCompressed<>,
                         RowCompressedParallel<>,
                         RowCompressedDynamic<>> AnnotatorDynamicNoSparseTypes;

typedef ::testing::Types<annotate::BinRelWTAnnotator,
                         annotate::BinRelWT_sdslAnnotator,
                         annotate::BRWTCompressed<>,
                         annotate::RainbowfishAnnotator,
                         annotate::RowFlatAnnotator,
                         annotate::ColumnCompressed<>> AnnotatorDumpTestTypes;


TYPED_TEST_CASE(AnnotatorTest, AnnotatorTypes);
TYPED_TEST_CASE(AnnotatorPresetTest, AnnotatorTypes);
TYPED_TEST_CASE(AnnotatorPreset2Test, AnnotatorTypes);
TYPED_TEST_CASE(AnnotatorPreset3Test, AnnotatorTypes);
TYPED_TEST_CASE(AnnotatorStaticTest, AnnotatorStaticTypes);
TYPED_TEST_CASE(AnnotatorStaticLargeTest, AnnotatorStaticLargeTypes);
TYPED_TEST_CASE(AnnotatorDynamicTest, AnnotatorDynamicTypes);
TYPED_TEST_CASE(AnnotatorDynamicNoSparseTest, AnnotatorDynamicNoSparseTypes);
TYPED_TEST_CASE(AnnotatorPresetDumpTest, AnnotatorDumpTestTypes);


std::set<std::string> convert_to_set(const std::vector<std::string> &vector);

#endif // __TEST_ANNOTATION_HPP__

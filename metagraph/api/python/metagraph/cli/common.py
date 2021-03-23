from enum import Enum


class AnnotationLabelsSource(Enum):
    SEQUENCE_HEADERS = 'sequence_headers'
    SEQUENCE_FILE_NAMES = 'sequence_file_names'

    def to_annotation_cmd_option(self):
        # TODO: check comparison
        if self == self.SEQUENCE_FILE_NAMES:
            return '--anno-filename'
        return '--anno-header'


class AnnotationFormats(Enum):
    # COLUMN = 'column' # TODO: need special case in the workflow
    ROW = 'row'
    BIN_REL_WT_SDSL = 'bin_rel_wt_sdsl' # TODO: friendlier name?
    BIN_REL_WT = 'bin_rel_wt'
    FLAT = 'flat'
    RBFISH = 'rbfish'
    BRWT = 'brwt'
    RELAXED_BRWT = 'relax.brwt'
    RB_BRWT = 'rb_brwt'
    #RELAXED_RB_BRWT = 'relax.rb_brwt' # not possible
    ROW_DIFF_BRWT = 'row_diff_brwt'
    RELAXED_ROW_DIFF_BRWT = 'relax.row_diff_brwt'

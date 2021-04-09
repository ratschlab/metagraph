import unittest
import time
import subprocess
import os


"""Formatting results"""


class TimeLoggingTestResult(unittest.TextTestResult):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dots = False
        self.__total_time = 0
        self.__num_successes = 0
        self.stream.write("\033[0;32;40m[----------]\033[0m {}\n".format("Run tests..."))

    def stopTestRun(self):
        self.stream.write("\033[0;32;40m[----------]\033[0m {}\n".format(""))
        self.stream.write("\033[0;32;40m[==========]\033[0m {}\n".format(
            "{} tests ran. ({:.2f} sec total)".format(self.testsRun, self.__total_time)
        ))

        self.stream.write("\033[0;32;40m[  PASSED  ]\033[0m {}\n".format(
            "{} tests.".format(self.__num_successes)
        ))

        # Report failures
        if len(self.failures):
            self.stream.write("\033[0;31;40m[  FAILED  ]\033[0m {} test(s), listed below:\n".format(
                len(self.failures)
            ))

        for failure in self.failures:
            self.stream.write("\033[0;31;40m[  FAILED  ]\033[0m {}\n".format(failure[0]))

        # Report errors
        if len(self.errors):
            self.stream.write("\033[0;31;40m[  ERRORS  ]\033[0m {} test(s), listed below:\n".format(
                len(self.errors)
            ))

        for error in self.errors:
            self.stream.write("\033[0;31;40m[  ERRORS  ]\033[0m {}\n".format(error[0]))

        super().stopTestRun()

    def startTest(self, test):
        name = self.getDescription(test)
        self.stream.write("\033[0;32;40m[ RUN      ]\033[0m {}\n".format(name))
        self._started_at = time.time()
        super().startTest(test)

    def addSuccess(self, test):
        elapsed = time.time() - self._started_at
        self.__total_time += elapsed
        self.__num_successes += 1
        name = self.getDescription(test)
        self.stream.write("\033[0;32;40m[       OK ]\033[0m {} ({:.03} sec)\n".format(name, elapsed))
        super().addSuccess(test)

    def addFailure(self, test, err):
        elapsed = time.time() - self._started_at
        self.__total_time += elapsed
        name = self.getDescription(test)
        self.stream.write("\033[0;31;40m[   FAIL   ]\033[0m {}\n".format(name))
        super().addFailure(test, err)
        self.stream.write("{}\n".format(self.failures[-1][1]))

    def addError(self, test, err):
        elapsed = time.time() - self._started_at
        self.__total_time += elapsed
        name = self.getDescription(test)
        self.stream.write("\033[0;31;40m[   ERROR  ]\033[0m {}\n".format(name))
        super().addError(test, err)
        self.stream.write("{}\n".format(self.errors[-1][1]))

    def printErrors(self):
        # all errors are printed in addFailure and addError
        pass


def get_test_class_name(cls, num, params_dict):
    # By default the generated class named includes either the "name"
    # parameter (if present), or the first string value.
    return "{}_{}_{}".format(cls.__name__, num, '_'.join(params_dict.values()))

graph_file_extension = {'succinct': '.dbg',
                        'bitmap': '.bitmapdbg',
                        'hash': '.orhashdbg',
                        'hashfast': '.hashfastdbg',
                        'hashstr': '.hashstrdbg'}

anno_file_extension = {'column': '.column.annodbg',
                       'row': '.row.annodbg',
                       'row_diff': '.row_diff.annodbg',
                       'row_sparse': '.row_sparse.annodbg',
                       'row_diff_brwt': '.row_diff_brwt.annodbg',
                       'row_diff_sparse': '.row_diff_sparse.annodbg',
                       'rb_brwt': '.rb_brwt.annodbg',
                       'brwt': '.brwt.annodbg',
                       'rbfish': '.rbfish.annodbg',
                       'flat': '.flat.annodbg'}

GRAPH_TYPES = [graph_type for graph_type, _ in graph_file_extension.items()]
ANNO_TYPES = [anno_type for anno_type, _ in anno_file_extension.items()]

def product(graph_types, anno_types):
    result  = []
    for graph in graph_types:
        for anno in anno_types:
            if graph == 'succinct' or not anno.startswith('row_diff'):
                result.append((graph, anno))
    return result

METAGRAPH = './metagraph'
NUM_THREADS = 4
def build_annotation(graph_filename, input_fasta, anno_repr, output_filename,
                     separate=False, no_fork_opt=False, no_anchor_opt=False):
    target_anno = anno_repr
    if anno_repr in {'rb_brwt', 'brwt', 'row_sparse'} or anno_repr.startswith('row_diff'):
        target_anno = anno_repr
        anno_repr = 'column'
    elif anno_repr in {'flat', 'rbfish'}:
        target_anno = anno_repr
        anno_repr = 'row'

    annotate_command = '{exe} annotate -p {num_threads} --anno-header -i {graph} \
            --anno-type {anno_repr} -o {outfile} {input}'.format(
        exe=METAGRAPH,
        num_threads=NUM_THREADS,
        graph=graph_filename,
        anno_repr=anno_repr,
        outfile=output_filename,
        input=input_fasta
    )
    res = subprocess.run([annotate_command], shell=True)
    assert(res.returncode == 0)

    if target_anno == anno_repr:
        return

    final_anno = target_anno
    if final_anno.startswith('row_diff'):
        target_anno = 'row_diff'

    annotate_command = '{exe} transform_anno -p {num_threads} \
            --anno-type {target_anno} -o {outfile} {input}'.format(
        exe=METAGRAPH,
        num_threads=NUM_THREADS,
        graph=graph_filename,
        target_anno=target_anno,
        outfile=output_filename,
        input=output_filename + anno_file_extension[anno_repr]
    )
    if target_anno == 'row_diff':
        annotate_command += ' -i ' + graph_filename

    if not no_fork_opt:
        print('-- Building RowDiff without fork optimization...')
        res = subprocess.run([annotate_command], shell=True)
        assert(res.returncode == 0)

    if target_anno == 'row_diff':
        without_input_anno = annotate_command.split(' ')
        without_input_anno.pop(-3)
        without_input_anno = ' '.join(without_input_anno)
        if not no_anchor_opt:
            if separate:
                print('-- Building RowDiff succ/pred...')
                res = subprocess.run(['echo \"\" | ' + without_input_anno + ' --row-diff-stage 1'], shell=True)
                assert(res.returncode == 0)
            res = subprocess.run([annotate_command + ' --row-diff-stage 1'], shell=True)
            assert(res.returncode == 0)
        if separate:
            print('-- Assigning anchors...')
            res = subprocess.run(['echo \"\" | ' + without_input_anno + ' --row-diff-stage 2'], shell=True)
            assert(res.returncode == 0)
        res = subprocess.run([annotate_command + ' --row-diff-stage 2'], shell=True)
        assert(res.returncode == 0)

        os.remove(output_filename + anno_file_extension[anno_repr])

        if final_anno != target_anno:
            annotate_command = f'{METAGRAPH} transform_anno --anno-type {final_anno} --greedy -o {output_filename} ' \
                               f'-i {graph_filename} -p {NUM_THREADS} {output_filename}.row_diff.annodbg'
            res = subprocess.run([annotate_command], shell=True)
            assert (res.returncode == 0)
            os.remove(output_filename + anno_file_extension['row_diff'])

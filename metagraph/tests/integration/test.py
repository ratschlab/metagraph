import glob
import unittest
import os
import time


"""Run all integration tests"""


def create_test_suite():
    test_file_strings = glob.glob(os.path.dirname(os.path.realpath(__file__)) + '/test_*.py')
    module_strings = [os.path.basename(test_filename)[:-3] for test_filename in test_file_strings]
    suites = [unittest.defaultTestLoader.loadTestsFromName(name) for name in module_strings]
    return unittest.TestSuite(suites)


class TimeLoggingTestResult(unittest.TextTestResult):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.stream.write("\033[0;32;40m[----------]\033[0m {}\n".format(""))

    def stopTestRun(self):
        self.stream.write("\033[0;32;40m[----------]\033[0m {}\n".format(""))
        super().stopTestRun()

    def startTest(self, test):
        name = self.getDescription(test)
        self.stream.write("\033[0;32;40m[ RUN      ]\033[0m {}\n".format(name))
        self._started_at = time.time()
        super().startTest(test)

    def addSuccess(self, test):
        elapsed = time.time() - self._started_at
        name = self.getDescription(test)
        self.stream.write("\033[0;32;40m[       OK ]\033[0m {} ({:.03} s)\n".format(name, elapsed))
        # super().addSuccess(test)

if __name__ == '__main__':
    text_runner = unittest.TextTestRunner(resultclass=TimeLoggingTestResult).run(create_test_suite())
    # unittest.main()

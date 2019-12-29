import glob
import unittest
import os
from helpers import TimeLoggingTestResult


"""Run all integration tests"""


def create_test_suite():
    test_file_strings = glob.glob(os.path.dirname(os.path.realpath(__file__)) + '/test_*.py')
    module_strings = [os.path.basename(test_filename)[:-3] for test_filename in test_file_strings]
    suites = [unittest.defaultTestLoader.loadTestsFromName(name) for name in module_strings]
    return unittest.TestSuite(suites)


if __name__ == '__main__':
    result = unittest.TextTestRunner(
        resultclass=TimeLoggingTestResult
    ).run(create_test_suite())

    exit(0 if result.wasSuccessful() else 1)

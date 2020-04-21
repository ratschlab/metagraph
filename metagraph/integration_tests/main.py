import glob
import unittest
import os
import fnmatch
import sys
import argparse
from helpers import TimeLoggingTestResult


"""Run all integration tests"""


def create_test_suite(filter_pattern="*"):
    test_file_strings = glob.glob(os.path.dirname(os.path.realpath(__file__)) + '/test_*.py')
    module_strings = [os.path.basename(test_filename)[:-3] for test_filename in test_file_strings]

    test_suite = unittest.TestSuite()
    for name in module_strings:
        for suite in unittest.defaultTestLoader.loadTestsFromName(name):
            try:
                for test in suite:
                    if fnmatch.fnmatchcase(test.__str__(), filter_pattern):
                        test_suite.addTest(test)
            except:
                test_suite.addTest(suite)

    return test_suite


if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser(description='Metagraph integration tests.')
        parser.add_argument('--test_filter', dest='filter', type=str, default="*",
                            help='filter test cases (default: run all)')
        args = parser.parse_args()

        result = unittest.TextTestRunner(
            resultclass=TimeLoggingTestResult
        ).run(create_test_suite(args.filter))

        exit(0 if result.wasSuccessful() else 1)

    except KeyboardInterrupt:
        sys.exit(130)

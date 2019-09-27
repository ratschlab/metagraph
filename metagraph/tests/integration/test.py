import glob
import unittest
import os


"""Run all integration tests"""


def create_test_suite():
    test_file_strings = glob.glob(os.path.dirname(os.path.realpath(__file__)) + '/test_*.py')
    module_strings = [os.path.basename(test_filename)[:-3] for test_filename in test_file_strings]
    suites = [unittest.defaultTestLoader.loadTestsFromName(name) for name in module_strings]
    return unittest.TestSuite(suites)


if __name__ == '__main__':
    text_runner = unittest.TextTestRunner().run(create_test_suite())
    # unittest.main()

import unittest
import time


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

        if len(self.failures):
            self.stream.write("\033[0;31;40m[  FAILED  ]\033[0m {} test(s), listed below:\n".format(
                len(self.failures)
            ))

        for failure in self.failures:
            self.stream.write("\033[0;31;40m[  FAILED  ]\033[0m {}\n".format(failure[0]))

        # All tracebacks have been printed in `addFailure`, so we don't need these
        for i in range(len(self.failures)):
            self.failures[i] = (self.failures[i][0], '')

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

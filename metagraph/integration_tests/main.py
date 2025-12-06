import glob
import unittest
import os
import fnmatch
import sys
import argparse
import concurrent.futures
import time
import signal
from helpers import TimeLoggingTestResult
from base import update_prefix


"""Run all integration tests"""


def matches_filter(test, filter_pattern):
    """Check if a test matches the filter pattern"""
    test_str = test.__str__()
    method_name = getattr(test, '_testMethodName', '')
    return (fnmatch.fnmatchcase(test_str, filter_pattern) or
            fnmatch.fnmatchcase(method_name, filter_pattern))

def load_tests_from_module(module_name, filter_pattern="*"):
    all_tests = []
    for suite in unittest.defaultTestLoader.loadTestsFromName(module_name):
        for test in suite:
            if matches_filter(test, filter_pattern):
                all_tests.append(test)
    return all_tests

def create_test_suite(filter_pattern="*"):
    """Create a test suite with filtered tests"""
    test_files = glob.glob(os.path.dirname(os.path.realpath(__file__)) + '/test_*.py')
    module_names = [os.path.basename(f)[:-3] for f in test_files]
    all_tests = []
    for module_name in module_names:
        all_tests += load_tests_from_module(module_name, filter_pattern)
    return unittest.TestSuite(all_tests)


def run_test_parallel(test_info, filter_pattern="*"):
    """Run a test file or chunk in parallel"""
    # Parse test info (either file path or chunk dict with indices)
    is_chunk = isinstance(test_info, dict)
    test_file = test_info['file'] if is_chunk else test_info
    module_name = os.path.basename(test_file)[:-3]

    if is_chunk:
        test_name = f"{module_name}_chunk_{test_info['chunk_id']}"
        start_idx, end_idx = test_info['start_idx'], test_info['end_idx']
    else:
        test_name = module_name

    # Load and filter tests
    all_tests = load_tests_from_module(module_name, filter_pattern)
    # Add tests (slice if chunk, all if file)
    tests_to_run = all_tests[start_idx:end_idx] if is_chunk else all_tests
    test_suite = unittest.TestSuite(tests_to_run)

    # Run with real-time output
    start_time = time.time()

    class RealTimeStream:
        def __init__(self): self.buffer = ""
        def write(self, text): print(text, end=''); self.buffer += text
        def flush(self): pass
        def getvalue(self): return self.buffer

    output_stream = RealTimeStream()
    result = unittest.TextTestRunner(
        resultclass=TimeLoggingTestResult,
        stream=output_stream,
        verbosity=1
    ).run(test_suite)

    duration = time.time() - start_time
    output_text = output_stream.getvalue()

    # Handle expected failures
    has_expected_failures = "expected failures=" in output_text and "OK" in output_text
    success = result.wasSuccessful() or has_expected_failures

    return {
        'test': test_name,
        'success': success,
        'tests_run': result.testsRun,
        'duration': duration,
        'output': output_text
    }

def run_tests_parallel(max_workers, filter_pattern="*"):
    """Run test files in parallel with chunking"""
    # Find all test files
    test_files = glob.glob(os.path.dirname(os.path.realpath(__file__)) + '/test_*.py')

    if not test_files:
        print("No test files found!")
        return False

    # Create chunks from test files, respecting class boundaries
    # Note: test_api uses server ports, so don't chunk it to avoid port collisions
    chunk_size = 20
    all_chunks = []

    for test_file in test_files:
        module_name = os.path.basename(test_file)[:-3]

        # Don't chunk test_api to avoid port collisions
        if module_name == 'test_api':
            all_chunks.append(test_file)
            continue

        # Load all tests and group by class
        all_tests = []
        class_ranges = {}
        try:
            for suite in unittest.defaultTestLoader.loadTestsFromName(module_name):
                for test in suite:
                    if matches_filter(test, filter_pattern):
                        class_name = test.__class__.__name__
                        if class_name not in class_ranges:
                            class_ranges[class_name] = len(all_tests)
                        all_tests.append(test)
        except Exception as e:
            print(f"Warning: Failed to load tests from {module_name}: {e}")
            continue

        if not all_tests:
            continue

        # Convert class_ranges to (start, count) format
        class_info = []
        class_names = list(class_ranges.keys())
        for i, class_name in enumerate(class_names):
            start_idx = class_ranges[class_name]
            end_idx = class_ranges[class_names[i + 1]] if i + 1 < len(class_names) else len(all_tests)
            count = end_idx - start_idx
            class_info.append((class_name, start_idx, count))

        # Split large classes into chunks
        for class_name, start_idx, count in class_info:
            for i in range(0, count, chunk_size):
                all_chunks.append({
                    'file': test_file,
                    'start_idx': start_idx + i,
                    'end_idx': min(start_idx + i + chunk_size, start_idx + count),
                    'chunk_id': len(all_chunks)
                })

    print(f"Running {len(all_chunks)} test chunks in parallel (num_workers={max_workers})")
    print(f"Filter pattern: {filter_pattern}")
    print("-" * 60)

    # Set up signal handling for graceful shutdown
    executor = None
    def signal_handler(signum, frame):
        print("\n\nInterrupted! Shutting down workers...")
        if executor:
            executor.shutdown(wait=False, cancel_futures=True)
        sys.exit(1)

    original_sigint = signal.signal(signal.SIGINT, signal_handler)
    original_sigterm = signal.signal(signal.SIGTERM, signal_handler)

    try:
        start_time = time.time()
        results = []

        executor = concurrent.futures.ProcessPoolExecutor(max_workers=max_workers)
        with executor:
            future_to_chunk = {executor.submit(run_test_parallel, chunk, filter_pattern): chunk
                              for chunk in all_chunks}

            for i, future in enumerate(concurrent.futures.as_completed(future_to_chunk), 1):
                chunk = future_to_chunk[future]
                try:
                    result = future.result()
                    results.append(result)

                    if result.get('tests_run', 0) > 0:
                        status = "✅ PASSED" if result['success'] else "❌ FAILED"
                        duration = f"({result['duration']:.1f}s)" if result['duration'] > 0 else ""
                        print(f"[{i}/{len(all_chunks)}] {result['test']:<25} {status} {duration}")
                        if not result['success'] and 'error' in result:
                            print(f"  Error: {result['error']}")
                except Exception as exc:
                    is_chunk = isinstance(chunk, dict)
                    test_name = f"{os.path.basename(chunk['file'] if is_chunk else chunk)[:-3]}"
                    if is_chunk:
                        test_name += f"_chunk_{chunk['chunk_id']}"
                    print(f"[{i}/{len(all_chunks)}] {test_name:<25} ❌ FAILED (Exception: {exc})")
                    results.append({'test': test_name, 'success': False, 'error': str(exc), 'duration': 0})

        # Print summary
        passed = sum(1 for r in results if r['success'])
        failed = len(results) - passed
        total_tests = sum(r.get('tests_run', 0) for r in results)
        duration = time.time() - start_time

        print("\n" + "=" * 60)
        print(f"FINAL SUMMARY")
        print("=" * 60)
        print(f"Total: {len([f for f in test_files if os.path.basename(f)[:-3] != 'test_api'])} test files, {len(results)} chunks, {total_tests} tests")
        print(f"Passed: {passed}")
        print(f"Failed: {failed}")
        print(f"Total duration: {duration:.1f}s")
        print("=" * 60)

        return failed == 0
    except KeyboardInterrupt:
        print("\n\nInterrupted by user. Cleaning up...")
        return False
    finally:
        signal.signal(signal.SIGINT, original_sigint)
        signal.signal(signal.SIGTERM, original_sigterm)


if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser(description='Metagraph integration tests.')
        parser.add_argument('--test_filter', dest='filter', type=str, default="*",
                            help='filter test cases (default: run all)')
        parser.add_argument('--gdb', dest='use_gdb', action='store_true',
                            help='run metagraph with gdb')
        parser.add_argument('--num_workers', dest='num_workers', type=int, default=None,
                            help='number of parallel workers (default: number of CPU cores)')
        args = parser.parse_args()

        # Default to number of CPU cores if not specified
        if args.num_workers is None:
            args.num_workers = os.cpu_count() or 1
        
        if args.use_gdb:
            update_prefix('gdb -ex run -ex bt -ex quit --args ')

        if args.num_workers > 1:
            # Run tests in parallel
            success = run_tests_parallel(args.num_workers, args.filter)
            exit(0 if success else 1)
        else:
            # Run tests sequentially (original behavior)
            result = unittest.TextTestRunner(
                resultclass=TimeLoggingTestResult
            ).run(create_test_suite(args.filter))
            exit(0 if result.wasSuccessful() else 1)

    except KeyboardInterrupt:
        sys.exit(130)

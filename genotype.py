#!/usr/bin/python3
"""
Multiplex genotyping script for interview coding challenges.
This probably would have been simpler with pandas,
but I chose to stick to python builtins.
"""
# imports
import sys
import argparse
import unittest

# StringIO import is different between python2.x and python3.x
try:
    # python 2.x: StringIO.StringIO
    from StringIO import StringIO
except ImportError:
    # python 3.x: io.StringIO
    from io import StringIO

# set script meta vars
__author__ = "Ryan Smith (ryan.smith.p@gmail.com)"
__version__ = "0.1.0"
__date__ = "2018-06-16"

def get_args(args):
    """Defines the command-line interface"""
    # init parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="author: " + __author__ + "\nversion: " + __version__
        )

    # add test argument
    parser.add_argument(
        '-t', '--test',
        action='store_true',
        required=False,
        help="Run unit tests"
        )

    # add input argument
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=False,
        help="Input file [stdin]"
        )

    # add output argument
    parser.add_argument(
        '-o', '--output',
        type=str,
        required=False,
        help="Output file [stdout]"
        )

    # parse the user arguments
    args = parser.parse_args(args)

    # send back the user input
    return args, parser

def set_io(args, parser):
    """Opens file handles depending on user arguments"""
    # if no input arg, check if part of pipe and if so, read stdin.
    if args.input is None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input = sys.stdin
    # otherwise, open input by filename
    else:
        args.input = open(args.input, 'r')
    # if no output specified, set to stdout
    if args.output is None:
        args.output = sys.stdout
    # otherwise, open output by filename
    else:
        args.output = open(args.output, 'w')
    return

def input_parser(input_file):
    """Generator function. Parses muliplexed genotyping experiments file"""
    # list to collect lines for each experiment
    experiment = []
    # iterate over lines in file
    for line in input_file:
        line = line.strip()
        # as long as line is not empty, append to experiment
        if line != '':
            experiment.append(line.split(","))
        # on empty lines, yield experiment list
        else:
            yield experiment
            experiment = []

    # make sure to yield final experiment (if no newline at last test)
    if experiment:
        yield experiment


def filter_mutants(mutant_sample_sets, normal_samples):
    """Filters any mutant from a test if called in normal tests"""
    # exclude normal sample_ids from individual mutant test sets
    return [test_set - normal_samples for test_set in mutant_sample_sets]


def get_single_mutants(filtered_mutants):
    """Returns a set of sample_ids uniquely identified as mutant"""
    # next(iter(x)) to get the value from a set of len 1
    #   WITHOUT leaving an empty set behind (with set.pop())
    return set([next(iter(x)) for x in filtered_mutants if len(x) == 1])


def check_unique(filtered_mutants, single_mutants):
    """Returns FALSE if any mutant samples were not uniquely identified"""
    # for each normal-filtered mutant test pool,
    #   remove any samples that were uniquely identified.
    #   if multiple samples still remain, they cannot be uniquely identified.
    # could report the ambiguous sample_ids and call the successful samples,
    #    but the instructions say that ALL samples must be mapped uniquely.
    for test_set in filtered_mutants:
        test_set = test_set - single_mutants
        if len(test_set) > 1:
            return False
    # return true if check passes
    return True


def check_consistent(filtered_mutants):
    """Checks if a filtered mutant callset has any nonunique samples"""
    # if any mutant test gets filtered to 0 elements,
    #   all of these test samples were identified in normal tests
    #   this can only happen with an erroneous NORMAL or MUTANT genotype call.
    for test_set in filtered_mutants:
        if len(test_set) == 0:
            # return False flag for failure
            return False

    # return true if check passes
    return True


def process_experiment(experiment):
    """Evaluates a single multiplexing experiment"""
    # stores normal sample ids
    normal_samples = set()
    # stores each mutant callset
    mutant_sample_sets = []
    # for each test pool,
    for test_set in experiment:
        # get genotype call from first item
        state = test_set[0]
        # all following items are sample_ids
        samples = test_set[1:]
        # collect mutant test pools in separate sets
        if state == "MUT":
            mutant_sample_sets.append(set(samples))
        # collect all normal test ids in a single set
        elif state == "NORM":
            normal_samples.update(samples)
    # for each mutant test set, remove samples found in normal pools
    filtered_mutants = filter_mutants(mutant_sample_sets, normal_samples)
    # get MUT tests that were narrowed down to a single mutant
    single_mutants = get_single_mutants(filtered_mutants)
    # check for uniqueness
    if not check_unique(filtered_mutants, single_mutants):
        return False, "NONUNIQUE"
    # check for consistency
    if not check_consistent(filtered_mutants):
        return False, "INCONSISTENT"
    # return True flag for success as well as the results
    return True, [single_mutants, normal_samples]


def output_results(success, result, output_file):
    """Writes the results of a genotyping experiment to output"""
    # write failure state if experiment failed
    if not success:
        output_file.write(result+"\n")
    elif success:
        # unpack results sets if success
        single_mutants, normal_samples = result
        # list to store all calls
        final_calls = []
        # get mutant and normal counts
        n_mut = len(single_mutants)
        n_norm = len(normal_samples)
        # append sample_ids with genotype call to final_calls
        for sample_id in single_mutants:
            final_calls.append((sample_id, "MUT"))
        for sample_id in normal_samples:
            final_calls.append((sample_id, "NORM"))

        # sort final_calls by sample_ids
        final_calls.sort(key=lambda x: int(x[0]))
        # write totals to output
        output_file.write("MUT COUNT: {0}\n".format(n_mut))
        output_file.write("NORM COUNT: {0}\n".format(n_norm))
        # write sample calls to output
        for line in final_calls:
            output_file.write(",".join(line)+"\n")

    # newline delimits each experiment
    output_file.write("\n")
    return


class TestArgs(unittest.TestCase):
    """Class for defining unit tests"""

    def test_get_args_testflag(self):
        """TEST get_args test FUNCTIONALITY"""
        # test flag
        args = get_args(["-t"])[0]
        self.assertTrue(args.test)

    def test_get_args_input(self):
        """TEST get_args input FUNCTIONALITY"""
        # only define input file
        args = get_args(["-i", "input"])[0]
        self.assertEqual(args.input, "input")
        self.assertEqual(args.output, None)

    def test_get_args_output(self):
        """TEST get_args output FUNCTIONALITY"""
        # only define output file
        args = get_args(["-o", "output"])[0]
        self.assertEqual(args.input, None)
        self.assertEqual(args.output, "output")

    def test_get_args_input_output(self):
        """TEST get_args input output FUNCTIONALITY"""
        # define both
        args = get_args(["-i", "input", "-o", "output"])[0]
        self.assertEqual(args.input, "input")
        self.assertEqual(args.output, "output")

    def test_get_args_no_args(self):
        """TEST get_args no args FUNCTIONALITY"""
        # define none
        args = get_args([])[0]
        self.assertEqual(args.input, None)
        self.assertEqual(args.output, None)

class TestUtilities(unittest.TestCase):
    """Class for defining utilities tests"""

    def setUp(self):
        """Init test module"""
        # expected filtered mutant sets for each test
        self.test_filtered_mutants = [
            [],
            [set(['12'])],
            [set([])],
            [set(['1', '0']), set(['1', '2'])]
        ]

        # mutant and normal sets for each test
        self.test_sample_sets = [
            [[], set(['1', '0', '2'])],
            [[set(['12', '110'])], set(['100', '110'])],
            [[set(['1', '2'])], set(['1', '0', '3', '2'])],
            [[set(['1', '0']), set(['1', '2'])], set([])]
        ]

        #single mutant sets for each test
        self.test_single_mutants = [set([]), set(['12']), set([]), set([])]
        # unique check result for each test
        self.test_uniques = [True, True, True, False]
        # consistency check result for each test
        self.test_consistents = [True, True, False, True]
        # test command line args

    def test_filter_mutants(self):
        """TEST filter_mutants FUNCTIONALITY"""
        # list to collect filtered mutants
        filtered_mutants = []
        # call filter_mutants on each test sample set
        for mutant_sample_sets, normal_samples in self.test_sample_sets:
            filtered_mutants.append(
                # sample_set is [mutant_sample_sets, normal_samples]
                filter_mutants(mutant_sample_sets, normal_samples)
            )
        # compare results to test_filtered_mutants
        self.assertEqual(filtered_mutants, self.test_filtered_mutants)

    def test_get_single_mutants(self):
        """TEST get_single_mutants FUNCTIONALITY"""
        # list to collect single mutants
        single_mutants = []
        # call get_single_mutants for
        for filtered_mutants in self.test_filtered_mutants:
            single_mutants.append(get_single_mutants(filtered_mutants))
        # compare results to test_single_mutants
        self.assertEqual(single_mutants, self.test_single_mutants)

    def test_check_unique(self):
        """TEST check_unique FUNCTIONALITY"""
        # list to collect unique mutants
        uniques = []
        # zip inputs together for looping
        test_input = zip(self.test_filtered_mutants, self.test_single_mutants)
        # call check_unique on each test input set
        for filtered_mutants, single_mutants in test_input:
            uniques.append(check_unique(filtered_mutants, single_mutants))
        # compare results to test_uniques
        self.assertEqual(uniques, self.test_uniques)

    def test_check_consistent(self):
        """TEST check_consistent FUNCTIONALITY"""
        # list collect consistency check output
        consistents = []
        # call check_consistent for each test
        for filtered_mutants in self.test_filtered_mutants:
            consistents.append(check_consistent(filtered_mutants))
        # compare results to test_consistentss
        self.assertEqual(consistents, self.test_consistents)

# TEST SUITE #
class TestGenotype(unittest.TestCase):
    """Class for defining higher-level genotyping functions"""

    def setUp(self):
        """Init test module"""
        # simulate input file for testing
        self.test_input_file = StringIO((
            "NORM,0,1\n"
            "NORM,1,2\n"
            "NORM,0,2\n"
            "\n"
            "NORM,100,110\n"
            "MUT,110,12\n"
            "\n"
            "NORM,0,1\n"
            "MUT,1,2\n"
            "NORM,1,3\n"
            "NORM,2,3\n"
            "\n"
            "MUT,0,1\n"
            "MUT,1,2\n"
            "\n"
            ))

        # expected output from test input
        self.test_experiments = [
            [
                ["NORM", "0", "1"],
                ["NORM", "1", "2"],
                ["NORM", "0", "2"]
            ],
            [
                ["NORM", "100", "110"],
                ["MUT", "110", "12"]
            ],
            [
                ["NORM", "0", "1"],
                ["MUT", "1", "2"],
                ["NORM", "1", "3"],
                ["NORM", "2", "3"]
            ],
            [
                ["MUT", "0", "1"],
                ["MUT", "1", "2"]
            ]
        ]

        # expected output of process_experiment for each test
        self.test_results = [
            (True, [set([]), set(['1', '0', '2'])]),
            (True, [set(['12']), set(['100', '110'])]),
            (False, 'INCONSISTENT'),
            (False, 'NONUNIQUE')
        ]

        # expected final output for testing
        self.test_output = (
            "MUT COUNT: 0\n"
            "NORM COUNT: 3\n"
            "0,NORM\n"
            "1,NORM\n"
            "2,NORM\n"
            "\n"
            "MUT COUNT: 1\n"
            "NORM COUNT: 2\n"
            "12,MUT\n"
            "100,NORM\n"
            "110,NORM\n"
            "\n"
            "INCONSISTENT\n"
            "\n"
            "NONUNIQUE\n"
            "\n"
        )

    def test_input_parser(self):
        """TEST input_parser FUNCTIONALITY"""
        # get all experiments from simulated input
        experiments = [x for x in input_parser(self.test_input_file)]
        # compare results to test_experiments
        self.assertEqual(experiments, self.test_experiments)

    def test_process_experiment(self):
        """TEST process_experiment FUNCTIONALITY"""
        # collect results for each test experiment
        results = []
        # call process_experiments for each test
        for experiment in self.test_experiments:
            results.append(process_experiment(experiment))
        # compare results to test_results
        self.assertEqual(results, self.test_results)

    def test_output_results(self):
        """TEST output_results FUNCTIONALITY"""
        # get StringIO simulated output file
        output_file = StringIO()
        # write results to StringIO
        for success, result in self.test_results:
            output_results(success, result, output_file)
        # compare output to test_output
        self.assertEqual(output_file.getvalue(), self.test_output)


def main():
    """Main function"""
    # get command line arguments
    args, parser = get_args(sys.argv[1:])
    # if test flag, run tests
    if args.test:
        # with tests in same file, unittest.main() still gets sys.argv
        # so, clear sys.argv before running unittest
        sys.argv = sys.argv[:1]
        unittest.main()
    # open file handles
    set_io(args, parser)
    # iterate over each experiment from input file using reader
    for experiment in input_parser(args.input):
        # get success flag and results
        success, results = process_experiment(experiment)
        # write the results to output
        output_results(success, results, args.output)


if __name__ == '__main__':
    main()

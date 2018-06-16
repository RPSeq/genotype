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
import StringIO

# set script meta vars
__author__ = "Ryan Smith (ryan.smith.p@gmail.com)"
__version__ = "0.1.0"
__date__ = "2018-06-13"

def get_args(args):
    """Defines the command-line interface"""
    # init parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="author: " + __author__ + "\nversion: " + __version__
        )

    # add input argument
    parser.add_argument(
        '-i', '--input',
        type=argparse.FileType('r'),
        required=False, default=None,
        help="Input file [stdin]"
        )

    # add output argument
    parser.add_argument(
        '-o', '--output',
        type=argparse.FileType('w'),
        required=False, default=sys.stdout,
        help="Output file [stdout]"
        )

    # parse the user arguments
    args = parser.parse_args(args)
    # if no input arg, check if part of pipe and if so, read stdin.
    if args.input is None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input = sys.stdin
    # send back the user input
    return args


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

    # we are finished with input
    input_file.close()
    # make sure to yield final experiment (if no newline at last test)
    if experiment:
        yield experiment


def filter_mutants(mutant_sample_sets, normal_samples):
    """Filters any mutant from a test if called in normal tests"""
    # exclude normal sample_ids from individual mutant test sets
    return [test - normal_samples for test in mutant_sample_sets]


def get_single_mutants(filtered_mutants):
    """Returns a set of sample_ids uniquely identified as mutant"""
    # next(iter(x)) to get the value from a set of len 1
    #   WITHOUT leaving an empty set behind.
    return set([next(iter(x)) for x in filtered_mutants if len(x) == 1])


def check_unique(normal_samples, single_mutants, all_samples):
    """Returns FALSE if any mutant samples were not uniquely identified"""
    # if single_mutants merged with normals doesn't cover all samples,
    #   at least one sample was never uniquely identified as mutant.
    # it would be possible to simply NOT call the ambiguous sample_ids,
    #    but the instructions say that ALL samples must be mapped uniquely.
    if len(normal_samples | single_mutants) != len(all_samples):
        # return False flag for failure
        return False
    # return true if check passes
    return True


def check_consistent(filtered_mutants):
    """Checks if a filtered mutant callset has any nonunique samples"""
    # if any mutant test gets filtered to 0 elements,
    #   all of these test samples were identified in normal tests
    #   this can only happen with an erroneous NORMAL or MUTANT genotype call.
    for test in filtered_mutants:
        if len(test) == 0:
            # return False flag for failure
            return False

    # return true if check passes
    return True


def process_experiment(experiment):
    """Evaluates a single multiplexing experiment"""
    # stores all sample ids
    all_samples = set()
    # stores normal sample ids
    normal_samples = set()
    # stores each mutant callset
    mutant_sample_sets = []
    # for each test pool,
    for test in experiment:
        # get genotype call from first item
        state = test[0]
        # all following items are sample_ids
        samples = test[1:]
        # add samples to all set
        all_samples.update(samples)
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
    if not check_unique(normal_samples, single_mutants, all_samples):
        return False, "NONUNIQUE"
    # check for consistency
    if not check_consistent(filtered_mutants):
        return False, "INCONSISTENT"
    # return True flag for success as well as the results
    return True, [single_mutants, normal_samples]


def writer(success, result, output_file):
    """Writes the results of a genotyping experiment to output"""
    # write failure state if experiment failed
    if not success:
        output_file.write(result+"\n")
    elif success:
        # unpack results sets if success
        single_mutants, normals = result
        # list to store all calls
        final_calls = []
        # get mutant and normal counts
        n_mut = len(single_mutants)
        n_norm = len(normals)
        # append sample_ids with genotype call to final_calls
        for sample_id in single_mutants:
            final_calls.append((sample_id, "MUT"))
        for sample_id in normals:
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


class Tester(unittest.TestCase):
    """Class for defining unit tests"""

    def setUp(self):
        """Init test module"""
        # simulate input file for testing
        self.test_input_file = StringIO.StringIO(
            (
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
            )
            )

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

    def test_input_parser(self):
        """TEST PARSER FUNCTIONALITY"""
        # get all experiments from simulated input
        experiments = [x for x in input_parser(self.test_input_file)]
        # check that output matches expectation
        self.assertEqual(experiments, self.test_experiments)

    def test_evaluate(self):
        """TEST PARSER FUNCTIONALITY"""
        #
        self.assertEqual(0, 0)

    def test_writer(self):
        pass


def main():
    """Main function"""
    # get command line arguments
    args = get_args(sys.argv[1:])
    # iterate over each experiment from input file using reader
    for experiment in input_parser(args.input):
        # get success flag and results
        success, results = process_experiment(experiment)
        # write the results to output
        writer(success, results, args.output)


if __name__ == '__main__':
    main()

#!/usr/bin/python3

# imports
import sys, argparse, unittest, StringIO, collections


# set script meta vars
__author__ = "Ryan Smith (ryan.smith.p@gmail.com)"
__version__ = "0.1.0"
__date__ = "2018-06-13"


#TODO test get_args
def getArgs(args):
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


def parser(inputFile):
    """Generator function. Parses muliplexed genotyping experiments file"""
    # list to collect lines for each experiment
    experiment = []
    # iterate over lines in file
    for line in inputFile:
        line = line.strip()
        # as long as line is not empty, append to experiment
        if line != '':
            experiment.append(line.split(","))
        # on empty lines, yield experiment list
        else:
            yield experiment
            experiment = []

    inputFile.close()

    # make sure to yield final experiment
    if experiment: yield experiment


def evaluate(experiment):
    """Evaluates a single multiplexing experiment"""
    # stores all sample ids
    all = set()
    # stores normal sample ids
    normals = set()
    # stores each mutant callset
    mutants = []

    # for each test,
    for test in experiment:
        # get MUT or NORM state at variant
        state = test[0]
        samples = test[1:]

        # add samples to all set
        all.update(samples)

        # collect mutant test sets separately
        if state == "MUT":
            mutants.append(set(samples))
        # collect normal test ids in a single set
        elif state == "NORM":
            normals.update(samples)

    # for each mutant test set, remove samples called in normal tests
    filteredMutants = [test - normals for test in mutants]

    # get MUT tests that can be narrowed down to a single mutant
    # next(iter(x)) is a trick to get the value from a set of len 1
    #   WITHOUT leaving an empty set behind.
    singleMutants = set([next(iter(x)) for x in filteredMutants if len(x) == 1])

    # if singleMutants merged with normals doesn't cover all samples,
    # at least one sample was never uniquely identified as mutant.
    if len(normals | singleMutants) != len(all):
        # return False flag for failure and the failure type
        return False, "NONUNIQUE"

    # if any mutant test gets filtered to 0 elements,
    # all of these test samples were identified in normal tests
    # this can only happen with an erroneous NORMAL or MUTANT genotype call.
    for test in filteredMutants:
        if len(test) == 0:
            # return False flag for failure and the failure type
            return False, "INCONSISTENT"

    # return True flag for success as well as the results
    return True, [singleMutants, normals]

def writer(success, result, outputFile):
    if not success:
        outputFile.write(result)
        outputFile.write("\n")

    elif success:
        # unpack results sets
        singleMutants, normals = result
        # list to store all calls
        all = []

        # get mutant and normal counts
        nMut = len(singleMutants)
        nNorm = len(normals)

        # append sampleIDs and corresponding results to all
        for sampleID in singleMutants:
            all.append((sampleID, "MUT"))
        for sampleID in normals:
            all.append((sampleID, "NORM"))

        # sort all by sampleIDs
        all.sort(key=lambda x: int(x[0]))

        # write totals to output
        outputFile.write("MUT COUNT: {0}\n".format(nMut))
        outputFile.write("NORM COUNT: {0}\n".format(nNorm))

        # write sample calls to output
        for line in all:
            outputFile.write(",".join(line)+"\n")

    # newline delimits each experiment
    outputFile.write("\n")


class Tester(unittest.TestCase):
    """Class for defining unit tests"""
    def setUp(self):
        # simulate input file for testing
        self.inputFile = StringIO.StringIO((
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
            "\n"))

        # expected output from test input
        self.experiments = [
            [["NORM","0","1"],
            ["NORM","1","2"],
            ["NORM","0","2"]],

            [["NORM","100","110"],
            ["MUT","110","12"]],

            [["NORM","0","1"],
            ["MUT","1","2"],
            ["NORM","1","3"],
            ["NORM","2","3"]],

            [["MUT","0","1"],
            ["MUT","1","2"]]
            ]

    def test_parser(self):
        """TEST PARSER FUNCTIONALITY"""
        # get all experiments from simulated input
        self.testExperiments = [x for x in parser(self.inputFile)]
        # check that output matches expectation
        self.assertEqual(self.testExperiments, self.experiments)

    # TODO implement this test
    def test_evaluate(self):
        """TEST PARSER FUNCTIONALITY"""
        self.assertEqual(0,0)

def main():
    # get command line arguments
    args = getArgs(sys.argv[1:])

    # iterate over each experiment from input file using reader
    for experiment in parser(args.input):
        # get success flag and results by evaluating experiment
        success, results = evaluate(experiment)
        # write the results to output
        writer(success, results, args.output)

if __name__ == '__main__':
    main()

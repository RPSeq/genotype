#!/usr/bin/python3

# imports
import sys, argparse, unittest, StringIO


# set script meta vars
__author__ = "Ryan Smith (ryan.smith.p@gmail.com)"
__version__ = "0.1.0"
__date__ = "2018-06-13 "


#TODO test get_args
def getArgs(args):
    """Defines the command-line interface"""
    # init parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="author: " + __author__ + "\nversion: " + __version__)

    # add input argument
    parser.add_argument(
        '-i',
        '--input',
        type=argparse.FileType('r'),
        required=False,
        default=None,
        help="Input file [stdin]")

    # add output argument
    parser.add_argument(
        '-o',
        '--output',
        type=argparse.FileType('w'),
        required=False,
        default=sys.stdout,
        help="Output file [stdout]")

    # parse the user arguments
    args = parser.parse_args(args)
    # if no input, check if part of pipe and if so, read stdin.
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

class TestReader(unittest.TestCase):
    """Unittest class for testing reader functionality."""
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

    def test_reader(self):
        """TEST READER FUNCTIONALITY"""
        # get all experiments from simulated input
        experiments = [x for x in parser(self.inputFile)]
        # check that output matches expectation
        self.assertEqual(experiments,self.experiments)



# TODO test process
def evaluate(experiment):
    """Evaluates a single multiplexing experiment"""
    # sets to store normal and mutant sample IDs
    all = set()
    normals = set()
    mutantSet = set()
    mutants = []

    # for each test,
    for test in experiment:
        # get MUT or NORM state at variant
        state = test[0]
        samples = test[1:]
        all.update(samples)
        # collect mutant test sets separately
        if state == "MUT":
            mutantSet.update(samples)
            mutants.append(set(samples))
        # collect normal test ids in a single set
        elif state == "NORM":
            normals.update(samples)

    # for each mutant test set, remove samples called in normal tests
    filteredMutants = [test - normals for test in mutants]


    # get tests narrowed down to a single mutant
    singleMutants = set([next(iter(x)) for x in filteredMutants if len(x) == 1])


    if len(singleMutants) < 1 and len(mutants) > 0:
        # print(singleMutants)
        print("NONUNIQUE")

    if len(normals | singleMutants) != len(all):
        print("NONUNIQUE")

    for x in filteredMutants:
        if len(x) == 0:
            print("INCONSISTENT")
            break

    print(len(singleMutants))
    print(len(normals))
    print(len(normals) + len(singleMutants))
    print(len(normals | singleMutants))
    print(len(all))
    # TODO implement unique mapping test
def main():
    # args = getArgs(sys.argv[1:])
    # for exp in reader(args.input):
    #     print("____")
    #     evaluate(exp)
    unittest.main()



if __name__ == '__main__':
    main()

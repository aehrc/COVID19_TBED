#!/bin/python

#------------------- Description & Notes --------------------#

'''
## Arguments
-f = FASTA files (.fa/.gz) / FASTQ files (.fa/.gz)
-g = GFF files (.gff)

-k = Kmer length
-w = Window size (bp)
-s = Step size (bp)
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import sys
import os
import argparse
from argparse import ArgumentParser
from argparse import ArgumentTypeError

# External imports

# Internal imports

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

class ArgParser(ArgumentParser):

    def __init__(self):
        ArgumentParser.__init__(self)

    def parse(self):
        pass

    def printArgs(self):
        args = self.parse_args()
        print(args)

    def isGTZeroInt(self, value):
        try:
            intValue = int(value)
            if (intValue <= 0):
                errMsg = "Invalid value. Must be > 0"
                raise ArgumentTypeError(errMsg)

        except (ArgumentTypeError, ValueError) as err:
            sys.exit(err)

        return intValue

    def isNumeric(self, value):
        try:
            floatValue = float(value)

        except (ValueError) as error:
            sys.exit(error)

        return floatValue

    def isFile(self, value):
        try:
            if (not os.path.exists(value) or os.path.isdir(value)):
                errMsg = "Invalid file:\t" + value
                raise ArgumentTypeError(errMsg)

        except ArgumentTypeError as err:
            sys.exit(err)

        return value

    def isDir(self, value):
        try:
            if (not os.path.exists(value) or not os.path.isdir(value)):
                errMsg = "Invalid dir:\t" + value
                raise ArgumentTypeError(errMsg)

        except ArgumentTypeError as err:
            sys.exit(err)

        return value

#------------------------------------------------------------#

class GenerateFastxArgParser(ArgParser):

    def __init__(self):
        ArgParser.__init__(self)

        ## Fixes an error that I really don't understand
        ## Related to the subclassing of ArgParser
        subparsers = self.add_subparsers(dest='cmd')
        subparsers._parser_class = ArgumentParser

        ## Create a subparsers
        self.createFixedSubparser(subparsers)
        self.createFeatureSubparser(subparsers)
        self.createHybridSubparser(subparsers)

    def createFixedSubparser(self, subparsers):
        p = subparsers.add_parser('fixed')
        p.add_argument("-w", help="Sequence record length (bp)",
            nargs=1, type=self.isGTZeroInt, required=True)
        p.add_argument("-i", help="ID \
            (Must match an ID in the FASTA/FASTQ files)",
            nargs='?')
        p.add_argument("-q", help="Write output as FASTQ reads",
            action='store_true')
        self.initArgs(p)

    def createFeatureSubparser(self, subparsers):
        p = subparsers.add_parser('feature')
        p.add_argument("-g", help="Genome annotation files (.gff)",
            nargs='+', type=self.isFile, required=True)
        p.add_argument("-i", help="ID \
            (Must match a feature type in the GFF file)",
            nargs='?')
        self.initArgs(p)

    def createHybridSubparser(self, subparsers):
        p = subparsers.add_parser('hybrid')
        p.add_argument("-s", help="ID \
            (Must match an ID in the FASTA/FASTQ files)",
            nargs='?')
        p.add_argument("-t", help="ID \
            (Must match an ID in the FASTA/FASTQ files)",
            nargs='?')
        self.initArgs(p)

    def initArgs(self, p):
        p.add_argument("-f", help="FASTA files (.fa/.gz)",
            nargs='+', type=self.isFile, required=True)
        p.add_argument("-o", help="Output file (.fa)",
            nargs=1, required=True)
        p.add_argument("-n", help="Number of sequence records",
            nargs=1, type=self.isGTZeroInt, required=True)

    def parse(self):
        try:
            args = self.parse_args()

            self.cmd      = args.cmd
            self.iFiles   = args.f
            self.oFile    = args.o[0]
            self.nSeqRecs = args.n[0]

        except (AttributeError, TypeError):
            self.print_help()
            sys.exit(2)

    def parseFixedArgs(self):
        args = self.parse_args()

        self.windowSize = args.w[0]
        self.fastxId    = None if args.i is None else args.i
        self.toFastq    = args.q

    def parseFeatureArgs(self):
        args = self.parse_args()

        self.annoFiles = args.g
        self.featureId = None if args.i is None else args.i

    def parseHybridArgs(self):
        args = self.parse_args()

        self.sId = None if args.s is None else args.s
        self.tId = None if args.t is None else args.t

#------------------------------------------------------------#

class CalculateKmerArgParser(ArgParser):

    def __init__(self):
        ArgParser.__init__(self)

        ## Fixes an error that I really don't understand
        ## Related to the subclassing of ArgParser
        subparsers = self.add_subparsers(dest='cmd')
        subparsers._parser_class = ArgumentParser

        ## Create a subparsers
        self.createCombinedSubparser(subparsers)
        self.createSplitSubparser(subparsers)
        self.createSequentialSubparser(subparsers)

    def createCombinedSubparser(self, subparsers):
        p = subparsers.add_parser('combined')
        self.initArgs(p)

    def createSplitSubparser(self, subparsers):
        p = subparsers.add_parser('split')
        self.initArgs(p)

    def createSequentialSubparser(self, subparsers):
        p = subparsers.add_parser('sequential')
        p.add_argument("-i", help="ID \
            (Must match an ID in the FASTA/FASTQ files)",
            nargs='?')
        p.add_argument("-w", help="Number of windows",
            nargs=1, type=self.isGTZeroInt, required=True)
        p.add_argument("-s", help="Number of steps",
            nargs=1, type=self.isGTZeroInt, required=True)
        self.initArgs(p)

    def initArgs(self, p):
        p.add_argument("-f", help="FASTA/FASTQ files (.fa/.fq)",
            nargs='+', type=self.isFile, required=True)
        p.add_argument("-k", help="Kmer length",
            nargs=1, type=self.isGTZeroInt, required=True)
        p.add_argument("-n", help="Include ambiguous bases (i.e., N's) \
            in the total count", action='store_false')
        p.add_argument("-o", help="Output file (.snappy.parquet)",
            nargs=1, required=True)

    def parse(self):
        try:
            args = self.parse_args()

            self.cmd        = args.cmd
            self.iFiles     = args.f
            self.kmerLength = args.k[0]
            self.includeN   = args.n
            self.oFile      = args.o[0]

        except (AttributeError, TypeError):
            self.print_help()
            sys.exit(2)

    def parseSequentialArgs(self):
        args = self.parse_args()

        self.fastxId  = None if args.i is None else args.i
        self.nWindows = args.w[0]
        self.nSteps   = args.s[0]

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------


#!/usr/bin/env python
# TITLE: script to find reg ex of interest
# author: Peter Thorpe September 2019.
# University of St Andrews


# imports
from Bio.Seq import Seq
from Bio import SeqIO
import time
import os
import errno
import logging
import logging.handlers
import sys
import re
from optparse import OptionParser


def find_reg(sequences, pattern):
    """func to iterate through the fasta and search for the
    reg expressions"""
    for seq_record in SeqIO.parse(sequences, "fasta"):
        # https://pythonforbiologists.com/regular-expressions
        # https://docs.python.org/2/howto/regex.html
        # eg GG(A|T)CC
        #GC(A|T|G|C)GC
        # if re.search(r"GG(A|T)CC", dna):
        # if re.search(r"GC[ATGC]GC", dna):
        # we dont want to match [^XYZ]
        # A question mark immediately following a character means that that character is optional – it can match either zero
        # or one times. So in the pattern GAT?C t
        # A plus sign immediately following a character or group means that the character or group must be present but can be repeated any number of times 
        # An asterisk immediately following a character or group means that the character or group is optional, but can also be repeated. In other words,
        #N-Arg_dibasic_convertase    [A-Z][RK|R]R[R[A-Z]]
        search_matches = re.search(pattern,
                            str(seq_record.seq.upper()))
        findallmatches = re.findall(pattern,
                            str(seq_record.seq.upper()))
        if search_matches:
            print('Match found: ', m.group())
            print(search_matches)
            print(search_matches.group())
            print(search_matches.start(), search_matches.end())
        else:
            print('No match')

        


usage = """Use as follows:

$ python reg_ex_finder.py -h
python reg_ex_finder.py.py -i in.fasta --reg expression -o outdata.txt

requires: biopython
"""

parser = OptionParser(usage=usage)

parser.add_option("-i",
                  dest="infasta",
                  default=os.path.join("test_data",
                                       "seq_of_interest.fasta"),
                  help="this is the fasta file with your sequences")
parser.add_option("-r", "--reg",
                  dest="reg",
                  default=r"[A-Z](RK|R)R(R[A-Z])",
                  help="reg ex of interests")
parser.add_option("-o", "--output",
                  dest="out_file",
                  default="tes_reg.txt",
                  help="Output filename ",
                  metavar="FILE")

# get the user options. TODO. arg parser instead
(options, args) = parser.parse_args()
infasta = options.infasta
reg = options.reg
outfile = open(options.out_file, "w")
logfile = "reg_ex_WARNINGS.log"
# Run as script
if __name__ == '__main__':
    start_time = time.time()
    # Set up logging
    logger = logging.getLogger('reg_ex_finder.py: %s'
                               % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)
    try:
        logstream = open(logfile, 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        logger.error("Could not open %s for logging" %
                     logfile)
        sys.exit(1)
    file_list = [infasta] #infasta]
    for user_file in file_list:
        if not os.path.isfile(user_file):
            print("cannot find file %s " % user_file)
    logger.info(sys.version_info)
    logger.info("Command-line: %s", ' '.join(sys.argv))
    logger.info("Starting testing: %s", time.asctime())
    # index the genome with biopython
    find_reg(infasta, reg)


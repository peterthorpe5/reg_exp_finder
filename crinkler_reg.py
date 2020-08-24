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
from collections import OrderedDict


def parse_db_file(database):
    """func to parse a tab sep database file.
    return two dictionaries
    db_name_to_exp, db_name_to_desc
    """
    db_name_to_exp = OrderedDict()
    db_name_to_desc = OrderedDict()
    f_open = open(database, "r")
    for line in f_open:
        if not line.strip():
            continue  # if the last line is blank
        if line.startswith("#"):
            continue  # comment line
        name, reg_exp, description = line.split()
        db_name_to_exp[name] = reg_exp
        db_name_to_desc[name] = description.rstrip()
    return db_name_to_exp, db_name_to_desc


def count_motifs(seq, pattern):
    chunks = re.split(pattern, seq)
    return(len(chunks)-1)


def find_reg(sequences, database, out, pattern, logger):
    """func to iterate through the fasta and search for the
    reg expressions"""
    f_out = open(out, "w")
    motifs = out + "_motifs_only"
    motifs_out = open(motifs, "w")
    if database:
        db_name_to_exp, db_name_to_desc = parse_db_file(database)
    for seq_record in SeqIO.parse(sequences, "fasta"):
        #print(seq_record.id)
        # https://pythonforbiologists.com/regular-expressions
        # https://docs.python.org/2/howto/regex.html
        # eg GG(A|T)CC
        # GC(A|T|G|C)GC
        # if re.search(r"GG(A|T)CC", dna):
        # if re.search(r"GC[ATGC]GC", dna):
        # we dont want to match [^XYZ]
        # A question mark immediately following a character means that that character is optional – it can match either zero
        # or one times. So in the pattern GAT?C t
        # A plus sign immediately following a character or group means that the character or group must be present but can be repeated any number of times 
        # An asterisk immediately following a character or group means that the character or group is optional, but can also be repeated. In other words,
        # N-Arg_dibasic_convertase    [A-Z][RK|R]R[R[A-Z]]
        seq = str(seq_record.seq.upper())
        search_matches = re.search(pattern, seq)
        findallmatches = re.findall(pattern, seq)
        all_matches = re.finditer(pattern, seq)
        for postivies in all_matches:
            base = postivies.group() 
            pos  = postivies.start() 
            info = (base + " found at position " + str(pos) + "-", postivies.end() )
            logger.info(info)
            print(info)
            SeqIO.write(seq_record, f_out, "fasta")
            seq_record.seq = seq_record.seq[int(pos): int(postivies.end())]
            SeqIO.write(seq_record, motifs_out, "fasta")
            
        if findallmatches:
            print(len(findallmatches))
        num_of_motifs = count_motifs(seq, pattern)
        print("we have %d number of motifs matches" % num_of_motifs)
        print("findallmatches", findallmatches)
        
        if search_matches:
            print('Match found: ', search_matches)
            #print(search_matches)
            #print(search_matches.group())
            #print(search_matches.start(), search_matches.end())

        else:
            print('No match')
        if database:
            for name, exp in db_name_to_exp.items():
                 search_matches = re.search(exp,
                            str(seq_record.seq.upper()))
                 if search_matches:
                     print("Match found: ", name,  search_matches.group(),
                           db_name_to_desc[name])
                     data = "\t".join([name,  str(search_matches.group()),
                                       db_name_to_desc[name]])
                

        


usage = """Use as follows:

$ python reg_ex_finder.py -h
python reg_ex_finder.py.py -i in.fasta --reg expression -o outdata.txt

requires: biopython
"""

parser = OptionParser(usage=usage)

parser.add_option("-i",
                  dest="infasta",
                  default=os.path.join("test_data",
                                       "phy_crnker.fasta"),
                  help="this is the fasta file with your sequences")
parser.add_option("-r", "--reg",
                  dest="reg",
                  default="^.{30,70}L[FY]LA[RK]",
                  help="reg ex of interests")

parser.add_option("-f", "--database",
                  dest="database",
                  default=False,
                  help="a tab separated text file with " +
                  " patterns of interets: " +
                  "pattern_name\treg_expr\tdescription")

parser.add_option("-o", "--output",
                  dest="out_file",
                  default="matches.fasta",
                  help="Output filename ",
                  metavar="FILE")

# get the user options. TODO. arg parser instead
(options, args) = parser.parse_args()
infasta = options.infasta
reg = options.reg
outfile = options.out_file
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
    find_reg(infasta, options.database, outfile, reg, logger)


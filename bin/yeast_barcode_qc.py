#!/usr/bin/env python

# TODO: WRITE SCRIPT TO DO BARCODE LEVEL QC ON THE YEAST .BED FILE

"""

"""
# standard library
import os
import sys
import argparse
# outside dependencies
import pandas as pd


def parse_args(args=None):
    # Description = "Quantify the 'hops' of each TF in the bam file"
    # Epilog = "Example usage: python add_read_group_and_tags.py \
    #     <input.bam> <output.bam> <genome.fasta> <genome.fasta.fai> \
    #     <id_length> <optional: insertion_length (default is 1)>"

    # parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    # parser.add_argument("bampath",
    #                      help="path to the input bam file")
    # parser.add_argument("paired_end",
    #                      help="Boolean. True if the bam file is created from paired end reads. May be either a string true/false or 1/0")
    # parser.add_argument("mapq_filter",
    #                      help = "minimum value above which to accept alignments")

    # return parser.parse_args(args)
    raise NotImplementedError


def a_function(param1, param2, param3):
    """Description

    :param param1:
    :param param2:
    :param param3:

    :return:

    """
    raise NotImplementedError

def main(args=None):
    args = parse_args(args)

    # Check inputs
    input_path_list = [args.param1, args.param2,]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" %input_path)

    # do stuff
    a_function(args.param1,
               bool(args.param2),
               int(args.param3))

    sys.exit(0)


if __name__ == "__main__":
    sys.exit(main())

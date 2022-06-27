#!/usr/bin/env python

# TODO: write docstring which includes an autoSQL specification of the custom
#       bed file

"""
written by: chase mateusiak, chasem@wustl.edu

Bed format v1.0
https://github.com/samtools/hts-specs/blob/master/BEDv1.pdf
"""
# standard library
import os
import sys
import argparse
# outside dependencies
import pandas as pd


def parse_args(args=None):
    Description = "Quantify the 'hops' of each TF in the bam file"
    Epilog = "Example usage: python add_read_group_and_tags.py \
        <input.bam> <output.bam> <genome.fasta> <genome.fasta.fai> \
        <id_length> <optional: insertion_length (default is 1)>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("bampath",
                         help="path to the input bam file")
    parser.add_argument("single_end",
                         help="Boolean. True if the bam file is created from \
                            single end reads. May be either a string true/false \
                                (case doesn't matter) or 1/0")
    parser.add_argument("mapq_filter",
                         help = "minimum value above which to accept alignments")

    return parser.parse_args(args)


def split_bed(bed_file, barcode_map):
    """

    """
    raise NotImplementedError()

def main(args=None):
    args = parse_args(args)

    # Check inputs
    input_path_list = [args.bampath]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" %input_path)

    # loop over the reads in the bam file and add the read group (header and tag)
    # and the XI and XZ tags
    split_bed(args.bed_file,
              args.barcode_map)

    sys.exit(0)


if __name__ == "__main__":
    sys.exit(main())

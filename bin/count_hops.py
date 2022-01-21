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
import pysam
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


def add_read_group_and_tags(bampath, single_end, mapq_filter):
    """Iterate through the bam file and extract to bed custom format reads which
    meet quality/alignment thresholds

    :param bampath: path to the sorted, indexed bam file with RG, XZ and XI tags
    :param paired_end: whether or not the library is from paired end sequencing
                       (True) or not (False)
    :param mapq_filter: exclude reads less than or equal to this value

    :return: none. Write the custom bed file to file as tsv. output name is
             the input name with .bam replaced by .bed

    """

    d = {
        'chrom':           [],
        'chromStart':      [],
        'chromEnd':        [],
        'barcode':         [],
        'strand':          [],
        'insert_seq':      [],
        'aln_flag':        [],
        'aln_mapq':        [],
    }

    outpath = os.path.splitext(os.path.basename(bampath))[0] + ".bed"

    input_bamfile = pysam.AlignmentFile(bampath, "rb")
    for read in input_bamfile.fetch():

        # determine if the read should be counted. In either case -- paired
        # or not -- throw out the read if the mapping quality is too low
        # TODO consider adding another condition on that the read must be
        # unique, using tags. however, likely the mapq is stringent enough
        count_read = False
        if single_end:
            # do not count read if unmapped, secondary, failing qc
            # or is supplementary alignment
            count_read = True if \
                sum([int(read.flag) & x for x in [0x4,0x100,0x200,0x800]]) == 0 \
                and int(read.mapping_quality) >= mapq_filter \
                else False
        else:
            # 83 and 99 are the flags for correctly paired R1 on reverse and
            # forward strand respectively
            count_read = True if int(read.flag) in [83,99] and \
                int(read.mapping_quality) >= mapq_filter \
                else False

        if count_read:
            d['chrom'].append(input_bamfile.get_reference_name(read.reference_id))
            d['chromStart'].append(read.get_tag("XI"))
            d['chromEnd'].append(int(read.get_tag("XI")) + len(read.get_tag("XZ")))
            d['barcode'].append(read.get_tag("RG"))
            d['strand'].append('-' if int(read.flag) & 0x10 else "+")
            d['insert_seq'].append(read.get_tag("XZ"))
            d['aln_flag'].append(read.flag)
            d['aln_mapq'].append(read.mapping_quality)

    input_bamfile.close()

    df = pd.DataFrame(data=d)

    df.to_csv(outpath, sep = "\t")

def main(args=None):
    args = parse_args(args)

    def parse_bool(single_end):
        single_end_uncase = single_end.lower()
        switcher = {
            "0": False,
            "1": True,
            "true": True,
            "false": False
        }

        if switcher.get(single_end_uncase, 1) is 1:
            raise ValueError("Invalid choice for argument single_end: %s. \
                Must be one of 0,1, true/True/TRUE, false/False/FALSE" %single_end)
        else:
            return switcher.get(single_end_uncase)

    single_end = parse_bool(args.single_end)


    # Check inputs
    input_path_list = [args.bampath]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" %input_path)

    # loop over the reads in the bam file and add the read group (header and tag)
    # and the XI and XZ tags
    add_read_group_and_tags(args.bampath,
                            single_end,
                            int(args.mapq_filter))

    sys.exit(0)


if __name__ == "__main__":
    sys.exit(main())

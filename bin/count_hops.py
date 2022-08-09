#!/usr/bin/env python


"""
written by: chase mateusiak, chasem@wustl.edu

Loop over alignments. Only read1 are recorded.

Reads are considering passing if they are primary alignments which are not
unmapped, secondary, failing sequencer/aligner QC or supplementary alignment.
If paired end, then the reads must be correctly paired.

Reads may also be filtered by a mapq threshold.

If require_exactly_length is set to True, then only reads which
align without soft clipping are retained.

Output are two bed files in a format specified below. One file has suffix
_passing.bed, the other _failing.bed, and represent non overlapping paritions
of read1.

Bed format v1.0
https://github.com/samtools/hts-specs/blob/master/BEDv1.pdf

Historic Note: This is a 'modified' ccf file which more closely follows the
bed specifications. Each entry is a hop -- no aggregation at this point.

The output follows a Bed6+3 format with the following fields:

chrom, chromStart, chromEnd, name(barcode),
score(aln_mapq), strand, reads*, insert_seq, aln_flag

* reads are set to 1 -- these are not aggregated at this stage, so every row
represents a single hop. This field exists in order to accomodate previous
scripts
"""

# standard library
import os
import sys
import argparse
from copy import deepcopy
import logging
# outside dependencies
import pysam
import pandas as pd

logging.basicConfig(
    level=logging.CRITICAL,
    format='[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
    datefmt='%H:%M:%S',
    stream=sys.stdout
)


BED_6_3_COLNAMES = ['chrom', 'chromStart', 'chromEnd',
                    'barcode', 'aln_mapq', 'strand',
                    'reads', 'insert_seq', 'aln_flag']

def parse_args(args=None):
    Description = "Extract reads which potentially describe transpositions " +\
    "from a bam file, transform to bed6+3 format with column headers " +\
         " ".join(BED_6_3_COLNAMES)
    Epilog = "Example usage: python add_read_group_and_tags.py \
        <input.bam> <output.bam> <genome.fasta> <genome.fasta.fai> \
        <id_length> <optional: insertion_length (default is 1)>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("bampath",
                         help="path to the input bam file")
    parser.add_argument("require_exact_length",
                        help="True to filter out any soft-clipped reads, False "+\
                        "otherwise. Default is False",
                        default='False')
    parser.add_argument("mapq_filter",
                         help = "minimum value above which to accept alignments")

    return parser.parse_args(args)


def add_read_group_and_tags(bampath, require_exact_length, mapq_filter, dry_run=False):
    """Iterate through the bam file and extract to bed custom format reads which
    meet quality/alignment thresholds

    :param bampath: path to the sorted, indexed bam file with RG, XZ and XI tags
    :param paired_end: whether or not the library is from paired end sequencing
                       (True) or not (False)
    :param mapq_filter: exclude reads less than or equal to this value
    :param require_exact_length: True to filter out any reads with soft clipping
    :param dry_run: default False. Set to True for testing. Purpose is to be
        able to set a breakpoint before the write statement. Setting to True
        prevents writing to file.

    :return: none. Write the custom bed file to file as tsv. output name is
             the input name with .bam replaced by .bed

    """

    # create dictionaries to hold the paritions of read1s
    passing_reads = {k:[] for k in BED_6_3_COLNAMES}
    failing_reads = deepcopy(passing_reads)

    # create output paths
    out_dict = {
        'passing': os.path.splitext(os.path.basename(bampath))[0] + "_passing.bed",
        'failing': os.path.splitext(os.path.basename(bampath))[0] + "_failing.bed"
    }

    # open the bame file and begin looping over reads
    input_bamfile = pysam.AlignmentFile(bampath, "rb")
    for read in input_bamfile.fetch():

        # only consider read 1
        if read.is_read1 or not read.is_paired:

            count_read = False

            # Below are further criteria to determine if read should be counted.
            # In either case -- paired or not -- throw out the read if the
            # mapping quality is too low

            if not read.is_paired:
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

            # if require_exact_length is true,
            # then do not allow soft-clipped reads
            if require_exact_length:
                if read.is_reverse:
                    count_read = True if read.cigartuples[-1][0] != 4 else False
                else:
                    count_read = True if read.cigartuples[0][0] != 4 else False

            # add read to correct dict -- this partitions read1s
            if count_read:
                try:
                    passing_reads['chrom']\
                        .append(input_bamfile\
                            .get_reference_name(read.reference_id))
                    passing_reads['chromStart'].append(read.get_tag("XI"))
                    passing_reads['chromEnd'].append(int(read.get_tag("XI")) +
                                                     len(read.get_tag("XZ")))
                    passing_reads['barcode'].append(read.get_tag("RG"))
                    passing_reads['aln_mapq'].append(read.mapping_quality)
                    passing_reads['strand']\
                        .append('-' if int(read.flag) & 0x10 else "+")
                    passing_reads['reads'].append(1)
                    passing_reads['insert_seq'].append(read.get_tag("XZ"))
                    passing_reads['aln_flag'].append(read.flag)
                except KeyError as e:
                    msg = "A expected key in the passing_reads dict does not match " +\
                        "the bed fields. If you see this error, you should open an " +\
                            "issue report on github. Please post this error:\n%s" %e
                    logging.critical(msg+" ", exc_info=(sys.exc_info()))
                    raise
            else:
                try:
                    failing_reads['chrom']\
                        .append(input_bamfile\
                            .get_reference_name(read.reference_id))
                    failing_reads['chromStart'].append(read.get_tag("XI"))
                    try:
                        failing_reads['chromEnd'].append(int(read.get_tag("XI")) +
                                                        len(read.get_tag("XZ")))
                    except ValueError as e:
                        print(e)
                        failing_reads['chromEnd'].append("*")
                    failing_reads['barcode'].append(read.get_tag("RG"))
                    failing_reads['aln_mapq'].append(read.mapping_quality)
                    failing_reads['strand']\
                        .append('-' if int(read.flag) & 0x10 else "+")
                    failing_reads['reads'].append(1)
                    failing_reads['insert_seq'].append(read.get_tag("XZ"))
                    failing_reads['aln_flag'].append(read.flag)
                except KeyError as e:
                    msg = "A expected key in the failing_reads dict does not match " +\
                        "the bed fields. If you see this error, you should open an " +\
                            "issue report on github. Please post this error:\n%s" %e
                    logging.critical(msg+" ", exc_info=(sys.exc_info()))
                    raise

    # close the bamfile
    input_bamfile.close()
    # convert read dicts to dataframes
    passing_reads_df = pd.DataFrame(data=passing_reads)
    failing_reads_df = pd.DataFrame(data=failing_reads)
    # write out
    if not dry_run:
        passing_reads_df.to_csv(out_dict.get('passing'),
                                sep = "\t",
                                index = False,
                                header = None)

        failing_reads_df.to_csv(out_dict.get('failing'),
                                sep = "\t",
                                index = False,
                                header = None)

def main(args=None):
    args = parse_args(args)

    def parse_bool(val):
        # translate from possible input to boolean
        switcher = {
            "0": False,
            "1": True,
            "true": True,
            "false": False
        }
        # throw error if cannot cast val to boolean based on switch statement
        if not isinstance(switcher.get(val.lower()),bool):
            raise ValueError("Invalid choice for argument single_end: %s. \
                Must be one of 0,1, true/True/TRUE, false/False/FALSE" %single_end)

        return switcher.get(val.lower())

    require_exact_length = parse_bool(args.require_exact_length)

    # Check inputs
    input_path_list = [args.bampath]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" %input_path)

    # loop over the reads in the bam file and add the read group (header and tag)
    # and the XI and XZ tags
    add_read_group_and_tags(args.bampath,
                            require_exact_length,
                            int(args.mapq_filter))

    sys.exit(0)


if __name__ == "__main__":
    sys.exit(main())

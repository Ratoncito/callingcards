#!/usr/bin/env python

"""
written by: chase mateusiak, chasem@wustl.edu

This script takes as input a bam file which has a barcode, probably added by
UMITools, of length n added to the end of each QNAME. The barcode is extracted
based on length (eg, if the barcode is 13 characters, then 13 characters are
extracted from the end of the QNAME string of each read). Currently, this
script loops over the bam file twice -- the first time, extracting the barcode
from each read, during which a unique set of barcodes are created. Then, a
write-able bamfile is opened, this unique set of barcodes is added as a RG
(read group) header, which makes parsing, eg splitting the bam into parts by
RG possible. Then, a second loop is performed during which the tags RG, XZ and
XI are added to each read line. RG is the read group (barcode), XZ is the
coordinate at which the transposon inserted on a given chromosome, and XI is the
sequence x bases upstream of the insertion site. If the read is unmapped,
XI and XZ are set to "*".
"""
# standard library
import os
import sys
import argparse
# outside dependencies
import pysam


def parse_args(args=None):
    Description = "Extract the barcode, added by UMItools, from each read id, \
                   Add the Read Group (RG) header, and add the appropriate \
                   RG tag to each alignment record, and add the XI and XZ tags \
                   which describe the insertion location and sequence at that \
                   location"
    Epilog = "Example usage: python add_read_group_and_tags.py \
        <input.bam> <output.bam> <genome.fasta> <genome.fasta.fai> \
        <id_length> <optional: insertion_length (default is 1)>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("bampath_in",
                         help="path to the input bam file")
    parser.add_argument("bampath_out",
                         help="path to the output bam file")
    parser.add_argument("genome_path",
                         help = "Path to the .fasta genome used in the \
                                 alignment")
    parser.add_argument("genome_index_path",
                         help = "path to the .fai produced by samtools \
                                 faidx from the genome .fasta")
    parser.add_argument("id_length",
                        help="length of the barcode in the bam ID line")
    parser.add_argument("insertion_length",
                        help = "The length of the transposase insertion. \
                                Default is 1 since pysam uses 0 based half \
                                    open intervals. [0,1) returns 1 base",
                        default = '1')

    return parser.parse_args(args)


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check the input to \
        add_read_group_to_bam.py -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check that the input bam \
            exists, is sorted, and is indexed -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def add_read_group_and_tags(bampath_in, bampath_out, genome_path,
                   genome_index_path, id_length, insertion_length):
    """
    :param bampath_in: path to the sorted, indexed bam file with UMI IDs
                       in the ID entry
    :param bampath_out: the same bam, but with the UMI ID added as a RG: tag
    :param genome_path:
    :param genome_index_path:
    :param id_length: length of the barcode which is appended, after an _,
                      to the read id

    Look at each read, extract the ID, add it as a tag (tag prefix RG) and
    write out to another bam. Output is the same bam, but with RG:<index>
    added to each line
    """

    # TODO: Right now, this loops over the file twice to create the header,
    #       and then to add the tags to each read. Need to figure out how to
    #       add the header without re-writing the whole bam so that only loops
    #       one time.

    # open files ---------------------------------------------------------------
    # including the index allows random access on disc
    genome = pysam.FastaFile(genome_path, genome_index_path)
    # open the input bam
    input_bamfile = pysam.AlignmentFile(bampath_in, "rb")

    # Get the set of unique barcodes in the input_bamfile ----------------------
    # extract current header
    new_header = input_bamfile.header.copy().to_dict()
    # instantiate a set object
    barcode_set = set()

    # loop through the reads in the input bam, extract unique barcodes, add to
    # barcode set
    for read in input_bamfile.fetch():
        # extract the ID from the ID line
        read_dict = read.to_dict()
        barcode = read_dict['name'][-id_length:]
        barcode_set.add(barcode)
    # Create new read group header. Note: this is used below in the tagged_bam
    new_header['RG'] = [{'ID': barcode} for barcode in barcode_set]

    # Add RG, XI and XZ tags to each line in the bam file ----------------------
    input_bamfile = pysam.AlignmentFile(bampath_in, "rb")
    # open a file to which to write. Note: the new header
    tagged_bam = pysam.AlignmentFile(bampath_out, "wb", header = new_header)
    for read in input_bamfile.fetch():

        # Extract XI and XZ tags -----------------------------------------------
        region_dict = dict()
        # (using the bitwise operator) check if the read is unmapped,
        # if so, set the region_dict start and end to *, indicating that there is
        # no alignment, and so there is no start and end region for the alignment
        if read.flag & 0x4:
            region_dict['start'] = "*"
            region_dict['end']   = "*"
        # if the bit flag 0x10 is set, the read reverse strand. Handle accordingly
        elif read.flag & 0x10:
            # A cigartuple looks like [(0,4), (2,2), (1,6),..,(4,68)] if read
            # is reverse complement. If it is forward, it would have the (4,68),
            # in this case, in the first position.
            # The first entry in the tuple is the cigar operation and the
            # second is the length. Note that pysam does order the tuples in the
            # reverse order from the sam cigar specs, so cigar 30M would be
            # (0,30). 4 is cigar S or BAM_CSOFT_CLIP. The list operation below
            # extracts the length of cigar operation 4 and returns a integer.
            # if 4 DNE, then soft_clip_length is 0.
            try:
                soft_clip_length = read.cigartuples[-1][1] \
                    if read.cigartuples[-1][0] == 4 \
                    else 0
            except TypeError:
                sys.exit("In bamfile %s, for read %s, cigar string %s \
                          is not parse-able" %
                (bampath_in, read.query_name, read.cigartuples))
            # The insertion point is at the end of the alignment
            insert = read.reference_end
            if insert is None:
                raise TypeError("failure to get read.reference_end from \
                                 read %s in bamfile %s"
                                %(read.query_name, bampath_in))
            # adjust insert for soft clipping
            # prevent from extending past end of chr
            insert = min(genome.get_reference_length(read.reference_name),
                         insert + soft_clip_length)
            # The span of the n bases preceding the insertion are:
            # [insert, insert + n)
            region_dict['start'] = insert
            # make sure the insert doesn't extend beyond the chr
            region_dict['end'] = min(genome.get_reference_length(read.reference_name),
                                     insert + insertion_length)
        # else, Read is in the forward orientation
        else:
            # see if clause for lengthy explanation. This examines the first
            # operation in the cigar string. If it is a soft clip (code 4),
            # the length of the soft clipping is stored. Else there is 0 soft
            # clipping
            try:
                soft_clip_length = read.cigartuples[0][1] \
                    if read.cigartuples[0][0] == 4 \
                    else 0
            except TypeError:
                sys.exit("In bamfile %s, for read %s, cigar string %s is \
                          not parse-able" %
                (bampath_in, read.query_name, read.cigartuples))
            # extract insert position
            insert = read.reference_start
            if insert is None:
                raise TypeError("failure to get read.reference_start from read \
                                 %s in bamfile %s"
                                %(read.query_name, bampath_in))
            # adjust insert for soft clipping
            # prevent from extending past end of chr
            insert = max(0,insert - soft_clip_length)
            # The span of the n bases preceding the insertions are:
            # [insert - n, insert). Take 'max' to prevent the 'start' value
            # from extending past the end of the chromosome
            region_dict['start'] = max(0,insert - insertion_length)
            region_dict['end']   = insert

        # Create the tag strings -----------------------------------------------
        tag_dict = dict()
        tag_dict['RG'] = read.query_name[-id_length:]
        tag_dict['XI'] = region_dict['start']
        # if the read was unmapped, set XZ to "*"
        if region_dict['start'] == "*" and region_dict['end'] == "*":
            tag_dict['XZ'] = "*"
        else:
            try:
                tag_dict['XZ'] = genome.fetch(read.reference_name,
                                    region_dict['start'],
                                    region_dict['end']).upper()
            except ValueError:
                sys.exit("In bamfile %s, for read %s, insert region %s:%s-%s \
                          is out of bounds" %
                        (bampath_in,
                        read.query_name,
                        read.reference_name,
                        region_dict['start'],
                        region_dict['end']))

        # Set tags -------------------------------------------------------------
        # the list comprehension outputs [None, None, None, ...]. Out catches
        # this so it isn't printed to std out
        # how else to do this?
        out = [read.set_tag(tag,tag_str) for tag, tag_str in tag_dict.items()]
        # Write to file --------------------------------------------------------
        tagged_bam.write(read)

    # Close files --------------------------------------------------------------
    genome.close()
    tagged_bam.close()
    input_bamfile.close()

    # Re-index the output ------------------------------------------------------
    # This is only here only to prevent the warning message:
    # bam timestamp and read timestamp are different
    # that you sometimes get when the bam is modified after the index is created
    pysam.index(bampath_out)

def main(args=None):
    args = parse_args(args)

    # Check inputs
    input_path_list = [args.bampath_in,
                       args.genome_path,
                       args.genome_index_path]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" %input_path)

    # loop over the reads in the bam file and add the read group (header and tag)
    # and the XI and XZ tags
    add_read_group_and_tags(args.bampath_in,
                            args.bampath_out,
                            args.genome_path,
                            args.genome_index_path,
                            int(args.id_length),
                            int(args.insertion_length))

    sys.exit(0)


if __name__ == "__main__":
    sys.exit(main())

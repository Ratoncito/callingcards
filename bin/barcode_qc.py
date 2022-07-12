#!/usr/bin/env python

# TODO better docstring, in particular, describe the 6x ? bed file format

"""
Given a 6x? bedfile detailing calling cards read insertions, split the
    name column into barcode components, tally in various ways (QC) and output
    a filtered bedfile with only those reads which pass barcode match standards.

    outputs the following files in the $PWD:
      - <input_filename>_<barcode_component>_tally.tsv
          for each barcode component, a tally of the unique elements in that
          component
      - <input_filename>_<barcode_component>_ppm.tsv
          a position probability matrix describing
          the percent probability of a given base in a given position given the
          multi-sequence alignment given a barcode component
      - <input_filename>_bc_fltr.bed
          a bed file in the same format as the input bedfile, but with a row
          filter applied such that only those records which conform to barcode
          AND insert_sequence specifications remain
"""
# standard library
import os
import sys
import argparse
import json
# outside dependencies
import pandas as pd


def parse_args(args=None):
    Description = "Examine hop barcodes, output some QC metrics, and a bed file \
        which matches the columns in the input bed, but is row filtered for only \
            those reads which meet barcode expectations"
    Epilog = "Example usage: python mammal_barcode_qc.py \
        <cc_format_bedfile.bed> <barcode_details.json>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("bed_path",
                         help="path to the input bam file")
    parser.add_argument("barcode_details",
                         help="A json which describes components of the barcode \
                            which occurs in the name column of the calling cards \
                                bed file format")
    return parser.parse_args(args)

def barcode_qc(bed_df,
               barcode_details,
               fltr_bed_output_name):
    """parse the calling cards bed file name column into components and perform
    This expects a dataframe with AT LEAST the columns (there may be more):
    ['chrom','chromStart','chromEnd','barcode', 'strand', 'insert_seq']

    :param bed_df:
    :param barcode_details:
    :param fltr_bed_output_name:

    """

    # fields by which to group and count transpositions. Not that in the case
    # of pooled barcodes, this occurs after splitting into individual barcode
    # tables
    coordinate_grouping_fields = ['chrom','chromStart','chromEnd','strand']
    # ccf fields and order
    ccf_col_order = ['chrom', 'chromStart', 'chromEnd', 'reads', 'strand']

    # Using the keys of the barcode_components dict, make a dictionary
    # with structure {component1: [], component2:[], ...} where the keys are
    # barcode components, and the values are empty lists
    barcode_vector_dict = { key:[] for key in barcode_details['components'].keys()}

    # iterate over the barcodes in the 'name' column of the cc bed file
    # and parse the barcode into the barcode components
    for barcode in bed_df.loc[:,'barcode']:
        [barcode_vector_dict[k].append(barcode[v[0]:v[1]]) \
        for k,v in barcode_details['indicies'].items()]

    # TODO rename the function below to reflect the fact that both barcodes
    # and insertion sites may be handled

    # create position probability matricies of the barcode components
    barcode_components_ppm(barcode_vector_dict)
    barcode_components_ppm({'insert_seq': bed_df.loc[:,'insert_seq']})

    # transform barcode_vector_dict to a dataframe
    barcode_component_df = pd.DataFrame(data=barcode_vector_dict)

    # TODO rename the function below to reflect the fact that both barcodes
    # and insertion sites may be handled

    # create tables which describe the variety of barcode components (eg,
    # group by the srt column and count how many of each unique srt seq
    # there is)
    counts_of_barcode_variety(barcode_component_df)
    counts_of_barcode_variety(pd.DataFrame({'insert_seq': bed_df.loc[:,'insert_seq']}))

    # In the following section, create a vector to filter out
    # barcodes which do not meet barcode expectation OR!! the insert sequence
    # expectation

    # initialize a boolean vector with length equal to the number of rows
    # in the bed dataframe
    barcode_fltr_vector = [True]*len(bed_df)
    # iterate over the parsed barcode vector dataframe
    for index,row in barcode_component_df.iterrows():
        # initialize variable keep_index to True
        keep_index = True
        # iterate over the barcode components
        for barcode_component in barcode_vector_dict.keys():
            # if the current barcode component is not in the list of expected
            # values for that component
            if row[barcode_component] not in \
                barcode_details['components'][barcode_component]:
                keep_index = False
            # if there is a specific insert seq on which to filter, also
            # exclude reads that don't match it
            try:
                if not bed_df.loc[index,'insert_seq'] in barcode_details['insert_seq']:
                    keep_index = False
            except KeyError:
                pass
        # if keep_index if false,
        if not keep_index:
            barcode_fltr_vector[index] = False

    # reduce the cc bed file to only those rows which have expected
    # barcodes, and write out
    col_list = list(bed_df.columns)
    try:
        # try to read the field tf_map into a dataframe. if there is a keyerror,
        # handle in the except block
        barcode_map = pd.DataFrame(barcode_details['tf_map'])
        # group by TF. If there is more than one field to group on, add to the
        # list
        barcode_grouping_fields = ["TF"]
        # row filter, join, transform mismatches to 'other'
        barcode_to_tf_df = pd.merge(barcode_component_df,barcode_map,
                        how = 'left',
                        on = list(barcode_map.columns[barcode_map.columns != "TF"]))
        col_list.append("TF")
        grouped_bed = pd.concat([bed_df, barcode_to_tf_df], axis=1)[col_list]\
                        .fillna({'TF':"other"})\
                        .groupby(barcode_grouping_fields)

        # output each grouped sheet with the name of the group in the filename
        for barcode_grouping_fields,fltr_df in grouped_bed:
            group_name = barcode_grouping_fields \
                if type(barcode_grouping_fields) is str \
                else "_".join(barcode_grouping_fields)
            fltr_df\
                .value_counts(coordinate_grouping_fields+['TF'])\
                .to_frame()\
                .sort_values(coordinate_grouping_fields)\
                .reset_index()\
                .rename(columns = {0:'reads'})\
                .reindex(ccf_col_order)\
                .to_csv(fltr_bed_output_name+"_"+group_name + "_bc_fltr.ccf",
                            sep = "\t",
                            header = None,
                            index = False)
    # if the tf_map slot is not in barcode_details, then just write out a single
    # table
    except KeyError:
        bed_df[barcode_fltr_vector]\
            .value_counts(coordinate_grouping_fields)\
            .to_frame()\
            .sort_values(coordinate_grouping_fields)\
            .reset_index()\
            .rename(columns = {0:'reads'})\
            .reindex(ccf_col_order)\
            .to_csv(fltr_bed_output_name+"_bc_fltr.ccf",
                                        sep = "\t",
                                        header = None,
                                        index = False)

def frequency_matrix(dna_list):
    """Tally the frequency of each base by position

    :param dna_list: a list of dna sequences of the same length

    this is hand coded because there isn't an easy way to do this in BioPython?
    In R, it is like so --

    library(Biostrings)
    start <- Sys.time()
    x = DNAStringSet(df$srt)
    freq_mat = consensusMatrix(x, baseOnly = TRUE)
    print(freq_mat)
    end <- Sys.time()
    print(paste0("Runtime: ", end-start))

    and takes .004 seconds. Could not find an easy function to do the same
    with biopython -- if one exists consider replace the function below with
    that
    """
    # CITE: https://hplgit.github.io/bioinf-py/doc/pub/html/main_bioinf.html
    n = len(dna_list[0])
    A = [0]*n
    T = [0]*n
    G = [0]*n
    C = [0]*n
    other = [0]*n
    for dna in dna_list:
        try:
            for index, base in enumerate(dna):
                if base == 'A':
                    A[index] += 1
                elif base == 'C':
                    C[index] += 1
                elif base == 'G':
                    G[index] += 1
                elif base == 'T':
                    T[index] += 1
                else:
                    other[index] +=1
        except IndexError:
            print("dna: %s; index: %s; n: %s" %(base, index, n))


    # make position frequency matrix
    df = pd.DataFrame({'A':A,'C':C,'G':G,'T':T,'other':other}).transpose()
    # transform to position probability matrix
    df = (df/df.sum(axis=0)*100)

    return(df)

def barcode_components_ppm(barcode_vector_dict):
    """Convert each barcode component into a position probability matrix.
    Writ the ppm for each component to file.

    :param barcode_vector_dict: a dictionary where the keys are the names of the
    components of the barcode, and the values are lists of those components
    which have been parsed out of the corresponding barcode
    """
    freq_table_dict = {k:frequency_matrix(v) for k,v in barcode_vector_dict.items()}

    for k,v in freq_table_dict.items():
        v.to_csv(k+"_ppm.tsv", sep="\t")

def counts_of_barcode_variety(barcode_component_df):
    """Group and tally by each of the barcode components, and write the output
    to file

    :param barcode_component_df: a dataframe which has columns corresponding to
    the expected components of the barcode in the name column of the calling
    cards bed file
    """
    # iterate over columns in the input df, group by that column, tally
    # and then write to file
    for barcode_component in barcode_component_df.columns:
        tally_series = barcode_component_df\
                            .groupby(by=barcode_component)\
                            .size()\
                            .sort_values(ascending=False)
        tally_series.name = "tally"
        tally_series.to_csv(barcode_component+"_tally.tsv", sep = "\t")

def main(args=None):
    args = parse_args(args)

    # Check inputs
    for input_path in [args.bed_path,
                       args.barcode_details]:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" %input_path)

    bed_colnames = ['chrom','chromStart','chromEnd','barcode','aln_mapq',
                    'strand', 'reads', 'insert_seq', 'aln_flag']

    bed_df = pd.read_csv(args.bed_path, sep = "\t", names=bed_colnames)

    with open(args.barcode_details) as f1:
        barcode_details = json.load(f1)

    if not list(barcode_details['components'].keys()).sort() == \
        list(barcode_details['indicies'].keys()).sort():
        raise ValueError('The keys in the barcode_component json and the ' +
            'barcode_component_indicies json are not the same')

    barcode_components_length_check_dict = { k:v[1]-v[0] for k,v in \
                                              barcode_details['indicies'].items()}

    for k,v in barcode_details['components'].items():
        for comp in v:
            if not len(comp) == barcode_components_length_check_dict[k]:
                raise ValueError("There exists a component in " +
                "barcode_components which is not the length described by the " +
                "barcode_details['indicies']: %s, %s"
                %(comp, barcode_components_length_check_dict[k]))

    fltr_bed_output_name = os.path.splitext(os.path.basename(args.bed_path))[0]

    # parse out the name column of the bed file into its components and
    # perform some QC/filtering
    barcode_qc(bed_df,
               barcode_details,
               fltr_bed_output_name)

    sys.exit(0)


if __name__ == "__main__":
    sys.exit(main())

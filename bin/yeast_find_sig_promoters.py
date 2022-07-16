#/usr/bin/env python

"""

:raises ValueError: If a sql command does not work, a value error is raised
:return:
:rtype: int
"""

# standard lib
from inspect import Attribute
import sys
import argparse
from os.path import exists
import logging
# third party
import pandas as pd
import sqlite3

logging.basicConfig(
    level=logging.CRITICAL,
    format='[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
    datefmt='%H:%M:%S',
    stream=sys.stdout
)

def create_db_table(con, df, tbl_details_dict):
    """create a table in the calling cards database

    :param con: psycopg connection to a sqlite database
    :type con: psycopg connection object
    :param df: a pandas dataframe to add to the database
    :type df: pandas dataframe
    :param tablename: name of the table in the database
    :type tablename: str
    :param table_colnames: names of the fields for this table in the database
    :type table_colnames: list/set of str
    :param index_name: name of hte index on the given columns of this table
    :type index_name: str
    :param index_col_string: field names to index
    :type index_col_string: list/set of str
    """
    # assign a cursor to the database at cur
    cur = con.cursor()
    # drop the table if it exists
    drop_sql = "DROP TABLE IF EXISTS %s" %tbl_details_dict.get('tablename')
    try:
        cur.execute(drop_sql)
    except sqlite3.OperationalError as e:
        msg = 'Could not drop table with sql: %s. Error: %s' %(drop_sql, e)
        logging.critical(msg+" ", exc_info=(sys.exc_info()))
        raise
    # create the table
    create_tbl_sql = "CREATE TABLE %s (%s, PRIMARY KEY(%s AUTOINCREMENT));" \
        %(tbl_details_dict.get('tablename'), ",".join(["id INTEGER"] +
          tbl_details_dict.get('colnames')), "id")
    try:
        cur.execute(create_tbl_sql)
    except sqlite3.OperationalError as e:
        cur.execute(drop_sql)
        msg = 'Could not create table with sql:\n"%s".\nTable has been dropped.\n \
            Error: %s.' %(create_tbl_sql, e)
        logging.critical(msg+" ", exc_info=(sys.exc_info()))
        raise

    # create index on table
    index_sql = 'CREATE INDEX %s ON %s %s' %(tbl_details_dict.get('index_name'),
                                             tbl_details_dict.get('tablename'),
                                          tbl_details_dict.get('index_col_string'))
    try:
        cur.execute(index_sql)
    except sqlite3.OperationalError as e:
        cur.execute(drop_sql)
        msg = "could not create index with sql:\n%s.\n"%index_sql +\
        "Table has been dropped.\nError: %s:" %e
        logging.critical(msg+" ", exc_info=(sys.exc_info()))
        raise

    df.to_sql(tbl_details_dict.get('tablename'),
              con=con,
              if_exists='append',
              index=False)



def create_hop_view(con, tablename, view_name):
    """Create view where hops are aggregated by chrom, pos and strand

    :param con: psycopg connection to a sqlite database
    :type con: psycopg connection object
    :param tablename: name of the table in the database
    :type tablename: str

    :raises ValueError: if database execution does not work
    """

    cur = con.cursor()

    drop_view_sql = 'DROP VIEW IF EXISTS "main"."%s";' %view_name

    cur.execute(drop_view_sql)

    create_view_sql = \
        'CREATE VIEW "%s" AS ' %view_name +\
        'SELECT chrom, chromStart,chromEnd, strand COUNT(*) as reads ' +\
        'FROM( SELECT p.chrom AS chrom, p.chromStart AS chromStart, ' +\
        'p.chromEnd AS chromEnd, p.strand as strand, '+\
        'x.chromStart as insert_pos ' +\
        'FROM %s as x ' %tablename +\
        'LEFT JOIN promoters as p ' +\
        'WHERE (x.chrom = p.chrom AND ' +\
        'x.strand = p.strand AND ' +\
        'x.chromStart BETWEEN P.chromStart AND p.chromEnd)) as t ' +\
        'GROUP BY chrom, chromStart, chromEnd, strand ' +\
        'ORDER BY chrom, chromStart ASC, chromEnd ASC, strand;'

    try:
        cur.execute(create_view_sql)
    except sqlite3.OperationalError as e:
        cur.execute(drop_view_sql)
        msg = "could not create view with sql:\n%s.\n"%create_view_sql +\
        "View has been dropped.\nError: %s:" %e
        logging.critical(msg+" ", exc_info=(sys.exc_info()))
        raise

def standardize_chr_names(df,chr_map_df, standard_chr_format):
    """_summary_

    :param df: _description_
    :type df: _type_
    :param chr_map_df: _description_
    :type chr_map_df: _type_
    :param standardized_chr_format: _description_
    :type standard_chr_format: _type_
    :return: _description_
    :rtype: _type_
    """

    # determine which name format the dataframe currently uses
    curr_chrom_format = -1
    for chr_format in chr_map_df.columns:
        if sum([True if x in chr_map_df[chr_format].unique() else False \
            for x in df['chrom'].unique()]) == len(df['chrom'].unique()):
            curr_chrom_format = chr_format
            break
    # raise an error if a format did not match
    if curr_chrom_format == -1:
        AttributeError("Chromosome names are not recognized in a recognized format")

    # if the current names are refseq, return
    if curr_chrom_format == standard_chr_format:
        return df
    # else, use the chr_map_df to map to the standardized_chr_format convention
    else:
        return pd.merge(df, chr_map_df[[curr_chrom_format,'refseq']],
                        how = 'left',
                        left_on = 'chrom',
                        right_on = curr_chrom_format)\
                        .drop(['chrom', curr_chrom_format], axis=1)\
                        .rename(columns = {'refseq':'chrom'})\
                        [df.columns]

def parse_args(args=None):
    Description = "create a sqlite database with background and expression data " +\
    "and calculate experimental significance. Optionally save the sqlite database."
    Epilog = "Example usage: yeast_find_sig_promoters.py \
        <input.bam> <output.bam> <genome.fasta> <genome.fasta.fai> \
        <id_length> <optional: insertion_length (default is 1)>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("promoter_bed_path",
                         help="path to to a bed file format promoter file")
    parser.add_argument("background_ccf_path",
                         help="path` to the background data table")
    parser.add_argument("experimental_ccf_path",
                         help="path to the experimental data table")
    parser.add_argument("chr_map",
                         help="a csv which maps between various "+\
                            "chromosome naming conventions, eg ucsc, ignomes " +\
                                "and refseq")
    parser.add_argument("standard_chr_format",
                        help="this must be one of the column names in "+\
                            "the chr_map. All chromosome identifiers will be "+\
                                "translated to this naming convention. "+\
                                    "default is 'refseq'", default="refseq")

    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    PROMOTER_BED_COLNAMES = ["chrom", "chromStart", "chromEnd", "name", "score", "strand"]
    CCF_COLNAMES = ['chrom', 'chromStart', 'chromEnd', 'reads', 'strand']
    INDEX_COL_STRING = '("chrom", "chromStart" ASC, "chromEnd" ASC, "strand");'

    con = sqlite3.connect(':memory:')
    promoter_dict = {
        'colnames': PROMOTER_BED_COLNAMES,
        'tablename': "promoters",
        'index_name': "promoter_index",
        'index_col_string': INDEX_COL_STRING
    }

    background_dict = {
        'colnames': CCF_COLNAMES,
        'tablename': "background",
        'hop_view': "background_hops",
        'index_name': "background_index",
        'index_col_string': INDEX_COL_STRING
    }

    experimental_dict = {
        'colnames': CCF_COLNAMES,
        'tablename': "experiment",
        'hop_view': "experiment_hops",
        'index_name': "experiment_index",
        'index_col_string': INDEX_COL_STRING
    }

    if not exists(args.chr_map):
        raise FileExistsError('File Not Found: %s' %args.chr_map)
    else:
        chr_map_df = pd.read_csv(args.chr_map)
        if not args.standard_chr_format in chr_map_df.columns:
            raise AttributeError('Standard chr format not in chr_map columns')

    # create promoter, background and experiment tables
    if not exists(args.promoter_table_path):
        raise FileExistsError("File Not Found: %s" %args.promoter_table_path)
    else:
        promoter_df = pd.read_csv(args.promoter_table_path,
                                  sep = "\t",
                                  names = promoter_dict.get('colnames'))
        promoter_df = standardize_chr_names(promoter_df,
                                            chr_map_df,
                                            args.standard_chr_format)
        create_db_table(con, promoter_df, promoter_dict)
        del promoter_df

    if not exists(args.background_data_path):
        raise FileExistsError("File Not Found: %s" %args.background_data_path)
    else:
        background_df = pd.read_csv(args.background_data_path,
                                     sep = "\t",
                                     names = background_dict.get('colnames'))
        background_df = standardize_chr_names(background_df,
                                            chr_map_df,
                                            args.standard_chr_format)
        create_db_table(con, background_df, background_dict)
        del background_df

    if not exists(args.experimental_data_path):
        raise FileExistsError("File Not Found: %s" %args.experimental_data_path)
    else:
        experimental_df = pd.read_csv(args.experimental_data_path,
                                  sep = "\t",
                                  names = experimental_dict.get('colnames')
        experimental_df = = standardize_chr_names(experimental_df,
                                                  chr_map_df,
                                                  args.standard_chr_format)
        create_db_table(con, experimental_df, experimental_dict)
        del experiment_df

    # create background hop view
    create_hop_view(con,
                    background_dict.tablename,
                    background_dict.view_name)
    # create experimental hop view
    create_hop_view(con,
                    experimental_dict.tablename,
                    experimental_dict.view_name)

    bg_df   = pd.read_sql_query("SELECT * FROM %s"
                                %background_dict.get('view_name'), con)
    expr_df = pd.read_sql_query("SELECT * FROM %s"
                                %experimental_dict.get('view_name'), con)

    # get nrows of the tables
    total_background_hops = nrow(bg_df)
    total_exp_hops        = nrow(expr_df)

    quant_df = bg_df %>%
        left_join(expr_df, by = c('chrom', 'chromStart', 'chromEnd')) %>%
        dplyr::rename(bg_hops = hops.x, exp_hops = hops.y) %>%
        filter(!is.na(exp_hops)) %>%
        mutate(poisson_pval =
           1-ppois(exp_hops, (bg_hops*(total_exp_hops/total_background_hops))))

    return quant_df

if __name__ == "__main__":
    sys.exit(main())



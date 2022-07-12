#/usr/bin/env python

"""

:raises ValueError: If a sql command does not work, a value error is raised
:return:
:rtype: int
"""

# standard lib
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

    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    CCF_COLNAMES = ['chrom', 'chromStart', 'chromEnd', 'reads', 'strand']
    INDEX_COL_STRING = '("chrom", "chromStart" ASC, "chromEnd" ASC, "strand");'

    con = sqlite3.connect(':memory:')
    promoter_dict = {
        'colnames': ["chrom", "chromStart", "chromEnd", "name", "score", "strand"],
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

    # create promoter, background and experiment tables
    if not exists(args.promoter_table_path):
        raise FileExistsError("File Not Found: %s" %args.promoter_table_path)
    else:
        promoter_df = pd.read_csv(args.promoter_table_path,
                                  sep = "\t",
                                  names = promoter_dict.get('colnames'))
        try:
            create_db_table(con, promoter_df, promoter_dict)
        except SystemExit(1):
            msg = 'Error'
            logging.critical(msg+" ", exc_info=(sys.exec_info()))
        del promoter_df

    if not exists(args.background_data_path):
        raise FileExistsError("File Not Found: %s" %args.background_data_path)
    else:
        background_df = pd.read_csv(args.background_data_path,
                                     sep = "\t",
                                     names = background_dict.get('colnames'))
        create_db_table(con, background_df, promoter_dict)
        del background_df

    if not exists(args.experimental_data_path):
        raise FileExistsError("File Not Found: %s" %args.experimental_data_path)
    else:
        experimental_df = pd.read_csv(args.experimental_data_path,
                                  sep = "\t",
                                  names = experimental_dict.get('colnames'))
        create_db_table(con, experimental_df, experimental_dict)
        del experiment_df

    # create hop view
    create_hop_view(con, background_dict.tablename, background_dict.view_name)
    create_hop_view(con, experimental_dict.tablename, experimental_dict.view_name)

    bg_df    = pd.read_sql_query("SELECT * FROM %s"
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



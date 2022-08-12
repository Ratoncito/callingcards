#!/usr/bin/env python

"""Calculate poisson and hypergeomtric p-values for a calling cards experiment.
Given the experiment qbed file, background qbed file, promoter definitions,
a csv which maps various chromosome naming conventions to one another (eg
refseq to ucsc), and a choice of which chromosome naming convention to use,
aggregate hops over the promoter regions and calculate the poisson and
hypergeometric pvalues. Manipulations of the data frames are largely performed
within sqlite, which may be done either in memory or saved to disc. Choose to
save the sqlite database to disc for easy re-definition of promoter regions.

Note: currently no 'lax' score (lax just 'relaxes' the frame in which hops are
  counted by some amount. In the yeast 3.0 scripts, it is 200bp) is not
  calculated.

Note: the strand of the hop is not considered for the background or
  experimental data -- if a hop occurs on either strand in a given promoter
  region, it is counted.

written by: chase mateusiak, chasem@wustl.edu, 202207

:raises AttributeError: Various custom AttributeErrors are raised
:raises FileExistsError: All cmd line param filepaths are checked for
    existence
:raises KeyError: If a key in one of the dictionaries describing attributes
    related to the promoter, background or experimental table DNE

:output: a statistics data frame which will have poisson_pval
    and hypergeom_pval columns with filename
    <experimental_qbed_basename>_stats.csv. Optionally, the sqlite database
    may be saved to disc. A cmd line argument directs the location.
"""

# standard lib
import sqlite3
from inspect import Attribute
import sys
import argparse
from os.path import exists,basename,splitext
import logging
# third party
import pandas as pd
from scipy.stats import poisson as pois
from scipy.stats import hypergeom

logging.basicConfig(
    level=logging.CRITICAL,
    format='[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
    datefmt='%H:%M:%S',
    stream=sys.stdout
)

def create_db_table(con, df, tbl_details_dict):
    """Create a table in the calling cards database

    :param con: sqlite3 connection object
    :type con: sqlite3 connection
    :param df: table to add to the database
    :type df: pandas DataFrame
    :param tbl_details_dict: a dictionary which has at least the attributes
    :type tbl_details_dict: dictionary
    """
    # assign a cursor to the database at cur
    cur = con.cursor()
    # drop the table if it exists
    drop_sql = "DROP TABLE IF EXISTS %s" %tbl_details_dict['tablename']
    try:
        cur.execute(drop_sql)
    except sqlite3.OperationalError as e:
        msg = 'Could not drop table with sql: %s. Error: %s' %(drop_sql, e)
        logging.critical(msg+" ", exc_info=(sys.exc_info()))
        raise
    # create the table
    create_tbl_sql = "CREATE TABLE %s (%s, PRIMARY KEY(%s AUTOINCREMENT));" \
        %(tbl_details_dict['tablename'], ",".join(["id INTEGER"] +
          list(df.columns)), "id")
    try:
        cur.execute(create_tbl_sql)
    except sqlite3.OperationalError as e:
        cur.execute(drop_sql)
        msg = 'Could not create table with sql:\n"%s".\nTable has been dropped.\n \
            Error: %s.' %(create_tbl_sql, e)
        logging.critical(msg+" ", exc_info=(sys.exc_info()))
        raise

    # create index on table
    index_sql = 'CREATE INDEX %s ON %s %s' %(tbl_details_dict['index_name'],
                                             tbl_details_dict['tablename'],
                                          tbl_details_dict['index_col_string'])
    try:
        cur.execute(index_sql)
    except sqlite3.OperationalError as e:
        cur.execute(drop_sql)
        msg = "could not create index with sql:\n%s.\n"%index_sql +\
        "Table has been dropped.\nError: %s:" %e
        logging.critical(msg+" ", exc_info=(sys.exc_info()))
        raise

    df.to_sql(tbl_details_dict['tablename'],
              con=con,
              if_exists='append',
              index=False)

def create_hop_view(con, tablename, view_name):
    """Create view where hops are aggregated by chrom, pos and strand

    :param con: sqlite3 connection to database
    :type con: sqlie3 connection object
    :param tablename: name of the table from which to create the view
    :type tablename: str
    :param view_name: name of the view
    :type view_name: str
    """

    cur = con.cursor()

    drop_view_sql = 'DROP VIEW IF EXISTS "main"."%s";' %view_name

    cur.execute(drop_view_sql)

    create_view_sql = 'CREATE VIEW %s AS '%view_name +\
               'SELECT chrom, chromStart,chromEnd, SUM(reads) as reads '+\
               'FROM(SELECT p.chrom AS chrom, p.chromStart AS chromStart, '+\
                            'p.chromEnd AS chromEnd, p.strand AS strand, '+\
                            'x.chromStart AS insert_pos, x.reads AS reads '+\
	                 'FROM %s as x '%tablename +\
	                 'LEFT JOIN promoters as p '+\
	                 'WHERE (x.chrom = p.chrom AND '+\
                        'x.chromStart BETWEEN p.chromStart AND p.chromEnd)) as t '+\
                'GROUP BY chrom, chromStart, chromEnd '+\
                'ORDER BY chrom, chromStart ASC, chromEnd ASC;'

    try:
        cur.execute(create_view_sql)
    except sqlite3.OperationalError as e:
        cur.execute(drop_view_sql)
        msg = "could not create view with sql:\n%s.\n"%create_view_sql +\
        "View has been dropped.\nError: %s:" %e
        logging.critical(msg+" ", exc_info=(sys.exc_info()))
        raise

def standardize_chr_names(df, chr_map_df, standard_chr_format):
    """translate the chrom column of a given df to a different chromosome naming
    convention. If the current chrom column naming convention matches the
    desired one, no change is made.

    :param df: dataframe with at least column chrom. all entries in chrom field
      must be in one of the columns of chr_map_df
    :type df: pandas DataFrame
    :param chr_map_df: fields correspond to chromosome naming conventions, eg
      maybe refseq,ucsc,ensembl
    :type chr_map_df: pandas DataFrame
    :param standardized_chr_format: a column of chr_map_df to which the chrom
      column of the input df will be translated
    :type standard_chr_format: pandas Dataframe

    :raises AttributeError: raised if the entries of the chrom column of the
      input dataframe are not entirely described in chr_map_df

    :return: the input dataframe with the chrom column translated to the
      standard_chr_format via the chr_map_df
    :rtype: pandas DataFrame
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
        raise AttributeError("Chromosome names are not "+\
            "recognized in a recognized format. Unique chr names which cause "+\
                " error are: %s" %",".join(df['chrom'].unique()))

    # if the current names are refseq, return
    if curr_chrom_format == standard_chr_format:
        return df
    # else, use the chr_map_df to map to the standardized_chr_format convention
    else:
        return pd.merge(df,
                        chr_map_df.reindex(columns=[curr_chrom_format,
                                                    standard_chr_format]),
                        how = 'left',
                        left_on = 'chrom',
                        right_on = curr_chrom_format)\
                        .drop(['chrom', curr_chrom_format], axis=1)\
                        .rename(columns = {standard_chr_format:'chrom'})\
                        [df.columns]

def parse_args(args=None):
    Description = "create a sqlite database with background and expression data " +\
    "and calculate experimental significance. Optionally save the sqlite database."
    Epilog = "Example usage: yeast_find_sig_promoters.py \
        <input.bam> <output.bam> <genome.fasta> <genome.fasta.fai> \
        <id_length> <optional: insertion_length (default is 1)>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("--promoter_bed_path", "-p",
                        type=str,
                        help="path to to a bed file format promoter file")
    parser.add_argument("--background_qbed_path", "-b",
                        type=str,
                        help="path` to the background data table")
    parser.add_argument("--experimental_qbed_path", "-e",
                        type=str,
                        help="path to the experimental data table")
    parser.add_argument("--chr_map", "-c",
                        type=str,
                        help="a csv which maps between various "+\
                            "chromosome naming conventions, eg ucsc, ignomes " +\
                                "and refseq")
    parser.add_argument("--standard_chr_format", "-s",
                        type=str,
                        help="this must be one of the column names in "+\
                            "the chr_map. All chromosome identifiers will be "+\
                                "translated to this naming convention. "+\
                                    "default is 'refseq'",
                        default="refseq")
    parser.add_argument("--sqlite_db_out", "-d",
                        type=str,
                        help="path to sqlite database output. "+\
                            "Default is ':memory:' for an in memory DB. " +\
                                "For instance, 'my_db.sqlite' would write to "+\
                                    "$PWD/my_db.sqlite.",
                        default=":memory:")
    parser.add_argument("--poisson_pseudocount", "-x",
                        type=float,
                        help="pseudocount to add to the poisson p-value "+\
                            "calculation. Default = 0.2",
                        default=0.2)

    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    PROMOTER_BED_COLNAMES = ["chrom", "chromStart", "chromEnd", "name", "score", "strand"]
    qbed_COLNAMES = ['chrom', 'chromStart', 'chromEnd', 'reads', 'strand']
    INDEX_COL_STRING = '("chrom", "chromStart" ASC, "chromEnd" ASC, "strand");'

    if args.sqlite_db_out != ":memory:" and\
         splitext(args.sqlite_db_out)[1] != '.sqlite':
        raise ValueError("sqlite must either be in memory -- ':memory:' -- or "+\
            "have the extension '.sqlite'.")

    con = sqlite3.connect(args.sqlite_db_out)

    promoter_dict = {
        'tablename': "promoters",
        'index_name': "promoter_index",
        'index_col_string': INDEX_COL_STRING
    }

    background_dict = {
        'tablename': "background",
        'hop_view': "background_hops",
        'index_name': "background_index",
        'index_col_string': INDEX_COL_STRING
    }

    experimental_dict = {
        'tablename': "experiment",
        'hop_view': "experiment_hops",
        'index_name': "experiment_index",
        'index_col_string': INDEX_COL_STRING
    }

    # raise error if poisson_pseudocount can't be converted to float, or
    # if the value is less than 0
    if(float(args.poisson_pseudocount) < 0):
        raise ValueError('poisson_pseudocount must be 0 or greater')

    poisson_pseudocount = float(args.poisson_pseudocount)

    if not exists(args.chr_map):
        raise FileExistsError('File Not Found: %s' %args.chr_map)
    else:
        chr_map_df = pd.read_csv(args.chr_map)
        if not args.standard_chr_format in chr_map_df.columns:
            raise AttributeError('Standard chr format not in chr_map columns')

    # create promoter, background and experiment tables
    if not exists(args.promoter_bed_path):
        raise FileExistsError("File Not Found: %s" %args.promoter_bed_path)
    else:
        promoter_df = pd.read_csv(args.promoter_bed_path,
                                  sep = "\t",
                                  names = PROMOTER_BED_COLNAMES)
        promoter_df = standardize_chr_names(promoter_df,
                                            chr_map_df,
                                            args.standard_chr_format)
        create_db_table(con, promoter_df, promoter_dict)
        del promoter_df

    if not exists(args.background_qbed_path):
        raise FileExistsError("File Not Found: %s" %args.background_qbed_path)
    else:
        background_df = pd.read_csv(args.background_qbed_path,
                                     sep = "\t",
                                     names = qbed_COLNAMES)
        background_df = standardize_chr_names(background_df,
                                            chr_map_df,
                                            args.standard_chr_format)
        create_db_table(con, background_df, background_dict)
        del background_df

    if not exists(args.experimental_qbed_path):
        raise FileExistsError("File Not Found: %s" %args.experimental_qbed_path)
    else:
        experimental_df = pd.read_csv(args.experimental_qbed_path,
                                  sep = "\t",
                                  names = qbed_COLNAMES)
        experimental_df = standardize_chr_names(experimental_df,
                                                  chr_map_df,
                                                  args.standard_chr_format)
        create_db_table(con, experimental_df, experimental_dict)
        del experimental_df

    # create background hop view
    create_hop_view(con,
                    background_dict['tablename'],
                    background_dict['hop_view'])
    # create experimental hop view
    create_hop_view(con,
                    experimental_dict['tablename'],
                    experimental_dict['hop_view'])

    bg_df   = pd.read_sql_query("SELECT * FROM %s"
                                %background_dict['hop_view'], con)
    expr_df = pd.read_sql_query("SELECT * FROM %s"
                                %experimental_dict['hop_view'], con)

    # get some hop numbers
    total_background_hops = len(bg_df)
    total_expr_hops       = len(expr_df)
    hop_ratio             = float(total_expr_hops) / float(total_background_hops)
    bg_plus_expr_hops     = total_background_hops + total_expr_hops

    # join background and experimental hop tables
    quant_df = pd.merge(bg_df, expr_df,
                        on = ['chrom', 'chromStart', 'chromEnd'])\
                 .rename(columns={'reads_x': 'bg_hops',
                                  'reads_y': 'expr_hops'})\
                 .dropna(axis = 'rows')

    # Calculate Statistics
    #usage: scistat.poisson.cdf(x,mu)
    # where X are the experimental hops,
    # and mu is the background hops * (total_expr_hops/total_background_hops)
    quant_df['poisson_pval'] = \
        [1-pois.cdf(quant_df.loc[index,'expr_hops'] + poisson_pseudocount,
                    (quant_df.loc[index, 'bg_hops'] * hop_ratio)+\
                        poisson_pseudocount)
        for index in range(len(quant_df))]

    # usage:
    # scistat.hypergeom.cdf(x,M,n,N)
    # where x is observed number of type I events
    # (white balls in draw) (experiment hops at locus)
    # M is total number of balls (total number of hops)
    # n is total number of white balls (total number of expeirment hops)
    # N is the number of balls drawn (total hops at a locus)
    quant_df['hypergeom_pval'] = \
        [1-hypergeom.cdf(quant_df.loc[index,'expr_hops']-1,
                         bg_plus_expr_hops,
                         total_expr_hops,
                         quant_df.loc[index,'expr_hops'] + \
                            quant_df.loc[index, 'bg_hops'])
        for index in range(len(quant_df))]

    # close db connection
    con.close()

    # write to CWD
    quant_output_path = splitext(basename(args.experimental_qbed_path))[0] + \
                            "_stats.csv"
    quant_df.to_csv(quant_output_path, index=False)

if __name__ == "__main__":
    sys.exit(main())

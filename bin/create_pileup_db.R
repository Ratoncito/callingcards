#!/usr/bin/env Rscript

library(callingCardsTools)
library(optparse)

main <- function(args) {

create_cc_database(args$pileup_path,
                   args$barcode_length,
                   args$background_data_path,
                   args$db_output_path)

} # end main()

parseArguments <- function() {
  # parse and return cmd line input

  option_list <- list(
    make_option(c("-p", "--pileup_path"),
                help = "path to the samtools mpileup output"),
    make_option(c("-b", "--barcode_length"),
                help = "length of the barcode which maps reads to TFs"),
    make_option(c("-o", "--db_output_path"),
                help = "path/name of the output sqlite database file"),
    make_option(c("-d", "--background_data_path"),
                help = paste0("path to the background data which will ",
                "be used to calculate enrichment")))

  args <- parse_args(
            OptionParser(option_list = option_list))

  return(args)
} # end parseAarguments

main(parseArguments()) # call main method

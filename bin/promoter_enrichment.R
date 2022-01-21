#!/usr/bin/env Rscript

library(callingCardsTools)
library(optparse)

main <- function(args) {

enrichment_df = calculateEnrichmentByHops(args$pileup_db,
                                          args$promoter_bedfile_path,
                                          args$barcode_tf_map_path,
                                          args$stranded)

write.csv(enrichment_df, paste0(args$out_name, "_enrichment.csv"))

} # end main()

parseArguments <- function() {
  # parse and return cmd line input

  option_list <- list(
    make_option(c("-p", "--pileup_db"),
                help = "path to the pileup sqlite database file"),
    make_option(c("-b", "--promoter_bedfile_path"),
                help = "a bed file of regions of interest"),
    make_option(c("-m", "--barcode_tf_map_path"),
                help = "path to a tsv which maps TFs to barcodes"),
    make_option(c("-s", "--stranded"),
                help = paste0("Boolean with values TRUE/FALSE. TRUE means ",
                "that only reads on the same strand of the feature are ",
                "counted. FALSE counts reads on either strand")),
    make_option(c("-o", "--out_name"),
                help = paste0("name of the output enrichment csv. ",
                "Enter only the basename, no extension")))

  args <- parse_args(
            OptionParser(option_list = option_list))

  return(args)
} # end parseArguments

main(parseArguments()) # call main method

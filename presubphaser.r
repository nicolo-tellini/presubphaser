#!/usr/bin/env Rscript
# Title: Homeologue Pairing via Shared BUSCO Genes
# Author: Nicolò T.
# Status: Draft
# Usage: Rscript presubphaser.r --input full_table.tsv --output table_hom_pairing [--chrom_regex SUPER]

# Options ----
options(warn = 1)
options(stringsAsFactors = FALSE)
options(scipen = 999)

# Argument parsing ----
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--input"),
              type    = "character",
              default = NULL,
              help    = "Path to BUSCO full_table.tsv file [required]",
              metavar = "FILE"),
  make_option(c("-o", "--output"),
              type    = "character",
              default = NULL,
              help    = "Output prefix (no extension) for matrix TSV and heatmap PDF [required]",
              metavar = "PREFIX"),
  make_option(c("-r", "--chrom_regex"),
              type    = "character",
              default = "SUPER",
              help    = "Regex pattern to filter chromosome/sequence names [default: %default]",
              metavar = "REGEX")
)

parser <- OptionParser(
  usage       = "%prog --input FILE --output PREFIX [--chrom_regex REGEX]",
  option_list = option_list,
  description = "Compute the matrix of single-copy BUSCO genes shared between chromosomes, as input for homeologue pairing."
)

opt <- parse_args(parser)

# Validate required arguments
if (is.null(opt$input)) {
  print_help(parser)
  stop("--input is required.", call. = FALSE)
}
if (is.null(opt$output)) {
  print_help(parser)
  stop("--output is required.", call. = FALSE)
}
if (!file.exists(opt$input)) {
  stop(paste("Input file not found:", opt$input), call. = FALSE)
}

table_file   <- opt$input
out_prefix   <- opt$output
chrom_regex  <- opt$chrom_regex


# Libraries ----
suppressPackageStartupMessages({
  library(gtools)
  library(pheatmap)
})

# Body ----
d <- read.table(table_file,
                sep              = "\t",
                header           = FALSE,
                comment.char     = "#",
                quote            = "",
                fill             = TRUE,
                stringsAsFactors = FALSE)

colnames(d)[1:3] <- c("busco", "status", "seq")

# Keep Complete and Duplicated hits, retain only busco + seq columns
d <- d[d$status %in% c("Complete", "Duplicated"), c("busco", "seq")]

# Filter to chromosomes matching the regex
d <- d[grepl(chrom_regex, d$seq), ]

# One entry per BUSCO per chromosome
d <- unique(d)

if (nrow(d) == 0) {
  stop(paste("No rows remaining after filtering with regex:", chrom_regex), call. = FALSE)
}

chroms <- mixedsort(unique(d$seq))

# Build presence/absence matrix and compute pairwise shared-BUSCO counts
tab    <- table(d$busco, d$seq)
pa     <- matrix(as.integer(tab > 0), nrow = nrow(tab), dimnames = dimnames(tab))
shared <- crossprod(pa)
shared <- shared[chroms, chroms, drop = FALSE]
#diag(shared) <- 0

# Write shared matrix
nome_matrice <- paste0(out_prefix, ".shared_matrix.tsv")
write.table(shared, file = nome_matrice, sep = "\t", quote = FALSE, col.names = NA)

# Build lower-triangle display matrix
shared_lower <- shared
shared_lower[upper.tri(shared_lower)] <- NA
diag(shared_lower) <- NA

col <- colorRampPalette(c("#F7F7F7", "gold2", "orange"))(100)
n   <- nrow(shared_lower)
num <- matrix("", n, n)
num[lower.tri(shared_lower, diag = FALSE)] <- shared_lower[lower.tri(shared_lower, diag = FALSE)]

# Write heatmap PDF
nome_heatmap <- paste0(out_prefix, ".heatmap.pdf")
pdf(nome_heatmap, width = 8, height = 8)
pheatmap(shared_lower,
         cluster_rows    = FALSE,
         cluster_cols    = FALSE,
         color           = col,
         na_col          = "white",
         display_numbers = num,
         border_color    = NA,
         number_color    = "black")
invisible(dev.off())

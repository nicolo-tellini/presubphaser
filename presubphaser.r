#!/usr/bin/env Rscript
# Title: Homeologue Pairing via Shared BUSCO Genes
# Author: Nicolò T.
# Status: Draft
# Usage: Rscript busco_heatmap.R --input full_table.tsv --output table_hom_pairing [--chrom_regex SUPER]

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
              metavar = "REGEX"),
  make_option(c("-k", "--clusters"),
              type    = "integer",
              default = NULL,
              help    = "Number of clusters (optional; if omitted only the heatmap and matrix are produced)",
              metavar = "INT")
)

parser <- OptionParser(
  usage       = "%prog --input FILE --output PREFIX [--chrom_regex REGEX] [--clusters K]",
  option_list = option_list,
  description = "Compute homeologue pairing across chromosomes via shared single-copy BUSCO genes."
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
n_clusters   <- opt$clusters


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

# Clustering (only if --clusters is provided) ----
if (!is.null(n_clusters)) {

  max_shared          <- max(shared[shared > 0])
  dist_mat            <- as.dist(1 - shared / max_shared)
  hc                  <- hclust(dist_mat, method = "average")  # UPGMA
  cluster_assignments <- cutree(hc, k = n_clusters)

  # Write cluster assignments
  clust_df   <- data.frame(
    chromosome = names(cluster_assignments),
    cluster    = cluster_assignments,
    stringsAsFactors = FALSE
  )
  clust_df   <- clust_df[order(clust_df$cluster, clust_df$chromosome), ]
  clust_file <- paste0(out_prefix, ".clusters.tsv")
  write.table(clust_df, file = clust_file, sep = "\t", quote = FALSE, row.names = FALSE)

  # Heatmap with cluster annotation
  anno_row <- data.frame(
    Cluster = factor(paste0("C", cluster_assignments)),
    row.names = names(cluster_assignments)
  )

  nome_heatmap_clust <- paste0(out_prefix, ".heatmap_clustered.pdf")
  pdf(nome_heatmap_clust, width = 16, height = 16)

  ord              <- order(cluster_assignments, names(cluster_assignments))
  shared_ord       <- shared[ord, ord]
  shared_lower_ord <- shared_ord
  shared_lower_ord[upper.tri(shared_lower_ord)] <- NA
  diag(shared_lower_ord) <- NA
  n2   <- nrow(shared_lower_ord)
  num2 <- matrix("", n2, n2)
  num2[lower.tri(shared_lower_ord, diag = FALSE)] <-
    shared_lower_ord[lower.tri(shared_lower_ord, diag = FALSE)]

  pheatmap(shared_lower_ord,
           cluster_rows    = FALSE,
           cluster_cols    = FALSE,
           color           = col,
           na_col          = "white",
           display_numbers = num2,
           border_color    = NA,
           number_color    = "black",
           annotation_row  = anno_row[rownames(shared_lower_ord), , drop = FALSE],
           annotation_col  = anno_row[colnames(shared_lower_ord), , drop = FALSE],
           main            = paste0("Chromosomes ordered by cluster (k=", n_clusters, ")"))
  invisible(dev.off())

  # Homeolog group summary
  idx      <- which(lower.tri(shared), arr.ind = TRUE)
  pairs_df <- data.frame(
    chr1   = rownames(shared)[idx[, 1]],
    chr2   = colnames(shared)[idx[, 2]],
    shared = shared[idx],
    stringsAsFactors = FALSE
  )

  lines <- sprintf("--- Homeolog groups (k = %d clusters) ---", n_clusters)
  for (g in sort(unique(cluster_assignments))) {
    members <- names(cluster_assignments[cluster_assignments == g])
    if (length(members) > 1) {
      sub_vals    <- pairs_df$shared[pairs_df$chr1 %in% members & pairs_df$chr2 %in% members]
      within_mean <- round(mean(sub_vals), 1)
    } else {
      within_mean <- NA
    }
    lines <- c(lines,
      sprintf("  Cluster %d [n=%d, within mean=%s]: %s",
              g, length(members),
              ifelse(is.na(within_mean), "NA", sprintf("%.1f", within_mean)),
              paste(members, collapse = ", ")))
  }
  lines <- c(lines, "===============================")

  summary_text <- paste(lines, collapse = "\n")
  cat("\n", summary_text, "\n\n")

  summary_file <- paste0(out_prefix, ".summary.txt")
  writeLines(summary_text, summary_file)

}

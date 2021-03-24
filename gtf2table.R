#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (!require("argparse")) install.packages("argparse")
  library(argparse)
  if (!require("data.table")) install.packages("data.table")
  library("data.table")
})

# create parser object
docstring <- "Converts a GTF file to a tab separated table with
ENSEMBL ids, gene names and gene biotype."
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action = "store_true", default = TRUE,
                    help = "Print extra output [default]")
parser$add_argument("-o", "--output", default = NULL,
                    help = "Output file name")
parser$add_argument("file", nargs = 1, help = "GTF file")

args <- parser$parse_args()
file <- args$file

# Check 
skip <- sum(sapply(readLines(con = file, n = 1000), function(x) {
  substr(x, 1, 2) == "##"
}))

# Read annotation file
if (args$verbose) cat(paste0("Reading GTF file...\n"))
annotation <- fread(file, sep = "\t", skip = skip, stringsAsFactors = F, showProgress = T)
if (args$verbose) cat(paste0("Read GTF file with ", nrow(annotation), " records...\n"))
annotation_gene <- annotation[annotation$V3 == "gene", ]
if (args$verbose) cat(paste0("Found ", nrow(annotation_gene), " genes...\n"))
if (args$verbose) cat(paste0("Collecting ENSEMBL IDs, HGNC symbols and gene biotype...\n"))
gene_info <- data.frame(do.call(rbind, lapply(annotation_gene$V9, function(x) {
  split_string <- unlist(strsplit(x, split = "; "))
  split_string <- gsub(pattern = "\\\"", replacement = "", split_string)
  bind_res <- do.call(rbind, strsplit(split_string, split = " "))
  conv_string <- bind_res[, 2]
  names(conv_string) <- bind_res[, 1]
  conv_string[c("gene_id", "gene_type", "gene_name")]
})), stringsAsFactors = F)

# Check for dupliacted genes
dup_genes <- sum(duplicated(gene_info$gene_name))
if (args$verbose) cat(paste0("Found ", dup_genes, " duplicated genes...\n"))

if (dup_genes > 0) {
  if (args$verbose) cat(paste0("Making genes unique...\n"))
  gene_info$gene_name <- make.unique(gene_info$gene_name)
}

# Export table
if (!is.null(args$output)) {
  outfile <- args$output
} else {
  outfile <- "genes.tsv"
}
if (args$verbose) cat(paste0("Writing table...\n"))
fwrite(gene_info, file = outfile, sep = "\t", col.names = T, row.names = F)

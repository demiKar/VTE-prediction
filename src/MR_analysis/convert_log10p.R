#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: convert_log10p.R <input_file> <output_file>")
}

input_file <- args[1]
output_file <- args[2]

# Read input file
df <- read.table(input_file, header = TRUE, sep = "", stringsAsFactors = FALSE)

# Check for LOG10P column
if (!"LOG10P" %in% colnames(df)) {
  stop("Input file must contain a 'LOG10P' column")
}

# Convert LOG10P to P
df$P <- 10^(-df$LOG10P)

# Write output
write.table(df, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE)


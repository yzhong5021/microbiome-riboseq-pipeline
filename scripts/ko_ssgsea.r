#!/usr/bin/env Rscript
library(GSVA)
library(data.table)
library(optparse)

option_list <- list(
  make_option(c("-t", "--tpm-file"),
              type="character",
              help="Path to the input TPM data file (2-column TSV: gene, tpm)",
              metavar="character"),
  make_option(c("-k", "--ko-file"),
              type="character",
              help="Path to the KEGG KO annotation file (2-column TSV: QUERY_NAME, name)",
              metavar="character"),
  make_option(c("-o", "--output-dir"),
              type="character",
              help="Path for the output directory for ssGSEA scores",
              metavar="character"),
  make_option(c("-f", "--output-file"),
	      type = "character",
	      help = "prefix of output file",
	      metavar="character"),
  make_option(c("-s", "--species-name"),
              type="character",
              help="Name of the species, used for output file naming and column headers",
              metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$`tpm-file`) || is.null(opt$`ko-file`) || is.null(opt$`output-dir`) || is.null(opt$`species-name`)){
  print_help(opt_parser)
  stop("All arguments (-t, -k, -o, -s) must be supplied.", call.=FALSE)
}

# loading raw tpm data

print(paste("Loading TPM data from:", opt$`tpm-file`))
tpm_data <- read.delim(opt$`tpm-file`, header = TRUE)

expression_matrix <- as.matrix(tpm_data[, "tpm", drop = FALSE])
rownames(expression_matrix) <- tpm_data$gene
colnames(expression_matrix) <- opt$`species-name`

print("Expression data prepared successfully:")
print(head(expression_matrix))

# loading and preparing gene sets

print(paste("Loading KEGG KO annotations from:", opt$`ko-file`))
ko_dt <- fread(
  opt$`ko-file`,
  col.names = c("QUERY_NAME", "name"),
  nThread = 8
)

print("File reading complete. Reshaping data into gene sets...")

ko_dt <- ko_dt[name != "" & !is.na(name)]

tidy_ko_data <- ko_dt[, .(name = unlist(strsplit(trimws(name), "[;,|]"))), by = QUERY_NAME]

ko_gene_sets <- split(tidy_ko_data$QUERY_NAME, tidy_ko_data$name)

ko_gene_sets <- lapply(ko_gene_sets, unique)

print(paste("Prepared", length(ko_gene_sets), "gene sets."))

print("Running ssGSEA analysis with GSVA...")

# filter out zero-count genes
is_expressed <- expression_matrix[, opt$`species-name`] > 0
expression_matrix_filtered <- expression_matrix[is_expressed, , drop = FALSE]

print(paste("Original number of genes:", nrow(expression_matrix)))
print(paste("Number of genes after filtering zero-tpm genes", nrow(expression_matrix_filtered)))

# run ssGSEA
ssgsea_scores <- gsva(
  expr = expression_matrix_filtered,
  gset.idx.list = ko_gene_sets,
  method = "ssgsea",
  min.sz = 4,
  parallel.sz = 8,
  verbose = TRUE,
  ssgsea.norm = TRUE
)

print("ssGSEA analysis complete.")

# saving results
ssgsea_df <- as.data.frame(ssgsea_scores)
ssgsea_df <- cbind(pathway = rownames(ssgsea_df), ssgsea_df)
rownames(ssgsea_df) <- NULL
colnames(ssgsea_df) <- c("pathway", "enrichment")

ssgsea_df <- ssgsea_df[order(ssgsea_df$enrichment, decreasing = TRUE), ]

output_file_path <- file.path(opt$`output-dir`, paste0(opt$`species-name`, "_", opt$`output-file`, "_gsea.tsv"))

write.table(ssgsea_df, file = output_file_path, sep = "\t", row.names = FALSE, quote = FALSE)

print(paste("Results saved to:", output_file_path))
print("Top 10 most active pathways/functions:")
print(head(ssgsea_df, 10))

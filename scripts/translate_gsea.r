#!/usr/bin/env Rscript

library(optparse)
library(GSVA)
library(data.table)

option_list <- list(
  make_option(c("-a", "--abundance-file"), type = "character", default = NULL,
              help = "Path to the input abundance file (e.g., VMS_B7_1_25_clean_genefamilies.tsv)", metavar = "character"),
  make_option(c("-g", "--go-map-file"), type = "character", default = NULL,
              help = "Path to the GO term to UniRef mapping file (e.g., map_go_uniref90.txt)", metavar = "character"),
  make_option(c("-n", "--go-name-file"), type = "character", default = NULL,
              help = "Path to the GO ID to GO name mapping file (e.g., map_go_name.tsv)", metavar = "character"),
  make_option(c("-o", "--output-file"), type = "character", default = "ssgsea_pathway_scores.tsv",
              help = "Path for the output SSGSEA pathway scores file [default: %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$`abundance-file`) || is.null(opt$`go-map-file`) || is.null(opt$`go-name-file`)) {
  print_help(opt_parser)
  stop("All input files (--abundance-file, --go-map-file, --go-name-file) must be specified.", call. = FALSE)
}

abundance_file <- opt$`abundance-file`
go_map_file <- opt$`go-map-file`
go_name_file <- opt$`go-name-file`
output_file <- opt$`output-file`


message("Reading abundance file...")
abundance_dt <- fread(abundance_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
setnames(abundance_dt, c("GeneFamily", "Abundance"))
message("Abundance file read. Initial rows:")
print(head(abundance_dt))

message("Filtering abundance data...")
abundance_dt <- abundance_dt[!grepl("\\|", GeneFamily)][GeneFamily != "UNMAPPED"]
message("Abundance data filtered. Filtered rows:")
print(head(abundance_dt))

expression_matrix <- as.matrix(abundance_dt[, .(Abundance)])
rownames(expression_matrix) <- abundance_dt$GeneFamily
colnames(expression_matrix) <- "MySample"

message("Reading and processing GO map file...")
go_map_lines <- readLines(go_map_file)
go_map_list <- strsplit(go_map_lines, "\\s+")
go_map_dt <- data.table(
  GO_Term = sapply(go_map_list, `[`, 1),
  GeneID = sapply(go_map_list, `[`, -1)
)
tidy_go_data <- go_map_dt[, .(GeneID = unlist(GeneID)), by = GO_Term]

message("GO map processed. Tidy GO data sample:")
print(head(tidy_go_data))

go_gene_sets <- split(tidy_go_data$GeneID, tidy_go_data$GO_Term)
go_gene_sets <- lapply(go_gene_sets, unique)

message(paste("Prepared", length(go_gene_sets), "GO term gene sets."))
message("Filtering expression matrix for expressed genes...")

is_expressed <- expression_matrix[, "MySample", drop = FALSE] > 0
expression_matrix_filtered <- expression_matrix[is_expressed, , drop = FALSE]

message(paste("Number of genes for GSEA (Abundance > 0):", nrow(expression_matrix_filtered)))




ssgsea_scores <- gsva(
  expr = expression_matrix_filtered,
  gset.idx.list = go_gene_sets,
  method = "ssgsea",
  min.sz = 5,
  max.sz = 500,
  parallel.sz = 32,
  verbose = TRUE,
  ssgsea.norm = TRUE
)

ssgsea_df <- as.data.frame(ssgsea_scores)
ssgsea_df <- cbind(GO_ID = rownames(ssgsea_df), ssgsea_df)
rownames(ssgsea_df) <- NULL

ssgsea_df <- as.data.table(ssgsea_scores, keep.rownames = "GO_ID")

go_name_map <- fread(go_name_file, header = FALSE, sep = "\t", col.names = c("GO_ID", "GO_Name"))
go_name_map[, GO_ID := trimws(GO_ID)]
ssgsea_df[, GO_ID := trimws(GO_ID)] # Also good to trim the other key, just in case

message("Merging scores with GO names...")
results_translated <- merge(ssgsea_df, go_name_map, by = "GO_ID", all.x = TRUE)

results_translated <- results_translated[, c("GO_ID", "GO_Name", "MySample")]

results_sorted <- results_translated[order(results_translated$MySample, decreasing = TRUE), ]

write.table(results_sorted, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

print(paste("Results saved to:", output_file))

#!/usr/bin/env Rscript

# Per-ORF DESeq2 interaction-model DTE (issue #168, Tier 2).
#
# Fits ~ condition + seq_type + condition:seq_type per ORF. The
# interaction coefficient is the per-ORF differential translation
# efficiency. Implementation mirrors the deltaTE module's "DESeq2
# interaction model" block but with:
#
#   - rows = ORFs (not genes)
#   - Ribo-seq numerator = per-ORF P-site count matrix
#     (`orf_psite_counts.tsv`, from issue #166).
#   - RNA-seq denominator = gene-level Salmon counts, joined to each ORF
#     via `orf_to_gene.tsv` (issue #167). ORFs with no host gene
#     (novel_u with no containing transcript) are dropped: they cannot
#     receive a denominator.
#
# Row-independence caveat: ORFs from the same gene share the gene-level
# RNA-seq denominator. After the join, those denominator rows are
# perfectly correlated. DESeq2 treats them as independent observations;
# this is a known limitation documented in docs/usage.md.

suppressPackageStartupMessages({
    library(DESeq2)
    library(data.table)
    library(dplyr)
    library(tibble)
    library(readr)
    library(purrr)
    library(ggplot2)
    library(BiocParallel)
})

is_valid_string <- function(x) !is.null(x) && nzchar(trimws(x))

parse_args <- function(x) {
    if (!nzchar(x) || x == "null") return(list())
    strsplit(x, " ?--")[[1]][-1] |>
        strsplit("\\\\s+", perl = TRUE) |>
        (\\(parts) set_names(map_chr(parts, 2), map_chr(parts, 1)))() |>
        discard(is.na) |>
        as.list()
}

################################################
## Parse Parameters                           ##
################################################

opt <- list(
    output_prefix       = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    ribo_count_file     = '$ribo_counts',
    rna_count_file      = '$rna_counts',
    orf_to_gene_file    = '$orf_to_gene',
    sample_file         = '$samplesheet',
    contrast_variable   = '$contrast_variable',
    reference_level     = '$reference',
    target_level        = '$target',
    sample_id_col       = "sample",
    seq_type_col        = "type",
    batch_col           = NULL,
    orf_id_col          = "orf_id",
    gene_id_col         = "gene_id",
    min_count           = as.integer(10),
    min_samples         = as.integer(3),
    alpha               = 0.05,
    fit_type            = "parametric",
    cores               = as.integer('$task.cpus')
)
opt_types <- map(opt, class)

args_opt <- parse_args('$task.ext.args')
for (ao in names(args_opt)) {
    if (!ao %in% names(opt)) stop(paste("Invalid option:", ao))
    if (!is.null(opt[[ao]])) args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
    opt[[ao]] <- args_opt[[ao]]
}

required_opts <- c("contrast_variable", "reference_level", "target_level", "output_prefix")
missing <- required_opts[!map_lgl(opt[required_opts], is_valid_string)]
if (length(missing) > 0) stop(paste("Missing required options:", paste(missing, collapse = ", ")))

walk(c("ribo_count_file", "rna_count_file", "orf_to_gene_file", "sample_file"), ~ {
    if (!is_valid_string(opt[[.x]]) || !file.exists(opt[[.x]])) {
        stop(paste0("Invalid or missing file: ", .x))
    }
})

################################################
## Setup Parallelization                      ##
################################################

if (opt\$cores > 1) register(MulticoreParam(opt\$cores))

################################################
## Read Inputs                                ##
################################################

ribo <- fread(opt\$ribo_count_file, data.table = FALSE)
rownames(ribo) <- ribo[[opt\$orf_id_col]]
ribo[[opt\$orf_id_col]] <- NULL
ribo <- mutate(ribo, across(everything(), as.integer))

rna_raw <- fread(opt\$rna_count_file, data.table = FALSE)
gene_id_col_idx <- which(colnames(rna_raw) == opt\$gene_id_col)
if (length(gene_id_col_idx) == 0) gene_id_col_idx <- 1
rownames(rna_raw) <- rna_raw[[gene_id_col_idx]]
# Drop gene_id, optional gene_name column (col 2)
non_sample_cols <- intersect(colnames(rna_raw), c(opt\$gene_id_col, "gene_name"))
rna <- rna_raw[, !colnames(rna_raw) %in% non_sample_cols, drop = FALSE]
rna <- mutate(rna, across(everything(), \\(x) as.integer(round(as.numeric(x)))))

orf_to_gene <- fread(opt\$orf_to_gene_file, data.table = FALSE)
required_o2g <- c(opt\$orf_id_col, opt\$gene_id_col)
if (!all(required_o2g %in% colnames(orf_to_gene))) {
    stop(paste("orf_to_gene.tsv must have columns:", paste(required_o2g, collapse = ", ")))
}
# Collapse to first gene per ORF (deterministic via file order).
orf_to_gene <- orf_to_gene |>
    distinct(across(all_of(opt\$orf_id_col)), .keep_all = TRUE) |>
    select(all_of(c(opt\$orf_id_col, opt\$gene_id_col)))

sample_sheet <- fread(opt\$sample_file, data.table = FALSE)

opt\$sample_id_col <- make.names(opt\$sample_id_col)
opt\$seq_type_col <- make.names(opt\$seq_type_col)
opt\$contrast_variable <- make.names(opt\$contrast_variable)
if (!is.null(opt\$batch_col)) opt\$batch_col <- make.names(opt\$batch_col)

required_cols <- c(opt\$sample_id_col, opt\$seq_type_col, opt\$contrast_variable)
missing_cols <- setdiff(required_cols, colnames(sample_sheet))
if (length(missing_cols) > 0) {
    stop(paste("Missing columns in sample sheet:", paste(missing_cols, collapse = ", ")))
}

sample_sheet <- sample_sheet |>
    distinct(across(all_of(opt\$sample_id_col)), .keep_all = TRUE) |>
    column_to_rownames(opt\$sample_id_col)

# Restrict samples to those present in both matrices.
ribo_samples <- intersect(colnames(ribo), rownames(sample_sheet))
rna_samples  <- intersect(colnames(rna),  rownames(sample_sheet))
keep_samples <- union(ribo_samples, rna_samples)
sample_sheet <- sample_sheet[keep_samples, , drop = FALSE]

sample_sheet <- sample_sheet |>
    mutate(
        across(all_of(opt\$contrast_variable), factor),
        across(all_of(opt\$seq_type_col), factor)
    )
if (!is.null(opt\$batch_col) && opt\$batch_col %in% colnames(sample_sheet)) {
    sample_sheet <- mutate(sample_sheet, across(all_of(opt\$batch_col), factor))
}

seq_type_values <- unique(as.character(sample_sheet[[opt\$seq_type_col]]))
ribo_type <- grep("ribo|rp|fp", seq_type_values, ignore.case = TRUE, value = TRUE)[1]
rna_type  <- grep("rna|mrna|total", seq_type_values, ignore.case = TRUE, value = TRUE)[1]
if (is.na(ribo_type) || is.na(rna_type)) {
    stop(paste("Cannot identify Ribo-seq/RNA-seq from seq_type column. Values:", paste(seq_type_values, collapse = ", ")))
}
cat("Seq types - Ribo:", ribo_type, "RNA:", rna_type, "\\n")

sample_sheet[[opt\$seq_type_col]] <- relevel(sample_sheet[[opt\$seq_type_col]], ref = rna_type)
sample_sheet[[opt\$contrast_variable]] <- relevel(sample_sheet[[opt\$contrast_variable]], ref = opt\$reference_level)

ribo_sample_ids <- rownames(sample_sheet)[sample_sheet[[opt\$seq_type_col]] == ribo_type]
rna_sample_ids  <- rownames(sample_sheet)[sample_sheet[[opt\$seq_type_col]] == rna_type]

ribo <- ribo[, ribo_sample_ids, drop = FALSE]
rna  <- rna[,  rna_sample_ids,  drop = FALSE]

################################################
## Build per-ORF combined count matrix        ##
################################################

# Keep only ORFs that map to a gene present in the RNA matrix.
gene_lookup <- setNames(orf_to_gene[[opt\$gene_id_col]], orf_to_gene[[opt\$orf_id_col]])
orf_ids <- rownames(ribo)
gene_ids_for_orfs <- gene_lookup[orf_ids]

keep <- !is.na(gene_ids_for_orfs) & gene_ids_for_orfs %in% rownames(rna)
n_dropped_no_gene <- sum(is.na(gene_ids_for_orfs))
n_dropped_no_rna  <- sum(!is.na(gene_ids_for_orfs) & !gene_ids_for_orfs %in% rownames(rna))
cat("Dropped", n_dropped_no_gene, "ORFs with no host gene (novel intergenic).\\n")
cat("Dropped", n_dropped_no_rna,  "ORFs whose host gene has no RNA-seq row.\\n")

orf_ids <- orf_ids[keep]
gene_ids_for_orfs <- gene_ids_for_orfs[keep]
ribo <- ribo[orf_ids, , drop = FALSE]

# Build the RNA block by replicating each gene's RNA row for every ORF
# that maps to it. This is the row-correlated step flagged in the
# usage docs.
rna_by_orf <- rna[gene_ids_for_orfs, , drop = FALSE]
rownames(rna_by_orf) <- orf_ids

# Min-count filter: keep ORFs with >= min_count P-sites in >= min_samples
# of the Ribo-seq samples. DESeq2 dispersion fits are unreliable on
# sparse rows (many zeros) and this filter is the empirical guard
# recommended in the spec.
pass_filter <- rowSums(ribo >= opt\$min_count) >= opt\$min_samples
n_filtered <- sum(!pass_filter)
cat("Filtered out", n_filtered, "low-count ORFs (min_count=",
    opt\$min_count, ", min_samples=", opt\$min_samples, ").\\n")
orf_ids <- orf_ids[pass_filter]
gene_ids_for_orfs <- gene_ids_for_orfs[pass_filter]
ribo <- ribo[orf_ids, , drop = FALSE]
rna_by_orf <- rna_by_orf[orf_ids, , drop = FALSE]

# Combined ORF x (ribo samples | rna samples) matrix.
colnames(ribo)       <- paste0(colnames(ribo), "__ribo")
colnames(rna_by_orf) <- paste0(colnames(rna_by_orf), "__rna")
combined <- cbind(ribo, rna_by_orf)

combined_coldata <- data.frame(
    row.names = colnames(combined),
    stringsAsFactors = FALSE
)
combined_coldata[[opt\$seq_type_col]] <- factor(
    c(rep(ribo_type, ncol(ribo)), rep(rna_type, ncol(rna_by_orf))),
    levels = c(rna_type, ribo_type)
)
combined_coldata[[opt\$contrast_variable]] <- factor(
    c(
        as.character(sample_sheet[ribo_sample_ids, opt\$contrast_variable]),
        as.character(sample_sheet[rna_sample_ids,  opt\$contrast_variable])
    ),
    levels = levels(sample_sheet[[opt\$contrast_variable]])
)
if (!is.null(opt\$batch_col) && opt\$batch_col %in% colnames(sample_sheet)) {
    combined_coldata[[opt\$batch_col]] <- factor(c(
        as.character(sample_sheet[ribo_sample_ids, opt\$batch_col]),
        as.character(sample_sheet[rna_sample_ids,  opt\$batch_col])
    ))
}

################################################
## DESeq2 interaction model                   ##
################################################

design_terms <- c(
    if (!is.null(opt\$batch_col) && opt\$batch_col %in% colnames(combined_coldata)) opt\$batch_col,
    opt\$contrast_variable,
    opt\$seq_type_col,
    paste0(opt\$contrast_variable, ":", opt\$seq_type_col)
)
design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))
cat("Design:", deparse(design_formula), "\\n")

dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(combined),
    colData   = combined_coldata,
    design    = design_formula
)
dds <- DESeq(dds, fitType = opt\$fit_type, parallel = (opt\$cores > 1))

result_names <- resultsNames(dds)
cat("Coefficients:", paste(result_names, collapse = ", "), "\\n")

interaction_coef <- c(
    grep(paste0(opt\$contrast_variable, ".*", opt\$seq_type_col), result_names, value = TRUE),
    grep(paste0(opt\$seq_type_col, ".*", opt\$contrast_variable), result_names, value = TRUE)
)[1]
if (is.na(interaction_coef)) stop("Could not find interaction coefficient")
cat("Interaction coefficient:", interaction_coef, "\\n")

res <- results(dds, name = interaction_coef, alpha = opt\$alpha)
res_df <- as.data.frame(res) |>
    rownames_to_column(opt\$orf_id_col) |>
    mutate(gene_id = gene_ids_for_orfs[match(get(opt\$orf_id_col), orf_ids)]) |>
    select(all_of(opt\$orf_id_col), gene_id, everything())

write_tsv(res_df, paste0(opt\$output_prefix, ".orf_dte.results.tsv"))

################################################
## Dispersion plot (QC)                       ##
################################################

dispersion_path <- paste0(opt\$output_prefix, ".orf_dte.dispersions.png")
tryCatch({
    png(dispersion_path, width = 720, height = 720, res = 100)
    plotDispEsts(dds)
    dev.off()
}, error = function(e) {
    message("plotDispEsts failed: ", conditionMessage(e))
    if (file.exists(dispersion_path)) file.remove(dispersion_path)
})

################################################
## Session info + versions                    ##
################################################

sink(paste0(opt\$output_prefix, ".orf_dte.R_sessionInfo.log"))
print(sessionInfo())
sink()

writeLines(c(
    paste0('"', '$task.process', '":'),
    paste0("    bioconductor-deseq2: ", as.character(packageVersion("DESeq2"))),
    paste0("    r-data.table: ", as.character(packageVersion("data.table")))
), "versions.yml")

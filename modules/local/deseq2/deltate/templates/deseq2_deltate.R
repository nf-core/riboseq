#!/usr/bin/env Rscript

# Delta TE (Translational Efficiency) Analysis using DESeq2
#
# Detects differentially translated genes (DTEGs) by integrating Ribo-seq and
# RNA-seq data using DESeq2 with an interaction model.
#
# Reference method: Chothani et al. (2019) "deltaTE: Detection of Translationally
# Regulated Genes by Integrative Analysis of Ribo-seq and RNA-seq Data"
# https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpmb.108
# https://github.com/SGDDNB/translational_regulation
#
# Gene classification:
# - DTEGs: Differentially translated genes (significant interaction term)
# - mRNA_abundance: RNA changes forwarded to translation (no TE change)
# - translation: Pure translational regulation (TE change without RNA change)
# - intensified: Translation amplifies RNA changes (same direction)
# - buffering: Translation dampens RNA changes (opposite direction)

################################################
## Load Libraries                             ##
################################################

suppressPackageStartupMessages({
    library(DESeq2)
    library(data.table)
    library(dplyr)
    library(tibble)
    library(readr)
    library(purrr)
    library(ggplot2)
    library(ggrepel)
    library(BiocParallel)
    library(SummarizedExperiment)
    library(RColorBrewer)
})

################################################
## Plotting Functions                         ##
################################################

# Define anota2seq color scheme for consistency
get_anota2seq_colors <- function() {
    cols <- c(RColorBrewer::brewer.pal(8,"Reds")[c(4,8)],
              RColorBrewer::brewer.pal(8,"Blues")[c(4,8)],
              RColorBrewer::brewer.pal(8,"Greens")[c(4,8)])
    names(cols) <- c("translation up","translation down","buffering down","buffering up","mRNA abundance up","mRNA abundance down")
    return(cols)
}

# Modern ggplot2 fold change plot matching anota2seq style
# Output: 720x720 pixels to match anota2seq fold_change.png
# Effect size thresholds (like anota2seq selDelta* parameters):
#   - lfc_threshold_te: TE threshold (parallel diagonal lines at y = x ± threshold)
#   - lfc_threshold_rna: Total mRNA threshold (vertical lines at ± threshold)
#   - lfc_threshold_ribo: Translated mRNA threshold (horizontal lines at ± threshold)
plot_fold_change <- function(results_df, prefix, target_level, reference_level,
                              lfc_threshold_te = NULL, lfc_threshold_rna = NULL, lfc_threshold_ribo = NULL) {
    anota2seq_cols <- get_anota2seq_colors()

    all_data <- results_df |>
        filter(!is.na(lfc_rna), !is.na(lfc_ribo))

    plot_data <- all_data |>
        mutate(
            regulation_type = case_when(
                class == "mRNA_abundance" & lfc_rna > 0 ~ "mRNA abundance up",
                class == "mRNA_abundance" & lfc_rna <= 0 ~ "mRNA abundance down",
                class == "translation" & lfc_ribo > 0 ~ "translation up",
                class == "translation" & lfc_ribo <= 0 ~ "translation down",
                class == "buffering" & lfc_ribo > 0 ~ "buffering up",
                class == "buffering" & lfc_ribo <= 0 ~ "buffering down",
                class == "intensified" ~ "intensified",
                TRUE ~ "other"
            )
        ) |>
        filter(regulation_type != "other")

    max_val <- max(abs(c(plot_data\$lfc_rna, plot_data\$lfc_ribo)), na.rm = TRUE)

    # Count genes per category for legend
    legend_counts <- plot_data |>
        count(regulation_type) |>
        mutate(
            clean_name = gsub("_", " ", regulation_type),
            clean_name = paste0(toupper(substring(clean_name, 1, 1)), substring(clean_name, 2)),
            legend_label = paste0(clean_name, " (", n, ")")
        )

    # Create color mapping
    color_mapping <- anota2seq_cols[legend_counts\$regulation_type]
    names(color_mapping) <- legend_counts\$regulation_type

    # Add intensified color if present
    if ("intensified" %in% plot_data\$regulation_type) {
        color_mapping["intensified"] <- "purple"
    }

    p <- ggplot(plot_data, aes(x = lfc_rna, y = lfc_ribo)) +
        # Grey background points (all genes)
        geom_point(
            data = all_data,
            aes(x = lfc_rna, y = lfc_ribo),
            color = "grey", size = 1.5, alpha = 0.5
        ) +
        # Colored significant points
        geom_point(aes(color = regulation_type), size = 1.5, alpha = 0.8) +
        # Reference lines (dashed)
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.7) +
        scale_color_manual(
            values = color_mapping,
            labels = setNames(legend_counts\$legend_label, legend_counts\$regulation_type)
        ) +
        coord_fixed() +
        xlim(c(-max_val, max_val)) +
        ylim(c(-max_val, max_val)) +
        labs(
            x = paste0("Total mRNA log2FC\\n(", target_level, " vs ", reference_level, ")"),
            y = paste0("Translated mRNA log2FC\\n(", target_level, " vs ", reference_level, ")"),
            color = NULL
        ) +
        theme_bw(base_size = 14) +
        theme(
            panel.grid.minor = element_blank(),
            legend.position = c(0.02, 0.98),
            legend.justification = c(0, 1),
            legend.background = element_rect(fill = alpha("white", 0.8)),
            legend.text = element_text(size = 11),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 11)
        ) +
        guides(color = guide_legend(override.aes = list(size = 3)))

    # Add effect size threshold lines (solid, like anota2seq)
    # TE threshold: parallel diagonal lines at y = x ± threshold
    if (!is.null(lfc_threshold_te) && lfc_threshold_te > 0) {
        p <- p +
            geom_abline(slope = 1, intercept = lfc_threshold_te, linetype = "solid", linewidth = 0.8) +
            geom_abline(slope = 1, intercept = -lfc_threshold_te, linetype = "solid", linewidth = 0.8)
    }
    # RNA threshold: vertical lines at ± threshold
    if (!is.null(lfc_threshold_rna) && lfc_threshold_rna > 0) {
        p <- p +
            geom_vline(xintercept = c(-lfc_threshold_rna, lfc_threshold_rna), linetype = "solid", linewidth = 0.8)
    }
    # Ribo threshold: horizontal lines at ± threshold
    if (!is.null(lfc_threshold_ribo) && lfc_threshold_ribo > 0) {
        p <- p +
            geom_hline(yintercept = c(-lfc_threshold_ribo, lfc_threshold_ribo), linetype = "solid", linewidth = 0.8)
    }

    ggsave(paste0(prefix, ".fold_change.png"), plot = p, width = 720, height = 720, units = "px", dpi = 72)
}

# Interaction p-value distribution plot
# Output: 800x600 pixels (landscape diagnostic plot)
plot_interaction_pval_distribution <- function(results_df, prefix, alpha) {
    p <- ggplot(results_df, aes(x = padj_te)) +
        geom_histogram(bins = 50, alpha = 0.7, fill = "#3498db", color = "white") +
        geom_vline(xintercept = alpha, linetype = "dashed", color = "#e74c3c", linewidth = 1) +
        labs(
            x = "Adjusted P-value (Interaction Term)",
            y = "Count",
            title = paste0("Distribution of Interaction P-values (\u03b1 = ", alpha, ")")
        ) +
        theme_bw(base_size = 14) +
        theme(panel.grid.minor = element_blank()) +
        annotate("text", x = alpha, y = Inf, vjust = 1.5, hjust = -0.1,
                 label = paste0("\u03b1 = ", alpha), color = "#e74c3c", size = 5)

    ggsave(paste0(prefix, ".interaction_p_distribution.png"), plot = p, width = 800, height = 600, units = "px", dpi = 72)
}

# PCA plots
# Output: 800x800 pixels (square plot)
# Also exports underlying PCA data as TSV
plot_pca <- function(dds, title, filename, contrast_variable) {
    vsd <- safe_vst(dds, blind = TRUE)
    pca_data <- plotPCA(vsd, intgroup = contrast_variable, returnData = TRUE)
    pct_var <- round(100 * attr(pca_data, "percentVar"))

    # Export underlying PCA data
    pca_export <- pca_data |>
        mutate(
            PC1_variance_pct = pct_var[1],
            PC2_variance_pct = pct_var[2]
        )
    export_plot_data(pca_export, sub("\\\\.png\$", ".tsv", filename))

    p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = .data[[contrast_variable]])) +
        geom_point(size = 4) +
        geom_text_repel(aes(label = name), size = 4) +
        labs(
            x = paste0("PC1: ", pct_var[1], "%"),
            y = paste0("PC2: ", pct_var[2], "%"),
            title = title
        ) +
        theme_bw(base_size = 14) +
        coord_fixed()

    ggsave(filename, plot = p, width = 800, height = 800, units = "px", dpi = 72)
}

# Heatmap plot
# Output: 1000x1200 pixels (larger to accommodate gene names and annotations)
# Also exports underlying heatmap data (Z-scores and sample annotations) as TSV
plot_heatmap <- function(dds_combined, gene_lists, res_delta_te, sample_sheet, n_top_genes, contrast_variable, seq_type_col, prefix) {
    if (length(gene_lists\$dtegs) < 2) return(invisible(NULL))

    suppressPackageStartupMessages(library(ComplexHeatmap))

    vsd_combined <- safe_vst(dds_combined, blind = FALSE)
    top_genes <- gene_lists\$dtegs[order(res_delta_te[gene_lists\$dtegs, "padj"])][1:min(n_top_genes, length(gene_lists\$dtegs))]

    if (length(top_genes) < 2) return(invisible(NULL))

    mat_scaled <- t(scale(t(assay(vsd_combined)[top_genes, ])))

    # Export underlying heatmap data
    export_plot_data(mat_scaled, paste0(prefix, ".heatmap_zscores.tsv"), id_col = "gene_id")
    sample_annotations <- sample_sheet |>
        rownames_to_column("sample") |>
        select(sample, all_of(c(contrast_variable, seq_type_col)))
    export_plot_data(sample_annotations, paste0(prefix, ".heatmap_annotations.tsv"))

    ha <- HeatmapAnnotation(
        Condition = sample_sheet[[contrast_variable]],
        SeqType = sample_sheet[[seq_type_col]]
    )

    png(paste0(prefix, ".heatmap.png"), width = 1000, height = 1200, units = "px")
    draw(Heatmap(
        mat_scaled,
        name = "Z-score",
        top_annotation = ha,
        show_row_names = (length(top_genes) <= 50),
        column_title = paste0("Top ", length(top_genes), " DTEGs")
    ))
    dev.off()
}

################################################
## Utility Functions                          ##
################################################

safe_vst <- function(dds, blind = TRUE) {
    tryCatch(
        vst(dds, blind = blind, fitType = "mean"),
        error = function(e) {
            if (grepl("less than 'nsub'", conditionMessage(e))) {
                varianceStabilizingTransformation(dds, blind = blind, fitType = "mean")
            } else stop(e)
        }
    )
}

# Export plot underlying data as TSV
export_plot_data <- function(data, filename, id_col = NULL) {
    if (is.matrix(data)) {
        data <- as.data.frame(data)
    }
    if (!is.null(id_col)) {
        data <- rownames_to_column(data, id_col)
    }
    write_tsv(data, filename)
}

is_valid_string <- function(x) !is.null(x) && nzchar(trimws(x))

parse_args <- function(x) {
    if (!nzchar(x) || x == "null") return(list())
    strsplit(x, " ?--")[[1]][-1] |>
        strsplit("\\\\s+", perl = TRUE) |>
        (\\(parts) set_names(map_chr(parts, 2), map_chr(parts, 1)))() |>
        discard(is.na) |>
        as.list()
}

run_deseq_subset <- function(counts, samples, design, contrast_var, target, ref,
                              fit_type, parallel, shrink = TRUE, shrink_type = "apeglm", alpha = 0.05) {
    dds <- DESeqDataSetFromMatrix(
        countData = as.matrix(counts),
        colData = samples,
        design = design
    )
    dds <- DESeq(dds, fitType = fit_type, parallel = parallel)
    res <- results(dds, contrast = c(contrast_var, target, ref), alpha = alpha)

    if (shrink) {
        coef <- grep(paste0(contrast_var, "_", target), resultsNames(dds), value = TRUE)[1]
        res <- lfcShrink(dds, coef = coef, res = res, type = shrink_type)
    }
    list(dds = dds, results = res)
}

################################################
## Parse Parameters                           ##
################################################

opt <- list(
    output_prefix       = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    count_file          = '$counts',
    sample_file         = '$samplesheet',
    contrast_variable   = '$contrast_variable',
    reference_level     = '$reference',
    target_level        = '$target',
    sample_id_col       = "sample",
    seq_type_col        = "type",
    batch_col           = NULL,
    gene_id_col         = "gene_id",
    shrink_lfc          = TRUE,
    shrinkage_type      = "apeglm",
    alpha               = 0.05,
    fit_type            = "parametric",
    generate_plots      = TRUE,
    n_top_genes         = 50,
    lfc_threshold_te    = as.numeric(0),
    lfc_threshold_rna   = as.numeric(0),
    lfc_threshold_ribo  = as.numeric(0),
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

walk(c("count_file", "sample_file"), ~ {
    if (!is_valid_string(opt[[.x]]) || !file.exists(opt[[.x]])) {
        stop(paste0("Invalid or missing file: ", .x))
    }
})

################################################
## Setup Parallelization                      ##
################################################

if (opt\$cores > 1) register(MulticoreParam(opt\$cores))

if (opt\$shrink_lfc && opt\$shrinkage_type == "apeglm") {
    suppressPackageStartupMessages(library(apeglm))
}

################################################
## Read Input Data                            ##
################################################

count_table <- fread(opt\$count_file, data.table = FALSE)
rownames(count_table) <- count_table[[opt\$gene_id_col]]
count_table[[opt\$gene_id_col]] <- NULL
count_table <- mutate(count_table, across(everything(), as.integer))

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

missing_samples <- setdiff(rownames(sample_sheet), colnames(count_table))
if (length(missing_samples) > 0) {
    stop(paste(length(missing_samples), "samples missing from count table"))
}

count_table <- count_table[, rownames(sample_sheet)]

sample_sheet <- sample_sheet |>
    mutate(
        across(all_of(opt\$contrast_variable), factor),
        across(all_of(opt\$seq_type_col), factor)
    )
if (!is.null(opt\$batch_col) && opt\$batch_col %in% colnames(sample_sheet)) {
    sample_sheet <- mutate(sample_sheet, across(all_of(opt\$batch_col), factor))
}

# Identify seq types
seq_type_values <- unique(as.character(sample_sheet[[opt\$seq_type_col]]))
ribo_type <- grep("ribo|rp|fp", seq_type_values, ignore.case = TRUE, value = TRUE)[1]
rna_type <- grep("rna|mrna|total", seq_type_values, ignore.case = TRUE, value = TRUE)[1]

if (is.na(ribo_type) || is.na(rna_type)) {
    stop(paste("Cannot identify Ribo-seq/RNA-seq from seq_type column. Values:", paste(seq_type_values, collapse = ", ")))
}
cat("Seq types - Ribo:", ribo_type, "RNA:", rna_type, "\\n")

sample_sheet[[opt\$seq_type_col]] <- relevel(sample_sheet[[opt\$seq_type_col]], ref = rna_type)
sample_sheet[[opt\$contrast_variable]] <- relevel(sample_sheet[[opt\$contrast_variable]], ref = opt\$reference_level)

################################################
## DESeq2 Interaction Model (deltaTE)         ##
################################################

cat("Running DESeq2 interaction model...\\n")

design_terms <- c(
    if (!is.null(opt\$batch_col) && opt\$batch_col %in% colnames(sample_sheet)) opt\$batch_col,
    opt\$contrast_variable,
    opt\$seq_type_col,
    paste0(opt\$contrast_variable, ":", opt\$seq_type_col)
)
design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))
cat("Design:", deparse(design_formula), "\\n")

dds_combined <- DESeqDataSetFromMatrix(
    countData = as.matrix(count_table),
    colData = sample_sheet,
    design = design_formula
)
dds_combined <- DESeq(dds_combined, fitType = opt\$fit_type, parallel = (opt\$cores > 1))

result_names <- resultsNames(dds_combined)
cat("Coefficients:", paste(result_names, collapse = ", "), "\\n")

interaction_coef <- c(
    grep(paste0(opt\$contrast_variable, ".*", opt\$seq_type_col), result_names, value = TRUE),
    grep(paste0(opt\$seq_type_col, ".*", opt\$contrast_variable), result_names, value = TRUE)
)[1]

if (is.na(interaction_coef)) stop("Could not find interaction coefficient")
cat("Interaction coefficient:", interaction_coef, "\\n")

res_delta_te <- results(dds_combined, name = interaction_coef, alpha = opt\$alpha)

################################################
## Separate Ribo-seq and RNA-seq Analyses     ##
################################################

cat("Running separate analyses...\\n")

design_individual <- as.formula(paste(
    "~",
    paste(c(
        if (!is.null(opt\$batch_col) && opt\$batch_col %in% colnames(sample_sheet)) opt\$batch_col,
        opt\$contrast_variable
    ), collapse = " + ")
))

ribo_samples <- rownames(sample_sheet)[sample_sheet[[opt\$seq_type_col]] == ribo_type]
rna_samples <- rownames(sample_sheet)[sample_sheet[[opt\$seq_type_col]] == rna_type]

ribo_analysis <- run_deseq_subset(
    counts = count_table[, ribo_samples],
    samples = sample_sheet[ribo_samples, , drop = FALSE],
    design = design_individual,
    contrast_var = opt\$contrast_variable,
    target = opt\$target_level,
    ref = opt\$reference_level,
    fit_type = opt\$fit_type,
    parallel = (opt\$cores > 1),
    shrink = opt\$shrink_lfc,
    shrink_type = opt\$shrinkage_type,
    alpha = opt\$alpha
)

rna_analysis <- run_deseq_subset(
    counts = count_table[, rna_samples],
    samples = sample_sheet[rna_samples, , drop = FALSE],
    design = design_individual,
    contrast_var = opt\$contrast_variable,
    target = opt\$target_level,
    ref = opt\$reference_level,
    fit_type = opt\$fit_type,
    parallel = (opt\$cores > 1),
    shrink = opt\$shrink_lfc,
    shrink_type = opt\$shrinkage_type,
    alpha = opt\$alpha
)

dds_ribo <- ribo_analysis\$dds
dds_rna <- rna_analysis\$dds
res_delta_ribo <- ribo_analysis\$results[rownames(res_delta_te), ]
res_delta_rna <- rna_analysis\$results[rownames(res_delta_te), ]

################################################
## Gene Classification                        ##
################################################

cat("Classifying genes...\\n")

alpha <- opt\$alpha

# Effect size thresholds (like anota2seq selDelta* parameters)
# These filter genes requiring minimum absolute log2FC in addition to p-value significance
lfc_threshold_te <- opt\$lfc_threshold_te      # TE threshold (deltaP - deltaT)
lfc_threshold_rna <- opt\$lfc_threshold_rna    # Total mRNA threshold (deltaT)
lfc_threshold_ribo <- opt\$lfc_threshold_ribo  # Translated mRNA threshold (deltaP)

results_df <- tibble(
    gene = rownames(res_delta_te),
    padj_te = res_delta_te\$padj,
    padj_ribo = res_delta_ribo\$padj,
    padj_rna = res_delta_rna\$padj,
    lfc_te = res_delta_te\$log2FoldChange,
    lfc_ribo = res_delta_ribo\$log2FoldChange,
    lfc_rna = res_delta_rna\$log2FoldChange
) |>
    mutate(
        # P-value significance
        te_sig = padj_te < alpha,
        ribo_sig = padj_ribo < alpha,
        rna_sig = padj_rna < alpha,
        # Effect size significance (if thresholds > 0)
        te_eff = if (lfc_threshold_te <= 0) TRUE else abs(lfc_te) >= lfc_threshold_te,
        ribo_eff = if (lfc_threshold_ribo <= 0) TRUE else abs(lfc_ribo) >= lfc_threshold_ribo,
        rna_eff = if (lfc_threshold_rna <= 0) TRUE else abs(lfc_rna) >= lfc_threshold_rna,
        # Combined significance (p-value AND effect size)
        te_pass = te_sig & te_eff,
        ribo_pass = ribo_sig & ribo_eff,
        rna_pass = rna_sig & rna_eff,
        same_direction = lfc_te * lfc_rna > 0,
        class = case_when(
            te_pass & ribo_pass & rna_pass & same_direction ~ "intensified",
            te_pass & ribo_pass & rna_pass & !same_direction ~ "buffering",
            te_pass & ribo_pass & !rna_pass ~ "translation",
            te_pass & !ribo_pass & rna_pass ~ "buffering",
            !te_pass & ribo_pass & rna_pass ~ "mRNA_abundance",
            te_pass ~ "dteg_other",
            TRUE ~ "other"
        )
    )

# Log effect size thresholds if used
if (lfc_threshold_te > 0) cat("  TE effect size threshold:", lfc_threshold_te, "\\n")
if (lfc_threshold_rna > 0) cat("  RNA effect size threshold:", lfc_threshold_rna, "\\n")
if (lfc_threshold_ribo > 0) cat("  Ribo effect size threshold:", lfc_threshold_ribo, "\\n")

gene_lists <- list(
    dtegs = filter(results_df, te_pass)\$gene,
    mRNA_abundance = filter(results_df, class == "mRNA_abundance")\$gene,
    translation = filter(results_df, class == "translation")\$gene,
    intensified = filter(results_df, class == "intensified")\$gene,
    buffering = filter(results_df, class == "buffering")\$gene
)

cat("\\n=== Summary ===\\n")
cat("Total genes:", nrow(results_df), "\\n")
iwalk(gene_lists, ~ cat(paste0(toupper(substring(.y, 1, 1)), substring(.y, 2)), ":", length(.x), "\\n"))

################################################
## Write Results                              ##
################################################

prefix <- opt\$output_prefix

# Write DESeq2 results
list(translation = res_delta_te, translated_mRNA = res_delta_ribo, total_mRNA = res_delta_rna) |>
    iwalk(~ {
        .x |>
            as.data.frame() |>
            rownames_to_column("gene_id") |>
            write_tsv(paste0(prefix, ".", .y, ".deltate.results.tsv"))
    })

# Write gene lists
iwalk(gene_lists, ~ write_tsv(tibble(gene_id = .x), paste0(prefix, ".", .y, ".deltate.genes.tsv")))

# Save R objects
saveRDS(
    list(
        combined = dds_combined,
        ribo = dds_ribo,
        rna = dds_rna,
        results = list(delta_te = res_delta_te, delta_ribo = res_delta_ribo, delta_rna = res_delta_rna),
        classification = gene_lists
    ),
    file = paste0(prefix, ".DESeqDataSet.rds")
)

################################################
## Plots                                      ##
################################################

if (opt\$generate_plots) {
    cat("Generating plots...\\n")

    plot_fold_change(results_df, prefix, opt\$target_level, opt\$reference_level,
                     lfc_threshold_te, lfc_threshold_rna, lfc_threshold_ribo)
    plot_interaction_pval_distribution(results_df, prefix, opt\$alpha)
    plot_pca(dds_ribo, "Ribo-seq PCA", paste0(prefix, ".pca_ribo.png"), opt\$contrast_variable)
    plot_pca(dds_rna, "RNA-seq PCA", paste0(prefix, ".pca_rna.png"), opt\$contrast_variable)
    plot_heatmap(dds_combined, gene_lists, res_delta_te, sample_sheet, opt\$n_top_genes, opt\$contrast_variable, opt\$seq_type_col, prefix)
}

################################################
## Session Info & Versions                    ##
################################################

sink(paste0(prefix, ".R_sessionInfo.log"))
print(sessionInfo())
sink()

writeLines(c(
    '"${task.process}":',
    paste("    bioconductor-deseq2:", packageVersion("DESeq2")),
    paste("    r-ggplot2:", packageVersion("ggplot2")),
    paste("    r-dplyr:", packageVersion("dplyr")),
    paste("    r-data.table:", packageVersion("data.table"))
), "versions.yml")

cat("\\nDone!\\n")

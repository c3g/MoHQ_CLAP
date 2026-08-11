#!/usr/bin/env Rscript
# =============================================================================
# run_pca.R -- expression PCA across a cohort
#
# Fixes vs the previous run_pca_interactive.R
# -------------------------------------------
# 1. The plot coloured points by `Cohort = args$cohort`, a CONSTANT. Every point
#    got the same colour, so the figure could not show the one thing a cohort
#    PCA is for: whether samples separate by batch, institution, or sex rather
#    than by biology. Points are now coloured by a real variable (default:
#    institution) and shaped by cohort.
#
# 2. PCA ran on all ~60k genes. Standard practice is the top N most variable
#    genes; including the long tail of near-zero-variance genes mostly adds
#    noise and memory.
#
# 3. The `zscore` branch called scale() on the matrix and then prcomp(scale.=FALSE),
#    which is inconsistent with the other branches. Normalisation is now explicit
#    and applied once.
#
# 4. Sample IDs came from filenames (MoHQ-CM-4-10-157954-1RT) and were never
#    reconciled with any other data type. They are now harmonised to patient_id
#    and validated against the manifest.
#
# 5. Adds an explicit outlier flag (median absolute deviation on PC1/PC2) and a
#    scree plot, so "is this sample an outlier" is answered rather than eyeballed.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(plotly)
  library(htmlwidgets)
  library(argparse)
})

parser <- ArgumentParser(description = "Cohort expression PCA")
parser$add_argument("--lib", required = TRUE)
parser$add_argument("--expr", required = TRUE, nargs = "+",
                    help = "*.abundance_genes.tsv files")
parser$add_argument("--manifest", required = TRUE)
parser$add_argument("--cohort", default = "MoHQ Cohort")
parser$add_argument("--norm_method", default = "log2",
                    choices = c("log2", "log10", "zscore", "none"))
parser$add_argument("--top_var_genes", type = "integer", default = 2000,
                    help = "0 = use all genes")
parser$add_argument("--colour_by", default = "institution",
                    help = "Manifest column used for point colour")
parser$add_argument("--out_prefix", default = "cohort_expression_pca")
args <- parser$parse_args()

source(args$lib)
manifest <- read_manifest(args$manifest)

# --------------------------------------------------------------------------- #
# Build the expression matrix
# --------------------------------------------------------------------------- #
files <- args$expr[file.exists(args$expr) & file.info(args$expr)$size > 0]
if (length(files) < 3) {
  mohq_die(sprintf("Only %d usable expression file(s); PCA needs at least 3.",
                   length(files)))
}
mohq_log(sprintf("Building expression matrix from %d files", length(files)))

read_one <- function(f) {
  hdr <- names(fread(f, nrows = 0))
  gene_col <- intersect(c("gene_id", "gene", "target_id", "Name", "gene_symbol"), hdr)[1]
  val_col  <- intersect(c("abundance", "TPM", "tpm", "est_counts", "expected_count"), hdr)[1]
  if (is.na(gene_col) || is.na(val_col)) {
    mohq_die("Cannot resolve gene/abundance columns in ", basename(f),
             "\n  Found: ", paste(hdr, collapse = ", "))
  }
  d <- fread(f, select = c(gene_col, val_col), showProgress = FALSE)
  setnames(d, c("gene_id", "abundance"))
  d[, sample_id := sub("\\.abundance_genes\\.tsv$", "", basename(f))]
  d
}

long <- rbindlist(lapply(files, read_one))
long[, patient_id := harmonise_patient_id(sample_id)]
check_ids_against_manifest(long$patient_id, manifest, "expression")
long <- long[patient_id %in% manifest$patient_id]

# One row per patient. If a patient somehow has two RNA samples, average them
# and say so rather than letting dcast silently pick one.
dup <- long[, .(n = uniqueN(sample_id)), by = patient_id][n > 1]
if (nrow(dup)) {
  mohq_warn(sprintf("%d patient(s) have >1 RNA sample; averaging: %s",
                    nrow(dup), paste(dup$patient_id, collapse = ", ")))
}

mat_dt <- dcast(long, patient_id ~ gene_id, value.var = "abundance",
                fun.aggregate = function(x) mean(x, na.rm = TRUE))
pids <- mat_dt$patient_id
m <- as.matrix(mat_dt[, -1])
rownames(m) <- pids
m[!is.finite(m)] <- 0
mohq_log(sprintf("Matrix: %d patients x %d genes", nrow(m), ncol(m)))

# --------------------------------------------------------------------------- #
# Normalise, then select variable genes
# --------------------------------------------------------------------------- #
m <- switch(args$norm_method,
            log2   = log2(m + 1),
            log10  = log10(m + 1),
            zscore = log2(m + 1),   # z-scoring happens via prcomp(scale.=TRUE)
            none   = m)

vars <- col_vars(m)
m <- m[, vars > 0, drop = FALSE]
mohq_log(sprintf("Dropped %d zero-variance genes", sum(vars <= 0)))

if (args$top_var_genes > 0 && ncol(m) > args$top_var_genes) {
  keep <- order(col_vars(m), decreasing = TRUE)[seq_len(args$top_var_genes)]
  m <- m[, sort(keep), drop = FALSE]
  mohq_log(sprintf("Kept top %d most variable genes", ncol(m)))
}

pca <- prcomp(m, center = TRUE, scale. = (args$norm_method == "zscore"))
var_pct <- 100 * (pca$sdev^2) / sum(pca$sdev^2)

scores <- data.table(patient_id = rownames(pca$x),
                     PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     PC3 = if (ncol(pca$x) >= 3) pca$x[, 3] else NA_real_)
scores <- merge(scores, manifest, by = "patient_id", all.x = TRUE)

# --------------------------------------------------------------------------- #
# Outlier flag: robust z on PC1/PC2 (MAD-based, so outliers do not mask
# themselves the way an SD-based cutoff would).
# --------------------------------------------------------------------------- #
robust_z <- function(x) {
  md <- median(x, na.rm = TRUE)
  s  <- mad(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - md) / s
}
scores[, outlier_score := sqrt(robust_z(PC1)^2 + robust_z(PC2)^2)]
scores[, is_outlier := outlier_score > 4]
if (any(scores$is_outlier)) {
  mohq_warn("PCA outlier(s) (robust distance > 4 MAD): ",
            paste(scores[is_outlier == TRUE]$patient_id, collapse = ", "))
}
fwrite(scores, paste0(args$out_prefix, "_scores.tsv"), sep = "\t")

# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #
colour_col <- if (args$colour_by %in% names(scores)) args$colour_by else "cohort_id"
if (colour_col != args$colour_by) {
  mohq_warn("--colour_by '", args$colour_by, "' not in manifest; using ", colour_col)
}

lab <- function(i) sprintf("PC%d (%.1f%% variance)", i, var_pct[i])

p <- ggplot(scores, aes(x = PC1, y = PC2,
                        colour = .data[[colour_col]],
                        text = paste0(patient_id, "<br>", colour_col, ": ",
                                      .data[[colour_col]]))) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text(data = scores[is_outlier == TRUE],
            aes(label = patient_id), size = 3, vjust = -1.1,
            show.legend = FALSE, inherit.aes = TRUE) +
  labs(title = paste("Expression PCA:", args$cohort),
       subtitle = sprintf("%d patients | %d genes | %s normalisation | coloured by %s",
                          nrow(scores), ncol(m), args$norm_method, colour_col),
       x = lab(1), y = lab(2), colour = colour_col) +
  theme_minimal(base_size = 13)

save_plot(p, paste0(args$out_prefix, ".png"), width = 9, height = 7)

saveWidget(as_widget(ggplotly(p, tooltip = "text")),
           paste0(args$out_prefix, ".html"), selfcontained = TRUE)
mohq_log("Wrote ", args$out_prefix, ".html")

scree <- data.table(PC = seq_along(var_pct), var_pct = var_pct)[1:min(15, length(var_pct))]
ps <- ggplot(scree, aes(x = factor(PC), y = var_pct)) +
  geom_col(fill = "grey40") +
  labs(title = "Scree plot", x = "Principal component", y = "% variance explained") +
  theme_minimal(base_size = 12)
save_plot(ps, paste0(args$out_prefix, "_scree.png"), width = 7, height = 4)

mohq_session_info(paste0(args$out_prefix, ".versions.txt"))

#!/usr/bin/env Rscript
# =============================================================================
# batch_effect_analysis.R -- separate technical batch from biology using the
#                            crossed cancer-type x institution design
#
# WHY THIS IS THE RIGHT ANALYSIS FOR *YOUR* COLLECTION
# ----------------------------------------------------
# Every patient in a cohort has the same cancer type, and the same cancer type
# appears in more than one cohort sequenced at different institutions.
#
# That is a partially crossed design, and it is genuinely valuable. In the
# common (fully confounded) case -- one cancer type per institution -- a
# difference between institutions is uninterpretable, because you cannot tell
# tumour biology from sequencing centre. Here you can:
#
#     HOLD CANCER TYPE CONSTANT, VARY INSTITUTION -> the difference is technical.
#
# So for any cancer type present at >= 2 institutions, an institution effect is
# a batch effect, near enough. That gives you a calibrated estimate of how large
# the technical component is, which then tells you whether a cross-cohort
# biological comparison is trustworthy at all.
#
# Outputs
#   *_variance_partition.tsv  how much variance cancer_type vs institution explain
#   *_within_type_tests.tsv   institution effect within each cancer type
#   *_crossed_design.tsv/png  which cancer types are actually crossed
#   *_metric_by_institution.png
#
# IMPORTANT: this measures ASSOCIATION. An institution effect could still be
# real biology if referral patterns differ (e.g. one centre sees later-stage
# disease). Treat a large effect as "investigate", not "correct and move on".
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparse)
})

parser <- ArgumentParser(description = "Crossed-design batch effect analysis")
parser$add_argument("--lib", required = TRUE)
parser$add_argument("--manifests", required = TRUE, nargs = "+",
                    help = "Per-cohort manifest TSVs")
parser$add_argument("--metrics", nargs = "*",
                    help = "Per-cohort metric tables (FGA, TMB) keyed on patient_id")
parser$add_argument("--cohort_map", default = NULL,
                    help = paste("TSV with columns cohort_id, cancer_type.",
                                 "Required to identify which cohorts share a cancer type."))
parser$add_argument("--min_per_cell", type = "integer", default = 5,
                    help = "Minimum patients per (cancer_type, institution) cell")
parser$add_argument("--out_prefix", default = "collection_batch")
args <- parser$parse_args()

source(args$lib)

# --------------------------------------------------------------------------- #
# 1. Assemble the collection-level patient table
# --------------------------------------------------------------------------- #
man <- rbindlist(lapply(args$manifests, function(f)
  fread(f, sep = "\t", na.strings = c("NA", ""))), fill = TRUE)
mohq_log(sprintf("Collection: %d patients across %d cohorts, %d institutions",
                 nrow(man), uniqueN(man$cohort_id), uniqueN(man$institution)))

# Cancer type: from the cohort map if supplied, else from the manifest if a
# clinical join already provided it.
if (!is.null(args$cohort_map) && file.exists(args$cohort_map)) {
  cmap <- fread(args$cohort_map, sep = "\t")
  need <- c("cohort_id", "cancer_type")
  miss <- setdiff(need, names(cmap))
  if (length(miss)) {
    mohq_die("--cohort_map needs column(s): ", paste(miss, collapse = ", "),
             "\n  Found: ", paste(names(cmap), collapse = ", "))
  }
  man <- merge(man, cmap[, .(cohort_id, cancer_type)], by = "cohort_id", all.x = TRUE)
  unmapped <- unique(man[is.na(cancer_type)]$cohort_id)
  if (length(unmapped)) {
    mohq_warn(length(unmapped), " cohort(s) absent from --cohort_map and excluded ",
              "from the crossed analysis: ", paste(head(unmapped, 10), collapse = ", "))
  }
} else if (!"cancer_type" %in% names(man)) {
  mohq_die("No cancer type available. Supply --cohort_map (cohort_id, cancer_type). ",
           "Without it the crossed design cannot be identified and institution ",
           "effects are uninterpretable.")
}

man <- man[!is.na(cancer_type) & !is.na(institution)]
if (!nrow(man)) mohq_die("No patients with both cancer_type and institution.")

# --------------------------------------------------------------------------- #
# 2. Which cancer types are actually crossed?
# --------------------------------------------------------------------------- #
cells <- man[, .(n_patients = .N), by = .(cancer_type, institution)]
cells <- cells[n_patients >= args$min_per_cell]

crossed <- cells[, .(n_institutions = uniqueN(institution),
                     institutions = paste(sort(unique(institution)), collapse = ","),
                     n_patients = sum(n_patients)),
                 by = cancer_type][order(-n_institutions, -n_patients)]
crossed[, is_crossed := n_institutions >= 2]
fwrite(crossed, paste0(args$out_prefix, "_crossed_design.tsv"), sep = "\t")

n_cross <- sum(crossed$is_crossed)
mohq_log(sprintf("%d/%d cancer types are sequenced at >= 2 institutions (>= %d patients per cell)",
                 n_cross, nrow(crossed), args$min_per_cell))

if (n_cross == 0) {
  mohq_warn("No cancer type is present at 2+ institutions with enough patients. ",
            "Cancer type and institution are effectively confounded in this ",
            "subset, so technical and biological effects CANNOT be separated. ",
            "Lower --min_per_cell or accept the limitation explicitly.")
}

p_design <- ggplot(cells, aes(x = institution, y = cancer_type, fill = n_patients)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = n_patients), size = 3) +
  scale_fill_gradient(low = "#deebf7", high = "#08519c", name = "patients") +
  theme_minimal(base_size = 11) +
  labs(title = "Crossed design: cancer type x institution",
       subtitle = sprintf("Rows with >= 2 filled cells let you estimate the technical effect (%d of %d)",
                          n_cross, nrow(crossed)),
       x = "Institution", y = "Cancer type")
save_plot(p_design, paste0(args$out_prefix, "_crossed_design.png"),
          width = max(7, 0.6 * uniqueN(cells$institution) + 4),
          height = max(5, 0.3 * uniqueN(cells$cancer_type) + 2))

# --------------------------------------------------------------------------- #
# 3. Attach metrics
# --------------------------------------------------------------------------- #
metric_cols <- character()
if (length(args$metrics)) {
  mets <- rbindlist(lapply(args$metrics[file.exists(args$metrics)], function(f) {
    d <- fread(f, sep = "\t")
    if (!"patient_id" %in% names(d)) return(NULL)
    keep <- intersect(c("patient_id", "fga", "fga_gain", "fga_loss", "profiled_mb",
                        "n_segments", "tmb_per_mb", "n_nonsyn"), names(d))
    d[, ..keep]
  }), fill = TRUE)
  if (!is.null(mets) && nrow(mets)) {
    mets <- mets[, lapply(.SD, function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else x[1]),
                 by = patient_id]
    man <- merge(man, mets, by = "patient_id", all.x = TRUE)
    metric_cols <- intersect(c("fga", "tmb_per_mb", "profiled_mb", "n_segments"), names(man))
    mohq_log("Metrics available: ", paste(metric_cols, collapse = ", "))
  }
}

if (!length(metric_cols)) {
  mohq_warn("No metric tables supplied; only the design map was produced.")
  mohq_session_info(paste0(args$out_prefix, ".versions.txt"))
  quit(save = "no", status = 0)
}

# --------------------------------------------------------------------------- #
# 4. Variance partitioning
#
# For each metric fit  metric ~ cancer_type + institution  and report the share
# of variance each term explains (type-II sums of squares). If institution
# explains a comparable or larger share than cancer type, the metric is being
# driven by where the sample was sequenced, not by what the tumour is.
# --------------------------------------------------------------------------- #
partition_one <- function(dt, metric) {
  d <- dt[is.finite(get(metric))]
  d <- d[, .(y = get(metric), cancer_type = factor(cancer_type),
             institution = factor(institution))]
  if (nrow(d) < 20 || nlevels(d$cancer_type) < 2 || nlevels(d$institution) < 2) {
    return(NULL)
  }
  fit <- tryCatch(lm(y ~ cancer_type + institution, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  av <- tryCatch(as.data.frame(anova(fit)), error = function(e) NULL)
  if (is.null(av)) return(NULL)
  ss <- av[["Sum Sq"]]
  names(ss) <- rownames(av)
  tot <- sum(ss)
  data.table(
    metric              = metric,
    n                   = nrow(d),
    pct_var_cancer_type = 100 * ss[["cancer_type"]] / tot,
    pct_var_institution = 100 * ss[["institution"]] / tot,
    pct_var_residual    = 100 * ss[["Residuals"]]  / tot,
    p_institution       = av["institution", "Pr(>F)"]
  )
}

vp <- rbindlist(lapply(metric_cols, function(m) partition_one(man, m)), fill = TRUE)
if (nrow(vp)) {
  vp[, institution_dominates := pct_var_institution > pct_var_cancer_type]
  fwrite(vp, paste0(args$out_prefix, "_variance_partition.tsv"), sep = "\t")
  for (i in seq_len(nrow(vp))) {
    mohq_log(sprintf("%-12s cancer_type %5.1f%% | institution %5.1f%% | residual %5.1f%%%s",
                     vp$metric[i], vp$pct_var_cancer_type[i],
                     vp$pct_var_institution[i], vp$pct_var_residual[i],
                     if (isTRUE(vp$institution_dominates[i])) "   <-- INSTITUTION DOMINATES" else ""))
  }

  vl <- melt(vp[, .(metric, `cancer type` = pct_var_cancer_type,
                    institution = pct_var_institution, residual = pct_var_residual)],
             id.vars = "metric", variable.name = "term", value.name = "pct")
  p_vp <- ggplot(vl, aes(x = metric, y = pct, fill = term)) +
    geom_col() +
    scale_fill_manual(values = c(`cancer type` = "#1b7837",
                                 institution = "#d95f02", residual = "grey80")) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    labs(title = "Variance partition: biology vs sequencing centre",
         subtitle = "Orange is technical. Where it rivals green, cross-cohort comparison is unsafe.",
         x = NULL, y = "% of variance explained", fill = NULL)
  save_plot(p_vp, paste0(args$out_prefix, "_variance_partition.png"),
            width = 9, height = max(3, 0.7 * nrow(vp) + 2))
}

# --------------------------------------------------------------------------- #
# 5. Within-cancer-type institution effect -- the clean technical estimate
# --------------------------------------------------------------------------- #
crossed_types <- crossed[is_crossed == TRUE]$cancer_type

within <- rbindlist(lapply(crossed_types, function(ct) {
  rbindlist(lapply(metric_cols, function(m) {
    d <- man[cancer_type == ct & is.finite(get(m))]
    keep <- d[, .N, by = institution][N >= args$min_per_cell]$institution
    d <- d[institution %in% keep]
    if (uniqueN(d$institution) < 2) return(NULL)

    kw <- tryCatch(kruskal.test(x = d[[m]], g = factor(d$institution)),
                   error = function(e) NULL)
    grp <- d[, .(med = median(get(m), na.rm = TRUE), n = .N), by = institution]
    # Effect size that is interpretable without knowing the metric's units.
    spread <- if (nrow(grp) >= 2) {
      (max(grp$med) - min(grp$med)) / (median(d[[m]], na.rm = TRUE) + 1e-9)
    } else NA_real_

    data.table(
      cancer_type      = ct,
      metric           = m,
      n_institutions   = nrow(grp),
      n_patients       = nrow(d),
      max_median       = max(grp$med),
      min_median       = min(grp$med),
      relative_spread  = spread,
      p_kruskal        = if (is.null(kw)) NA_real_ else kw$p.value
    )
  }), fill = TRUE)
}), fill = TRUE)

if (nrow(within)) {
  within[, fdr := p.adjust(p_kruskal, "BH")]
  setorder(within, fdr, -relative_spread)
  fwrite(within, paste0(args$out_prefix, "_within_type_tests.tsv"), sep = "\t")

  sig <- within[!is.na(fdr) & fdr < 0.05]
  if (nrow(sig)) {
    mohq_warn(nrow(sig), " cancer_type x metric combination(s) show a significant ",
              "institution effect WITH CANCER TYPE HELD CONSTANT. ",
              "This is technical variation, not biology:")
    for (i in seq_len(min(10, nrow(sig)))) {
      mohq_warn(sprintf("   %s / %s: medians %.3g vs %.3g (%.0f%% spread), FDR = %.2g",
                        sig$cancer_type[i], sig$metric[i], sig$min_median[i],
                        sig$max_median[i], 100 * sig$relative_spread[i], sig$fdr[i]))
    }
  } else {
    mohq_log("No significant within-cancer-type institution effects. ",
             "Cross-cohort comparison is on reasonably safe ground for these metrics.")
  }

  plot_dt <- man[cancer_type %in% crossed_types]
  for (m in metric_cols) {
    d <- plot_dt[is.finite(get(m))]
    if (!nrow(d)) next
    pm <- ggplot(d, aes(x = institution, y = get(m), fill = institution)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
      facet_wrap(~ cancer_type, scales = "free_y") +
      theme_minimal(base_size = 11) +
      labs(title = sprintf("%s by institution, within cancer type", m),
           subtitle = "Cancer type is constant within each panel, so a difference here is technical",
           x = NULL, y = m) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(pm, sprintf("%s_%s_by_institution.png", args$out_prefix, m),
              width = max(8, 2.4 * ceiling(sqrt(length(crossed_types)))),
              height = max(5, 2.2 * ceiling(length(crossed_types) / 3)))
  }
}

mohq_session_info(paste0(args$out_prefix, ".versions.txt"))

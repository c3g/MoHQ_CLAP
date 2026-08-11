#!/usr/bin/env Rscript
# =============================================================================
# collection_rollup.R -- collection-level view across all ~65 cohorts
#
# Consumes per-cohort SUMMARY tables only (manifests, completeness matrices,
# Key_metrics.csv), never raw MAFs or VCFs. That is what makes a 4,500-patient
# rollup cheap: the inputs are kilobytes per cohort, so this runs in one small
# task no matter how large the collection grows.
#
# What it answers:
#   * where is the collection incomplete, by cohort and by asset
#   * has the assay drifted over the years / across institutions
#   * which cohorts are outliers on QC and should be looked at before analysis
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparse)
})

parser <- ArgumentParser(description = "Collection-level rollup")
parser$add_argument("--lib", required = TRUE)
parser$add_argument("--manifests", required = TRUE, nargs = "+")
parser$add_argument("--completeness", nargs = "*")
parser$add_argument("--cohort_map", default = NULL,
                    help = "TSV: cohort_id, cancer_type [, accrual_year]")
parser$add_argument("--collection", default = "MoHQ")
parser$add_argument("--out_prefix", default = "collection")
args <- parser$parse_args()

source(args$lib)

# --------------------------------------------------------------------------- #
# 1. Collection inventory
# --------------------------------------------------------------------------- #
man <- rbindlist(lapply(args$manifests, function(f)
  fread(f, sep = "\t", na.strings = c("NA", ""))), fill = TRUE)

if (!is.null(args$cohort_map) && file.exists(args$cohort_map)) {
  cmap <- fread(args$cohort_map, sep = "\t")
  join_cols <- intersect(c("cancer_type", "accrual_year", "sequencing_platform"),
                         names(cmap))
  if (length(join_cols)) {
    man <- merge(man, cmap[, c("cohort_id", join_cols), with = FALSE],
                 by = "cohort_id", all.x = TRUE)
  }
}

mohq_log(sprintf("Collection %s: %d patients | %d cohorts | %d institutions",
                 args$collection, nrow(man), uniqueN(man$cohort_id),
                 uniqueN(man$institution)))

cohort_summary <- man[, .(
  n_patients        = .N,
  n_analysis_ready  = sum(analysis_ready == "yes", na.rm = TRUE),
  pct_ready         = round(100 * mean(analysis_ready == "yes", na.rm = TRUE), 1),
  n_multi_tumour    = sum(multi_tumour == "yes", na.rm = TRUE),
  n_with_rna        = sum(!is.na(expression_genes)),
  pct_with_rna      = round(100 * mean(!is.na(expression_genes)), 1),
  n_with_sv         = sum(!is.na(gripss_somatic)),
  n_with_cpsr       = sum(!is.na(cpsr_zip))
), by = .(institution, cohort_id,
          cancer_type = if ("cancer_type" %in% names(man)) cancer_type else NA_character_)]
setorder(cohort_summary, pct_ready, -n_patients)
fwrite(cohort_summary, paste0(args$out_prefix, "_cohort_summary.tsv"), sep = "\t")

mohq_log(sprintf("Analysis-ready: %d/%d patients (%.1f%%)",
                 sum(man$analysis_ready == "yes", na.rm = TRUE), nrow(man),
                 100 * mean(man$analysis_ready == "yes", na.rm = TRUE)))

worst <- head(cohort_summary[pct_ready < 100], 15)
if (nrow(worst)) {
  mohq_log("Least-complete cohorts:")
  for (i in seq_len(nrow(worst))) {
    mohq_log(sprintf("   %-16s %3d patients, %5.1f%% ready",
                     worst$cohort_id[i], worst$n_patients[i], worst$pct_ready[i]))
  }
}

# --- cohort size distribution: which cohorts can support which analyses ----- #
# This is decision-relevant. Mutual exclusivity and driver discovery need
# n >= ~30-50; signature FITTING works far lower than de-novo extraction.
cohort_summary[, power_tier := fifelse(n_analysis_ready >= 50, "full (drivers, exclusivity)",
                              fifelse(n_analysis_ready >= 20, "moderate (signatures, recurrence)",
                              fifelse(n_analysis_ready >= 10, "limited (descriptive only)",
                                      "underpowered")))]
tier_counts <- cohort_summary[, .(n_cohorts = .N, n_patients = sum(n_analysis_ready)),
                              by = power_tier][order(-n_patients)]
fwrite(tier_counts, paste0(args$out_prefix, "_power_tiers.tsv"), sep = "\t")
mohq_log("Cohorts by achievable analysis depth:")
for (i in seq_len(nrow(tier_counts))) {
  mohq_log(sprintf("   %-36s %3d cohorts, %5d patients",
                   tier_counts$power_tier[i], tier_counts$n_cohorts[i],
                   tier_counts$n_patients[i]))
}

p_size <- ggplot(cohort_summary, aes(x = reorder(cohort_id, n_analysis_ready),
                                     y = n_analysis_ready, fill = institution)) +
  geom_col() +
  geom_hline(yintercept = c(20, 50), linetype = "dashed", colour = "grey40") +
  coord_flip() +
  theme_minimal(base_size = 9) +
  labs(title = sprintf("Analysis-ready patients per cohort: %s", args$collection),
       subtitle = "Dashed lines: n = 20 (signatures, recurrence) and n = 50 (driver discovery, exclusivity)",
       x = NULL, y = "Analysis-ready patients")
save_plot(p_size, paste0(args$out_prefix, "_cohort_sizes.png"),
          width = 10, height = max(5, 0.18 * nrow(cohort_summary) + 2))

# --------------------------------------------------------------------------- #
# 2. Completeness across the collection
# --------------------------------------------------------------------------- #
if (length(args$completeness)) {
  comp <- rbindlist(lapply(args$completeness[file.exists(args$completeness)],
                           function(f) fread(f, sep = "\t")), fill = TRUE)
  id_cols <- intersect(c("cohort_id", "patient_id", "analysis_ready"), names(comp))
  asset_cols <- setdiff(names(comp), id_cols)

  long <- melt(comp, id.vars = id_cols, measure.vars = asset_cols,
               variable.name = "asset", value.name = "present")
  long[, present := present == "yes"]

  by_asset <- long[, .(pct = round(100 * mean(present), 1), n = .N), by = asset]
  setorder(by_asset, pct)
  fwrite(by_asset, paste0(args$out_prefix, "_completeness_by_asset.tsv"), sep = "\t")

  agg <- long[, .(pct = 100 * mean(present)), by = .(cohort_id, asset)]
  ord_c <- agg[, .(m = mean(pct)), by = cohort_id][order(m)]$cohort_id
  ord_a <- by_asset$asset
  agg[, cohort_id := factor(cohort_id, levels = ord_c)]
  agg[, asset := factor(asset, levels = ord_a)]

  p_comp <- ggplot(agg, aes(x = asset, y = cohort_id, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_gradient(low = "#d6d6d6", high = "#1b7837", limits = c(0, 100),
                        name = "% present") +
    theme_minimal(base_size = 9) +
    labs(title = sprintf("Collection completeness: %s", args$collection),
         subtitle = sprintf("%d cohorts x %d assets | cohorts and assets ordered least- to most-complete",
                            uniqueN(agg$cohort_id), length(asset_cols)),
         x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7), panel.grid = element_blank())
  save_plot(p_comp, paste0(args$out_prefix, "_completeness.png"),
            width = max(9, 0.3 * length(asset_cols) + 4),
            height = max(6, 0.18 * uniqueN(agg$cohort_id) + 2))
}

# --------------------------------------------------------------------------- #
# 3. Assay drift from Key_metrics.csv
#
# Every patient directory has one. Concatenating them across 4,500 patients and
# several years is the cheapest possible early-warning system for a changed
# capture kit, a bad run, or a coverage shift -- all of which would otherwise
# surface later as fake biology in a cohort figure.
#
# VERIFY: the exact column names in Key_metrics.csv are not documented here, so
# every numeric column is picked up generically and the useful ones identified
# by name matching. Inspect *_qc_metrics_available.tsv after the first run.
# --------------------------------------------------------------------------- #
km_files <- man[!is.na(key_metrics)]$key_metrics
km_files <- km_files[file.exists(km_files)]

if (length(km_files)) {
  mohq_log(sprintf("Reading %d Key_metrics.csv files", length(km_files)))

  read_km <- function(f) {
    d <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
    if (is.null(d) || !nrow(d)) return(NULL)
    # These files are small and their layout varies; melt to long form so we do
    # not depend on a fixed schema.
    d[, .src := basename(dirname(f))]
    melt(d, id.vars = ".src", variable.name = "metric", value.name = "value",
         variable.factor = FALSE)
  }

  km <- rbindlist(lapply(km_files, read_km), fill = TRUE)
  if (!is.null(km) && nrow(km)) {
    km[, value_num := suppressWarnings(as.numeric(gsub("[,%]", "", as.character(value))))]
    km <- km[!is.na(value_num)]
    km[, patient_id := harmonise_patient_id(.src)]
    km <- merge(km, man[, .(patient_id, institution, cohort_id)],
                by = "patient_id", all.x = TRUE)

    avail <- km[, .(n_patients = uniqueN(patient_id),
                    median = median(value_num, na.rm = TRUE)), by = metric]
    setorder(avail, -n_patients)
    fwrite(avail, paste0(args$out_prefix, "_qc_metrics_available.tsv"), sep = "\t")
    mohq_log("QC metrics found: ", paste(head(avail$metric, 12), collapse = ", "))

    # Focus on the metrics that actually indicate drift.
    focus <- avail[grepl("cover|depth|dup|insert|contam|error|q30|map",
                         metric, ignore.case = TRUE)]$metric
    focus <- head(focus, 8)

    if (length(focus)) {
      kf <- km[metric %in% focus]
      fwrite(kf[, .(patient_id, cohort_id, institution, metric, value_num)],
             paste0(args$out_prefix, "_qc_long.tsv"), sep = "\t")

      p_qc <- ggplot(kf, aes(x = institution, y = value_num, fill = institution)) +
        geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
        facet_wrap(~ metric, scales = "free_y") +
        theme_minimal(base_size = 10) +
        labs(title = sprintf("Sequencing QC by institution: %s", args$collection),
             subtitle = "Systematic offsets here propagate into every downstream cohort comparison",
             x = NULL, y = NULL) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
      save_plot(p_qc, paste0(args$out_prefix, "_qc_by_institution.png"),
                width = 12, height = max(5, 2.4 * ceiling(length(focus) / 3)))

      # Per-cohort QC outliers: robust z against the collection median.
      out_tbl <- kf[, {
        md <- median(value_num, na.rm = TRUE); s <- mad(value_num, na.rm = TRUE)
        .(cohort_median = median(value_num, na.rm = TRUE),
          robust_z = if (is.finite(s) && s > 0)
            (median(value_num, na.rm = TRUE) - md) / s else 0)
      }, by = .(metric, cohort_id)]
      out_tbl <- out_tbl[abs(robust_z) > 3][order(-abs(robust_z))]
      if (nrow(out_tbl)) {
        fwrite(out_tbl, paste0(args$out_prefix, "_qc_outlier_cohorts.tsv"), sep = "\t")
        mohq_warn(nrow(out_tbl), " cohort x metric combination(s) are QC outliers ",
                  "(robust z > 3); see ", args$out_prefix, "_qc_outlier_cohorts.tsv")
      }
    }
  }
} else {
  mohq_warn("No Key_metrics.csv files resolved; skipping the QC drift analysis.")
}

mohq_session_info(paste0(args$out_prefix, ".versions.txt"))

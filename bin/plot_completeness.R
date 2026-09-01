#!/usr/bin/env Rscript
# =============================================================================
# plot_completeness.R -- turn the completeness matrix into a reviewable figure
#
# This is the "checking for completeness of the project" task, promoted from a
# side script to a first-class, versioned pipeline output. Because it runs from
# the same manifest the analyses consume, the completeness figure can never
# disagree with what was actually analysed.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparse)
})

parser <- ArgumentParser(description = "Cohort data-completeness heatmap")
parser$add_argument("--lib", required = TRUE)
parser$add_argument("--completeness", required = TRUE)
parser$add_argument("--manifest", required = TRUE)
parser$add_argument("--cohort", default = "MoHQ")
parser$add_argument("--max_tile_rows", type = "integer", default = 150,
                    help = "Above this many patients, aggregate the heatmap to cohort level")
parser$add_argument("--drop_never_present", type = "logical", default = TRUE,
                    help = paste("Exclude assets absent for EVERY patient from the",
                                 "heatmap and list them in a table instead. These are",
                                 "almost always a harvest-scope choice rather than a",
                                 "per-patient gap, and they crowd out the variation",
                                 "the grid exists to show. (default TRUE)"))
parser$add_argument("--out_prefix", default = "completeness")
args <- parser$parse_args()

source(args$lib)

comp <- fread(args$completeness, sep = "\t")
man  <- read_manifest(args$manifest)

id_cols <- intersect(c("cohort_id", "patient_id", "analysis_ready"), names(comp))
asset_cols <- setdiff(names(comp), id_cols)

long <- melt(comp, id.vars = id_cols, measure.vars = asset_cols,
             variable.name = "asset", value.name = "present")
long[, present := present == "yes"]
long <- merge(long, man[, .(patient_id, institution)], by = "patient_id", all.x = TRUE)

# Per-asset completeness, ordered worst-first: the actionable list.
by_asset <- long[, .(n_present = sum(present), n_total = .N), by = asset]
by_asset[, pct := 100 * n_present / n_total]
setorder(by_asset, pct, asset)
fwrite(by_asset, paste0(args$out_prefix, "_by_asset.tsv"), sep = "\t")

by_patient <- long[, .(n_present = sum(present), n_total = .N), by = .(patient_id, institution)]
by_patient[, pct := 100 * n_present / n_total]
setorder(by_patient, pct)
fwrite(by_patient, paste0(args$out_prefix, "_by_patient.tsv"), sep = "\t")

mohq_log(sprintf("Completeness: %.1f%% of %d patient x asset cells populated",
                 100 * mean(long$present), uniqueN(long$patient_id)))
worst <- by_asset[pct < 100][seq_len(min(8, .N))]
if (nrow(worst)) {
  mohq_log("Least complete assets: ",
           paste(sprintf("%s (%.0f%%)", worst$asset, worst$pct), collapse = ", "))
}

# --------------------------------------------------------------------------- #
# THREE STATES, NOT TWO.
#
# This matrix describes the LOCAL harvested tree, not the object store. Those
# are different questions and the plot did not say which it was answering.
#
# harvest_juno.py transfers a chosen set of file types (--sets, default "core").
# Every other asset is absent locally by design, and the heatmap rendered all of
# them as plain "missing" -- so roughly half the grid was grey and looked like a
# catastrophic delivery failure. It was a download scope.
#
# An asset absent for EVERY patient is almost never a per-patient gap; it is a
# type that was never harvested. Shown as its own state so the eye can separate
# "we did not fetch this" from "this patient is missing something".
#
# The remaining risk is the opposite error: a genuinely collection-wide missing
# asset now looks like a harvest choice. So the label says "check", not
# "not harvested", and the count is reported in the log either way.
# --------------------------------------------------------------------------- #
never_present <- by_asset[n_present == 0]$asset
long[, state := fifelse(asset %in% never_present, "absent for all",
                fifelse(present, "present", "missing"))]
if (length(never_present)) {
  mohq_log(sprintf(paste0("%d asset(s) are absent for EVERY patient: %s.\n",
                          "  Most likely they were not harvested (see harvest_juno.py ",
                          "--sets) rather than missing from the delivery. Confirm ",
                          "against the object store before reporting them as gaps."),
                   length(never_present),
                   paste(utils::head(as.character(never_present), 10), collapse = ", ")))
}

# --------------------------------------------------------------------------- #
# DROP THE NEVER-PRESENT ASSETS FROM THE TILES.
#
# Giving them their own colour made them interpretable but not readable: with
# ~24 of ~60 assets never harvested, a third of the grid is a solid block that
# carries one bit of information ("we did not fetch these") and crowds out the
# per-patient variation the figure exists to show.
#
# So they come out of the plot and go into a table instead -- same information,
# proportionate space. The heatmap is then about what differs BETWEEN patients,
# which is the only thing a per-patient grid can tell you.
#
# --drop_never_present false restores the old behaviour, for the case where you
# genuinely want to see the full expected asset list.
# --------------------------------------------------------------------------- #
if (isTRUE(args$drop_never_present) && length(never_present)) {
  fwrite(by_asset[n_present == 0, .(asset, n_present, n_total, pct)],
         paste0(args$out_prefix, "_assets_never_present.tsv"), sep = "\t")
  long      <- long[!asset %in% never_present]
  by_asset  <- by_asset[n_present > 0]
  asset_cols <- setdiff(asset_cols, as.character(never_present))
  mohq_log(sprintf("Heatmap shows %d asset(s) with at least one patient; %d never-present asset(s) listed in %s_assets_never_present.tsv",
                   length(asset_cols), length(never_present), args$out_prefix))
  if (!nrow(long))
    mohq_die("Every asset is absent for every patient -- nothing to plot. ",
             "Check the harvest before interpreting this cohort.")
}

long[, asset := factor(asset, levels = by_asset$asset)]
long[, patient_id := factor(patient_id, levels = by_patient$patient_id)]

n_pat <- uniqueN(long$patient_id)
w <- max(8, 0.32 * length(asset_cols) + 3)

# SCALE GUARD: one tile row per patient is unreadable beyond ~150 patients, and
# at 4,500 the naive height (0.22 in/patient) would be a 990-inch image that
# ggsave refuses or that no one can open. Above the threshold we aggregate to
# cohort level; the per-patient detail is still in the TSVs.
if (n_pat <= args$max_tile_rows) {
  p <- ggplot(long, aes(x = asset, y = patient_id, fill = state)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    facet_grid(institution ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c(present = "#1b7837",
                                 missing = "#d6d6d6",
                                 `absent for all` = "#fde0c5"),
                      breaks = c("present", "missing", "absent for all"),
                      labels = c("present", "missing for this patient",
                                 "absent for ALL patients - check if harvested")) +
    theme_minimal(base_size = 10) +
    labs(title = sprintf("Data completeness: %s", args$cohort),
         subtitle = sprintf(paste0("%d patients x %d assets | least- to most-complete | ",
                                   "scope: files present in the LOCAL harvested tree, ",
                                   "not the object store"),
                            n_pat, length(asset_cols)),
         x = NULL, y = NULL, fill = NULL) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
          axis.text.y = element_text(size = 7),
          panel.grid = element_blank(), legend.position = "top",
          strip.text.y = element_text(angle = 0, face = "bold"))
  h <- max(5, 0.22 * n_pat + 2)
} else {
  mohq_log(sprintf("%d patients > --max_tile_rows (%d); aggregating heatmap to cohort level.",
                   n_pat, args$max_tile_rows))
  agg <- long[, .(pct = 100 * mean(present), n = uniqueN(patient_id)),
              by = .(cohort_id, asset)]
  ord <- agg[, .(m = mean(pct)), by = cohort_id][order(m)]$cohort_id
  agg[, cohort_id := factor(cohort_id, levels = ord)]

  p <- ggplot(agg, aes(x = asset, y = cohort_id, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient(low = "#d6d6d6", high = "#1b7837", limits = c(0, 100),
                        name = "% present") +
    theme_minimal(base_size = 10) +
    labs(title = sprintf("Data completeness: %s", args$cohort),
         subtitle = sprintf("%d patients in %d cohorts x %d assets | per-patient detail in the TSVs",
                            n_pat, uniqueN(long$cohort_id), length(asset_cols)),
         x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
          axis.text.y = element_text(size = 7),
          panel.grid = element_blank(), legend.position = "top")
  h <- max(5, 0.22 * uniqueN(agg$cohort_id) + 2)

  # The worst individual patients still need to be nameable.
  worst_pat <- head(by_patient[order(pct)], 50)
  fwrite(worst_pat, paste0(args$out_prefix, "_worst_patients.tsv"), sep = "\t")
  mohq_log("Wrote the 50 least-complete patients to ",
           args$out_prefix, "_worst_patients.tsv")
}

save_plot(p, paste0(args$out_prefix, "_heatmap.png"), width = w, height = h)

pa <- ggplot(by_asset, aes(x = reorder(asset, pct), y = pct)) +
  geom_col(fill = "grey35") +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey60") +
  coord_flip() +
  theme_minimal(base_size = 11) +
  labs(title = "Completeness by asset", x = NULL, y = "% of patients with the file")
save_plot(pa, paste0(args$out_prefix, "_by_asset.png"), width = 8, height = w)

mohq_session_info(paste0(args$out_prefix, ".versions.txt"))

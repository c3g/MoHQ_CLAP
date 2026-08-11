#!/usr/bin/env Rscript
# =============================================================================
# run_fga_burden.R -- Fraction of Genome Altered per patient
#
# Answers the "# Check if there is a better way than to hardcode it??" comment
# in the original script. Yes: do not use a fixed genome size at all.
#
# The old denominator was a constant (2875 Mb here, 2744 Mb in the comparative
# script -- two values for one metric). A constant denominator assumes every
# sample was profiled across the whole genome. When that is not true, samples
# with a smaller callable footprint look artificially quiet, and because
# callable footprint tracks coverage and capture batch, the metric acquires a
# batch effect that looks like biology.
#
# FGA is now altered length / length actually profiled FOR THAT SAMPLE
# (the cBioPortal definition). `profiled_mb` is reported alongside so the
# footprint is visible and QC-able instead of being absorbed into the number.
#
# Thresholds now come from the same --amp/--del parameters as every other
# script, rather than a second hardcoded +/-0.2.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(argparse)
})

parser <- ArgumentParser(description = "Per-patient FGA burden")
parser$add_argument("--lib", required = TRUE)
parser$add_argument("--seg", required = TRUE)
parser$add_argument("--manifest", required = TRUE)
parser$add_argument("--cohort", default = "MoHQ Cohort")
parser$add_argument("--amp", type = "double", default = 0.58)
parser$add_argument("--del", type = "double", default = -0.58)
parser$add_argument("--min_profiled_mb", type = "double", default = 1000)
parser$add_argument("--out_prefix", default = "cohort_cnv_burden")
args <- parser$parse_args()

source(args$lib)
manifest <- read_manifest(args$manifest)

seg <- read_seg(args$seg, manifest = manifest)
fga <- compute_fga(seg, args$amp, args$del,
                   autosomes_only = TRUE,
                   min_profiled_mb = args$min_profiled_mb)

# Attach institution so batch structure is visible in the burden plot -- with
# multiple contributing institutions this is the first place a systematic
# difference in CNV calling will show up.
fga <- merge(fga, manifest[, .(patient_id, institution, cohort_id)],
             by = "patient_id", all.x = TRUE)

fwrite(fga, paste0(args$out_prefix, ".tsv"), sep = "\t")
mohq_log(sprintf("FGA: median %.1f%%, IQR %.1f-%.1f%%",
                 100 * median(fga$fga), 100 * quantile(fga$fga, .25),
                 100 * quantile(fga$fga, .75)))

# --- implausible values ----------------------------------------------------- #
#
# FGA at or near 1.0 means every profiled segment in the autosome was called
# altered. That is not a biological state; it is a normalisation failure -- a
# log2 baseline offset, or a purity/ploidy fit that placed the whole genome off
# neutral. The first real run had two patients at exactly 100% loss, and those
# two also drove the near-100% loss bands in the genome-wide frequency plot.
#
# Flagged rather than dropped: which patients are affected is information, and
# silently removing them would hide a systematic problem.
fga[, implausible := fga >= 0.99]
n_bad <- sum(fga$implausible)
if (n_bad) {
  mohq_warn(sprintf(paste0(
    "%d patient(s) have FGA >= 99%%, which is not biologically plausible and ",
    "indicates a copy-number normalisation failure: %s\n",
    "  They are flagged in the output table and marked on the plot. They also ",
    "inflate the genome-wide frequency plot, so check that figure against the ",
    "cohort with these patients excluded before interpreting it."),
    n_bad, paste(fga$patient_id[fga$implausible], collapse = ", ")))
}
# Zero is a legitimate result for copy-number-quiet tumours, but a large block
# of exactly zero is worth a look rather than an assumption.
n_zero <- sum(fga$fga == 0)
if (n_zero > 0.5 * nrow(fga)) {
  mohq_warn(sprintf(paste0(
    "%d of %d patients (%.0f%%) have FGA of exactly zero. Some tumour types are ",
    "genuinely copy-number quiet, but check that segments were actually loaded ",
    "for these patients and that the +/-%.2f log2 thresholds suit this assay."),
    n_zero, nrow(fga), 100 * n_zero / nrow(fga), args$amp))
}

# --- stacked gain/loss waterfall ------------------------------------------- #
long <- melt(fga[, .(patient_id, institution, low_coverage, implausible,
                     Gain = 100 * fga_gain, Loss = 100 * fga_loss)],
             id.vars = c("patient_id", "institution", "low_coverage", "implausible"),
             variable.name = "type", value.name = "pct")

order_by <- fga[order(-fga)]$patient_id
long[, patient_id := factor(patient_id, levels = order_by)]

p <- ggplot(long, aes(x = patient_id, y = pct, fill = type)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.1) +
  scale_fill_manual(values = MOHQ_PALETTE[c("Gain", "Loss")]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     labels = function(x) paste0(x, "%")) +
  theme_minimal(base_size = 12) +
  labs(
    title = sprintf("Somatic copy-number burden: %s (n = %d)",
                    args$cohort, nrow(fga)),
    subtitle = sprintf("FGA = altered autosomal length / profiled autosomal length, per patient | gain >= %.2f, loss <= %.2f log2",
                       args$amp, args$del),
    x = "Patient", y = "Fraction of Genome Altered"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "top", legend.title = element_blank(),
    plot.title = element_text(face = "bold")
  )

# Mark patients whose profiled footprint is too small to compare fairly, and
# patients whose FGA is not biologically possible. Two different problems, so
# two different marks rather than one ambiguous asterisk.
caps <- character()
if (any(fga$low_coverage)) {
  flagged <- fga[low_coverage == TRUE]$patient_id
  p <- p + annotate("text", x = flagged, y = 0, label = "!", vjust = 1.6,
                    colour = "red", fontface = "bold", size = 4)
  caps <- c(caps, sprintf("! = profiled over < %.0f Mb of autosome; interpret with caution",
                          args$min_profiled_mb))
}
if (any(fga$implausible)) {
  bad <- fga[implausible == TRUE]$patient_id
  p <- p + annotate("text", x = bad, y = 100, label = "✖", vjust = -0.2,
                    colour = "#b91c1c", fontface = "bold", size = 4)
  caps <- c(caps, paste("X = FGA >= 99%, not biologically plausible.",
                        "Copy-number normalisation failed for this patient;",
                        "exclude it before interpreting the genome-wide frequency plot."))
}
if (length(caps)) p <- p + labs(caption = paste(caps, collapse = "\n"))

save_plot(p, paste0(args$out_prefix, ".png"), width = 11, height = 5.5)

# --- FGA by institution: an explicit batch-effect check -------------------- #
if (uniqueN(fga$institution) > 1) {
  pb <- ggplot(fga, aes(x = institution, y = 100 * fga)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
    theme_minimal(base_size = 13) +
    labs(title = "FGA by contributing institution",
         subtitle = "A systematic shift here is a batch effect until proven otherwise",
         x = "Institution", y = "Fraction of Genome Altered (%)")
  save_plot(pb, paste0(args$out_prefix, "_by_institution.png"), width = 6, height = 5)

  if (uniqueN(fga$institution) >= 2) {
    # `fga ~ ...` resolves to the fga COLUMN of the fga data.table; explicit
    # vectors avoid depending on that coincidence.
    kw <- kruskal.test(x = fga$fga, g = factor(fga$institution))
    mohq_log(sprintf("Kruskal-Wallis FGA ~ institution: p = %.3g", kw$p.value))
    fwrite(data.table(test = "kruskal_fga_by_institution",
                      statistic = unname(kw$statistic), p_value = kw$p.value),
           paste0(args$out_prefix, "_institution_test.tsv"), sep = "\t")
  }
}

mohq_session_info(paste0(args$out_prefix, ".versions.txt"))

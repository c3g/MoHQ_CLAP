#!/usr/bin/env Rscript
# =============================================================================
# run_cnv_frequency.R -- genome-wide CNV gain/loss frequency across a cohort
#
# PRIMARY FIX vs the previous version
# -----------------------------------
# The old script binned segments with
#     mutate(bin = floor(start / bin_size) * bin_size)
# which places each segment in its START bin only. A segment spanning 50 Mb
# contributed to one 1 Mb bin instead of fifty. Since large segments carry the
# bulk of the copy-number signal, the plot systematically under-reported the
# broad events it was built to display -- and did so silently, because the
# output still looked like a plausible frequency plot.
#
# Segments are now expanded across every bin they overlap
# (see expand_segments_to_bins in mohq_common.R).
#
# Secondary fixes:
#   * frequency is now a PERCENTAGE of patients profiled at that bin, not a raw
#     count, so bins with different denominators are comparable
#   * columns resolved by name, not by blind positional renaming
#   * chromosomes ordered numerically; sex chromosomes optional
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(argparse)
})

parser <- ArgumentParser(description = "Cohort genome-wide CNV frequency")
parser$add_argument("--lib", required = TRUE)
parser$add_argument("--seg", required = TRUE)
parser$add_argument("--manifest", required = TRUE)
parser$add_argument("--cohort", default = "MoHQ Cohort")
parser$add_argument("--amp", type = "double", default = 0.58)
parser$add_argument("--del", type = "double", default = -0.58)
parser$add_argument("--bin_size", type = "double", default = 1e6)
parser$add_argument("--include_sex_chr", action = "store_true",
                    help = "Include chrX/chrY (off by default: without sex-matched normalisation these bins are dominated by patient sex, not by tumour biology)")
parser$add_argument("--out_prefix", default = "cohort_cnv_frequency")
args <- parser$parse_args()

source(args$lib)
manifest <- read_manifest(args$manifest)

seg <- read_seg(args$seg, manifest = manifest)
seg <- classify_segments(seg, args$amp, args$del)

chr_levels <- if (args$include_sex_chr) c(AUTOSOMES, "X", "Y") else AUTOSOMES
seg <- seg[chr %in% chr_levels]
if (!nrow(seg)) mohq_die("No segments on the requested chromosomes.")

n_patients <- uniqueN(seg$patient_id)
mohq_log(sprintf("Cohort %s: %d patients", args$cohort, n_patients))

# --- denominator: how many patients are actually profiled at each bin? ------ #
# Using the cohort size as a flat denominator would understate frequency in
# regions where some patients have no segment coverage at all.
covered <- expand_segments_to_bins(seg, args$bin_size)
denom <- covered[, .(n_profiled = uniqueN(patient_id)), by = .(chr, bin)]

altered <- covered[type != "Neutral"]
freq <- altered[, .(n_patients = uniqueN(patient_id)), by = .(chr, bin, type)]
freq <- merge(freq, denom, by = c("chr", "bin"), all.x = TRUE)
freq[, pct := 100 * n_patients / n_profiled]
# Losses point left of the axis.
freq[type == "Loss", pct := -pct]
freq[, chr := factor(chr, levels = chr_levels)]

fwrite(freq[order(chr, bin)], paste0(args$out_prefix, ".tsv"), sep = "\t")

# --- pad each chromosome to its true length -------------------------------- #
#
# facet_grid(space = "free_y") allocates vertical space by the DATA RANGE in
# each facet, not by chromosome length. Chromosomes with alterations spread
# across many bins therefore got taller panels than chromosomes with few, and
# the first real run drew chr6 and chr16 larger than chr1 -- which invites
# exactly the wrong reading of a genome-wide figure.
#
# Adding one zero-valued row at each chromosome's first and last bin fixes the
# range without changing any plotted value.
pad <- rbindlist(lapply(levels(freq$chr), function(cc) {
  len <- GRCH38_CHROM_LEN[[as.character(cc)]]
  if (is.null(len)) return(NULL)
  data.table(chr = cc, bin = c(0, floor(len / args$bin_size) * args$bin_size),
             type = "Gain", n_patients = 0, n_profiled = NA_integer_, pct = 0)
}), fill = TRUE)
if (nrow(pad)) {
  pad[, chr := factor(chr, levels = chr_levels)]
  freq <- rbind(freq, pad, fill = TRUE)
} else {
  mohq_warn("No chromosome lengths matched; panel heights will follow the data range.")
}

# --- plot ------------------------------------------------------------------ #
p <- ggplot(freq, aes(y = bin, x = pct, fill = type)) +
  geom_col(width = args$bin_size, orientation = "y") +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
  facet_grid(chr ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(labels = function(x) paste0(abs(x), "%")) +
  scale_fill_manual(values = MOHQ_PALETTE[c("Gain", "Loss")]) +
  theme_classic(base_size = 13) +
  labs(
    title = sprintf("CNV landscape: %s (n = %d)", args$cohort, n_patients),
    subtitle = sprintf("%.0f kb bins | gain >= %.2f, loss <= %.2f log2 | %s",
                       args$bin_size / 1e3, args$amp, args$del,
                       if (args$include_sex_chr) "incl. sex chromosomes" else "autosomes only"),
    x = "Patients altered (%)", y = "Chromosome"
  ) +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90", linetype = "dashed"),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    strip.placement = "outside",
    panel.spacing = unit(0, "lines"),
    panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.4),
    legend.position = "top", legend.title = element_blank()
  )

save_plot(p, paste0(args$out_prefix, ".png"), width = 8, height = 14)
mohq_session_info(paste0(args$out_prefix, ".versions.txt"))

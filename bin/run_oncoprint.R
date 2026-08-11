#!/usr/bin/env Rscript
# =============================================================================
# run_oncoprint.R -- cohort oncoplot, TMB, mutual exclusivity, pathway summary
#
# Fixes vs the previous run_oncoprints.R
# --------------------------------------
# 1. RNA MAFs are no longer swept in. `list.files(pattern="\\.maf$")` matched
#    both *_D.*.maf (somatic DNA) and *_R.*.maf (RNA-derived). Those are
#    different assays; the RNA files are also 10-20x larger, which is where
#    most of the memory pressure came from. File selection is now the
#    manifest's job and load_somatic_mafs() refuses RNA files outright.
#
# 2. The PCGR driver filter compared a column to the string "TRUE". fread
#    parses those columns as logical, so `oncogene == "TRUE"` was FALSE for
#    every row and every CNV was silently discarded -- an empty cnTable, no
#    error, an oncoplot with no copy-number track. See read_pcgr_cna().
#
# 3. maftools expects cnTable columns Gene / Sample_name / CN; the old code
#    passed Hugo_Symbol / Tumor_Sample_Barcode / Variant_Classification.
#
# 4. Variants are filtered to non-synonymous classes. Without this the "top 20
#    mutated genes" is largely TTN / MUC16 / OBSCN -- long genes, not drivers.
#
# 5. Sample barcodes are canonicalised to patient_id, so one patient with two
#    tumour samples counts once rather than twice in every frequency.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(maftools)
  library(argparse)
})

parser <- ArgumentParser(description = "Cohort oncoplot and mutation summaries")
parser$add_argument("--lib", required = TRUE)
parser$add_argument("--mafs", required = TRUE, nargs = "+", help = "Somatic (_D) MAF files")
parser$add_argument("--cna", nargs = "*",
                    help = "PCGR *.cna_segments.tsv.gz files")
parser$add_argument("--manifest", required = TRUE)
parser$add_argument("--cohort", default = "MoHQ Cohort")
parser$add_argument("--amp", type = "double", default = 0.58)
parser$add_argument("--del", type = "double", default = -0.58)
parser$add_argument("--top", type = "integer", default = 25)
parser$add_argument("--callable_mb", type = "double", default = 30,
                    help = paste("Territory the TMB NUMERATOR was counted over, in Mb.",
                                 "This script counts NON-SYNONYMOUS (coding) variants, so the",
                                 "denominator is the CODING footprint (~30-40 Mb) even for WGS.",
                                 "Use ~2800 only if you also count non-coding variants."))
parser$add_argument("--exclude_flags", type = "logical", default = TRUE,
                    help = "Drop FLAGS artefact-prone genes before ranking (default TRUE)")
parser$add_argument("--extra_flags", default = "",
                    help = "Comma-separated extra gene symbols to exclude")
parser$add_argument("--gene_panel", default = "",
                    help = paste("Optional file of HGNC symbols, one per line. If given,",
                                 "the oncoplot is restricted to these genes -- an allowlist",
                                 "(e.g. COSMIC Cancer Gene Census) rather than a blocklist."))
parser$add_argument("--keep_synonymous", action = "store_true")
parser$add_argument("--out_prefix", default = "cohort_oncoprint")
args <- parser$parse_args()

source(args$lib)
manifest <- read_manifest(args$manifest)

# --------------------------------------------------------------------------- #
# Load
# --------------------------------------------------------------------------- #
maf_dt <- load_somatic_mafs(args$mafs, manifest = manifest,
                            nonsyn_only = !args$keep_synonymous)

# --------------------------------------------------------------------------- #
# FLAGS filter.
#
# The unfiltered ranking is written to disk first, so the decision is auditable
# and reversible rather than invisible. Then the artefact-prone genes are
# dropped, because without this the top-25 list is mucins, WASH/NBPF/GOLGA
# family members and FLG -- the difficult regions of the reference genome, not
# the biology of the cohort.
# --------------------------------------------------------------------------- #
rank_all <- maf_dt[, .(n_patients = uniqueN(patient_id), n_variants = .N),
                   by = Hugo_Symbol][order(-n_patients)]
rank_all[, is_flags := Hugo_Symbol %in% FLAGS_GENES]
fwrite(rank_all, paste0(args$out_prefix, "_gene_rank_unfiltered.tsv"), sep = "\t")

if (args$exclude_flags) {
  n_before <- nrow(maf_dt)
  maf_dt <- drop_flags(maf_dt, extra = if (nzchar(args$extra_flags))
                                         trimws(strsplit(args$extra_flags, ",")[[1]])
                                       else character())
  mohq_log(sprintf("FLAGS filter removed %s of %s variants (%.1f%%)",
                   format(n_before - nrow(maf_dt), big.mark = ","),
                   format(n_before, big.mark = ","),
                   100 * (n_before - nrow(maf_dt)) / max(n_before, 1)))
  if (!nrow(maf_dt))
    mohq_die("Every variant was removed by the FLAGS filter. Re-run with --exclude_flags false and inspect ",
             paste0(args$out_prefix, "_gene_rank_unfiltered.tsv"))
} else {
  mohq_warn("FLAGS filter is OFF. Expect the top genes to be dominated by ",
            "recurrently mis-called artefact genes (mucins, WASH/NBPF/GOLGA, FLG).")
}

# An allowlist, if supplied, is applied AFTER the blocklist and supersedes it in
# practice: it answers "which genes am I making a claim about?" rather than the
# open-ended "which genes do I not believe?".
if (nzchar(args$gene_panel)) {
  maf_dt <- restrict_to_panel(maf_dt, args$gene_panel)
} else {
  mohq_warn("No --gene_panel supplied. The oncoplot ranks ALL remaining genes, ",
            "so the top of the list is whatever survived the artefact filter. ",
            "For a figure anyone will interpret, restrict to a curated list ",
            "(COSMIC Cancer Gene Census or an OncoKB panel).")
}

cn_dt <- read_pcgr_cna(args$cna, args$amp, args$del,
                       driver_only = TRUE, manifest = manifest)

# Clinical annotation drives the oncoplot side bars.
clin <- manifest[patient_id %in% unique(maf_dt$patient_id),
                 .(Tumor_Sample_Barcode = patient_id, institution, cohort_id,
                   multi_tumour)]

maf_obj <- read.maf(
  maf = as.data.frame(maf_dt),
  cnTable = if (nrow(cn_dt)) as.data.frame(cn_dt) else NULL,
  clinicalData = as.data.frame(clin),
  verbose = FALSE
)

# --------------------------------------------------------------------------- #
# Oncoplot
# --------------------------------------------------------------------------- #
n_patients <- uniqueN(maf_dt$patient_id)

png(paste0(args$out_prefix, ".png"), width = 14, height = 9, units = "in", res = 300)
oncoplot(maf = maf_obj, top = args$top,
         clinicalFeatures = intersect(c("institution", "cohort_id"), names(clin)),
         sortByAnnotation = TRUE, draw_titv = TRUE,
         titleText = sprintf("%s  (n = %d patients)", args$cohort, n_patients))
dev.off()
mohq_log("Wrote ", args$out_prefix, ".png")

# --------------------------------------------------------------------------- #
# Gene summary + TMB
# --------------------------------------------------------------------------- #
gs <- as.data.table(getGeneSummary(maf_obj))
fwrite(gs, paste0(args$out_prefix, "_gene_summary.tsv"), sep = "\t")

tmb <- compute_tmb(maf_dt, callable_mb = args$callable_mb)
tmb <- merge(tmb, manifest[, .(patient_id, institution)], by = "patient_id", all.x = TRUE)
fwrite(tmb, paste0(args$out_prefix, "_tmb.tsv"), sep = "\t")
mohq_log(sprintf("TMB: median %.2f mut/Mb (assuming %.0f Mb callable)",
                 median(tmb$tmb_per_mb), args$callable_mb))

# Denominator sanity check. The numerator here is CODING variants, so a
# whole-genome denominator understates TMB by roughly two orders of magnitude.
# The first real run used 2800 Mb and produced 0.01-0.05 mut/Mb; the same data
# over the coding footprint gives 0.9-4.7, which is a believable range.
if (!args$keep_synonymous && args$callable_mb > 100) {
  mohq_warn(sprintf(paste0(
    "callable_mb = %.0f with a NON-SYNONYMOUS numerator. Coding variants divided ",
    "by a whole-genome territory understates TMB by ~%.0fx. Use ~30-40 Mb, or ",
    "pass --keep_synonymous to count genome-wide."), args$callable_mb, args$callable_mb / 35))
}

# Linear scale, NOT log.
#
# geom_col draws each bar from zero, and log10(0) is -Inf, so on a log scale
# ggplot draws every bar from the top of the panel downwards. The first real run
# produced exactly that: 20 bars apparently pinned at 1.0 with a small notch,
# which reads as "every patient has the same TMB". The values here span less than
# one order of magnitude, so a linear axis loses nothing.
p_tmb <- ggplot2::ggplot(tmb, ggplot2::aes(x = reorder(patient_id, -tmb_per_mb),
                                           y = tmb_per_mb, fill = institution)) +
  ggplot2::geom_col() +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::labs(title = paste("Tumour mutational burden:", args$cohort),
                subtitle = sprintf("non-synonymous variants / %.0f Mb (median %.2f/Mb, n = %d)",
                                   args$callable_mb, median(tmb$tmb_per_mb), nrow(tmb)),
                x = "Patient", y = "Mutations per Mb") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
                 panel.grid.major.x = ggplot2::element_blank())
save_plot(p_tmb, paste0(args$out_prefix, "_tmb.png"), width = 11, height = 5)

# --------------------------------------------------------------------------- #
# Cohort-scale analyses that are impossible per-patient
# --------------------------------------------------------------------------- #

# Mutual exclusivity / co-occurrence. Needs a reasonable cohort size to mean
# anything, so it is skipped rather than reported spuriously on tiny cohorts.
# Your CM-4 / CQ-34 / IQ-33 cohorts are ~5 patients each, so this will skip
# until cohorts are pooled -- which is itself an argument for pooling.
n_samples <- n_patients
if (n_samples >= 15) {
  png(paste0(args$out_prefix, "_somatic_interactions.png"),
      width = 9, height = 8, units = "in", res = 300)
  si <- tryCatch(
    somaticInteractions(maf = maf_obj, top = args$top, pvalue = c(0.05, 0.01)),
    error = function(e) { mohq_warn("somaticInteractions failed: ", e$message); NULL })
  dev.off()
  if (!is.null(si)) {
    fwrite(as.data.table(si), paste0(args$out_prefix, "_somatic_interactions.tsv"), sep = "\t")
  }
} else {
  mohq_warn(sprintf("Only %d patients; skipping mutual-exclusivity testing ",
                    n_samples),
            "(underpowered -- any 'significant' pair would be noise).")
}

# Oncogenic pathway involvement (Sanchez-Vega curated pathways).
# maftools renamed OncogenicPathways() -> pathways() around v2.12, so try both
# rather than depending on one version.
pw <- tryCatch({
  res <- if (exists("pathways")) {
    pathways(maf = maf_obj, plotType = "treemap")
  } else {
    OncogenicPathways(maf = maf_obj)
  }
  fwrite(as.data.table(res), paste0(args$out_prefix, "_pathways.tsv"), sep = "\t")
  res
}, error = function(e) {
  mohq_warn("Oncogenic pathway analysis unavailable: ", e$message); NULL
})

# Per-patient variant classification summary, useful for QC drift over years.
png(paste0(args$out_prefix, "_summary.png"), width = 12, height = 7,
    units = "in", res = 300)
plotmafSummary(maf = maf_obj, addStat = "median", dashboard = TRUE)
dev.off()

mohq_session_info(paste0(args$out_prefix, ".versions.txt"))

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
parser$add_argument("--oncodrive", type = "logical", default = TRUE,
                    help = paste("Rank genes by POSITIONAL CLUSTERING as well as by",
                                 "frequency (maftools::oncodrive, the OncodriveCLUST",
                                 "method). Frequency alone favours long genes; this",
                                 "asks whether mutations pile up at specific residues.",
                                 "(default TRUE)"))
parser$add_argument("--oncodrive_min_mut", type = "integer", default = 5,
                    help = "Minimum mutations in a gene before it is tested (default 5)")
parser$add_argument("--gtf", default = "",
                    help = paste("GRCh38 GTF. When supplied, the gene ranking gains a",
                                 "coding length and a mutations-per-kb column, which",
                                 "removes the length bias that makes long genes",
                                 "outrank drivers on raw frequency."))
parser$add_argument("--out_prefix", default = "cohort_oncoprint")
args <- parser$parse_args()

source(args$lib)
manifest <- read_manifest(args$manifest)

# --------------------------------------------------------------------------- #
# Load
# --------------------------------------------------------------------------- #
# ASK FOR THE PROTEIN-CHANGE COLUMN.
#
# load_somatic_mafs() keeps MAF_MIN_COLS plus whatever extra_cols requests, and
# discards the rest -- sensible, since these MAFs carry 116 columns and only a
# handful are used. But oncodrive needs an amino-acid position, so without
# asking for it here the clustering test reports "no amino-acid change column
# found" against files that have HGVSp_Short sitting in column 37.
#
# The column is requested by NAME and intersected against each file's header,
# so a MAF that genuinely lacks it still loads -- the clustering test then
# skips with a clear message rather than the loader failing.
maf_dt <- load_somatic_mafs(args$mafs, manifest = manifest,
                            nonsyn_only = !args$keep_synonymous,
                            extra_cols  = c("HGVSp_Short", "HGVSp",
                                            "Protein_position", "Amino_acids"))

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
# is_flags must use the SAME rule as drop_flags(), which matches the named list
# OR the family patterns. Testing only the named list marked AGAP4, RGPD8 and
# every GOLGA/NBPF member as FALSE while still removing them -- so this table,
# whose whole purpose is to show what the filter did, disagreed with the filter.
rank_all[, is_flags := Hugo_Symbol %in% FLAGS_GENES |
                       grepl(paste(FLAGS_PATTERNS, collapse = "|"), Hugo_Symbol)]
# LENGTH-NORMALISED RANKING.
#
# The available background correction. Positional clustering (oncodrive) needs
# a protein-change column, which these MAFs do not carry -- vcf2maf run with
# --inhibit-vep does not emit HGVSp_Short. Coding length needs only a GTF, and
# it targets the bias directly: divide by how big a target the gene is.
#
# Reported ALONGSIDE the raw count, never instead of it. Dividing by length
# gives short genes a small-denominator advantage, so a two-mutation 0.5 kb
# gene can top the list on noise. min_kb keeps those out of the ranking rather
# than pretending the ratio is stable there.
cds <- if (nzchar(args$gtf)) cds_length_kb(args$gtf) else NULL
if (!is.null(cds)) {
  rank_all <- merge(rank_all, cds, by.x = "Hugo_Symbol", by.y = "gene", all.x = TRUE)
  rank_all[, mut_per_kb := fifelse(is.na(cds_kb) | cds_kb <= 0,
                                   NA_real_, n_variants / cds_kb)]
  min_kb <- 0.5
  rank_all[, rank_per_kb := fifelse(!is.na(mut_per_kb) & cds_kb >= min_kb,
                                    frank(-mut_per_kb, ties.method = "first"),
                                    NA_real_)]
  setorder(rank_all, -n_patients)

  n_missing <- sum(is.na(rank_all$cds_kb))
  if (n_missing)
    mohq_warn(n_missing, " gene symbol(s) had no CDS in the GTF and carry no ",
              "length-normalised rank. Usually a symbol-version mismatch ",
              "between the MAF's annotation and the GTF release.")

  top_raw <- head(rank_all[order(-n_patients)]$Hugo_Symbol, 6)
  top_norm <- head(rank_all[!is.na(rank_per_kb)][order(rank_per_kb)]$Hugo_Symbol, 6)
  mohq_log("top by raw count      : ", paste(top_raw, collapse = ", "))
  mohq_log("top by mutations/kb   : ", paste(top_norm, collapse = ", "))
} else if (nzchar(args$gtf)) {
  mohq_warn("--gtf supplied but no coding lengths could be read from it.")
}

fwrite(rank_all, paste0(args$out_prefix, "_gene_rank_unfiltered.tsv"), sep = "\t")

# KEEP THE UNFILTERED VARIANTS FOR THE CLUSTERING TEST.
#
# oncodrive must see the SAME input the frequency ranking saw, or it answers
# nothing. Run it downstream of drop_flags() and it can only re-rank genes the
# blocklist already approved: it cannot demote one that was removed (pointless)
# and cannot rescue one removed in error (harmful). Under --gene_panel it would
# merely re-rank the panel.
#
# The entire claim being made is "frequency puts RP1 first, clustering puts
# BRAF first, on identical input, with no blocklist involved". That claim only
# holds if the clustering test runs BEFORE any filtering.
#
# Cost: a second read.maf over the full variant set. On a hypermutated cohort
# that is real memory, so it is built only when the test is switched on.
maf_dt_unfiltered <- if (isTRUE(args$oncodrive)) copy(maf_dt) else NULL

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

# --------------------------------------------------------------------------- #
# RANK BY CLUSTERING, NOT ONLY BY FREQUENCY.
#
# The problem this addresses: a long gene is a bigger target, so it collects
# more passenger mutations and rises up a frequency ranking without being a
# driver. Filtering the offenders by name does not converge -- on MoHQ-CM-6,
# after two rounds of artefact filtering, BRAF still ranked third behind RP1
# and MGAM2, both of which are simply large.
#
# oncodrive (the OncodriveCLUST method) asks a different question: are this
# gene's mutations CONCENTRATED at particular residues, more than chance? A
# passenger process scatters mutations along the coding sequence; positive
# selection piles them on the same few positions. Gene length cancels out,
# because the test is about distribution within the gene rather than count.
#
# WHAT IT WILL AND WILL NOT FIND -- state this when reporting it.
#   Finds:    hotspot oncogenes. BRAF V600, NRAS Q61, KRAS G12, PIK3CA E545.
#   Misses:   tumour suppressors. CDKN2A, NF1, TP53 are inactivated by
#             truncating mutations spread across the gene, which is exactly
#             the pattern this test calls "not clustered".
#
# So it is a COMPLEMENT to the frequency ranking, not a replacement, and it is
# not a full background-rate model. dNdScv and MutSigCV model the expected
# mutation rate per gene from length and sequence context and would rank both
# classes correctly; neither is in the container, and adding them means a
# rebuild. Worth doing -- and worth saying that this is the cheaper stand-in.
# --------------------------------------------------------------------------- #
if (isTRUE(args$oncodrive) && !is.null(maf_dt_unfiltered)) {
  # On the UNFILTERED variants -- see the note at the copy() above. This is the
  # independent ranking, not a re-ranking of what the blocklist allowed.
  mohq_log(sprintf("oncodrive: testing %s genes on the UNFILTERED variant set (%s variants)",
                   format(uniqueN(maf_dt_unfiltered$Hugo_Symbol), big.mark = ","),
                   format(nrow(maf_dt_unfiltered), big.mark = ",")))
  maf_obj_unfiltered <- tryCatch(
    read.maf(maf = as.data.frame(maf_dt_unfiltered), verbose = FALSE),
    error = function(e) { mohq_warn("read.maf on unfiltered set failed: ", e$message); NULL })

  # NAME THE AMINO-ACID COLUMN EXPLICITLY.
  #
  # oncodrive scores positional clustering, so it needs the protein change per
  # variant. Left to auto-detect it looks for names like 'AAChange' and dies
  # with "argument is of length zero" when it finds none -- vcf2maf writes
  # HGVSp_Short. Resolve it here and fail with a readable message instead.
  aa_col <- intersect(c("HGVSp_Short", "AAChange", "Protein_Change",
                        "amino_acid_change", "HGVSp"),
                      names(maf_dt_unfiltered))[1]

  # BACKGROUND: predefined, not cohort-derived, and that is a real limitation.
  #
  # oncodrive normally estimates the background clustering score from SYNONYMOUS
  # mutations -- they are under no positional selection, so they show what
  # clustering looks like by chance in this data. This pipeline loads
  # non-synonymous variants only (51,706 of 19,486,243 here), so there are none
  # and maftools falls back to constants (mean 0.279, SD 0.13) derived from
  # other cohorts.
  #
  # The ranking is still informative -- it is the same statistic applied to
  # every gene -- but the p-values are calibrated against someone else's data,
  # not yours. Treat the ORDER as the result and the significance as indicative.
  # Loading synonymous variants too would fix it and would mean reading ~19.5M
  # variants per cohort instead of ~52k.
  has_syn <- FALSE
  # is.na, not is.null: intersect(...)[1] on an empty intersection returns
  # NA_character_, so the is.null() guard passed and AACol = NA went into
  # maftools, which fell back to its default 'AAChange' and died looking for a
  # column nobody had. A missing value that is not NULL is the classic R trap.
  if (length(aa_col) == 0 || is.na(aa_col)) {
    mohq_warn("No amino-acid change column found (looked for HGVSp_Short, ",
              "AAChange, Protein_Change). Skipping the clustering test; ",
              "available columns: ", paste(head(names(maf_dt_unfiltered), 25),
                                           collapse = ", "))
    od <- NULL
  } else {
    mohq_log("oncodrive: using '", aa_col, "' for protein position; ",
             "background = predefined constants (no synonymous variants loaded)")
    od <- if (is.null(maf_obj_unfiltered)) NULL else tryCatch(
      oncodrive(maf = maf_obj_unfiltered, AACol = aa_col,
                minMut = args$oncodrive_min_mut, pvalMethod = "zscore",
                bgEstimate = has_syn),
      error = function(e) { mohq_warn("oncodrive failed: ", e$message); NULL })
  }

  if (!is.null(od) && nrow(od)) {
    od <- as.data.table(od)
    fwrite(od, paste0(args$out_prefix, "_oncodrive.tsv"), sep = "\t")

    sig <- od[!is.na(fdr) & fdr < 0.1]
    mohq_log(sprintf("oncodrive: %d gene(s) tested, %d with FDR < 0.1%s",
                     nrow(od), nrow(sig),
                     if (nrow(sig)) paste0(": ", paste(head(sig[order(fdr)]$Hugo_Symbol, 8),
                                                       collapse = ", ")) else ""))

    # Compare the two rankings explicitly. A gene high on frequency and absent
    # from the clustering result is the length-bias case; a gene high on
    # clustering and modest on frequency is the one frequency was hiding.
    # Both rankings come from rank_all, which is the UNFILTERED frequency table
    # written before drop_flags() -- so the two columns are computed on exactly
    # the same variants and the comparison means something. is_flags is carried
    # through as a third opinion: it shows what the blocklist would have done.
    cmp <- merge(rank_all[, .(Hugo_Symbol, n_patients, n_variants, is_flags)],
                 od[, .(Hugo_Symbol, clusters, fdr)],
                 by = "Hugo_Symbol", all.x = TRUE)
    setorder(cmp, -n_patients)
    cmp[, rank_frequency := seq_len(.N)]
    setorder(cmp, fdr, na.last = TRUE)
    cmp[, rank_clustering := fifelse(is.na(fdr), NA_integer_, seq_len(.N))]
    setorder(cmp, rank_frequency)
    fwrite(cmp, paste0(args$out_prefix, "_rank_comparison.tsv"), sep = "\t")

    top_freq <- head(cmp[order(rank_frequency)]$Hugo_Symbol, 5)
    top_clus <- head(cmp[!is.na(rank_clustering)][order(rank_clustering)]$Hugo_Symbol, 5)
    mohq_log("top 5 by frequency : ", paste(top_freq, collapse = ", "))
    mohq_log("top 5 by clustering: ", paste(top_clus, collapse = ", "))

    if (nrow(sig)) {
      png(paste0(args$out_prefix, "_oncodrive.png"),
          width = 8, height = 6, units = "in", res = 300)
      tryCatch(plotOncodrive(res = od, fdrCutOff = 0.1, useFraction = TRUE),
               error = function(e) mohq_warn("plotOncodrive failed: ", e$message))
      dev.off()
    } else {
      mohq_warn("oncodrive found no gene with FDR < 0.1. With a small cohort ",
                "that is expected: clustering needs enough mutations per gene ",
                "to be detectable. The table is still written.")
    }
  }
}

# Oncogenic pathway involvement (Sanchez-Vega curated pathways).
# maftools renamed OncogenicPathways() -> pathways() around v2.12, so try both
# rather than depending on one version.
# CAPTURE THE PLOT, not just the table.
#
# Both pathways() and OncogenicPathways() DRAW as a side effect and return the
# table. With no device open the drawing went to a default Rplots.pdf in the
# work directory -- so `_pathways.tsv` was written, `_pathways.png` never was,
# and the report showed an explanatory note above a figure that did not exist.
# The `pdf 2` line in the run log was that stray device closing.
pw <- tryCatch({
  png(paste0(args$out_prefix, "_pathways.png"), width = 9, height = 7,
      units = "in", res = 300)
  res <- tryCatch(
    if (exists("pathways")) pathways(maf = maf_obj, plotType = "treemap")
    else                    OncogenicPathways(maf = maf_obj),
    finally = dev.off())
  fwrite(as.data.table(res), paste0(args$out_prefix, "_pathways.tsv"), sep = "\t")
  res
}, error = function(e) {
  mohq_warn("Oncogenic pathway analysis unavailable: ", e$message); NULL
})

# An empty PNG is worse than none: the report would show a blank panel with no
# indication anything went wrong. Remove it so absent() reports it honestly.
pw_png <- paste0(args$out_prefix, "_pathways.png")
if (file.exists(pw_png) && file.info(pw_png)$size < 1000) {
  unlink(pw_png)
  mohq_warn("Pathway plot came out empty and was removed.")
}

# Per-patient variant classification summary, useful for QC drift over years.
png(paste0(args$out_prefix, "_summary.png"), width = 12, height = 7,
    units = "in", res = 300)
plotmafSummary(maf = maf_obj, addStat = "median", dashboard = TRUE)
dev.off()

mohq_session_info(paste0(args$out_prefix, ".versions.txt"))

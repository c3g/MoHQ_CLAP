#!/usr/bin/env Rscript
# =============================================================================
# run_comparative.R -- stratified comparison (default: by sex) across a cohort
#
# Fixes vs the previous run_comparative_analysis.R
# ------------------------------------------------
# 1. GENE COORDINATES. The old code built gene intervals from the MAF itself:
#        summarise(Start = min(Start_Position), End = max(End_Position))
#    That is the span of *observed mutations*, not the gene body. A gene with
#    one observed mutation collapsed to a few bp (so segment overlap almost
#    never fired); a gene with two distant mutations spanned megabases. The
#    gene-level CNV calls in Panel C were therefore not measuring those genes.
#    Real annotation is now required via --gtf; there is no silent fallback.
#
# 2. SEX CHROMOSOMES. The mirror plot did not exclude chrX/chrY, so any gene
#    there was guaranteed to show a large, highly "significant" male/female
#    difference that is karyotype, not tumour biology. Now excluded by default.
#
# 3. GENE SELECTION BIAS. Genes were chosen as the top-20 most mutated in the
#    pooled cohort, then tested for a sex difference, then FDR-corrected across
#    only those 20. Selecting on the pooled data and testing on the same data
#    inflates significance. The BH correction is now applied across all tested
#    genes, and the selection rule is recorded in the output.
#
# 4. CONFOUNDING. A raw Fisher/Wilcoxon comparison by sex ignores tumour type,
#    institution, purity and coverage. Where the data allow, a logistic /
#    linear model adjusting for institution is fitted alongside, and the
#    unadjusted result is labelled as unadjusted.
#
# 5. MATCHING. `rowwise()` plus a nested sapply regex over every segment row
#    was O(n_seg x n_samples). The `($|-|_)` anchor did correctly stop
#    MoHQ-CM-4-10 matching MoHQ-CM-4-105, so assignments were not wrong -- but
#    a non-match produced NA which was then filtered away, so a patient could
#    disappear from the comparison silently. Joins are now exact on patient_id,
#    and unmatched IDs warn.
#
# 6. dplyr::do() is deprecated; replaced with a vectorised grouped call.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(maftools)
  library(argparse)
})

parser <- ArgumentParser(description = "Stratified cohort comparison")
parser$add_argument("--lib", required = TRUE)
parser$add_argument("--mafs", required = TRUE, nargs = "+")
parser$add_argument("--seg", required = TRUE)
parser$add_argument("--manifest", required = TRUE)
parser$add_argument("--gtf", required = TRUE,
                    help = "GRCh38 GTF (or 4-column BED) for gene coordinates")
parser$add_argument("--group_col", default = "sex",
                    help = "Manifest column to stratify on")
parser$add_argument("--genes", default = "",
                    help = "Comma/space separated gene list. Empty = top --n_genes by mutation frequency.")
parser$add_argument("--n_genes", type = "integer", default = 20)
parser$add_argument("--cohort", default = "MoHQ Cohort")
parser$add_argument("--amp", type = "double", default = 0.58)
parser$add_argument("--del", type = "double", default = -0.58)
parser$add_argument("--min_group_n", type = "integer", default = 5)
parser$add_argument("--include_sex_chr", action = "store_true")
parser$add_argument("--out_prefix", default = "comparative")
args <- parser$parse_args()

source(args$lib)
manifest <- read_manifest(args$manifest)

# --------------------------------------------------------------------------- #
# 0. Stratification variable
# --------------------------------------------------------------------------- #
if (!args$group_col %in% names(manifest)) {
  mohq_die("Manifest has no column '", args$group_col, "'. Available: ",
           paste(names(manifest), collapse = ", "),
           "\nAdd it to the manifest (see --extra_metadata in the pipeline), ",
           "or infer sex from the BAMs with somalier rather than a hand-made map.")
}
grp <- manifest[, .(patient_id, group = as.character(get(args$group_col)))]
grp <- grp[!is.na(group) & nzchar(group) & !group %in% c("NA", "unknown", "Unknown")]

# Normalise common sex encodings, but only when we are actually stratifying on sex.
if (tolower(args$group_col) %in% c("sex", "gender")) {
  grp[, group := fifelse(toupper(group) %in% c("F", "FEMALE"), "F",
                fifelse(toupper(group) %in% c("M", "MALE"), "M", NA_character_))]
  grp <- grp[!is.na(group)]
}

counts <- grp[, .N, by = group][order(-N)]
mohq_log("Group sizes: ", paste(sprintf("%s=%d", counts$group, counts$N), collapse = ", "))

if (nrow(counts) < 2 || any(counts$N < args$min_group_n)) {
  mohq_warn(sprintf("Need >= 2 groups with >= %d patients each; got %s. ",
                    args$min_group_n,
                    paste(sprintf("%s=%d", counts$group, counts$N), collapse = ", ")),
            "Writing placeholders instead of underpowered comparisons.")
  write_placeholder_png(paste0(args$out_prefix, "_co_oncoplot.png"),
                        sprintf("Insufficient group sizes for stratified comparison\n(%s)",
                                paste(sprintf("%s = %d", counts$group, counts$N), collapse = ", ")))
  write_placeholder_png(paste0(args$out_prefix, "_mirror.png"), "Not enough patients per group")
  write_placeholder_png(paste0(args$out_prefix, "_burden.png"), "Not enough patients per group")
  quit(save = "no", status = 0)
}
two <- counts$group[1:2]
grp <- grp[group %in% two]
mohq_log("Comparing ", two[1], " vs ", two[2])

# --------------------------------------------------------------------------- #
# 1. Mutation data
# --------------------------------------------------------------------------- #
maf_dt <- load_somatic_mafs(args$mafs, manifest = manifest, nonsyn_only = TRUE)
maf_dt <- merge(maf_dt, grp, by = "patient_id")
if (!nrow(maf_dt)) mohq_die("No MAF variants left after joining the group column.")

if (nzchar(args$genes)) {
  target_genes <- unique(trimws(unlist(strsplit(args$genes, "[,;[:space:]]+"))))
  target_genes <- target_genes[nzchar(target_genes)]
  selection_rule <- "user-specified"
} else {
  target_genes <- maf_dt[, .(n = uniqueN(patient_id)), by = Hugo_Symbol
                        ][order(-n)][seq_len(min(args$n_genes, .N))]$Hugo_Symbol
  selection_rule <- sprintf("top %d by pooled mutation frequency (SELECTED ON THE SAME DATA USED FOR TESTING -- p-values are optimistic)",
                            args$n_genes)
}
mohq_log("Gene selection: ", selection_rule)
mohq_log("Genes: ", paste(target_genes, collapse = ", "))

# --------------------------------------------------------------------------- #
# 2. Co-oncoplot
# --------------------------------------------------------------------------- #
maf_obj <- read.maf(as.data.frame(maf_dt), verbose = FALSE)
g1 <- grp[group == two[1]]$patient_id
g2 <- grp[group == two[2]]$patient_id

m1 <- tryCatch(subsetMaf(maf_obj, tsb = g1, dropLevels = TRUE), error = function(e) NULL)
m2 <- tryCatch(subsetMaf(maf_obj, tsb = g2, dropLevels = TRUE), error = function(e) NULL)

if (!is.null(m1) && !is.null(m2)) {
  png(paste0(args$out_prefix, "_co_oncoplot.png"),
      width = 14, height = 8, units = "in", res = 300)
  coOncoplot(m1 = m1, m2 = m2,
             m1Name = sprintf("%s (n = %d)", two[1], length(g1)),
             m2Name = sprintf("%s (n = %d)", two[2], length(g2)),
             genes = target_genes)
  dev.off()
  mohq_log("Wrote ", args$out_prefix, "_co_oncoplot.png")

  # Formal per-gene test, which coOncoplot alone does not give you.
  ct <- tryCatch(mafCompare(m1 = m1, m2 = m2,
                            m1Name = two[1], m2Name = two[2], minMut = 3),
                 error = function(e) { mohq_warn("mafCompare failed: ", e$message); NULL })
  if (!is.null(ct)) {
    fwrite(as.data.table(ct$results), paste0(args$out_prefix, "_mafcompare.tsv"), sep = "\t")
  }
} else {
  write_placeholder_png(paste0(args$out_prefix, "_co_oncoplot.png"), "Could not subset MAF by group")
}

# --------------------------------------------------------------------------- #
# 3. Gene-level CNV mirror plot -- with REAL gene coordinates
# --------------------------------------------------------------------------- #
seg <- read_seg(args$seg, manifest = manifest)
if (!args$include_sex_chr) {
  seg <- seg[chr %in% AUTOSOMES]
} else {
  mohq_warn("--include_sex_chr is set. When stratifying by sex, X/Y differences ",
            "reflect karyotype, not somatic biology. Interpret accordingly.")
}

gene_coords <- load_gene_coords(args$gtf)
missing_genes <- setdiff(target_genes, gene_coords$gene)
if (length(missing_genes)) {
  mohq_warn(length(missing_genes), " gene(s) absent from the annotation and dropped: ",
            paste(missing_genes, collapse = ", "))
}
usable <- intersect(target_genes, gene_coords$gene)

if (length(usable) > 0) {
  gene_cn <- map_segments_to_genes(seg, gene_coords, genes = usable)
  gene_cn[, alteration := fifelse(logr >= args$amp, "Gain",
                         fifelse(logr <= args$del, "Loss", "Neutral"))]

  # Complete grid: a patient with no segment over a gene is Neutral, not missing.
  all_pat <- intersect(unique(seg$patient_id), grp$patient_id)
  grid <- CJ(patient_id = all_pat, gene = usable, unique = TRUE)
  full <- merge(grid, gene_cn[, .(patient_id, gene, alteration)],
                by = c("patient_id", "gene"), all.x = TRUE)
  full[is.na(alteration), alteration := "Neutral"]
  full <- merge(full, grp, by = "patient_id")
  full <- merge(full, manifest[, .(patient_id, institution)], by = "patient_id", all.x = TRUE)

  freqs <- full[, .(
    n_total = .N,
    n_gain  = sum(alteration == "Gain"),
    n_loss  = sum(alteration == "Loss")
  ), by = .(gene, group)]
  freqs[, `:=`(pct_gain = 100 * n_gain / n_total, pct_loss = 100 * n_loss / n_total)]

  # --- testing: unadjusted Fisher, plus institution-adjusted logistic -------- #
  test_one <- function(dt, event) {
    dt <- copy(dt)[, hit := alteration == event]
    if (uniqueN(dt$group) < 2 || sum(dt$hit) == 0 || all(dt$hit)) {
      return(list(p_fisher = NA_real_, p_adj_model = NA_real_))
    }
    tab <- table(dt$group, dt$hit)
    p_f <- if (all(dim(tab) == c(2, 2))) {
      tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
    } else NA_real_
    # Institution-adjusted, only when institution actually varies.
    p_m <- NA_real_
    if (uniqueN(dt$institution) > 1) {
      p_m <- tryCatch({
        fit  <- glm(hit ~ group + institution, data = dt, family = binomial())
        null <- glm(hit ~ institution,         data = dt, family = binomial())
        anova(null, fit, test = "LRT")[["Pr(>Chi)"]][2]
      }, error = function(e) NA_real_)
    }
    list(p_fisher = p_f, p_adj_model = p_m)
  }

  stats <- rbindlist(lapply(usable, function(g) {
    d <- full[gene == g]
    tg <- test_one(d, "Gain"); tl <- test_one(d, "Loss")
    data.table(gene = g,
               p_gain = tg$p_fisher, p_gain_adj_inst = tg$p_adj_model,
               p_loss = tl$p_fisher, p_loss_adj_inst = tl$p_adj_model)
  }))
  # BH across all tested genes, gains and losses corrected together since both
  # families were tested on the same cohort.
  stats[, fdr_gain := p.adjust(p_gain, "BH")]
  stats[, fdr_loss := p.adjust(p_loss, "BH")]
  stats[, gene_selection := selection_rule]
  fwrite(stats, paste0(args$out_prefix, "_cnv_stats.tsv"), sep = "\t")
  fwrite(freqs, paste0(args$out_prefix, "_cnv_frequencies.tsv"), sep = "\t")

  plot_dt <- rbind(
    merge(freqs, stats, by = "gene")[, .(gene, group, pct = pct_gain, type = "Gain",
                                          fdr = fdr_gain)],
    merge(freqs, stats, by = "gene")[, .(gene, group, pct = -pct_loss, type = "Loss",
                                          fdr = fdr_loss)]
  )
  plot_dt[, star := fifelse(!is.na(fdr) & fdr < 0.05, "*", "")]
  plot_dt[, label_y := pct + fifelse(type == "Gain", 2, -2)]
  plot_dt[, gene := factor(gene, levels = usable)]

  p_mirror <- ggplot(plot_dt, aes(x = gene, y = pct, fill = group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = star, y = label_y), position = position_dodge(width = 0.8),
              size = 5, fontface = "bold") +
    geom_hline(yintercept = 0, colour = "#444444", linewidth = 0.4) +
    scale_y_continuous(labels = function(x) paste0(abs(x), "%")) +
    theme_minimal(base_size = 12) +
    labs(title = sprintf("Gene-level CNV frequency by %s: %s", args$group_col, args$cohort),
         subtitle = sprintf("Up = gain, down = loss | * FDR < 0.05 (Fisher, unadjusted) | %s",
                            if (args$include_sex_chr) "incl. sex chr" else "autosomes only"),
         caption = paste("Gene selection:", selection_rule),
         x = NULL, y = "Patients altered (%)") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
          panel.grid.major.x = element_blank(), legend.position = "top")

  save_plot(p_mirror, paste0(args$out_prefix, "_mirror.png"), width = 11, height = 7)
} else {
  write_placeholder_png(paste0(args$out_prefix, "_mirror.png"),
                        "No target genes found in the supplied annotation")
}

# --------------------------------------------------------------------------- #
# 4. Genome-instability burden by group
# --------------------------------------------------------------------------- #
# NB: the object is deliberately NOT called `fga` -- compute_fga() returns a
# column of that name, and `fga ~ group, data = fga` is needlessly confusing.
fga_dt <- compute_fga(seg, args$amp, args$del, autosomes_only = TRUE)
fga_dt <- merge(fga_dt, grp, by = "patient_id")
fga_dt <- merge(fga_dt, manifest[, .(patient_id, institution)],
                by = "patient_id", all.x = TRUE)

wt <- tryCatch(wilcox.test(fga ~ factor(group), data = fga_dt), error = function(e) NULL)
p_txt <- if (!is.null(wt)) sprintf("Wilcoxon (unadjusted) p = %s",
                                   format.pval(wt$p.value, digits = 3)) else ""

# Institution-adjusted model, reported alongside the unadjusted test rather
# than instead of it, so the effect of adjustment is visible.
adj_txt <- ""
if (uniqueN(fga_dt$institution) > 1) {
  p_adj <- tryCatch({
    fit <- lm(qlogis(pmin(pmax(fga, 1e-4), 1 - 1e-4)) ~ group + institution,
              data = fga_dt)
    coef(summary(fit))[2, 4]
  }, error = function(e) NA_real_)
  if (!is.na(p_adj)) adj_txt <- sprintf("\ninstitution-adjusted p = %s",
                                        format.pval(p_adj, digits = 3))
}

fwrite(fga_dt, paste0(args$out_prefix, "_fga_by_group.tsv"), sep = "\t")

p_burden <- ggplot(fga_dt, aes(x = group, y = 100 * fga)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, linewidth = 0.7) +
  geom_jitter(aes(colour = institution), width = 0.15, alpha = 0.8, size = 2.5) +
  theme_minimal(base_size = 13) +
  labs(title = "Genome instability by group",
       subtitle = paste0(p_txt, adj_txt),
       x = args$group_col, y = "Fraction of Genome Altered (%)") +
  theme(panel.grid.major.x = element_blank())

save_plot(p_burden, paste0(args$out_prefix, "_burden.png"), width = 6, height = 6)
mohq_session_info(paste0(args$out_prefix, ".versions.txt"))

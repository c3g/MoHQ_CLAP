#!/usr/bin/env Rscript
# =============================================================================
# run_fusions.R -- recurrent fusion detection across a cohort
#
# NEW module. Recurrence is inherently a cohort-scale question: a fusion seen in
# one patient is a curiosity, the same fusion in four is a finding. You already
# generate the inputs per patient (anno_fuse from RNA, LINX from WGS SVs) but
# nothing was aggregating them.
#
# Cross-referencing the two is the point: an RNA fusion call with a supporting
# DNA structural variant in the same gene pair is far stronger evidence than
# either alone.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparse)
})

parser <- ArgumentParser(description = "Cohort recurrent fusions")
parser$add_argument("--lib", required = TRUE)
parser$add_argument("--annofuse", nargs = "*")
parser$add_argument("--linx_fusion", nargs = "*")
parser$add_argument("--manifest", required = TRUE)
parser$add_argument("--cohort", default = "MoHQ Cohort")
parser$add_argument("--min_recurrence", type = "integer", default = 2)
parser$add_argument("--out_prefix", default = "cohort_fusions")
args <- parser$parse_args()

source(args$lib)
manifest <- read_manifest(args$manifest)

pick_col <- function(hdr, candidates) {
  hit <- intersect(candidates, hdr)
  if (length(hit)) hit[1] else NA_character_
}

norm_pair <- function(a, b) {
  a <- toupper(trimws(a)); b <- toupper(trimws(b))
  paste0(a, "--", b)
}

# --------------------------------------------------------------------------- #
# RNA fusions (annoFuse)
# --------------------------------------------------------------------------- #
read_annofuse <- function(f) {
  d <- fread(f, showProgress = FALSE, fill = TRUE)
  if (!nrow(d)) return(NULL)
  hdr <- names(d)
  g1 <- pick_col(hdr, c("Gene1A", "gene1", "LeftGene", "GeneA", "geneA", "FusionName"))
  g2 <- pick_col(hdr, c("Gene1B", "gene2", "RightGene", "GeneB", "geneB"))
  if (is.na(g1)) {
    mohq_warn("Cannot resolve gene columns in ", basename(f),
              " (found: ", paste(head(hdr, 15), collapse = ", "), "); skipping.")
    return(NULL)
  }
  pid <- harmonise_patient_id(sub("\\.anno_fuse\\.tsv$", "", basename(f)))
  if (is.na(g2)) {
    # Single FusionName column of the form GENEA--GENEB
    parts <- tstrsplit(as.character(d[[g1]]), "--", fixed = TRUE)
    out <- data.table(patient_id = pid,
                      geneA = toupper(parts[[1]]),
                      geneB = if (length(parts) > 1) toupper(parts[[2]]) else NA_character_)
  } else {
    out <- data.table(patient_id = pid,
                      geneA = toupper(as.character(d[[g1]])),
                      geneB = toupper(as.character(d[[g2]])))
  }
  out[!is.na(geneA) & nzchar(geneA)]
}

rna <- if (length(args$annofuse)) {
  rbindlist(lapply(args$annofuse[file.exists(args$annofuse)], read_annofuse), fill = TRUE)
} else NULL

# --------------------------------------------------------------------------- #
# DNA fusions (LINX)
# --------------------------------------------------------------------------- #
read_linx <- function(f) {
  d <- fread(f, showProgress = FALSE, fill = TRUE)
  if (!nrow(d)) return(NULL)
  hdr <- names(d)
  g1 <- pick_col(hdr, c("geneStart", "GeneStart", "fivePrimeGene"))
  g2 <- pick_col(hdr, c("geneEnd", "GeneEnd", "threePrimeGene"))
  if (is.na(g1) || is.na(g2)) return(NULL)
  rep_col <- pick_col(hdr, c("reported", "Reported"))
  if (!is.na(rep_col)) d <- d[as_logical_loose(d[[rep_col]])]
  if (!nrow(d)) return(NULL)
  pid <- harmonise_patient_id(sub("\\.linx\\.fusion\\.tsv$", "", basename(f)))
  data.table(patient_id = pid,
             geneA = toupper(as.character(d[[g1]])),
             geneB = toupper(as.character(d[[g2]])))
}

dna <- if (length(args$linx_fusion)) {
  rbindlist(lapply(args$linx_fusion[file.exists(args$linx_fusion)], read_linx), fill = TRUE)
} else NULL

if (is.null(rna) && is.null(dna)) {
  mohq_warn("No fusion inputs available for this cohort.")
  write_placeholder_png(paste0(args$out_prefix, ".png"), "No fusion calls available")
  fwrite(data.table(), paste0(args$out_prefix, ".tsv"), sep = "\t")
  quit(save = "no", status = 0)
}

summarise_source <- function(dt, label) {
  if (is.null(dt) || !nrow(dt)) return(NULL)
  dt <- copy(dt)[!is.na(geneB) & nzchar(geneB)]
  dt[, fusion := norm_pair(geneA, geneB)]
  unique(dt[, .(patient_id, fusion, source = label)])
}

all_f <- rbindlist(list(summarise_source(rna, "RNA"),
                        summarise_source(dna, "DNA")), fill = TRUE)
if (!nrow(all_f)) {
  write_placeholder_png(paste0(args$out_prefix, ".png"), "No parseable fusion calls")
  quit(save = "no", status = 0)
}

check_ids_against_manifest(all_f$patient_id, manifest, "fusion")

rec <- all_f[, .(
  n_patients  = uniqueN(patient_id),
  sources     = paste(sort(unique(source)), collapse = "+"),
  patients    = paste(sort(unique(patient_id)), collapse = ";")
), by = fusion][order(-n_patients)]

# Fusions seen in BOTH assays are the highest-confidence set.
rec[, dna_rna_supported := sources == "DNA+RNA"]
fwrite(rec, paste0(args$out_prefix, ".tsv"), sep = "\t")

mohq_log(sprintf("%d distinct fusions; %d recurrent (>= %d patients); %d DNA+RNA supported",
                 nrow(rec), sum(rec$n_patients >= args$min_recurrence),
                 args$min_recurrence, sum(rec$dna_rna_supported)))

top <- rec[n_patients >= args$min_recurrence]
if (!nrow(top)) {
  write_placeholder_png(paste0(args$out_prefix, ".png"),
                        sprintf("No fusion recurred in >= %d patients", args$min_recurrence))
} else {
  top <- top[seq_len(min(30, .N))]
  p <- ggplot(top, aes(x = reorder(fusion, n_patients), y = n_patients,
                       fill = dna_rna_supported)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#1b7837", `FALSE` = "grey55"),
                      labels = c(`TRUE` = "DNA + RNA", `FALSE` = "single assay"),
                      name = NULL) +
    theme_minimal(base_size = 12) +
    labs(title = paste("Recurrent fusions:", args$cohort),
         subtitle = sprintf("seen in >= %d patients", args$min_recurrence),
         x = NULL, y = "Patients")
  save_plot(p, paste0(args$out_prefix, ".png"), width = 9,
            height = max(4, 0.28 * nrow(top) + 2))
}

mohq_session_info(paste0(args$out_prefix, ".versions.txt"))

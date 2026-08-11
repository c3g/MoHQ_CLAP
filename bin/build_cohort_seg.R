#!/usr/bin/env Rscript
# =============================================================================
# build_cohort_seg.R -- assemble a cohort-level .seg from per-patient CNVkit VCFs
#
# Why this exists: previously the cohort .seg was produced by hand, outside
# Nextflow, and dropped into `../results/all_cohorts_analysis_ready/...`. That
# breaks provenance (nothing records how it was made), breaks `-resume`, and
# makes the cohort silently stale whenever a patient is reprocessed. Building it
# inside the pipeline from the manifest fixes all three.
#
# Sample IDs are canonicalised to patient_id here, once, so every downstream
# join is exact.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(argparse)
})

parser <- ArgumentParser(description = "Build a cohort .seg from CNVkit VCFs")
parser$add_argument("--lib", required = TRUE, help = "Path to mohq_common.R")
parser$add_argument("--vcfs", required = TRUE, nargs = "+", help = "CNVkit VCF(s)")
parser$add_argument("--manifest", required = TRUE)
parser$add_argument("--cohort", required = TRUE)
parser$add_argument("--out", default = NULL)
parser$add_argument("--min_markers", type = "integer", default = 0,
                    help = "Drop segments supported by fewer probes than this")
args <- parser$parse_args()

source(args$lib)                      # defines %||% among other helpers
manifest <- read_manifest(args$manifest)
outfile <- args$out %||% paste0(args$cohort, ".seg")

# INFO keys that CNVkit / other callers use for the log2 ratio, in priority order.
LOGR_KEYS <- c("FOLD_CHANGE_LOG", "LOG2", "SEGMEAN", "LOG2RATIO")

parse_info <- function(info, key) {
  pat <- paste0("(^|;)", key, "=([^;]+)")
  m <- regmatches(info, regexpr(pat, info))
  out <- rep(NA_character_, length(info))
  hit <- nzchar(m)
  out[which(hit)] <- sub(pat, "\\2", m[hit])
  out
}

read_cnvkit_vcf <- function(path) {
  # data.table reads .gz natively; skip the ## header block.
  v <- fread(path, sep = "\t", header = FALSE, skip = "#CHROM", quote = "",
             showProgress = FALSE)
  if (!nrow(v)) return(NULL)
  setnames(v, 1:8, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))

  end <- suppressWarnings(as.numeric(parse_info(v$INFO, "END")))
  if (all(is.na(end))) {
    mohq_die("No END= in INFO for ", basename(path),
             " -- this does not look like a CNVkit segment VCF.")
  }

  logr <- NULL
  used_key <- NA_character_
  for (k in LOGR_KEYS) {
    cand <- suppressWarnings(as.numeric(parse_info(v$INFO, k)))
    if (!all(is.na(cand))) { logr <- cand; used_key <- k; break }
  }
  if (is.null(logr)) {
    mohq_die("Could not find a log2-ratio INFO key in ", basename(path),
             ". Tried: ", paste(LOGR_KEYS, collapse = ", "),
             "\n  First INFO field: ", substr(v$INFO[1], 1, 200),
             "\n  Add the correct key to LOGR_KEYS in build_cohort_seg.R.")
  }
  probes <- suppressWarnings(as.numeric(parse_info(v$INFO, "PROBES")))

  # Patient ID comes from the filename: <patient>.cnvkit.vcf.gz
  pid <- harmonise_patient_id(sub("\\.cnvkit\\.vcf\\.gz$", "", basename(path)))

  data.table(
    sample     = pid,
    patient_id = pid,
    chr        = normalise_chr(v$CHROM),
    start      = as.numeric(v$POS),
    end        = end,
    markers    = probes,
    logr       = logr,
    logr_key   = used_key
  )
}

mohq_log("Reading ", length(args$vcfs), " CNVkit VCF(s) for cohort ", args$cohort)
segs <- rbindlist(lapply(args$vcfs, read_cnvkit_vcf), fill = TRUE)
if (!nrow(segs)) mohq_die("No segments parsed from any VCF.")

keys <- unique(segs$logr_key)
if (length(keys) > 1) {
  mohq_warn("Different log2 INFO keys across patients (", paste(keys, collapse = ", "),
            ") -- values may not be on the same scale. Investigate before trusting FGA.")
}
segs[, logr_key := NULL]

n0 <- nrow(segs)
segs <- segs[!is.na(logr) & !is.na(start) & !is.na(end) & end > start]
if (args$min_markers > 0 && !all(is.na(segs$markers))) {
  segs <- segs[is.na(markers) | markers >= args$min_markers]
}
mohq_log(sprintf("Retained %s/%s segments", format(nrow(segs), big.mark = ","),
                 format(n0, big.mark = ",")))

check_ids_against_manifest(segs$patient_id, manifest, "CNVkit")

# Standard IGV .seg column order and names.
out <- segs[, .(Sample = patient_id, Chromosome = chr, Start = start,
                End = end, Num_Probes = markers, Segment_Mean = logr)]
setorder(out, Sample, Chromosome, Start)
fwrite(out, outfile, sep = "\t", na = "NA", quote = FALSE)
mohq_log("Wrote ", outfile, " (", uniqueN(out$Sample), " patients, ",
         format(nrow(out), big.mark = ","), " segments)")

# GISTIC marker file, derived from the same segments so the two always agree.
markers <- unique(rbindlist(list(
  out[, .(Chromosome, Position = Start)],
  out[, .(Chromosome, Position = End)]
)))
setorder(markers, Chromosome, Position)
markers[, Marker_Name := paste0("mk_", Chromosome, "_", Position)]
fwrite(markers[, .(Marker_Name, Chromosome, Position)],
       paste0(args$cohort, ".markers.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
mohq_log("Wrote ", args$cohort, ".markers.txt (", nrow(markers), " markers)")

mohq_session_info(paste0(args$cohort, ".build_seg.versions.txt"))

#!/usr/bin/env Rscript
# =============================================================================
# mohq_common.R -- shared utilities for MoHQ cohort-scale analysis
#
# Every analysis script sources this file (staged by Nextflow and passed via
# --lib). Putting the shared logic here means the FGA definition, the segment
# threshold, and the sample-ID rules exist exactly once. Previously each script
# reimplemented them and they had drifted apart -- e.g. the genome denominator
# was 2875 Mb in one script and 2744 Mb in another, producing two different
# FGA values for the same sample in the same report.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# --------------------------------------------------------------------------- #
# Logging
# --------------------------------------------------------------------------- #

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

mohq_log <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))
}

mohq_warn <- function(...) {
  message(sprintf("[%s] WARNING: %s", format(Sys.time(), "%H:%M:%S"),
                  paste0(..., collapse = "")))
}

mohq_die <- function(...) {
  stop(paste0(..., collapse = ""), call. = FALSE)
}

# --------------------------------------------------------------------------- #
# Sample identity
# --------------------------------------------------------------------------- #

#' Reduce any MoHQ identifier to its canonical patient_id.
#'
#' This replaces the previous fuzzy prefix matching
#' (`grepl(paste0("^", s, "($|-|_)"), barcode)` inside a `rowwise()` loop).
#'
#' To be fair to that approach: the `($|-|_)` anchor DOES correctly prevent
#' MoHQ-CM-4-10 from matching MoHQ-CM-4-105, so it was not producing wrong
#' patient assignments. The problems were:
#'   * O(n_segments x n_samples) with rowwise() -- very slow at cohort scale
#'   * it took `[...][1]`, so a duplicated metadata row silently won at random
#'   * a non-match produced NA, which was then filtered out, so a patient could
#'     vanish from an analysis with no message at all
#'   * `s` was interpolated into a regex unescaped (harmless for these IDs, but
#'     it means the matcher is only accidentally correct)
#'
#' Reducing every identifier to a canonical form and joining exactly removes
#' all four, and makes an unmatched ID a loud warning rather than a silent drop.
#'
#' Handles every identifier form observed in the delivery tree:
#'   MoHQ-CM-4-10-157954-1DT   (full sample ID, numeric sample_num)
#'   MoHQ-CQ-34-01-RCC01n-1DN  (full sample ID, alphanumeric sample_num)
#'   MoHQ-CM-4-10_D            (PCGR barcode)
#'   MoHQ-CM-4-10              (already canonical)
harmonise_patient_id <- function(x) {
  x <- trimws(as.character(x))
  # Drop trailing file-ish decorations first.
  x <- sub("\\.(variants|sorted|bam|vcf|maf|tsv)$", "", x)
  # Full sample ID -> first four hyphen tokens.
  x <- sub("^([A-Za-z]+-[A-Za-z]+-[A-Za-z0-9]+-[A-Za-z0-9]+)-[A-Za-z0-9]+-[0-9]+(DN|DT|RN|RT)$",
           "\\1", x)
  # PCGR barcode suffixes.
  x <- sub("_(D|R|T|N|DNA|RNA)$", "", x)
  x
}

#' Assert that harmonisation produced IDs that exist in the manifest.
#' Fails loudly rather than silently dropping samples, which is what the old
#' `inner_join` on fuzzy matches did.
check_ids_against_manifest <- function(ids, manifest, what = "input") {
  ids <- unique(ids)
  known <- unique(manifest$patient_id)
  unmatched <- setdiff(ids, known)
  if (length(unmatched)) {
    mohq_warn(sprintf(
      "%d/%d %s IDs are not in the manifest and will be dropped: %s",
      length(unmatched), length(ids), what,
      paste(utils::head(unmatched, 10), collapse = ", ")))
  }
  if (length(unmatched) == length(ids)) {
    mohq_die("No ", what, " IDs matched the manifest at all. ",
             "This almost always means an ID-format change upstream. ",
             "Observed: ", paste(utils::head(ids, 3), collapse = ", "),
             " | expected e.g.: ", paste(utils::head(known, 3), collapse = ", "))
  }
  invisible(intersect(ids, known))
}

read_manifest <- function(path) {
  m <- fread(path, sep = "\t", header = TRUE, na.strings = c("NA", ""))
  req <- c("patient_id", "cohort_id", "institution")
  miss <- setdiff(req, names(m))
  if (length(miss)) mohq_die("Manifest is missing column(s): ", paste(miss, collapse = ", "))
  m[]
}

# --------------------------------------------------------------------------- #
# Segment (.seg) handling
# --------------------------------------------------------------------------- #

# Canonical internal names. The old code did
#   colnames(seg) <- c("sample","chr","start","end","markers","logr")
# which renames by POSITION and silently mislabels every column if the file has
# a different column order or an extra column. We resolve by name instead.
SEG_COL_PATTERNS <- list(
  sample  = "^(sample|sample_?id|tumor_?sample_?barcode|id)$",
  chr     = "^(chrom|chr|chromosome|seqnames)$",
  start   = "^(start|loc\\.start|segment_?start|start_?position)$",
  end     = "^(end|loc\\.end|segment_?end|end_?position)$",
  markers = "^(num_?probes|num\\.mark|markers|n_?probes|probes|num_?mark)$",
  logr    = "^(seg\\.mean|seg_?mean|log_?r|log2|log2_?ratio|mean_?log2_?ratio|segment_?mean)$"
)

#' Read a .seg file, resolving columns by name with an explicit failure mode.
read_seg <- function(path, manifest = NULL) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    mohq_die("Segment file missing or empty: ", path)
  }
  raw <- fread(path, sep = "\t", header = TRUE, na.strings = c("NA", "", "."))
  found <- names(raw)
  norm  <- tolower(gsub("[^A-Za-z0-9._]", "_", found))

  mapping <- vapply(SEG_COL_PATTERNS, function(pat) {
    hit <- which(grepl(pat, norm))
    if (length(hit) == 0) NA_character_ else found[hit[1]]
  }, character(1))

  hard_required <- c("sample", "chr", "start", "end", "logr")
  missing <- hard_required[is.na(mapping[hard_required])]
  if (length(missing)) {
    mohq_die(
      "Could not resolve required .seg column(s): ", paste(missing, collapse = ", "),
      "\n  File:    ", path,
      "\n  Columns: ", paste(found, collapse = ", "),
      "\n  Add the actual column name to SEG_COL_PATTERNS in mohq_common.R.")
  }

  seg <- data.table(
    sample  = as.character(raw[[mapping[["sample"]]]]),
    chr     = as.character(raw[[mapping[["chr"]]]]),
    start   = as.numeric(raw[[mapping[["start"]]]]),
    end     = as.numeric(raw[[mapping[["end"]]]]),
    markers = if (is.na(mapping[["markers"]])) NA_real_ else as.numeric(raw[[mapping[["markers"]]]]),
    logr    = as.numeric(raw[[mapping[["logr"]]]])
  )

  seg[, chr := normalise_chr(chr)]
  seg[, patient_id := harmonise_patient_id(sample)]

  n0 <- nrow(seg)
  seg <- seg[!is.na(start) & !is.na(end) & !is.na(logr) & end > start]
  if (nrow(seg) < n0) {
    mohq_warn(sprintf("Dropped %d of %d segments with missing coords/logr or end <= start.",
                      n0 - nrow(seg), n0))
  }
  if (!is.null(manifest)) check_ids_against_manifest(seg$patient_id, manifest, "segment")
  mohq_log(sprintf("Loaded %s segments for %d patients from %s",
                   format(nrow(seg), big.mark = ","),
                   uniqueN(seg$patient_id), basename(path)))
  seg[]
}

normalise_chr <- function(x) {
  x <- sub("^chr", "", as.character(x), ignore.case = TRUE)
  x <- sub("^23$", "X", x)
  x <- sub("^24$", "Y", x)
  x <- sub("^(MT|M)$", "MT", x)
  x
}

AUTOSOMES <- as.character(1:22)

#' Classify segments as Gain / Loss / Neutral using a single shared threshold.
classify_segments <- function(seg, amp_threshold, del_threshold) {
  stopifnot(amp_threshold > 0, del_threshold < 0)
  seg <- copy(seg)
  seg[, type := fifelse(logr >= amp_threshold, "Gain",
               fifelse(logr <= del_threshold, "Loss", "Neutral"))]
  seg[]
}

# --------------------------------------------------------------------------- #
# Fraction of Genome Altered
# --------------------------------------------------------------------------- #

#' Fraction of Genome Altered, per sample.
#'
#' FIX: the previous implementation divided altered length by a hardcoded
#' constant (2875 Mb in run_fga_burden.R, 2744 Mb in run_comparative_analysis.R).
#' Two problems with that:
#'   1. Two different constants gave two different FGA values for the same
#'      sample in the same report.
#'   2. A fixed denominator assumes every sample was profiled over the whole
#'      genome. Samples with sparser CNV calling then look artificially quiet,
#'      which correlates with coverage -- i.e. it manufactures a batch effect.
#'
#' The denominator here is the total length actually profiled FOR THAT SAMPLE
#' (the cBioPortal definition), so FGA is comparable across samples with
#' different callable footprints. `profiled_mb` is returned so you can see and
#' QC that footprint rather than having it silently absorbed into the metric.
compute_fga <- function(seg, amp_threshold, del_threshold, autosomes_only = TRUE,
                        min_profiled_mb = 1000) {
  d <- copy(as.data.table(seg))
  if (autosomes_only) d <- d[chr %in% AUTOSOMES]
  if (!nrow(d)) mohq_die("No segments left after autosome filtering -- check chromosome naming.")

  d[, len := end - start]
  d[, `:=`(is_gain = logr >= amp_threshold,
           is_loss = logr <= del_threshold)]

  res <- d[, .(
    n_segments  = .N,
    profiled_mb = sum(len) / 1e6,
    gain_mb     = sum(len[is_gain]) / 1e6,
    loss_mb     = sum(len[is_loss]) / 1e6
  ), by = .(patient_id)]

  res[, altered_mb := gain_mb + loss_mb]
  res[, `:=`(
    fga      = altered_mb / profiled_mb,
    fga_gain = gain_mb    / profiled_mb,
    fga_loss = loss_mb    / profiled_mb
  )]

  # QC: a sample profiled over far less than the autosomal genome is not
  # comparable to the rest of the cohort. Surface it, do not hide it.
  thin <- res[profiled_mb < min_profiled_mb]
  if (nrow(thin)) {
    mohq_warn(sprintf(
      "%d sample(s) profiled over < %.0f Mb of autosome and may not be comparable: %s",
      nrow(thin), min_profiled_mb,
      paste(sprintf("%s (%.0f Mb)", thin$patient_id, thin$profiled_mb), collapse = ", ")))
  }
  res[, low_coverage := profiled_mb < min_profiled_mb]
  res[order(-fga)]
}

# --------------------------------------------------------------------------- #
# Genome binning
# --------------------------------------------------------------------------- #

#' Expand each segment across EVERY bin it overlaps.
#'
#' FIX: the previous run_frequency.R computed
#'     mutate(bin = floor(start / bin_size) * bin_size)
#' which assigns a segment to its START bin only. A 50 Mb deletion therefore
#' contributed to exactly one 1 Mb bin instead of fifty. Because large segments
#' carry most of the real copy-number signal, the genome-wide frequency plot
#' systematically under-represented exactly what it was meant to show.
expand_segments_to_bins <- function(seg, bin_size = 1e6, max_rows = 5e7) {
  d <- as.data.table(seg)
  d <- d[end > start]
  b0 <- floor(d$start / bin_size)
  b1 <- floor((d$end - 1) / bin_size)
  n  <- as.integer(b1 - b0 + 1L)

  total <- sum(as.numeric(n))
  if (total > max_rows) {
    mohq_die(sprintf(
      "Binning would create %.2g rows (bin_size = %g). Increase --bin_size.",
      total, bin_size))
  }

  # rep.int + sequence() is vectorised and works on every R >= 3.x.
  idx <- rep.int(b0, n) + (sequence(n) - 1L)

  out <- data.table(
    patient_id = rep.int(d$patient_id, n),
    chr        = rep.int(d$chr,        n),
    bin        = idx * bin_size,
    type       = rep.int(d$type,       n)
  )
  mohq_log(sprintf("Expanded %s segments into %s bins of %.0f kb",
                   format(nrow(d), big.mark = ","),
                   format(nrow(out), big.mark = ","), bin_size / 1e3))
  out[]
}

# --------------------------------------------------------------------------- #
# Gene coordinates
# --------------------------------------------------------------------------- #

#' Load real gene coordinates from a GTF (optionally gzipped) or a BED/TSV.
#'
#' FIX: run_comparative_analysis.R derived "gene coordinates" as
#'     min(Start_Position) / max(End_Position) of that gene's mutations in the MAF.
#' That is the span of observed mutations, not the gene body. For a gene with a
#' single observed mutation the interval collapses to a few base pairs, so the
#' segment overlap almost never fires; for a gene with two distant mutations it
#' can span megabases. Either way the gene-level CNV calls were not measuring
#' the gene. Real annotation is required -- there is no safe fallback, so this
#' function errors rather than degrading silently.
load_gene_coords <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
    mohq_die("A gene annotation file is required (--gtf). Supply the same GRCh38 ",
             "GTF used by the RNA pipeline, or a 4-column BED (chr,start,end,gene).")
  }
  is_gtf <- grepl("\\.gtf(\\.gz)?$", path, ignore.case = TRUE)

  if (is_gtf) {
    g <- fread(path, sep = "\t", header = FALSE, quote = "", skip = "\t",
               col.names = c("chr", "src", "feature", "start", "end",
                             "score", "strand", "frame", "attr"))
    g <- g[feature == "gene"]
    if (!nrow(g)) mohq_die("No 'gene' features found in GTF: ", path)
    g[, gene := sub('.*gene_name "([^"]+)".*', "\\1", attr)]
    # Ensembl GTFs without gene_name fall back to gene_id.
    no_name <- !grepl("gene_name", g$attr)
    if (any(no_name)) {
      g[no_name, gene := sub('.*gene_id "([^"]+)".*', "\\1", attr)]
    }
    g <- g[, .(chr, start, end, gene)]
  } else {
    g <- fread(path, header = TRUE)
    if (ncol(g) < 4) {
      g <- fread(path, header = FALSE)
      if (ncol(g) < 4) mohq_die("Gene BED needs >= 4 columns (chr,start,end,gene): ", path)
    }
    setnames(g, 1:4, c("chr", "start", "end", "gene"))
    g <- g[, .(chr, start, end, gene)]
  }

  g[, chr := normalise_chr(chr)]
  g <- g[!is.na(gene) & nzchar(gene) & !is.na(start) & !is.na(end)]
  # Collapse duplicate gene symbols (PAR regions, patches) to their widest span
  # on their most common chromosome.
  g <- g[, .(start = min(start), end = max(end)),
         by = .(gene, chr)][, .SD[which.max(end - start)], by = gene]
  mohq_log(sprintf("Loaded coordinates for %s genes from %s",
                   format(nrow(g), big.mark = ","), basename(path)))
  setkey(g, chr, start, end)
  g[]
}

#' Assign a copy-number call to each (patient, gene) via true interval overlap.
#' Uses data.table::foverlaps -- vectorised, unlike the previous per-gene for()
#' loop wrapped around a rowwise() filter.
map_segments_to_genes <- function(seg, gene_coords, genes = NULL) {
  gc2 <- if (is.null(genes)) copy(gene_coords) else gene_coords[gene %in% genes]
  if (!nrow(gc2)) mohq_die("None of the requested genes are present in the annotation.")

  s <- copy(as.data.table(seg))[, .(patient_id, chr, start, end, logr)]
  setkey(s, chr, start, end)
  setkey(gc2, chr, start, end)

  ov <- foverlaps(gc2, s, type = "any", nomatch = NULL)
  if (!nrow(ov)) mohq_die("No segment/gene overlaps found -- check chromosome naming on both sides.")

  # Where several segments hit one gene, keep the most extreme call. This
  # matches the previous intent but is done with a proper grouped max.
  ov[, absr := abs(logr)]
  res <- ov[order(-absr)][, .SD[1], by = .(patient_id, gene)]
  res[, .(patient_id, gene, chr, logr)]
}

# --------------------------------------------------------------------------- #
# MAF handling
# --------------------------------------------------------------------------- #

MAF_MIN_COLS <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
                  "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification",
                  "Variant_Type", "Tumor_Sample_Barcode")

# Non-synonymous classes -- maftools' default definition.
# --------------------------------------------------------------------------- #
# FLAGS -- FrequentLy mutAted GeneS (Shyr et al., BMC Med Genomics 2014)
#
# These genes appear near the top of the recurrence ranking in essentially every
# exome and genome cohort ever sequenced, and almost never because of biology.
# They are long, repetitive, and sit in segmental duplications, so short reads
# mismap and callers emit false positives.
#
# The first real MoHQ-HM-19 oncoplot had 25 of 25 top genes from this list --
# mucins, WASH and NBPF and GOLGA family members, FLG, HRNR, AHNAK2 -- and not
# one recognisable driver. The figure was not wrong; it was measuring the
# reference genome's difficult regions.
#
# Excluding them is standard practice, but it is a JUDGEMENT, not a fact: a real
# driver in this list would be removed too. So it is a parameter (exclude_flags),
# the removed genes are reported, and the unfiltered ranking is still written to
# disk next to the filtered one.
# --------------------------------------------------------------------------- #
FLAGS_GENES <- c(
  "TTN","MUC16","OBSCN","AHNAK2","SYNE1","FLG","MUC5B","DNAH5","GPR98","FAT3",
  "PKHD1L1","FAT4","DNAH2","CDH23","DNAH3","MYH13","DNAH8","DNAH1","DNAH9",
  "PKD1","MDN1","RNF213","RYR1","DNAH10","USH2A","DNAH17","DNAH11","HMCN1",
  "MUC17","ZNF469","MUC2","FSIP2","MUC12","MUC4","MUC5AC","MUC6","MUC19",
  "HRNR","PLEC","SYNE2","NEB","CSMD1","CSMD3","LRP1B","XIRP2","VPS13D",
  "AHNAK","MACF1","SPTA1","PCLO","RYR2","RYR3","ABCA13","APOB","COL6A3",
  "DST","DYNC1H1","FCGBP","HYDIN","LAMA5","MYO18B","NBEA","SACS","SZT2",
  "TENM1","TRIO","UBR4","VCAN","ZFHX4","EPPK1","PRB2","LILRB3","TUBB8B",
  # segmental-duplication families seen at the top of the first MoHQ oncoplot
  "WASH6P","WASHC1","WASH1","WASH2P","WASH3P","WASH4P","WASH5P","WASH7P",
  "NBPF1","NBPF3","NBPF8","NBPF9","NBPF10","NBPF12","NBPF14","NBPF15","NBPF20",
  "GOLGA6B","GOLGA6L1","GOLGA6L2","GOLGA6L4","GOLGA6L5P","GOLGA6L6","GOLGA6L9",
  "GOLGA8A","GOLGA8B","GOLGA8K","GOLGA8M","GOLGA8N","GOLGA8O",
  "FOXD4L1","FOXD4L3","FOXD4L4","FOXD4L5","FOXD4L6","KLF18","RAMEF10","PRIM2",
  "ANKRD20A1","ANKRD20A2","ANKRD20A3","ANKRD20A4","TRD20A4P","PDE4DIP","NOTCH2NL"
)

# A list of individual symbols cannot win this.
#
# After the first FLAGS pass removed the mucins and the WASH/GOLGA/NBPF members
# that were in the list, the next 25 genes were the SAME PROBLEM one tier down:
# other members of the same families (NBPF11, NPIPA5, NPIPB11, GOLGA6L10,
# ANKRD36C, USP17L22, POTEH, RAMEF15), plus olfactory receptors and other
# tandem-duplicated clusters. Adding those 25 by name would produce a third
# tier.
#
# These families are large, and their members are numbered. Match the family.
#
# I also learned to read symbols from the DATA, not from a plot: "TRD20A4P" went
# into the list above from a clipped axis label. The real symbol is ANKRD20A4P,
# so the entry never matched anything and the gene stayed at the top.
FLAGS_PATTERNS <- c(
  "^ANKRD20A",     # segmental duplication family
  "^ANKRD36",
  "^NBPF[0-9]",    # neuroblastoma breakpoint family
  "^NPIP[AB]",     # nuclear pore interacting protein family
  "^GOLGA[0-9]",   # golgins in segmental duplications
  "^WASH[0-9]",    # WASH family and pseudogenes
  "^USP17L",       # tandem repeat array
  "^POTE[A-Z]",
  "^RAMEF",
  "^FAM90A",
  "^TBC1D3",
  "^OR[0-9]+[A-Z][0-9]",   # olfactory receptors: large, polymorphic, low mappability
  "^LILR[AB][0-9]",        # LILR cluster
  "^MUC[0-9]",             # mucins
  "^HLA-",                 # HLA: high polymorphism causes systematic miscalls
  "^IG[HKL][VDJC]",        # immunoglobulin loci -- somatically rearranged
  "^TR[ABGD][VDJC]",       # T-cell receptor loci -- likewise
  "^PRAMEF",
  "^KRTAP",                # keratin-associated, tandem arrays
  "^ZNF7[0-9][0-9]",       # some large ZNF clusters
  "^DEFB[0-9]",
  "^NUTM2",
  "-AS[0-9]$", "^LINC[0-9]", "^LOC[0-9]"   # non-coding / provisional symbols
)

#' Drop artefact-prone genes, by explicit symbol AND by family pattern.
#'
#' Reports what went, because a filter that removes a third of the data
#' silently is not a filter, it is a hidden assumption.
drop_flags <- function(dt, extra = character(), patterns = FLAGS_PATTERNS,
                       report = TRUE) {
  flags <- union(FLAGS_GENES, extra)
  hit <- dt$Hugo_Symbol %in% flags
  if (length(patterns)) {
    rx <- paste(patterns, collapse = "|")
    hit <- hit | grepl(rx, dt$Hugo_Symbol)
  }
  if (report && any(hit)) {
    tab <- sort(table(dt$Hugo_Symbol[hit]), decreasing = TRUE)
    mohq_log(sprintf("FLAGS filter: removed %s variant(s) in %d gene(s); most frequent: %s",
                     format(sum(hit), big.mark = ","), length(tab),
                     paste(utils::head(names(tab), 8), collapse = ", ")))
  }
  dt[!hit]
}

#' Restrict to a curated gene list.
#'
#' The defensible way to build a cohort oncoplot. A blocklist answers "which
#' genes do I not believe?", which is open-ended. An allowlist answers "which
#' genes am I making a claim about?", which is finite and citable -- the COSMIC
#' Cancer Gene Census or an OncoKB panel, whichever the group already uses.
#'
#' Supply a file with one HGNC symbol per line (lines starting with # ignored).
restrict_to_panel <- function(dt, panel_file) {
  if (is.null(panel_file) || !nzchar(panel_file) || !file.exists(panel_file)) return(dt)
  genes <- readLines(panel_file, warn = FALSE)
  genes <- trimws(genes[!grepl("^\\s*#", genes) & nzchar(trimws(genes))])
  before <- nrow(dt)
  out <- dt[Hugo_Symbol %in% genes]
  mohq_log(sprintf("Gene panel %s: %d symbols; kept %s of %s variants (%.1f%%), %d genes present",
                   basename(panel_file), length(genes),
                   format(nrow(out), big.mark = ","), format(before, big.mark = ","),
                   100 * nrow(out) / max(before, 1), uniqueN(out$Hugo_Symbol)))
  if (!nrow(out))
    mohq_die("No variant fell in the supplied gene panel. Check the symbols match ",
             "the MAF's Hugo_Symbol (HGNC), and that the file is one symbol per line.")
  out
}

# GRCh38 primary autosome lengths (bp). Used so a genome-wide plot allocates
# vertical space by real chromosome length rather than by how many bins happen
# to carry data -- otherwise chr6 can appear larger than chr1.
GRCH38_CHROM_LEN <- c(
  `1`=248956422, `2`=242193529, `3`=198295559, `4`=190214555, `5`=181538259,
  `6`=170805979, `7`=159345973, `8`=145138636, `9`=138394717, `10`=133797422,
  `11`=135086622, `12`=133275309, `13`=114364328, `14`=107043718, `15`=101991189,
  `16`=90338345,  `17`=83257441,  `18`=80373285,  `19`=58617616,  `20`=64444167,
  `21`=46709983,  `22`=50818468,  X=156040895,   Y=57227415)

NONSYN_CLASSES <- c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
                    "Missense_Mutation", "Nonsense_Mutation", "Splice_Site",
                    "Translation_Start_Site", "Nonstop_Mutation")

#' Read and concatenate somatic MAFs, keyed on patient_id.
#'
#' `files` must be DNA (_D) MAFs only. The old scripts globbed `*.maf`, which
#' also matched the `_R` RNA MAFs -- those are RNA-derived calls, are 10-20x
#' larger, and do not belong in a somatic DNA oncoplot. File selection is now
#' the manifest's job, so this function just validates what it is handed.
load_somatic_mafs <- function(files, manifest = NULL, nonsyn_only = TRUE,
                              extra_cols = character()) {
  files <- files[file.exists(files) & file.info(files)$size > 0]
  if (!length(files)) mohq_die("No readable MAF files supplied.")

  bad <- grep("_R\\.|_R_", basename(files), value = TRUE)
  if (length(bad)) {
    mohq_die("RNA-derived MAF(s) passed to a somatic loader: ",
             paste(bad, collapse = ", "),
             "\nPass only _D MAFs (manifest column `somatic_maf`).")
  }

  want <- unique(c(MAF_MIN_COLS, extra_cols))
  dt <- rbindlist(lapply(files, function(f) {
    hdr <- names(fread(f, skip = "Hugo_Symbol", nrows = 0))
    sel <- intersect(want, hdr)
    missing_req <- setdiff(MAF_MIN_COLS, hdr)
    if (length(missing_req)) {
      mohq_die("MAF ", basename(f), " lacks required column(s): ",
               paste(missing_req, collapse = ", "))
    }
    d <- fread(f, skip = "Hugo_Symbol", select = sel, fill = TRUE,
               showProgress = FALSE)
    d[, source_file := basename(f)]
    d
  }), fill = TRUE)

  dt[, Chromosome := normalise_chr(Chromosome)]
  dt[, patient_id := harmonise_patient_id(Tumor_Sample_Barcode)]
  # Downstream everything is per patient, so make the barcode canonical too.
  dt[, Tumor_Sample_Barcode := patient_id]

  # ---------------------------------------------------------------------- #
  # GUARD: the MAFs shipped inside the MoHQ collection are NOT gene-annotated
  # -- Hugo_Symbol is "Unknown" for essentially every row. Those files must be
  # regenerated from the PCGR VCFs with vcf2maf + VEP (the VCF2MAF process does
  # this) before any mutation analysis.
  #
  # Without this check the downstream `Hugo_Symbol != "Unknown"` filter would
  # silently discard every variant and hand maftools an empty table -- an
  # oncoplot with no genes and no error. Fail loudly instead.
  # ---------------------------------------------------------------------- #
  # JUDGE BY DISTINCT GENE SYMBOLS, NOT BY THE 'Unknown' FRACTION.
  #
  # This check previously failed above 50% Unknown, which is WRONG for whole
  # genome data: most somatic variants in WGS are intergenic, and "Unknown" is
  # the correct value for them. A correctly annotated MoHQ cohort comes in
  # around 55% Unknown with several thousand distinct genes, and the old
  # threshold rejected it -- good data, refused, with a message confidently
  # blaming the input.
  #
  # The failure actually worth catching is the delivered MoHQ MAFs, where
  # Hugo_Symbol is "Unknown" for EVERY row. That shows up as zero distinct
  # symbols, not as a fraction near one half.
  #
  # Same reasoning and same thresholds as the guard in modules/local/vcf2maf.nf;
  # keep the two in step.
  unknown_frac <- mean(is.na(dt$Hugo_Symbol) | dt$Hugo_Symbol %in% c("Unknown", "", "."))
  n_genes <- length(unique(dt$Hugo_Symbol[
    !is.na(dt$Hugo_Symbol) & !dt$Hugo_Symbol %in% c("Unknown", "", ".")]))

  min_genes <- as.integer(Sys.getenv("MOHQ_MIN_ANNOTATED_GENES", "200"))
  max_frac  <- as.numeric(Sys.getenv("MOHQ_MAX_UNKNOWN_FRAC", "0.95"))

  if (n_genes < min_genes || unknown_frac > max_frac) {
    mohq_die(sprintf(
      paste0("These MAFs are not gene-annotated: only %d distinct gene symbol(s) ",
             "across %s variants (%.1f%% 'Unknown').\n",
             "  Expected several thousand distinct genes. The MAFs distributed in ",
             "the MoHQ collection\n",
             "  are unannotated and must be regenerated from the PCGR VCFs. Point ",
             "--mafs at the\n",
             "  VCF2MAF output (derived_maf_dir), not at the delivered ",
             "`somatic_maf`.\n",
             "  Files given: %s"),
      n_genes, format(nrow(dt), big.mark = ","), 100 * unknown_frac,
      paste(basename(files), collapse = ", ")))
  }
  mohq_log(sprintf("Annotation check: %s distinct genes, %.1f%% Unknown (intergenic variants are expected to be Unknown)",
                   format(n_genes, big.mark = ","), 100 * unknown_frac))

  n_all <- nrow(dt)
  if (nonsyn_only) {
    dt <- dt[Variant_Classification %in% NONSYN_CLASSES]
    mohq_log(sprintf("Kept %s/%s non-synonymous variants",
                     format(nrow(dt), big.mark = ","), format(n_all, big.mark = ",")))
  }
  dt <- dt[!is.na(Hugo_Symbol) & nzchar(Hugo_Symbol) & Hugo_Symbol != "Unknown"]

  if (!is.null(manifest)) check_ids_against_manifest(dt$patient_id, manifest, "MAF")
  mohq_log(sprintf("Loaded %s variants across %d patients from %d MAF file(s)",
                   format(nrow(dt), big.mark = ","),
                   uniqueN(dt$patient_id), length(files)))
  dt[]
}

#' Tumour mutational burden per patient.
#' `callable_mb` defaults to the standard WES-ish 30 Mb; for WGS pass ~2800.
compute_tmb <- function(maf_dt, callable_mb = 30) {
  tmb <- maf_dt[, .(n_nonsyn = .N), by = patient_id]
  tmb[, tmb_per_mb := n_nonsyn / callable_mb]
  tmb[order(-tmb_per_mb)]
}

# --------------------------------------------------------------------------- #
# PCGR CNA gene table
# --------------------------------------------------------------------------- #

#' Read PCGR per-gene CNA segment tables into a maftools-compatible cnTable.
#'
#' Two fixes over the previous process_pcgr():
#'
#'  1. `filter(oncogene == "TRUE" | tumor_suppressor == "TRUE")` compared a
#'     column to the STRING "TRUE". data.table::fread parses those columns as
#'     logical, so the comparison was FALSE for every row and every CNV was
#'     silently dropped -- an empty cnTable with no error. We coerce explicitly.
#'
#'  2. maftools' `cnTable` expects columns Gene / Sample_name / CN. The old code
#'     supplied Hugo_Symbol / Tumor_Sample_Barcode / Variant_Classification, so
#'     even a populated table was not read as intended. Verify against your
#'     maftools version -- `?read.maf` documents the expected names.
#'
#' Column names differ between PCGR releases, so they are resolved by pattern
#' with an explicit error listing what was actually found.
PCGR_CNA_PATTERNS <- list(
  gene   = "^(symbol|gene_symbol|gene|hugo_symbol)$",
  sample = "^(sample_id|sample|tumor_sample_barcode)$",
  logr   = "^(log_?r|logr_?mean|log2_?ratio|segment_?mean|mean_?log2_?ratio)$",
  onco   = "^(oncogene)$",
  tsg    = "^(tumor_?suppressor|tumour_?suppressor)$"
)

as_logical_loose <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  xx <- tolower(trimws(as.character(x)))
  xx %in% c("true", "t", "yes", "y", "1")
}

read_pcgr_cna <- function(files, amp_threshold, del_threshold,
                          driver_only = TRUE, manifest = NULL) {
  files <- files[file.exists(files) & file.info(files)$size > 0]
  if (!length(files)) {
    mohq_warn("No PCGR CNA files supplied; continuing without copy-number annotation.")
    return(data.table(Gene = character(), Sample_name = character(), CN = character()))
  }

  out <- rbindlist(lapply(files, function(f) {
    d <- fread(f, fill = TRUE, showProgress = FALSE)
    norm <- tolower(gsub("[^A-Za-z0-9_]", "_", names(d)))
    pick <- function(pat) {
      hit <- which(grepl(pat, norm))
      if (length(hit)) names(d)[hit[1]] else NA_character_
    }
    cg <- pick(PCGR_CNA_PATTERNS$gene)
    cs <- pick(PCGR_CNA_PATTERNS$sample)
    cl <- pick(PCGR_CNA_PATTERNS$logr)
    if (any(is.na(c(cg, cl)))) {
      mohq_die("Cannot resolve PCGR CNA columns in ", basename(f),
               "\n  Found: ", paste(names(d), collapse = ", "),
               "\n  Update PCGR_CNA_PATTERNS in mohq_common.R.")
    }
    co <- pick(PCGR_CNA_PATTERNS$onco)
    ct <- pick(PCGR_CNA_PATTERNS$tsg)

    res <- data.table(
      Gene        = as.character(d[[cg]]),
      Sample_name = if (is.na(cs)) harmonise_patient_id(basename(f)) else as.character(d[[cs]]),
      logr        = suppressWarnings(as.numeric(d[[cl]])),
      is_onco     = if (is.na(co)) FALSE else as_logical_loose(d[[co]]),
      is_tsg      = if (is.na(ct)) FALSE else as_logical_loose(d[[ct]])
    )
    res[!is.na(logr)]
  }), fill = TRUE)

  n_all <- nrow(out)
  out <- out[logr >= amp_threshold | logr <= del_threshold]
  if (driver_only) {
    n_pre <- nrow(out)
    keep <- out[is_onco | is_tsg]
    if (!nrow(keep)) {
      mohq_warn("Driver filter removed all CNVs (", n_pre, " -> 0). ",
                "The oncogene/tumor_suppressor columns are likely absent or ",
                "differently encoded in this PCGR version; keeping all CNVs instead.")
    } else {
      out <- keep
    }
  }

  out[, CN := fifelse(logr >= amp_threshold, "Amp", "Del")]
  out[, Sample_name := harmonise_patient_id(Sample_name)]
  out <- unique(out[, .(Gene, Sample_name, CN)])
  out <- out[!is.na(Gene) & nzchar(Gene)]

  if (!is.null(manifest)) check_ids_against_manifest(out$Sample_name, manifest, "PCGR CNA")
  mohq_log(sprintf("PCGR CNA: %s gene-level calls retained (from %s rows, %d file(s))",
                   format(nrow(out), big.mark = ","),
                   format(n_all, big.mark = ","), length(files)))
  out[]
}

# --------------------------------------------------------------------------- #
# Numeric helpers
# --------------------------------------------------------------------------- #

#' Column variances without pulling in matrixStats.
#' The old code used `apply(expr_matrix, 2, var)` on a data.frame, which coerces
#' the whole thing per column and is very slow on a ~60k-gene matrix.
col_vars <- function(m) {
  n <- nrow(m)
  if (n < 2) return(rep(0, ncol(m)))
  mu <- colMeans(m, na.rm = TRUE)
  colSums((m - rep(mu, each = n))^2, na.rm = TRUE) / (n - 1)
}

# --------------------------------------------------------------------------- #
# Plot helpers
# --------------------------------------------------------------------------- #

MOHQ_PALETTE <- c(Gain = "#D55E00", Loss = "#0072B2",
                  Amp  = "#D55E00", Del  = "#0072B2",
                  F    = "#CC79A7", M    = "#009E73")

save_plot <- function(plot, file, width, height, dpi = 300) {
  ggplot2::ggsave(file, plot = plot, width = width, height = height,
                  dpi = dpi, limitsize = FALSE)
  mohq_log("Wrote ", file)
}

#' Write an empty-but-valid placeholder so the report never breaks on a
#' missing panel. Replaces the previous `ifEmpty(file("no_oncoprint.png"))`
#' pattern in main.nf, which referenced files that did not exist on disk.
write_placeholder_png <- function(file, message_text) {
  p <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = message_text, size = 5) +
    ggplot2::theme_void()
  ggplot2::ggsave(file, p, width = 8, height = 3, dpi = 150)
  mohq_warn("Wrote placeholder: ", file, " (", message_text, ")")
}

mohq_session_info <- function(file = "versions.txt") {
  si <- utils::sessionInfo()
  pkgs <- sort(vapply(si$otherPkgs, function(p) paste0(p$Package, " ", p$Version), character(1)))
  writeLines(c(paste("R", getRversion()), pkgs), file)
  invisible(file)
}

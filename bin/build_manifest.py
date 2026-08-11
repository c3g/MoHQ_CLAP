#!/usr/bin/env python3
"""
build_manifest.py -- walk a MoHQ delivery tree and emit a canonical sample manifest
                     plus a data-completeness matrix.

This replaces the fuzzy prefix-regex sample matching that was previously scattered
across the R scripts. Every downstream script joins on `patient_id`, which is
derived here, once, deterministically.

Tree layout assumed (verified against the MoHQ delivery tree):

    <root>/
      <INSTITUTION>/            e.g. CM, CQ, IQ
        <COHORT>/               e.g. MoHQ-CM-4, MoHQ-CM-60, MoHQ-CQ-34
          <PATIENT>/            e.g. MoHQ-CM-4-10, MoHQ-CM-60-2
            alignment/ expression/ parameters/ raw_cnv/ raw_data/
            reports/{pcgr/} svariants/{linx/} variants/{caller_vcfs/}
            Key_metrics.csv log.csv Methods.html Readme.html Warnings.html

Patient dirs are detected by *content* (presence of known subdirectories), not by
depth, so the walker still works if an institution adds or removes a nesting level.

ID grammar (derived empirically from the delivery tree -- note the deliberate
tolerance for alphanumeric sample numbers):

    MoHQ - <INST> - <cohort_num> - <patient_num> - <sample_num> - <aliquot><TYPE>
    \________________________________________/    \___________/   \_____/\____/
                  patient_id                       alphanumeric!    int   DN|DT|RN|RT

    cohort_id  = first 3 hyphen tokens   (MoHQ-CM-4,  MoHQ-CM-60)
    patient_id = first 4 hyphen tokens   (MoHQ-CM-4-10, MoHQ-CM-60-2)

Real-world quirks this handles explicitly (all observed in the tree):
  * sample_num is NOT always numeric        -> RCC01, RCC01n
  * normal and tumour have DIFFERENT sample_nums -> 105-305328-1DN vs 105-305329-2DT
  * aliquot index does not always start at 1 -> MoHQ-CM-4-105 has 2DT and no 1DT
  * PCGR filenames use three different conventions across pipeline versions:
        <patient>_D.acmg.grch38.maf
        <patient>_D.pcgr_acmg.grch38.maf
        <patient>.acmg.grch38.maf          <-- no _D/_R suffix, AMBIGUOUS, flagged
  * some patients have several TumourPair runs (.1.ini ... .4.ini)
  * entire assay classes can be missing per cohort (CQ has no RNA at all)

Outputs (all TSV, written to --outdir):
    <cohort>.manifest.tsv     one row per patient, absolute paths per asset
    <cohort>.completeness.tsv one row per patient, yes/no per asset
    <cohort>.ambiguities.tsv  every asset where >1 file matched, with the choice made
    <cohort>.samples.tsv      one row per *sample* (DN/DT/RT), for QC joins

Stdlib only -- no third-party deps, so it runs anywhere.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import json
import os
import re
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

# --------------------------------------------------------------------------- #
# Scale note (4,500 patients / ~65 cohorts on Lustre)
#
# The first version of this script called Path.glob() once per asset pattern:
# 51 patterns x 4,500 patients = ~230,000 directory listings. On a local disk
# that is ~30 s. On Lustre /scratch, where each readdir is a network metadata
# round trip, it is tens of minutes of pure metadata traffic -- the access
# pattern Alliance support flags as abusive.
#
# Now: ONE os.scandir per patient subdirectory (11 per patient), cached in
# memory, with every pattern matched via fnmatch against that listing. That is
# a 4.6x reduction in directory listings.
#
# On top of that, --cache stores a per-patient signature built from the mtime
# of each subdirectory. A directory's mtime changes when files are added or
# removed, which is exactly what a manifest records, so an unchanged patient is
# detected with ~11 cheap stat() calls and skipped entirely -- no readdir at
# all. Incremental runs therefore touch only new or modified patients.
# --------------------------------------------------------------------------- #

# Every directory we ever need to list. Asset patterns are validated against
# this at import time so a new pattern cannot silently fall outside the scan.
SCAN_DIRS: tuple[str, ...] = (
    "",                       # patient root: Key_metrics.csv, log.csv, *.html
    "alignment",
    "expression",
    "parameters",
    "raw_cnv",
    "raw_data",
    # Some deliveries put RNA alignments and RNA variant calls under
    # rna_variants/ rather than alignment/ and variants/. Missed entirely until
    # a comparison with the object-store audit script surfaced it; every RNA
    # asset was reported absent for those patients.
    "rna_variants",
    "reports",
    "reports/pcgr",
    "svariants",
    "svariants/linx",
    "variants",
    "variants/caller_vcfs",
)

# --------------------------------------------------------------------------- #
# ID parsing
# --------------------------------------------------------------------------- #

# sample_num is [A-Za-z0-9]+ on purpose: CQ uses RCC01 / RCC01n, IQ uses
# 0016225, GC uses BL01, MU uses FF1.
#
# seqtype is [A-Z]{0,2}[DR][NT] rather than a fixed DN|DT|RN|RT list, because
# real IDs carry a tissue-qualifier prefix that a fixed list silently rejects:
#     MoHQ-GC-31-576-BL01-1FDN   <-- 'FDN', not 'DN'
#     MoHQ-MU-8-2-FF1-1DT        <-- plain 'DT'
# With the old fixed list, every GC sample named ...1FDN failed to parse and was
# dropped from discover_samples() with no error -- the patient would appear to
# have no normal sample at all.
SAMPLE_RE = re.compile(
    r"^(?P<collection>[A-Za-z]+)"
    r"-(?P<institution>[A-Za-z]+)"
    r"-(?P<cohort_num>[A-Za-z0-9]+)"
    r"-(?P<patient_num>[A-Za-z0-9]+)"
    r"-(?P<sample_num>[A-Za-z0-9]+)"
    r"-(?P<aliquot>\d+)"
    r"(?P<seqtype>[A-Z]{0,2}[DR][NT])$"
)

PATIENT_RE = re.compile(
    r"^(?P<collection>[A-Za-z]+)"
    r"-(?P<institution>[A-Za-z]+)"
    r"-(?P<cohort_num>[A-Za-z0-9]+)"
    r"-(?P<patient_num>[A-Za-z0-9]+)$"
)

# Subdirectories that mark a directory as a patient deliverable.
PATIENT_MARKER_DIRS = {"alignment", "variants", "reports", "raw_data", "parameters"}

# Filename suffixes that carry a full sample ID we can harvest.
SAMPLE_BEARING_SUFFIXES = (
    ".bam", ".bam.bai", ".variants.bam",
    ".abundance_genes.tsv", ".abundance_transcripts.tsv",
    ".anno_fuse.tsv",
)


def parse_sample_id(sid: str) -> dict | None:
    """Parse a full sample ID. Returns None if it does not match the grammar."""
    m = SAMPLE_RE.match(sid)
    if not m:
        return None
    d = m.groupdict()
    d["sample_id"] = sid
    d["patient_id"] = "-".join(sid.split("-")[:4])
    d["cohort_id"] = "-".join(sid.split("-")[:3])
    d["aliquot"] = int(d["aliquot"])
    # Derive from the D/R and N/T characters, not from the whole token, so a
    # qualifier prefix like 'FDN' still classifies correctly.
    st = d["seqtype"]
    d["molecule"] = "DNA" if "D" in st else "RNA"
    d["tissue"] = "normal" if st.endswith("N") else "tumour"
    # Anything before the D/R is a tissue qualifier (F = FFPE/fresh-frozen etc.)
    d["seqtype_qualifier"] = st[:-2] or ""
    return d


# --------------------------------------------------------------------------- #
# Genome build
#
# GRCh37 and GRCh38 coordinates are NOT interchangeable, so mixing builds inside
# one cohort silently corrupts every coordinate-based analysis -- the .seg file,
# the CNV frequency plot, the gene overlap. Nothing errors; the numbers are just
# wrong. Hence recording it per patient.
#
# *** FILENAMES ARE NOT TRUSTWORTHY FOR THIS. ***
# Every object in the MoHQ collection is named "...grch38..." (27,268 of them,
# zero mentioning grch37), yet a delivery report states MoHQ-CM-1-1 was ALIGNED
# against GRCh37. So the name asserts something the data may contradict.
#
# The value below is therefore labelled `genome_build_from_filename` and is a
# CLAIM, not a fact. Verify with bin/check_genome_build.sh, which reads the
# actual chromosome lengths out of the BAM and VCF headers, and pass the result
# back via --build-table. When a verified build is supplied it always wins.
# --------------------------------------------------------------------------- #
BUILD_PATTERNS = (
    (re.compile(r"\bgrch38\b|\bhg38\b", re.I), "GRCh38"),
    (re.compile(r"\bgrch37\b|\bhg19\b", re.I), "GRCh37"),
)


def detect_genome_build(paths: Iterable[str]) -> str:
    """Infer the assembly from the resolved filenames. 'unknown' if silent."""
    found: set[str] = set()
    for p in paths:
        if not p or p == "NA":
            continue
        base = os.path.basename(p)
        for rx, build in BUILD_PATTERNS:
            if rx.search(base):
                found.add(build)
    if len(found) == 1:
        return found.pop()
    if len(found) > 1:
        return "MIXED:" + "+".join(sorted(found))
    return "unknown"


# --------------------------------------------------------------------------- #
# Processing provenance
#
# GenPipes writes its full config beside each patient, recording every module
# version used. This is a BATCH VARIABLE and belongs in the manifest.
#
# Confirmed in MoHQ-CM-1: patients were processed in at least two rounds --
# PCGR 0.9.2 for some, PCGR 1.0.3 for others, WITHIN THE SAME COHORT. PCGR
# renamed output columns between major versions, and differences between the
# two groups are technical, not biological.
# --------------------------------------------------------------------------- #
MODULE_RE = re.compile(r"^\s*module_(\w+)\s*=\s*(\S+)", re.M)


def read_ini_versions(ini_path: str) -> dict[str, str]:
    """Extract module versions from a GenPipes .ini. Cheap: these are small."""
    if not ini_path or ini_path == "NA":
        return {}
    try:
        with open(ini_path, "r", errors="replace") as fh:
            text = fh.read(400_000)          # configs are small; cap defensively
    except OSError:
        return {}
    out = {}
    for m in MODULE_RE.finditer(text):
        out[m.group(1).lower()] = m.group(2)
    return out


def patient_id_of(name: str) -> str | None:
    """Return the patient_id embedded in a directory or file name, else None."""
    stem = name.split(".")[0]
    if PATIENT_RE.match(stem):
        return stem
    parsed = parse_sample_id(stem)
    if parsed:
        return parsed["patient_id"]
    # e.g. "MoHQ-CM-4-10_D" from PCGR barcodes
    stem2 = re.sub(r"_(D|R|T|N|DNA|RNA)$", "", stem)
    if PATIENT_RE.match(stem2):
        return stem2
    return None


# --------------------------------------------------------------------------- #
# Asset catalogue
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Asset:
    """One logical data product we expect to find for a patient.

    patterns are glob patterns relative to the patient directory, tried in
    priority order. `{p}` is substituted with the patient_id. The FIRST pattern
    that yields any match wins; within that pattern, matches are sorted and
    `prefer`/`avoid` substrings break ties deterministically.
    """
    name: str
    patterns: tuple[str, ...]
    required: bool = False          # counted in the "core completeness" score
    prefer: tuple[str, ...] = ()    # substrings that win a tie
    avoid: tuple[str, ...] = ()     # substrings that lose a tie
    note: str = ""


# NOTE ON MAF GLOBS ---------------------------------------------------------
# `_D` = variants called from tumour/normal DNA (what you want for a somatic
# cohort oncoplot). `_R` = variants called from RNA -- these files are 10-20x
# larger and are NOT somatic DNA calls. The old `list.files(pattern="\\.maf$")`
# swept both into the same oncoplot. They are separated here, deliberately.
ASSETS: tuple[Asset, ...] = (
    # ---- somatic DNA (the analytic backbone) ----
    #
    # THE PCGR VCF IS THE REAL MUTATION INPUT, NOT THE MAF.
    # The MAFs distributed in the collection are not gene-annotated: Hugo_Symbol
    # is "Unknown" throughout. They must be regenerated from these VCFs with
    # vcf2maf + VEP (the VCF2MAF pipeline process). `somatic_maf` below is kept
    # only for provenance and is deliberately NOT `required`.
    Asset("pcgr_vcf",
          ("reports/pcgr/{p}_D.*.vcf.gz", "reports/pcgr/{p}_D*.vcf.gz"),
          required=False,      # see mutation_vcf fallback below
          avoid=("missing_topups",),
          note="PCGR somatic VCF -- preferred input to vcf2maf."),
    Asset("somatic_maf",
          ("reports/pcgr/{p}_D.*.maf", "reports/pcgr/{p}_D*.maf"),
          required=False,
          avoid=("missing_topups",),
          note="Collection MAF: UNANNOTATED (Hugo_Symbol='Unknown'). Provenance only."),
    # NOT required: no analysis script reads it. Marking it required would
    # exclude patients for lacking a file we never open. Harvested anyway --
    # it is the natural input for any future tier / actionability analysis.
    Asset("somatic_tiers",
          ("reports/pcgr/{p}_D.*.snvs_indels.tiers.tsv", "reports/pcgr/{p}_D*tiers.tsv"),
          required=False,
          avoid=("missing_topups",),
          note="PCGR tier calls. Not consumed by the current analyses."),
    Asset("cna_segments",
          ("reports/pcgr/{p}_D.*.cna_segments.tsv.gz", "reports/pcgr/{p}_D*cna_segments.tsv.gz"),
          required=True,
          avoid=("missing_topups",)),
    # NOT independently required: a patient with a PCGR VCF does not need it.
    # What IS required is `mutation_vcf` -- the PCGR VCF if present, otherwise
    # this. That check happens after asset resolution, so a patient is only
    # excluded when BOTH are absent.
    #
    # Still needed for: the ensemble fallback, RUN_PCGR gap-fill input, and
    # mutational signatures. Harvest it with  --sets core signatures.
    Asset("ensemble_somatic_vcf",
          ("variants/{p}.ensemble.somatic.vt.annot.vcf.gz",),
          required=False,
          note="fallback mutation source, gap-fill input, signature extraction"),
    Asset("ensemble_germline_vcf",
          ("variants/{p}.ensemble.germline.vt.annot.vcf.gz",)),

    # ---- ambiguous legacy naming: no _D / _R suffix ----
    # Observed on MoHQ-CM-4-49. These are promoted into the corresponding
    # somatic_* slot only if that slot is otherwise empty, and always flagged.
    Asset("unsuffixed_maf",
          ("reports/pcgr/{p}.acmg*.maf", "reports/pcgr/{p}.pcgr_acmg*.maf"),
          note="LEGACY: no _D/_R suffix, cannot tell DNA from RNA. Verify manually."),
    Asset("unsuffixed_tiers",
          ("reports/pcgr/{p}.acmg*.snvs_indels.tiers.tsv",
           "reports/pcgr/{p}.pcgr_acmg*.snvs_indels.tiers.tsv"),
          note="LEGACY: no _D/_R suffix."),
    Asset("unsuffixed_cna",
          ("reports/pcgr/{p}.acmg*.cna_segments.tsv.gz",
           "reports/pcgr/{p}.pcgr_acmg*.cna_segments.tsv.gz"),
          note="LEGACY: no _D/_R suffix."),

    # ---- RNA-derived variants (kept separate, NOT for somatic oncoplots) ----
    Asset("rna_maf",
          ("reports/pcgr/{p}_R.*.maf", "reports/pcgr/{p}_R*.maf")),
    Asset("rna_tiers",
          ("reports/pcgr/{p}_R.*.snvs_indels.tiers.tsv", "reports/pcgr/{p}_R*tiers.tsv")),
    Asset("rna_hc_vcf",
          ("variants/*RT.hc.vt.annot.filt.vcf.gz",)),

    # ---- copy number ----
    Asset("cnvkit_vcf", ("raw_cnv/{p}.cnvkit.vcf.gz",), required=True),

    # ---- expression ----
    Asset("expression_genes",
          ("expression/*RT.abundance_genes.tsv",),
          note="kallisto/RSEM gene abundances; absent for whole cohorts (e.g. CQ)"),
    Asset("expression_transcripts", ("expression/*RT.abundance_transcripts.tsv",)),

    # ---- structural variants ----
    Asset("gripss_somatic", ("svariants/*.gripss.filtered.somatic.vcf.gz",)),
    Asset("gripss_germline", ("svariants/*.gripss.filtered.germline.vcf.gz",)),
    Asset("gridss_vcf", ("svariants/*.gridss.vcf.gz", "svariants/*.gridss")),
    Asset("driver_catalog_somatic", ("svariants/*.driver.catalog.somatic.tsv",)),
    Asset("driver_catalog_germline", ("svariants/*.driver.catalog.germline.tsv",)),
    Asset("linx_fusion", ("svariants/linx/*.linx.fusion.tsv",)),
    Asset("linx_drivers", ("svariants/linx/*.linx.driver.catalog.tsv",)),
    Asset("linx_svs", ("svariants/linx/*.linx.svs.tsv",)),
    Asset("linx_clusters", ("svariants/linx/*.linx.clusters.tsv",)),
    Asset("linx_breakend", ("svariants/linx/*.linx.breakend.tsv",)),
    Asset("circos_png", ("svariants/*.circos.png",)),
    Asset("purple_ensemble", ("variants/{p}.purple_ensemble.zip",),
          note="PURPLE purity/ploidy -- unlock for tumour-purity QC"),
    Asset("purple_sv", ("svariants/{p}.purple_sv.zip",)),

    # ---- fusions (RNA) ----
    Asset("anno_fuse", ("reports/*RT.anno_fuse.tsv",)),

    # ---- germline predisposition ----
    Asset("cpsr_zip", ("reports/{p}_D.cpsr.zip",)),

    # ---- QC / provenance ----
    Asset("key_metrics", ("Key_metrics.csv",), required=True),
    Asset("multiqc_dna", ("reports/{p}_D.multiqc.html",)),
    Asset("multiqc_rna", ("reports/{p}_R.multiqc.html",)),
    Asset("pcgr_html_dna", ("reports/{p}_D.pcgr.html", "reports/{p}.pcgr.html")),
    Asset("pcgr_html_rna", ("reports/{p}_R.pcgr.html",)),
    Asset("log_csv", ("log.csv",)),
    Asset("warnings_html", ("Warnings.html",)),
    Asset("methods_html", ("Methods.html",)),
    Asset("tumourpair_ini", ("parameters/{p}.TumourPair*.ini",)),
    Asset("rna_ini", ("parameters/{p}.RNA.*.ini",)),
)

ASSET_BY_NAME = {a.name: a for a in ASSETS}


def _validate_asset_dirs() -> None:
    """Every asset pattern must live in a directory we actually scan."""
    bad = []
    for a in ASSETS:
        for p in a.patterns:
            sub = p.rpartition("/")[0]
            if sub not in SCAN_DIRS:
                bad.append(f"{a.name}: '{p}' -> '{sub}' not in SCAN_DIRS")
    if bad:
        raise RuntimeError("Asset patterns outside SCAN_DIRS:\n  " + "\n  ".join(bad))


_validate_asset_dirs()


# --------------------------------------------------------------------------- #
# Resolution
# --------------------------------------------------------------------------- #

@dataclass
class Ambiguity:
    patient_id: str
    asset: str
    chosen: str
    alternatives: list[str] = field(default_factory=list)
    reason: str = ""


def _rank(name: str, asset: Asset) -> tuple:
    """Deterministic tie-break: avoid-substrings last, prefer-substrings first."""
    avoided = any(a in name for a in asset.avoid)
    preferred = any(p in name for p in asset.prefer)
    return (avoided, not preferred, name)


Listing = dict          # {subdir_relpath: [filename, ...]}


def scan_patient(pdir: Path) -> Listing:
    """One os.scandir per subdirectory. This is the ONLY place we read dirs."""
    out: Listing = {}
    for sub in SCAN_DIRS:
        d = pdir / sub if sub else pdir
        try:
            with os.scandir(d) as it:
                out[sub] = [e.name for e in it]
        except (FileNotFoundError, NotADirectoryError, PermissionError):
            out[sub] = []
    return out


def patient_signature(pdir: Path) -> list:
    """Cheap change-detector: mtime of each scanned subdirectory.

    A directory's mtime changes when entries are added or removed, which is
    precisely what the manifest records (presence of files, not their content).
    ~11 stat() calls, no readdir -- far cheaper than rescanning on Lustre.
    """
    sig = []
    for sub in SCAN_DIRS:
        d = pdir / sub if sub else pdir
        try:
            sig.append([sub, d.stat().st_mtime_ns])
        except (FileNotFoundError, NotADirectoryError, PermissionError):
            sig.append([sub, None])
    return sig


def resolve_asset(listing: Listing, pdir: Path, patient_id: str, asset: Asset,
                  ambiguities: list[Ambiguity]) -> str | None:
    for pattern in asset.patterns:
        sub, _, fpat = pattern.rpartition("/")
        names = listing.get(sub, ())
        if not names:
            continue
        hits = fnmatch.filter(names, fpat.format(p=patient_id))
        if not hits:
            continue
        hits.sort(key=lambda n: _rank(n, asset))
        base = pdir / sub if sub else pdir
        chosen = str(base / hits[0])
        if len(hits) > 1:
            ambiguities.append(Ambiguity(
                patient_id=patient_id,
                asset=asset.name,
                chosen=chosen,
                alternatives=[str(base / h) for h in hits[1:]],
                reason=f"{len(hits)} files matched '{pattern}'",
            ))
        return chosen
    return None


def discover_samples(listing: Listing) -> list[dict]:
    """Harvest every full sample ID visible in this patient directory."""
    seen: dict[str, dict] = {}
    for sub in ("alignment", "expression", "raw_data", "reports", "svariants", "variants"):
        for name in listing.get(sub, ()):
            stem = name
            for suf in sorted(SAMPLE_BEARING_SUFFIXES, key=len, reverse=True):
                if stem.endswith(suf):
                    stem = stem[: -len(suf)]
                    break
            else:
                stem = stem.split(".")[0]
            parsed = parse_sample_id(stem)
            if parsed and parsed["sample_id"] not in seen:
                seen[parsed["sample_id"]] = parsed
    return sorted(seen.values(), key=lambda s: (s["seqtype"], s["aliquot"], s["sample_id"]))


def find_patient_dirs(root: Path) -> list[Path]:
    """Detect patient directories by content, not by depth."""
    out: list[Path] = []
    for dirpath, dirnames, _filenames in os.walk(root):
        p = Path(dirpath)
        if PATIENT_RE.match(p.name) and (set(dirnames) & PATIENT_MARKER_DIRS):
            out.append(p)
            dirnames[:] = []          # do not descend into a patient dir
    return sorted(out)


def read_build_table(path: Path) -> dict[str, str]:
    """patient_id -> verified build, from check_genome_build.sh output.

    Prefers the VCF-derived build: the cohort analyses work on variant and
    segment coordinates, so what the CALLS are on matters more than what the
    alignment was on if the two ever differ.
    """
    out: dict[str, str] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pid = (row.get("patient_id") or "").strip()
            if not pid:
                continue
            for key in ("build_from_vcf", "build_from_bam"):
                v = (row.get(key) or "").strip()
                if v in ("GRCh37", "GRCh38"):
                    out[pid] = v
                    break
    return out


def build_rows(patient_dirs: Iterable[Path],
               cache: dict | None = None,
               build_table: dict[str, str] | None = None
               ) -> tuple[list[dict], list[dict], list[Ambiguity], dict]:
    manifest_rows: list[dict] = []
    sample_rows: list[dict] = []
    ambiguities: list[Ambiguity] = []
    new_cache: dict = {}
    n_cached = 0

    for pdir in patient_dirs:
        patient_id = pdir.name
        m = PATIENT_RE.match(patient_id)
        if not m:
            print(f"[warn] skipping unparseable patient dir: {pdir}", file=sys.stderr)
            continue
        cohort_id = "-".join(patient_id.split("-")[:3])

        # --- cache check: ~11 stat() calls, no readdir ---------------------- #
        key = str(pdir)
        sig = patient_signature(pdir)
        if cache is not None:
            hit = cache.get(key)
            if hit and hit.get("sig") == sig:
                manifest_rows.append(hit["row"])
                sample_rows.extend(hit["samples"])
                ambiguities.extend(Ambiguity(**a) for a in hit["ambig"])
                new_cache[key] = hit
                n_cached += 1
                continue

        listing = scan_patient(pdir)
        pat_ambig: list[Ambiguity] = []
        samples = discover_samples(listing)
        tumour_dna = [s for s in samples if s["seqtype"] == "DT"]
        normal_dna = [s for s in samples if s["seqtype"] == "DN"]
        tumour_rna = [s for s in samples if s["seqtype"] == "RT"]

        # Primary tumour = lowest aliquot index. Do NOT assume aliquot 1 exists
        # (MoHQ-CM-4-105 has a 2DT and no 1DT).
        primary_dt = tumour_dna[0]["sample_id"] if tumour_dna else ""
        primary_dn = normal_dna[0]["sample_id"] if normal_dna else ""
        primary_rt = tumour_rna[0]["sample_id"] if tumour_rna else ""

        row: dict[str, str] = {
            "collection": m.group("collection"),
            "institution": m.group("institution"),
            "cohort_id": cohort_id,
            "patient_id": patient_id,
            "patient_dir": str(pdir),
            "dna_normal_sample": primary_dn,
            "dna_tumour_sample": primary_dt,
            "rna_tumour_sample": primary_rt,
            "n_dna_tumour_samples": str(len(tumour_dna)),
            "n_dna_normal_samples": str(len(normal_dna)),
            "n_rna_tumour_samples": str(len(tumour_rna)),
            "multi_tumour": "yes" if len(tumour_dna) > 1 else "no",
            "all_samples": ";".join(s["sample_id"] for s in samples),
        }

        for asset in ASSETS:
            hit = resolve_asset(listing, pdir, patient_id, asset, pat_ambig)
            row[asset.name] = hit if hit else "NA"

        # Legacy promotion: some deliveries (observed: MoHQ-CM-4-49) predate the
        # _D/_R suffix convention. Promote the unsuffixed file into the somatic
        # slot ONLY when that slot is empty, and always record an ambiguity so
        # the choice is auditable rather than silent.
        for legacy, target in (("unsuffixed_maf", "somatic_maf"),
                               ("unsuffixed_tiers", "somatic_tiers"),
                               ("unsuffixed_cna", "cna_segments")):
            if row[legacy] != "NA" and row[target] == "NA":
                row[target] = row[legacy]
                pat_ambig.append(Ambiguity(
                    patient_id=patient_id, asset=target,
                    chosen=row[legacy], alternatives=[],
                    reason="LEGACY naming: no _D/_R suffix; assumed DNA. VERIFY MANUALLY.",
                ))

        # ------------------------------------------------------------------ #
        # Mutation input, with a RECORDED fallback.
        #
        # Older deliveries (those using the `_D.acmg.` MAF naming) shipped no
        # PCGR VCF at all -- 6/20 patients in the sample tree, i.e. ~30%. Since
        # the collection MAFs are unannotated, requiring the PCGR VCF would
        # silently drop every one of those patients from all mutation analyses.
        #
        # They do all have an ensemble somatic VCF, which vcf2maf can consume.
        # BUT THE TWO ARE NOT EQUIVALENT: the PCGR VCF is post-PCGR filtering
        # and tiering, while the ensemble VCF is the caller union. Variant
        # counts will differ systematically, which biases TMB and mutation
        # frequency -- and because the source tracks delivery version, that bias
        # correlates with time and possibly institution. Exactly the kind of
        # batch effect this pipeline is built to detect.
        #
        # So: use the fallback, but RECORD which source each patient used, and
        # never compare counts across sources without adjusting for it.
        # ------------------------------------------------------------------ #
        # pcgr_vcf_status describes what CAN be done, independently of what the
        # pipeline is configured to do:
        #   delivered   -- PCGR VCF shipped with the data
        #   regenerable -- no PCGR VCF, but an ensemble VCF exists, so RUN_PCGR
        #                  can rebuild one. NOT a reason to exclude the patient.
        #   unavailable -- neither; nothing can be done
        if row["pcgr_vcf"] != "NA":
            row["pcgr_vcf_status"] = "delivered"
            row["mutation_vcf"] = row["pcgr_vcf"]
            row["mutation_vcf_source"] = "pcgr_delivered"
        elif row["ensemble_somatic_vcf"] != "NA":
            row["pcgr_vcf_status"] = "regenerable"
            row["mutation_vcf"] = row["ensemble_somatic_vcf"]
            # The manifest is built BEFORE the pipeline runs, so it cannot know
            # whether gap-fill will be enabled. This records the input that is
            # available; main.nf writes the source ACTUALLY used to
            # mutation_provenance.tsv once routing has happened.
            row["mutation_vcf_source"] = "ensemble_pending_gapfill"
        else:
            row["pcgr_vcf_status"] = "unavailable"
            row["mutation_vcf"] = "NA"
            row["mutation_vcf_source"] = "none"

        # Processing round: PCGR / R / genome module versions from the ini.
        vers = read_ini_versions(row.get("tumourpair_ini", "NA"))
        row["pcgr_version"] = vers.get("pcgr", "unknown")
        row["genpipes_version"] = vers.get("genpipes") or vers.get("mugqic_pipelines", "unknown")
        row["r_version"] = vers.get("r", "unknown")

        # Genome build. The filename value is a CLAIM; a verified value from
        # check_genome_build.sh (header-derived) overrides it when available.
        claimed = detect_genome_build(
            [row["pcgr_vcf"], row["somatic_maf"], row["somatic_tiers"],
             row["cna_segments"]])
        row["genome_build_from_filename"] = claimed

        verified = (build_table or {}).get(patient_id)
        if verified:
            row["genome_build"] = verified
            row["genome_build_source"] = "verified_from_headers"
            if verified != claimed and claimed != "unknown":
                ambiguities.append(Ambiguity(
                    patient_id=patient_id, asset="genome_build",
                    chosen=verified, alternatives=[claimed],
                    reason=f"FILENAME CLAIMS {claimed} BUT HEADERS SAY {verified}. "
                           f"Headers win. Report this to the data providers.",
                ))
        else:
            row["genome_build"] = claimed
            row["genome_build_source"] = "filename_claim_UNVERIFIED"

        # A patient is analysis-ready if every `required` asset is present.
        missing_core = [a.name for a in ASSETS if a.required and row[a.name] == "NA"]
        if row["mutation_vcf"] == "NA":
            missing_core.append("mutation_vcf")
        row["analysis_ready"] = "yes" if not missing_core else "no"
        row["missing_core"] = ";".join(missing_core) if missing_core else "NA"

        manifest_rows.append(row)

        pat_samples = [{
            "cohort_id": cohort_id,
            "patient_id": patient_id,
            "sample_id": s["sample_id"],
            "sample_num": s["sample_num"],
            "aliquot": str(s["aliquot"]),
            "seqtype": s["seqtype"],
            "molecule": s["molecule"],
            "tissue": s["tissue"],
            "is_primary_tumour": "yes" if s["sample_id"] == primary_dt else "no",
        } for s in samples]
        sample_rows.extend(pat_samples)
        ambiguities.extend(pat_ambig)

        new_cache[key] = {
            "sig": sig,
            "row": row,
            "samples": pat_samples,
            "ambig": [a.__dict__ for a in pat_ambig],
        }

    if cache is not None and n_cached:
        print(f"[manifest] {n_cached} patient(s) unchanged since last run "
              f"(skipped the filesystem scan)")
    return manifest_rows, sample_rows, ambiguities, new_cache


# --------------------------------------------------------------------------- #
# Writers
# --------------------------------------------------------------------------- #

def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                           extrasaction="ignore", restval="NA")
        w.writeheader()
        w.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, type=Path,
                    help="Directory to walk (collection root, institution dir, or a single cohort dir)")
    ap.add_argument("--cohort-id", default=None,
                    help="Emit only this cohort (e.g. MoHQ-CM-4). Default: all found.")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--prefix", default=None,
                    help="Output filename prefix. Default: --cohort-id, else 'all_cohorts'.")
    ap.add_argument("--cache", type=Path, default=None,
                    help="JSON cache. Unchanged patients (by subdirectory mtime) are "
                         "skipped without reading their directories. Safe to delete.")
    ap.add_argument("--build-table", type=Path, default=None,
                    help="TSV from bin/check_genome_build.sh. Genome build read from "
                         "BAM/VCF headers OVERRIDES the filename claim, which the MoHQ "
                         "collection is known to get wrong.")
    ap.add_argument("--path-root", type=Path, default=None,
                    help="Prefix to record in the manifest instead of the resolved "
                         "--root. Only needed if --root is a COPY of the collection "
                         "rather than a symlink to it.")
    args = ap.parse_args(argv)

    if not args.root.is_dir():
        ap.error(f"--root is not a directory: {args.root}")

    # Record ABSOLUTE paths in the manifest.
    #
    # Every asset path is built as `root / patient / subdir / file`, so a relative
    # --root produces relative paths in the TSV. Nextflow then resolves those
    # against launchDir and you get, from a real run:
    #     No such file or directory: /home/mslavova/pipeline/MoHQ-HM-19/...
    # ...because Nextflow stages the cohort directory into the task work dir as a
    # symlink named after the cohort, so --root arrives as the bare name.
    #
    # .resolve() follows that symlink back to the real collection, which is what
    # we want recorded: the manifest should point at the data, not at a work
    # directory that -resume or a cleanup may remove.
    scan_root = args.root.resolve()
    emit_root = (args.path_root.resolve() if args.path_root else scan_root)
    args.root = scan_root

    if args.path_root is None and any(
            part in {"work", ".nextflow"} for part in scan_root.parts):
        print(f"[warn] --root resolves inside a work directory: {scan_root}\n"
              f"       Manifest paths will point there and may vanish on cleanup.\n"
              f"       Pass --path-root <real collection dir> if that is not intended.",
              file=sys.stderr)

    t0 = time.time()
    cache = None
    if args.cache:
        if args.cache.exists():
            try:
                cache = json.loads(args.cache.read_text())
                print(f"[manifest] loaded cache with {len(cache)} patient(s) "
                      f"from {args.cache}")
            except (json.JSONDecodeError, OSError) as e:
                print(f"[warn] cache unreadable ({e}); rebuilding from scratch",
                      file=sys.stderr)
                cache = {}
        else:
            cache = {}

    pdirs = find_patient_dirs(args.root)
    if args.cohort_id:
        pdirs = [p for p in pdirs
                 if "-".join(p.name.split("-")[:3]) == args.cohort_id]
    if not pdirs:
        print(f"[error] no patient directories found under {args.root}"
              + (f" for cohort {args.cohort_id}" if args.cohort_id else ""),
              file=sys.stderr)
        return 2

    build_table = None
    if args.build_table:
        if not args.build_table.exists():
            ap.error(f"--build-table not found: {args.build_table}")
        build_table = read_build_table(args.build_table)
        print(f"[manifest] loaded verified genome build for {len(build_table)} patient(s)")

    manifest_rows, sample_rows, ambiguities, new_cache = build_rows(
        pdirs, cache, build_table)

    if emit_root != scan_root:
        path_cols = ["patient_dir", "mutation_vcf"] + [a.name for a in ASSETS]
        old, new = str(scan_root), str(emit_root)
        for r in manifest_rows:
            for c in path_cols:
                if r.get(c, "NA").startswith(old):
                    r[c] = new + r[c][len(old):]
        for a in ambiguities:
            if a.chosen.startswith(old):
                a.chosen = new + a.chosen[len(old):]
            a.alternatives = [new + x[len(old):] if x.startswith(old) else x
                              for x in a.alternatives]
        print(f"[manifest] paths recorded under {new}")

    # A relative path in the manifest is a silent time bomb: the file exists when
    # this script runs and does not exist when Nextflow resolves it against
    # launchDir. Refuse to write one.
    bad = [(r["patient_id"], c, r[c])
           for r in manifest_rows
           for c in (["mutation_vcf"] + [a.name for a in ASSETS])
           if r.get(c, "NA") != "NA" and not r[c].startswith("/")]
    if bad:
        print(f"[error] {len(bad)} manifest path(s) are not absolute, e.g. "
              f"{bad[0][0]} {bad[0][1]}={bad[0][2]}", file=sys.stderr)
        return 3

    if args.cache:
        args.cache.parent.mkdir(parents=True, exist_ok=True)
        tmp = args.cache.with_suffix(args.cache.suffix + ".tmp")
        tmp.write_text(json.dumps(new_cache))
        tmp.replace(args.cache)        # atomic: never leave a half-written cache

    prefix = args.prefix or args.cohort_id or "all_cohorts"
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    manifest_fields = (
        ["collection", "institution", "cohort_id", "patient_id", "patient_dir",
         "dna_normal_sample", "dna_tumour_sample", "rna_tumour_sample",
         "n_dna_tumour_samples", "n_dna_normal_samples", "n_rna_tumour_samples",
         "multi_tumour", "analysis_ready", "missing_core",
         "genome_build", "genome_build_source", "genome_build_from_filename",
         "pcgr_version", "genpipes_version", "r_version",
         "pcgr_vcf_status", "mutation_vcf", "mutation_vcf_source", "all_samples"]
        + [a.name for a in ASSETS]
    )
    write_tsv(outdir / f"{prefix}.manifest.tsv", manifest_rows, manifest_fields)

    completeness_rows = []
    for r in manifest_rows:
        cr = {k: r[k] for k in ("cohort_id", "patient_id", "analysis_ready")}
        for a in ASSETS:
            cr[a.name] = "yes" if r[a.name] != "NA" else "no"
        completeness_rows.append(cr)
    write_tsv(outdir / f"{prefix}.completeness.tsv", completeness_rows,
              ["cohort_id", "patient_id", "analysis_ready"] + [a.name for a in ASSETS])

    write_tsv(outdir / f"{prefix}.ambiguities.tsv",
              [{"patient_id": a.patient_id, "asset": a.asset, "chosen": a.chosen,
                "n_alternatives": str(len(a.alternatives)),
                "alternatives": ";".join(a.alternatives) or "NA",
                "reason": a.reason} for a in ambiguities],
              ["patient_id", "asset", "chosen", "n_alternatives", "alternatives", "reason"])

    write_tsv(outdir / f"{prefix}.samples.tsv", sample_rows,
              ["cohort_id", "patient_id", "sample_id", "sample_num", "aliquot",
               "seqtype", "molecule", "tissue", "is_primary_tumour"])

    # --- console summary --------------------------------------------------- #
    # Kept aggregate: at 4,500 patients, listing individuals is unreadable.
    n = len(manifest_rows)
    ready = sum(1 for r in manifest_rows if r["analysis_ready"] == "yes")
    n_cohorts = len({r["cohort_id"] for r in manifest_rows})
    print(f"[manifest] {n} patients across {n_cohorts} cohort(s); "
          f"{ready}/{n} analysis-ready  [{time.time() - t0:.1f}s]")

    # Separate "some patients are missing this" from "no patient has it".
    # A 0/n asset almost always means it was never harvested, not that the
    # collection lacks it -- listing those together with genuine per-patient
    # gaps buries the signal in noise.
    absent_everywhere = []
    for a in ASSETS:
        have = sum(1 for r in manifest_rows if r[a.name] != "NA")
        if have == 0:
            absent_everywhere.append(a.name)
        elif have < n:
            flag = "  <-- REQUIRED" if a.required else ""
            print(f"  {a.name:<26} {have:>5}/{n}{flag}")

    if absent_everywhere:
        req = [x for x in absent_everywhere if ASSET_BY_NAME[x].required]
        print(f"  ({len(absent_everywhere)} asset(s) absent for ALL patients -- "
              f"normally means they were not harvested, not that they are missing "
              f"upstream:")
        print("     " + ", ".join(absent_everywhere[:12])
              + (" ..." if len(absent_everywhere) > 12 else "") + ")")
        if req:
            print(f"  !! {len(req)} of these are REQUIRED: {', '.join(req)}")
            print( "     Add the relevant --sets to harvest_juno.py and re-harvest.")
        if "tumourpair_ini" in absent_everywhere:
            print( "  !! tumourpair_ini absent -> pcgr_version cannot be determined.")
            print( "     Re-harvest to pick up parameters/*.ini (tiny; rclone skips")
            print( "     everything already downloaded).")

    if ambiguities:
        by_asset: dict[str, int] = {}
        for a in ambiguities:
            by_asset[a.asset] = by_asset.get(a.asset, 0) + 1
        print(f"[manifest] {len(ambiguities)} ambiguous asset(s) across "
              f"{len({a.patient_id for a in ambiguities})} patient(s); "
              f"see {prefix}.ambiguities.tsv")
        for asset, cnt in sorted(by_asset.items(), key=lambda kv: -kv[1])[:8]:
            print(f"  {asset:<26} {cnt:>5}")
        legacy = sum(1 for a in ambiguities if "LEGACY" in a.reason)
        if legacy:
            print(f"  !! {legacy} legacy unsuffixed file(s) assumed to be DNA -- VERIFY")

    # Genome build. Mixing GRCh37 and GRCh38 in one cohort silently corrupts
    # every coordinate-based analysis, so this is reported prominently.
    builds: dict[str, int] = {}
    for r in manifest_rows:
        builds[r["genome_build"]] = builds.get(r["genome_build"], 0) + 1
    n_verified = sum(1 for r in manifest_rows
                     if r["genome_build_source"] == "verified_from_headers")
    print(f"[manifest] genome build ({n_verified}/{n} verified from file headers):")
    for b, cnt in sorted(builds.items(), key=lambda kv: -kv[1]):
        print(f"  {b:<26} {cnt:>5}/{n}")
    if n_verified < n:
        print(f"  !! {n - n_verified} patient(s) rely on the FILENAME only, which this")
        print( "     collection is known to get wrong (files named grch38 whose")
        print( "     alignment was GRCh37). Verify with:")
        print( "       bash bin/check_genome_build.sh <INST/COHORT> 999 > builds.tsv")
        print( "       python3 bin/build_manifest.py ... --build-table builds.tsv")
    conflicts = [a for a in ambiguities if a.asset == "genome_build"]
    if conflicts:
        print(f"  !! {len(conflicts)} patient(s) where the filename CONTRADICTS the headers:")
        for a in conflicts[:5]:
            print(f"       {a.patient_id}: name says {a.alternatives[0]}, headers say {a.chosen}")
    real_builds = {b for b in builds if b not in ("unknown",) and not b.startswith("MIXED")}
    if len(real_builds) > 1 or any(b.startswith("MIXED") for b in builds):
        print("  !! MORE THAN ONE GENOME BUILD IN THIS SET.")
        print("     GRCh37 and GRCh38 coordinates are not interchangeable. A .seg file,")
        print("     CNV frequency plot or gene overlap built across both is wrong, and")
        print("     nothing will error. Split by `genome_build` and analyse separately,")
        print("     or lift over to a single build before running the pipeline.")

    # Processing round. Confirmed to vary WITHIN a cohort (MoHQ-CM-1 has both
    # PCGR 0.9.2 and 1.0.3), so this is a covariate, not a footnote.
    pv: dict[str, int] = {}
    for r in manifest_rows:
        pv[r["pcgr_version"]] = pv.get(r["pcgr_version"], 0) + 1
    if len(pv) > 1 or "unknown" not in pv:
        print("[manifest] PCGR version:")
        for v, cnt in sorted(pv.items(), key=lambda kv: -kv[1]):
            print(f"  {v:<26} {cnt:>5}/{n}")
    if len([v for v in pv if v != "unknown"]) > 1:
        print("  !! MORE THAN ONE PCGR VERSION -- the cohort was processed in")
        print("     separate rounds. Output columns and calling behaviour differ")
        print("     between PCGR majors, so differences between these groups are")
        print("     TECHNICAL. Adjust for `pcgr_version` or restrict to one round.")

        # Which versions do the patients NEEDING regeneration belong to? This
        # decides what RUN_PCGR should load for each of them -- regenerating a
        # 1.4.1 patient with 1.0.3 would introduce a difference that was not
        # there before.
        cross: dict[tuple[str, str], int] = {}
        for r in manifest_rows:
            cross[(r["pcgr_vcf_status"], r["pcgr_version"])] = \
                cross.get((r["pcgr_vcf_status"], r["pcgr_version"]), 0) + 1
        print("\n  PCGR version x VCF availability (drives the gap-fill choice):")
        print(f"    {'status':<14}{'pcgr version':<28}{'n':>4}")
        for (status, ver), cnt in sorted(cross.items()):
            print(f"    {status:<14}{ver:<28}{cnt:>4}")
        regen_vers = {v for (s, v) in cross if s == "regenerable" and v != "unknown"}
        if len(regen_vers) > 1:
            print(f"    -> patients needing regeneration span {len(regen_vers)} versions.")
            print( "       Leave params.pcgr_module unset so RUN_PCGR uses each patient's")
            print( "       OWN recorded version, rather than forcing one for all.")

    # PCGR VCF availability. Patients without one are NOT excluded -- they are
    # gap-fill candidates.
    st: dict[str, int] = {}
    for r in manifest_rows:
        st[r["pcgr_vcf_status"]] = st.get(r["pcgr_vcf_status"], 0) + 1
    print("[manifest] PCGR VCF availability:")
    for k, cnt in sorted(st.items(), key=lambda kv: -kv[1]):
        label = {"delivered":   "delivered with the data",
                 "regenerable": "MISSING -> RUN_PCGR can rebuild it (not excluded)",
                 "unavailable": "no PCGR and no ensemble VCF -> excluded"}.get(k, k)
        print(f"  {k:<14} {cnt:>5}/{n}   {label}")

    if st.get("regenerable"):
        print(f"\n  {st['regenerable']} patient(s) need their PCGR VCF regenerated.")
        print( "    run_pcgr_gapfill = true   -> RUN_PCGR rebuilds them from the ensemble")
        print( "                                 VCF; every patient then follows one path")
        print( "    run_pcgr_gapfill = false  -> they use the ensemble VCF directly. That is")
        print( "                                 the caller UNION, not PCGR-filtered calls, so")
        print( "                                 their variant counts are not comparable with")
        print( "                                 the delivered patients. Adjust for")
        print( "                                 mutation_vcf_source, or enable gap-fill.")
        print( "    Either way these patients ARE analysed -- see mutation_provenance.tsv")
        print( "    in the results for the source actually used per patient.")

    multi = [r["patient_id"] for r in manifest_rows if r["multi_tumour"] == "yes"]
    if multi:
        shown = ", ".join(multi[:5]) + (f", +{len(multi) - 5} more" if len(multi) > 5 else "")
        print(f"[manifest] {len(multi)} patient(s) with >1 tumour sample; only the "
              f"lowest-aliquot tumour is used downstream ({shown})")

    not_ready = [r for r in manifest_rows if r["analysis_ready"] != "yes"]
    if not_ready:
        reasons: dict[str, int] = {}
        for r in not_ready:
            reasons[r["missing_core"]] = reasons.get(r["missing_core"], 0) + 1
        print(f"[manifest] {len(not_ready)} patient(s) excluded from analyses:")
        for reason, cnt in sorted(reasons.items(), key=lambda kv: -kv[1])[:8]:
            print(f"  missing {reason:<40} {cnt:>5} patient(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

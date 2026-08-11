#!/usr/bin/env python3
"""
harvest_juno.py -- pull ONLY the cohort-analysis files out of the Juno object
                   store, preserving the collection's directory layout.

Run this on Cardinal (Juno is reachable from there). It writes a tree that
build_manifest.py can walk unchanged, because the layout is preserved exactly:

    <dest>/CM/MoHQ-CM-3/MoHQ-CM-3-1/reports/pcgr/...
                                    /variants/...
                                    /raw_cnv/...

WHY HARVEST RATHER THAN MOUNT OR STREAM
---------------------------------------
Measured against the 20-patient sample tree: 99.7% of the collection is BAM and
FASTQ, which no cohort-level analysis reads. The files this pipeline actually
consumes are ~98 MB per patient:

    full collection   ~2.7 PB   (4,500 patients x ~607 GB)
    harvest           ~431 GB   (4,500 patients x ~98 MB)
    reduction         ~6,300x

So a one-time selective copy is small, and after that every analysis is local,
resumable and fast.

IMPROVEMENTS OVER THE PER-FILE `rclone cat | ssh` LOOP
------------------------------------------------------
The original transfer script did, per file:
    1. one SSH round trip to test -f the destination
    2. one `rclone cat` piped through a second SSH
That is two SSH connections per file, strictly serial. At ~10 harvested files
per patient x 4,500 patients that is ~90,000 sequential SSH round trips.

Here instead:
    * ONE recursive listing per cohort (cached to disk, reused across runs)
    * pattern matching done locally against that listing
    * ONE `rclone copy --files-from` per cohort, which transfers in parallel
      (--transfers) and does its own size/mtime skip logic, so resuming is
      free and no destination-existence probing is needed
    * if the destination is on a different host, an rclone SFTP remote replaces
      the ssh pipe entirely, keeping the parallelism

Requires: rclone on PATH. No Python dependencies beyond the stdlib.
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import os
import subprocess
import sys
import time
from pathlib import Path

# --------------------------------------------------------------------------- #
# What to harvest.
#
# Deliberately a SEPARATE, simpler list from build_manifest.ASSETS: these are
# glob patterns matched against object keys relative to the cohort prefix, and
# they must stay readable because they decide what lands on scratch.
#
# NOTE the MAF situation: the collection's *_D.*.maf files are NOT gene
# annotated (Hugo_Symbol = "Unknown"), so the *VCF* is the real input and gets
# converted locally by the VCF2MAF pipeline step. The unannotated MAF is
# harvested anyway -- it is small, and keeping it makes the provenance auditable.
# --------------------------------------------------------------------------- #
HARVEST_SETS: dict[str, tuple[str, ...]] = {
    # Everything the cohort pipeline reads. ~98 MB/patient.
    "core": (
        "*/reports/pcgr/*_D*.vcf.gz",              # THE mutation input -> vcf2maf
        "*/reports/pcgr/*_D*.vcf.gz.tbi",
        "*/reports/pcgr/*_D*.maf",                 # unannotated; provenance only
        "*/reports/pcgr/*_D*snvs_indels.tiers.tsv",
        "*/reports/pcgr/*_D*cna_segments.tsv.gz",
        "*/raw_cnv/*.cnvkit.vcf.gz",
        "*/raw_cnv/*.cnvkit.vcf.gz.tbi",
        "*/expression/*.abundance_genes.tsv",
        "*/reports/*.anno_fuse.tsv",
        "*/svariants/linx/*.linx.fusion.tsv",
        "*/svariants/linx/*.linx.driver.catalog.tsv",
        "*/svariants/*.driver.catalog.somatic.tsv",
        "*/Key_metrics.csv",
        "*/log.csv",
        "*/Warnings.html",
        # GenPipes configs: tiny (~46 KB) but they are the ONLY record of which
        # PCGR / R / genome versions produced each patient. Without them the
        # manifest cannot populate pcgr_version, and processing-round becomes an
        # invisible batch variable.
        "*/parameters/*.ini",
    ),
    # Adds the ensemble somatic VCFs needed for mutational signatures.
    # Roughly doubles the footprint (~50 MB/patient more).
    "signatures": (
        "*/variants/*.ensemble.somatic.vt.annot.vcf.gz",
        "*/variants/*.ensemble.somatic.vt.annot.vcf.gz.tbi",
    ),
    # Structural-variant detail for LINX/GRIPSS cohort work.
    "sv": (
        "*/svariants/*.gripss.filtered.somatic.vcf.gz",
        "*/svariants/linx/*.linx.svs.tsv",
        "*/svariants/linx/*.linx.clusters.tsv",
        "*/svariants/linx/*.linx.breakend.tsv",
    ),
    # Germline predisposition.
    "germline": (
        "*/reports/*_D.cpsr.zip",
        "*/variants/*.ensemble.germline.vt.annot.vcf.gz",
    ),
    # Expression at transcript level.
    "transcripts": (
        "*/expression/*.abundance_transcripts.tsv",
    ),
}


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, check=False, text=True, **kw)


def rclone_available() -> None:
    p = run(["rclone", "version"], capture_output=True)
    if p.returncode != 0:
        sys.exit("rclone not found on PATH. Load it or install per "
                 "https://rclone.org/install/")
    print(f"[harvest] {p.stdout.splitlines()[0]}")


def list_cohort(remote: str, cohort_path: str, cache: Path | None,
                refresh: bool) -> list[dict]:
    """One recursive listing of a cohort prefix. Cached to disk.

    Uses `rclone lsjson --recursive --files-only`, which returns Path/Size/
    ModTime in a single paginated LIST sequence -- far cheaper than probing
    each expected file, and it also gives the sizes needed to report the
    transfer volume before committing to it.
    """
    if cache and cache.exists() and not refresh:
        age_h = (time.time() - cache.stat().st_mtime) / 3600
        print(f"[harvest] reusing cached listing ({age_h:.1f} h old): {cache}")
        return json.loads(cache.read_text())

    src = f"{remote}/{cohort_path}" if cohort_path else remote
    print(f"[harvest] listing {src} ...", flush=True)
    t0 = time.time()
    # --fast-list issues far fewer LIST requests on S3-compatible backends at
    # the cost of holding the listing in memory, which is the right trade here.
    p = run(["rclone", "lsjson", "--recursive", "--files-only", "--fast-list", src],
            capture_output=True)
    if p.returncode != 0:
        sys.exit(f"rclone lsjson failed:\n{p.stderr}")
    entries = json.loads(p.stdout or "[]")
    print(f"[harvest] {len(entries):,} objects listed in {time.time() - t0:.1f}s")

    if cache:
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(entries))
    return entries


def select(entries: list[dict], patterns: tuple[str, ...]) -> list[dict]:
    out, seen = [], set()
    for e in entries:
        path = e.get("Path", "")
        if path in seen:
            continue
        for pat in patterns:
            if fnmatch.fnmatch(path, pat):
                out.append(e)
                seen.add(path)
                break
    return out


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PB"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--remote", required=True,
                    help="e.g. juno:d5f8b8e8e3e2442f81573b2f0951013b:MOH-Q")
    ap.add_argument("--cohort", required=True,
                    help="Path under the remote, e.g. CM/MoHQ-CM-3")
    ap.add_argument("--dest", required=True,
                    help="Destination root. Local path, or an rclone remote "
                         "such as narval-sftp:/scratch/slavova/MoHQ")
    ap.add_argument("--sets", nargs="+", default=["core"],
                    choices=sorted(HARVEST_SETS) + ["all"],
                    help="Which harvest sets to pull (default: core)")
    ap.add_argument("--transfers", type=int, default=16)
    ap.add_argument("--checkers", type=int, default=16)
    ap.add_argument("--cache-dir", type=Path, default=Path(".harvest_cache"))
    ap.add_argument("--refresh-listing", action="store_true")
    ap.add_argument("--dry-run", action="store_true",
                    help="List and size the selection without transferring")
    ap.add_argument("--verify", action="store_true",
                    help="After transfer, list the destination and confirm every "
                         "selected object arrived at the right size. Strongly "
                         "recommended for cross-cluster (Cardinal->Narval) copies.")
    ap.add_argument("--extra-rclone", nargs="*", default=[],
                    help="Extra flags passed through to rclone copy")
    args = ap.parse_args(argv)

    rclone_available()

    sets = sorted(HARVEST_SETS) if "all" in args.sets else args.sets
    patterns: tuple[str, ...] = tuple(p for s in sets for p in HARVEST_SETS[s])
    print(f"[harvest] sets: {', '.join(sets)}  ({len(patterns)} patterns)")

    cohort_name = args.cohort.rstrip("/").split("/")[-1]
    cache = args.cache_dir / f"{cohort_name}.lsjson.json"
    entries = list_cohort(args.remote, args.cohort, cache, args.refresh_listing)
    if not entries:
        print(f"[harvest] nothing listed under {args.cohort}", file=sys.stderr)
        return 2

    chosen = select(entries, patterns)
    total_all = sum(e.get("Size", 0) for e in entries)
    total_sel = sum(e.get("Size", 0) for e in chosen)
    patients = {e["Path"].split("/")[0] for e in chosen if "/" in e["Path"]}

    print(f"[harvest] {cohort_name}: {len(patients)} patients")
    print(f"[harvest] selected {len(chosen):,} / {len(entries):,} objects")
    print(f"[harvest] volume    {human(total_sel)} / {human(total_all)}"
          f"  ({100 * total_sel / max(total_all, 1):.3f}%"
          f", {total_all / max(total_sel, 1):.0f}x smaller)")

    # Per-asset breakdown: makes it obvious when a pattern silently matches
    # nothing because of an upstream naming change.
    # Index files (.tbi/.bai) are inherently optional -- nothing in the pipeline
    # does random access, so their absence is normal and should not read as an
    # error next to a genuinely missing data file.
    print("[harvest] by pattern:")
    for pat in patterns:
        hits = [e for e in chosen if fnmatch.fnmatch(e["Path"], pat)]
        optional = pat.endswith((".tbi", ".bai", ".md5"))
        if hits:
            flag = ""
        elif optional:
            flag = "   (none delivered -- fine, indexes are not used)"
        else:
            flag = "   <-- MATCHED NOTHING"
        print(f"    {pat:<52} {len(hits):>5} files"
              f" {human(sum(h.get('Size', 0) for h in hits)):>10}{flag}")

    if args.dry_run:
        print("[harvest] --dry-run: stopping before transfer")
        return 0

    # rclone copy --files-from: one parallel transfer, with rclone's own
    # size/mtime skip logic providing resume for free.
    listfile = args.cache_dir / f"{cohort_name}.files-from.txt"
    listfile.parent.mkdir(parents=True, exist_ok=True)
    listfile.write_text("\n".join(e["Path"] for e in chosen) + "\n")

    src = f"{args.remote}/{args.cohort}"
    dst = f"{args.dest.rstrip('/')}/{args.cohort}"
    cmd = ["rclone", "copy", src, dst,
           "--files-from", str(listfile),
           "--transfers", str(args.transfers),
           "--checkers", str(args.checkers),
           "--retries", "5",
           "--low-level-retries", "20",
           "--stats", "30s",
           "--stats-one-line",
           "--progress"] + args.extra_rclone

    print(f"[harvest] {' '.join(cmd)}", flush=True)
    t0 = time.time()
    p = run(cmd)
    dt = time.time() - t0
    if p.returncode != 0:
        print(f"[harvest] rclone copy exited {p.returncode}", file=sys.stderr)
        return p.returncode

    print(f"[harvest] {cohort_name} done in {dt / 60:.1f} min "
          f"({human(total_sel / max(dt, 1))}/s)")

    # ----------------------------------------------------------------- #
    # Verification.
    #
    # rclone compares size and modtime as it goes, but a cross-cluster SFTP
    # copy has more ways to end up with a short file than a local one, and a
    # truncated VCF or MAF fails LATER and confusingly -- halfway through VEP,
    # or as a mysteriously small variant count. Confirming the landing state
    # explicitly is cheap next to that.
    # ----------------------------------------------------------------- #
    if args.verify:
        print("[verify] listing destination ...", flush=True)
        vp = run(["rclone", "lsjson", "--recursive", "--files-only", "--fast-list", dst],
                 capture_output=True)
        if vp.returncode != 0:
            print(f"[verify] could not list destination:\n{vp.stderr}", file=sys.stderr)
            return 1
        got = {e["Path"]: e.get("Size", 0) for e in json.loads(vp.stdout or "[]")}

        missing, wrong = [], []
        for e in chosen:
            p, want = e["Path"], e.get("Size", 0)
            if p not in got:
                missing.append(p)
            elif got[p] != want:
                wrong.append((p, want, got[p]))

        receipt = args.cache_dir / f"{cohort_name}.receipt.tsv"
        with receipt.open("w") as fh:
            fh.write("path\tsource_size\tdest_size\tstatus\n")
            for e in chosen:
                p, want = e["Path"], e.get("Size", 0)
                have = got.get(p)
                status = ("MISSING" if have is None
                          else "SIZE_MISMATCH" if have != want else "ok")
                fh.write(f"{p}\t{want}\t{have if have is not None else 'NA'}\t{status}\n")

        if missing or wrong:
            print(f"[verify] FAILED: {len(missing)} missing, {len(wrong)} wrong size",
                  file=sys.stderr)
            for p in missing[:5]:
                print(f"  missing        {p}", file=sys.stderr)
            for p, w, g in wrong[:5]:
                print(f"  size mismatch  {p}  expected {w}, got {g}", file=sys.stderr)
            print(f"[verify] full receipt: {receipt}", file=sys.stderr)
            print("[verify] re-run the same command; rclone re-copies only the bad files.",
                  file=sys.stderr)
            return 1

        print(f"[verify] all {len(chosen):,} objects present at the expected size")
        print(f"[verify] receipt: {receipt}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

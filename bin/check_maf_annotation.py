#!/usr/bin/env python3
"""
Is the MAF gene-annotated, or is Hugo_Symbol "Unknown" throughout?

The MoHQ delivered MAFs were produced by a GenPipes version that did not run the
annotation step, so Hugo_Symbol is "Unknown" for every row. Any oncoplot, gene
summary or mutual-exclusivity analysis built on them is empty or meaningless --
and, critically, produces NO ERROR. This script measures how far that goes.

WHAT COUNTS AS ANNOTATED
------------------------
Not the fraction of "Unknown". These are whole-genome somatic calls, and most
somatic variants in WGS are intergenic, where "Unknown" is the CORRECT value. A
healthy WGS MAF here runs ~50-60% Unknown.

The signal is the number of DISTINCT gene symbols:

    annotated      thousands of distinct symbols
    unannotated    zero, or a handful

So that is what this reports, with the fraction alongside for context.

USAGE
-----
Local (already-harvested tree):

    python3 bin/check_maf_annotation.py --root /home/you/MoHQ/collection

Remote, WITHOUT downloading -- streams each MAF over rclone and discards it:

    python3 bin/check_maf_annotation.py \\
        --remote juno:PROJECT:MOH-Q --sample 3

    python3 bin/check_maf_annotation.py \\
        --remote juno:PROJECT:MOH-Q --cohorts HM/MoHQ-HM-19 CM/MoHQ-CM-1 --sample 5

--sample N looks at N patients per cohort, which is plenty to establish whether
a cohort is uniformly affected. Use --sample 0 for every patient.

RNA MAFs (_R) are reported separately and never mixed with somatic (_D) ones.
"""

from __future__ import annotations

import argparse
import gzip
import io
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

# vcf2maf/PCGR write Hugo_Symbol first. Verified against the delivered files
# rather than assumed -- if the header ever moves, we locate it by name below.
UNKNOWN = {"Unknown", "", "NA", "."}

PATIENT_RE = re.compile(r"(MoHQ-[A-Za-z]+-[A-Za-z0-9]+-[A-Za-z0-9]+)")
COHORT_RE = re.compile(r"(MoHQ-[A-Za-z]+-[A-Za-z0-9]+)")


class MafStats:
    __slots__ = ("path", "rows", "unknown", "genes", "molecule", "error")

    def __init__(self, path, molecule):
        self.path = path
        self.molecule = molecule
        self.rows = 0
        self.unknown = 0
        self.genes: set[str] = set()
        self.error: str | None = None

    @property
    def frac(self) -> float:
        return self.unknown / self.rows if self.rows else 0.0

    @property
    def verdict(self) -> str:
        if self.error:
            return "ERROR"
        if self.rows == 0:
            return "EMPTY"
        if len(self.genes) == 0:
            return "UNANNOTATED"
        if len(self.genes) < 50:
            return "SUSPECT"
        return "annotated"


def scan_stream(fh, stats: MafStats, max_rows: int = 0) -> None:
    """Read a MAF line by line. Never loads the file into memory."""
    col = None
    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if col is None:
            # Locate Hugo_Symbol by NAME, not position. Cheap, and it means a
            # column reorder shows up as an error rather than as nonsense.
            if "Hugo_Symbol" in fields:
                col = fields.index("Hugo_Symbol")
                continue
            stats.error = "no Hugo_Symbol column in header"
            return
        if col >= len(fields):
            continue
        sym = fields[col]
        stats.rows += 1
        if sym in UNKNOWN:
            stats.unknown += 1
        else:
            stats.genes.add(sym)
        if max_rows and stats.rows >= max_rows:
            return


def scan_local(path: Path, max_rows: int) -> MafStats:
    molecule = "RNA" if re.search(r"_R\b|_R\.", path.name) else "DNA"
    st = MafStats(str(path), molecule)
    try:
        if path.suffix == ".gz":
            with gzip.open(path, "rt", errors="replace") as fh:
                scan_stream(fh, st, max_rows)
        else:
            with open(path, "r", errors="replace") as fh:
                scan_stream(fh, st, max_rows)
    except OSError as e:
        st.error = str(e)
    return st


def scan_remote(remote_path: str, max_rows: int) -> MafStats:
    """Stream via `rclone cat` -- nothing is written to disk."""
    molecule = "RNA" if re.search(r"_R\.", remote_path) else "DNA"
    st = MafStats(remote_path, molecule)
    proc = subprocess.Popen(["rclone", "cat", remote_path],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        stream = io.TextIOWrapper(proc.stdout, errors="replace")
        if remote_path.endswith(".gz"):
            stream = io.TextIOWrapper(gzip.GzipFile(fileobj=proc.stdout),
                                      errors="replace")
        scan_stream(stream, st, max_rows)
    except Exception as e:                       # noqa: BLE001 - report, don't crash
        st.error = f"{type(e).__name__}: {e}"
    finally:
        # We usually stop early; killing rclone avoids transferring the rest.
        if proc.poll() is None:
            proc.kill()
        proc.wait()
    return st


def key_of(path: str) -> tuple[str, str]:
    c = COHORT_RE.search(path)
    p = PATIENT_RE.search(path)
    return (c.group(1) if c else "?", p.group(1) if p else "?")


def rclone_list(remote: str, sub: str) -> list[str]:
    cmd = ["rclone", "lsf", "--recursive", "--files-only",
           "--include", "*.maf", "--include", "*.maf.gz",
           f"{remote}/{sub}" if sub else remote]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        print(f"[warn] rclone lsf failed for {sub}: {out.stderr.strip()[:200]}",
              file=sys.stderr)
        return []
    base = f"{remote}/{sub}" if sub else remote
    return [f"{base}/{line}" for line in out.stdout.splitlines() if line.strip()]


def report(results: list[MafStats], show_all: bool) -> int:
    by_cohort: dict[str, list[MafStats]] = defaultdict(list)
    for st in results:
        by_cohort[key_of(st.path)[0]].append(st)

    print()
    print(f"{'cohort':<18}{'patient':<20}{'mol':<5}{'rows':>9}{'unknown':>9}"
          f"{'frac':>7}{'genes':>8}  verdict")
    print("-" * 88)

    for cohort in sorted(by_cohort):
        for st in sorted(by_cohort[cohort], key=lambda s: s.path):
            if not show_all and st.verdict == "annotated" and st.molecule == "RNA":
                continue
            _, patient = key_of(st.path)
            print(f"{cohort:<18}{patient:<20}{st.molecule:<5}{st.rows:>9,}"
                  f"{st.unknown:>9,}{st.frac:>7.2f}{len(st.genes):>8,}  {st.verdict}"
                  + (f"  [{st.error}]" if st.error else ""))

    print()
    print("=" * 88)
    print("SUMMARY BY COHORT (DNA/somatic MAFs only)")
    print("=" * 88)
    print(f"{'cohort':<18}{'files':>7}{'annotated':>11}{'unannotated':>13}"
          f"{'suspect':>9}{'median genes':>14}")
    print("-" * 88)

    n_unann = 0
    for cohort in sorted(by_cohort):
        dna = [s for s in by_cohort[cohort] if s.molecule == "DNA"]
        if not dna:
            continue
        ann = sum(1 for s in dna if s.verdict == "annotated")
        una = sum(1 for s in dna if s.verdict == "UNANNOTATED")
        sus = sum(1 for s in dna if s.verdict == "SUSPECT")
        n_unann += una
        counts = sorted(len(s.genes) for s in dna)
        med = counts[len(counts) // 2] if counts else 0
        print(f"{cohort:<18}{len(dna):>7}{ann:>11}{una:>13}{sus:>9}{med:>14,}")

    print()
    if n_unann:
        print(f"{n_unann} somatic MAF(s) have ZERO distinct gene symbols.")
        print("Those cannot support any gene-level analysis. Regenerate them from")
        print("the PCGR VCFs (the pipeline's VCF2MAF step does this).")
    else:
        print("Every somatic MAF scanned carries real gene symbols.")
    print()
    print("Reminder: a high 'frac' is NORMAL for whole-genome data -- most somatic")
    print("variants are intergenic. Judge by 'genes', not by 'frac'.")
    return 1 if n_unann else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--root", type=Path, help="Local harvested collection")
    src.add_argument("--remote", help="rclone remote, e.g. juno:PROJECT:MOH-Q")
    ap.add_argument("--cohorts", nargs="*", default=None,
                    help="Remote mode: cohort subpaths, e.g. HM/MoHQ-HM-19. "
                         "Default: every cohort found.")
    ap.add_argument("--sample", type=int, default=3,
                    help="Patients per cohort (0 = all). Default 3.")
    ap.add_argument("--max-rows", type=int, default=200000,
                    help="Stop after N data rows per file (0 = whole file). "
                         "200k is far more than enough to count distinct genes.")
    ap.add_argument("--all", action="store_true",
                    help="Also list annotated RNA MAFs in the per-file table.")
    args = ap.parse_args(argv)

    results: list[MafStats] = []

    if args.root:
        if not args.root.is_dir():
            ap.error(f"--root is not a directory: {args.root}")
        mafs = sorted(list(args.root.rglob("*.maf")) + list(args.root.rglob("*.maf.gz")))
        if not mafs:
            print(f"[error] no .maf files under {args.root}", file=sys.stderr)
            return 2
        by_cohort: dict[str, list[Path]] = defaultdict(list)
        for m in mafs:
            by_cohort[key_of(str(m))[0]].append(m)
        for cohort, files in sorted(by_cohort.items()):
            per_patient: dict[str, list[Path]] = defaultdict(list)
            for f in files:
                per_patient[key_of(str(f))[1]].append(f)
            chosen = sorted(per_patient)
            if args.sample:
                chosen = chosen[:args.sample]
            print(f"[scan] {cohort}: {len(files)} MAF(s), "
                  f"scanning {len(chosen)} patient(s)", file=sys.stderr)
            for p in chosen:
                for f in sorted(per_patient[p]):
                    results.append(scan_local(f, args.max_rows))
    else:
        cohorts = args.cohorts
        if not cohorts:
            out = subprocess.run(["rclone", "lsf", "--dirs-only", args.remote],
                                 capture_output=True, text=True)
            insts = [d.strip("/") for d in out.stdout.split()]
            cohorts = []
            for inst in insts:
                sub = subprocess.run(["rclone", "lsf", "--dirs-only",
                                      f"{args.remote}/{inst}"],
                                     capture_output=True, text=True)
                cohorts += [f"{inst}/{c.strip('/')}" for c in sub.stdout.split()]
            print(f"[scan] found {len(cohorts)} cohort(s)", file=sys.stderr)

        for sub in cohorts:
            files = rclone_list(args.remote, sub)
            per_patient: dict[str, list[str]] = defaultdict(list)
            for f in files:
                per_patient[key_of(f)[1]].append(f)
            chosen = sorted(per_patient)[:args.sample] if args.sample else sorted(per_patient)
            print(f"[scan] {sub}: {len(files)} MAF(s), sampling {len(chosen)} patient(s)",
                  file=sys.stderr)
            for p in chosen:
                for f in per_patient[p]:
                    results.append(scan_remote(f, args.max_rows))

    if not results:
        print("[error] nothing scanned", file=sys.stderr)
        return 2
    return report(results, args.all)


if __name__ == "__main__":
    sys.exit(main())

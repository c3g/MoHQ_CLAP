#!/bin/bash
# ===========================================================================
# run_tests.sh -- pre-flight checks. Run on Narval before touching real data.
#
#   ./tests/run_tests.sh                 # uses ./pipeline.sif
#   ./tests/run_tests.sh /path/to.sif
#
# Three layers, cheapest first:
#   1. Python  -- manifest builder, on a synthetic MoHQ tree
#   2. R       -- shared library unit tests, ~5 s
#   3. Nextflow-- stub run, exercises the whole DAG without touching data
#
# If all three pass, the remaining risk is in your real file formats
# (PCGR column names, CNVkit INFO keys), which only real data can settle.
# ===========================================================================
set -uo pipefail

SIF="${1:-pipeline.sif}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"

pass=0; fail=0
hdr()  { printf '\n\033[1m=== %s ===\033[0m\n' "$1"; }
good() { printf '  \033[32mOK\033[0m    %s\n' "$1"; pass=$((pass+1)); }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; fail=$((fail+1)); }

# --------------------------------------------------------------------------
hdr "1. Python: manifest builder"

if python3 -c "import sys; sys.path.insert(0,'bin'); import build_manifest" 2>/dev/null; then
    good "build_manifest.py imports (all asset patterns validated)"
else
    bad "build_manifest.py failed to import"
    python3 -c "import sys; sys.path.insert(0,'bin'); import build_manifest"
fi

# Build a small synthetic tree covering the naming quirks that exist in the
# real collection, and check the manifest resolves all of them.
TREE=$(mktemp -d)
trap 'rm -rf "$TREE"' EXIT
python3 - "$TREE" <<'PY'
import os, sys
root = sys.argv[1]
# Three PCGR naming conventions + an alphanumeric sample number + a 2DT-only
# patient -- every quirk found in the real tree.
cases = [
    ("CM/MoHQ-CM-4/MoHQ-CM-4-10",  "MoHQ-CM-4-10",  "_D.acmg.grch38",      "157954-1DT"),
    ("CM/MoHQ-CM-4/MoHQ-CM-4-105", "MoHQ-CM-4-105", "_D.pcgr_acmg.grch38", "305329-2DT"),
    ("CQ/MoHQ-CQ-34/MoHQ-CQ-34-01","MoHQ-CQ-34-01", "_D.pcgr_acmg.grch38", "RCC01-1DT"),
]
for path, pid, conv, dt in cases:
    base = os.path.join(root, path)
    for d in ("reports/pcgr","variants","raw_cnv","expression","alignment"):
        os.makedirs(os.path.join(base,d), exist_ok=True)
    for f in [f"reports/pcgr/{pid}{conv}.vcf.gz",
              f"reports/pcgr/{pid}{conv}.maf",
              f"reports/pcgr/{pid}{conv}.snvs_indels.tiers.tsv",
              f"reports/pcgr/{pid}{conv}.cna_segments.tsv.gz",
              f"variants/{pid}.ensemble.somatic.vt.annot.vcf.gz",
              f"raw_cnv/{pid}.cnvkit.vcf.gz",
              "Key_metrics.csv",
              f"alignment/{pid}-{dt}.bam"]:
        open(os.path.join(base,f),"w").close()
PY

OUT=$(mktemp -d)
if python3 bin/build_manifest.py --root "$TREE" --outdir "$OUT" --prefix t >/dev/null 2>&1; then
    n_ready=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $h["analysis_ready"]=="yes"' "$OUT/t.manifest.tsv" | wc -l)
    [ "$n_ready" -eq 3 ] && good "3/3 patients resolved across all PCGR naming conventions" \
                         || bad "only $n_ready/3 patients analysis-ready"
    grep -q "RCC01" "$OUT/t.samples.tsv" && good "alphanumeric sample IDs (RCC01) parsed" \
                                         || bad "alphanumeric sample IDs not parsed"
    grep -q "2DT" "$OUT/t.samples.tsv" && good "aliquot-2-only tumour handled" \
                                       || bad "2DT sample not found"
else
    bad "build_manifest.py errored on the synthetic tree"
fi
rm -rf "$OUT"

# Cache correctness: a removed file must change the manifest.
CACHE=$(mktemp); OUT=$(mktemp -d)
python3 bin/build_manifest.py --root "$TREE" --outdir "$OUT" --prefix a --cache "$CACHE" >/dev/null 2>&1
rm -f "$TREE"/CM/MoHQ-CM-4/MoHQ-CM-4-10/raw_cnv/*.cnvkit.vcf.gz
python3 bin/build_manifest.py --root "$TREE" --outdir "$OUT" --prefix b --cache "$CACHE" >/dev/null 2>&1
if ! diff -q "$OUT/a.manifest.tsv" "$OUT/b.manifest.tsv" >/dev/null 2>&1; then
    good "manifest cache detects a removed file"
else
    bad "manifest cache did NOT notice a removed file"
fi
rm -rf "$OUT" "$CACHE"

# Paths in the manifest must be ABSOLUTE, and must survive the way Nextflow
# actually calls this: it stages the cohort directory into the task work dir as
# a symlink named after the cohort, so --root arrives as a bare relative name.
# A relative path here passes every check inside the task and then fails in the
# workflow, where file() resolves it against launchDir instead:
#     No such file or directory: <launchDir>/MoHQ-HM-19/MoHQ-HM-19-10/variants/...
STAGE=$(mktemp -d); OUT=$(mktemp -d)
ln -sfn "$TREE/CM/MoHQ-CM-4" "$STAGE/MoHQ-CM-4"
if ( cd "$STAGE" && python3 "$OLDPWD/bin/build_manifest.py" \
        --root MoHQ-CM-4 --cohort-id MoHQ-CM-4 --outdir "$OUT" --prefix s ) >/dev/null 2>&1
then
    n_rel=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
                        {v=$h["pcgr_vcf"]; if (v!="NA" && v !~ /^\//) c++}
                        END{print c+0}' "$OUT/s.manifest.tsv")
    n_ok=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
                       {v=$h["pcgr_vcf"]; if (v ~ /^\//) c++} END{print c+0}' "$OUT/s.manifest.tsv")
    if [ "$n_rel" -eq 0 ] && [ "$n_ok" -gt 0 ]; then
        good "manifest records absolute paths through a staged symlink"
    else
        bad "manifest has $n_rel relative path(s) - will break in the workflow"
    fi
else
    bad "build_manifest.py errored on a staged-symlink root"
fi
rm -rf "$STAGE" "$OUT"

# argparse builds a PYTHON command string from these R calls, and it has no
# representation for an empty vector: `default = character()` renders as the
# literal `default=)`, a SyntaxError. R then sees multi-line python output where
# it expected one JSON line and dies with the useless message
#     the condition has length > 1
# Omit the default instead -- it becomes NULL, and length(NULL) == 0, so every
# `if (length(x))` guard behaves the same.
EMPTY_DEF=$(grep -n "default *= *\(character\|numeric\|integer\|logical\)()" bin/*.R || true)
if [ -z "$EMPTY_DEF" ]; then
    good "no empty-vector argparse defaults (invalid python codegen)"
else
    bad "empty-vector argparse default(s) will emit invalid python:"
    printf '%s\n' "$EMPTY_DEF" | sed 's|^|          |'
fi

# The FLAGS list must actually cover the genes that dominated the first real
# oncoplot. A list that misses them is worse than no list: it looks like the
# problem was handled.
FIRST_ONCOPLOT_TOP="MUC12 TRD20A4P MUC4 WASH6P NBPF3 EPPK1 GOLGA6L4 GOLGA6L6 NBPF10 \
RAMEF10 HRNR MUC16 NBPF14 WASHC1 AHNAK2 FLG GOLGA6B GOLGA8K KLF18 LILRB3 MUC5AC \
NBPF1 PRB2 FOXD4L1 TUBB8B"
MISSED=""
for g in $FIRST_ONCOPLOT_TOP; do
    grep -q "\"$g\"" bin/mohq_common.R || MISSED="$MISSED $g"
done
if [ -z "$MISSED" ]; then
    good "FLAGS list covers all 25 top genes from the first real oncoplot"
else
    bad "FLAGS list misses:$MISSED"
fi

# callable_mb must match the TMB numerator. run_oncoprint.R counts coding
# variants, so a whole-genome denominator understates TMB ~90x.
CM=$(grep -oE "^callable_mb: *[0-9]+" params/cardinal_test.yaml | grep -oE "[0-9]+" || echo "")
if [ -n "$CM" ] && [ "$CM" -gt 100 ]; then
    bad "params/cardinal_test.yaml has callable_mb=$CM with a coding numerator (expect ~30-40)"
else
    good "callable_mb (${CM:-unset}) is consistent with a coding TMB numerator"
fi

# --------------------------------------------------------------------------
hdr "2. R: shared library unit tests"

RUNNER=""
if [ -f "$SIF" ] && command -v apptainer >/dev/null 2>&1; then
    RUNNER="apptainer exec $SIF Rscript"
    echo "  using container: $SIF"
elif command -v Rscript >/dev/null 2>&1; then
    RUNNER="Rscript"
    echo "  using system Rscript (container not found at $SIF)"
else
    bad "no R available - build pipeline.sif or module load r"
fi

if [ -n "$RUNNER" ]; then
    if $RUNNER tests/test_mohq_common.R bin/mohq_common.R; then
        good "all shared-library tests passed"
    else
        bad "shared-library tests FAILED (see above)"
    fi
fi

# --------------------------------------------------------------------------
hdr "3. Nextflow: stub run"

if command -v nextflow >/dev/null 2>&1; then
    if nextflow run . -profile test -stub-run \
         --cohort_dir "$TREE/CM/MoHQ-CM-4" --cohort_name MoHQ-CM-4 \
         --run_collection false --run_report false \
         -ansi-log false >/tmp/stub.log 2>&1; then
        good "stub run completed - workflow wiring is valid"
    else
        bad "stub run failed - see /tmp/stub.log"
        tail -25 /tmp/stub.log
    fi
    rm -rf work .nextflow* 2>/dev/null
else
    echo "  nextflow not on PATH - skipping (module load nextflow)"
fi

# --------------------------------------------------------------------------
printf '\n%s\n' "--------------------------------------------------------------"
printf '  %d passed, %d failed\n' "$pass" "$fail"
if [ "$fail" -gt 0 ]; then
    printf '\n  Send me the output above and I will fix it.\n'
    exit 1
fi
printf '  Pre-flight clean. Next: one real cohort, with -resume.\n'

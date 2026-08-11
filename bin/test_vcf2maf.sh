#!/bin/bash
# ===========================================================================
# test_vcf2maf.sh -- prove the VEP setup works on ONE real patient
#
#   bash bin/test_vcf2maf.sh CM/MoHQ-CM-1/MoHQ-CM-1-1
#   bash bin/test_vcf2maf.sh /path/to/local.vcf.gz          # local file
#
# "VEP is installed" is not the goal; "we can turn a PCGR VCF into a MAF with
# real gene symbols" is. This does exactly what the pipeline's VCF2MAF process
# does -- same container, same flags -- on a single patient, then reports
# whether Hugo_Symbol is actually populated.
#
# Run this BEFORE launching a cohort. A whole run that ends in an empty
# oncoplot costs hours; this costs a couple of minutes.
# ===========================================================================
set -uo pipefail

REMOTE="${MOHQ_REMOTE:-juno:d5f8b8e8e3e2442f81573b2f0951013b:MOH-Q}"
VEP_ROOT="${VEP_ROOT:-$HOME/MoHQ/VEP}"
SIF="${VEP_SIF:-$VEP_ROOT/vep_115.sif}"
CACHE="$VEP_ROOT/vep_data"
VEP_RELEASE="${VEP_RELEASE:-115}"
ASSEMBLY="${ASSEMBLY:-GRCh38}"
V2M="${VCF2MAF_DIR:-$HOME/tools/vcf2maf}/vcf2maf.pl"
TARGET="${1:?usage: test_vcf2maf.sh <INST/COHORT/PATIENT | local.vcf.gz>}"

ok()  { printf '  \033[32mOK\033[0m    %s\n' "$1"; }
bad() { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; }
hdr() { printf '\n\033[1m=== %s ===\033[0m\n' "$1"; }

module load StdEnv/2023 apptainer/1.4.5 2>/dev/null || module load apptainer 2>/dev/null

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

# --------------------------------------------------------------------------
hdr "1. prerequisites"
for f in "$SIF" "$V2M"; do
    [ -e "$f" ] && ok "$f" || { bad "missing: $f"; exit 1; }
done
CACHE_DIR="$CACHE/homo_sapiens/${VEP_RELEASE}_${ASSEMBLY}"
[ -d "$CACHE_DIR" ] && ok "cache: $CACHE_DIR" || { bad "missing cache: $CACHE_DIR"; exit 1; }
FASTA=$(ls "$CACHE_DIR"/*.fa.gz 2>/dev/null | head -1)
[ -n "$FASTA" ] && ok "fasta: $(basename "$FASTA")" || bad "no .fa.gz in cache dir"

# --------------------------------------------------------------------------
hdr "2. fetch the VCF"
if [ -f "$TARGET" ]; then
    cp "$TARGET" "$TMP/in.vcf.gz"; PID=$(basename "$TARGET" .vcf.gz)
    ok "local file: $TARGET"
else
    PID=$(basename "$TARGET")
    vcf=$(rclone lsf "$REMOTE/$TARGET/reports/pcgr/" 2>/dev/null \
          | grep -E '_D.*\.vcf\.gz$' | head -1)
    [ -z "$vcf" ] && { bad "no _D VCF under $TARGET/reports/pcgr/"; exit 1; }
    rclone copyto "$REMOTE/$TARGET/reports/pcgr/$vcf" "$TMP/in.vcf.gz" 2>/dev/null \
        || { bad "download failed"; exit 1; }
    ok "$vcf  ($(du -h "$TMP/in.vcf.gz" | cut -f1))"
fi

# --------------------------------------------------------------------------
hdr "3. run vcf2maf (same steps as the pipeline)"
# Identical normalisation to modules/local/vcf2maf.nf: PCGR writes chr1-style
# names, the cache and FASTA use Ensembl-style 1.
gunzip -c "$TMP/in.vcf.gz" | sed -e 's/^chr//' -e 's/^M\t/MT\t/' > "$TMP/in.vcf"
n_in=$(grep -vc '^#' "$TMP/in.vcf" || echo 0)
ok "$n_in variants in"

echo "  running VEP -- a few minutes for a whole-genome VCF ..."
if apptainer exec --bind "$VEP_ROOT:$VEP_ROOT" --bind "$TMP:$TMP" "$SIF" \
        perl "$V2M" \
        --input-vcf "$TMP/in.vcf" \
        --output-maf "$TMP/out.maf" \
        --tumor-id "$PID" \
        --ncbi-build "$ASSEMBLY" \
        --vep-data "$CACHE" \
        --cache-version "$VEP_RELEASE" \
        --ref-fasta "$FASTA" \
        --vep-forks 4 \
        --filter-vcf 0 \
        > "$TMP/log" 2>&1; then
    ok "vcf2maf completed"
else
    bad "vcf2maf failed -- last 25 lines:"; tail -25 "$TMP/log"; exit 1
fi

# --------------------------------------------------------------------------
hdr "4. is Hugo_Symbol actually populated?"
[ -s "$TMP/out.maf" ] || { bad "no MAF produced"; exit 1; }

python3 - "$TMP/out.maf" <<'PY'
import csv, sys, collections
path = sys.argv[1]
rows = []
with open(path) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        rows.append(line)
        break
    rdr = csv.DictReader([rows[0]] + fh.readlines(), delimiter='\t')
    data = list(rdr)

if not data:
    print("  \033[31mFAIL\033[0m  MAF has a header but no variants"); sys.exit(1)

sym = [ (r.get('Hugo_Symbol') or '').strip() for r in data ]
unk = sum(1 for s in sym if s in ('', 'Unknown', '.'))
frac = unk / len(sym)

print(f"  variants out       : {len(sym):,}")
print(f"  Hugo_Symbol Unknown: {unk:,}  ({frac:.1%})")

top = collections.Counter(s for s in sym if s not in ('', 'Unknown', '.')).most_common(10)
print("  most-mutated genes :", ", ".join(f"{g}({n})" for g, n in top) or "none")

if frac > 0.5:
    print("\n  \033[31mFAIL\033[0m  Annotation did not work -- this is the same state as")
    print("        the collection MAFs. Check that --cache-version matches the")
    print("        VEP release in the container.")
    sys.exit(1)
print("\n  \033[32mPASS\033[0m  Real gene symbols. The VEP setup is working end to end.")
PY
rc=$?

# --------------------------------------------------------------------------
hdr "5. sanity check on the result"
cat <<'EOF'
  Look at the gene list above before trusting it:
    * plausible cancer genes (TP53, KRAS, PIK3CA...) -> good
    * only TTN / MUC16 / OBSCN -> those are long genes, not drivers. Not a
      failure: the pipeline filters to non-synonymous variants, which this
      raw test does not.
    * gene names that look like transcript IDs -> cache/version mismatch
EOF

exit $rc

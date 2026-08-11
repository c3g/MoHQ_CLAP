#!/bin/bash
# ===========================================================================
# check_genome_build.sh -- determine the ACTUAL reference build from file
#                          headers, and compare it with what the filename says.
#
#   bash bin/check_genome_build.sh CM/MoHQ-CM-1            # first 5 patients
#   bash bin/check_genome_build.sh CM/MoHQ-CM-1 20         # first 20
#   bash bin/check_genome_build.sh CM/MoHQ-CM-1 20 > build_audit.tsv
#
# WHY THIS EXISTS
# ---------------
# Every file in the collection is named "...grch38..." (27,268 of them; zero
# mention grch37). But a delivery report states that MoHQ-CM-1-1 was ALIGNED
# against GRCh37. If that is true, the filenames are asserting something false,
# and any build inference from filenames -- including the one in
# build_manifest.py -- is unsafe.
#
# So we ask the data instead. Two independent, hard-to-fake signals:
#
#   1. BAM @SQ headers: chromosome LENGTHS.
#        chr1  GRCh37 = 249,250,621     GRCh38 = 248,956,422
#        chr2  GRCh37 = 243,199,373     GRCh38 = 242,193,529
#      A length cannot be renamed. This is the ground truth for alignment.
#
#   2. VCF ##contig=<ID=...,length=...> and ##reference= headers.
#      Tells you what the VARIANT CALLS are on, which is what actually matters
#      for the cohort analyses -- and which may differ from the BAM if a
#      liftover happened somewhere in between.
#
# Only the first few MB of each object is fetched, never the whole BAM.
#
# Output is TSV, so it can be fed straight into build_manifest.py:
#     bash bin/check_genome_build.sh CM/MoHQ-CM-1 999 > builds.tsv
#     python3 bin/build_manifest.py --root ... --build-table builds.tsv
# ===========================================================================
set -uo pipefail

REMOTE="${MOHQ_REMOTE:-juno:d5f8b8e8e3e2442f81573b2f0951013b:MOH-Q}"
COHORT="${1:?usage: check_genome_build.sh <INSTITUTION/COHORT> [n_patients]}"
N="${2:-5}"
BYTES="${HEADER_BYTES:-4000000}"     # 4 MB is ample for a BAM header

command -v samtools >/dev/null 2>&1 || {
    echo "samtools not found. module load samtools/1.22.1" >&2; exit 1; }

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

# --- length -> build -------------------------------------------------------
build_from_len() {
    case "$1" in
        249250621) echo GRCh37 ;;
        248956422) echo GRCh38 ;;
        *)         echo unknown ;;
    esac
}

# --- BAM: read the header via a truncated fetch ----------------------------
build_from_bam() {
    local obj="$1"
    rclone cat --count "$BYTES" "$obj" > "$TMP/h.bam" 2>/dev/null
    # samtools will complain about the truncated stream AFTER emitting the
    # header, which is all we need.
    local hdr; hdr=$(samtools view -H "$TMP/h.bam" 2>/dev/null)
    [ -z "$hdr" ] && { echo "unreadable|"; return; }

    local len; len=$(awk '$1=="@SQ"{
        name=""; ln="";
        for(i=2;i<=NF;i++){
            if($i ~ /^SN:/){name=substr($i,4)}
            if($i ~ /^LN:/){ln=substr($i,4)}
        }
        if(name=="1"||name=="chr1"){print ln; exit}
    }' <<<"$hdr")

    # @PG / @SQ UR: often record the reference path outright.
    local refpath; refpath=$(grep -oE 'UR:[^[:space:]]+' <<<"$hdr" | head -1 | cut -c4-)
    [ -z "$refpath" ] && refpath=$(awk '$1=="@PG"' <<<"$hdr" \
        | grep -oE '/[^[:space:]]*\.(fa|fasta)(\.gz)?' | head -1)

    echo "$(build_from_len "${len:-0}")|${refpath}"
}

# --- VCF: read ##contig / ##reference --------------------------------------
build_from_vcf() {
    local obj="$1"
    local hdr; hdr=$(rclone cat --count 400000 "$obj" 2>/dev/null \
                     | zcat 2>/dev/null | head -3000)
    [ -z "$hdr" ] && { echo "unreadable|"; return; }

    local len; len=$(grep -oE '^##contig=<ID=(chr)?1,length=[0-9]+' <<<"$hdr" \
                     | head -1 | grep -oE '[0-9]+$')
    local ref; ref=$(grep -m1 '^##reference=' <<<"$hdr" | cut -d= -f2-)
    echo "$(build_from_len "${len:-0}")|${ref}"
}

# --------------------------------------------------------------------------
mapfile -t PATIENTS < <(rclone lsf --dirs-only "$REMOTE/$COHORT" 2>/dev/null \
                        | sed 's:/$::' | head -n "$N")
[ ${#PATIENTS[@]} -eq 0 ] && { echo "no patients under $REMOTE/$COHORT" >&2; exit 2; }

printf 'patient_id\tbuild_from_bam\tbuild_from_vcf\tbuild_from_filename\tagreement\tbam_reference\n'

mismatch=0; checked=0
for p in "${PATIENTS[@]}"; do
    # tumour DNA BAM
    bam=$(rclone lsf "$REMOTE/$COHORT/$p/alignment/" 2>/dev/null \
          | grep -E 'DT\.bam$' | head -1)
    [ -z "$bam" ] && bam=$(rclone lsf "$REMOTE/$COHORT/$p/alignment/" 2>/dev/null \
          | grep -E 'DN\.bam$' | head -1)

    # somatic VCF (what the analyses actually consume)
    vcf=$(rclone lsf "$REMOTE/$COHORT/$p/variants/" 2>/dev/null \
          | grep -E 'ensemble\.somatic.*\.vcf\.gz$' | head -1)
    [ -z "$vcf" ] && vcf=$(rclone lsf "$REMOTE/$COHORT/$p/reports/pcgr/" 2>/dev/null \
          | grep -E '_D.*\.vcf\.gz$' | head -1)

    bam_res="no_bam|"; vcf_res="no_vcf|"
    [ -n "$bam" ] && bam_res=$(build_from_bam "$REMOTE/$COHORT/$p/alignment/$bam")
    if [ -n "$vcf" ]; then
        if [[ "$vcf" == *.vcf.gz ]] && rclone lsf "$REMOTE/$COHORT/$p/variants/" 2>/dev/null | grep -q "^${vcf}$"; then
            vcf_res=$(build_from_vcf "$REMOTE/$COHORT/$p/variants/$vcf")
        else
            vcf_res=$(build_from_vcf "$REMOTE/$COHORT/$p/reports/pcgr/$vcf")
        fi
    fi

    b_bam="${bam_res%%|*}"; b_ref="${bam_res#*|}"
    b_vcf="${vcf_res%%|*}"

    # what the FILENAMES claim
    b_name=$(rclone lsf "$REMOTE/$COHORT/$p/reports/pcgr/" 2>/dev/null \
             | grep -oiE 'grch3[78]|hg(19|38)' | head -1 \
             | sed -e 's/^grch37$/GRCh37/I' -e 's/^grch38$/GRCh38/I' \
                   -e 's/^hg19$/GRCh37/I'  -e 's/^hg38$/GRCh38/I')
    [ -z "$b_name" ] && b_name="none"

    real=$(printf '%s\n%s\n' "$b_bam" "$b_vcf" | grep -E '^GRCh3[78]$' | sort -u | tr '\n' ',' | sed 's/,$//')
    if [ -z "$real" ]; then
        agree="undetermined"
    elif [[ "$real" == *,* ]]; then
        agree="BAM_VCF_DISAGREE"; mismatch=$((mismatch+1))
    elif [ "$b_name" = "none" ]; then
        agree="ok_no_filename_claim"
    elif [ "$real" = "$b_name" ]; then
        agree="ok"
    else
        agree="FILENAME_WRONG"; mismatch=$((mismatch+1))
    fi
    checked=$((checked+1))

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$p" "$b_bam" "$b_vcf" "$b_name" "$agree" "$b_ref"
done

{
  echo
  echo "checked $checked patient(s); $mismatch with a problem"
  if [ "$mismatch" -gt 0 ]; then
    cat <<'EOF'

  FILENAME_WRONG   the filename asserts a build the data contradicts.
                   Trust the header, not the name. Report this to the data
                   providers -- anyone else using this collection is exposed
                   to the same error.

  BAM_VCF_DISAGREE alignment and variant calls are on different builds. This
                   can be legitimate (a liftover step) but must be deliberate.
                   For the cohort analyses, what matters is the VCF/segment
                   build, since that is what the coordinates come from.

  Next: run this across all cohorts, save as TSV, and feed it to the manifest:
      python3 bin/build_manifest.py --root ... --build-table builds.tsv
EOF
  fi
} >&2

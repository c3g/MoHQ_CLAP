#!/bin/bash
# ===========================================================================
# check_tool_versions.sh -- which PCGR / VEP / GenPipes versions produced
#                           each patient's files?
#
#   bash bin/check_tool_versions.sh CM/MoHQ-CM-1/MoHQ-CM-1-1     # one patient
#   bash bin/check_tool_versions.sh CM/MoHQ-CM-1 --cohort 10     # sample a cohort
#
# WHY THIS MATTERS BEYOND CURIOSITY
# ---------------------------------
# A multiyear collection is processed in rounds, and tool versions change
# between them. That has two concrete consequences:
#
#  1. FILE FORMATS SHIFT. PCGR renamed columns between major versions -- which
#     is why read_pcgr_cna() in mohq_common.R resolves columns by pattern
#     rather than by fixed name. Knowing the version spread tells you how much
#     variation to expect.
#
#  2. VERSION IS A BATCH VARIABLE. If half the collection was called with PCGR
#     1.x and half with 2.x, differences between those halves are technical.
#     Same logic as `mutation_vcf_source` in the manifest: record it, and never
#     compare across it without adjusting.
#
# Three sources, most informative first:
#   parameters/*.ini   GenPipes config -- records EVERY module version used
#   VCF header         ##PCGR_version, ##VEP, ##source lines
#   PCGR HTML report   version in the page header (largest download; last resort)
# ===========================================================================
set -uo pipefail

REMOTE="${MOHQ_REMOTE:-juno:d5f8b8e8e3e2442f81573b2f0951013b:MOH-Q}"
TARGET="${1:?usage: check_tool_versions.sh <INST/COHORT/PATIENT | INST/COHORT --cohort [n]>}"
MODE="${2:-}"
N="${3:-10}"

hdr() { printf '\n\033[1m=== %s ===\033[0m\n' "$1"; }

# --------------------------------------------------------------------------
inspect_patient() {
    local path="$1" label="$2"

    # ---- 1. GenPipes ini: the richest source ----------------------------- #
    local ini
    ini=$(rclone lsf "$REMOTE/$path/parameters/" 2>/dev/null \
          | grep -iE '\.ini$' | grep -i tumourpair | head -1)
    if [ -n "$ini" ]; then
        local content
        content=$(rclone cat "$REMOTE/$path/parameters/$ini" 2>/dev/null)
        # GenPipes inis list module versions as  module_x = mugqic/tool/version
        local pcgr vep gp
        pcgr=$(grep -iE '^\s*module_pcgr\s*=' <<<"$content" | head -1 | sed 's/.*=\s*//')
        vep=$(grep -iE '^\s*module_vep\s*=|^\s*module_ensembl' <<<"$content" | head -1 | sed 's/.*=\s*//')
        # FIX: the previous version grepped any x.y.z near the word "genpipes",
        # which matched unrelated strings and reported nonsense like "10.12.28"
        # (not a real GenPipes version). Take the module path instead.
        gp=$(grep -iE '^\s*module_(genpipes|mugqic_pipelines)\s*=' <<<"$content" \
             | head -1 | sed 's/.*=\s*//')
        [ -z "$gp" ] && gp=$(grep -oE 'mugqic/(genpipes|mugqic_pipelines)/[0-9][^[:space:]]*' \
                             <<<"$content" | head -1)
        printf '%s\tini\t%s\t%s\t%s\t%s\n' "$label" "${pcgr:-?}" "${vep:-?}" "${gp:-?}" "$ini"
    else
        printf '%s\tini\t-\t-\t-\tno TumourPair ini\n' "$label"
    fi

    # ---- 2. VCF header ---------------------------------------------------- #
    local vcf
    vcf=$(rclone lsf "$REMOTE/$path/reports/pcgr/" 2>/dev/null \
          | grep -E '_D.*\.vcf\.gz$' | head -1)
    if [ -n "$vcf" ]; then
        local h
        h=$(rclone cat --count 300000 "$REMOTE/$path/reports/pcgr/$vcf" 2>/dev/null \
            | zcat 2>/dev/null | head -400)
        # Grep broadly rather than assuming exact key names -- they differ
        # between PCGR releases.
        local vlines
        vlines=$(grep -iE '^##(PCGR|source|VEP|tool|pipeline)' <<<"$h" \
                 | cut -c1-150 | tr '\n' ' | ')
        printf '%s\tvcf\t%s\n' "$label" "${vlines:-no version lines in header}"
    fi
}

# --------------------------------------------------------------------------
if [ "$MODE" = "--cohort" ]; then
    hdr "sampling $N patients from $TARGET"
    printf 'patient\tsource\tpcgr\tvep\tgenpipes\tdetail\n'
    for p in $(rclone lsf --dirs-only "$REMOTE/$TARGET" 2>/dev/null \
               | sed 's:/$::' | head -n "$N"); do
        inspect_patient "$TARGET/$p" "$p"
    done
    cat >&2 <<'EOF'

  Read the `ini` rows: if the pcgr column is not identical across patients, the
  collection was processed in more than one round. Record that as a covariate
  before comparing anything across those groups.
EOF
else
    hdr "$TARGET"
    inspect_patient "$TARGET" "$(basename "$TARGET")"

    hdr "full VCF header (first 60 meta lines)"
    vcf=$(rclone lsf "$REMOTE/$TARGET/reports/pcgr/" 2>/dev/null \
          | grep -E '_D.*\.vcf\.gz$' | head -1)
    if [ -n "$vcf" ]; then
        rclone cat --count 300000 "$REMOTE/$TARGET/reports/pcgr/$vcf" 2>/dev/null \
            | zcat 2>/dev/null | grep '^##' | head -60
    else
        echo "(no _D VCF found)"
    fi

    hdr "GenPipes ini: module versions"
    ini=$(rclone lsf "$REMOTE/$TARGET/parameters/" 2>/dev/null | grep -i '\.ini$' | head -1)
    if [ -n "$ini" ]; then
        echo "from: $ini"
        rclone cat "$REMOTE/$TARGET/parameters/$ini" 2>/dev/null \
            | grep -iE '^\s*module_|^\s*genome|^\s*assembly|^\s*source' \
            | sed 's/^\s*//' | sort -u | head -40
    else
        echo "(no ini found)"
    fi
fi

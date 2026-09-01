#!/bin/bash
# ===========================================================================
# check_csq.sh -- do the PCGR VCFs already carry VEP annotation (a CSQ field)?
#
#   bash bin/check_csq.sh HM/MoHQ-HM-19          # one cohort, first 5 patients
#   bash bin/check_csq.sh HM/MoHQ-HM-19 20       # first 20 patients
#   bash bin/check_csq.sh MU/MoHQ-MU-16 20 --external 20231122_vepvcfs
#
# Why bother: if CSQ is present, vcf2maf --inhibit-vep reuses it and you skip
# the ~30 GB VEP cache and minutes-per-patient VEP runs entirely. If it is
# absent, you need the full VEP install. That is a big enough difference to be
# worth 30 seconds of checking.
#
# Why a script rather than one command: one file is not representative. We have
# already been caught out once -- in the 20-patient sample tree, every patient
# with the older `_D.acmg.` MAF naming had NO PCGR VCF at all, yet
# MoHQ-HM-19-1 has the old MAF naming AND a PCGR VCF. Conventions vary across
# cohorts, so sample several patients before deciding.
#
# --external: CHECK THE BATCH DIRECTORY TOO, OR THE ANSWER IS HALF AN ANSWER
# --------------------------------------------------------------------------
# Some cohorts deliver PCGR VCFs in a dated batch directory beside the patient
# folders instead of inside them. In MoHQ-MU-16, `20231122_vepvcfs/` supplies
# 35 of the 86 delivered PCGR VCFs -- for those patients it is the ONLY one --
# and they are named <patient_id>.pcgr_acmg.grch38.vcf.gz, so the in-patient
# `_D` pattern below does not match them.
#
# Without --external this script would scan only the in-patient files, report
# "all annotated, keep inhibit_vep: true", and say nothing about 40% of the
# delivered cohort. Those batch files are ALSO an older vintage -- 32 of them
# were produced by PCGR 0.9.2 -- so they are the likelier group to differ, and
# they were the group not being looked at.
#
# If the two groups disagree, that is the MIXED verdict: run VEP for everyone.
# An annotated half and a re-annotated half is precisely the systematic
# difference the batch-effect analysis exists to detect.
#
# Only the first ~200 KB of each object is fetched, not the whole VCF.
# ===========================================================================
set -uo pipefail

REMOTE="${MOHQ_REMOTE:-juno:d5f8b8e8e3e2442f81573b2f0951013b:MOH-Q}"
COHORT="${1:?usage: check_csq.sh <INSTITUTION/COHORT> [n_patients] [--external DIR]}"
shift
N=5
EXT=""
# n_patients is optional and positional, so accept it only if it looks like a
# number -- otherwise `check_csq.sh COHORT --external dir` would silently set
# N="--external" and sample zero patients.
if [ $# -gt 0 ] && [[ "$1" =~ ^[0-9]+$ ]]; then N="$1"; shift; fi
while [ $# -gt 0 ]; do
    case "$1" in
        --external) EXT="${2:?--external needs a directory name}"; shift 2 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

echo "cohort : $COHORT"
echo "sample : first $N patient(s)"
[ -n "$EXT" ] && echo "batch  : $EXT"
echo

# Returns 0 = CSQ present, 1 = absent. Prints the INFO keys it did see on
# stderr when absent, so a "no" is explicable rather than bare.
csq_of() {
    local url="$1" hdr
    hdr=$(rclone cat --count 200000 "$url" 2>/dev/null | zcat 2>/dev/null | head -500)
    if grep -q '^##INFO=<ID=CSQ' <<<"$hdr"; then return 0; fi
    grep -oE '^##INFO=<ID=[A-Za-z_]+' <<<"$hdr" \
        | sed 's/.*ID=//' | tr '\n' ',' | cut -c1-60 >&2
    return 1
}

mapfile -t PATIENTS < <(rclone lsf --dirs-only "$REMOTE/$COHORT" 2>/dev/null \
                        | sed 's:/$::' | head -n "$N")
if [ ${#PATIENTS[@]} -eq 0 ]; then
    echo "No patient directories under $REMOTE/$COHORT" >&2
    exit 2
fi

# One listing of the batch directory, reused for every patient below. Without
# it a patient whose only PCGR VCF lives there is reported as "no _D VCF",
# which reads as missing data when the file simply sits elsewhere.
declare -A IN_BATCH=()
if [ -n "$EXT" ]; then
    while read -r f; do
        [ -n "$f" ] && IN_BATCH["${f%%.*}"]="$f"
    done < <(rclone lsf "$REMOTE/$COHORT/$EXT/" 2>/dev/null | grep -E '\.vcf\.gz$')
    echo "  $EXT holds ${#IN_BATCH[@]} PCGR VCF(s)"
    echo
fi

csq=0; nocsq=0; novcf=0
ext_csq=0; ext_nocsq=0

for p in "${PATIENTS[@]}"; do
    vcf=$(rclone lsf "$REMOTE/$COHORT/$p/reports/pcgr/" 2>/dev/null \
          | grep -E '_D.*\.vcf\.gz$' | head -1)

    if [ -z "$vcf" ]; then
        if [ -n "${IN_BATCH[$p]:-}" ]; then
            printf '  %-24s \033[36mbatch dir only\033[0m  %s\n' "$p" "${IN_BATCH[$p]}"
        else
            printf '  %-24s \033[33mno _D VCF\033[0m\n' "$p"
            novcf=$((novcf+1))
        fi
        continue
    fi

    # 2>&1 >/dev/null captures only what csq_of wrote to stderr (the INFO keys
    # it did find) while keeping its exit status. One fetch, not two.
    if keys=$(csq_of "$REMOTE/$COHORT/$p/reports/pcgr/$vcf" 2>&1 >/dev/null); then
        printf '  %-24s \033[32mCSQ present\033[0m   %s\n' "$p" "$vcf"
        csq=$((csq+1))
    else
        printf '  %-24s \033[31mno CSQ\033[0m        INFO keys: %s\n' \
               "$p" "${keys:-none found}"
        nocsq=$((nocsq+1))
    fi
done

# --- the batch directory, sampled independently ----------------------------
if [ -n "$EXT" ] && [ ${#IN_BATCH[@]} -gt 0 ]; then
    echo
    echo "  --- $EXT (older vintage, sampled separately) ---"
    i=0
    for pid in $(printf '%s\n' "${!IN_BATCH[@]}" | sort); do
        [ "$i" -ge "$N" ] && break
        i=$((i+1))
        f="${IN_BATCH[$pid]}"
        if keys=$(csq_of "$REMOTE/$COHORT/$EXT/$f" 2>&1 >/dev/null); then
            printf '  %-24s \033[32mCSQ present\033[0m   %s\n' "$pid" "$f"
            ext_csq=$((ext_csq+1))
        else
            printf '  %-24s \033[31mno CSQ\033[0m        INFO keys: %s\n' \
                   "$pid" "${keys:-none found}"
            ext_nocsq=$((ext_nocsq+1))
        fi
    done
fi

echo
echo "  in-patient   CSQ present : $csq"
echo "  in-patient   no CSQ      : $nocsq"
echo "  in-patient   no _D VCF   : $novcf"
if [ -n "$EXT" ]; then
    echo "  $EXT  CSQ present : $ext_csq"
    echo "  $EXT  no CSQ      : $ext_nocsq"
fi
echo

# The two groups are judged TOGETHER. A cohort where the in-patient files carry
# CSQ and the batch files do not is the mixed case, not a pass.
csq=$((csq + ext_csq)); nocsq=$((nocsq + ext_nocsq))

if [ "$csq" -gt 0 ] && [ "$nocsq" -eq 0 ] && [ "$novcf" -eq 0 ]; then
    cat <<'EOF'
  => All sampled VCFs are already VEP-annotated.
     Set  inhibit_vep: true  and skip the VEP install entirely.
     Still verify on one patient that the resulting MAF has real gene symbols.
EOF
elif [ "$csq" -gt 0 ]; then
    cat <<'EOF'
  => MIXED. Some VCFs carry annotation, some do not.
     Use inhibit_vep: false (run VEP) so every patient is treated identically.
     Mixing annotated and re-annotated calls in one cohort introduces exactly
     the kind of systematic difference the batch-effect analysis is meant to
     detect -- do not do it accidentally.
EOF
    if [ -n "$EXT" ] && [ "$ext_nocsq" -gt 0 ] && [ "$nocsq" -eq "$ext_nocsq" ]; then
        cat <<EOF

     NOTE: every unannotated file found is in $EXT, and the in-patient files
     are all annotated. That is a vintage split, not a scatter -- the batch
     directory is older. It is still the mixed case and still means VEP for
     everyone, but it tells you the cause.
EOF
    fi
else
    cat <<'EOF'
  => No annotation present. You need the full VEP install:
       bash bin/setup_vep.sh --install
EOF
fi

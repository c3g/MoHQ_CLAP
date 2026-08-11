#!/bin/bash
# ===========================================================================
# check_csq.sh -- do the PCGR VCFs already carry VEP annotation (a CSQ field)?
#
#   bash bin/check_csq.sh HM/MoHQ-HM-19          # one cohort, first 5 patients
#   bash bin/check_csq.sh HM/MoHQ-HM-19 20       # first 20 patients
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
# Only the first ~200 KB of each object is fetched, not the whole VCF.
# ===========================================================================
set -uo pipefail

REMOTE="${MOHQ_REMOTE:-juno:d5f8b8e8e3e2442f81573b2f0951013b:MOH-Q}"
COHORT="${1:?usage: check_csq.sh <INSTITUTION/COHORT> [n_patients]}"
N="${2:-5}"

echo "cohort : $COHORT"
echo "sample : first $N patient(s)"
echo

mapfile -t PATIENTS < <(rclone lsf --dirs-only "$REMOTE/$COHORT" 2>/dev/null \
                        | sed 's:/$::' | head -n "$N")
if [ ${#PATIENTS[@]} -eq 0 ]; then
    echo "No patient directories under $REMOTE/$COHORT" >&2
    exit 2
fi

csq=0; nocsq=0; novcf=0

for p in "${PATIENTS[@]}"; do
    vcf=$(rclone lsf "$REMOTE/$COHORT/$p/reports/pcgr/" 2>/dev/null \
          | grep -E '_D.*\.vcf\.gz$' | head -1)

    if [ -z "$vcf" ]; then
        printf '  %-24s \033[33mno _D VCF\033[0m\n' "$p"
        novcf=$((novcf+1)); continue
    fi

    # --count limits the transfer to the header region; zcat will complain
    # about the truncated stream, which is expected and harmless.
    hdr=$(rclone cat --count 200000 "$REMOTE/$COHORT/$p/reports/pcgr/$vcf" 2>/dev/null \
          | zcat 2>/dev/null | head -500)

    if grep -q '^##INFO=<ID=CSQ' <<<"$hdr"; then
        printf '  %-24s \033[32mCSQ present\033[0m   %s\n' "$p" "$vcf"
        csq=$((csq+1))
    else
        # PCGR sometimes writes annotation under other keys; report what we saw
        # rather than a bare "no", so the reason is visible.
        keys=$(grep -oE '^##INFO=<ID=[A-Za-z_]+' <<<"$hdr" \
               | sed 's/.*ID=//' | tr '\n' ',' | cut -c1-60)
        printf '  %-24s \033[31mno CSQ\033[0m        INFO keys: %s\n' "$p" "${keys:-none found}"
        nocsq=$((nocsq+1))
    fi
done

echo
echo "  CSQ present : $csq"
echo "  no CSQ      : $nocsq"
echo "  no _D VCF   : $novcf"
echo
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
else
    cat <<'EOF'
  => No annotation present. You need the full VEP install:
       bash bin/setup_vep.sh --install
EOF
fi

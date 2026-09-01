#!/bin/bash
# ===========================================================================
# why_unknown_pcgr.sh -- why does a patient have pcgr_version = "unknown"?
#
#   bash bin/why_unknown_pcgr.sh [manifest.tsv] [n]
#
# pcgr_version is read from `module_pcgr = ...` in the patient's GenPipes
# TumourPair ini. "unknown" therefore means one of three quite different
# things, and they imply different actions:
#
#   A. no ini was harvested        -> harvest */parameters/*.ini and rebuild
#   B. an ini exists, no module_pcgr line
#                                  -> PCGR WAS NEVER CONFIGURED for this patient
#   C. an ini exists with a differently-spelled key
#                                  -> the manifest parser needs widening
#
# B is the interesting one. If the patients with no PCGR VCF are also the
# patients whose ini never mentions PCGR, then the output was not "produced and
# later deleted" -- that step never ran. Same visible symptom, different
# conversation with the data team, and it changes what regenerating means:
# not restoring a lost file, but producing one that never existed.
# ===========================================================================
set -uo pipefail
MAN="${1:-$(ls audit/*.manifest.tsv results/*/manifest/*.manifest.tsv 2>/dev/null | head -1)}"
N="${2:-6}"
[ -f "$MAN" ] || { echo "usage: $0 <manifest.tsv> [n]"; exit 2; }

echo "manifest: $MAN"
echo

awk -F'\t' -v n="$N" '
NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
$h["pcgr_vcf_status"]=="regenerable" && ($h["pcgr_version"]=="unknown" || $h["pcgr_version"]=="NA") {
    print $h["patient_id"] "\t" $h["tumourpair_ini"]
}' "$MAN" | head -n "$N" | while IFS=$'\t' read -r pid ini; do

    printf '\033[1m%s\033[0m\n' "$pid"
    if [ "$ini" = "NA" ] || [ -z "$ini" ]; then
        echo "    A. no TumourPair ini in the manifest"
        echo "       -> harvest */parameters/*.ini for this cohort and rebuild the manifest"
    elif [ ! -f "$ini" ]; then
        echo "    A. ini recorded but not on disk: $ini"
        echo "       -> harvest */parameters/*.ini and rebuild the manifest"
    else
        echo "    ini: $ini"
        if grep -qiE '^\s*module_pcgr\s*=' "$ini"; then
            echo "    C. module_pcgr IS present -- the parser missed it:"
            grep -iE '^\s*module_pcgr\s*=' "$ini" | sed 's/^/         /'
        else
            echo "    B. NO module_pcgr line in this ini."
            if grep -qiE '^\[report_pcgr\]|pcgr' "$ini"; then
                echo "       but PCGR is mentioned elsewhere:"
                grep -inE '\[report_pcgr\]|pcgr' "$ini" | head -4 | sed 's/^/         /'
            else
                echo "       and PCGR is not mentioned anywhere in the file."
                echo "       -> this patient's run had no PCGR step. The missing VCF was"
                echo "          never produced, not deleted."
            fi
            echo "    module_ lines that ARE present (first 6):"
            grep -iE '^\s*module_' "$ini" | head -6 | sed 's/^/         /'
        fi
    fi
    echo
done

cat <<'EOF'
------------------------------------------------------------------
If most say B: the 27 patients were never PCGR-processed. Worth telling
your supervisor -- it contradicts "produced then deleted", and it means
--pcgr_module_fallback is choosing a version for a run that never happened
rather than restoring one. Still legitimate, but it is a decision to state,
not a restoration to assume.

If most say A: just harvest the inis. Cheap, and it removes the guess:
    python3 bin/harvest_juno.py --remote juno:$MOHQ_PROJECT:MOH-Q \
        --cohort HM/MoHQ-HM-19 --dest /home/$USER/MoHQ/collection \
        --sets core --refresh-listing
    rm -rf .manifest_cache
EOF

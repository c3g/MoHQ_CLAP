#!/bin/bash
# ===========================================================================
# harvest_all.sh -- harvest every MoHQ cohort from Juno. Run on Cardinal.
#
#   ./harvest_all.sh                 # harvest to local scratch on Cardinal
#   ./harvest_all.sh --dry-run       # size it up first, transfer nothing
#
# WHY NOT THE `rclone cat | ssh` LOOP
# -----------------------------------
# The original transfer script opened two SSH connections per file (one to
# test -f the destination, one to receive the stream) and ran strictly serially.
# At ~10 harvested files x 4,500 patients that is ~90,000 sequential round
# trips. harvest_juno.py instead does one listing and one parallel
# `rclone copy --files-from` per cohort.
#
# CROSS-HOST NOTE
# ---------------
# If Juno is only reachable from Cardinal but you compute on Narval, do NOT
# pipe through ssh. Define an SFTP remote and let rclone keep its parallelism:
#
#   # ~/.config/rclone/rclone.conf  (on Cardinal)
#   [narval]
#   type = sftp
#   host = narval.alliancecan.ca
#   user = slavova
#   key_file = ~/.ssh/id_rsa
#
# then set DEST="narval:/home/slavova/scratch/MoHQ/collection".
#
# Better still, if you have an allocation on Cardinal: harvest to Cardinal
# scratch and run the pipeline there. 431 GB is small, and it removes the
# cross-cluster hop entirely.
# ===========================================================================
set -euo pipefail

# --- configure -------------------------------------------------------------
# The project ID below is the public C3G OpenStack project used in the SD4H
# rclone documentation, not a private credential. It still belongs in config
# rather than committed source, so override it via the environment.
REMOTE="${MOHQ_REMOTE:-juno:d5f8b8e8e3e2442f81573b2f0951013b:MOH-Q}"
DEST="${MOHQ_DEST:-$HOME/scratch/MoHQ/collection}"
SETS="${MOHQ_SETS:-core}"          # core | signatures | sv | germline | all
TRANSFERS="${MOHQ_TRANSFERS:-16}"
CACHE_DIR="${MOHQ_CACHE_DIR:-$HOME/.mohq_harvest_cache}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG="${MOHQ_LOG:-$HOME/mohq_harvest_$(date +%Y%m%d_%H%M%S).log}"

echo "remote : $REMOTE"
echo "dest   : $DEST"
echo "sets   : $SETS"
echo "log    : $LOG"

# --- discover institutions and cohorts, rather than hardcoding them --------
# One listing at each level; cheap, and it means a newly added cohort is picked
# up automatically instead of being silently missed.
mapfile -t INSTITUTIONS < <(rclone lsf --dirs-only "$REMOTE" | sed 's:/$::')
echo "institutions: ${INSTITUTIONS[*]}"

COHORTS=()
for inst in "${INSTITUTIONS[@]}"; do
    while IFS= read -r c; do
        [ -n "$c" ] && COHORTS+=("$inst/${c%/}")
    done < <(rclone lsf --dirs-only "$REMOTE/$inst")
done
echo "found ${#COHORTS[@]} cohorts"

# --- harvest ---------------------------------------------------------------
failed=()
for cohort in "${COHORTS[@]}"; do
    echo "=== $cohort ===" | tee -a "$LOG"
    if python3 "$HERE/harvest_juno.py" \
            --remote    "$REMOTE" \
            --cohort    "$cohort" \
            --dest      "$DEST" \
            --sets      $SETS \
            --transfers "$TRANSFERS" \
            --cache-dir "$CACHE_DIR" \
            "$@" 2>&1 | tee -a "$LOG"; then
        echo "OK   $cohort" | tee -a "$LOG"
    else
        echo "FAIL $cohort" | tee -a "$LOG"
        failed+=("$cohort")
    fi
done

echo
echo "=== summary ===" | tee -a "$LOG"
echo "harvested: $(( ${#COHORTS[@]} - ${#failed[@]} )) / ${#COHORTS[@]}" | tee -a "$LOG"
if [ ${#failed[@]} -gt 0 ]; then
    echo "failed cohorts:" | tee -a "$LOG"
    printf '  %s\n' "${failed[@]}" | tee -a "$LOG"
    echo "Re-run harvest_all.sh; completed transfers are skipped by rclone." | tee -a "$LOG"
    exit 1
fi

echo "Total on disk: $(du -sh "$DEST" 2>/dev/null | cut -f1)" | tee -a "$LOG"
echo
echo "Next:"
echo "  python3 $HERE/build_manifest.py --root $DEST --outdir audit/ --cache .cache.json"

#!/bin/bash
# ===========================================================================
# audit_all_builds.sh -- genome-build audit across every institution and cohort
#
#   bash bin/audit_all_builds.sh            > build_audit.tsv    # 3/cohort
#   bash bin/audit_all_builds.sh 10         > build_audit.tsv    # 10/cohort
#
# Spot-check evidence so far is unanimous GRCh38 -- including MoHQ-MU-8-2,
# which the readset report claimed was GRCh37. This runs the same header check
# across the whole collection so the conclusion rests on more than two cohorts,
# and leaves a dated artefact you can attach to a bug report about the
# `Reference` column.
#
# Sampling a few patients per cohort is deliberate: a build change would apply
# to a whole processing round, not to scattered individuals. Raise the count if
# you want stronger coverage.
# ===========================================================================
set -uo pipefail

REMOTE="${MOHQ_REMOTE:-juno:d5f8b8e8e3e2442f81573b2f0951013b:MOH-Q}"
N="${1:-3}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

command -v samtools >/dev/null || { echo "module load samtools/1.22.1" >&2; exit 1; }

echo "auditing $REMOTE ($N patients per cohort)" >&2

printed_header=0
for inst in $(rclone lsf --dirs-only "$REMOTE" 2>/dev/null | sed 's:/$::'); do
    for coh in $(rclone lsf --dirs-only "$REMOTE/$inst" 2>/dev/null | sed 's:/$::'); do
        echo "  -> $inst/$coh" >&2
        out=$(bash "$HERE/check_genome_build.sh" "$inst/$coh" "$N" 2>/dev/null)
        if [ "$printed_header" -eq 0 ]; then
            printf 'institution\tcohort\t%s\n' "$(head -1 <<<"$out")"
            printed_header=1
        fi
        tail -n +2 <<<"$out" | grep -v '^$' | while IFS= read -r line; do
            printf '%s\t%s\t%s\n' "$inst" "$coh" "$line"
        done
    done
done

#!/bin/bash
# ===========================================================================
# check_cardinal.sh -- report what Cardinal gives us, so the pipeline can be
#                      configured from facts instead of guesses.
#
#   bash bin/check_cardinal.sh
#
# Paste the output back and I will fill in the `cardinal` profile for you.
# Reads only; changes nothing.
# ===========================================================================

line() { printf '\n\033[1m--- %s ---\033[0m\n' "$1"; }

line "host"
hostname
cat /etc/os-release 2>/dev/null | grep -E '^(NAME|VERSION)=' || uname -a

line "scheduler"
if command -v sinfo >/dev/null 2>&1; then
    echo "SLURM: yes"
    echo "partitions:"
    sinfo -o "%20P %5a %10l %6D %s" 2>/dev/null | head -12
    echo "your accounts:"
    sacctmgr -nP show assoc user="$USER" format=Account,Partition 2>/dev/null | sort -u | head
    [ -z "$(sacctmgr -nP show assoc user=$USER 2>/dev/null)" ] && echo "  (none reported - you may not need --account here)"
else
    echo "SLURM: no (sinfo not found)"
    command -v qsub >/dev/null && echo "PBS/Torque present" || echo "no obvious scheduler - may be a standalone server"
fi

line "container runtime"
for c in apptainer singularity docker podman; do
    printf '%-12s %s\n' "$c" "$(command -v $c || echo -)"
done
if command -v module >/dev/null 2>&1; then
    echo "modules matching apptainer/singularity:"
    module -t avail 2>&1 | grep -iE 'apptainer|singularity' | head
fi

line "nextflow / java"
printf '%-12s %s\n' nextflow "$(command -v nextflow || echo -)"
printf '%-12s %s\n' java "$(command -v java || echo -)"
java -version 2>&1 | head -1
command -v module >/dev/null 2>&1 && { echo "modules matching nextflow/java:"; module -t avail 2>&1 | grep -iE 'nextflow|^java' | head; }

line "R (only needed if you skip the container)"
printf '%-12s %s\n' Rscript "$(command -v Rscript || echo -)"
command -v module >/dev/null 2>&1 && { echo "modules matching r/:"; module -t avail 2>&1 | grep -iE '^r/|^r-' | head -5; }

line "tools the pipeline shells out to"
for c in rclone perl python3 tabix bcftools; do
    printf '%-12s %s\n' "$c" "$(command -v $c || echo -)"
done
python3 --version 2>/dev/null

line "vcf2maf + VEP (needed for the mutation analyses)"
for p in "$HOME/vcf2maf/vcf2maf.pl" "$HOME/scratch/tools/vcf2maf/vcf2maf.pl" \
         "$HOME/tools/vcf2maf/vcf2maf.pl"; do
    [ -f "$p" ] && echo "found vcf2maf: $p"
done
command -v vep >/dev/null && echo "vep on PATH: $(command -v vep)"
echo "(if neither is here, they need copying from Narval or reinstalling)"

line "storage and quota"
for v in HOME SCRATCH PROJECT; do
    printf '%-10s %s\n' "\$$v" "${!v:-(unset)}"
done
echo
df -h "$HOME" 2>/dev/null | tail -1
[ -n "${SCRATCH:-}" ] && df -h "$SCRATCH" 2>/dev/null | tail -1
command -v diskusage_report >/dev/null && { echo; diskusage_report 2>/dev/null | head -12; }
command -v quota >/dev/null && quota -s 2>/dev/null | head -5

line "juno reachability"
if command -v rclone >/dev/null 2>&1; then
    if rclone lsf --max-depth 1 juno:d5f8b8e8e3e2442f81573b2f0951013b:MOH-Q 2>&1 | head -5; then
        echo "(if you see CM/ CQ/ IQ/ above, Juno works here)"
    fi
else
    echo "rclone not installed - see notes"
fi

line "done"
echo "Paste everything above back and I will write the cardinal profile."

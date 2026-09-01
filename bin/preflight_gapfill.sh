#!/bin/bash
# ===========================================================================
# preflight_gapfill.sh -- ten minutes of checks before an unattended PCGR run
#
#   bash bin/preflight_gapfill.sh params/cardinal_test.yaml
#
# RUN_PCGR has never executed. Everything it depends on -- the PCGR modules,
# the reference bundle, the per-patient version recorded in the manifest -- is
# unverified. Each of these fails in a way that would waste the whole night,
# and each takes seconds to check now.
# ===========================================================================
set -uo pipefail
PARAMS="${1:-params/cardinal_test.yaml}"
ok=0; bad=0
good() { printf '  \033[32mOK\033[0m    %s\n' "$1"; ok=$((ok+1)); }
fail() { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; bad=$((bad+1)); }
warn() { printf '  \033[33mWARN\033[0m  %s\n' "$1"; }
hdr()  { printf '\n\033[1m=== %s ===\033[0m\n' "$1"; }

hdr "1. PCGR versions this cohort needs"
# THE MANIFEST MUST BELONG TO THE COHORT BEING CHECKED.
#
# This used to be `ls audit/*.manifest.tsv results/*/manifest/*.manifest.tsv
# | head -1`, which returns whichever path sorts first -- in practice
# audit/all_cohorts.manifest.tsv, regardless of the cohort in the params file.
# Running it for MoHQ-MU-16 therefore reported MoHQ-HM-19's PCGR versions,
# passed all nine checks, and printed "safe to launch overnight". The real run
# then died on mugqic/pcgr/0.9.2, a version HM-19 does not use and which
# provides no pcgr executable.
#
# A preflight that validates the wrong subject is worse than none: it converts
# an unknown into a false assurance. Resolve the cohort explicitly and fail if
# its manifest is absent, rather than silently checking something else.
COHORT=$(grep -oE '^cohort_name: *"?[^"#]+' "$PARAMS" 2>/dev/null \
         | sed 's/.*: *//;s/"//g;s/ *$//')
OUTDIR=$(grep -oE '^outdir: *"?[^"#]+' "$PARAMS" 2>/dev/null \
         | sed 's/.*: *//;s/"//g;s/ *$//')
if [ -z "$COHORT" ]; then
    fail "cannot read cohort_name from $PARAMS -- cannot verify the right cohort"
    MAN=""
else
    echo "        cohort under test: $COHORT"
    MAN=$(ls "${OUTDIR:-results}/$COHORT/manifest/$COHORT.manifest.tsv" \
             "audit/$COHORT.manifest.tsv" 2>/dev/null | head -1)
fi

if [ -z "$MAN" ]; then
    fail "no manifest for $COHORT. Run the pipeline once (BUILD_MANIFEST is enough)"
    echo "        before the preflight can tell you which PCGR versions it needs."
    echo "        Manifests present, none of them this cohort's:"
    ls audit/*.manifest.tsv "${OUTDIR:-results}"/*/manifest/*.manifest.tsv 2>/dev/null \
        | sed 's/^/          /' | head -5
else
    good "manifest: $MAN"
    VERSIONS=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
                           $h["pcgr_version"]!="NA" && $h["pcgr_version"]!="" {print $h["pcgr_version"]}' \
               "$MAN" | sort -u)
    if [ -z "$VERSIONS" ]; then
        fail "no pcgr_version recorded. RUN_PCGR cannot pick a version per patient."
        echo "        Either harvest */parameters/*.ini, or set --pcgr_module explicitly."
    else
        echo "        versions needed:"; printf '          %s\n' $VERSIONS
        N_UNKNOWN=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
                                 $h["pcgr_vcf_status"]=="regenerable" &&
                                 ($h["pcgr_version"]=="unknown" || $h["pcgr_version"]=="NA")' \
                    "$MAN" | wc -l)
        for v in $VERSIONS; do
            # "unknown" is a manifest sentinel meaning "no usable ini", not a
            # module name. Handled by --pcgr_module_fallback, not by loading it.
            if [ "$v" = "unknown" ]; then
                if [ "$N_UNKNOWN" -gt 0 ]; then
                    warn "$N_UNKNOWN patient(s) have NO recorded PCGR version (no usable ini)."
                    echo "        They cannot be gap-filled unless you choose a version for them:"
                    echo "            --pcgr_module_fallback mugqic/pcgr/1.0.3"
                    echo "        Others keep their own version. The choice is recorded in"
                    echo "        mutation_provenance.tsv as pcgr_regenerated_fallback."
                fi
                continue
            fi
            if module load "$v" 2>/dev/null; then
                good "module loads: $v"; module unload "$v" 2>/dev/null
            else
                fail "module NOT loadable: $v   <-- every patient on this version will fail"
            fi
        done
    fi
    # How many patients will actually go through PCGR?
    N=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
                    $h["pcgr_vcf_status"]=="regenerable"' "$MAN" | wc -l)
    echo "        patients to gap-fill: $N"
    [ "$N" -gt 0 ] || warn "nothing to gap-fill -- the run will just redo the analyses"
fi

hdr "1b. Is the pcgr EXECUTABLE actually reachable?"
# Loading the module is not the same as the command existing, and RUN_PCGR runs
# on the HOST (container = null) precisely so the module is usable. The first
# gap-fill attempt failed with exit 127 on every patient because the process
# still inherited the global container.
# EVERY version, not just the first.
#
# This tested only `head -1` of the version list. MoHQ-MU-16 needs
# mugqic/pcgr/0.9.2, which LOADS cleanly and provides no pcgr executable --
# so the module check in section 1 passes and the patients still die with
# exit 127. Whichever version happened to sort first was the only one proved.
for v in ${VERSIONS:-}; do
    [ "$v" = "unknown" ] && continue
    if module load "$v" 2>/dev/null; then
        if command -v pcgr >/dev/null 2>&1; then
            good "pcgr on PATH after $v: $(command -v pcgr)"
        else
            # COUNT FIRST, THEN JUDGE.
            #
            # VERSIONS lists every version recorded in the manifest, including
            # those belonging to patients whose PCGR VCF was DELIVERED. PCGR is
            # never invoked for those patients, so a broken module they merely
            # happen to be labelled with cannot fail anything.
            #
            # Calling fail() before counting made MoHQ-MU-16 report "1 BLOCKING
            # problem" for mugqic/pcgr/0.9.2 while printing "patients affected:
            # 0" directly underneath. A preflight that cries wolf gets ignored,
            # and then a real blocker gets ignored with it -- the same failure
            # mode as the earlier bug that validated the wrong cohort, just
            # pointing the other way.
            n_v=$(awk -F'\t' -v want="$v" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
                    $h["pcgr_vcf_status"]=="regenerable" && $h["pcgr_version"]==want' \
                  "$MAN" 2>/dev/null | wc -l)
            if [ "${n_v:-0}" -gt 0 ]; then
                fail "$v loads but provides NO 'pcgr' executable -- exit 127 for ${n_v} patient(s)"
                echo "        remedy: map them to a working version with --pcgr_module_map"
                echo "                e.g. --pcgr_module_map '$v=mugqic/pcgr/1.0.3'"
                echo "        (recorded in mutation_provenance.tsv, so the substitution is visible)"
            else
                warn "$v provides no 'pcgr' executable, but 0 patients need it regenerated"
                echo "        Every patient on $v already has a DELIVERED PCGR VCF, so PCGR"
                echo "        is never invoked for them. Not blocking."
            fi
        fi
        module unload "$v" 2>/dev/null
    fi
done
grep -q "container = null" modules/local/pcgr.nf \
    && good "RUN_PCGR sets container = null (runs on the host, so modules work)" \
    || fail "RUN_PCGR has NO container override -- module loads on the host, script runs in pipeline.sif -> exit 127"

hdr "2. PCGR reference bundle"
BUNDLE=$(grep -oE "pcgr_bundle *= *'[^']+'" nextflow.config | head -1 | cut -d"'" -f2)
BUNDLE=$(grep -oE "^pcgr_bundle: *\"?[^\"]+\"?" "$PARAMS" 2>/dev/null | sed 's/.*: *//;s/"//g' || echo "$BUNDLE")
if [ -d "$BUNDLE" ]; then
    good "bundle exists: $BUNDLE"
    ASM=$(grep -oE "^genome_build: *\"?[A-Za-z0-9]+" "$PARAMS" 2>/dev/null | sed 's/.*: *//;s/"//' )
    ASM_LC=$(echo "${ASM:-GRCh38}" | tr 'A-Z' 'a-z')
    [ -d "$BUNDLE/data/$ASM_LC" ] && good "assembly data present: $BUNDLE/data/$ASM_LC" \
                                  || fail "no $ASM_LC subdirectory under the bundle"
else
    fail "bundle NOT found: $BUNDLE"
fi

hdr "3. Disk headroom"
DEST=$(grep -oE "^outdir: *\"?[^\"]+" "$PARAMS" 2>/dev/null | sed 's/.*: *//;s/"//')
AVAIL=$(df -BG "${DEST:-$HOME}" 2>/dev/null | awk 'NR==2{gsub("G","",$4); print $4}')
if [ -n "$AVAIL" ]; then
    [ "$AVAIL" -gt 50 ] && good "${AVAIL} GB free" \
                        || fail "only ${AVAIL} GB free -- PCGR writes several GB per patient"
else
    warn "could not determine free space"
fi

hdr "4. Workflow still wires up"
if command -v nextflow >/dev/null 2>&1; then
    # --max_cpus/--max_memory belt-and-braces: conf/base.config already clamps
    # under -stub-run, but if that file is stale the stub fails with
    #     Process requirement exceeds available CPUs -- req: 4; avail: 1
    # which looks like a gap-fill problem and is not one.
    if nextflow run . -profile cardinal_local -stub-run -params-file "$PARAMS" \
         --run_pcgr_gapfill true --pcgr_module_fallback mugqic/pcgr/1.0.3 \
         --max_cpus 1 --max_memory 4.GB >/tmp/pf_stub.log 2>&1; then
        good "stub run passes with --run_pcgr_gapfill true"
        grep -oE "RUN_PCGR +\| *[0-9]+ of [0-9]+" /tmp/pf_stub.log | tail -1 | sed 's/^/        /'
    else
        fail "stub run FAILED -- fix before launching."
        # Show the actual cause. `tail` alone usually catches only Nextflow's
        # trailing noise ("Join mismatch...") and hides the real message.
        grep -m2 -A6 -E "^ERROR|Parameter validation failed|needs a regenerated" \
            /tmp/pf_stub.log | sed 's/^/        /' || tail -15 /tmp/pf_stub.log | sed 's/^/        /'
    fi
else
    warn "nextflow not on PATH (source bin/load_cardinal_modules.sh)"
fi

hdr "5. One real patient, to prove PCGR runs at all"
# Unquoted heredoc, so the command below names the cohort ACTUALLY under test.
# It used to hardcode params/cardinal_test.yaml and MoHQ-HM-19 regardless of
# $PARAMS -- advice about a different cohort, which is how the wrong-manifest
# bug started.
cat <<EOF
  This is the check that matters most, and the only one this script cannot do
  for you quickly -- it needs a real PCGR execution (~1 h).

  If you have an hour before sleeping, run ONE patient first:

      nextflow run . -profile cardinal,vcf2maf_host \\
          -params-file $PARAMS \\
          --run_pcgr_gapfill true --cohort_name $COHORT \\
          -c conf/overnight.config -resume

  and watch until the first RUN_PCGR task completes. If it does, the remaining
  patients differ only in input. If you would rather not wait, launch anyway --
  conf/overnight.config keeps a failure from killing the whole run.
EOF

printf '\n%s\n' "----------------------------------------------------------------"
if [ "$bad" -eq 0 ]; then
    echo "  $ok checks passed, nothing blocking. Safe to launch overnight."
else
    echo "  $bad BLOCKING problem(s). Fix these first -- each one would fail"
    echo "  every gap-fill task and waste the night."
fi

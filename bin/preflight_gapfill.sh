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
MAN=$(ls audit/*.manifest.tsv results/*/manifest/*.manifest.tsv 2>/dev/null | head -1)
if [ -z "$MAN" ]; then
    fail "no manifest found -- run the pipeline once without gapfill first"
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
FIRSTV=$(printf '%s\n' $VERSIONS | grep -v '^unknown$' | head -1)
if [ -n "${FIRSTV:-}" ] && module load "$FIRSTV" 2>/dev/null; then
    if command -v pcgr >/dev/null 2>&1; then
        good "pcgr on PATH after loading $FIRSTV: $(command -v pcgr)"
    else
        fail "module $FIRSTV loads but 'pcgr' is NOT on PATH -- exit 127 on every patient"
        echo "        executables the module does provide:"
        module show "$FIRSTV" 2>&1 | grep -iE "PATH|prepend" | head -3 | sed 's/^/          /'
    fi
    module unload "$FIRSTV" 2>/dev/null
fi
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
cat <<'EOF'
  This is the check that matters most, and the only one this script cannot do
  for you quickly -- it needs a real PCGR execution (~1 h).

  If you have an hour before sleeping, run ONE patient first:

      nextflow run . -profile cardinal,vcf2maf_host \
          -params-file params/cardinal_test.yaml \
          --run_pcgr_gapfill true --cohort_name MoHQ-HM-19 \
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

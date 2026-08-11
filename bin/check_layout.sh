#!/bin/bash
# ===========================================================================
# check_layout.sh -- is every file where the pipeline expects it?
#
#   cd ~/pipeline && bash bin/check_layout.sh
#
# Nextflow is strict about two things and silent about a third:
#   * main.nf and nextflow.config must sit at the project root
#   * include statements resolve literal paths -- modules/local/x.nf, so a file
#     at modules/x.nf is simply not found
#   * scripts in bin/ are added to PATH, but ONLY run if they are executable.
#     A non-executable bin/ script fails at task run time with "command not
#     found", long after you would rather have known.
# ===========================================================================
set -uo pipefail

ok()   { printf '  \033[32m ok \033[0m %s\n' "$1"; }
miss() { printf '  \033[31mMISS\033[0m %s\n' "$1"; MISSING=$((MISSING+1)); }
warn() { printf '  \033[33mfix \033[0m %s\n' "$1"; }
MISSING=0

ROOT_FILES="main.nf nextflow.config pipeline.def"
CONF="conf/base.config"
MODULES="modules/local/manifest.nf modules/local/analysis.nf modules/local/collection.nf
         modules/local/vcf2maf.nf modules/local/pcgr.nf"
ASSETS="assets/report_template.Rmd assets/NO_COHORT_MAP"
PARAMS="params/cardinal_test.yaml"
BIN_PY="bin/build_manifest.py bin/harvest_juno.py"
BIN_R="bin/mohq_common.R bin/build_cohort_seg.R bin/run_oncoprint.R
       bin/run_cnv_frequency.R bin/run_fga_burden.R bin/run_pca.R
       bin/run_comparative.R bin/run_fusions.R bin/plot_completeness.R
       bin/collection_rollup.R bin/batch_effect_analysis.R"
BIN_SH="bin/load_cardinal_modules.sh bin/setup_vep.sh bin/harvest_all.sh
        bin/check_genome_build.sh bin/check_tool_versions.sh bin/check_csq.sh
        bin/audit_all_builds.sh bin/test_vcf2maf.sh bin/check_cardinal.sh"

printf '\n\033[1m=== project root (%s) ===\033[0m\n' "$PWD"
for f in $ROOT_FILES; do [ -f "$f" ] && ok "$f" || miss "$f"; done

printf '\n\033[1m=== conf/ ===\033[0m\n'
for f in $CONF; do [ -f "$f" ] && ok "$f" || miss "$f"; done

printf '\n\033[1m=== modules/local/  (note: TWO levels) ===\033[0m\n'
for f in $MODULES; do [ -f "$f" ] && ok "$f" || miss "$f"; done
if [ -d modules ] && ! [ -d modules/local ]; then
    warn "modules/ exists but modules/local/ does not -- move the .nf files down one level:"
    warn "    mkdir -p modules/local && mv modules/*.nf modules/local/"
fi

printf '\n\033[1m=== assets/ and params/ ===\033[0m\n'
for f in $ASSETS $PARAMS; do [ -f "$f" ] && ok "$f" || miss "$f"; done

printf '\n\033[1m=== bin/  (Nextflow puts this on PATH) ===\033[0m\n'
for f in $BIN_PY $BIN_R $BIN_SH; do [ -f "$f" ] && ok "$f" || miss "$f"; done

printf '\n\033[1m=== executable bits ===\033[0m\n'
NOEXEC=0
for f in $BIN_PY $BIN_R $BIN_SH; do
    [ -f "$f" ] && [ ! -x "$f" ] && { warn "not executable: $f"; NOEXEC=$((NOEXEC+1)); }
done
if [ "$NOEXEC" -gt 0 ]; then
    warn "fix all of them at once:   chmod +x bin/*"
else
    ok "all bin/ scripts are executable"
fi

printf '\n\033[1m=== optional (nice to have) ===\033[0m\n'
for f in tests/run_tests.sh tests/test_mohq_common.R \
         assets/samplesheet_example.csv assets/cohort_map_example.tsv \
         params/example.yaml START_HERE.md; do
    [ -f "$f" ] && ok "$f" || printf '  \033[90m  - \033[0m %s (optional)\n' "$f"
done

# ---------------------------------------------------------------------------
# FRESHNESS
#
# Existence is not enough. Every file below is copied between machines by hand,
# and a stale copy passes every existence check while failing at run time --
# usually several minutes in, with an error that points somewhere else.
#
# Each entry is a string that only the CURRENT version of that file contains.
# If the marker is missing, the file on disk is older than the one it should be.
# ---------------------------------------------------------------------------
printf '\n%s\n' "=== freshness (is each file the current version?) ==="
STALE=0
check_marker() {
    local file="$1" marker="$2" what="$3"
    if [ ! -f "$file" ]; then
        printf '  \033[31mMISS\033[0m  %s\n' "$file"; STALE=$((STALE+1)); return
    fi
    if grep -q -- "$marker" "$file" 2>/dev/null; then
        ok "$file"
    else
        printf '  \033[31mSTALE\033[0m %s  (missing: %s)\n' "$file" "$what"
        STALE=$((STALE+1))
    fi
}

check_marker modules/local/analysis.nf   "exclude_flags"        "FLAGS filter passed to ONCOPRINT"
check_marker modules/local/vcf2maf.nf    "min_annotated_genes"  "gene-count annotation guard"
check_marker bin/mohq_common.R           "FLAGS_GENES"          "FLAGS gene list"
check_marker bin/run_oncoprint.R         "drop_flags"           "FLAGS filtering"
check_marker bin/run_fga_burden.R        "implausible"          "FGA >= 99% flag"
check_marker bin/run_cnv_frequency.R     "GRCH38_CHROM_LEN"     "chromosome-length panel scaling"
check_marker bin/build_manifest.py       "path-root"            "absolute-path manifest fix"
check_marker main.nf                     "PLOT_COMPLETENESS.out.plots" "completeness figure staged to report"
check_marker nextflow.config             "vcf2maf_host"         "host-modules vcf2maf profile"
check_marker assets/report_template.Rmd  "Eligible patients"    "per-analysis denominators"
check_marker bin/plot_completeness.R     "absent for all"       "three-state heatmap (harvest scope)"
check_marker bin/build_manifest.py       "rna_variants"         "rna_variants/ directory scanned"

if [ "$STALE" -gt 0 ]; then
    printf '\n  \033[31m%d file(s) are STALE or missing.\033[0m Re-copy them before running.\n' "$STALE"
    printf '  A stale file passes every existence check and then fails mid-run.\n'
fi

printf '\n%s\n' "------------------------------------------------------------"
if [ "$MISSING" -eq 0 ]; then
    echo "  Layout is complete."
    [ "$NOEXEC" -gt 0 ] && echo "  Run: chmod +x bin/*" \
                        || echo "  Next: nextflow run . -profile cardinal_local -stub-run ..."
else
    echo "  $MISSING required file(s) missing -- see MISS lines above."
fi

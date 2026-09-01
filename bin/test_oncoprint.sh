#!/bin/bash
# ===========================================================================
# test_oncoprint.sh -- run run_oncoprint.R standalone, outside Nextflow
#
#   bash bin/test_oncoprint.sh MoHQ-CM-6
#   bash bin/test_oncoprint.sh MoHQ-CM-6 assets/panel_hartwig.txt
#
# WHY OUTSIDE NEXTFLOW
# --------------------
# Iterating on an R script through the pipeline means: edit, launch, wait for
# the scheduler, then hunt for the error in a work directory. Standalone, the
# same script runs against the same inputs in seconds and errors land in the
# terminal.
#
# Outputs go to a scratch directory so nothing here can overwrite the real
# results in ~/MoHQ/results. Compare, then re-run the pipeline once it is right.
# ===========================================================================
set -uo pipefail

COHORT="${1:?usage: test_oncoprint.sh <COHORT_ID> [gene_panel.txt]}"
PANEL="${2:-}"

# Coding lengths for the length-normalised ranking. Override with GTF=... if
# you want a different Ensembl release; 104 is closest to the VEP 105 that
# annotated these calls, so the gene symbols line up.
GTF="${GTF:-/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl104.gtf}"
[[ -f "$GTF" ]] || { echo "note: GTF not found, skipping length normalisation: $GTF" >&2; GTF=""; }

RES="$HOME/MoHQ/results/$COHORT"
MAFDIR="$HOME/MoHQ/derived_mafs/$COHORT"
OUT="$HOME/MoHQ/oncoprint_test/$COHORT"
SIF="$HOME/pipeline/pipeline.sif"

MANIFEST=$(ls "$RES"/manifest/*.manifest.annotated.tsv \
              "$RES"/manifest/*.manifest.tsv 2>/dev/null | head -1)
[[ -n "$MANIFEST" ]] || { echo "no manifest under $RES/manifest/" >&2; exit 1; }

mapfile -t MAFS < <(ls "$MAFDIR"/*.maf 2>/dev/null)
[[ ${#MAFS[@]} -gt 0 ]] || { echo "no MAFs under $MAFDIR" >&2; exit 1; }

# cna_segments paths live in the manifest; the module passes them as one arg.
# Resolve them the same way so the copy-number overlay is exercised too --
# testing without it would leave the cnTable join untested, and that join is
# where a VEP-release mismatch would show up.
CNA=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
                  $h["cna_segments"]!="NA" && $h["cna_segments"]!="" {print $h["cna_segments"]}' \
      "$MANIFEST" | tr '\n' ' ')

mkdir -p "$OUT"
echo "cohort   : $COHORT"
echo "manifest : $MANIFEST"
echo "mafs     : ${#MAFS[@]}"
echo "cna      : $(wc -w <<<"$CNA") file(s)"
echo "panel    : ${PANEL:-none (all genes)}"
echo "out      : $OUT"
echo

cd "$OUT" || exit 1

# --bind /cvmfs, or the container cannot see the GTF, the PCGR bundle or the
# VEP cache -- all of which live there. Without it R's file.exists() is FALSE
# for a file that is plainly on disk, and the step reports "no coding lengths
# could be read" about a path you just verified by hand.
BIND=""
[[ -d /cvmfs ]] && BIND="--bind /cvmfs:/cvmfs:ro"

apptainer exec --cleanenv $BIND "$SIF" Rscript "$HOME/pipeline/bin/run_oncoprint.R" \
    --lib        "$HOME/pipeline/bin/mohq_common.R" \
    --mafs       "${MAFS[@]}" \
    ${CNA:+--cna $CNA} \
    --manifest   "$MANIFEST" \
    --cohort     "$COHORT" \
    --callable_mb 30 \
    --exclude_flags TRUE \
    --oncodrive     TRUE \
    --oncodrive_min_mut 5 \
    ${GTF:+--gtf "$GTF"} \
    ${PANEL:+--gene_panel "$PANEL"} \
    --out_prefix "${COHORT}_test"

rc=$?
echo
echo "exit: $rc"
[ $rc -ne 0 ] && exit $rc

# Which amino-acid column do these MAFs actually carry? oncodrive needs one,
# and knowing the answer takes a second -- far cheaper than reading a stack
# trace from inside maftools.
echo "=== amino-acid column present in the MAFs ==="
# grep -v '^#' first: a MAF opens with '#version 2.4', so head -1 returns the
# comment and not the header. That made this report NONE FOUND against a file
# whose column 37 is HGVSp_Short -- a check that contradicted the run it was
# meant to explain.
grep -v '^#' "${MAFS[0]}" | head -1 | tr '\t' '\n' \
  | grep -nE '^(HGVSp_Short|HGVSp|AAChange|Protein_Change|amino_acid_change)$' \
  || echo "  NONE FOUND -- oncodrive cannot score positional clustering"

echo
echo "=== outputs ==="
ls -la "$OUT" | grep -E '\.(png|tsv)$'

echo
echo "=== raw count vs mutations-per-kb, top 20 by raw count ==="
RANK="$OUT/${COHORT}_test_gene_rank_unfiltered.tsv"
if [ -f "$RANK" ] && head -1 "$RANK" | grep -q mut_per_kb; then
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;
                      printf "%-14s %7s %7s %8s %9s %8s\n",
                             "gene","pts","vars","cds_kb","mut/kb","flags"; next}
                NR<=21 {printf "%-14s %7s %7s %8.1f %9.2f %8s\n",
                        $h["Hugo_Symbol"], $h["n_patients"], $h["n_variants"],
                        $h["cds_kb"]+0, $h["mut_per_kb"]+0, $h["is_flags"]}' "$RANK"
    echo
    echo "=== top 20 by mutations-per-kb (>= 0.5 kb CDS) ==="
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
                $h["rank_per_kb"] != "NA" && $h["rank_per_kb"] != "" {
                  printf "%s\t%s\t%s\t%s\n",
                         $h["rank_per_kb"], $h["Hugo_Symbol"],
                         $h["mut_per_kb"], $h["n_patients"]}' "$RANK" \
      | sort -n | head -20 \
      | awk -F'\t' 'BEGIN{printf "%-6s %-14s %9s %6s\n","rank","gene","mut/kb","pts"}
                    {printf "%-6s %-14s %9.2f %6s\n",$1,$2,$3,$4}'
else
    echo "  (no length-normalised ranking -- was --gtf supplied?)"
fi

echo
echo "=== frequency vs clustering ==="
CMP="$OUT/${COHORT}_test_rank_comparison.tsv"
if [ -f "$CMP" ]; then
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;
                      printf "%-14s %6s %6s %8s %8s\n",
                             "gene","freq","clust","fdr","is_flags"; next}
                NR<=16 {printf "%-14s %6s %6s %8s %8s\n",
                        $h["Hugo_Symbol"], $h["rank_frequency"],
                        $h["rank_clustering"], $h["fdr"], $h["is_flags"]}' "$CMP"
else
    echo "  (none -- oncodrive needs a protein-change column these MAFs lack)"
fi

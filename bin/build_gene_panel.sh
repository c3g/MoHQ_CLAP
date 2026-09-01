#!/bin/bash
# ===========================================================================
# build_gene_panel.sh -- an allow-list of driver genes, from data already here
#
#   bash bin/build_gene_panel.sh ~/MoHQ/collection > assets/panel_hartwig.txt
#   bash bin/build_gene_panel.sh ~/MoHQ/collection/CM/MoHQ-CM-6 > panel_cm6.txt
#
# WHY THIS EXISTS
# ---------------
# Ranking genes by mutation frequency in whole-genome data returns long genes,
# because a long gene is a bigger target for passenger mutations. Filtering the
# offenders one family at a time does not converge: each pass removes the
# current tier and the next appears. MoHQ-CM-6 was the fourth demonstration --
# BRAF ranked third, behind RP1 and MGAM2.
#
# An allow-list converges by construction: it names what you want rather than
# everything you do not.
#
# COSMIC's Cancer Gene Census needs registration and cannot be redistributed.
# But the collection ALREADY SHIPS a driver gene panel, implicitly: PURPLE and
# LINX are run against Hartwig's DriverGenePanel, and every
# *.driver.catalog.somatic.tsv names the genes from it that were called in that
# patient. Union those across a cohort -- or the whole collection -- and you
# have a defensible, locally derived panel with no download and no licence.
#
# WHAT THIS PANEL IS AND IS NOT
# -----------------------------
# It is the set of driver genes Hartwig looked for AND found somewhere in the
# input. It is therefore biased toward what these tumour types actually carry,
# which is usually what you want for a cohort of this collection -- and it is
# NOT the complete Hartwig panel, because a gene never altered in any patient
# here never appears in any catalogue. Build it over the WHOLE collection
# rather than one cohort if you want the broader list.
# ===========================================================================
set -uo pipefail

ROOT="${1:?usage: build_gene_panel.sh <collection or cohort dir> [> panel.txt]}"
[[ -d "$ROOT" ]] || { echo "not a directory: $ROOT" >&2; exit 1; }

# Hartwig writes the symbol in a column called `gene`. Older releases used
# `Gene`. Resolve by NAME rather than position -- column order is not stable
# across releases, and picking column 1 silently produced chromosome names.
extract() {
    awk -F'\t' '
      FNR == 1 {
        col = 0
        for (i = 1; i <= NF; i++)
          if (tolower($i) == "gene" || tolower($i) == "hugo_symbol") col = i
        next
      }
      col && $col != "" { print $col }
    ' "$@"
}

mapfile -t FILES < <(find "$ROOT" -name '*driver.catalog.somatic.tsv' -o \
                                  -name '*.linx.driver.catalog.tsv' 2>/dev/null)

if [ ${#FILES[@]} -eq 0 ]; then
    echo "No driver catalogue files under $ROOT." >&2
    echo "Harvest them with:  --sets core   (they are in the core set)" >&2
    exit 2
fi

echo "# Driver gene panel derived from ${#FILES[@]} Hartwig driver catalogue(s)"
echo "# source: $ROOT"
echo "# built:  $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "#"
echo "# One HGNC symbol per line. Consumed by run_oncoprint.R --gene_panel,"
echo "# which keeps only variants whose Hugo_Symbol appears here."
echo "#"
echo "# This is the set of driver genes Hartwig SEARCHED FOR and FOUND in this"
echo "# input -- not the full Hartwig panel, and not a statement that every gene"
echo "# listed is a driver in every tumour type."

extract "${FILES[@]}" | sort -u | grep -vE '^(gene|Gene|NA|)$'

# Counts to stderr so they do not land in the panel file itself.
n=$(extract "${FILES[@]}" | sort -u | grep -cvE '^(gene|Gene|NA|)$')
echo "[panel] ${#FILES[@]} catalogue(s) -> $n distinct gene symbol(s)" >&2
if [ "$n" -lt 50 ]; then
    echo "[panel] WARNING: only $n genes. That is small for a driver panel." >&2
    echo "        Build over the whole collection rather than one cohort," >&2
    echo "        or check the catalogues actually have a 'gene' column." >&2
fi

#!/bin/bash
# ===========================================================================
# load_cardinal_modules.sh -- everything the pipeline needs, from modules.
#
#   source bin/load_cardinal_modules.sh
#
# Must be SOURCED, not executed, or the modules load into a subshell and vanish.
#
# Correction to what I told you earlier: nextflow IS available on Cardinal
# (13 versions). My check script grepped `module -t avail` and missed it. You
# do not need to install anything by hand.
# ===========================================================================

if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
    echo "ERROR: source this file, don't run it:" >&2
    echo "    source bin/load_cardinal_modules.sh" >&2
    exit 1
fi

# Nextflow. 26.04.4 is the default but is very new; 25.04.6 is a safer choice
# for a pipeline you are debugging. Both satisfy the >=23.04 requirement.
module load nextflow/25.04.6 2>/dev/null || module load nextflow

# Containers.
module load apptainer/1.3.5 2>/dev/null || module load apptainer

# Perl for vcf2maf.pl. The MUGQIC build is the one GenPipes uses.
module load mugqic/perl/5.40.1 2>/dev/null || module load perl

# htslib gives tabix/bgzip, which some steps expect on PATH.
module load mugqic/htslib/1.23 2>/dev/null || module load htslib 2>/dev/null || true

echo "loaded:"
for m in nextflow apptainer perl htslib; do
    printf '  %-10s %s\n' "$m" "$(command -v $m 2>/dev/null || echo '(not on PATH)')"
done
nextflow -version 2>&1 | grep -i version | head -2

cat <<'EOF'

If you are running WITHOUT the container (-profile cardinal,no_container),
also load an R stack. Two candidates on Cardinal:

    module load r-bundle-bioconductor/3.21          # Alliance build
    module load mugqic/R_Bioconductor/4.5.1_3.21    # MUGQIC build (used by GenPipes)

Check which one already has what the scripts need:

    Rscript -e 'for (p in c("data.table","dplyr","ggplot2","tidyr","argparse",
                            "plotly","htmlwidgets","maftools","rmarkdown","DT"))
                  cat(sprintf("%-14s %s\n", p, requireNamespace(p, quietly=TRUE)))'

Anything FALSE installs to a personal library:
    mkdir -p ~/R && Rscript -e 'install.packages("argparse", lib="~/R")'
    export R_LIBS_USER=~/R
EOF

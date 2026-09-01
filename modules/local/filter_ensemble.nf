// ===========================================================================
// modules/local/filter_ensemble.nf -- the step between the ensemble VCF and PCGR
//
// WHY THIS EXISTS
// ---------------
// RUN_PCGR used to be fed the raw ensemble somatic VCF, on the stated premise
// that this was what GenPipes gave PCGR. It was not. GenPipes' `report_pcgr`
// consumes the FILTERED ensemble VCF, and the difference is not cosmetic:
//
//   ensemble.somatic.vt.annot.vcf.gz          423,764 variants   (raw union)
//   ...annot.2caller.flt.vcf.gz                23,604 variants   (PCGR input)
//
// Ninety-four percent of the calls are removed before PCGR ever sees them.
// Regenerating from the raw file produced VCFs with ~365,000 variants where
// the delivered ones have ~24,000 -- a fifteen-fold difference that correlates
// exactly with which patients happened to lose files. That is the batch
// artefact this pipeline exists to prevent, manufactured by the process meant
// to prevent it.
//
// The chain was recovered from two pieces of primary evidence rather than
// inferred: the [filter_ensemble] section of the delivered runs' config trace,
// and the ##bcftools_viewCommand line preserved in the delivered PCGR VCF's
// own header.
//
//   1. format2pcgr.py  -- adds TAL/TDP/TVAF/NDP/NVAF as INFO tags, and drops
//                         variants supported by fewer than `call_filter`
//                         callers (MOH used 2)              -> .2caller.vcf.gz
//   2. bcftools view -i'TDP>=10 && TVAF>=0.05 &&
//                       NDP>=10 && NVAF<=0.05'              -> .2caller.flt.vcf.gz
//   3. PCGR
//
// WHY NOT REIMPLEMENT THE TAG MATH
// --------------------------------
// It looks trivial -- TVAF is just alt/(ref+alt). It is not: format2pcgr.py
// branches on the FIRST caller in the CALLERS field and reads a different
// FORMAT field for each. varscan2 uses DP4, strelka2 indels use TAR/TIR, GATK
// and everything else use AD, and there is a negative-AD fallback that reads DP
// instead. A hand-written version agreeing with varscan2 records would disagree
// silently on every other caller, and produce VAFs that look entirely
// reasonable. So we run MOH's own script, unchanged, from CVMFS.
// ===========================================================================

process FILTER_ENSEMBLE {
    tag   { patient_id }
    label 'process_low'

    // Host execution + environment modules, for the same reason as RUN_PCGR:
    // format2pcgr.py needs cyvcf2, which lives in the mugqic python module, not
    // in pipeline.sif. A `module` directive with an inherited container gives
    // exit 127 -- the mistake already made twice in this pipeline.
    container = null
    module    { [params.python_module,
                 params.mugqic_tools_module,
                 params.bcftools_module].findAll { it }.join(':') }

    // Keep failures inspectable (see the long note in pcgr.nf).
    scratch = false

    // Cheap to recompute but awkward to chase down; persisting it also lets you
    // diff a regenerated PCGR input against a delivered one directly.
    // enabled: !workflow.stubRun -- see the note in vcf2maf.nf. A touched
    // 0-byte VCF published here would be silently handed to PCGR.
    publishDir "${params.filtered_vcf_dir}/${cohort_id}",
               mode: 'copy', overwrite: false, enabled: !workflow.stubRun

    input:
    tuple val(cohort_id), val(patient_id), path(ensemble_vcf), val(tumour_sample)

    output:
    tuple val(cohort_id), val(patient_id),
          path("${patient_id}.ensemble.somatic.vt.annot.2caller.flt.vcf.gz"), emit: vcf
    path "${patient_id}.filter_ensemble.log",                                 emit: log

    script:
    def min_callers = params.filter_ensemble_min_callers
    def expr        = params.filter_ensemble_expr
    def out2c       = "${patient_id}.ensemble.somatic.vt.annot.2caller.vcf.gz"
    def outflt      = "${patient_id}.ensemble.somatic.vt.annot.2caller.flt.vcf.gz"
    """
    set -euo pipefail
    exec > >(tee ${patient_id}.filter_ensemble.log) 2>&1

    # ------------------------------------------------------------------ #
    # WHICH COLUMN IS THE TUMOUR
    #
    # format2pcgr.py's --tumor_name is documented as "ordinal location of
    # tumor in vcf", but the code does sample_list.index(args.tumor_name):
    # it wants the sample NAME. Passing 0 or 1 raises ValueError.
    #
    # Getting this backwards is the worst available failure: the script would
    # happily annotate the NORMAL's depth and VAF as TDP/TVAF, the bcftools
    # filter would still pass a plausible number of variants, and every
    # downstream figure would be built on the wrong sample with nothing
    # visibly wrong. So resolve it explicitly and refuse to guess.
    #
    # In this collection tumour is index 0 and normal index 1 -- the reverse
    # of the usual convention -- which is exactly why it is not hardcoded.
    # ------------------------------------------------------------------ #
    bcftools query -l ${ensemble_vcf} > .samples.txt
    echo "[filter] samples in VCF:"; sed 's/^/           /' .samples.txt

    TUMOUR="${tumour_sample}"
    if [ "\$TUMOUR" = "NA" ] || [ -z "\$TUMOUR" ]; then
        # The manifest had no dna_tumour_sample. Fall back to the naming
        # convention (...DT = DNA tumour, ...DN = DNA normal), but only if it
        # identifies EXACTLY one sample.
        TUMOUR=\$(grep -E 'DT\$' .samples.txt || true)
        n=\$(printf '%s' "\$TUMOUR" | grep -c . || true)
        if [ "\$n" -ne 1 ]; then
            echo "ERROR: cannot identify the tumour sample for ${patient_id}." >&2
            echo "  manifest dna_tumour_sample : ${tumour_sample}" >&2
            echo "  samples matching /DT\$/     : \${n}" >&2
            echo "Set dna_tumour_sample in the manifest; do not let this be guessed." >&2
            exit 1
        fi
        echo "[filter] tumour resolved by naming convention: \$TUMOUR"
    else
        if ! grep -qx "\$TUMOUR" .samples.txt; then
            echo "ERROR: manifest names tumour sample '\$TUMOUR' for ${patient_id}," >&2
            echo "but the VCF contains only:" >&2
            sed 's/^/  /' .samples.txt >&2
            exit 1
        fi
        echo "[filter] tumour from manifest: \$TUMOUR"
    fi

    # ------------------------------------------------------------------ #
    # 1. Annotate + require >= N callers.
    # ------------------------------------------------------------------ #
    n_in=\$(bcftools view -H ${ensemble_vcf} | wc -l)

    # Invoked THROUGH python3, never executed directly.
    #
    # The script's shebang is "#! python3" with CRLF line endings, so the
    # kernel looks for an interpreter literally named "python3\r" and fails:
    #     bad interpreter: No such file or directory
    # It is on PATH via the mugqic_tools module but is not runnable as a
    # command. Python itself reads CRLF source fine, so the interpreter call
    # sidesteps both problems.
    SCRIPT=\$(command -v format2pcgr.py 2>/dev/null || true)
    if [ -z "\$SCRIPT" ]; then
        SCRIPT="\${MUGQIC_TOOLS_HOME:-}/python-tools/format2pcgr.py"
    fi
    if [ ! -f "\$SCRIPT" ]; then
        echo "ERROR: format2pcgr.py not found. Is ${params.mugqic_tools_module} loaded?" >&2
        exit 1
    fi
    echo "[filter] using \$SCRIPT"

    python3 "\$SCRIPT" \\
        --input_file  ${ensemble_vcf} \\
        --output_file ${out2c} \\
        --variant_type somatic \\
        --filter      ${min_callers} \\
        --tumor_name  "\$TUMOUR"

    n_2c=\$(bcftools view -H ${out2c} | wc -l)

    # The five tags are the entire point of this step; PCGR's --tumor_dp_tag
    # etc. reference them by name and silently filter nothing if they are
    # absent. That silence is what made the original bug invisible for a full
    # run, so assert here instead of discovering it three processes later.
    # Header written once, then grepped -- NOT `bcftools view -h | grep -q`.
    # grep -q exits on first match, the writer gets SIGPIPE, and pipefail turns
    # a successful match into a failed pipeline. It happens to work here because
    # a VCF header fits in the pipe buffer, but the identical pattern in
    # RUN_PCGR (on a 98 MB file) reported every output as missing its CSQ field.
    bcftools view -h ${out2c} > .hdr.txt
    missing=""
    for t in TAL TDP TVAF NDP NVAF; do
        grep -q "^##INFO=<ID=\$t," .hdr.txt || missing="\$missing \$t"
    done
    if [ -n "\$missing" ]; then
        echo "ERROR: format2pcgr.py did not add:\$missing" >&2
        exit 1
    fi

    # ------------------------------------------------------------------ #
    # 2. Depth / VAF filter -- verbatim from [filter_ensemble]
    #    somatic_filter_options, and confirmed against the
    #    ##bcftools_viewCommand recorded in the delivered PCGR VCFs.
    # ------------------------------------------------------------------ #
    bcftools view -Oz -i'${expr}' -o ${outflt} ${out2c}
    tabix -p vcf ${outflt}

    n_flt=\$(bcftools view -H ${outflt} | wc -l)
    printf '[filter] %s: %s raw -> %s (>=%s callers) -> %s (depth/VAF)\\n' \\
        "${patient_id}" "\$n_in" "\$n_2c" "${min_callers}" "\$n_flt"

    # An empty result is not a valid answer. It means the CALLERS field was
    # absent, or the tags were computed from the wrong sample column, or the
    # thresholds do not suit this data -- all of which must stop the run rather
    # than remove a patient from the cohort without saying so.
    if [ "\$n_flt" -lt 1 ]; then
        echo "ERROR: filtering left 0 variants for ${patient_id}." >&2
        exit 1
    fi

    # Delivered patients in the cohorts seen so far land near ~24k here. An
    # order of magnitude out CAN mean the filters did not apply as they did in
    # the original runs -- but it can equally be real: melanoma, MSI and POLE
    # tumours legitimately carry hundreds of thousands of somatic variants.
    # So warn, never fail, and say both possibilities out loud.
    if [ "\$n_flt" -gt 200000 ]; then
        echo "WARNING: \$n_flt variants survived filtering for ${patient_id}." >&2
        echo "  Cohorts seen so far land near ~24k. Hypermutated tumour types" >&2
        echo "  (melanoma, MSI, POLE) reach this range legitimately -- check the" >&2
        echo "  tumour type before treating it as a filtering failure. Tags:" >&2

        # SIGPIPE-SAFE, and this is why it matters here.
        #
        # The previous form was:
        #     bcftools view -H ... | head -1 | cut -f8 | tr ... | grep ...
        # `head -1` closes the pipe after one line, bcftools is killed with
        # SIGPIPE, and under `set -euo pipefail` that 141 becomes the task's
        # exit status. So the diagnostic written FOR hypermutated samples
        # killed exactly the hypermutated samples -- and only them, which is
        # why it went unnoticed until a melanoma cohort crossed the threshold.
        #
        # pipefail off for this one line: the awk exit status is what matters,
        # and a diagnostic must never be able to fail the step it describes.
        set +o pipefail
        bcftools view -H ${outflt} 2>/dev/null \\
          | awk 'NR==1 {
                   n = split(\$8, a, ";")
                   for (i = 1; i <= n; i++)
                     if (a[i] ~ /^(TAL|TDP|TVAF|NDP|NVAF)=/) print "  " a[i]
                   exit
                 }' >&2 || true
        set -o pipefail
    fi
    """

    stub:
    """
    touch ${patient_id}.ensemble.somatic.vt.annot.2caller.flt.vcf.gz
    touch ${patient_id}.filter_ensemble.log
    """
}

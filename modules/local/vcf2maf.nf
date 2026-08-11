// ===========================================================================
// modules/local/vcf2maf.nf -- regenerate gene-annotated MAFs from PCGR VCFs
//
// WHY THIS EXISTS
// ---------------
// The MAFs distributed inside the MoHQ collection are NOT gene-annotated:
// Hugo_Symbol is "Unknown" for essentially every row. Any oncoplot, gene
// summary, TMB-by-gene or mutual-exclusivity analysis built on them is
// meaningless. The annotation has to be regenerated from the PCGR VCFs.
//
// This is your run_vcf2maf_array.sh, moved into the pipeline. What that buys:
//
//  * -resume. The array script recomputes everything on every invocation.
//    VEP is the most expensive step in the whole pipeline; at 4,500 patients
//    that is the difference between hours and days on a re-run.
//  * No hardcoded --array bounds. `#SBATCH --array=0-215  # Adjust this number`
//    is a manual step that silently truncates a cohort if you forget: with 216
//    in the array and 230 VCFs, 14 patients are dropped with no error. Nextflow
//    derives the count from the channel.
//  * Per-patient failure isolation and retry, rather than one array element
//    failing unnoticed among hundreds of log files.
//  * Published to a persistent directory, so the expensive output is reused
//    across pipeline versions rather than living in a scratch folder.
// ===========================================================================

process VCF2MAF {
    tag   { "${patient_id}" }
    label 'process_vep'

    // Run INSIDE the VEP image, not pipeline.sif. pipeline.sif has R; this step
    // needs perl + vep. Running here means `vep` is already on PATH, so no
    // --vep-path wrapper is needed and there is no nested-container problem
    // (apptainer inside apptainer generally does not work).
    container { workflow.stubRun ? null : (params.vep_container ?: params.container) }

    // Persistent, content-addressed by patient: the single most expensive
    // artefact in the pipeline should survive work-directory cleanup.
    // enabled: !workflow.stubRun -- STUB OUTPUT MUST NEVER BE PUBLISHED.
    //
    // The stub writes a one-row MAF containing a fake TP53 call. With
    // `overwrite: false`, that file then wins forever: COLLECT_DERIVED_MAFS
    // finds it, skips conversion, and the patient enters every downstream
    // analysis with a single invented mutation.
    //
    // MoHQ-HM-19-1 sat in derived_mafs/ with 1 variant while its cohort-mates
    // had ~9,000-38,000, put there by a `-stub-run` validation and reused by
    // every real run afterwards. `-stub-run` is advertised in this repo's
    // README as the safe thing to do after any edit; it was quietly corrupting
    // the most expensive artefact the pipeline produces.
    publishDir "${params.derived_maf_dir}/${cohort_id}",
               mode: 'copy', overwrite: false, enabled: !workflow.stubRun

    input:
    tuple val(cohort_id), val(patient_id), path(vcf)
    // val, NOT path, for both references.
    //
    // `path ref_fasta` symlinks ONLY the .fa.gz into the task directory, leaving
    // its .fai/.gzi indexes behind -- and samtools faidx needs them adjacent.
    // These are large, read-only reference files on a shared filesystem; staging
    // buys nothing and breaks index adjacency. Pass absolute paths instead.
    val vep_cache       // unused when params.inhibit_vep is true
    val ref_fasta       // ALWAYS required: vcf2maf uses samtools faidx for
                        // reference alleles whether or not VEP runs

    output:
    tuple val(cohort_id), val(patient_id), path("${patient_id}.maf"), emit: maf
    path "${patient_id}.vcf2maf.log",                                 emit: log

    script:
    // --filter-vcf / --max-filter-ac were REMOVED from vcf2maf. They existed to
    // flag common variants against an ExAC VCF; ExAC is superseded by gnomAD and
    // the options went with it. Older versions took `--filter-vcf 0` to disable
    // the behaviour, which is what this used to pass -- and a current vcf2maf
    // rejects the flag outright:
    //     Unknown option: filter-vcf
    //
    // Nothing replaces it, because nothing needs to: these are somatic calls
    // that PCGR (or the ensemble caller) has already filtered, and common-variant
    // removal is not this step's job. If a future vcf2maf reinstates it, set
    // params.vcf2maf_filter_vcf and restore the flag here.
    if (params.vcf2maf_filter_vcf) {
        error "params.vcf2maf_filter_vcf is set, but current vcf2maf has no " +
              "--filter-vcf option (removed with ExAC support). Unset it."
    }

    // --------------------------------------------------------------------- //
    // TWO MODES.
    //
    // inhibit_vep = true  (STRONGLY PREFERRED where it works)
    //   PCGR runs VEP itself, so its output VCF already carries a CSQ INFO
    //   field with the full annotation. `--inhibit-vep` makes vcf2maf parse
    //   that existing annotation instead of re-running VEP. This:
    //     * removes the VEP cache requirement entirely (~30 GB not copied)
    //     * turns a multi-minute per-patient VEP run into seconds
    //     * uses exactly the annotation PCGR used, so your MAFs agree with
    //       the PCGR HTML reports rather than being a second opinion
    //   Verify first that the CSQ field is present:
    //     zcat <pcgr_vcf> | grep -m1 '^##INFO=<ID=CSQ'
    //
    // inhibit_vep = false
    //   Re-annotate from scratch. Needs --vep_cache and --ref_fasta.
    // --------------------------------------------------------------------- //
    // --vep-path is only needed when vep is NOT already on PATH. Running inside
    // the VEP image it is, so leave vep_path unset in that case.
    // --ref-fasta is passed in BOTH modes. vcf2maf calls samtools faidx on it to
    // read reference alleles, independently of VEP. Omitting it under
    // --inhibit-vep does not disable the requirement; it just falls back to
    // vcf2maf's hardcoded default, which is a GRCh37 path:
    //     ERROR: Provided --ref-fasta is missing or empty:
    //            ~/.vep/homo_sapiens/112_GRCh37/Homo_sapiens.GRCh37...fa.gz
    // Note --ncbi-build does NOT influence that default.
    def annotation_args = params.inhibit_vep
        ? "--inhibit-vep --ref-fasta ${ref_fasta}"
        : "--vep-data ${vep_cache} --vep-forks ${task.cpus} --ref-fasta ${ref_fasta}" +
          " --cache-version ${params.vep_cache_version}" +
          (params.vep_path ? " --vep-path ${params.vep_path}" : '')

    // vcf2maf shells out to `samtools faidx` to read reference alleles. The
    // Ensembl VEP image does NOT ship the samtools binary -- VEP itself uses the
    // Bio::DB::HTS perl bindings instead -- so inside that container you get:
    //     ERROR: Please install samtools on your PATH, or specify --samtools-exec
    // Point params.samtools_exec at a binary visible INSIDE the container
    // (a CVMFS path works if /cvmfs is bind-mounted; check before assuming).
    // Both checks are UNCONDITIONAL in vcf2maf -- they run before the
    // --inhibit-vep branch, so turning VEP off does not avoid either:
    //     295: ( $samtools ) or die "ERROR: Please install samtools ..."
    //     299: ( $tabix )    or die "ERROR: Please install tabix ..."
    // (liftOver at 308 is only reached with --remap-chain, which we never pass.)
    //
    // The VEP image ships htslib under /opt/vep/src/htslib, which usually
    // contains tabix and bgzip but NOT samtools. Set these to paths valid
    // INSIDE the container; /cvmfs paths need `-B /cvmfs` in the profile.
    def samtools_arg = params.samtools_exec ? "--samtools-exec ${params.samtools_exec}" : ''
    def tabix_arg    = params.tabix_exec    ? "--tabix-exec ${params.tabix_exec}"       : ''
    """
    set -euo pipefail

    # PCGR VCFs use 'chr1'-style names; VEP caches and reference FASTAs use
    # Ensembl-style '1'. Preserved from the original array script.
    gunzip -c ${vcf} \\
        | sed -e 's/^chr//' -e 's/^M\\t/MT\\t/' \\
        > ${patient_id}.tmp.vcf

    if [ "${params.inhibit_vep}" = "true" ]; then
        if ! grep -q '^##INFO=<ID=CSQ' ${patient_id}.tmp.vcf; then
            echo "ERROR: inhibit_vep is set but this VCF has no CSQ INFO field," >&2
            echo "so there is no existing VEP annotation to parse." >&2
            echo "Set inhibit_vep: false and supply vep_cache + ref_fasta." >&2
            exit 1
        fi
    fi

    # --tumor-id is set to the canonical patient_id, so the MAF barcode needs
    # no downstream harmonisation at all.
    perl ${params.vcf2maf_path} \\
        --input-vcf   ${patient_id}.tmp.vcf \\
        --output-maf  ${patient_id}.maf \\
        --tumor-id    ${patient_id} \\
        --ncbi-build  ${params.genome_build} \\
        ${annotation_args} ${samtools_arg} ${tabix_arg} \\
        2>&1 | tee ${patient_id}.vcf2maf.log

    rm -f ${patient_id}.tmp.vcf

    # Did annotation actually happen?
    #
    # The failure this guards against is the original MoHQ problem: MAFs where
    # Hugo_Symbol is "Unknown" for EVERY row, which silently produce an empty
    # oncoplot rather than an error.
    #
    # An earlier version failed above 50% Unknown. That is wrong for whole-genome
    # data: most somatic variants in WGS are intergenic, and vcf2maf correctly
    # writes Unknown for those. A real WGS patient here came out at 54% Unknown
    # with 23,605 variants -- perfectly healthy, and rejected.
    #
    # The discriminating signal is the number of DISTINCT gene symbols. Failed
    # annotation gives 0; working annotation on WGS gives thousands. The
    # fraction is still reported, and still checked, but only at a level that
    # cannot be reached by ordinary intergenic variants.
    total=\$(awk -F'\\t' '!/^#/ && \$1!="Hugo_Symbol"' ${patient_id}.maf | wc -l)
    if [ "\$total" -gt 0 ]; then
        unknown=\$(awk -F'\\t' '!/^#/ && \$1=="Unknown"' ${patient_id}.maf | wc -l)
        genes=\$(awk -F'\\t' '!/^#/ && \$1!="Hugo_Symbol" && \$1!="Unknown" && \$1!="" {print \$1}' \\
                 ${patient_id}.maf | sort -u | wc -l)
        frac=\$(awk -v u="\$unknown" -v t="\$total" 'BEGIN{printf "%.3f", u/t}')
        echo "[vcf2maf] ${patient_id}: \$total variants, \$unknown Unknown (\$frac), \$genes distinct genes"

        if [ "\$genes" -lt ${params.min_annotated_genes} ]; then
            echo "ERROR: only \$genes distinct gene symbols (need >= ${params.min_annotated_genes})." >&2
            echo "Annotation did not happen. With --inhibit-vep, check the input VCF" >&2
            echo "really carries a CSQ INFO field; otherwise check the VEP cache." >&2
            exit 1
        fi
        if awk -v f="\$frac" 'BEGIN{exit !(f > ${params.max_unknown_frac})}'; then
            echo "ERROR: \$frac of rows have Hugo_Symbol=Unknown (limit ${params.max_unknown_frac})." >&2
            exit 1
        fi
    fi
    """

    stub:
    """
    printf 'Hugo_Symbol\\tTumor_Sample_Barcode\\nTP53\\t${patient_id}\\n' > ${patient_id}.maf
    touch ${patient_id}.vcf2maf.log
    """
}


// ---------------------------------------------------------------------------
// Reuse already-converted MAFs instead of re-running VEP.
//
// VCF2MAF is by far the most expensive step. If a MAF already exists in
// params.derived_maf_dir from a previous run, take it and skip conversion
// entirely -- this survives even a cleared Nextflow work directory.
// ---------------------------------------------------------------------------
process COLLECT_DERIVED_MAFS {
    tag   { cohort_id }
    label 'process_single'

    input:
    tuple val(cohort_id), val(patient_ids)

    output:
    tuple val(cohort_id), path("existing_mafs.txt"), emit: existing

    script:
    """
    touch existing_mafs.txt
    for p in ${patient_ids.join(' ')}; do
        f="${params.derived_maf_dir}/${cohort_id}/\${p}.maf"
        if [ -s "\$f" ]; then echo "\$f" >> existing_mafs.txt; fi
    done
    echo "[collect] \$(wc -l < existing_mafs.txt) of ${patient_ids.size()} MAFs already converted" >&2
    """

    stub:
    """
    touch existing_mafs.txt
    """
}

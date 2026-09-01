// ===========================================================================
// modules/local/manifest.nf
//
// Ingest: walk a cohort's delivery tree, build the canonical sample manifest
// and the data-completeness matrix, and assemble the cohort .seg.
//
// This replaces the previous hand-staging into
//   ../results/all_cohorts_analysis_ready/<cohort>/
// which sat outside Nextflow entirely, so nothing recorded how those files
// were produced and `-resume` could not tell when they changed.
// ===========================================================================

process BUILD_MANIFEST {
    tag   { cohort_id }
    label 'process_single'

    // enabled: !workflow.stubRun -- a -stub-run writes placeholder files;
    // publishing them overwrites real results. A stub manifest (2 fake
    // patients) replaced MoHQ-MU-16's real 102-patient manifest this way,
    // and the preflight then validated the cohort against it.
    publishDir "${params.outdir}/${cohort_id}/manifest", mode: 'copy', enabled: !workflow.stubRun

    input:
    tuple val(cohort_id), path(cohort_dir)

    output:
    tuple val(cohort_id), path("${cohort_id}.manifest.tsv"),     emit: manifest
    tuple val(cohort_id), path("${cohort_id}.completeness.tsv"), emit: completeness
    tuple val(cohort_id), path("${cohort_id}.samples.tsv"),      emit: samples
    tuple val(cohort_id), path("${cohort_id}.ambiguities.tsv"),  emit: ambiguities
    path "versions.yml",                                          emit: versions

    script:
    // Clinical metadata is joined separately by MERGE_METADATA so that the
    // manifest builder stays a pure function of the delivery tree.
    //
    // The cache lives OUTSIDE the task work directory on purpose: it must
    // survive across runs to be useful. It is keyed on subdirectory mtimes, so
    // a stale cache cannot produce a wrong manifest -- only a slower one if
    // deleted. Nextflow's own -resume still gates whether this task runs at all;
    // the cache helps when the task DOES rerun (new patients in a cohort).
    def cache_arg = params.use_manifest_cache
        ? "--cache ${params.manifest_cache_dir}/${cohort_id}.json"
        : ''

    // PCGR VCFs delivered outside the patient folders (see nextflow.config).
    // A RELATIVE value resolves inside the staged cohort directory, which is
    // where such batch directories live -- MoHQ-MU-16/20231122_vepvcfs. An
    // absolute path is passed through untouched, for a batch kept elsewhere.
    //
    // Nextflow stages cohort_dir as a whole, so the batch directory comes with
    // it and needs no separate input channel.
    def ext_arg = ''
    if (params.external_pcgr_dir) {
        def p = params.external_pcgr_dir.toString()
        ext_arg = p.startsWith('/') ? "--external-pcgr-dir ${p}"
                                    : "--external-pcgr-dir ${cohort_dir}/${p}"
    }
    """
    # Drop the MATLAB Compiler Runtime from the library path.
    #
    # GISTIC ships MCR v83, whose bin/glnxa64 contains libexpat.so.1.5.0 from
    # ~2012. The container's %environment put that first on LD_LIBRARY_PATH for
    # every process, so python3 loaded the ancient expat and died with
    #     symbol lookup error: undefined symbol: XML_SetHashSalt
    #
    # This has to happen INSIDE the container: apptainer sources %environment on
    # startup, which is after Nextflow's `beforeScript` has already run on the
    # host. R is unaffected, so only this python step needs it.
    #
    # Rebuilding with the current pipeline.def makes this redundant -- there,
    # MCR is scoped to the GISTIC process alone.
    unset LD_LIBRARY_PATH

    mkdir -p ${params.manifest_cache_dir}

    build_manifest.py \\
        --root      ${cohort_dir} \\
        --cohort-id ${cohort_id} \\
        --outdir    . \\
        --prefix    ${cohort_id} ${cache_arg} ${ext_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    // An EMPTY manifest would make every downstream channel empty, so the stub
    // run would validate syntax and nothing else. Emit two realistic patients
    // instead -- one with a delivered PCGR VCF, one 'regenerable' -- so the
    // joins, filters, the genome-build guard and the gap-fill routing are all
    // actually exercised.
    //
    // Asset paths must be UNIQUE PER PATIENT. An earlier version pointed every
    // asset at /dev/null, which exists (so checkIfExists passes) but means that
    // any process collecting files from several patients stages two inputs both
    // named "null":
    //     input file name collision -- multiple input files named: null
    // So we create per-patient dummy files here and reference them by absolute
    // path. They live in this task's work directory, which persists as long as
    // work/ does -- exactly as long as the stub run needs them.
    """
    D=\$PWD/stub_assets
    mkdir -p \$D
    for p in ${cohort_id}-1 ${cohort_id}-2; do
        for a in pcgrvcf maf tiers cna ensvcf cnvkit expr annofuse linxfus linxdrv keymetrics ini; do
            : > \$D/\${p}.\${a}
        done
    done

    cat > ${cohort_id}.manifest.tsv <<TSV
collection	institution	cohort_id	patient_id	patient_dir	dna_normal_sample	dna_tumour_sample	rna_tumour_sample	n_dna_tumour_samples	n_dna_normal_samples	n_rna_tumour_samples	multi_tumour	analysis_ready	missing_core	genome_build	genome_build_source	genome_build_from_filename	pcgr_version	genpipes_version	r_version	pcgr_vcf_status	mutation_vcf	mutation_vcf_source	all_samples	pcgr_vcf	somatic_maf	somatic_tiers	cna_segments	ensemble_somatic_vcf	cnvkit_vcf	expression_genes	anno_fuse	linx_fusion	key_metrics	tumourpair_ini
MoHQ	${cohort_id.split('-')[1]}	${cohort_id}	${cohort_id}-1	\$PWD	NA	NA	NA	1	1	1	no	yes	NA	${params.genome_build}	stub	${params.genome_build}	mugqic/pcgr/1.0.3	NA	NA	delivered	\$D/${cohort_id}-1.pcgrvcf	pcgr_delivered	NA	\$D/${cohort_id}-1.pcgrvcf	\$D/${cohort_id}-1.maf	\$D/${cohort_id}-1.tiers	\$D/${cohort_id}-1.cna	\$D/${cohort_id}-1.ensvcf	\$D/${cohort_id}-1.cnvkit	\$D/${cohort_id}-1.expr	\$D/${cohort_id}-1.annofuse	\$D/${cohort_id}-1.linxfus	\$D/${cohort_id}-1.keymetrics	\$D/${cohort_id}-1.ini
MoHQ	${cohort_id.split('-')[1]}	${cohort_id}	${cohort_id}-2	\$PWD	NA	NA	NA	1	1	1	no	yes	NA	${params.genome_build}	stub	${params.genome_build}	mugqic/pcgr/1.4.1	NA	NA	regenerable	\$D/${cohort_id}-2.ensvcf	ensemble_pending_gapfill	NA	NA	\$D/${cohort_id}-2.maf	\$D/${cohort_id}-2.tiers	\$D/${cohort_id}-2.cna	\$D/${cohort_id}-2.ensvcf	\$D/${cohort_id}-2.cnvkit	\$D/${cohort_id}-2.expr	\$D/${cohort_id}-2.annofuse	\$D/${cohort_id}-2.linxfus	\$D/${cohort_id}-2.keymetrics	\$D/${cohort_id}-2.ini
TSV

    printf 'cohort_id\\tpatient_id\\tanalysis_ready\\tcna_segments\\n%s\\t%s-1\\tyes\\tyes\\n%s\\t%s-2\\tyes\\tyes\\n' \\
        "${cohort_id}" "${cohort_id}" "${cohort_id}" "${cohort_id}" > ${cohort_id}.completeness.tsv

    printf 'cohort_id\\tpatient_id\\tsample_id\\tsample_num\\taliquot\\tseqtype\\tmolecule\\ttissue\\tis_primary_tumour\\n' \\
        > ${cohort_id}.samples.tsv
    printf 'patient_id\\tasset\\tchosen\\tn_alternatives\\talternatives\\treason\\n' \\
        > ${cohort_id}.ambiguities.tsv
    touch versions.yml
    """
}


process MERGE_METADATA {
    tag   { cohort_id }
    label 'process_single'

    // enabled: !workflow.stubRun -- a -stub-run writes placeholder files;
    // publishing them overwrites real results. A stub manifest (2 fake
    // patients) replaced MoHQ-MU-16's real 102-patient manifest this way,
    // and the preflight then validated the cohort against it.
    publishDir "${params.outdir}/${cohort_id}/manifest", mode: 'copy', enabled: !workflow.stubRun

    input:
    tuple val(cohort_id), path(manifest), path(extra_meta)

    output:
    tuple val(cohort_id), path("${cohort_id}.manifest.annotated.tsv"), emit: manifest

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(data.table))

    man <- fread("${manifest}", sep = "\\t")
    ext <- fread("${extra_meta}")

    # The join key is whatever column in the extra metadata holds patient IDs.
    # We do NOT prefix-match: an exact join or a loud failure.
    key <- intersect(c("patient_id", "Patient", "Sample", "sample_id", "ID"), names(ext))[1]
    if (is.na(key)) {
        stop("Extra metadata has no recognisable patient column. Found: ",
             paste(names(ext), collapse = ", "))
    }
    setnames(ext, key, "join_key")

    # Reduce to canonical patient IDs using the same rule as everything else.
    ext[, patient_id := sub("_(D|R|T|N)\$", "", trimws(as.character(join_key)))]
    ext[, patient_id := sub("^([A-Za-z]+-[A-Za-z]+-[A-Za-z0-9]+-[A-Za-z0-9]+)-[A-Za-z0-9]+-[0-9]+(DN|DT|RN|RT)\$",
                            "\\\\1", patient_id)]
    ext[, join_key := NULL]

    matched <- sum(man\$patient_id %in% ext\$patient_id)
    message(sprintf("[metadata] %d/%d manifest patients matched extra metadata",
                    matched, nrow(man)))
    if (matched == 0) {
        stop("No patient IDs matched between the manifest and ${extra_meta}. ",
             "Manifest e.g.: ", paste(head(man\$patient_id, 3), collapse = ", "))
    }

    out <- merge(man, ext, by = "patient_id", all.x = TRUE)
    fwrite(out, "${cohort_id}.manifest.annotated.tsv", sep = "\\t", na = "NA", quote = FALSE)
    """

    stub:
    """
    touch ${cohort_id}.manifest.annotated.tsv
    """
}


process BUILD_COHORT_SEG {
    tag   { cohort_id }
    label 'process_medium'

    // enabled: !workflow.stubRun -- a -stub-run writes placeholder files;
    // publishing them overwrites real results. A stub manifest (2 fake
    // patients) replaced MoHQ-MU-16's real 102-patient manifest this way,
    // and the preflight then validated the cohort against it.
    publishDir "${params.outdir}/${cohort_id}/cnv", mode: 'copy', enabled: !workflow.stubRun

    input:
    tuple val(cohort_id), path(manifest), path(vcfs)
    path rlib

    output:
    tuple val(cohort_id), path("${cohort_id}.seg"),         emit: seg
    tuple val(cohort_id), path("${cohort_id}.markers.txt"), emit: markers
    path "*.versions.txt",                                   emit: versions

    script:
    """
    build_cohort_seg.R \\
        --lib      ${rlib} \\
        --vcfs     ${vcfs} \\
        --manifest ${manifest} \\
        --cohort   ${cohort_id} \\
        --min_markers ${params.min_seg_markers}
    """

    stub:
    """
    touch ${cohort_id}.seg ${cohort_id}.markers.txt ${cohort_id}_build_cohort_seg.versions.txt
    """
}


process PLOT_COMPLETENESS {
    tag   { cohort_id }
    label 'process_single'

    // enabled: !workflow.stubRun -- a -stub-run writes placeholder files;
    // publishing them overwrites real results. A stub manifest (2 fake
    // patients) replaced MoHQ-MU-16's real 102-patient manifest this way,
    // and the preflight then validated the cohort against it.
    publishDir "${params.outdir}/${cohort_id}/completeness", mode: 'copy', enabled: !workflow.stubRun

    input:
    tuple val(cohort_id), path(completeness), path(manifest)
    path rlib

    output:
    tuple val(cohort_id), path("*.png"), emit: plots
    tuple val(cohort_id), path("*.tsv"), emit: tables
    path "*.versions.txt",               emit: versions

    script:
    """
    plot_completeness.R \\
        --lib           ${rlib} \\
        --completeness  ${completeness} \\
        --manifest      ${manifest} \\
        --cohort        ${cohort_id} \\
        --max_tile_rows ${params.completeness_max_tile_rows} \\
        --drop_never_present ${params.completeness_drop_never_present} \\
        --out_prefix    ${cohort_id}_completeness
    """

    stub:
    """
    touch ${cohort_id}_completeness.png ${cohort_id}_completeness.tsv ${cohort_id}_plot_completeness.versions.txt
    """
}

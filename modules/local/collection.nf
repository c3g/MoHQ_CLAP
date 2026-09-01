// ===========================================================================
// modules/local/collection.nf -- collection tier (all ~65 cohorts together)
//
// These processes consume per-cohort SUMMARY tables only, never raw MAFs or
// VCFs. That is deliberate: it keeps the collection tier cheap and constant-cost
// as the collection grows from 4,500 patients to whatever comes next.
// ===========================================================================

process COLLECTION_ROLLUP {
    tag   'collection'
    label 'process_medium'

    // enabled: !workflow.stubRun -- a -stub-run writes placeholder files;
    // publishing them overwrites real results. A stub manifest (2 fake
    // patients) replaced MoHQ-MU-16's real 102-patient manifest this way,
    // and the preflight then validated the cohort against it.
    publishDir "${params.outdir}/_collection", mode: 'copy', enabled: !workflow.stubRun

    input:
    path manifests,    stageAs: 'manifests/*'
    path completeness, stageAs: 'completeness/*'
    path rlib
    path cohort_map

    output:
    path "*.png", emit: plots
    path "*.tsv", emit: tables
    path "*.versions.txt", emit: versions

    script:
    def cmap = cohort_map.name != 'NO_COHORT_MAP' ? "--cohort_map ${cohort_map}" : ''
    """
    collection_rollup.R \\
        --lib          ${rlib} \\
        --manifests    manifests/* \\
        --completeness completeness/* ${cmap} \\
        --collection   ${params.collection_name} \\
        --out_prefix   ${params.collection_name}_collection
    """

    stub:
    """
    touch ${params.collection_name}_collection.png ${params.collection_name}_collection.tsv
    touch collection_rollup.versions.txt
    """
}


process BATCH_EFFECT_ANALYSIS {
    tag   'collection'
    label 'process_medium'

    // enabled: !workflow.stubRun -- a -stub-run writes placeholder files;
    // publishing them overwrites real results. A stub manifest (2 fake
    // patients) replaced MoHQ-MU-16's real 102-patient manifest this way,
    // and the preflight then validated the cohort against it.
    publishDir "${params.outdir}/_collection/batch_effects", mode: 'copy', enabled: !workflow.stubRun

    input:
    path manifests, stageAs: 'manifests/*'
    path metrics,   stageAs: 'metrics/*'
    path rlib
    path cohort_map

    output:
    path "*.png", emit: plots
    path "*.tsv", emit: tables
    path "*.versions.txt", emit: versions

    script:
    // The whole point of this process is the crossed cancer_type x institution
    // design, which cannot be identified without the cohort map.
    """
    batch_effect_analysis.R \\
        --lib          ${rlib} \\
        --manifests    manifests/* \\
        --metrics      metrics/* \\
        --cohort_map   ${cohort_map} \\
        --min_per_cell ${params.batch_min_per_cell} \\
        --out_prefix   ${params.collection_name}_batch
    """

    stub:
    """
    touch ${params.collection_name}_batch.png ${params.collection_name}_batch.tsv
    touch batch_effect_analysis.versions.txt
    """
}

// ===========================================================================
// modules/local/analysis.nf -- cohort-scale analysis processes
//
// Design notes on what changed from the original main.nf:
//
//  * Inputs are FILE LISTS, not directories. Passing `path maf_dir` meant
//    Nextflow hashed only the directory name for -resume, so adding or
//    changing a MAF inside it did not invalidate the cached task and you got
//    silently stale results. Staging individual files fixes that.
//
//  * No `when:` directives. The original had both `when: params.run_x` on the
//    process AND `if (params.run_x)` in the workflow. Conditional logic lives
//    in the workflow only -- one place, easier to reason about.
//
//  * Every process declares a `label`, so resources come from conf/base.config
//    instead of every task inheriting the same global 4h / 32 GB.
//
//  * `stub:` blocks everywhere, so `-stub-run` can validate the whole DAG in
//    seconds without touching real data. Use this on every edit.
// ===========================================================================

process ONCOPRINT {
    tag   { cohort_id }
    label 'process_high_memory'

    // enabled: !workflow.stubRun -- a -stub-run writes placeholder files;
    // publishing them overwrites real results. A stub manifest (2 fake
    // patients) replaced MoHQ-MU-16's real 102-patient manifest this way,
    // and the preflight then validated the cohort against it.
    publishDir "${params.outdir}/${cohort_id}/mutations", mode: 'copy', enabled: !workflow.stubRun

    input:
    tuple val(cohort_id), path(manifest), path(mafs), path(cna)
    path rlib
    // Staged ONLY so that editing it invalidates the -resume cache.
    //
    // Nextflow hashes a task from its script text, its declared path inputs and
    // its container. Scripts in bin/ are put on PATH but are NOT inputs, so
    // editing run_oncoprint.R alone leaves the cached result valid and -resume
    // silently returns the OLD figure. That is the most dangerous kind of stale:
    // the run succeeds and the output looks current.
    path rscript
    path gene_panel

    output:
    tuple val(cohort_id), path("*_oncoprint.png"),  emit: oncoprint
    tuple val(cohort_id), path("*_tmb.png"),        emit: tmb_plot
    tuple val(cohort_id), path("*.tsv"),            emit: tables
    tuple val(cohort_id), path("*.png"),            emit: all_plots
    path "*.versions.txt",                          emit: versions

    script:
    def cna_arg = cna ? "--cna ${cna}" : ''
    """
    run_oncoprint.R \\
        --lib         ${rlib} \\
        --mafs        ${mafs} ${cna_arg} \\
        --manifest    ${manifest} \\
        --cohort      ${cohort_id} \\
        --amp         ${params.amp_threshold} \\
        --del         ${params.del_threshold} \\
        --top         ${params.oncoplot_top_genes} \\
        --callable_mb ${params.callable_mb} \\
        --exclude_flags ${params.exclude_flags} \\
        --extra_flags "${params.extra_flags}" \\
        ${params.gene_panel ? "--gene_panel ${file(params.gene_panel).name}" : ''} \\
        ${params.gtf ? "--gtf ${params.gtf}" : ''} \\
        --oncodrive         ${params.run_oncodrive} \\
        --oncodrive_min_mut ${params.oncodrive_min_mut} \\
        --out_prefix  ${cohort_id}_oncoprint
    """

    stub:
    """
    touch ${cohort_id}_oncoprint.png ${cohort_id}_oncoprint_tmb.png \\
          ${cohort_id}_oncoprint.tsv ${cohort_id}_oncoprint.versions.txt
    """
}


process CNV_FREQUENCY {
    tag   { cohort_id }
    label 'process_medium'

    // enabled: !workflow.stubRun -- a -stub-run writes placeholder files;
    // publishing them overwrites real results. A stub manifest (2 fake
    // patients) replaced MoHQ-MU-16's real 102-patient manifest this way,
    // and the preflight then validated the cohort against it.
    publishDir "${params.outdir}/${cohort_id}/cnv", mode: 'copy', enabled: !workflow.stubRun

    input:
    tuple val(cohort_id), path(manifest), path(seg)
    path rlib

    output:
    tuple val(cohort_id), path("*_frequency.png"), emit: plot
    tuple val(cohort_id), path("*_frequency.tsv"), emit: table
    path "*.versions.txt",                          emit: versions

    script:
    def sexchr = params.include_sex_chr ? '--include_sex_chr' : ''
    """
    run_cnv_frequency.R \\
        --lib        ${rlib} \\
        --seg        ${seg} \\
        --manifest   ${manifest} \\
        --cohort     ${cohort_id} \\
        --amp        ${params.amp_threshold} \\
        --del        ${params.del_threshold} \\
        --bin_size   ${params.cnv_bin_size} ${sexchr} \\
        --out_prefix ${cohort_id}_cnv_frequency
    """

    stub:
    """
    touch ${cohort_id}_cnv_frequency.png ${cohort_id}_cnv_frequency.tsv ${cohort_id}_cnv_frequency.versions.txt
    """
}


process CNV_BURDEN {
    tag   { cohort_id }
    label 'process_low'

    // enabled: !workflow.stubRun -- a -stub-run writes placeholder files;
    // publishing them overwrites real results. A stub manifest (2 fake
    // patients) replaced MoHQ-MU-16's real 102-patient manifest this way,
    // and the preflight then validated the cohort against it.
    publishDir "${params.outdir}/${cohort_id}/cnv", mode: 'copy', enabled: !workflow.stubRun

    input:
    tuple val(cohort_id), path(manifest), path(seg)
    path rlib

    output:
    tuple val(cohort_id), path("*_burden.png"), emit: plot
    tuple val(cohort_id), path("*.tsv"),        emit: tables
    tuple val(cohort_id), path("*.png"),        emit: all_plots
    path "*.versions.txt",                      emit: versions

    script:
    """
    run_fga_burden.R \\
        --lib             ${rlib} \\
        --seg             ${seg} \\
        --manifest        ${manifest} \\
        --cohort          ${cohort_id} \\
        --amp             ${params.amp_threshold} \\
        --del             ${params.del_threshold} \\
        --min_profiled_mb ${params.min_profiled_mb} \\
        --out_prefix      ${cohort_id}_cnv_burden
    """

    stub:
    """
    touch ${cohort_id}_cnv_burden.png ${cohort_id}_cnv_burden.tsv ${cohort_id}_cnv_burden.versions.txt
    """
}


process EXPRESSION_PCA {
    tag   { cohort_id }
    label 'process_high_memory'

    // enabled: !workflow.stubRun -- a -stub-run writes placeholder files;
    // publishing them overwrites real results. A stub manifest (2 fake
    // patients) replaced MoHQ-MU-16's real 102-patient manifest this way,
    // and the preflight then validated the cohort against it.
    publishDir "${params.outdir}/${cohort_id}/expression", mode: 'copy', enabled: !workflow.stubRun

    input:
    tuple val(cohort_id), path(manifest), path(expr)
    path rlib

    output:
    tuple val(cohort_id), path("*_pca.html"), emit: html
    tuple val(cohort_id), path("*.png"),      emit: plots
    tuple val(cohort_id), path("*_scores.tsv"), emit: scores
    path "*.versions.txt",                    emit: versions

    script:
    """
    run_pca.R \\
        --lib           ${rlib} \\
        --expr          ${expr} \\
        --manifest      ${manifest} \\
        --cohort        ${cohort_id} \\
        --norm_method   ${params.pca_norm_method} \\
        --top_var_genes ${params.pca_top_var_genes} \\
        --colour_by     ${params.pca_colour_by} \\
        --out_prefix    ${cohort_id}_expression_pca
    """

    stub:
    """
    touch ${cohort_id}_expression_pca.html ${cohort_id}_expression_pca.png \\
          ${cohort_id}_expression_pca_scores.tsv ${cohort_id}_expression_pca.versions.txt
    """
}


process COMPARATIVE {
    tag   { cohort_id }
    label 'process_medium'

    // enabled: !workflow.stubRun -- a -stub-run writes placeholder files;
    // publishing them overwrites real results. A stub manifest (2 fake
    // patients) replaced MoHQ-MU-16's real 102-patient manifest this way,
    // and the preflight then validated the cohort against it.
    publishDir "${params.outdir}/${cohort_id}/comparative", mode: 'copy', enabled: !workflow.stubRun

    input:
    tuple val(cohort_id), path(manifest), path(mafs), path(seg)
    path rlib
    path gtf

    output:
    tuple val(cohort_id), path("*.png"), emit: plots
    tuple val(cohort_id), path("*.tsv"), emit: tables, optional: true
    path "*.versions.txt",               emit: versions

    script:
    def sexchr = params.include_sex_chr ? '--include_sex_chr' : ''
    def genes  = params.target_genes ?: ''
    """
    run_comparative.R \\
        --lib          ${rlib} \\
        --mafs         ${mafs} \\
        --seg          ${seg} \\
        --manifest     ${manifest} \\
        --gtf          ${gtf} \\
        --group_col    ${params.comparative_group_col} \\
        --genes        "${genes}" \\
        --n_genes      ${params.comparative_n_genes} \\
        --cohort       ${cohort_id} \\
        --amp          ${params.amp_threshold} \\
        --del          ${params.del_threshold} \\
        --min_group_n  ${params.comparative_min_group_n} ${sexchr} \\
        --exclude_flags ${params.exclude_flags} \\
        --min_profiled_mb ${params.min_profiled_mb} \\
        --out_prefix   ${cohort_id}_comparative
    """

    stub:
    """
    touch ${cohort_id}_comparative_mirror.png ${cohort_id}_comparative.versions.txt
    """
}


process RECURRENT_FUSIONS {
    tag   { cohort_id }
    label 'process_low'

    // enabled: !workflow.stubRun -- a -stub-run writes placeholder files;
    // publishing them overwrites real results. A stub manifest (2 fake
    // patients) replaced MoHQ-MU-16's real 102-patient manifest this way,
    // and the preflight then validated the cohort against it.
    publishDir "${params.outdir}/${cohort_id}/fusions", mode: 'copy', enabled: !workflow.stubRun

    input:
    tuple val(cohort_id), path(manifest), path(annofuse), path(linx)
    path rlib

    output:
    tuple val(cohort_id), path("*.png"), emit: plot
    tuple val(cohort_id), path("*.tsv"), emit: table
    path "*.versions.txt",               emit: versions

    script:
    def af = annofuse ? "--annofuse ${annofuse}"      : ''
    def lx = linx     ? "--linx_fusion ${linx}"       : ''
    """
    run_fusions.R \\
        --lib             ${rlib} ${af} ${lx} \\
        --manifest        ${manifest} \\
        --cohort          ${cohort_id} \\
        --min_recurrence  ${params.min_fusion_recurrence} \\
        --out_prefix      ${cohort_id}_fusions
    """

    stub:
    """
    touch ${cohort_id}_fusions.png ${cohort_id}_fusions.tsv ${cohort_id}_recurrent_fusions.versions.txt
    """
}


process GISTIC {
    tag   { cohort_id }
    label 'process_high_memory'

    // enabled: !workflow.stubRun -- a -stub-run writes placeholder files;
    // publishing them overwrites real results. A stub manifest (2 fake
    // patients) replaced MoHQ-MU-16's real 102-patient manifest this way,
    // and the preflight then validated the cohort against it.
    publishDir "${params.outdir}/${cohort_id}/gistic", mode: 'copy', enabled: !workflow.stubRun

    input:
    tuple val(cohort_id), path(seg), path(markers)
    path refgene

    output:
    tuple val(cohort_id), path("gistic_out/*"), emit: results
    path "gistic.log",                          emit: log

    script:
    // -broad 1 / -armpeel 1: the original set both to 0, which switches OFF
    // arm-level and broad-event analysis. Broad events are most of what GISTIC
    // is useful for at cohort scale, so they are on by default here.
    """
    mkdir -p gistic_out mcr_cache
    export MCR_CACHE_ROOT=\$PWD/mcr_cache

    # GISTIC is the ONLY step that needs the MATLAB Compiler Runtime on the
    # library path. It is scoped here rather than set container-wide because
    # MCR v83 bundles 2014-vintage system libraries (libexpat and friends) that
    # break anything else linking the modern ones -- python3 in particular.
    export LD_LIBRARY_PATH=\${MCR_LD_PATH}\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}

    gp_gistic2_from_seg \\
        -b        gistic_out \\
        -seg      ${seg} \\
        -mk       ${markers} \\
        -refgene  ${refgene} \\
        -genegistic 1 \\
        -smallmem   1 \\
        -broad      ${params.gistic_broad} \\
        -armpeel    ${params.gistic_armpeel} \\
        -conf       ${params.gistic_conf} \\
        -ta         ${params.gistic_ta} \\
        -td         ${params.gistic_td} \\
        -scent      ${params.gistic_scent} \\
        -gcm        extreme \\
        -savegene   1 \\
        > gistic.log 2>&1

    # `[ -s glob ]` does not expand a multi-match glob reliably; count instead.
    # GISTIC exits 0 on several internal failures, so the output must be checked
    # explicitly or the pipeline will happily carry on with an empty directory.
    n_lesions=\$(ls gistic_out/all_lesions.conf_*.txt 2>/dev/null | wc -l)
    if [ "\$n_lesions" -eq 0 ]; then
        echo "GISTIC produced no lesion table -- see gistic.log" >&2
        tail -50 gistic.log >&2
        exit 1
    fi
    """

    stub:
    """
    mkdir -p gistic_out && touch gistic_out/all_lesions.txt gistic.log
    """
}


process COHORT_REPORT {
    tag   { cohort_id }
    label 'process_low'

    // enabled: !workflow.stubRun -- a -stub-run writes placeholder files;
    // publishing them overwrites real results. A stub manifest (2 fake
    // patients) replaced MoHQ-MU-16's real 102-patient manifest this way,
    // and the preflight then validated the cohort against it.
    publishDir "${params.outdir}/${cohort_id}/report", mode: 'copy', enabled: !workflow.stubRun

    input:
    tuple val(cohort_id), path(manifest), path(completeness), path(figures, stageAs: 'figures/*')
    path rmd

    output:
    tuple val(cohort_id), path("${cohort_id}_report.html"), emit: report

    script:
    // The original passed nine separate positional `path` inputs and padded
    // absent ones with `ifEmpty(file("no_oncoprint.png"))` -- files that did
    // not exist on disk, so staging failed. Here every figure that exists is
    // collected into one directory and the template discovers them at render
    // time, so a missing panel degrades to a note instead of a crash.
    // The parameter record is written HERE, not passed through a channel.
    //
    // The channel version silently did not stage -- the report kept saying
    // "_No parameter record staged._" while the version files staged fine, so
    // the mechanism worked for one input and not the other and there was no
    // error either way. Writing the file in the process removes the failure
    // mode entirely: if the process ran, the record exists.
    //
    // Only parameters that change a NUMBER are recorded. Paths and toggles are
    // noise in a results document.
    // RECORD THE TOGGLES, NOT ONLY THE THRESHOLDS.
    //
    // This list was numeric settings only, so every run_* switch resolved to
    // its default in the report. Two consequences, both silent:
    //
    //   1. absent() decides between "you turned this off" and "this failed"
    //      by reading its toggle. With the toggle unrecorded it always chose
    //      the alarming wording, for all seven panels.
    //   2. The eligibility table counts gap-filled patients as having a PCGR
    //      VCF only when run_pcgr_gapfill is true. Unrecorded, it read false,
    //      so MoHQ-HM-19 reported 20 eligible for mutations after a gap-fill
    //      run that produced 44 MAFs.
    //
    // params[k] already includes command-line overrides, so `--run_pcgr_gapfill
    // true` is captured once the key is listed here.
    def reported = ['genome_build','callable_mb','amp_threshold','del_threshold',
                    'cnv_bin_size','include_sex_chr','min_seg_markers',
                    'min_profiled_mb','oncoplot_top_genes','exclude_flags',
                    'extra_flags','pca_norm_method','pca_top_var_genes',
                    'pca_colour_by','min_fusion_recurrence','inhibit_vep',
                    'vep_cache_version','min_annotated_genes','max_unknown_frac',
                    // toggles read by absent() and by the eligibility table
                    'run_pcgr_gapfill','run_oncoprint','run_cnv_frequency',
                    'run_cnv_burden','run_pca','run_fusions','run_comparative',
                    'run_completeness_plot','run_gistic','run_collection',
                    'comparative_group_col','comparative_min_group_n',
                    'comparative_n_genes','pcgr_module_fallback','extra_metadata',
                    // Which genes a figure was built on, and how they were ranked.
                    // Without these the report cannot say whether an oncoplot was
                    // restricted to a panel or ranked over everything -- the first
                    // question anyone asks of a filtered gene list.
                    'gene_panel','gtf','run_oncodrive','oncodrive_min_mut',
                    'completeness_drop_never_present']
    def rows = reported.collect { k ->
        "${k}\t${params.containsKey(k) && params[k] != null ? params[k] : 'unset'}"
    }.join('\n')
    """
    cp ${rmd} ./report.Rmd
    mkdir -p figures

    cat > figures/run_params.tsv <<'PARAMS_EOF'
param\tvalue
${rows}
PARAMS_EOF

    Rscript -e "rmarkdown::render('report.Rmd', \\
        params = list(cohort = '${cohort_id}', \\
                      manifest = '${manifest}', \\
                      completeness = '${completeness}', \\
                      figdir = 'figures', \\
                      run_params = 'figures/run_params.tsv'), \\
        output_file = '${cohort_id}_report.html')"
    """

    stub:
    """
    touch ${cohort_id}_report.html
    """
}

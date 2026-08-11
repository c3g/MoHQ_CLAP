#!/usr/bin/env nextflow
// ===========================================================================
// MoHQ cohort-scale analysis pipeline
//
// Usage
//   nextflow run . -profile narval -params-file params/MoHQ-CM-4.yaml -resume
//   nextflow run . -profile narval --input samplesheet.csv -resume   # all cohorts
//
// Always validate structural edits with a stub run first -- it exercises the
// whole DAG in seconds without touching data:
//   nextflow run . -profile test,stub -stub-run
//
// ---------------------------------------------------------------------------
// WHAT CHANGED, AND WHY
// ---------------------------------------------------------------------------
// 1. SAMPLESHEET-DRIVEN, MULTI-COHORT. Previously one invocation handled one
//    cohort via a per-cohort YAML. Now a samplesheet of cohort roots drives a
//    single execution: cohorts run in parallel, share one resume history, and
//    produce one provenance record.
//
// 2. FILES, NOT DIRECTORIES. The old processes took `path maf_dir`. Nextflow
//    hashes a directory input by NAME, so -resume happily reused cached tasks
//    after the contents changed. Every process now receives explicit file
//    lists resolved from the manifest.
//
// 3. NO PLACEHOLDER FILES. The old `ifEmpty(file("no_oncoprint.png"))` idiom
//    referenced files that do not exist; Nextflow fails when staging them.
//    Optional inputs are handled with `remainder: true` joins and empty lists.
//
// 4. NO `when:` DIRECTIVES. Toggles are evaluated once, here, rather than
//    being duplicated as both `if (params.run_x)` and `when: params.run_x`.
//
// 5. INGEST IS IN THE PIPELINE. The manifest, the completeness matrix and the
//    cohort .seg are all built here rather than hand-staged into
//    ../results/all_cohorts_analysis_ready/, which had no provenance.
// ===========================================================================

nextflow.enable.dsl = 2

include { BUILD_MANIFEST; MERGE_METADATA; BUILD_COHORT_SEG;
          PLOT_COMPLETENESS }                       from './modules/local/manifest.nf'
include { ONCOPRINT; CNV_FREQUENCY; CNV_BURDEN; EXPRESSION_PCA;
          COMPARATIVE; RECURRENT_FUSIONS; GISTIC;
          COHORT_REPORT }                           from './modules/local/analysis.nf'
include { COLLECTION_ROLLUP;
          BATCH_EFFECT_ANALYSIS }                   from './modules/local/collection.nf'
include { VCF2MAF }                                 from './modules/local/vcf2maf.nf'
include { RUN_PCGR }                                from './modules/local/pcgr.nf'
include { FILTER_ENSEMBLE }                         from './modules/local/filter_ensemble.nf'

// --------------------------------------------------------------------------- //
// Parameter validation -- fail at second zero, not after an hour of queueing
// --------------------------------------------------------------------------- //
def validateParams() {
    def errors = []

    if (!params.input && !params.cohort_dir) {
        errors << "Provide --input <samplesheet.csv> or --cohort_dir <path> (+ --cohort_name)."
    }
    if (params.cohort_dir && !params.cohort_name) {
        errors << "--cohort_dir requires --cohort_name."
    }
    if (params.amp_threshold <= 0) errors << "--amp_threshold must be > 0 (got ${params.amp_threshold})."
    if (params.del_threshold >= 0) errors << "--del_threshold must be < 0 (got ${params.del_threshold})."
    // Reference-file existence checks are skipped under -stub-run for the same
    // reason: the point of a stub is to validate WIRING before the environment
    // is complete.
    if (params.run_comparative && !params.gtf && !workflow.stubRun) {
        errors << "--run_comparative needs --gtf. Gene coordinates cannot be inferred " +
                  "from the MAF: doing so measures the span of observed mutations, not the gene."
    }
    if (params.gtf && !file(params.gtf).exists() && !workflow.stubRun) {
        errors << "--gtf not found: ${params.gtf}"
    }
    if (params.run_gistic && !params.refgene && !workflow.stubRun) {
        errors << "--run_gistic needs --refgene (the GISTIC .mat reference)."
    }
    if (params.run_gistic && params.refgene && !file(params.refgene).exists() && !workflow.stubRun) {
        errors << "--refgene not found: ${params.refgene}"
    }
    if (params.run_gistic && params.seg_file && !params.marker_file) {
        errors << "--run_gistic with a supplied --seg_file also needs --marker_file. " +
                  "Leave --seg_file unset to have the pipeline build both from CNVkit VCFs."
    }
    if (params.extra_metadata && !file(params.extra_metadata).exists()) {
        errors << "--extra_metadata not found: ${params.extra_metadata}"
    }

    // vcf2maf is not optional: the collection MAFs are unannotated, so without
    // this every mutation analysis silently has nothing to work with.
    // A stub run executes no task scripts, so requiring real tool paths would
    // block the one check that is meant to be runnable before anything is set
    // up. Structural validation below still applies.
    def needs_mafs = (params.run_oncoprint || params.run_comparative) && !workflow.stubRun
    if (needs_mafs) {
        // vcf2maf.pl is always needed. The VEP cache and reference FASTA are
        // needed only when re-annotating from scratch -- with inhibit_vep we
        // parse the CSQ field PCGR already wrote into the VCF.
        def required = params.inhibit_vep ? ['vcf2maf_path']
                                          : ['vcf2maf_path', 'vep_cache', 'ref_fasta']
        required.each { k ->
            if (!params[k]) {
                errors << "--${k} is required for mutation analyses. The MAFs in the " +
                          "MoHQ collection are NOT gene-annotated (Hugo_Symbol='Unknown'); " +
                          "they are regenerated from the PCGR VCFs." +
                          (params.inhibit_vep ? '' :
                           " Set --inhibit_vep true to reuse PCGR's own VEP annotation " +
                           "and avoid needing a VEP cache at all.")
            } else if (!file(params[k]).exists()) {
                errors << "--${k} not found: ${params[k]}"
            }
        }
    }
    if (params.run_comparative && !params.extra_metadata &&
        !(params.comparative_group_col in ['institution', 'cohort_id', 'multi_tumour'])) {
        errors << "--comparative_group_col '${params.comparative_group_col}' is not derivable " +
                  "from the delivery tree. Supply --extra_metadata with that column, or " +
                  "group by institution/cohort_id."
    }

    // The vcf2maf_host profile drops the VEP container and runs vcf2maf with
    // host modules. That is only coherent with --inhibit-vep, so the profile
    // sets params.inhibit_vep = true.
    //
    // But a -params-file OVERRIDES params set in a config profile. So if the
    // params file still says `inhibit_vep: false`, the profile's process
    // settings apply (no container, modules loaded) while the param does not --
    // and vcf2maf tries to run VEP on the bare host, failing with something
    // completely unrelated-looking:
    //     ERROR: Cannot find VEP script under: /home/you/miniconda3/bin
    // Catch the contradiction here instead.
    if (workflow.profile.tokenize(',').contains('vcf2maf_host') && !params.inhibit_vep) {
        errors << "-profile vcf2maf_host requires inhibit_vep: true, but it is false. " +
                  "A -params-file overrides profile params, so set `inhibit_vep: true` " +
                  "IN YOUR PARAMS FILE (not the profile). Without VEP in a container " +
                  "there is no vep for vcf2maf to call."
    }

    if (errors) {
        // `error` (not `exit 1`): it raises properly through Nextflow and is
        // not deprecated the way a bare exit in a workflow script is.
        error "Parameter validation failed:\n  - " + errors.join("\n  - ")
    }
}

// --------------------------------------------------------------------------- //
// Helper: collect one manifest column into a per-cohort file list.
//
// Every cohort appears in the output even when it has zero files for that
// column (CQ, for instance, has no RNA at all), so downstream joins never drop
// a cohort silently.
// --------------------------------------------------------------------------- //
def perCohortFiles(ch_rows, ch_keys, String column) {
    def ch_files = ch_rows
        .filter { cid, row -> row[column] && row[column] != 'NA' && row[column].trim() }
        .map    { cid, row -> tuple(cid, file(row[column], checkIfExists: true)) }
        .groupTuple()

    return ch_keys
        .join(ch_files, remainder: true)
        .map { cid, _k, files -> tuple(cid, files ?: []) }
}


workflow {

    validateParams()

    // ---------------------------------------------------------------------- //
    // 1. Cohort inputs
    // ---------------------------------------------------------------------- //
    ch_cohorts = params.input
        ? Channel.fromPath(params.input, checkIfExists: true)
              .splitCsv(header: true)
              .map { row ->
                  if (!row.cohort_id || !row.cohort_dir) {
                      error "Samplesheet needs columns 'cohort_id' and 'cohort_dir'. Got: ${row.keySet()}"
                  }
                  tuple(row.cohort_id, file(row.cohort_dir, checkIfExists: true))
              }
        : Channel.of(tuple(params.cohort_name, file(params.cohort_dir, checkIfExists: true)))

    ch_rlib = Channel.value(file("${projectDir}/bin/mohq_common.R", checkIfExists: true))

    // ---------------------------------------------------------------------- //
    // 2. Ingest: manifest + completeness
    // ---------------------------------------------------------------------- //
    BUILD_MANIFEST(ch_cohorts)

    // Optionally enrich with clinical metadata (sex, tumour type, outcome...).
    if (params.extra_metadata) {
        ch_extra = Channel.value(file(params.extra_metadata, checkIfExists: true))
        MERGE_METADATA(BUILD_MANIFEST.out.manifest.combine(ch_extra))
        ch_manifest = MERGE_METADATA.out.manifest
    } else {
        ch_manifest = BUILD_MANIFEST.out.manifest
    }

    // The MANIFEST is required -- it is how every later step finds its files.
    // The completeness FIGURE is just a report built from it, and can be off.
    if (params.run_completeness_plot) {
        PLOT_COMPLETENESS(
            BUILD_MANIFEST.out.completeness.join(ch_manifest),
            ch_rlib
        )
    }

    // ---------------------------------------------------------------------- //
    // 3. Resolve per-cohort file lists from the manifest
    //
    //    splitCsv(elem: 1) parses the manifest file while preserving the
    //    cohort_id in element 0.
    // ---------------------------------------------------------------------- //
    ch_rows = ch_manifest
        .splitCsv(header: true, sep: '\t', elem: 1)
        .map { cid, row -> tuple(cid, row) }

    // Only patients with the required assets feed the analyses. Excluding them
    // here -- once, visibly -- beats each R script silently dropping rows.
    ch_ready = ch_rows.filter { cid, row -> row.analysis_ready == 'yes' }

    // ---------------------------------------------------------------------- //
    // Genome-build guard.
    //
    // GRCh37 and GRCh38 coordinates are not interchangeable. A cohort .seg
    // built across both, or a gene overlap against a single GTF, would be
    // wrong -- and nothing downstream would raise an error. Fail here instead.
    // ---------------------------------------------------------------------- //
    ch_ready
        .map    { cid, row -> tuple(cid, row.genome_build ?: 'unknown') }
        .unique()
        .groupTuple()
        .map { cid, blds ->
            def real = blds.findAll { it != 'unknown' }
            if (real.unique().size() > 1) {
                error """
                Cohort ${cid} mixes genome builds: ${real.unique().join(', ')}.
                GRCh37 and GRCh38 coordinates cannot be combined -- the cohort
                .seg, the CNV frequency plot and every gene overlap would be
                silently wrong.
                Split the cohort by the manifest's `genome_build` column and run
                each separately, or lift over to one build first.
                """.stripIndent()
            }
            if (real.isEmpty()) {
                log.warn "Cohort ${cid}: genome build could not be determined from " +
                         "filenames. Assuming --genome_build ${params.genome_build}."
            } else if (real[0] != params.genome_build) {
                log.warn "Cohort ${cid}: files look like ${real[0]} but " +
                         "--genome_build is ${params.genome_build}. vcf2maf will be " +
                         "told ${params.genome_build}; fix this before trusting output."
            }
            return tuple(cid, real.isEmpty() ? params.genome_build : real[0])
        }
        .set { ch_build }

    ch_keys = ch_manifest.map { cid, _m -> tuple(cid, 'k') }

    // ---------------------------------------------------------------------- //
    // 3b. Regenerate annotated MAFs from the PCGR VCFs.
    //
    // The MAFs shipped in the collection are NOT gene-annotated -- Hugo_Symbol
    // is "Unknown" throughout -- so they cannot drive any mutation analysis.
    // vcf2maf + VEP runs here, per patient, cached by -resume and published to
    // a persistent directory so this expensive step happens once per patient
    // for the life of the project rather than once per pipeline run.
    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //
    // Mutation VCF, in priority order:
    //   1. the delivered PCGR VCF
    //   2. a PCGR VCF we regenerated earlier (regenerated_pcgr_dir)
    //   3. run PCGR now on the ensemble VCF   [run_pcgr_gapfill]
    //   4. fall back to the raw ensemble VCF  [only if gapfill is off]
    //
    // 4 is a LAST RESORT and is recorded in the manifest as
    // `mutation_vcf_source = ensemble_fallback`: the ensemble VCF is the caller
    // union while the PCGR VCF is filtered and tiered, so variant counts differ
    // systematically. Since which patients lost their PCGR VCF is an accident
    // of data management, that difference would masquerade as biology.
    // ---------------------------------------------------------------------- //
    ch_have_pcgr = ch_ready
        .filter { cid, row -> row.pcgr_vcf && row.pcgr_vcf != 'NA' }
        .map    { cid, row -> tuple(cid, row.patient_id,
                                    file(row.pcgr_vcf, checkIfExists: true)) }

    ch_need_pcgr = ch_ready
        .filter { cid, row -> (!row.pcgr_vcf || row.pcgr_vcf == 'NA') &&
                              row.ensemble_somatic_vcf && row.ensemble_somatic_vcf != 'NA' }
        .map    { cid, row ->
            // Reuse an earlier regeneration if one exists.
            def prev = file("${params.regenerated_pcgr_dir}/${cid}/" +
                            "${row.patient_id}_D.pcgr_acmg.${params.genome_build.toLowerCase()}.vcf.gz")
            // Each patient is rebuilt with the PCGR version ITS OWN run used,
            // recorded in its GenPipes ini. A cohort can span versions
            // (HM-19: 1.0.3 and 1.4.1), and forcing one on everyone would
            // create a difference that did not previously exist.
            // Version resolution, most specific first:
            //   1. params.pcgr_module      -- force one version on everyone
            //   2. the patient's OWN version, from its GenPipes ini
            //   3. params.pcgr_module_fallback -- only when the ini is missing
            //
            // 3 exists because some patients have no usable ini and record
            // pcgr_version = "unknown". Previously that was a hard error, which
            // is defensible but leaves the whole gap-fill blocked by a handful
            // of patients. The fallback is OPT-IN and recorded, so the choice
            // is visible rather than assumed: those patients get a version they
            // were not originally processed with, which is a real (small) batch
            // variable and belongs in mutation_provenance.tsv.
            def own = (row.pcgr_version && row.pcgr_version != 'unknown'
                       && row.pcgr_version != 'NA') ? row.pcgr_version : null
            def mod = params.pcgr_module ?: (own ?: params.pcgr_module_fallback)
            def src = params.pcgr_module ? 'forced'
                    : (own ? 'own_ini' : (params.pcgr_module_fallback ? 'fallback' : 'none'))
            // dna_tumour_sample travels with the patient: FILTER_ENSEMBLE needs
            // to tell PCGR's depth/VAF tags apart from the normal's, and
            // resolving that from a name is safer than from a column position.
            tuple(cid, row.patient_id,
                  file(row.ensemble_somatic_vcf, checkIfExists: true),
                  prev.exists() ? prev : null, mod, src,
                  row.dna_tumour_sample ?: 'NA')
        }

    ch_reused = ch_need_pcgr.filter { c, p, e, prev, m, s, t -> prev != null }
                            .map    { c, p, e, prev, m, s, t -> tuple(c, p, prev) }
    ch_torun  = ch_need_pcgr.filter { c, p, e, prev, m, s, t -> prev == null }
                            .map    { c, p, e, prev, m, s, t ->
                                if (!m) {
                                    error "Patient ${p} needs a regenerated PCGR VCF but its " +
                                          "PCGR version is unknown (no usable GenPipes ini).\n" +
                                          "  Choose one of:\n" +
                                          "    --pcgr_module_fallback mugqic/pcgr/1.0.3   " +
                                          "(only the unknown patients; others keep their own version)\n" +
                                          "    --pcgr_module mugqic/pcgr/1.0.3            " +
                                          "(force one version on EVERY patient)\n" +
                                          "  Or harvest the missing */parameters/*.ini and rebuild the manifest."
                                }
                                tuple(c, p, e, m, t)
                            }

    // Say plainly how many patients are being processed with a version they
    // were not originally run with -- easy to miss, and it is a batch variable.
    ch_need_pcgr.filter { c, p, e, prev, m, s, t -> prev == null && s == 'fallback' }
                .count()
                .subscribe { nfb ->
                    if (nfb > 0) log.warn """
                    ${nfb} patient(s) have no recorded PCGR version and will be
                    regenerated with the FALLBACK ${params.pcgr_module_fallback}.
                    That is not necessarily the version their original run used.
                    Recorded as pcgr_regenerated_fallback in mutation_provenance.tsv;
                    treat it as a covariate, or harvest their parameters/*.ini and rerun.
                    """.stripIndent()
                }

    if (params.run_pcgr_gapfill) {
        // PCGR is fed the FILTERED ensemble VCF, not the raw one -- GenPipes'
        // report_pcgr step consumes .2caller.flt.vcf.gz. Feeding it the raw
        // union produced ~365k-variant VCFs against ~24k in the delivered ones.
        FILTER_ENSEMBLE(ch_torun.map { c, p, e, m, t -> tuple(c, p, e, t) })

        // Re-attach each patient's PCGR module, which FILTER_ENSEMBLE does not
        // carry. Joined on (cohort, patient) rather than zipped: channel order
        // is not guaranteed, and a silent mis-pairing here would rebuild
        // patients with another patient's PCGR version.
        ch_pcgr_in = FILTER_ENSEMBLE.out.vcf
            .map  { c, p, v -> tuple([c, p], v) }
            .join( ch_torun.map { c, p, e, m, t -> tuple([c, p], m) } )
            .map  { key, v, m -> tuple(key[0], key[1], v, m) }

        RUN_PCGR(ch_pcgr_in)
        ch_vcfs = ch_have_pcgr.mix(ch_reused).mix(RUN_PCGR.out.vcf)

        // Provenance: the manifest is written before routing, so it cannot know
        // which source each patient ended up using. Record it here.
        ch_have_pcgr.map { c, p, v -> "${c}\t${p}\tpcgr_delivered" }
            .mix(ch_reused.map      { c, p, v -> "${c}\t${p}\tpcgr_regenerated_earlier" })
            .mix(RUN_PCGR.out.vcf.map { c, p, v -> "${c}\t${p}\tpcgr_regenerated_now" })
            .collectFile(name: 'mutation_provenance.tsv',
                         storeDir: "${params.outdir}/_provenance",
                         seed: "cohort_id\tpatient_id\tmutation_vcf_source\n",
                         sort: true, newLine: true)
    } else {
        // Gapfill disabled: use the ensemble VCF directly, having warned.
        ch_fallback = ch_ready
            .filter { cid, row -> (!row.pcgr_vcf || row.pcgr_vcf == 'NA') &&
                                  row.mutation_vcf && row.mutation_vcf != 'NA' }
            .map    { cid, row -> tuple(cid, row.patient_id,
                                        file(row.mutation_vcf, checkIfExists: true)) }
        ch_vcfs = ch_have_pcgr.mix(ch_reused).mix(ch_fallback)

        ch_have_pcgr.map { c, p, v -> "${c}\t${p}\tpcgr_delivered" }
            .mix(ch_reused.map   { c, p, v -> "${c}\t${p}\tpcgr_regenerated_earlier" })
            .mix(ch_fallback.map { c, p, v -> "${c}\t${p}\tensemble_fallback" })
            .collectFile(name: 'mutation_provenance.tsv',
                         storeDir: "${params.outdir}/_provenance",
                         seed: "cohort_id\tpatient_id\tmutation_vcf_source\n",
                         sort: true, newLine: true)

        ch_fallback.count().subscribe { nfb ->
            if (nfb > 0) {
                log.warn """
                ${nfb} patient(s) have no PCGR VCF and will use their ENSEMBLE
                somatic VCF instead. That is the caller union, not PCGR-filtered
                calls, so their variant counts are not comparable with the rest.
                Either adjust for `mutation_vcf_source`, or set
                --run_pcgr_gapfill true to regenerate their PCGR VCFs properly.
                """.stripIndent()
            }
        }

        // --inhibit-vep parses a CSQ field that only a VEP-annotated VCF has.
        // PCGR writes one; the raw ensemble VCF does not. Feeding those in
        // anyway means each fails the CSQ guard inside VCF2MAF, and since one
        // failed task terminates the run, the patients that WOULD have worked
        // never execute.
        //
        // So exclude them here rather than letting them fail. This is a routing
        // decision, not error suppression: there is no annotation to parse, and
        // the alternative (run VEP on them) is exactly what run_pcgr_gapfill is
        // for. They stay recorded in mutation_provenance.tsv above.
        if (params.inhibit_vep) {
            ch_vcfs = ch_have_pcgr.mix(ch_reused)
            ch_fallback.count().subscribe { nfb ->
                if (nfb > 0) {
                    log.warn """
                    inhibit_vep is ON, so the ${nfb} ensemble-VCF patient(s) are
                    EXCLUDED from VCF2MAF: their VCFs carry no CSQ field, so there
                    is nothing for --inhibit-vep to parse. Downstream analyses run
                    on the PCGR-annotated patients only.
                    To include them, set --run_pcgr_gapfill true (regenerates real
                    PCGR VCFs), or --inhibit_vep false with a working VEP setup.
                    """.stripIndent()
                }
            }
        }
    }

    // With inhibit_vep the PCGR VCF's existing CSQ annotation is parsed, so no
    // VEP cache or reference FASTA is needed. Stage placeholders in that case.
    //
    // The two placeholders must have DIFFERENT NAMES. Both pointed at
    // assets/NO_COHORT_MAP, so Nextflow staged two inputs into the same task
    // with the same filename and refused:
    //     input file name collision -- multiple input files named: NO_COHORT_MAP
    // Same class of bug as the /dev/null stub earlier: a placeholder is still a
    // staged file, and staging is by basename.
    ch_vep_cache = params.inhibit_vep
        ? Channel.value('NO_VEP_CACHE')
        : Channel.value(file(params.vep_cache, checkIfExists: true).toAbsolutePath().toString())

    // ref_fasta is required in BOTH modes -- see modules/local/vcf2maf.nf.
    if (!params.ref_fasta) {
        error "--ref_fasta is required. vcf2maf uses samtools faidx on it for " +
              "reference alleles even with --inhibit-vep; without it vcf2maf " +
              "silently falls back to a GRCh37 default path."
    }
    def _fa = file(params.ref_fasta, checkIfExists: true)
    // samtools faidx needs the indexes NEXT TO the fasta. Checked here because
    // the failure otherwise appears 20 tasks in, as a samtools error.
    def _fai = file("${_fa}.fai")
    def _gzi = file("${_fa}.gzi")
    if (!_fai.exists() || (_fa.name.endsWith('.gz') && !_gzi.exists())) {
        log.warn "ref_fasta index missing next to ${_fa} " +
                 "(.fai${_fa.name.endsWith('.gz') ? ' / .gzi' : ''}). " +
                 "Create with: samtools faidx ${_fa}"
    }
    ch_ref_fasta = Channel.value(_fa.toAbsolutePath().toString())

    VCF2MAF(ch_vcfs, ch_vep_cache, ch_ref_fasta)

    ch_mafs = VCF2MAF.out.maf
        .map { cid, pid, maf -> tuple(cid, maf) }
        .groupTuple()

    ch_cna       = perCohortFiles(ch_ready, ch_keys, 'cna_segments')
    ch_cnvkit    = perCohortFiles(ch_ready, ch_keys, 'cnvkit_vcf')
    ch_expr      = perCohortFiles(ch_rows,  ch_keys, 'expression_genes')
    ch_annofuse  = perCohortFiles(ch_rows,  ch_keys, 'anno_fuse')
    ch_linxfus   = perCohortFiles(ch_rows,  ch_keys, 'linx_fusion')

    // ---------------------------------------------------------------------- //
    // 4. Cohort .seg (built here, not hand-staged)
    // ---------------------------------------------------------------------- //
    if (params.seg_file) {
        // Escape hatch: use a pre-existing .seg rather than rebuilding.
        ch_seg = ch_manifest.map { cid, _m ->
            tuple(cid, file(params.seg_file, checkIfExists: true)) }
        ch_markers = params.marker_file
            ? ch_manifest.map { cid, _m -> tuple(cid, file(params.marker_file, checkIfExists: true)) }
            : Channel.empty()
    } else {
        BUILD_COHORT_SEG(
            ch_manifest.join(ch_cnvkit).filter { cid, m, vcfs -> vcfs.size() > 0 },
            ch_rlib
        )
        ch_seg     = BUILD_COHORT_SEG.out.seg
        ch_markers = BUILD_COHORT_SEG.out.markers
    }

    // ---------------------------------------------------------------------- //
    // 5. Analyses
    // ---------------------------------------------------------------------- //
    ch_figures = Channel.empty()

    // PLOT_COMPLETENESS runs earlier, in the ingest section, and its output was
    // never mixed in here -- so the report's Data completeness section always
    // reported "No completeness matrix was produced" even on runs where the
    // process had succeeded. The figure existed in results/; it just never
    // reached the report.
    if (params.run_completeness_plot) {
        ch_figures = ch_figures.mix(PLOT_COMPLETENESS.out.plots)
                               .mix(PLOT_COMPLETENESS.out.tables)
    }

    if (params.run_oncoprint) {
        // run_oncoprint.R is staged as an input purely so that editing it
        // invalidates the -resume cache. Nextflow does not hash bin/ scripts.
        ch_onco_script = Channel.value(
            file("${projectDir}/bin/run_oncoprint.R", checkIfExists: true))
        ch_panel = params.gene_panel
            ? Channel.value(file(params.gene_panel, checkIfExists: true))
            : Channel.value(file("${projectDir}/assets/NO_GENE_PANEL"))

        ONCOPRINT(
            ch_manifest.join(ch_mafs).join(ch_cna)
                       .filter { cid, m, mafs, cna -> mafs.size() > 0 },
            ch_rlib,
            ch_onco_script,
            ch_panel
        )
        ch_figures = ch_figures.mix(ONCOPRINT.out.all_plots)
    }

    if (params.run_cnv_frequency) {
        CNV_FREQUENCY(ch_manifest.join(ch_seg), ch_rlib)
        ch_figures = ch_figures.mix(CNV_FREQUENCY.out.plot)
    }

    if (params.run_cnv_burden) {
        CNV_BURDEN(ch_manifest.join(ch_seg), ch_rlib)
        ch_figures = ch_figures.mix(CNV_BURDEN.out.all_plots)
    }

    if (params.run_pca) {
        // A cohort with fewer than three RNA samples cannot yield a PCA.
        // Filtering here means CQ is skipped cleanly instead of erroring.
        EXPRESSION_PCA(
            ch_manifest.join(ch_expr).filter { cid, m, expr ->
                if (expr.size() < 3) {
                    log.warn "Cohort ${cid}: only ${expr.size()} RNA sample(s); skipping PCA."
                    return false
                }
                return true
            },
            ch_rlib
        )
        ch_figures = ch_figures.mix(EXPRESSION_PCA.out.plots)
    }

    if (params.run_comparative) {
        ch_gtf = Channel.value(file(params.gtf, checkIfExists: true))
        COMPARATIVE(
            ch_manifest.join(ch_mafs).join(ch_seg)
                       .filter { cid, m, mafs, seg -> mafs.size() > 0 },
            ch_rlib, ch_gtf
        )
        ch_figures = ch_figures.mix(COMPARATIVE.out.plots)
    }

    if (params.run_fusions) {
        RECURRENT_FUSIONS(
            ch_manifest.join(ch_annofuse).join(ch_linxfus)
                       .filter { cid, m, af, lx -> af.size() > 0 || lx.size() > 0 },
            ch_rlib
        )
        ch_figures = ch_figures.mix(RECURRENT_FUSIONS.out.plot)
    }

    if (params.run_gistic) {
        ch_refgene = Channel.value(file(params.refgene, checkIfExists: true))
        GISTIC(ch_seg.join(ch_markers), ch_refgene)
    }

    // ---------------------------------------------------------------------- //
    // 6. Collection tier -- all cohorts together
    //
    // Consumes summary tables only, so cost is flat in collection size. The
    // batch-effect analysis exploits the fact that the SAME cancer type appears
    // at MULTIPLE institutions: holding cancer type constant, any residual
    // institution effect is technical. That separation is only possible because
    // the design is crossed rather than confounded.
    // ---------------------------------------------------------------------- //
    if (params.run_collection) {
        ch_cohort_map = params.cohort_map
            ? Channel.value(file(params.cohort_map, checkIfExists: true))
            : Channel.value(file("${projectDir}/assets/NO_COHORT_MAP"))

        COLLECTION_ROLLUP(
            ch_manifest.map { cid, m -> m }.collect(),
            BUILD_MANIFEST.out.completeness.map { cid, c -> c }.collect(),
            ch_rlib,
            ch_cohort_map
        )

        // Per-patient metric tables feed the variance partition. Only the ones
        // that were actually produced are mixed in.
        ch_metrics = Channel.empty()
        if (params.run_cnv_burden) {
            ch_metrics = ch_metrics.mix(
                CNV_BURDEN.out.tables.map { cid, t -> t }.flatten()
                    .filter { it.name.endsWith('_cnv_burden.tsv') })
        }
        if (params.run_oncoprint) {
            ch_metrics = ch_metrics.mix(
                ONCOPRINT.out.tables.map { cid, t -> t }.flatten()
                    .filter { it.name.endsWith('_tmb.tsv') })
        }

        if (params.cohort_map) {
            BATCH_EFFECT_ANALYSIS(
                ch_manifest.map { cid, m -> m }.collect(),
                ch_metrics.collect().ifEmpty([]),
                ch_rlib,
                ch_cohort_map
            )
        } else {
            log.warn "--cohort_map not supplied: skipping the batch-effect analysis. " +
                     "Without a cohort_id -> cancer_type map the crossed design cannot " +
                     "be identified, and institution effects cannot be told apart from biology."
        }
    }

    // ---------------------------------------------------------------------- //
    // 7. Per-cohort report -- collect whatever figures were actually produced
    // ---------------------------------------------------------------------- //
    if (params.run_report) {
        ch_rmd = Channel.value(file("${projectDir}/assets/report_template.Rmd",
                                    checkIfExists: true))

        // Version files from every process that emitted one, so the Provenance
        // section stops reading "No version files staged".
        ch_versions = Channel.empty()
        if (params.run_oncoprint)     ch_versions = ch_versions.mix(ONCOPRINT.out.versions)
        if (params.run_cnv_frequency) ch_versions = ch_versions.mix(CNV_FREQUENCY.out.versions)
        if (params.run_cnv_burden)    ch_versions = ch_versions.mix(CNV_BURDEN.out.versions)
        if (params.run_pca)           ch_versions = ch_versions.mix(EXPRESSION_PCA.out.versions)
        if (params.run_fusions)       ch_versions = ch_versions.mix(RECURRENT_FUSIONS.out.versions)
        if (params.run_completeness_plot)
                                      ch_versions = ch_versions.mix(PLOT_COMPLETENESS.out.versions)

        ch_all_figs = ch_figures
            .map { cid, files -> tuple(cid, files instanceof List ? files : [files]) }
            .groupTuple()
            .map { cid, nested -> tuple(cid, nested.flatten()) }
            // Attach the version files and the parameter record to every cohort.
            .combine(ch_versions.collect().ifEmpty([]).toList())
            .map { cid, figs, vers -> tuple(cid, (figs + vers.flatten()).unique()) }

        COHORT_REPORT(
            ch_manifest
                .join(BUILD_MANIFEST.out.completeness)
                .join(ch_all_figs),
            ch_rmd
        )
    }
}


workflow.onComplete {
    log.info """
    ---------------------------------------------------------------
    Pipeline   : ${workflow.manifest.name} ${workflow.manifest.version}
    Completed  : ${workflow.complete}
    Duration   : ${workflow.duration}
    Success    : ${workflow.success}
    Results    : ${params.outdir}
    Work dir   : ${workflow.workDir}
    Command    : ${workflow.commandLine}
    ---------------------------------------------------------------
    """.stripIndent()
}

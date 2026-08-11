// ===========================================================================
// modules/local/pcgr.nf -- regenerate MISSING PCGR VCFs
//
// WHY
// ---
// Some patients never received a PCGR VCF: it was produced, then deleted, and
// is not recoverable. In MoHQ-HM-19 that is 25 of 47 patients. Meanwhile the
// delivered MAFs are unannotated collection-wide, so those 25 have no usable
// mutation data at all.
//
// The tempting shortcut is to run vcf2maf on their ensemble somatic VCF
// instead. Do not do that silently. The PCGR VCF is post-filtering; the
// ensemble VCF is the raw caller union. Variant COUNTS differ systematically
// between the two, which biases TMB and every mutation frequency -- and the
// split correlates with which patients happened to lose files, i.e. an
// accident of data management masquerading as biology.
//
// THE INPUT IS THE *FILTERED* ENSEMBLE VCF, NOT THE RAW ONE.
//
// This module previously fed PCGR the raw ensemble VCF, asserting in this very
// comment that doing so reproduced the delivered path. It did not. GenPipes'
// report_pcgr step consumes .ensemble.somatic.vt.annot.2caller.flt.vcf.gz,
// which has already been through TWO filters (>=2 callers, then depth/VAF).
// The scale of what was being skipped:
//
//     raw ensemble (delivered patient -1) : 423,764 variants
//     PCGR input after filtering          :  23,604 variants
//
// So the shortcut this comment warned against was, in effect, what the module
// itself was doing -- with the added irony of committing exactly the batch
// artefact it was written to prevent. Regenerated patients came out with
// ~365,000 variants against ~24,000 for delivered ones.
//
// FILTER_ENSEMBLE (modules/local/filter_ensemble.nf) now performs the missing
// step, so RUN_PCGR receives the same file GenPipes gave PCGR. Do not wire the
// raw ensemble VCF back into this process.
//
// Everything above was recovered from primary evidence -- the [filter_ensemble]
// section of the delivered config trace and the ##bcftools_viewCommand line in
// the delivered PCGR VCF headers -- after being wrong when inferred from the
// delivered files alone.
//
// This process runs ONLY for patients missing a PCGR VCF. Patients who have
// one are untouched -- we never re-derive data that already exists.
//
// PCGR is an environment module on Cardinal (mugqic/pcgr/1.0.3), not a
// container, so the `module` directive is used rather than `container`.
// ===========================================================================

process RUN_PCGR {
    tag   { "${patient_id} [${pcgr_module}]" }
    label 'process_pcgr'

    // PER-PATIENT PCGR VERSION.
    //
    // A cohort can span several PCGR versions -- MoHQ-HM-19 has 1.0.3 (37
    // patients) and 1.4.1 (12). Each patient's GenPipes ini records the version
    // ITS run used, so we regenerate with that one. Forcing a single version on
    // everyone would introduce a difference that was not previously there:
    // a 1.4.1 patient rebuilt with 1.0.3 stops matching its own cohort-mates.
    //
    // params.pcgr_module, if set, overrides this for every patient. Leave it
    // null to use each patient's recorded version.
    // htslib joins pcgr on the module list so `tabix` is available for the
    // index fallback below. The pcgr conda env may ship it, but relying on that
    // is the kind of assumption that produced exit 127 earlier.
    module { params.htslib_module ? "${pcgr_module}:${params.htslib_module}" : pcgr_module }

    // MUST run on the HOST, not in a container.
    //
    // nextflow.config sets a global default container (pipeline.sif). The
    // `module` directive loads on the host, but the script then executes INSIDE
    // that container -- where pcgr does not exist. Result: exit 127, "command
    // not found", on every patient, and because 127 looks like any other
    // failure the retry logic cheerfully repeats it.
    //
    // Same mistake as VCF2MAF, fixed there with container = null in the
    // vcf2maf_host profile. A `module` directive and a container are almost
    // always contradictory: if a process uses modules, it runs on the host.
    container = null

    // scratch = false, deliberately, overriding the cardinal profile.
    //
    // With scratch = true the task runs in a node-local temp directory and
    // Nextflow's wrapper ends with:
    //     nxf_unstage() { if [[ ${nxf_main_ret:=0} == 0 ]]; then ...copy out... }
    //     on_exit()     { rm -rf $NXF_SCRATCH }
    // Outputs are copied back ONLY on success, and the scratch directory is
    // removed unconditionally. So when PCGR fails, its log and its partial
    // out/ directory are deleted before anyone can read them -- the work dir
    // contains nothing but .command.* files, .command.err is empty, and there
    // is no error message anywhere on disk.
    //
    // That cost an evening: three separate attempts to find a failure message
    // that was being destroyed microseconds after it was written.
    //
    // PCGR's I/O is modest and this is the step most likely to fail, so running
    // it directly in the work directory is the right trade.
    scratch = false

    // Persist outside the work directory: this is expensive to regenerate and
    // should survive a work/ cleanup, exactly like the converted MAFs.
    // enabled: !workflow.stubRun -- same trap as VCF2MAF. The stub `touch`es a
    // 0-byte VCF; published with overwrite:false it becomes the permanent
    // "regenerated" file for that patient, is found by existingRegeneratedPcgr(),
    // and skips PCGR forever while contributing nothing.
    publishDir "${params.regenerated_pcgr_dir}/${cohort_id}",
               mode: 'copy', overwrite: false, enabled: !workflow.stubRun

    input:
    tuple val(cohort_id), val(patient_id), path(ensemble_vcf), val(pcgr_module)

    output:
    tuple val(cohort_id), val(patient_id),
          path("${patient_id}_D.pcgr_acmg.${params.genome_build.toLowerCase()}.vcf.gz"),
          emit: vcf, optional: true
    path "${patient_id}.pcgr.log", emit: log

    script:
    def assembly = params.genome_build.toLowerCase()   // pcgr wants grch38
    def site     = params.pcgr_tumor_site
    def extra    = params.pcgr_extra_args ?: ''
    """
    set -euo pipefail

    # DO NOT LET THE LAUNCHING SHELL'S PYTHON REACH IN HERE.
    #
    # PCGR ships its own python. If the shell that started Nextflow had a
    # python environment module loaded (mugqic/python/3.10.4, say, for
    # FILTER_ENSEMBLE's format2pcgr.py), PYTHONHOME/PYTHONPATH propagate
    # through Nextflow into every SLURM job, and PCGR's interpreter tries to
    # load a foreign standard library:
    #
    #     Fatal Python error: init_sys_streams: can't initialize sys standard streams
    #     ImportError: cannot import name 'open_code' from 'io' (unknown location)
    #     -> SIGABRT, exit 134, on every patient
    #
    # The failure has nothing to do with PCGR and everything to do with who
    # typed `module load` before `nextflow run`. That is far too fragile to
    # leave as a documented caveat, so scrub it here.
    #
    # Only variables pointing INTO a foreign python tree are removed --
    # blanket-unsetting PYTHONPATH would break the pcgr module if a future
    # version uses it for its own libraries.
    case "\${PYTHONHOME:-}" in
        */software/python/*) echo "[pcgr] dropping inherited PYTHONHOME=\$PYTHONHOME"
                             unset PYTHONHOME ;;
    esac
    if [ -n "\${PYTHONPATH:-}" ]; then
        CLEAN=\$(printf '%s' "\$PYTHONPATH" | tr ':' '\\n' \\
                | grep -v '/software/python/' | paste -sd: -)
        if [ "\$CLEAN" != "\$PYTHONPATH" ]; then
            echo "[pcgr] pruning inherited PYTHONPATH"
            if [ -n "\$CLEAN" ]; then export PYTHONPATH="\$CLEAN"; else unset PYTHONPATH; fi
        fi
    fi

    mkdir -p out

    # PCGR requires the tabix index NEXT TO the input VCF:
    #     ERROR - Tabix file (i.e. '.gz.tbi') is not present for the bgzipped
    #             VCF input file (...). Please make sure your input VCF is
    #             properly compressed and indexed (bgzip + tabix)
    #
    # Nextflow stages the VCF as a symlink but NOT its index -- the .tbi sits
    # beside the original in the collection and never reaches the work
    # directory. The harvest does download it (*.vcf.gz.tbi, in the
    # "signatures" set), so rather than plumb a second input through the
    # channel, link it if present and build it if not. Indexing a somatic VCF
    # takes seconds.
    if [ ! -e "${ensemble_vcf}.tbi" ]; then
        src=\$(readlink -f ${ensemble_vcf})
        if [ -e "\$src.tbi" ]; then
            ln -sf "\$src.tbi" "${ensemble_vcf}.tbi"
            echo "[pcgr] linked existing index: \$src.tbi"
        else
            echo "[pcgr] no .tbi found; building one with tabix"
            tabix -p vcf ${ensemble_vcf}
        fi
    fi

    # The depth/AF tag names in pcgr_extra_args (TAL/TDP/TVAF/NDP/NVAF) are
    # copied verbatim from the delivered runs' [report_pcgr] section. They are
    # specific to the GenPipes ensemble VCF -- with different tag names PCGR
    # silently filters on fields that are not there.
    #
    # --vcf2maf is deliberately NOT passed. The delivered runs used it, and it
    # is what produced the unannotated MAFs. We want the VCF only; the MAF comes
    # from the VCF2MAF process with a working VEP cache.

    pcgr \\
        --input_vcf        ${ensemble_vcf} \\
        --pcgr_dir         ${params.pcgr_bundle} \\
        --output_dir       out \\
        --genome_assembly  ${assembly} \\
        --sample_id        ${patient_id}_D \\
        --tumor_site       ${site} \\
        --assay            ${params.pcgr_assay} \\
        --force_overwrite \\
        ${extra} \\
        2>&1 | tee ${patient_id}.pcgr.log

    # `| tee` NOT `> file`.
    #
    # The redirect sent everything to a file that then had to be found and read
    # by hand. When the first real gap-fill failed, the log was invisible in
    # .command.out and .command.err was empty, so there was no error message
    # anywhere obvious -- twenty minutes hunting for something that had been
    # captured all along and hidden.
    #
    # With pipefail set, \${PIPESTATUS[0]} is not needed: the pipeline's status
    # is pcgr's. But say so explicitly rather than relying on it.
    pcgr_status=\${PIPESTATUS[0]}
    if [ "\$pcgr_status" -ne 0 ]; then
        echo "PCGR exited \$pcgr_status for ${patient_id}. Last 40 lines above." >&2
        exit "\$pcgr_status"
    fi

    # PICK THE RIGHT OUTPUT FILE BY NAME, NOT BY find | head -1.
    #
    # PCGR 1.0.3 writes THREE VCFs, and `find ... | head -1` returns whichever
    # the filesystem lists first -- an arbitrary choice that silently picked the
    # wrong one:
    #
    #   <sample>.pcgr_acmg.grch38.mp_input.vcf.gz   286,430 variants, NO CSQ
    #       ^ the unannotated input written for mutational-signature analysis
    #   <sample>.pcgr_acmg.grch38.pass.vcf.gz       365,430 variants, has CSQ
    #   <sample>.pcgr_acmg.grch38.vcf.gz            365,430 variants, has CSQ
    #       ^ THIS one, and it already matches the delivered naming convention
    #
    # Picking mp_input produced a VCF with no CSQ field, which then failed two
    # steps later in VCF2MAF with "inhibit_vep is set but this VCF has no CSQ",
    # pointing at the wrong process entirely. It also made the regenerated
    # variant counts look inexplicable next to the delivered ones.
    expected="out/${patient_id}_D.pcgr_acmg.${assembly}.vcf.gz"
    if [ -s "\$expected" ]; then
        found="\$expected"
    else
        # Fall back only to files that could plausibly be the annotated output.
        found=\$(find out -name '*.vcf.gz' \\
                    ! -name '*germline*' ! -name '*mp_input*' ! -name '*.pass.*' \\
                    | sort | head -1)
    fi
    if [ -z "\$found" ]; then
        echo "PCGR produced no annotated VCF. Files in out/:" >&2
        ls -la out/ >&2
        tail -40 ${patient_id}.pcgr.log >&2
        exit 1
    fi
    echo "[pcgr] selected output: \$found"

    # The whole point of this file is its CSQ annotation -- VCF2MAF parses it
    # with --inhibit-vep. Check HERE, so a wrong pick fails in the process that
    # made the choice rather than surfacing later as someone else's error.
    # NOTE ON `set +o pipefail` BELOW -- it is not laziness.
    #
    # `zcat big.vcf.gz | grep -qm1 PATTERN` does the right thing in an
    # interactive shell and the WRONG thing here. grep -q exits the moment it
    # matches; zcat then writes into a closed pipe, dies of SIGPIPE with status
    # 141, and `pipefail` promotes that to the pipeline's status. A SUCCESSFUL
    # match therefore reports failure.
    #
    # This cost a full round of 24 PCGR runs: every output was reported as
    # having no CSQ field, including files verified by hand to have one. The
    # guard added to catch a silent error became a silent error itself.
    #
    # Small outputs hide this -- if the writer finishes before grep exits (the
    # pipe buffer absorbs it) there is no SIGPIPE. It only bites on large
    # files, i.e. exactly the real ones.
    csq_present() {
        local f="\$1" rc
        set +o pipefail
        zcat "\$f" 2>/dev/null | grep -qm1 '^##INFO=<ID=CSQ'
        rc=\$?
        set -o pipefail
        return \$rc
    }

    if ! csq_present "\$found"; then
        echo "PCGR output \$found has NO CSQ INFO field." >&2
        echo "That is the annotation VCF2MAF needs. Files produced:" >&2
        for v in out/*.vcf.gz; do
            if csq_present "\$v"; then s=yes; else s=NO; fi
            printf '  %-60s CSQ:%s\\n' "\$(basename \$v)" "\$s" >&2
        done
        exit 1
    fi

    # EXISTENCE IS NOT ENOUGH.
    #
    # The first gap-fill attempt published a 0-byte VCF for MoHQ-HM-19-2: PCGR
    # exited 0, wrote an empty file, `find` found it, and this step copied and
    # published it as though it were a result. An empty VCF downstream means a
    # patient silently contributes no variants to the oncoplot -- exactly the
    # class of failure this pipeline exists to prevent, reintroduced by me in
    # the one process that had never run.
    #
    # So check the file is non-empty AND actually contains variant records.
    if [ ! -s "\$found" ]; then
        echo "PCGR wrote an EMPTY VCF for ${patient_id}: \$found (0 bytes)." >&2
        echo "PCGR exited 0, so this is a silent failure. Last 40 log lines:" >&2
        tail -40 ${patient_id}.pcgr.log >&2
        exit 1
    fi
    n_var=\$(zcat "\$found" 2>/dev/null | grep -vc '^#' || true)
    if [ "\${n_var:-0}" -lt 1 ]; then
        echo "PCGR VCF for ${patient_id} has 0 variant records (header only)." >&2
        echo "Check the depth/AF tag names in pcgr_extra_args match this VCF:" >&2
        echo "  \$(zcat ${ensemble_vcf} 2>/dev/null | grep -m1 '^##FORMAT' || true)" >&2
        tail -40 ${patient_id}.pcgr.log >&2
        exit 1
    fi
    echo "[pcgr] ${patient_id}: \$n_var variant records"

    cp "\$found" ${patient_id}_D.pcgr_acmg.${assembly}.vcf.gz
    """

    stub:
    """
    touch ${patient_id}_D.pcgr_acmg.${params.genome_build.toLowerCase()}.vcf.gz
    touch ${patient_id}.pcgr.log
    """
}


// ---------------------------------------------------------------------------
// Reuse a previously regenerated VCF rather than re-running PCGR.
// PCGR is slow; this keeps the work when the Nextflow work dir is cleared.
// ---------------------------------------------------------------------------
def existingRegeneratedPcgr(cohort_id, patient_id) {
    def f = file("${params.regenerated_pcgr_dir}/${cohort_id}/" +
                 "${patient_id}_D.pcgr_acmg.${params.genome_build.toLowerCase()}.vcf.gz")
    return f.exists() ? f : null
}

# MoHQ cohort-scale analysis pipeline

A Nextflow pipeline that takes the delivered MoHQ files for one or many cohorts
and produces cohort-level results: an oncoplot, copy-number frequency and
burden, expression PCA, recurrent fusions, and a completeness audit.

It exists because the collection could not be analysed as delivered. Every MAF
in the collection has `Hugo_Symbol = "Unknown"` for every row, so any gene-level
analysis built on them returns an empty result **with no error**. The pipeline
regenerates annotated MAFs from the PCGR VCFs and records, per patient, exactly
which source and which tool versions were used.

---

## Contents

**New here?** Read [`docs/TUTORIAL.md`](docs/TUTORIAL.md) instead. It walks one
real cohort from harvest to report — every command, every number, and what to
check before moving to the next step. This README is the reference you come back
to once you know the shape of a run.

- [Quick start](#quick-start)
- [How the pipeline works](#how-the-pipeline-works)
- [What each file does](#what-each-file-does)
  - [Pipeline entry points](#pipeline-entry-points)
  - [Nextflow modules](#nextflow-modules-modulelocal)
  - [Analysis scripts called by the pipeline](#analysis-scripts-called-by-the-pipeline-bin)
  - [Utilities you run by hand](#utilities-you-run-by-hand-bin)
- [Parameters](#parameters)
- [Profiles](#profiles)
- [Outputs](#outputs)
- [Testing](#testing)
- [Known constraints](#known-constraints)

---

## Quick start

The setup that currently works on Cardinal:

```bash
source bin/load_cardinal_modules.sh
export NXF_WORK=/home/$USER/MoHQ/work

# 1. Validate the whole workflow structure. Seconds, touches no data.
nextflow run . -profile cardinal_local -stub-run \
    -params-file params/cardinal_test.yaml

# 2. Real run, one cohort, submitted to Slurm.
nextflow run . -profile cardinal,vcf2maf_host \
    -params-file params/cardinal_test.yaml -resume
```

Run step 1 after **any** edit to `main.nf` or a module. It catches wiring errors
in about four seconds, before they cost a real run.

Getting data in the first place:

```bash
export MOHQ_PROJECT=<openstack project id>

# See what would transfer -- moves nothing
python3 bin/harvest_juno.py --remote juno:$MOHQ_PROJECT:MOH-Q \
    --cohort HM/MoHQ-HM-19 --dest /home/$USER/MoHQ/collection --dry-run

# Then without --dry-run, adding --verify to confirm sizes after transfer
```

---

## How the pipeline works

```
  delivery tree (per cohort)
        |
        v
  BUILD_MANIFEST ......... walk the tree; emit manifest, completeness,
        |                   samples and ambiguities tables
        |
        +--> PLOT_COMPLETENESS ...... which patients are missing which files
        |
        +--> route each patient to a mutation source:
        |       PCGR VCF delivered  -> use it
        |       regenerated earlier -> reuse it
        |       neither             -> FILTER_ENSEMBLE -> RUN_PCGR
        |                                             [run_pcgr_gapfill]
        |                              or exclude    [inhibit_vep]
        |
        v
  VCF2MAF ................ regenerate a gene-annotated MAF per patient,
        |                   published to derived_maf_dir and reused later
        |
        +--> BUILD_COHORT_SEG ...... CNVkit VCFs -> one cohort .seg
        |
        v
  ONCOPRINT  CNV_FREQUENCY  CNV_BURDEN  EXPRESSION_PCA  RECURRENT_FUSIONS
  COMPARATIVE [opt]  GISTIC [opt]
        |
        v
  COHORT_REPORT .......... one self-contained HTML per cohort
        |
        v
  COLLECTION_ROLLUP + BATCH_EFFECT_ANALYSIS   [run_collection, many cohorts]
```

Two design decisions worth knowing:

**Analysis steps receive file lists, not directories.** Nextflow identifies a
task by hashing its inputs, and for a directory only the *name* is hashed.
Passing directories meant adding or changing a file inside one did not
invalidate the cached result, so `-resume` could silently return stale output.

**Converted MAFs are published outside the work directory**, to
`derived_maf_dir`. Regenerating them is the most expensive step, so they are
reused across runs and survive work-directory cleanup.

---

## What each file does

### Pipeline entry points

| File | Purpose |
|---|---|
| `main.nf` | The workflow. Reads the samplesheet or single-cohort parameters, validates them, routes each patient to a mutation source, and calls the modules. |
| `nextflow.config` | All parameters and their defaults, plus the execution profiles. |
| `conf/base.config` | Resource requests per process label, with `check_max()` clamping to `max_cpus` / `max_memory` / `max_time`. Clamps hard under `-stub-run`. |
| `pipeline.def` | Apptainer definition for the R analysis image. Pins the base image, R version and a dated package snapshot so a rebuild reproduces the same versions. |

### Nextflow modules (`modules/local/`)

| Module | Processes | What they do |
|---|---|---|
| `manifest.nf` | `BUILD_MANIFEST`, `MERGE_METADATA`, `BUILD_COHORT_SEG`, `PLOT_COMPLETENESS` | Ingest. Turn a delivery tree into tables, join optional clinical metadata, assemble the cohort `.seg`, plot completeness. |
| `vcf2maf.nf` | `VCF2MAF`, `COLLECT_DERIVED_MAFS` | Regenerate gene-annotated MAFs. `COLLECT_DERIVED_MAFS` reuses any MAF already present in `derived_maf_dir` instead of recomputing it. |
| `filter_ensemble.nf` | `FILTER_ENSEMBLE` | The step between the ensemble VCF and PCGR. Runs MOH's own `format2pcgr.py` to add the depth/VAF INFO tags and keep variants seen by ≥2 callers, then applies the depth/VAF filter. Without it PCGR receives ~15× too many variants. |
| `pcgr.nf` | `RUN_PCGR` | Gap-fill: regenerate a PCGR VCF for patients whose delivered one is missing, from the **filtered** VCF above, using **that patient's own** recorded PCGR version. |
| `analysis.nf` | `ONCOPRINT`, `CNV_FREQUENCY`, `CNV_BURDEN`, `EXPRESSION_PCA`, `COMPARATIVE`, `RECURRENT_FUSIONS`, `GISTIC`, `COHORT_REPORT` | The cohort-level analyses and the per-cohort HTML report. |
| `collection.nf` | `COLLECTION_ROLLUP`, `BATCH_EFFECT_ANALYSIS` | Cross-cohort tier. Needs several cohorts and `cohort_map`. |

### Analysis scripts called by the pipeline (`bin/`)

You do not run these directly; Nextflow does. Each takes `--lib bin/mohq_common.R`.

| Script | Called by | What it does |
|---|---|---|
| `mohq_common.R` | all R scripts | Shared library. Segment-to-bin expansion, FGA with a per-sample denominator, MAF loading with an annotation guard, PCGR CNA column resolution, patient-ID harmonisation. **Most of the correctness lives here.** |
| `build_manifest.py` | `BUILD_MANIFEST` | Walks a delivery tree and emits `manifest`, `completeness`, `samples` and `ambiguities` TSVs. Handles all three PCGR naming conventions, alphanumeric sample numbers and tissue-qualified sample types. Caches on directory mtimes. |
| `build_cohort_seg.R` | `BUILD_COHORT_SEG` | Reads per-patient CNVkit VCFs and writes one cohort `.seg` plus a markers file. |
| `plot_completeness.R` | `PLOT_COMPLETENESS` | Heatmap of which assets each patient has, plus a summary table. |
| `run_oncoprint.R` | `ONCOPRINT` | maftools oncoplot, TMB, gene summary. Joins the PCGR CNA table onto the MAF by gene symbol. |
| `run_cnv_frequency.R` | `CNV_FREQUENCY` | Genome-wide gain/loss frequency across the cohort, binned. |
| `run_fga_burden.R` | `CNV_BURDEN` | Fraction of genome altered per patient, with a per-sample profiled-length denominator. |
| `run_pca.R` | `EXPRESSION_PCA` | PCA of expression, interactive HTML, coloured by a real variable. |
| `run_fusions.R` | `RECURRENT_FUSIONS` | Recurrent fusions from Anno-FUSE and LINX outputs. |
| `run_comparative.R` | `COMPARATIVE` | Stratified comparison between two groups. Needs a GTF for gene coordinates. |
| `collection_rollup.R` | `COLLECTION_ROLLUP` | Completeness and QC drift across all cohorts. |
| `batch_effect_analysis.R` | `BATCH_EFFECT_ANALYSIS` | Holds cancer type constant and varies institution, so technical effect can be measured rather than assumed. Needs `cohort_map`. |

### Utilities you run by hand (`bin/`)

| Script | What it answers |
|---|---|
| `harvest_juno.py` | Copy only the files the analyses read from the object store. `--dry-run` first, `--verify` to confirm sizes. Roughly a 6,000-fold reduction versus the full collection. |
| `harvest_all.sh` | `harvest_juno.py` across every cohort. |
| `check_maf_annotation.py` | Are the delivered MAFs annotated? Judges by **distinct gene symbols**, not by the fraction of `Unknown` — see [Known constraints](#known-constraints). Works locally or by streaming from the object store. |
| `check_genome_build.sh` | Which genome build, read from BAM and VCF **headers** rather than filenames. |
| `audit_all_builds.sh` | The above across every cohort. |
| `check_tool_versions.sh` | Which PCGR, VEP and GenPipes versions produced each patient's files. Version is a batch variable. |
| `check_csq.sh` | Do the PCGR VCFs carry a `CSQ` field? Decides whether `inhibit_vep` is usable. |
| `check_cardinal.sh` | What the cluster provides: scheduler, partitions, modules, disk. |
| `check_layout.sh` | Are all expected pipeline files present and executable? Run after copying files between machines. |
| `load_cardinal_modules.sh` | `source` this to load Nextflow, Java and Apptainer. |
| `setup_vep.sh` | Check for or install a VEP cache and `vcf2maf`. Only needed when `inhibit_vep = false`. |
| `test_vcf2maf.sh` | Run vcf2maf on one patient outside Nextflow, for debugging. |

---

## Parameters

Set these in a YAML file and pass `-params-file`, or on the command line as
`--name value`. A `-params-file` **overrides** anything set in a profile.

### Input and output

| Parameter | Default | Meaning |
|---|---|---|
| `input` | `null` | Samplesheet CSV with columns `cohort_id,cohort_dir`. Use this for many cohorts. |
| `cohort_dir` | `null` | Single-cohort mode: path to one cohort's delivery tree. |
| `cohort_name` | `null` | Single-cohort mode: the cohort ID. |
| `outdir` | `./results` | Where results are published. |
| `extra_metadata` | `null` | Optional TSV/CSV of clinical metadata, joined on patient ID. |
| `collection_name` | `MoHQ` | Label used in collection-level outputs. |

### Mutation source and annotation

| Parameter | Default | Meaning |
|---|---|---|
| `inhibit_vep` | `true` | Parse the `CSQ` field PCGR already wrote instead of re-running VEP. Seconds per patient instead of tens of minutes, and gene symbols match the PCGR reports. Requires a delivered PCGR VCF. |
| `ref_fasta` | `null` | **Required in both modes.** vcf2maf calls `samtools faidx` on it regardless of whether VEP runs. Needs `.fai` (and `.gzi` if bgzipped) beside it. |
| `vep_cache` | `null` | VEP cache directory. Only when `inhibit_vep = false`. |
| `vep_cache_version` | `115` | Cache release. Match the release PCGR used, or gene symbols may disagree between the MAF and the CNA table. |
| `vep_container` | `null` | Image with `vep` on `PATH`. Only when `inhibit_vep = false`. |
| `vep_path` | `null` | Leave unset when running inside `vep_container`. |
| `samtools_exec` | `null` | Path to `samtools` **inside** whatever container VCF2MAF runs in. Unnecessary with `-profile vcf2maf_host`. |
| `tabix_exec` | `null` | Same, for `tabix`. vcf2maf hard-requires both. |
| `derived_maf_dir` | `<project>/derived_mafs` | Where regenerated MAFs are published and reused. Keep off any purged filesystem. |
| `min_annotated_genes` | `200` | Reject a regenerated MAF with fewer distinct gene symbols than this. |
| `max_unknown_frac` | `0.95` | Reject above this fraction of `Unknown`. Deliberately high — see [Known constraints](#known-constraints). |
| `genome_build` | `GRCh38` | Passed to vcf2maf as `--ncbi-build`. |

### PCGR gap-fill

| Parameter | Default | Meaning |
|---|---|---|
| `run_pcgr_gapfill` | `false` | Regenerate PCGR VCFs for patients whose delivered one is missing. Expensive; plan it. |
| `pcgr_module` | `null` | Force one PCGR version. Leave null to use each patient's own recorded version. |
| `pcgr_bundle` | CVMFS path | PCGR reference bundle. |
| `pcgr_tumor_site` | `0` | `0` = unspecified. |
| `pcgr_assay` | `WGS` | Sets the TMB denominator. `WGS` ≈ 2800 Mb, exome ≈ 30 Mb. |
| `pcgr_extra_args` | see config | Copied verbatim from the original GenPipes ini so regenerated files match delivered ones. |
| `regenerated_pcgr_dir` | `<project>/regenerated_pcgr` | Where gap-filled VCFs are kept and reused. |

### Which analyses run

| Parameter | Default | Needs |
|---|---|---|
| `run_completeness_plot` | `true` | — |
| `run_oncoprint` | `true` | annotated MAFs |
| `run_cnv_frequency` | `true` | cohort `.seg` |
| `run_cnv_burden` | `true` | cohort `.seg` |
| `run_pca` | `true` | expression files |
| `run_fusions` | `true` | Anno-FUSE / LINX outputs |
| `run_report` | `true` | — |
| `run_comparative` | `false` | `gtf` **and** a group column |
| `run_gistic` | `false` | `refgene` (.mat) |
| `run_collection` | `false` | several cohorts **and** `cohort_map` |

### Analysis settings

| Parameter | Default | Meaning |
|---|---|---|
| `amp_threshold` | `0.58` | log2 ratio above which a segment is a gain. |
| `del_threshold` | `-0.58` | log2 ratio below which a segment is a loss. |
| `cnv_bin_size` | `1000000` | Bin width in bp (1 Mb) for the CNV frequency plot. See below. |

**On `cnv_bin_size`.** Copy-number calls are *segments* — "this stretch of
chromosome 7, from 12.4 Mb to 61.9 Mb, is gained" — and every patient's segments
start and end in different places, so they cannot be counted directly. Binning
lays a fixed ruler over the genome, here one mark every 1 Mb, and asks of each
bin: how many patients are gained here? That gives one comparable number per
position, which is what the frequency plot draws.

1 Mb is a **convention**, not a value derived from this data — it is what TCGA
and most cohort CNV figures use, and it was chosen on that basis. The trade-off:
smaller bins show finer detail but produce noisier plots and more rows (the
binning step refuses to run past 50 M rows); larger bins smooth real focal
events away. 
    
If the cohort's biology is focal (e.g. *ERBB2* amplification in breast), 1 Mb may
be too coarse; drop to 100 kb and compare. That is a judgement call for the
analysis, not a fixed property of the pipeline, which is why it is a parameter.
| `include_sex_chr` | `false` | Include chrX/chrY. Leave off for sex-stratified comparisons, or karyotype dominates. |
| `min_seg_markers` | `0` | Drop segments with fewer supporting markers. |
| `min_profiled_mb` | `1000` | Flag samples whose profiled genome is smaller than this. |
| `oncoplot_top_genes` | `25` | Genes shown in the oncoplot. |
| `callable_mb` | `30` | **TMB denominator, in Mb. Must match what the numerator counts.** `run_oncoprint.R` counts non-synonymous (coding) variants, so the denominator is the coding footprint (~30–40 Mb) **even for WGS**. Use ~2800 only if you also count non-coding variants. |
| `pca_norm_method` | `log2` | `log2`, `log10`, `zscore` or `none`. |
| `pca_top_var_genes` | `2000` | Most variable genes used; `0` = all. |
| `pca_colour_by` | `institution` | Must be a real varying column, not a constant. |
| `comparative_group_col` | `institution` | Column to stratify by. |
| `comparative_n_genes` | `20` | Genes in the comparison. |
| `comparative_min_group_n` | `5` | Minimum patients per group. |
| `target_genes` | `null` | Comma-separated list; null = data-driven. |
| `min_fusion_recurrence` | `2` | Patients a fusion must appear in. |
| `batch_min_per_cell` | `5` | Minimum patients per (cancer type, institution) cell. |
| `gistic_*` | see config | GISTIC options. `broad` and `armpeel` default to on. |

### Execution

| Parameter | Default | Meaning |
|---|---|---|
| `container` | `null` | Path to `pipeline.sif`. |
| `slurm_account` | `null` | Slurm account, if not set by the profile. |
| `max_cpus` | `32` | Upper clamp on any task. |
| `max_memory` | `250.GB` | Upper clamp. |
| `max_time` | `48.h` | Upper clamp. |
| `manifest_cache_dir` | `<project>/.manifest_cache` | Manifest cache location. Safe to delete. |
| `use_manifest_cache` | `true` | Cache is keyed on directory mtimes, so a stale cache cannot produce a wrong manifest — only a slower run. |
| `completeness_max_tile_rows` | `150` | Above this, the completeness heatmap switches to a summarised form. |

---

## Profiles

Combine with commas: `-profile cardinal,vcf2maf_host`.

| Profile | Effect |
|---|---|
| `cardinal` | Slurm executor, account `def-c3g`, high-memory tasks to `cnode-hm-`, Apptainer enabled. |
| `cardinal_local` | Same environment, local executor. For a login node or inside `salloc`. |
| `vcf2maf_host` | Runs VCF2MAF on the host with `mugqic/samtools` and `mugqic/htslib` instead of inside the VEP image, and forces `inhibit_vep = true`. Sets that step to 1 CPU, which is what parsing an existing annotation needs. |
| `narval` | Alliance cluster settings. |
| `no_container` | Use the cluster's R module instead of `pipeline.sif`. |
| `apptainer` / `docker` | Container engine only. |
| `test` | Small resources for a quick check. |
| `debug` | Keeps work directories and prints more. |

---

## Outputs

```
results/<cohort_id>/
    manifest/      manifest, completeness, samples, ambiguities TSVs
    completeness/  heatmap and summary table
    mutations/     oncoplot, TMB plot, gene summary tables
    cnv/           cohort .seg, frequency plot, FGA burden
    expression/    PCA HTML, scores
    fusions/       recurrence plot and table
    report/        <cohort_id>_report.html
results/_provenance/
    mutation_provenance.tsv    which mutation source each patient used
```

`mutation_provenance.tsv` records which mutation source each patient used:
`pcgr_delivered`, `pcgr_regenerated_now`, or `pcgr_regenerated_earlier`.

Historically this column mattered because the two groups were **not
comparable** — regenerated patients were built from the raw ensemble VCF, which
is the unfiltered caller union. That is fixed: `FILTER_ENSEMBLE` now applies the
same two filters GenPipes applied before PCGR, so all patients follow one path.
Verified on a delivered patient, where the reconstruction reproduces the
delivered file exactly:

| stage | variants |
|---|---|
| raw ensemble somatic VCF | 423,764 |
| after `>= 2 callers` | 59,670 |
| after depth/VAF filter | **23,604** |
| delivered PCGR VCF (ground truth) | **23,604** |

Across MoHQ-HM-19 the per-patient MAF counts now form one continuous
distribution (≈9,700–38,700) instead of two groups differing fifteen-fold.

Keep checking this column anyway: it is the audit trail, and if a future change
breaks the equivalence the counts will separate again.

---

## Testing

```bash
./tests/run_tests.sh pipeline.sif
```

Three layers: Python checks on the manifest builder, R unit tests on the shared
library, and a Nextflow stub run. The stub run is the cheapest useful check —
run it after any workflow edit.

---

## Known constraints

**No delivered MAF in the collection is annotated.** Confirmed across all
cohorts and institutions. The pipeline regenerates them; it does not fix the
delivered files.

**About a third of patients have no delivered PCGR VCF.** With
`inhibit_vep = true` they are excluded from VCF2MAF, because there is no `CSQ`
field to parse. `run_pcgr_gapfill` handles them properly, at real compute cost.

**Annotation is judged by distinct gene symbols, not by the `Unknown`
fraction.** These are whole-genome calls, and most somatic variants are
intergenic, where `Unknown` is correct. A healthy patient here runs about 54%
`Unknown` with roughly 4,000 distinct genes. A threshold near 50% rejects valid
data — this was found the hard way, and it is why `max_unknown_frac` is 0.95
and `min_annotated_genes` exists.

**Genome build is read from file headers, never from filenames or reports.**
The readset report states GRCh37 for samples whose headers say GRCh38. The
headers are right.

**PCGR version varies within cohorts.** Recorded per patient and reused during
gap-fill, so reprocessing does not change a patient's classification. Treat it
as a covariate.

**`callable_mb` must match the TMB numerator, not the assay.** The numerator
here is non-synonymous (coding) variants, so the denominator is the coding
footprint (~30–40 Mb) even for whole-genome data. The first real run used 2800
and produced 0.01–0.05 mut/Mb — about 90× too low; the same variants over 30 Mb
give 0.9–4.7, which is a believable range.

**The top genes are filtered.** `exclude_flags` removes FLAGS genes (mucins,
WASH/NBPF/GOLGA families, FLG, HRNR, AHNAK2 and similar) before ranking. Without
it, 25 of the top 25 genes in the first real oncoplot were artefacts of
mismapping in repetitive regions and not one recognisable driver appeared. The
unfiltered ranking is still written to `*_gene_rank_unfiltered.tsv`, so the
decision is auditable.

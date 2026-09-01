# Tutorial — one cohort from start to finish

A complete worked example on **MoHQ-CM-6** (melanoma, 39 patients). Every command
and every number below is from a real run; where output is shown, that is what
was printed.

The README is the reference showing what each parameter does, what each file is. This
is a complete tutorial showing what you should type, in what order, and **what to check before
moving on**.

**Time:** about an hour of wall-clock for a 39-patient cohort, most of it
unattended. Add 20–30 minutes the first time for the clone and the container
build in step 0.

**Source:** <https://github.com/c3g/MoHQ_CLAP>

---

## Step 0 — get the pipeline and build the container

Only needed once per machine.

### Clone

```bash
cd ~
git clone https://github.com/c3g/MoHQ_CLAP.git pipeline
cd pipeline
```

The repository is public: <https://github.com/c3g/MoHQ_CLAP>. It contains the
workflow and the scripts.

### Load the cluster tools

```bash
source bin/load_cardinal_modules.sh     # nextflow, java, apptainer
```

Check what the cluster actually provides before assuming:

```bash
bash bin/check_cardinal.sh
```

### Build the analysis image

The R environment lives in an Apptainer image built from `pipeline.def`:

```bash
apptainer build pipeline.sif pipeline.def
```

Ten to twenty minutes, and it needs internet access so do it on a login node,
not a compute node. It is built `FROM rocker/r-ver:4.4.1` with R packages
installed from a **dated** Posit Package Manager snapshot.

The definition ends with a `%test` block that loads every R package.

Confirm it works:

```bash
apptainer exec pipeline.sif Rscript -e 'library(maftools); library(data.table); cat("ok\n")'
```

Then point the params file at it:

```yaml
container: "/home/<you>/pipeline/pipeline.sif"
```

### Configure rclone for the object store

Harvesting reads from Juno over S3. You need a remote called `juno` in
`~/.config/rclone/rclone.conf`, and that file holds secret keys:

```bash
chmod 600 ~/.config/rclone/rclone.conf
rclone lsd juno:$MOHQ_PROJECT:MOH-Q | head     # should list the institutions
```

### Check the layout

```bash
bash bin/check_layout.sh
```

Confirms every expected script is present and executable. It is worth running after
any clone or copy between machines, since a non-executable `bin/` script fails
at the point Nextflow tries to run it.

---

## Before you start a run

```bash
cd ~/pipeline
source bin/load_cardinal_modules.sh
export MOHQ_PROJECT=<the OpenStack project id>
export NXF_WORK=/home/$USER/MoHQ/work
```

---

## Step 1 — harvest

Look before you transfer. `--dry-run` lists and counts, moves nothing:

```bash
python3 bin/harvest_juno.py \
    --remote juno:$MOHQ_PROJECT:MOH-Q \
    --cohort CM/MoHQ-CM-6 \
    --dest   /home/$USER/MoHQ/collection/ \
    --sets   core signatures \
    --dry-run
```

```
[harvest] MoHQ-CM-6: 39 patients
[harvest] selected 744 / 4,793 objects
[harvest] volume    21.1 GB / 27.9 TB  (0.074%, 1354x smaller)
[harvest] by pattern:
    */reports/pcgr/*_D*.vcf.gz                    25 files     4.5 GB
    */raw_cnv/*.cnvkit.vcf.gz                     37 files     1.5 MB
    */expression/*.abundance_genes.tsv            39 files    60.1 MB
    */variants/*.ensemble.somatic.vt.annot.vcf.gz 38 files     8.8 GB
    ...
```

**Read the per-pattern counts — this is the audit before the audit.** Here,
25 of 39 patients have a delivered PCGR VCF and 38 have an ensemble VCF. That
already tells you 13 patients will need gap-fill and one can't be helped.

### Why `--sets core signatures` and not just `core`

`core` does **not** include `*.ensemble.somatic.vt.annot.vcf.gz`, which is the
input to `FILTER_ENSEMBLE`. Harvest `core` alone and patients without a
delivered PCGR VCF are recorded as `unavailable` rather than `regenerable`,
gap-fill silently has nothing to do, and they drop out of every mutation
analysis while the run reports success.

Then the real transfer:

```bash
python3 bin/harvest_juno.py --remote juno:$MOHQ_PROJECT:MOH-Q \
    --cohort CM/MoHQ-CM-6 --dest /home/$USER/MoHQ/collection/ \
    --sets core signatures --verify
```

`--verify` re-lists the destination and checks every object arrived at the
expected size. It writes a receipt to `.harvest_cache/`.

---

## Step 2 — two checks that take seconds and save hours

### Is there a batch directory?

Some cohorts ship PCGR VCFs *beside* the patient folders rather than inside
them. MoHQ-MU-16 keeps 35 of its 86 there.

```bash
ls ~/MoHQ/collection/CM/MoHQ-CM-6/ | grep -v '^MoHQ-CM-6-'
```

Nothing for CM-6. If it prints a directory of `*.vcf.gz`, add it as
`external_pcgr_dir` in your params — otherwise those patients look regenerable
and you re-run PCGR to produce files the project already shipped.

### Do the PCGR VCFs carry CSQ?

`inhibit_vep: true` tells vcf2maf to reuse PCGR's own VEP annotation instead of
re-running VEP — seconds per patient instead of minutes. It is only valid if the
annotation is there.

```bash
cd ~/MoHQ/collection/CM/MoHQ-CM-6
for f in */reports/pcgr/*_D*.vcf.gz; do
  zcat "$f" 2>/dev/null | awk '/^##INFO=<ID=CSQ/{print "CSQ";f=1;exit}
                               /^#CHROM/{exit} END{if(!f)print "NONE"}'
done | sort | uniq -c
cd ~/pipeline
```

```
     25 CSQ
```

All present, so `inhibit_vep: true` stands. **Mixed output means
`inhibit_vep: false`** — half-annotated and half-re-annotated in one cohort is
exactly the systematic difference the batch-effect analysis exists to detect.

---

## Step 3 — the manifest, and the preflight

The manifest is the pipeline's source of truth. Build it standalone first; it
takes seconds and tells you what the run will do.

```bash
python3 bin/build_manifest.py \
  --root      ~/MoHQ/collection/CM/MoHQ-CM-6 \
  --cohort-id MoHQ-CM-6 \
  --outdir    audit
```

```
[manifest] 39 patients across 1 cohort(s); 37/39 analysis-ready
  pcgr_vcf              25/39
  cna_segments          38/39  <-- REQUIRED
  cnvkit_vcf            37/39  <-- REQUIRED
[manifest] PCGR version:
  mugqic/pcgr/1.0.3     26/39
  mugqic/pcgr/1.4.1     12/39
  unknown                1/39
  !! MORE THAN ONE PCGR VERSION
[manifest] PCGR VCF availability:
  delivered    25/39
  regenerable  13/39
  unavailable   1/39
[manifest] 2 patient(s) excluded from analyses
```

**Three things to read here.**

*Do the numbers add up?* 25 + 13 + 1 = 39. If they don't, something is
mis-parsed.

*How many are analysis-ready?* 37 of 39. A patient missing any **required**
asset is excluded from **every** analysis, not just the one that needs it —
here, two patients lack `cnvkit_vcf`.

*How many PCGR versions?* Three. Output columns and calling behaviour differ
between PCGR majors, so differences between those groups are technical. Set
`pca_colour_by: "pcgr_version"` so you can see it.

Then the preflight, which checks the things that would waste an unattended run:

```bash
bash bin/preflight_gapfill.sh params/mohq_cm_6.yaml
```

It verifies the manifest belongs to **this** cohort, that every PCGR version the
gap-fill needs both loads *and* provides an executable, that the reference
bundle is present, that there is disk headroom, and that the workflow still
wires up under `-stub-run`.

```
  9 checks passed, nothing blocking. Safe to launch overnight.
```

---

## Step 4 — the params file

Copy an existing one and change what is a property of the *cohort*. Everything
that is a property of the *cluster* stays as it is.

```yaml
cohort_dir:  "/home/mslavova/MoHQ/collection/CM/MoHQ-CM-6"
cohort_name: "MoHQ-CM-6"
outdir:      "/home/mslavova/MoHQ/results"
container:   "/home/mslavova/pipeline/pipeline.sif"

inhibit_vep: true            # verified in step 2
callable_mb: 30              # CODING footprint -- see below

pca_colour_by: "pcgr_version"

run_oncoprint: true
run_cnv_frequency: true
run_cnv_burden: true
run_pca: true
run_fusions: true
run_report: true
```

**`callable_mb` is the parameter most likely to be silently wrong.** The
denominator must match what the numerator counts. `run_oncoprint.R` counts
non-synonymous **coding** variants, so the denominator is the coding footprint
(~30 Mb) even though this is whole-genome sequencing. 

---

## Step 5 — run it

```bash
nextflow run . -profile cardinal,vcf2maf_host \
    -params-file params/mohq_cm_6.yaml \
    -c conf/overnight.config \
    --run_pcgr_gapfill true -resume
```

```
[54/9ac3db] FILTER_ENSEMBLE (MoHQ-CM-6-3P3S1)              | 13 of 13 ✔
[9f/5c2cb9] RUN_PCGR (MoHQ-CM-6-3P3S1 [mugqic/pcgr/1.0.3]) | 13 of 13 ✔
[5f/ed6ae1] VCF2MAF (MoHQ-CM-6-17)                         | 37 of 37 ✔
[dc/7345cb] ONCOPRINT (MoHQ-CM-6)                          | 1 of 1 ✔
[95/036b64] COHORT_REPORT (MoHQ-CM-6)                      | 1 of 1 ✔

  status        : completed
  tasks ok      : 54
  tasks failed  : 0
  tasks ignored : 0
```

45 minutes, 17.5 CPU hours.

**Check the counts against the manifest.** `FILTER_ENSEMBLE 13 of 13` matches
the 13 regenerable patients; `VCF2MAF 37 of 37` matches analysis-ready. A
smaller number means patients were dropped.

**`tasks ignored` must be 0.** An ignored task is a patient or a figure missing
from the results with the run still reporting success.

---

## Step 6 — read the results

```
results/MoHQ-CM-6/
├── manifest/     the manifest, completeness matrix, ambiguities
├── mutations/    oncoplot, TMB, gene rankings, pathways, interactions
├── cnv/          genome-wide frequency, FGA burden, the cohort .seg
├── expression/   PCA and scree plot
├── fusions/      recurrent fusions
└── report/       MoHQ-CM-6_report.html   <- start here
```

The HTML report is self-contained with figures embedded and no external files so it
can be sent to a researcher as one attachment.

**Open it and read the eligibility table first.** It says how many patients each
analysis could use and why, because the denominators differ between
panels. Then each figure carries a note on how to read it and what would make it
misleading.

Two sanity checks worth doing every time:

```bash
# variant counts should form ONE distribution across mutation sources.
# Two groups differing ~15x means the wrong VCF was fed to PCGR.
cut -f3 ~/MoHQ/results/_provenance/mutation_provenance.tsv | sort | uniq -c

# do the top genes look like this tumour type?
head -12 ~/MoHQ/results/MoHQ-CM-6/mutations/MoHQ-CM-6_oncoprint_gene_summary.tsv | column -t
```

For CM-6, **BRAF appears in 19 of 37 patients** — about 50%, which is the
published melanoma rate. Recovering the canonical driver of a tumour type the
pipeline was never tuned for is the strongest evidence the chain is correct.

---

## Step 7 — iterating on one analysis

Do **not** re-run the pipeline to tweak a figure. Run the script standalone:

```bash
bash bin/test_oncoprint.sh MoHQ-CM-6
bash bin/test_oncoprint.sh MoHQ-CM-6 assets/panel.txt   # with an allow-list
```

It writes to `~/MoHQ/oncoprint_test/`, so it cannot touch real results, and
errors land in your terminal instead of a work directory.

---

## When it goes wrong

Real failures from real runs, with what they actually mean.

### `exit status (141)` — intermittent, some patients pass

SIGPIPE (128 + 13). A pipeline where the reader exits early (`head -1`,
`grep -q`) kills the writer, and `set -o pipefail` turns that into a task
failure. **Intermittent is the tell**: whether the signal fires depends on
whether the writer filled the pipe buffer first, which depends on file size.

### `exit status (3)` from BUILD_MANIFEST

`[error] N manifest path(s) are not absolute`. A relative path in the manifest
resolves against whatever directory the reader happens to be in. The message
names the offending patient and column. Usually an `external_pcgr_dir` given
relative; make it absolute.

### `status : FAILED (nothing ran)`

Everything was cached and nothing needed doing — the run is fine. Distinguish it
from the real failure by `tasks cached`: if that is non-zero, you are up to date.

### A figure or table missing from the report, no error

Two usual causes. The output was never staged into the report (check the
`ch_figures` mixes in `main.nf`), or an R chunk built a widget that was not the
**last** expression in its chunk. Only the last is auto-printed, so anything
after a `dt_show()` silently discards it.

### `no coding lengths could be read` about a file you can `ls`

The container cannot see it. `/cvmfs` needs binding into the image
(`--bind /cvmfs:/cvmfs:ro`); without it `file.exists()` is FALSE inside the
container for a path that is plainly on disk. The PCGR bundle and VEP cache live
there too.

### `-resume` re-runs something you did not change

Scripts in `bin/` are put on `PATH`, not passed as inputs, so editing one does
not always invalidate the cache. The analysis processes stage their R script as
a declared input specifically to force invalidation. If a task re-runs
unexpectedly, an upstream file moved. For example gap-filled PCGR VCFs being
read from the published directory instead of the task work directory on the
second run. It settles after one cycle.

---

## What to do next

- **A second cohort.** Most of what the audit found came from running an
  unfamiliar cohort and looking at what changed.
- **`bin/check_genome_build.sh`** — every patient's build currently comes from
  the *filename*, and this collection is documented to get that wrong.
- **The collection tier** (`run_collection: true` plus a `cohort_map.tsv`) —
  needs the same cancer type at two or more institutions to separate technical
  from biological effects. See `docs/BREAST_COLLECTION_RUNBOOK.md`.

#!/bin/bash
# ===========================================================================
# setup_vep.sh -- install VEP (container + cache + FASTA) on Cardinal
#
#   bash bin/setup_vep.sh                 # check what's needed, download nothing
#   bash bin/setup_vep.sh --install       # actually download (~30 GB, slow)
#
# Based on your colleague's recipe, with four changes that matter at scale:
#
#  1. VERSION PINNED. `docker://ensemblorg/ensembl-vep` (no tag) pulls whatever
#     is newest today. Rebuild in six months and you get a different VEP, which
#     for a multiyear project means annotations that silently stop matching.
#     Worse: VEP and its cache MUST be the same release. An unpinned image plus
#     a pinned cache is a mismatch waiting to happen.
#     Default here is 115, matching the cache path in your Narval scripts
#     (.../homo_sapiens/115_GRCh38/...), so new results stay comparable to old.
#
#  2. RESUMABLE + IDEMPOTENT. The cache is ~25-30 GB. Re-running does not
#     re-download what is already there.
#
#  3. BIND MOUNTS HANDLED. Apptainer auto-mounts $HOME but not arbitrary paths.
#     Installing to /project or /scratch without --bind gives confusing
#     "cannot write" errors from inside the container.
#
#  4. VERIFIED. Ends by annotating a real test VCF. An install that downloaded
#     but cannot annotate is worse than no install, because you find out
#     halfway through a cohort run.
#
# RUN THIS IN A JOB OR tmux, not bare on the login node -- it is a long download.
# ===========================================================================
set -uo pipefail

# --- configuration ---------------------------------------------------------
VEP_RELEASE="${VEP_RELEASE:-115}"
# Space-separated. GRCh38 only -- and that is now VERIFIED, not assumed:
#
#   * 27,268 objects named grch38, zero named grch37
#   * BAM headers checked on CM/MoHQ-CM-1 (5 patients) and MU/MoHQ-MU-8:
#     all GRCh38, by both chr1 length (248,956,422) and reference path
#     (/cvmfs/soft.mugqic/.../Homo_sapiens.GRCh38/genome/bwa_index/...)
#   * MoHQ-MU-8-2 is a sample the readset report CLAIMED was GRCh37; its header
#     says GRCh38. The report's `Reference` column is unreliable.
#
# Do NOT add GRCh37 "just in case" -- it is a separate ~15-25 GB download and,
# on current evidence, would never be used. Re-check with
# bin/check_genome_build.sh if a new institution is added to the collection.
ASSEMBLY="${ASSEMBLY:-GRCh38}"
SPECIES="${SPECIES:-homo_sapiens}"
VEP_ROOT="${VEP_ROOT:-$HOME/MoHQ/VEP}"
SIF="${VEP_SIF:-$VEP_ROOT/vep_${VEP_RELEASE}.sif}"
CACHE="$VEP_ROOT/vep_data"

DO_INSTALL=0
[ "${1:-}" = "--install" ] && DO_INSTALL=1

ok()   { printf '  \033[32mOK\033[0m    %s\n' "$1"; }
warn() { printf '  \033[33mTODO\033[0m  %s\n' "$1"; }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; }
hdr()  { printf '\n\033[1m=== %s ===\033[0m\n' "$1"; }

echo "VEP release : $VEP_RELEASE"
echo "assembly    : $ASSEMBLY"
echo "install to  : $VEP_ROOT"
echo "mode        : $([ $DO_INSTALL -eq 1 ] && echo INSTALL || echo 'check only (--install to download)')"

module load StdEnv/2023 apptainer/1.4.5 2>/dev/null || module load apptainer 2>/dev/null
command -v apptainer >/dev/null || { bad "apptainer not available"; exit 1; }

mkdir -p "$VEP_ROOT" "$CACHE"

# --------------------------------------------------------------------------
hdr "1. VEP container"
if [ -s "$SIF" ]; then
    ok "image present: $SIF"
else
    warn "image missing: $SIF"
    if [ $DO_INSTALL -eq 1 ]; then
        echo "  pulling docker://ensemblorg/ensembl-vep:release_${VEP_RELEASE}.0 ..."
        if ! apptainer pull --name "$SIF" \
                "docker://ensemblorg/ensembl-vep:release_${VEP_RELEASE}.0"; then
            bad "pull failed"
            cat <<EOF

  If this failed on NETWORK, SD4H may block Docker Hub. Options:
    a) pull it on Narval and copy the .sif across (one file, ~1-2 GB)
    b) ask SD4H support whether a registry mirror is available
    c) check whether mugqic/pcgr already ships a usable VEP:
         module load mugqic/pcgr/2.1.2 && which vep

  If it failed because the TAG does not exist, list what does:
    apptainer pull --name /tmp/t.sif docker://ensemblorg/ensembl-vep:release_${VEP_RELEASE}.0
  and adjust VEP_RELEASE.
EOF
            exit 1
        fi
        ok "pulled $SIF"
    fi
fi

# Bind the install root explicitly: apptainer auto-mounts \$HOME but not
# /project or /scratch, and a missing bind looks like a permissions bug.
BIND="--bind $VEP_ROOT:$VEP_ROOT"

# --------------------------------------------------------------------------
hdr "2. cache + reference FASTA"
declare -A FASTA_FOR
for asm in $ASSEMBLY; do
    CACHE_DIR="$CACHE/$SPECIES/${VEP_RELEASE}_${asm}"
    if [ -d "$CACHE_DIR" ] && [ -n "$(ls -A "$CACHE_DIR" 2>/dev/null)" ]; then
        ok "$asm cache present: $CACHE_DIR ($(du -sh "$CACHE_DIR" 2>/dev/null | cut -f1))"
    else
        warn "$asm cache missing: $CACHE_DIR  (~15-30 GB download)"
        if [ $DO_INSTALL -eq 1 ] && [ -s "$SIF" ]; then
            echo "  running INSTALL.pl -a cf for $asm ..."
            # -a cf = auto-install Cache and Fasta.
            # --NO_UPDATE stops it offering to upgrade past the pinned image,
            # which would break the release/cache version match.
            apptainer exec $BIND "$SIF" INSTALL.pl \
                --CACHEDIR "$CACHE" \
                --AUTO cf \
                --SPECIES "$SPECIES" \
                --ASSEMBLY "$asm" \
                --CACHE_VERSION "$VEP_RELEASE" \
                --NO_UPDATE --NO_TEST || { bad "INSTALL.pl failed for $asm"; exit 1; }
            ok "$asm cache installed"
        fi
    fi
    f=$(ls "$CACHE_DIR"/*.fa.gz 2>/dev/null | head -1)
    FASTA_FOR[$asm]="$f"
    [ -n "$f" ] && ok "$asm FASTA: $f" || warn "$asm: no .fa.gz under $CACHE_DIR"
done

# Primary assembly = the first listed; used for the params summary below.
PRIMARY=$(echo $ASSEMBLY | awk '{print $1}')
CACHE_DIR="$CACHE/$SPECIES/${VEP_RELEASE}_${PRIMARY}"
FASTA="${FASTA_FOR[$PRIMARY]:-}"

# --------------------------------------------------------------------------
hdr "3. vep_wrapper for vcf2maf"
# vcf2maf shells out to a `vep` EXECUTABLE. With a containerised VEP it needs a
# small wrapper on disk, which is what --vep-path points at.
WRAPPER_DIR="$VEP_ROOT/bin"
mkdir -p "$WRAPPER_DIR"
cat > "$WRAPPER_DIR/vep" <<EOF
#!/bin/bash
# Auto-generated by setup_vep.sh -- lets vcf2maf call containerised VEP.
exec apptainer exec --bind $VEP_ROOT:$VEP_ROOT \\
     ${SIF} vep "\$@"
EOF
chmod +x "$WRAPPER_DIR/vep"
ok "wrapper: $WRAPPER_DIR/vep   (use as --vep-path $WRAPPER_DIR)"

# --------------------------------------------------------------------------
hdr "4. vcf2maf"
V2M="${VCF2MAF_DIR:-$HOME/tools/vcf2maf}"
if [ -f "$V2M/vcf2maf.pl" ]; then
    ok "vcf2maf: $V2M/vcf2maf.pl"
else
    warn "vcf2maf missing at $V2M/vcf2maf.pl"
    if [ $DO_INSTALL -eq 1 ]; then
        mkdir -p "$(dirname "$V2M")"
        git clone --depth 1 https://github.com/mskcc/vcf2maf.git "$V2M" \
            && ok "cloned vcf2maf" \
            || bad "clone failed (network?) - copy it from Narval instead"
    fi
fi

# --------------------------------------------------------------------------
hdr "5. verify it actually annotates"
if [ -s "$SIF" ] && [ -d "$CACHE_DIR" ] && [ -n "$(ls -A "$CACHE_DIR" 2>/dev/null)" ]; then
    TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
    cat > "$TMP/test.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=17>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
17	7676154	.	G	A	.	.	.
EOF
    if apptainer exec $BIND --bind "$TMP:$TMP" "$SIF" vep \
            --input_file "$TMP/test.vcf" --output_file "$TMP/out.txt" \
            --cache --dir_cache "$CACHE" --offline \
            --assembly "$PRIMARY" --cache_version "$VEP_RELEASE" \
            --symbol --force_overwrite >/dev/null 2>"$TMP/err"; then
        if grep -q 'TP53' "$TMP/out.txt" 2>/dev/null; then
            ok "annotated a known TP53 position correctly"
        else
            warn "VEP ran but produced no TP53 symbol - check --symbol / cache version"
        fi
    else
        bad "VEP failed to run:"; tail -8 "$TMP/err"
    fi
else
    warn "skipping verification until the image and cache are installed"
fi

# --------------------------------------------------------------------------
hdr "add to your params file"
cat <<EOF
inhibit_vep:  false
vcf2maf_path: "$V2M/vcf2maf.pl"
vep_cache:    "$CACHE"
vep_path:     "$WRAPPER_DIR"
ref_fasta:    "${FASTA:-<set once the FASTA is downloaded>}"
vep_container: "$SIF"
EOF

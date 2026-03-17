#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-./ld_vcf}"
mkdir -p "$OUTDIR"

# Primary: IGSR/EBI (GRCh38 phased biallelic SNV+INDEL; 2019-03-12 release)
BASE1="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL"

# Fallback mirror: UCSC (same filenames; handy if EBI is flaky)
BASE2="https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes"

# Helper: download URL -> output file (resume supported)
dl () {
  local url="$1"
  local out="$2"
  if [[ -f "$out" ]]; then
    echo "Already exists: $(basename "$out")"
    return 0
  fi
  echo "Downloading: $url"
  wget -c -O "$out" "$url"
}

# Helper: choose a working base by checking chr1 VCF existence
pick_base () {
  local testfile="ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
  if wget --spider -q "${BASE1}/${testfile}"; then
    echo "$BASE1"
  elif wget --spider -q "${BASE2}/${testfile}"; then
    echo "$BASE2"
  else
    echo "ERROR: Could not find ${testfile} at either base URL." >&2
    echo "Tried:" >&2
    echo "  ${BASE1}/${testfile}" >&2
    echo "  ${BASE2}/${testfile}" >&2
    exit 1
  fi
}

BASE="$(pick_base)"
echo "Using base URL: $BASE"
echo "Output dir: $OUTDIR"
echo

for CHR in $(seq 1 22); do
  VCF="ALL.chr${CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
  TBI="${VCF}.tbi"

  echo "=== Chromosome ${CHR} ==="
  dl "${BASE}/${VCF}" "${OUTDIR}/${VCF}"
  dl "${BASE}/${TBI}" "${OUTDIR}/${TBI}"
  echo
done

echo "Done."

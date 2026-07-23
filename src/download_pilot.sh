#!/usr/bin/env bash
# Resumable, idempotent pilot download: chr21+chr22 for 4 high-coverage archaics + masks + moderns + ancestral.
# Safe to re-run: skips files whose local size already matches the remote Content-Length.
set -u
BASE="https://cdna.eva.mpg.de/neandertal"
ROOT="C:/Users/benne/Desktop/Archaic Genomics Pipeline/superarchaic-neanderthal/data"
RAW="$ROOT/raw"; MASK="$ROOT/masks"; MOD="$ROOT/modern"; ANC="$ROOT/ancestral"
mkdir -p "$RAW" "$MASK" "$MOD" "$ANC"
LOG="$ROOT/../results/logs/download_pilot.log"
mkdir -p "$(dirname "$LOG")"
say(){ echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

# fetch <url> <outfile>  — resume; skip when local size == remote size.
fetch(){
  local url="$1" out="$2"
  local remote local
  remote=$(curl -sIL --max-time 60 "$url" 2>/dev/null | grep -i '^content-length' | tail -1 | tr -d '\r' | awk '{print $2}')
  if [[ -f "$out" && -n "${remote:-}" ]]; then
    local=$(stat -c %s "$out" 2>/dev/null || echo 0)
    if [[ "$local" == "$remote" ]]; then say "SKIP (complete) $(basename "$out")"; return 0; fi
  fi
  say "GET  $(basename "$out")  (remote=${remote:-?} bytes)"
  curl -fL -C - --retry 5 --retry-delay 5 --max-time 5400 -o "$out" "$url" \
    && say "DONE $(basename "$out")" \
    || say "FAIL $(basename "$out")  <-- $url"
}

say "=== PILOT DOWNLOAD START (chr21, chr22) ==="

for C in 21 22; do
  # --- Archaic all-sites VCFs ---
  fetch "$BASE/altai/AltaiNeandertal/VCF/AltaiNea.hg19_1000g.$C.mod.vcf.gz" "$RAW/AltaiNea.$C.vcf.gz"
  fetch "$BASE/altai/Denisovan/DenisovaPinky.hg19_1000g.$C.mod.vcf.gz"       "$RAW/Denisova.$C.vcf.gz"
  fetch "$BASE/Vindija/VCF/Vindija33.19/chr${C}_mq25_mapab100.vcf.gz"        "$RAW/Vindija33.19.$C.vcf.gz"
  fetch "$BASE/Chagyrskaya/VCF/chr${C}.noRB.vcf.gz"                          "$RAW/Chagyrskaya.$C.vcf.gz"
  # --- Callability masks (BED of callable regions) ---
  fetch "$BASE/Vindija/FilterBed/Altai/chr${C}_mask.bed.gz"        "$MASK/Altai.chr$C.mask.bed.gz"
  fetch "$BASE/Vindija/FilterBed/Denisova/chr${C}_mask.bed.gz"     "$MASK/Denisova.chr$C.mask.bed.gz"
  fetch "$BASE/Vindija/FilterBed/Vindija33.19/chr${C}_mask.bed.gz" "$MASK/Vindija33.19.chr$C.mask.bed.gz"
  fetch "$BASE/Chagyrskaya/FilterBed/chr${C}_mask.bed.gz"          "$MASK/Chagyrskaya.chr$C.mask.bed.gz"
done

# --- Modern humans: 1000G phase3 chr21 (chr22 already present in ../vcf) ---
fetch "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
      "$MOD/ALL.chr21.phase3.vcf.gz"

# --- Ancestral (Ensembl EPO 6-primate, GRCh37) — whole-genome tarball (non-fatal if slow) ---
fetch "http://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2" \
      "$ANC/homo_sapiens_ancestor_GRCh37_e71.tar.bz2"

say "=== PILOT DOWNLOAD END ==="
say "Local sizes:"; du -h "$RAW"/*.vcf.gz "$MASK"/*.bed.gz "$MOD"/*.vcf.gz 2>/dev/null | tee -a "$LOG"

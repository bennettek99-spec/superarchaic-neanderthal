#!/usr/bin/env bash
# Laptop-friendly staggered acquisition: process chromosomes one at a time, and
# within a chromosome one genome at a time — download the big VCF, extract its small
# variant cache, VERIFY the cache, then DELETE the VCF. Peak disk = a single VCF
# (~6 GB worst case, chr1/2), never the ~250 GB full download.
#
# Safe & resumable: a genome is skipped once its cache exists; a raw VCF is deleted
# ONLY after its cache is confirmed non-empty; interrupted downloads resume (curl -C -).
#
# Usage:
#   bash src/stagger.sh 20 19 18            # process these chromosomes, in order
#   bash src/stagger.sh $(seq 1 22)         # whole autosome set (deletes as it goes)
#   KEEP_VCF=1 bash src/stagger.sh 20       # keep VCFs (debug); default deletes them
set -u
BASE="https://cdna.eva.mpg.de/neandertal"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW="$ROOT/data/raw"; MASK="$ROOT/data/masks"; MOD="$ROOT/data/modern"; CACHE="$ROOT/data/cache"
PY=/c/Users/benne/AppData/Local/Python/pythoncore-3.14-64/python.exe
mkdir -p "$RAW" "$MASK" "$MOD" "$CACHE"
LOG="$ROOT/results/logs/stagger.log"; mkdir -p "$(dirname "$LOG")"
say(){ echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
KEEP_VCF="${KEEP_VCF:-0}"

# genome -> "url_dir/filename_template"   ({c} = chromosome)
declare -A URL=(
  [AltaiNea]="altai/AltaiNeandertal/VCF/AltaiNea.hg19_1000g.{c}.mod.vcf.gz"
  [Denisova]="altai/Denisovan/DenisovaPinky.hg19_1000g.{c}.mod.vcf.gz"
  [Vindija33.19]="Vindija/VCF/Vindija33.19/chr{c}_mq25_mapab100.vcf.gz"
  [Chagyrskaya]="Chagyrskaya/VCF/chr{c}.noRB.vcf.gz"
)
declare -A MASKURL=(
  [Altai]="Vindija/FilterBed/Altai/chr{c}_mask.bed.gz"
  [Denisova]="Vindija/FilterBed/Denisova/chr{c}_mask.bed.gz"
  [Vindija33.19]="Vindija/FilterBed/Vindija33.19/chr{c}_mask.bed.gz"
  [Chagyrskaya]="Chagyrskaya/FilterBed/chr{c}_mask.bed.gz"
)
ORDER=(Denisova AltaiNea Vindija33.19 Chagyrskaya)

fetch(){ # url out  — resume; abort+retry a stalled connection (no 33-min hangs)
  local url="$1" out="$2" rem loc
  rem=$(curl -sIL --max-time 60 "$url" 2>/dev/null | grep -i '^content-length' | tail -1 | tr -d '\r' | awk '{print $2}')
  if [[ -f "$out" && -n "${rem:-}" ]]; then
    loc=$(stat -c %s "$out" 2>/dev/null || echo 0)
    [[ "$loc" == "$rem" ]] && return 0
  fi
  curl -fL -C - --retry 8 --retry-delay 5 --retry-all-errors \
    --connect-timeout 30 --speed-limit 50000 --speed-time 60 \
    --max-time 7200 -o "$out" "$url"
}

cache_ok(){ [[ -f "$1" ]] && [[ "$(zcat "$1" 2>/dev/null | head -2 | wc -l)" -ge 1 ]]; }

for C in "$@"; do
  t0=$(date +%s)
  say "===== chr$C START ====="
  # --- archaics: download -> extract -> delete, one genome at a time ---
  for G in "${ORDER[@]}"; do
    cfile="$CACHE/$G.chr$C.variants.tsv.gz"
    if cache_ok "$cfile"; then say "  chr$C $G: cache present, skip"; continue; fi
    vcf="$RAW/$G.$C.vcf.gz"
    tmpl="${URL[$G]}"; url="$BASE/${tmpl/\{c\}/$C}"
    say "  chr$C $G: downloading $(basename "$url")"
    if ! fetch "$url" "$vcf" || ! gzip -t "$vcf" 2>/dev/null; then
      say "  chr$C $G: DOWNLOAD/gzip FAILED (partial kept for resume)"; continue; fi
    say "  chr$C $G: extracting cache"
    PYTHONIOENCODING=utf-8 "$PY" "$ROOT/src/extract_variants.py" --chroms "$C" --genomes "$G" >>"$LOG" 2>&1
    if cache_ok "$cfile"; then
      [[ "$KEEP_VCF" == "1" ]] || { rm -f "$vcf"; say "  chr$C $G: cached OK, deleted VCF ($(basename "$vcf"))"; }
    else
      say "  chr$C $G: EXTRACT FAILED — keeping VCF for inspection"
    fi
  done
  # --- masks (tiny, keep) ---
  for M in Altai Denisova Vindija33.19 Chagyrskaya; do
    out="$MASK/$M.chr$C.mask.bed.gz"
    [[ -f "$out" ]] || { tmpl="${MASKURL[$M]}"; fetch "$BASE/${tmpl/\{c\}/$C}" "$out" && say "  chr$C mask $M ok"; }
  done
  # --- modern 1000G: download -> extract -> delete ---
  mcache="$CACHE/1000G.chr$C.freq.tsv.gz"
  if cache_ok "$mcache"; then
    say "  chr$C 1000G: cache present, skip"
  else
    mvcf="$MOD/ALL.chr$C.phase3.vcf.gz"
    rm -f "$CACHE/1000G.chr$C.freq.tsv.tmp.gz"   # clear any stale partial cache
    murl="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$C.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    if fetch "$murl" "$mvcf" && gzip -t "$mvcf" 2>/dev/null; then
      say "  chr$C 1000G: extracting cache"
      PYTHONIOENCODING=utf-8 "$PY" "$ROOT/src/extract_modern.py" --chroms "$C" --vcf "$mvcf" >>"$LOG" 2>&1
      if cache_ok "$mcache"; then
        [[ "$KEEP_VCF" == "1" ]] || { rm -f "$mvcf"; say "  chr$C 1000G: cached OK, deleted VCF"; }
      else
        say "  chr$C 1000G: EXTRACT FAILED"
      fi
    else
      say "  chr$C 1000G: DOWNLOAD/gzip FAILED (partial kept for resume)"
    fi
  fi
  dt=$(( $(date +%s) - t0 ))
  disk=$(du -sh "$ROOT/data" 2>/dev/null | awk '{print $1}')
  say "===== chr$C DONE in $((dt/60))m$((dt%60))s | data/ now $disk ====="
done
say "ALL REQUESTED CHROMOSOMES PROCESSED. Caches in data/cache/ (small); raw VCFs deleted."

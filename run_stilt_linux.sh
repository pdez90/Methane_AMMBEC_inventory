#!/usr/bin/env bash
# run_stilt_linux.sh — run prepared STILT receptor directories with a Linux hycs_std ------------
# For METHANE_STILT_ENGINE=external: 03_footprints.R writes <FOOT_DIR>/<id>/{CONTROL,SETUP.CFG,
# ASCDATA.CFG,*.ASC}; this script runs each one that has no particle output yet and leaves
# PARTICLE_STILT.DAT.gz behind, which `Rscript scripts/03_footprints.R --collect` turns into foot.rds.
#
#   bash run_stilt_linux.sh <FOOT_DIR> <hycs_std> [n_parallel] [max_runs] [path_from] [path_to]
#
# path_from/path_to rewrite the met-file prefix in CONTROL when the folder is mounted at a
# different path on the Linux machine (e.g. /Users/priyanka/Downloads -> $HOME/mnt).
# Each receptor runs in a scratch copy (TMPDIR), so a killed run never leaves partial output.
set -u
FOOT=${1:?FOOT_DIR}; EXE=${2:?hycs_std}; NP=${3:-4}; MAX=${4:-100000}; FROM=${5:-}; TO=${6:-}
run1() {
  d=$1; id=$(basename "$d"); w=$(mktemp -d); cp "$d"/CONTROL "$d"/SETUP.CFG "$d"/ASCDATA.CFG "$w"/
  cp "$d"/*.ASC "$w"/ 2>/dev/null
  [ -n "$FROM" ] && sed -i "s#$FROM#$TO#" "$w"/CONTROL
  cp "$EXE" "$w"/hycs_std
  (cd "$w" && ./hycs_std > stilt.log 2>&1)
  if [ -f "$w"/PARTICLE_STILT.DAT ] && [ "$(wc -l < "$w"/PARTICLE_STILT.DAT)" -gt 10 ]; then
    gzip -c "$w"/PARTICLE_STILT.DAT > "$d"/PARTICLE_STILT.DAT.gz.tmp && mv "$d"/PARTICLE_STILT.DAT.gz.tmp "$d"/PARTICLE_STILT.DAT.gz
    cp "$w"/stilt.log "$d"/stilt_linux.log; echo "ok $id"
  else echo "FAIL $id: $(tail -2 "$w"/stilt.log | tr '\n' ' ')"; cp "$w"/stilt.log "$d"/stilt_linux.log; fi
  rm -rf "$w"
}
export -f run1; export EXE FROM TO
pending=$(for d in "$FOOT"/*/; do [ -f "$d/CONTROL" ] && [ ! -f "$d/foot.rds" ] && [ ! -f "$d/PARTICLE_STILT.DAT.gz" ] && echo "${d%/}"; done)
echo "pending: $(echo "$pending" | grep -c . )"
echo "$pending" | grep . | head -n "$MAX" | xargs -P "$NP" -I{} bash -c 'run1 "$@"' _ {}

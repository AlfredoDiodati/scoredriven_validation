#!/bin/bash
# LC_ALL is fixed so the numbers in the provenance file do not change decimal
# separator with the shell that launched it.
export LC_ALL=C
# The whole 1000 x 1000 experiment on one machine.
#
# bin/abm_system_simulate runs configurations one after another inside a single
# process, which is 147 core-hours for the full design. The machine finishes
# more runs per minute with several of those processes at once, and the measured
# peak on the machine docs/ABM_SYSTEM_SIMULATION.md describes is eight: past
# that the runs compete for the same memory system and throughput falls. This
# splits the design into contiguous shards, one process each, and waits.
#
# Nothing here decides how results are stored. Each process writes compressed
# .npz archives through abm_system_write_batch, ten replications to an archive,
# which is the same layout applications/abm_system_convert_rdata.c produces and
# applications/abm_system_fit_qvarma.c reads.
#
# Resuming is re-running. An archive already on disk is skipped whole, so an
# interrupted run, a re-submitted job and a deliberate extension to more
# replications are the same command.
#
# Usage:
#
#   ./applications/abm_system_simulate_all.sh [SHARDS [FIRST LAST [N_MC]]]
#
# Defaults: 8 shards over configurations 1 to 1000, 1000 replications each.
# A pilot is the same script over a few configurations with a few replications.
#
# ABM_SYSTEM_OUTPUT_DIR chooses where the archives go and defaults to
# dataset/abm_system. The script refuses to start if that directory already
# holds subdirectories some other writer produced, because the fitting stage
# reads every subdirectory it finds and cannot tell two datasets apart.
#
# Writes out/abm_system_simulate_all_provenance.txt and one log per shard under
# out/abm_system_simulate/.
set -u

BIN=./bin/abm_system_simulate
DESIGN=dataset/abm_system_design.csv

# The design decides how many configurations there are, so that number is read
# from it rather than repeated here.
n_design=0
[ -f "$DESIGN" ] && n_design=$(( $(wc -l < "$DESIGN") - 1 ))

SHARDS="${1:-8}"
FIRST="${2:-1}"
LAST="${3:-$n_design}"
N_MC="${4:-1000}"
OUTPUT_DIR="${ABM_SYSTEM_OUTPUT_DIR:-dataset/abm_system}"
REPORT_DIR=out/abm_system_simulate
STAMP=out/abm_system_simulate_all_provenance.txt

[ -x "$BIN" ] || { echo "$BIN is not built: run make bin/abm_system_simulate" >&2; exit 1; }
[ -f "$DESIGN" ] || { echo "$DESIGN is missing: run make app-abm_system_design" >&2; exit 1; }
[ -x model/dsk_sfc/dsk_SFC ] || { echo "the model is not built: run make model" >&2; exit 1; }

if [ "$LAST" -lt "$FIRST" ] || [ "$FIRST" -lt 1 ]; then
  echo "configurations $FIRST to $LAST are not a range" >&2; exit 1
fi
if [ "$SHARDS" -lt 1 ]; then echo "at least one shard" >&2; exit 1; fi

# Every subdirectory of the output directory is a configuration to the fitting
# stage. One left there by the older .Rdata route would be fitted alongside
# these and there would be nothing in the results to say so.
if [ -d "$OUTPUT_DIR" ]; then
  foreign=$(find "$OUTPUT_DIR" -mindepth 1 -maxdepth 1 -type d ! -name 'cop_*' | head -5)
  if [ -n "$foreign" ]; then
    echo "$OUTPUT_DIR already holds directories this experiment did not write:" >&2
    echo "$foreign" | sed 's/^/  /' >&2
    echo >&2
    echo "The fitting stage reads every subdirectory it finds and cannot tell the" >&2
    echo "two datasets apart. Move the old one aside, or point this run elsewhere:" >&2
    echo >&2
    echo "  mkdir -p dataset/abm_system_rdata && mv $OUTPUT_DIR/EstimationSeries* dataset/abm_system_rdata/" >&2
    echo "  ABM_SYSTEM_OUTPUT_DIR=dataset/abm_system_design_runs $0 $*" >&2
    exit 1
  fi
fi

mkdir -p "$REPORT_DIR" "$OUTPUT_DIR"

# The design a stored archive was produced under, recorded beside the archives
# rather than only in the provenance file. The row index is the identity every
# archive, fit and confidence set entry is named by, so a design redrawn with
# --force renumbers all of them, and a later run would append new-parameter data
# into existing cop_NNNN directories with nothing to say so.
design_md5=$(md5sum "$DESIGN" | cut -d' ' -f1)
fingerprint="$OUTPUT_DIR/.design_md5"
if [ -f "$fingerprint" ]; then
  stored=$(cat "$fingerprint")
  if [ "$stored" != "$design_md5" ]; then
    echo "$OUTPUT_DIR was filled under a different parameter design." >&2
    echo "  stored   $stored" >&2
    echo "  $DESIGN  $design_md5" >&2
    echo >&2
    echo "Configuration numbers identify rows of the design, so adding to this" >&2
    echo "directory under a redrawn design would mix two meanings of cop_NNNN." >&2
    echo "Restore the design this directory was built with, or point this run" >&2
    echo "elsewhere with ABM_SYSTEM_OUTPUT_DIR." >&2
    exit 1
  fi
else
  echo "$design_md5" > "$fingerprint"
fi

# Archives already present, so the throughput reported at the end counts what
# this launch produced rather than what was on disk before it started.
archives_before=$(find "$OUTPUT_DIR" -name 'batch_*.npz' | wc -l)

total=$((LAST - FIRST + 1))
if [ "$SHARDS" -gt "$total" ]; then SHARDS="$total"; fi
per=$((total / SHARDS))
remainder=$((total % SHARDS))

START_EPOCH=$(date +%s)
{
  # Named for what this launch actually covers, not for the full design, so a
  # pilot's provenance file does not describe a run that did not happen.
  echo "$total x $N_MC DSK experiment, launched $(date -d @$START_EPOCH '+%F %T')"
  echo
  echo "configurations   $FIRST to $LAST ($total)"
  echo "replications     $N_MC per configuration"
  echo "shards           $SHARDS concurrent processes"
  echo "output           $OUTPUT_DIR, compressed .npz, 10 replications per archive"
  echo
  echo "model            $(readlink -f model/dsk_sfc/dsk_SFC), built $(stat -c %y model/dsk_sfc/dsk_SFC)"
  echo "driver           $(readlink -f $BIN), built $(stat -c %y $BIN)"
  echo "design           $DESIGN  md5 $design_md5"
  echo "base parameters  model/dsk_sfc/dsk_sfc_inputs.json  md5 $(md5sum model/dsk_sfc/dsk_sfc_inputs.json | cut -d' ' -f1)"
  echo "host             $(uname -srm), $(nproc) hardware threads"
  echo "scratch          ${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
  echo
  echo "shard ranges"
} > "$STAMP"

pids=()
next="$FIRST"
for shard in $(seq 1 "$SHARDS"); do
  width="$per"
  if [ "$shard" -le "$remainder" ]; then width=$((per + 1)); fi
  shard_first="$next"
  shard_last=$((next + width - 1))
  next=$((shard_last + 1))

  log="$REPORT_DIR/shard_$(printf '%02d' "$shard").log"
  printf "  shard %2d  configurations %4d to %4d  %s\n" \
         "$shard" "$shard_first" "$shard_last" "$log" >> "$STAMP"

  ABM_SYSTEM_OUTPUT_DIR="$OUTPUT_DIR" "$BIN" "$shard_first" "$shard_last" "$N_MC" > "$log" 2>&1 &
  pids+=($!)
done

# A run of this length outlives any reasonable idle timeout, and this machine is
# configured to suspend on a closed lid. A suspend mid-run loses only the batch
# in flight, but it also stops the clock this file reports.
inhibit=""
if command -v systemd-inhibit >/dev/null 2>&1; then
  systemd-inhibit --what=handle-lid-switch:sleep:idle:shutdown \
                  --who="abm_system_simulate_all" \
                  --why="$total configurations x $N_MC replications" \
                  --mode=block sleep infinity &
  inhibit=$!
fi

failed=0
for pid in "${pids[@]}"; do
  wait "$pid" || failed=$((failed + 1))
done
[ -n "$inhibit" ] && kill "$inhibit" 2>/dev/null

END_EPOCH=$(date +%s)
ELAPSED=$((END_EPOCH - START_EPOCH))
archives=$(find "$OUTPUT_DIR" -name 'batch_*.npz' | wc -l)
written=$((archives - archives_before))
{
  echo
  echo "finished $(date -d @$END_EPOCH '+%F %T')"
  printf "wall clock       %d s (%.2f h)\n" "$ELAPSED" "$(awk -v s=$ELAPSED 'BEGIN{print s/3600}')"
  echo "shards failing   $failed of $SHARDS"
  echo "archives on disk $archives"
  echo "archives written $written by this launch, $archives_before already present"
  echo "stored size      $(du -sh "$OUTPUT_DIR" | cut -f1)"
  if [ "$ELAPSED" -gt 0 ] && [ "$written" -gt 0 ]; then
    # An upper bound: every archive is assumed full at ten replications, which
    # a batch with a failed run is not. out/abm_system_simulate_manifest.txt
    # holds the completed and failed counts per configuration.
    printf "throughput       %.1f runs per minute over this launch, at most\n" \
           "$(awk -v a=$written -v s=$ELAPSED 'BEGIN{print a*10*60/s}')"
  fi
  echo
  echo "completion per configuration is in out/abm_system_simulate_manifest.txt,"
  echo "and the reason for every rejected replication in $REPORT_DIR."
} >> "$STAMP"

cat "$STAMP"
exit $((failed > 0))

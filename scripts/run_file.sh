#!/usr/bin/env sh
set -eu

PROCESSORS="${1:-4}"
INPUT="${2:-examples/input_small.txt}"
OUTPUT="${3:-build/result.txt}"
MPIRUN="${MPIRUN:-}"

if [ -z "$MPIRUN" ]; then
  if [ -x /opt/homebrew/bin/mpirun ]; then
    MPIRUN=/opt/homebrew/bin/mpirun
  else
    MPIRUN=mpirun
  fi
fi

make all
"$MPIRUN" -np "$PROCESSORS" ./bin/parallel_fft_poly --input "$INPUT" --output "$OUTPUT" --print

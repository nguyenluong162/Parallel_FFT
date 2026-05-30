#!/usr/bin/env sh
set -eu

PROCESSORS="${1:-4}"
DEGREE="${2:-16}"
MPIRUN="${MPIRUN:-}"

if [ -z "$MPIRUN" ]; then
  if [ -x /opt/homebrew/bin/mpirun ]; then
    MPIRUN=/opt/homebrew/bin/mpirun
  else
    MPIRUN=mpirun
  fi
fi

make all
"$MPIRUN" -np "$PROCESSORS" ./bin/parallel_fft_poly --demo "$DEGREE" --print

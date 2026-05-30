#!/usr/bin/env sh
set -eu

find include src tests scripts -type f \( -name '*.h' -o -name '*.c' -o -name '*.sh' \) -print0 | xargs -0 wc -l

#!/bin/bash
set -e

BUNDLE_DIR="$1"
OUTPUT_LIB="$2"
AR="$3"

cd "$BUNDLE_DIR"

# Create the combined archive with all object files
find . -name '*.o' -print0 | xargs -0 "$AR" rcs "$OUTPUT_LIB"

echo "Created bundled library: $OUTPUT_LIB"

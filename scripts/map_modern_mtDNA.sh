#!/bin/bash
# Usage: bash script.sh reference read_path output_path

# Variables to change:
threads=32
remove_temp=false

set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 reference.fa read_path output_path output_name" >&2
    exit 1
fi

reference="$1"
read_path="$2"
output_path="$3"
output_name="$4"

#R1 and #R2?

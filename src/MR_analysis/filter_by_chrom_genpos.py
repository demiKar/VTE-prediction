import sys

if len(sys.argv) != 3:
    print("Usage: python filter_by_first_column.py <reference_file> <target_file>")
    sys.exit(1)

reference_file = sys.argv[1]
target_file = sys.argv[2]

# Load CHROM_GENPOS values from reference (first column only)
ref_ids = set()
with open(reference_file) as ref:
    for line in ref:
        if line.strip():
            ref_ids.add(line.strip().split()[0])

# Print lines from target if first column is in ref_ids
with open(target_file) as target:
    for line in target:
        if line.strip() and line.strip().split()[0] in ref_ids:
            print(line, end="")


import sys

distance_threshold = 1000000

selected = {}  # SNP -> BP

with open(sys.argv[1]) as fp:
	fp.readline()
	for line in fp:

		fields = line.strip().split()

		if len(fields) != 0:

			snp = fields[2]
			bp = int(fields[3])
			pval = fields[4]

			# Check distance from all previously selected SNPs
			too_close = any(abs(bp - other_bp) < distance_threshold for other_bp in selected.values())
			if not too_close:
				selected[snp] = bp
				print(f"{snp} {pval}")


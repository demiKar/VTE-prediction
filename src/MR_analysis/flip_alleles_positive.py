import sys

with open(sys.argv[1]) as fp:
	header = fp.readline()
	print(header, end="")

	for line in fp:
		cols = line.strip().split()
		beta = float(cols[6])
		if beta < 0:
			cols[6] = str(-beta)				   # Flip beta
			cols[3], cols[4] = cols[4], cols[3]	# Swap alleles
			af = float(cols[5])
			cols[5] = str(1 - af)				  # Flip allele freq
		print(" ".join(cols))


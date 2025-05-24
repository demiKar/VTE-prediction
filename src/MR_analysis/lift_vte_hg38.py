import sys

chrpos_oldnew = {}

with open(sys.argv[1]) as fp: #hglft_genome_161ba9_3081c0.bed
	for line in fp:
		line = line.rstrip().split()

		chromish, posish = line[3].split(":")

		chrom = chromish.replace("chr", "")
		pos = posish.split("-")[0]

		old_chrompos = chrom + "\t" + pos

		new_chrompos = line[0].replace("chr", "") + "\t" + line[2]

		chrpos_oldnew[old_chrompos] = new_chrompos



with open(sys.argv[2]) as fp: #vte_meta_sorted.txt
	print("CHROM_GENPOS\t" + fp.readline().rstrip())
	for line in fp:
		line = line.rstrip().split()

		chrompos = line[1] + "\t" + line[2]

		toprint = ""
		if chrompos in chrpos_oldnew:
			newchr, newpos = chrpos_oldnew[chrompos].split()
			toprint = newchr + "_" + newpos + "\t" + line[0] + "\t" + chrpos_oldnew[chrompos]
		
			rest = "\t".join(line[3:])
			print(toprint + "\t" + rest)






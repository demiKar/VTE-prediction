import sys

cp_as = {} #chromPos to alleles1and2
with open(sys.argv[1]) as fp: #example_files/STable2.txt
	fp.readline()
	for line in fp:
	        line = line.rstrip().split()
                chromPos = line[1] + "\t" + line[2]
	        cp_as[chromPos] = line[3] + "\t" + line[4]

#positions to use
chrPos = set()
with open(sys.argv[2]) as fp: #example_files/chrPos_toKeep.txt
	for line in fp:
		chrPos.add(line.rstrip())

#parse more args
chr_ind = int(sys.argv[4])
pos_ind = int(sys.argv[5])
a1_ind = int(sys.argv[6])
a2_ind = int(sys.argv[7])
freq_ind = int(sys.argv[8])
beta_ind = int(sys.argv[9])
se_ind = int(sys.argv[10])
p_ind = int(sys.argv[11])

def flip(allele):
        if(allele == "A"):
            return "T"
        elif allele == "T":
            return "A"
        elif allele == "C":
            return "G"
        elif allele == "G":
            return "C"
        else:
            return "ERROR_non-SNP"

with open(sys.argv[3]) as fp: #whatever full GWAS results e.g., example_files/LDL_EUR_formatted.txt
        print("CHROM\tPOS\tOTHER_allele\tEFFECT_allele\tFREQ_EFFECT_allele\tBETA\tSE\tP")
	for line in fp:
	        old_l = line.rstrip()
    	        line = old_l.split()

                chromPos = line[chr_ind] + "\t" + line[pos_ind]

		if chromPos in chrPos:
	       	        a1 = line[a1_ind].upper()
	                a2 = line[a2_ind].upper()
	                a1_f = flip(a1)
	                a2_f = flip(a2)

	                #combine with tabs
	                a1_a2 = a1 + "\t" + a2
	                a2_a1 = a2 + "\t" + a1
	                a1_a2_f = a2_f + "\t" + a1_f
	                a2_a1_f = a2_f + "\t" + a1_f

		        toprint = ""
		        if(cp_as[chromPos] == a1_a2): #alleles same
			        for i in [chr_ind, pos_ind, a1_ind, a2_ind, freq_ind, beta_ind, se_ind, p_ind]:
	     	                        if i == a1_ind:
	                                       	toprint += a1 + "\t"
	                                elif i == a2_ind:
	                                       	toprint += a2  + "\t"
	                                else:
	                                	toprint += line[i] + "\t"
		                print(toprint.rstrip())
	                elif (cp_as[chromPos] == a1_a2_f): #alleles wrong strand
			        for i in [chr_ind, pos_ind, a1_ind, a2_ind, freq_ind, beta_ind, se_ind, p_ind]:
	                        	if i == a1_ind:
	                                        toprint += a1_f + "\t"
	                                elif i == a2_ind:
	                                        toprint += a2_f  + "\t"
	                                else:
	                                        toprint += line[i] + "\t"
		                print(toprint.rstrip())
	                elif(cp_as[chromPos] == a2_a1): #alleles flipped
			        for i in [chr_ind, pos_ind, a1_ind, a2_ind, freq_ind, beta_ind, se_ind, p_ind]:
	                                if i == a1_ind:
	                                        toprint += a2 + "\t"
	                                elif i == a2_ind:
	                                        toprint += a1 + "\t"
	                                elif i == beta_ind:
			                        toprint += str((float(line[beta_ind]) * -1.00)) + "\t"
			                elif(i == freq_ind):
						if line[freq_ind] != "NA":
				               		toprint += str((1.00 - float(line[freq_ind]))) + "\t"
						else:
				               		toprint += "NA\t"
	                                else:
	                                        toprint += line[i] + "\t"
			        print(toprint.rstrip())
	                elif(cp_as[chromPos] == a2_a1_f): #alleles flipped and wrong strand
			        for i in [chr_ind, pos_ind, a1_ind, a2_ind, freq_ind, beta_ind, se_ind, p_ind]:
	    		                if(i == a1_ind):
				                toprint += a2_f + "\t"
			                elif(i == a2_ind):
				                toprint += a1_f + "\t"
			                elif(i == beta_ind):
				                toprint += str((float(line[beta_ind]) * -1.00)) + "\t"
			               	elif(i == freq_ind):
				                toprint += str((1.00 - float(line[freq_ind]))) + "\t"
			                else:
				                toprint += line[i] + "\t"
			        print(toprint.rstrip())
			else:
	                        print("ERROR:wrong_alleles_betwen_height_and_CAD")

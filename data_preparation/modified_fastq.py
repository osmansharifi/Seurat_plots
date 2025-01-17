import csv

#Random sequences that do NOT exist in the mouse genome
mecp2_WT = "TCACCTAAGGGCCCCGGGCTTTATTACTCTCCCCTACGCTAATCGGGTCCTAGTGGGCAATGTCTGCGGCGTACTGA\
TAGAAGTTACAGGTGCGTGGTCCCGAGAGCAATCGCAAAACACGCAGCATTGTAAGGAAGGACAATCGAACAAGCGCGTAGATATAATA\
ACCACGGGCTTGGGCTGCTGGCGTTCGCTACGTTGCTAGGGCCGTTGTGCCACTAGTACTATGACCTTCCCTCGGGTCTGTAGTGTGGA\
GGCGACGACGGGAAGTGACTAAAGAACTGTTTCGATCTTACGGATAGCCTACCCGTCGCCTCCAAGCATAAATATGGGGCGAGGGAAAG\
CGTCCCGCCGCGGTACGAGACACTAGGGAACTTAAGGAAAATGGCAGGGCACTCTCTTAGCTGTATTCAACAGAGTCTATCCGTCCGGC\
GACCTATCGTAACGTGGGTCGCAGAGTAACCATGCATAAATATGCCTGCGGTAACTTCTGACGTTTA"
mecp2_MUT = "ATGTTTCTAGAATATATTTACGAAAAGTATGACGCCCTGTAGTCTCGCTTCTCCGCAACTAAACATACGGGACAGA\
TTGTATGTGGCTAAGTTTACACTGACACCTGCCCCCCCCTCTGGTCGCCCCCGGAGGGAATACGGTCAGCACACCACTGGGCTCATGTG\
AAAGAGCGTTAGCAATTACCCAAATTCGGTGTTACGACTTCAGACTCAAGGAAACAGCTGGCCACGGAAACATGTCGCGATTTAGGCGG\
CCGGTGTTTATCCCGTTCCGGTGTTGTAGATGAGACTGCGCACACACAAAAACCCCTCGACCGGGATTACTCGCTACCCCTCAGTTGTG\
GATGGGTCCCCCAGTGATGCGCCTAGCTTGGAGTCGGGAGCGCGGCCGAGAACGCCTCTTCTTGACTTGTTGAACGGTCACTTAAGACA\
CAACAGTTCGGTACATACTCTGTTTTGGGTGGCCGGTCAATCCTATGACTGGAACCAAATGGCTCAAG"
WT_r1_addon = (mecp2_WT[10:160]) # dont start at 0
WT_r2_addon = (mecp2_WT[340:490])
MUT_r1_addon = (mecp2_MUT[10:160])
MUT_r2_addon = (mecp2_MUT[340:490])

#spoofing function
def spoofreads(bc, WT_r1_addon):




#Read barcodes and glue them to its corresponding random sequence
with open('/Users/osman/Documents/GitHub/snRNA-seq-pipeline/\
data_preparation/fake_barcode_counts.tsv','r') as file:
	barcodes = csv.reader(file, delimiter="\t")
	header = next(barcodes)
	for bc, wt, mut in barcodes:
		wt = int(wt)
		mut = int(mut)
		for i in range(wt):
			print(i)
			#spoofreads(barcode, WT_r1_addon)
		#for i in range(mut):
			#print(i)
			#spoofreads(bc, MUT_r1_addon)
		
			
#create spoof read functions


'''
# example barcode count file
barcode	WT_count	MUT_count
WWWWWWWWWWWWW	3	0
XXXXXXXXXXXXX	0	2
YYYYYYYYYYYYY	3	0
ZZZZZZZZZZZZZ	0	2

# example fastq reads
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:
FFFF: @A00351:327:H2MGWDSXY:1:2375:29722:14888 1:N:0:GAAGGAAC
CTAACTTCAGTATCTGATGAGCGTGGTTTCTTATATGGGGAGGCTGGCTTCGGCTGGACAAGGCTGAGGTGACTGCAGCT
CCTTCTGGAGAGCTTCGGTTTCTCGTATTTGTCTTCTCTACCTTCGCTCGTCCAGCTGTCTTTCGATTTGC +
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
FFFFFFFFFFFFFFFFFFF,FFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFF
@A00351:327:H2MGWDSXY:1:2501:28583:32424 1:N:0:GAAGGAAC
ATAACGCTCTAACCGACATAGGAGCTTTTCTTATATGGGCAGGAAAGAAAATATTTCTTTGCTTTTGTTCCATATAGCAT
GTTTACTTTTACTTTATTTATGCATAGTTGAGGCTATTGATTGGCATTAATTACTTAAAAGCAATGAGTTA


'''





'''
- Input will be 1 file containing barcode counts
- Glue barcode to the correct addon and append it to the r1 and r2 fastqs x number of times
'''
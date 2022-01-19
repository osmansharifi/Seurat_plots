#!/usr/bin/env python
import csv
import sys
import gzip

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
r1_addon = print(mecp2_WT[:150])
r2_addon = print(mecp2_WT[350:])

with open('/Users/osman/Documents/GitHub/snRNA-seq-pipeline/\
data_preparation/fake_barcode_counts.csv','r') as csvfile:
	readit = csv.reader(csvfile)
	for line in readit:
		print(line) # shows barcode file contents


#read fastq files
input_file = ('/Users/osman/Documents/GitHub/snRNA-seq-pipeline/data_preparation/fastq_subset/190106Csub_S9_L001_R1_001.fastq.gz')

if input_file.endswith('.gz'):
	fp = gzip.open(input_file, 'rt')
else:
	fp = open(input_file, 'r')

while True: 
	description = fp.readline()
	sequence = fp.readline() 
	extra1 = fp.readline()
	extra2 = fp.readline()
	#if description == '':
		#break
	#description = description.replace('@', '>', 1)
	print(description + sequence)
	
#fastq_file = open(('/Users/osman/Documents/GitHub/snRNA-seq-pipeline/data_preparation/fastq_subset/190106Csub_S9_L001_R1_001.fastq.gz'), 'r')
#print(fastq_file)
'''
def readFastq(filename):
	fastq_open = open(filename, "rU")     
	for line in fastq_open:
    	header1 = line.rstrip()
        sequences = next(fastq_open).rstrip()
        header2 = next(fastq_open).rstrip()
        qualities = next(fastq_open).rstrip()
        if len(seq) == 0:
        	break
        return sequences, qualities
    fastq_open.close()

seqs, quals = readFastq('/Users/osman/Documents/GitHub/snRNA-seq-pipeline/data_preparation/fastq_subset/190106Csub_S9_L001_R1_001.fastq.gz')
print (seqs)
'''
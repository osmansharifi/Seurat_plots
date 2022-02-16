import argparse
import os
import sys

parser = argparse.ArgumentParser(description='Setup for alleler')
parser.add_argument('fasta', type=str, metavar='<fasta>',
	help='path to Mecp2 fasta file')
parser.add_argument('fq1gz', type=str, metavar='<fq1>',
	help='path to fastq file 1 (gzipped)')
parser.add_argument('fq2gz', type=str, metavar='<fq2>',
	help='path to fastq file 2 (gzipped)')
parser.add_argument('target', type=str, metavar='<target>',
	help='target directory')
arg = parser.parse_args()

fasta = os.path.abspath(arg.fasta)
fq1 = os.path.abspath(arg.fq1gz)
fq2 = os.path.abspath(arg.fq2gz)
target = os.path.abspath(arg.target)

makefile = os.path.abspath('Makefile')
sam2fa = os.path.abspath('sam2fa.pl')
samtrim = os.path.abspath('samtrim.pl')
alleler = os.path.abspath('alleler.py')

if os.path.isdir(target):
	print(f'Directory "{target}" already exists. Remove or change target.')
	sys.exit(1)

os.system(f'mkdir {target}')
os.system(f'ln -s {fasta} {target}/sequence')
os.system(f'ln -s {fq1} {target}/reads1.fq.gz')
os.system(f'ln -s {fq2} {target}/reads2.fq.gz')
os.system(f'ln -s {makefile} {target}/Makefile')
os.system(f'ln -s {sam2fa} {target}/sam2fa.pl')
os.system(f'ln -s {samtrim} {target}/samtrim.pl')
os.system(f'ln -s {alleler} {target}/alleler.py')

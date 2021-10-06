#!/bin/env python3

import argparse
import os
import subprocess
from multiprocessing import Pool

## Command line processing ##

parser = argparse.ArgumentParser(
	description='Local wrapper for cell ranger',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--mem', required=False, type=int, default=62,
	metavar='<int>', help='memory in GB')
parser.add_argument('--cpu', required=False, type=int, default=30,
	metavar='<int>', help='number of threads')
parser.add_argument('--bin', required=False, type=str,
	default='/share/biocore/software/bin/cellranger',
	metavar='<path>', help='path to cell ranger executable directory')
parser.add_argument('--tx', required=False, type=str,
	default='/share/genomes/cellranger_genomes/refdata-cellranger-mm10-3.0.0',
	metavar='<path>', help='path to transcriptome')
parser.add_argument('--fastq', required=True, type=str,
	metavar='<path>', help='path to fastq files')
parser.add_argument('--project', required=True, type=str,
	metavar='<path>', help='path to project directory')
parser.add_argument('--sample', required=True, type=str,
	metavar='<str>', help='name of project or name of file')
parser.add_argument('--job', required=False, type=int, default=1,
	metavar='<int>', help='number of jobs to run concurrently')
parser.add_argument('--test', action='store_true', help='test mode, does not actually run cellranger')
parser.add_argument('--quiet', action='store_true', help='echo parameters')
arg = parser.parse_args()

## File checking ##

if not os.path.exists(arg.bin):
	raise IOError('cellranger executable not found', arg.bin)
if not os.path.exists(arg.tx):
	raise IOError('transcriptome file not found', arg.tx)
if not os.path.exists(arg.fastq):
	raise IOError('fastq directory not found', arg.fastq)
if not os.path.exists(arg.project):
	raise IOError('project directory not found', arg.project)

base = 'time {} count --transcriptome={} --fastqs={} --localcores={} --localmem={}'.format(
	arg.bin, arg.tx, arg.fastq, arg.cpu, arg.mem)

## Verbosity ##

if not arg.quiet:
	print('Memory:', arg.mem)
	print('Threads:', arg.cpu)
	print('Executable:', arg.bin)
	print('Transcriptome:', arg.tx)
	print('Fastq files:', arg.fastq)
	print('Project dir:', arg.project)
	print('Sample/file name:', arg.sample)
	print('Path:', os.environ['PATH'])
	#print('Host:', subprocess.getoutput('hostname'))
	print('User:', os.environ['USER'])

## Threading ##

def worker(name):
	cmd = base + ' --id={} --sample={}'.format(name, name)
	if not arg.quiet:
		print(cmd)
	if not arg.test:
		if os.path.exists(name + '/' + 'outs/metrics_summary.csv'):
			print(name + ' exists, skipping')
		else:
			os.system(cmd)

## Main ##

if os.path.exists(arg.sample):
	files = []
	fp = open(arg.sample, 'r')
	for line in fp:
		files.append(line.rstrip())
	p = Pool(arg.job)
	p.map(worker, files)
else:
	cmd = base + ' --id={} --sample={}'.format(arg.sample, arg.sample)
	if not arg.quiet:
		print(cmd)
	if not arg.test:
		os.system(cmd)




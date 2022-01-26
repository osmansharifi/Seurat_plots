#!/usr/bin/env python3

import argparse
import gzip
import os
import sys

def parse_ab_blast(filename):
	
	fp = open(filename)

	qid, sid = None, None
	qbeg, qend, sbeg, send = None, None, None, None
	score, pct, qstr, sstr = None, None, None, None
	data = False
	
	while True:
		line = fp.readline()
		if line == '':
			break
		elif line.startswith('Query='):
			f = line.split()
			qid, sid = f[1], None
			qbeg, qend, sbeg, send = None, None, None, None
			score, pct, qstr, sstr = None, None, None, None
			data = False
		elif line.startswith('Parameters:'):
			if data:
				yield qid, sid, qbeg, qend, sbeg, send, score, pct, qstr, sstr
				qbeg, qend, sbeg, send = None, None, None, None
				score, pct, qstr, sstr = None, None, None, None
				data = False
		elif line.startswith('>'):
			if data:
				yield qid, sid, qbeg, qend, sbeg, send, score, pct, qstr, sstr
				qbeg, qend, sbeg, send = None, None, None, None
				score, pct, qstr, sstr = None, None, None, None
				data = False
			f = line.split()
			sid = f[0][1:]
		elif line.startswith(' Score ='):
			if data:
				yield qid, sid, qbeg, qend, sbeg, send, score, pct, qstr, sstr
				qbeg, qend, sbeg, send = None, None, None, None
				score, pct, qstr, sstr = None, None, None, None
				data = False
			data = True
			f = line.split()
			score = int(f[2])
			line2 = fp.readline()
			f = line2.split()
			ids = f[2]
			match, total = ids.split('/')
			pct = int(match)/int(total)
			qstr, sstr = '', ''
		elif line.startswith('Query:'):
			f = line.split()
			qstr += f[2]
			if qbeg == None: qbeg = f[1]
			qend = f[3]
		elif line.startswith('Sbjct:'):
			f = line.split()
			sstr += f[2]
			if sbeg == None: sbeg = f[1]
			send = f[3]

	fp.close()


#######
# CLI #
#######

parser = argparse.ArgumentParser(description='Mecp2 allele identifier')
parser.add_argument('blast', type=str, metavar='<blast>',
	help='path to blast report')
arg = parser.parse_args()


#############
# Main Loop #
#############

found = {}
mass = {}
for qid, sid, qbeg, qend, sbeg, send, score, pct, qstr, sstr \
		in parse_ab_blast(arg.blast):
	sbeg = int(sbeg)
	send = int(send)
	bcs = qid.split('.')
	bc = bcs[0]
	umi = bcs[2]
	if 105 >= sbeg and 105 <= send:
		if bc not in found: found[bc] = {}
		if umi not in found[bc]: found[bc][umi] = {'A': 0, 'T': 0}
		if   qstr[105-sbeg] == 'A': found[bc][umi]['A'] += 1
		elif qstr[105-sbeg] == 'T': found[bc][umi]['T'] += 1
		else: sys.stderr.write(f'oddly, found {qstr[105-sbeg]}\n')
	else:
		if bc not in mass: mass[bc] = {}
		if umi not in mass[bc]: mass[bc][umi] = 0
		mass[bc][umi] += 1

for bc in sorted(found):
	for umi in sorted(found[bc]):
		count = mass[bc] if bc in mass else 0
		A = found[bc][umi]['A']
		T = found[bc][umi]['T']
		m = mass[bc][umi] if bc in mass and umi in mass[bc] else 0
		print(bc, umi, A, T, m, sep='\t')




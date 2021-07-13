import korflib
import os
import sys

assert(len(sys.argv) == 2)
def wrap(seq, step=80):
	for i in range(0, len(seq), step):
		print(seq[i:i+step])
for id, seq in korflib.read_fasta(sys.argv[1]):
	if id.startswith("chromX"):
		beg = 10
		end = 30
		mut = 16
		mx = seq[:beg] + "N"*(end-beg) + seq[end:]
		print(f">{id}")
		wrap(mx)
		print(">A")
		wrap(seq[beg:end])
		print(">B")
		wrap(seq[beg:mut] + "t" + seq[mut+1:end])
	else:
		print(f">{id}")
		wrap(seq)
	
import korflib
import sys

#Make sure command line has 2 arguments
assert(len(sys.argv) == 2)

#Make the length of each sequence line 60 bases long
def wrap(seq, step=60):
	for i in range(0, len(seq), step):
		print(seq[i:i+step])
		
#Read in fasta file and change specific sequences into Ns
for id, seq in korflib.read_fasta(sys.argv[1]):
	if id.startswith("X"):
		beg = 74036400
		end = 74085786		#exon beg is 74085691 
		mut = 74085586		#the T at this site was changed to A
		mx = seq[:beg] + "N"*(end-beg) + seq[end:]
		print(f">{id}")
		wrap(mx)
		print(">Mecp2_e1 dna:chromosome chromosome:GRCm38:Mecp2_e1:1:49386:1 REF")
		wrap(seq[beg:end])
		print(">Mecp2_e2 dna:chromosome chromosome:GRCm38:Mecp2_e2:1:49386:1 REF")
		wrap(seq[beg:mut-1] + "a" + seq[mut:end])
	else:
		print(f">{id}")
		wrap(seq)
	
import random
from random import choice

#create 2 random 500bp sequences for the fake mecp2_WT and mecp2_MUT chromosomes

def String(length):
	DNA=""
	random.seed(400)
	for count in range(length):
		DNA+=choice("ACGT")
	return DNA
print(String(500))
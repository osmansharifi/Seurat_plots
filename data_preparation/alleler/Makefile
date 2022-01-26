# Makefile for finding Mecp2 alleles

# Inputs - use aliases to actual files

SEQ = sequence
FQ = reads1.fq.gz reads2.fq.gz

# Alignment Settings

BLASTOPS = M=1 N=-3 S=40 nogap kap wordmask=dust novalidctxok warnings notes
BWAOPS = 

# Main rules

default:
	make $(SEQ).alleler

clean:
	rm -f $(SEQ).*

# BWA Setup

$(SEQ).bwt: $(SEQ)
	bwa index -a is $(SEQ)

# BLAST Setup

$(SEQ).xnt: $(SEQ)
	xdformat -n $(SEQ)

# BWA Alignment

$(SEQ).sam: $(SEQ).bwt $(FQ)
	bwa mem $(SEQ) $(FQ) $(BWAOPS) | perl samtrim.pl > $(SEQ).sam

# Convert aligned sam to fasta

$(SEQ).aln.fa: $(SEQ).sam
	perl sam2fa.pl $(SEQ).sam > $(SEQ).aln.fa

# BLAST Alignment

$(SEQ).blast: $(SEQ).aln.fa $(SEQ).xnt
	ab-blastn $(SEQ) $(SEQ).aln.fa $(BLASTOPS) > $(SEQ).blast

# Allele Counts

$(SEQ).alleler: $(SEQ).blast
	python3 alleler.py $(SEQ).blast > $(SEQ).alleler


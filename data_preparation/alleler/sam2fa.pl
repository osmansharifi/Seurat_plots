use strict;
use warnings;

my %bct;
while (<>) {
	my @f = split;
	my $seq = $f[9];
	my $barcode = substr($seq, 0, 16);
	my $umi = substr($seq, 16, 10);
	my $read = substr($seq, 26);
	$bct{$barcode}++;
	print ">$barcode.$bct{$barcode}.$umi\n$read\n";
}

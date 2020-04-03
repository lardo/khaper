use strict;
use warnings;

die "perl $0 <filt.al> <out.al> <max_cov:30>\n" if(@ARGV==0);

my($infile,$outfile,$max_cov) = @ARGV;
$max_cov ||= 30;

my %hash;
my %hash_cov;

open FI, $infile or die $!;
open FO,">$outfile" or die $!;
my $pre_pb = "";
while (<FI>) {
	my($pb, $contig, $pb_len, $contig_len) = split;
	$hash_cov{$contig} += $pb_len;
	# die "$pb, $contig, $pb_len, $contig_len\n";
	next if ($hash_cov{$contig} > $contig_len*$max_cov);
	$hash{$contig} .= "$pb ";
}
while (my ($ref,$pac) = each %hash){
	print FO ">$ref\n";
	my @pac_array = split(/\s+/,$pac);
	foreach my $pac_content (@pac_array){
		print FO "$pac_content\n";
	}
}

close FI;
close FO;

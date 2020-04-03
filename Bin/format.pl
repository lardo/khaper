use strict;
use warnings;

die "perl $0 <kmer.histo> <outFile> <number of kmer_kind:from *.stats District>
This script is used for formating the data for kmer analysis.
example:
perl $0 17mer.freq kmer.xls 1000000\n" if(@ARGV==0);

my($file,$out,$node)= @ARGV;

open FI,"$file" or die $!;
open FO,">$out";
my $line = 0;
while(<FI>){
	$line++;
	chomp;
	my($count,$num) = split;
	my $percent = $num/$node;
	print FO "$count\t$percent\t$num\n";
	last if($line>=800);
}
close FI;
close FO;

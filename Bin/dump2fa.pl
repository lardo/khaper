#!/usr/bin/perl
use warnings;
use strict;
use FindBin qw($Bin);
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path getcwd);

BEGIN { 
    push (@INC,"$Bin"); 
} 
use SeqIO;

die "perl $0 <kmer.dump> <min_occ> <max_occ>" if(@ARGV!=3);

my($kmer_file,$min_occ,$max_occ) = @ARGV;
$max_occ ||= 0;

my $in_hdl = myOpen($kmer_file);
my $index = 1;
while (<$in_hdl>) {
	chomp;
	my($kmer,$freq) = split;
    next if(!defined $freq);
	next if($freq<$min_occ);
	my $kmer_r = revCom($kmer);

	if($max_occ<=0){
		print ">read_$index\n$kmer\n";
		print ">read_$index\n$kmer_r\n";
	}else{
		print ">read_$index\n$kmer\n" if($freq<=$max_occ);
		print ">read_$index\n$kmer_r\n" if($freq<=$max_occ);
	}
	$index++;
}
close $in_hdl;

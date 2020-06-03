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


die "perl $0 <infile>\n" if(@ARGV==0);

my($infile) = @ARGV;
my %hdl_hash;

my $in_hdl = myOpen($infile);
while (my $id = <$in_hdl>) {

	my $seq = <$in_hdl>;
	chomp $id;
	chomp $seq;

	my @array = split(/\s+/,$id);
	$id = "@array[0..$#array-1]";
	my $info = $array[-1];
	my($out_prefix,$len1,$len2) = split(/:/,$info);
	
	if(!defined $len2){
		my $out_file;
		$out_file = "$out_prefix.single.fasta.gz";
		if(!exists $hdl_hash{$out_file}){
			my $out_hdl = myOpen("|gzip>$out_file");
			$hdl_hash{$out_file} = $out_hdl;
		}
		my $out_hdl = $hdl_hash{$out_file};
		print $out_hdl "$id\n$seq\n";
	}

	if(defined $len2){
		my($read1,$read2);

		$read1 = "$out_prefix.1.fasta.gz";
		$read2 = "$out_prefix.2.fasta.gz";

		if(!exists $hdl_hash{$read1}){
			my $out_hdl1 = myOpen("|gzip>$read1");
			my $out_hdl2 = myOpen("|gzip>$read2");
			$hdl_hash{$read1} = $out_hdl1;
			$hdl_hash{$read2} = $out_hdl2;
		}

		my $out_hdl1 = $hdl_hash{$read1};
		my $out_hdl2 = $hdl_hash{$read2};

		my $seq1 = substr($seq,0,$len1);
		my $seq2 = substr($seq,-$len2,$len2);
		
		print $out_hdl1 "$id/1\n$seq1\n";
		print $out_hdl2 "$id/2\n$seq2\n";
	}

}
close $in_hdl;

foreach my $file(keys %hdl_hash)
{
	my $hdl = $hdl_hash{$file};
	close $hdl;
}
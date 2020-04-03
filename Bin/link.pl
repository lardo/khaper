use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path getcwd);

BEGIN {
    push (@INC,"$Bin");
}
use SeqIO;

die "perl $0 <contig.fa> <short.lst> <graph.info> <outfile> <simi:0>\n" if(@ARGV==0);

my($contig, $short_file, $graph, $outfile, $SIMI) = @ARGV;
$SIMI ||= 0;

# Record short contig
# ==================================================================================================
# |
my %short_hash;
my $in_hdl = myOpen($short_file);
while (my $id=<$in_hdl>)
{
	chomp $id;
	$short_hash{$id}=1;
}
close $in_hdl;
# |
# ==================================================================================================


# Load contigs
# ==================================================================================================
# |
my %seq_hash;

my($id, $seq);
$in_hdl = myOpen($contig);
while (getSeq($in_hdl, \$id, \$seq)!=-1) {
	next if(exists $short_hash{$id});
	$seq_hash{$id} = $seq;
}
close $in_hdl;
# |
# ==================================================================================================


# Output Scaffold
# ==================================================================================================
# |
my $index = 0;

open FO, ">$outfile";
open LOG,">$outfile.log";
open FO2,">$outfile.graph";
$in_hdl = myOpen($graph);

while (<$in_hdl>)
{
	chomp;

	$index++;
	my $id = "scaffold_$index";

	my $scaffold = "";
	my $ovl = 0;
	my $index2 = 0;
	my @links = split(/\s+/,$_);

	my @array = ();
	for(my $i=0; $i<=$#links; $i++)
	{
		my($contig, $strand, $len, $dist) = split(/,/,$links[$i]);

		my $seq = $seq_hash{$contig};
		$seq = revCom($seq) if($strand eq "-");
		if($scaffold eq "")
		{
			$scaffold = $seq;
			print LOG "$id|$index2\t$contig\t$strand\t0\n";
		}else{
			my ($merge, $simi) = &link($scaffold, $seq, $ovl);
			if($merge eq ""){
				print FO ">$id|$index2\n$scaffold\n";
				print FO2 ">$id|$index2\n@array\n";
				@array = ();

				$scaffold = $seq;
				$index2++;
				print LOG "$id|$index2\t$contig\t$strand\t0\n";
			}else{
				$scaffold = $merge;
				print LOG "$id|$index2\t$contig\t$simi\n";
			}
		}
		push(@array, $links[$i]);

		$ovl = $dist;
		delete $seq_hash{$contig};
	}
	print FO ">$id|$index2\n$scaffold\n";
	print FO2 ">$id|$index2\n@array\n";
}

while (my($id, $seq) = each %seq_hash) {
	$index++;

	my $id = "scaffold_$index";

	print FO ">$id|0\n$seq\n";
}

close $in_hdl;
close FO;
close FO2;
close LOG;
# |
# ==================================================================================================


# 进行连接
# ==================================================================================================
# |
sub link{
	my($seq1, $seq2, $ovl) = @_;
	my $len1 = length($seq1);
	my $len2 = length($seq2);

	my $part1 = substr($seq1, -$ovl, $ovl);
	my $part2 = substr($seq2, 0, 2*$ovl);

	my $jump = 1;
	$jump = 50 if($SIMI>0.1);
	my($simi, $anchor) = align($part1, $seq2, 13, $jump);


	return "" if($simi<$SIMI);
	return "" if($anchor eq "-1 -1");

	my($pos1, $pos2) = split(/\s+/, $anchor);

	my $merge = substr($seq1,0,$len1-$ovl+$pos1);
	$merge .= substr($seq2, $pos2);

	return ($merge, $simi);
}
# |
# ==================================================================================================

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

die "perl $0 <contig.fa> <graph.info> <outfile>\n" if(@ARGV==0);

my($contig, $graph, $outfile) = @ARGV;

# Load contigs
# ==================================================================================================
# |
my %seq_hash;

my($id, $seq);
my $in_hdl = myOpen($contig);
while (getSeq($in_hdl, \$id, \$seq)!=-1) {
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
$in_hdl = myOpen($graph);

while (<$in_hdl>)
{
	chomp;

	$index++;
	my $id = "scaffold_$index";
	
	my $scaffold;
	my @links = split(/\s+/,$_);
	foreach my $link(@links)
	{
		my($contig, $strand, $len, $dist) = split(/,/,$link);
		$dist = 101 if($dist<0);

		my $seq = $seq_hash{$contig};
		$seq = revCom($seq) if($strand eq "-");
		my $n_seq = "N"x$dist;
		$scaffold .= $seq;
		$scaffold .= $n_seq;

		delete $seq_hash{$contig};
	}
	print FO ">$id\n$scaffold\n";
}

while (my($id, $seq) = each %seq_hash) {
	$index++;

	my $id = "scaffold_$index";
	
	print FO ">$id\n$seq\n";
}

close $in_hdl;
close FO;
# |
# ==================================================================================================

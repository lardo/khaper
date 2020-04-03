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

die "perl $0 <graph.info> <outfile> <readfile1> <readfile2> ..<readfilen>\n" if(@ARGV==0);

my($graph, $outfile, @inputs) = @ARGV;


my $in_hdl;
my $ou_hdl;
my $ou2_hdl;

my %id_hash;
my %seq_hash;

# Record the graph's reads
# ==================================================================================================
# |
$in_hdl = myOpen($graph);
while (my $id = <$in_hdl>)
{
	my $nodes = <$in_hdl>;
	chomp $nodes;
	my @nodes = split(/\s+/, $nodes);
	foreach my $node(@nodes)
	{
		my($id, $strand) = split(/,/, $node);
		$id_hash{$id} = $strand;
	}
}
close $in_hdl;

foreach my $file(@inputs)
{
	my($id, $seq);
	$in_hdl = myOpen($file);
	while (getSeq($in_hdl,\$id,\$seq)!=-1) {
		next if(!exists $id_hash{$id});
		$seq_hash{$id} = $seq;
	}
	close $in_hdl;
}
# |
# ==================================================================================================


# Link the sequence for graph
# ==================================================================================================
# |
$in_hdl = myOpen($graph);
$ou_hdl = myOpen(">$outfile");
$ou2_hdl = myOpen(">$outfile.log");
while (my $id = <$in_hdl>)
{
	chomp $id;
	my $nodes = <$in_hdl>;
	chomp $nodes;
	my @nodes = split(/\s+/, $nodes);

	my $index = 0;
	my $pre_id;
	my $backbone = "";
	my @ctg_array = ();

	for(my $i=0; $i<@nodes; $i++)
	{
		my($ctg, $strand) = split(/,/, $nodes[$i]);
		my $seq = $seq_hash{$ctg};
		$seq = revCom($seq) if($strand eq "-");


		push(@ctg_array, $ctg);
		if($backbone eq "")
		{
			$backbone = $seq;
		}else{
			my($simi, $anchor, $time) = (0, "-1 -1", 0);
			while ($simi<0.1 && $time<4)
			{
				$time++;
				my $k = 17-2*$time;
				($simi, $anchor) = align($backbone, $seq, $k, 5);
			}

			if($simi<0.1)
			{
				debug("$pre_id vs $ctg: $simi");
				print $ou_hdl "$id|$index\n$backbone\n";
				print $ou2_hdl "$id|$index\n@ctg_array\n";

				@ctg_array = ();
				$pre_id  = $ctg;
				$backbone = "$seq";
				$index++;
				next;
			}

			my $len1 = length($backbone);
			my $len2 = length($seq);

			my($pos1, $pos2) = split(/\s+/, $anchor);
			my $tail1 = $len1-$pos1;
			my $tail2 = $len2-$pos2;

            # connect with different situation
            if($pos1<$pos2 && $tail1<$tail2){
                $backbone = $seq;
            }elsif($tail1>$tail2 && $i==$#nodes){
                $backbone = $backbone;
            }else{
                $backbone  = substr($backbone, 0, $pos1);
    			$backbone .= substr($seq, $pos2);
            }
		}
		$pre_id  = $ctg;
	}
	print $ou_hdl "$id|$index\n$backbone\n";
	print $ou2_hdl "$id|$index\n@ctg_array\n";
}
close $in_hdl;
close $ou_hdl;
close $ou2_hdl;
# |
# ==================================================================================================

use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path getcwd);
use List::Util qw(max min);

BEGIN {
    push (@INC,"$Bin");
}
use SeqIO;

die "perl $0 <graph.info> <node.pos> <outpefix>  \n" if(@ARGV==0);

my($ingraph, $nodefile, $outfile) = @ARGV;

my($in_hdl, $ou_hdl, $ou2_hdl);

# 记录backbone
# ==============================================================================
# |
debug("record backbone..");

my %bac_num_hash;
my %ctg_bac_hash;

my %ctg_idx_hash;
my %idx_ctg_hash;

$in_hdl = myOpen($ingraph);
while (my $line=<$in_hdl>) 
{
	# backbone
	chomp $line;
	$line =~ /^>(\S+)/;
	my $backbone = $1;

	# contigs
	$line = <$in_hdl>;
	chomp $line;
	my @nodes = split(/\s+/, $line);
	my $num = @nodes;

	# ctg2backbone
	for(my $i=0; $i<$num; $i++)
	{
		my $node = $nodes[$i];
		my($ctg, $strand, $pos) = split(/,/, $node);
		
		# 将contig转成数字，方便排序
		$ctg_idx_hash{"$ctg,$strand"} = $i;
		$idx_ctg_hash{"$backbone,$i"} = "$ctg,$strand";

		# contig在backbone上的情况
		$ctg_bac_hash{$ctg} = "$backbone $strand $pos";
	}

	$bac_num_hash{$backbone} = $num;
}
close $in_hdl;
# |
# ==============================================================================






# 对pacbio reads 进行归类
# ==============================================================================
# |
debug("class pacbio reads..");

$in_hdl = myOpen($nodefile);
$ou_hdl = myOpen(">$outfile.pb");
$ou2_hdl= myOpen(">$outfile.mb");
while (my $line=<$in_hdl>) 
{
	# pacbio
	chomp $line;
	$line =~ /^>(\S+)/;
	my $id = $1;
	
	# contigs
	$line = <$in_hdl>;
	chomp $line;
	my @nodes = split(/\s+/, $line);
	
	# target backbone
	my @array = targetBac(@nodes);
	my($backbone, $strand, $num) = split(/\s+/, $array[0]);

	# output the pacbio read #
	my $head = "$id,$strand\t$backbone\t";
	@nodes = revPB(@nodes) if($strand eq "-");

	my $info = "";
	my @matrix = ();
	foreach my $node (@nodes){
		my($ctg, $strand, $pos) = split(/,/, $node);
		next if(!exists $ctg_bac_hash{$ctg});
		
		my($strand2, $backbone2, $pos2, $index) = getCtgInfo($ctg);
		next if($backbone2 ne $backbone);
		next if($strand ne $strand2);

		if(isMatch($pos, $pos2, \@matrix)!=0){
			$info .= "$index\t";
		}
	}
	$num = scalar @matrix;
	print $ou_hdl "$head\t$num\t$info\n";

	# output the maybe connection
	if(@array>1){
		foreach(@array)
		{
			my($backbone, $strand, $num) = split;
			print $ou2_hdl "$backbone,$strand\t";
		}
		print $ou2_hdl "\n";
	}
}
close $in_hdl;
close $ou_hdl;
close $ou2_hdl;

`sort -k 2,2 -k 4,4n -k 3,3nr -S 1g -o $outfile.pb $outfile.pb`;
# |
# ==============================================================================


# 对backbone进行连接
# ==============================================================================
# |
my %idx_pb_hash;

my $pre_bac = "";
$in_hdl = myOpen("$outfile.pb");
$ou_hdl = myOpen(">$outfile.lk");
while (<$in_hdl>) {
	chomp;
	my($pb, $backbone, $num, @indexs) = split;
	outBackbone() if($backbone ne $pre_bac);
	for(my $i=0; $i<$num; $i++)
	{
		my $idx = $indexs[$i];
		$idx_pb_hash{$idx}{$pb} = $indexs[-1];
	}
	$pre_bac = $backbone;
}
outBackbone();
close $in_hdl;
close $ou_hdl;

sub outBackbone
{
	return if($pre_bac eq "");

	my $num = $bac_num_hash{$pre_bac};
	print $ou_hdl ">$pre_bac\n";
	for(my $idx=0; $idx<$num; $idx++)
	{
		my($max, $best_pb) = (-1, "");
		while (my($pb, $idx2) = each %{$idx_pb_hash{$idx}}) 
		{
			next if($idx2<$max);
			$max = $idx2;
			$best_pb = $pb;
		}
		my $ctg = $idx_ctg_hash{"$pre_bac,$idx"};
		print $ou_hdl "$ctg\t";
		
		if($max > $idx){
			print $ou_hdl "$best_pb\t";
			$idx = $max-1;
		}
	}
	print $ou_hdl "\n";
	%idx_pb_hash = ();
}
# |
# ==============================================================================




# 节点能否与backbone一致
# ==============================================================================
# |
sub isMatch
{
	my($x, $y, $arr_ref) = @_;
	if(@{$arr_ref}==0)
	{
		push(@{$arr_ref}, "$x,$y");
		return 1;
	}
	# 斜率
	my($x2, $y2) = split(/,/, ${$arr_ref}[-1]);
	return 0 if($y2==$y);
	
	my $slope = ($x2-$x)/($y2-$y);
	return 0 if($slope<0.9 or $slope>1.1);

	push(@{$arr_ref}, "$x,$y");
	return 1;
}
# |
# ==============================================================================



# 判断pb所属的backbone
# ==============================================================================
# |
sub targetBac
{
	my @nodes = @_;

	my $num = @nodes;
	
	# ctg 2 backbone #
	my %bac_hash;
	for(my $i=0; $i<$num; $i++)
	{
		my $node = $nodes[$i];
		my($ctg, $strand, $pos) = split(/,/, $node);

		my($strand2, $backbone, $pos2, $index) = getCtgInfo($ctg);
		next if($index == -1);

		my $strand_ = "+";
		$strand_ = "-" if($strand ne $strand2);
		$bac_hash{"$backbone $strand_"}++;
	}

	# get the backbone #
	my @array;
	while (my($backbone, $num) = each %bac_hash) 
	{
		push(@array, "$backbone $num");
	}
	@array = sort by_num(@array);
}


sub getCtgInfo
{
	my $ctg = shift;
	my($strand, $backbone, $pos, $index) = ("", "", -1, -1);
	if(exists $ctg_bac_hash{$ctg}){
		($backbone, $strand, $pos) = split(/\s+/, $ctg_bac_hash{$ctg});
		$index = $ctg_idx_hash{"$ctg,$strand"};
	}
	return($strand, $backbone, $pos, $index);
}


sub by_num
{
	my($id1, $strand1, $num1) = split(/\s+/, $a);
	my($id2, $strand2, $num2) = split(/\s+/, $b);
	return $num2<=>$num1;
}
# |
# ==============================================================================




# 反转节点
# ==============================================================================
# |
sub revPB
{
	my @array = @_;
	@array = reverse @array;
	my($id_, $strand_, $pos_) = split(/,/, $array[0]);
	for(my $i=0; $i<=$#array; $i++)
	{
		my($id, $strand, $pos) = split(/,/, $array[$i]);
		$strand = ($strand eq "+") ? "-":"+";
		$pos = $pos_-$pos;
		$array[$i] = "$id,$strand,$pos";
	}
	return @array;
}
# |
# ==============================================================================


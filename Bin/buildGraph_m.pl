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

die "perl $0 <node.pos>  <out.graph> <min_occ:2>\n" if(@ARGV==0);

my($infile, $outfile, $min_occ) = @ARGV;
$min_occ ||= 2;

my($in_hdl, $ou_hdl);

# 记录连接关系
# ==============================================================================
# |
my %nodes; # 记录节点，方便遍历
my(%left_hash, %right_hash);

$in_hdl = myOpen($infile);
while (<$in_hdl>) 
{
	chomp;
	next if(/^>/);
	my @nodes = split;
	my $num = @nodes;

	# 记录两两连接关系
	for(my $i=0; $i<$num; $i++)
	{
		# 左节点
		my $left_node = $nodes[$i];
		my($id1, $strand1, $pos1) = split(/,/, $left_node);

		$left_node = "$id1,$strand1";
		$nodes{$id1} = 1;

		# 右节点
		for(my $j=$i+1; $j<$num; $j++)
		{
			my $right_node = $nodes[$j];
			my($id2, $strand2, $pos2) = split(/,/, $right_node);
			
			next if($id1 eq $id2);

			my $dist = $pos2-$pos1;
			next if($dist<300);
			$right_node = "$id2,$strand2";
			$right_hash{$left_node}{$right_node} .= "$dist\t";
			$left_hash{$right_node}{$left_node}  .= "$dist\t";

			my $left_node_r  = revLink($right_node);
			my $right_node_r = revLink($left_node);

			$left_hash{$right_node_r}{$left_node_r}  .= "$dist\t";
			$right_hash{$left_node_r}{$right_node_r} .= "$dist\t";
		}
	}
}
close $in_hdl;
# |
# ==============================================================================



# 进行连接
# ==============================================================================
# |
my @scaffold;

my $index = 0;
my %parse_hash;	# 存储遍历的数据
open FO, ">$outfile" or die $!;
foreach my $ctg (keys %nodes)
{
	@scaffold = ();
	next if(exists $parse_hash{$ctg});

	# 向右连接
	my $next = "$ctg,+";
	while ($next ne "") {		
		my($contig, $strand) = split(/,/, $next);

		$parse_hash{$contig} = 1;
		push(@scaffold, $next);
		
		$next = getNext($next, \%right_hash);
	}

	# 向左连接
	@scaffold = revScf(@scaffold);
	$next = $scaffold[-1];
	$next = getNext($next, \%right_hash);
	while ($next ne "") {
		my($contig, $strand) = split(/,/, $next);

		$parse_hash{$contig} = 1;
		push(@scaffold, $next);
		
		$next = getNext($next, \%right_hash);
	}

	# 输出框架
	$index++;
	print FO ">backbone_$index\n";

	my $pre_dist = 0;
	my $num = scalar(@scaffold);
	for(my $i=0; $i<$num; $i++){
		my $link1 = $scaffold[$i];
		my($dist_, $num_) = (0, 0);
		if($i>0){
			my $link2 = $scaffold[$i-1];
			($num_, $dist_) = getLinkInfo($link2, $link1, \%right_hash);
			$pre_dist += $dist_;
		}
		print FO "$link1,$pre_dist\t";
	}
	print FO "\n";
}
close FO;

sub revScf
{
	my @scf = @_;
	@scf = reverse(@scf);
	for(my $i=0; $i<@scf; $i++)
	{
		my($id, $strand) = split(/,/, $scf[$i]);
		my $strand_r = ($strand eq "+")? "-":"+";
		$scf[$i] = "$id,$strand_r";
	}
	return @scf;
}
# |
# ==============================================================================


# 取得下一个节点
# ==============================================================================
# |
sub getNext{
	my($link) = @_;
	my $next_link = "";

	my @link2s = keys %{$right_hash{$link}};
	my @array = sortLinks($link, \@link2s, \%right_hash, $min_occ);
	
	return "" if(@array==0);

	# 判定这个节点是否被支持
	my $end_link = (split(/\s+/,$array[-1]))[0];
	for(my $i=0; $i<@array; $i++)
	{
		my($link2, $num, $dist, $error_num) = split(/\s+/, $array[$i]);
		my($contig, $strand) = split(/,/, $link2);

		# next if($error_num>2);
		next if(exists $parse_hash{$contig});
		if($link2 ne $end_link){
			next if(!exists $right_hash{$link2}{$end_link});
		}

		if($link2 eq $end_link){
			next if($dist<2000);
		}

		# 左边最远距离结点
		my @link3s = keys %{$left_hash{$link2}};
		my @array3 = sortLinks($link2, \@link3s, \%left_hash, $min_occ);
		next if(@array3==0);
		my $link3 = (split(/\s+/,$array3[-1]))[0];
		if($link3 ne $link){
			next if(!exists $right_hash{$link3}{$link});
		}

		$next_link = $link2;
		last;
	}
	return $next_link;
}


sub sortLinks{
	my($start, $links_arr, $hash_ref, $min_occ) = @_;
	$min_occ ||= 2;

	my @array;
	foreach my $link(@{$links_arr})
	{
		my($num, $dist, $error_num) = getLinkInfo($start, $link, $hash_ref);
		next if($num<$min_occ);
		push(@array, "$link $num $dist $error_num");
	}
	@array = sort by_link_dist @array;	
	return @array;
}

sub getLinkInfo
{
	my($link1, $link2, $hash_ref) = @_;
	my $info = ${$hash_ref}{$link1}{$link2};
	
	my @dist_arr = split(/\s+/, $info);
	@dist_arr = sort @dist_arr;
	
	my $num = scalar @dist_arr;
	my $dist = $dist_arr[$num/2];

	my $error_num = 0;
	foreach my $dist1(@dist_arr){
		my $min = min($dist1, $dist);
		next if($min==0);

		my $distant = abs($dist-$dist1);
		my $perc = $distant/$min;
		$error_num++ if($perc>0.15 && $distant>100);
	}
	return($num, $dist, $error_num);
}

sub by_link_dist{
	my($link1, $num1, $dist1, $error_num1) = split(/\s+/, $a);
	my($link2, $num2, $dist2, $error_num2) = split(/\s+/, $b);
	return $dist1<=>$dist2;
}
# |
# ==============================================================================



# 反转节点
# ==============================================================================
# |
sub revLink{
	my $link = shift;
	my($contig, $strand) = split(/,/, $link);
	if($strand eq "+"){
		$strand = "-";
	}else{
		$strand = "+";
	}
	return "$contig,$strand";
}
# |
# ==============================================================================


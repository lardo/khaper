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

die "perl $0 <filt.al> <out.graph>\n" if(@ARGV==0);

my($infile, $outfile) = @ARGV;

# 记录被包含的序列
# ==============================================================================
# |
my %short_hash;

my $in_hdl = myOpen($infile);
while (<$in_hdl>) {
	chomp;
	my($q_id, $t_id, $q_len, $t_len,
		$strand, $ovl, $q_pos, $t_pos, $score) = split;

	if($q_len==$t_len && $q_len<=$ovl)
	{
		$short_hash{$q_id} = 1 if($q_id lt $t_id);
		$short_hash{$t_id} = 1 if($t_id lt $q_id);
		next;
	}

	$short_hash{$q_id} = 1 if($q_len<=$ovl);
	$short_hash{$t_id} = 1 if($t_len<=$ovl);
}
close $in_hdl;
# |
# ==============================================================================


# 记录连接关系
# ==============================================================================
# |
my %nodes; # 记录节点，方便遍历
my(%left_hash, %right_hash);

$in_hdl = myOpen($infile);
while (<$in_hdl>)
{
	chomp;
	my($q_id, $t_id, $q_len, $t_len,
		$strand, $ovl, $q_pos, $t_pos, $score) = split;

	next if(exists $short_hash{$q_id});
	next if(exists $short_hash{$t_id});

	# record the node #
	my $left_node  = "$q_id,$strand,$q_len";
	my $right_node = "$t_id,+,$t_len";
	if($q_pos<$t_pos)
	{
		$left_node  = "$t_id,+,$t_len";
		$right_node = "$q_id,$strand,$q_len";
	}
	$right_hash{$left_node}{$right_node} = "$ovl $score";
	$left_hash{$right_node}{$left_node}  = "$ovl $score";

	# reverse the node #
	my $left_node_r  = revLink($right_node);
	my $right_node_r = revLink($left_node);

	$left_hash{$right_node_r}{$left_node_r}  = "$ovl $score";
	$right_hash{$left_node_r}{$right_node_r} = "$ovl $score";

	# record nodes #
	$nodes{$left_node}  = 1;
	$nodes{$right_node} = 1;
	$nodes{$left_node_r}  = 1;
	$nodes{$right_node_r} = 1;
}
close $in_hdl;
# |
# ==============================================================================


# 进行连接
# ==============================================================================
# |
my %parse_hash;	# 存储遍历的数据
open FO, ">$outfile" or die $!;
foreach my $link (keys %nodes)
{
	my @scaffold;

	my($contig, $strand, $len) = split(/,/, $link);
	next if(exists $parse_hash{$contig});

	my $next = $link;
	# 向右连接
	while ($next ne "") {
		my($contig, $strand, $len) = split(/,/, $next);
		push(@scaffold, $next);
		$parse_hash{$contig} = 1;

		$next = getNext($next, \%right_hash);
	}

	# 向左连接
	$next = getNext($link, \%left_hash);
	while ($next ne "") {
		my($contig, $strand, $len) = split(/,/, $next);
		@scaffold = ($next, @scaffold);
		$parse_hash{$contig} = 1;

		$next = getNext($next, \%left_hash);
	}

	# 输出
	for(my $i=0; $i<=$#scaffold; $i++)
	{
		my $link1 = $scaffold[$i];
		my $ovl = 0;
		if($i<$#scaffold)
		{
			my $link2 = $scaffold[$i+1];
			my($ovl_, $score, $extend) = getLinkInfo($link1, $link2, \%right_hash);
			$ovl = $ovl_;

			# 去掉被包含的序列 #
			my @others = keys %{$right_hash{$link1}};
			foreach my $other(@others)
			{
				my($id, $strand, $len) = split(/,/, $other);
				next if(exists $parse_hash{$id});

				if(exists $left_hash{$other}{$link2} or
						exists $right_hash{$other}{$link2})
				{
					$parse_hash{$id} =1;
					$short_hash{$id} =1;
				}
			}
		}
		print FO "$link1,$ovl\t";
	}
	print FO "\n";
}
close FO;
# |
# ==============================================================================

# 输出要过滤的序列ID
# ==============================================================================
# |
open FO,">short.lst";
foreach my $id(keys %short_hash)
{
	print FO "$id\n";
}
close FO;
# |
# ==============================================================================



# 取得两个节点的连接长度、分数、连完后的长度
# ==============================================================================
# |
sub getLinkInfo
{
	my($link1, $link2, $hash_ref) = @_;
	my $info = ${$hash_ref}{$link1}{$link2};

	my($id1, $strand1, $len1) = split(/,/, $link1);
	my($id2, $strand2, $len2) = split(/,/, $link2);

	my($ovl, $score) = split(/\s+/, $info);
	my $extend = $len1+$len2-$ovl;

	return($ovl, $score, $extend);
}
# |
# ==============================================================================


# 按连接的分数排序
# ==============================================================================
# |
sub by_link_score{
	my($link1, $ovl1, $score1, $extend1) = split(/\s+/, $a);
	my($link2, $ovl2, $score2, $extend2) = split(/\s+/, $b);

	return $score2<=>$score1;
}

sub by_link_ovl{
	my($link1, $ovl1, $score1, $extend1) = split(/\s+/, $a);
	my($link2, $ovl2, $score2, $extend2) = split(/\s+/, $b);

	return $ovl2<=>$ovl1;
}
# |
# ==============================================================================



# 对节点进行排序
# ==============================================================================
# |
sub sortLinks{
	my($start, $links_arr, $hash_ref) = @_;
	my @array;
	foreach my $link(@{$links_arr})
	{
		my($ovl, $score, $extend) = getLinkInfo($start, $link, $hash_ref);
		push(@array, "$link $ovl $score $extend");
	}
	# 从大肠肝菌的数据来看
	# 用overlap来连接和判断准确度最高
	@array = sort by_link_ovl @array;
	return @array;
}
# |
# ==============================================================================



# 取得下一个节点
# ==============================================================================
# |
sub getNext{
	my($link, $hash_ref) = @_;
	my $next_link = "";

	my @link2s = keys %{${$hash_ref}{$link}};
	my @array = sortLinks($link, \@link2s, $hash_ref);

	foreach my $next (@array){
		my($link2, $ovl, $score, $extend) = split(/\s+/, $next);
		my($contig, $strand, $len) = split(/,/, $link2);
		next if(exists $parse_hash{$contig});

		if(isBest($link, $link2))
		{
			$next_link = $link2;
			last;
		}
	}
	return $next_link;
}
# |
# ==============================================================================


# 是否是end_link的最好节点
# ==============================================================================
# |
sub isBest{
	my($start_link, $end_link) = @_;

	my $hash_ref = \%right_hash;
	$hash_ref = \%left_hash if(exists $left_hash{$end_link}{$start_link});

	my @link2s = keys %{${$hash_ref}{$end_link}};
	return 1 if(@link2s==1);

	my @array = sortLinks($end_link, \@link2s, $hash_ref);
	my($link2, $ovl, $score, $extend) = split(/\s+/, $array[0]);

	return 1 if($link2 eq $start_link);
	return 1 if(exists $right_hash{$link2}{$start_link});
	return 1 if(exists $left_hash{$link2}{$start_link});

	return 0;
}
# |
# ==============================================================================


# 反转节点
# ==============================================================================
# |
sub revLink{
	my $link = shift;
	my($contig, $strand, $len) = split(/,/, $link);
	if($strand eq "+"){
		$strand = "-";
	}else{
		$strand = "+";
	}
	return "$contig,$strand,$len";
}
# |
# ==============================================================================

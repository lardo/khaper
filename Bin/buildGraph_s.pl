use strict;
use warnings;

die "perl $0 <link.info> <out.graph> <min_conn:2> <min_dist:10> <strategy:1>\n" if(@ARGV==0);

my($infile, $outfile, $min_conn, $min_dist, $strategy) = @ARGV;
$min_conn ||= 2;
$min_dist ||= 10;
$strategy ||= 1;

# 记录连接关系
# ==============================================================================
# |
my %nodes; # 记录节点，方便遍历

my(%left_hash, %right_hash);
open FI, $infile or die $!;
while (<FI>) {
	my($pre_pb, $link1, $link2, $dist) = split;

    next if($dist<$min_dist);

	$right_hash{$link1}{$link2} .= "$dist ";
	$left_hash{$link2}{$link1} .= "$dist ";

	$nodes{$link1} = 1;
	$nodes{$link2} = 1;

	# record the other strand chrome
	$link1 = revLink($link1);
	$link2 = revLink($link2);

	$right_hash{$link2}{$link1} .= "$dist ";
	$left_hash{$link1}{$link2} .= "$dist ";

	$nodes{$link1} = 1;
	$nodes{$link2} = 1;
}
close FI;

# 清除不符合条件的连接
cleanLink(\%left_hash);
cleanLink(\%right_hash);
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
		my $distant = 0;
		if($i<$#scaffold)
		{
			my $link2 = $scaffold[$i+1];
			my($dist, $num) = getLinkInfo($link1, $link2, \%right_hash);
			$distant = $dist;
		}
		print FO "$link1,$distant\t";
	}
	print FO "\n";
}
close FO;
# |
# ==============================================================================


# 按连接的条数排序
# ==============================================================================
# |
sub by_link_num{
	my($link1, $num1, $dist1) = split(/\s+/, $a);
	my($link2, $num2, $dist2) = split(/\s+/, $b);

	return $num2<=>$num1 if($num1 != $num2);
	return $dist1<=>$dist2;
}

sub by_link_dist{
	my($link1, $num1, $dist1) = split(/\s+/, $a);
	my($link2, $num2, $dist2) = split(/\s+/, $b);

	return $dist1<=>$dist2;
}
# |
# ==============================================================================


# 取得两个节点的连接距离和支持数
# ==============================================================================
# |
sub getLinkInfo
{
	my($link1, $link2, $hash_ref) = @_;
	my $dist_info = ${$hash_ref}{$link1}{$link2};
	my @dists = split(/\s+/,$dist_info);
	@dists = sort {return $a<=>$b} @dists;
	my $num = scalar @dists;
	my $mediant = $dists[$num/2];

	return($mediant, $num);
}
# |
# ==============================================================================


# 取得下一个节点
# ==============================================================================
# |
sub getNexts{
	my($link, $hash_ref, $sort) = @_;
	# sort
	# 1: by link num
	# 2: by link distant
	$sort ||= 1;

	my @array;
	my @link2s = keys %{${$hash_ref}{$link}};
	foreach my $link2(@link2s)
	{
		my($mediant, $num) = getLinkInfo($link, $link2, $hash_ref);
		push(@array, "$link2 $num $mediant");
	}
	@array = sort by_link_num @array if($sort == 1);
	@array = sort by_link_num @array if($sort == 2);
	return @array;
}

sub getNext{
	my($link, $hash_ref) = @_;
	my $next_link = "";

	my @array = getNexts($link, $hash_ref);

	foreach my $next (@array){
		my($link2, $num, $dist) = split(/\s+/, $next);
		my($contig, $strand, $len) = split(/,/, $link2);

		next if(exists $parse_hash{$contig});
		next if(!isBest($link, $link2));

		$next_link = $link2;
		last;
	}
	return $next_link;
}
# |
# ==============================================================================

# 判断是否能连接
# ==============================================================================
# |
sub isBest{
	my($start_link, $end_link) = @_;

    return 1 if($strategy !=1);

	my $hash_ref = \%right_hash;
	$hash_ref = \%left_hash if(exists $left_hash{$end_link}{$start_link});

	my @array = getNexts($end_link, $hash_ref, 2);
	my($link2, $num, $dist) = split(/\s+/, $array[-1]);

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


# 清除不满足条件的连接
# ==============================================================================
# |
sub cleanLink{
	my $hash_ref = shift;
	foreach my $link1(%{$hash_ref})
	{
		my @link2s = keys %{${$hash_ref}{$link1}};
		foreach my $link2(@link2s)
		{
			my $dist_info = ${$hash_ref}{$link1}{$link2};
			my @dists = split(/\s+/,$dist_info);
			if(@dists<$min_conn){
				delete ${$hash_ref}{$link1}{$link2};
			}
		}
	}
}
# |
# ==============================================================================

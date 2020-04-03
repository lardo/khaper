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

die "perl $0 <filt.al> <outfile>\n" if(@ARGV==0);

my($infile, $outfile) = @ARGV;

# 将read和contig转化成编号，节省内存
# ==============================================================================
# |
my @pbc_arr;
my @ctg_arr;

my %pbc2idx_hash;
my %ctg2idx_hash;

my $ctg_len = 0;
my $pbc_len = 0;
my $in_hdl = myOpen($infile);
while (<$in_hdl>) 
{
	chomp;
	my($q_id, $t_id, $q_len, $t_len,
		$strand, $ovl, $q_pos, $t_pos, $score) = split;

	if(!exists $pbc2idx_hash{$q_id}){
		push(@pbc_arr, $q_id);
		$pbc2idx_hash{$q_id} = $#pbc_arr;
		$pbc_len += $q_len;
	}

	if(!exists $ctg2idx_hash{$t_id}){
		push(@ctg_arr, $t_id);
		$ctg2idx_hash{$t_id} = $#ctg_arr;
		$ctg_len += $t_len;
	}
}
close $in_hdl;
debug("total pacbio len: $pbc_len");
debug("total contig len: $ctg_len");
# |
# ==============================================================================



# 记录contig对应的pacbio
# 记录pacbio对应的contig
# ==============================================================================
# |

my %ctg_pbc_hash;
my %pbc_ctg_hash;
$in_hdl = myOpen($infile);
while (<$in_hdl>) 
{
	chomp;
	my($q_id, $t_id, $q_len, $t_len,
		$strand, $ovl, $q_pos, $t_pos, $score) = split;

	my $pbc_idx = $pbc2idx_hash{$q_id};
	my $ctg_idx = $ctg2idx_hash{$t_id};

	$ctg_pbc_hash{$ctg_idx}{$pbc_idx} = $q_len;
	$pbc_ctg_hash{$pbc_idx}{$ctg_idx}++;
}
close $in_hdl;
# |
# ==============================================================================



# 输出每个tag对应的pacbio
# ==============================================================================
# |
debug("Class reads...");
my %ctg_hash;
my $max_pbc_num=20;
my $min_ctg_num=10;

my $ou_hdl = myOpen(">$outfile");
foreach my $ctg_idx (keys %ctg_pbc_hash)
{
	next if(exists $ctg_hash{$ctg_idx});
	$ctg_hash{$ctg_idx} = 1;

	my @pacbio_arr = keys %{$ctg_pbc_hash{$ctg_idx}};
	
	# get the related pacbio
	my @array;
	foreach my $pbc_idx (@pacbio_arr)
	{
		my $len =  $ctg_pbc_hash{$ctg_idx}{$pbc_idx};
		push(@array, "$pbc_idx $len");
	}
	@array = sort by_len(@array);

	# get related contig
	my $pb_num = 0;
	my %temp_hash;
	foreach(@array)
	{
		$pb_num++;
		my($pbc_idx, $len) = split;
		print $ou_hdl "$ctg_arr[$ctg_idx]\t$pbc_arr[$pbc_idx]\t$len\n";

		my @ctg_idxs = keys %{$pbc_ctg_hash{$pbc_idx}};
		foreach(@ctg_idxs){
			$temp_hash{$_}++;
		}
		last if($pb_num>$max_pbc_num);
	}

	while (my($ctg_idx, $num) = each %temp_hash) 
	{
		next if(exists $ctg_hash{$ctg_idx});
		next if($num<$min_ctg_num);
		$ctg_hash{$ctg_idx} = $num;
		# debug("log $ctg_idx");
	}
}
close $ou_hdl;
# |
# ==============================================================================



# 排序
# ==============================================================================
# |
sub by_len
{
	my($id1, $len1) = split(/\s+/, $a);
	my($id2, $len2) = split(/\s+/, $b);

	return $len2<=>$len1;
}
# |
# ==============================================================================


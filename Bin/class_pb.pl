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

die "perl $0 <class.tab> <pacbio.fa> <outdir:Split>\n" if(@ARGV==0);

my($classfile, $readfile, $outdir) = @ARGV;
$outdir ||= "Split";

my $in_hdl;
my $ou_hdl;
my $outfile = "temp.fa";
# 对pacbio的read进行聚类
# ==============================================================================
# |
debug("Loading class file...");
my %pb_tag_hash;

my $index = 0;
my $pre_tag = "";
my @array = ();
$in_hdl = myOpen("$classfile");
while (<$in_hdl>) 
{
	chomp;
	my($tag, $pb, $len) = split;
	record() if($pre_tag ne $tag);
	push(@array, $pb);
	$pre_tag = $tag;
}
close $in_hdl;
record();
# |
# ==============================================================================



# 对pacbio数据进行格式转换
# ==============================================================================
# |
debug("Loading pacbio file...");

my($id, $seq);
$in_hdl = myOpen($readfile);
$ou_hdl = myOpen(">$outfile");
while (getSeq($in_hdl, \$id, \$seq)!=-1) 
{
	next if(!exists $pb_tag_hash{$id});
	my $tags = 	$pb_tag_hash{$id};
	my @tag_arr = split(/\s+/, $tags);

	my $len = length($seq);
	foreach my $tag(@tag_arr)
	{
		print $ou_hdl "$tag\t$len\t$id\t$seq\n";
	}
}
close $in_hdl;
close $ou_hdl;


`sort -k 1,1n -k 2,2nr -S 5G -o $outfile $outfile`;
# |
# ==============================================================================



# 对每个文件进行拆分
# ==============================================================================
# |
debug("Split File...");
$outfile = abs_path($outfile);

mkdir($outdir);
chdir($outdir);

$index = -1;
@array = ();
$pre_tag = -1;
my $base  = 0;

my $max_base = 30000000;
$in_hdl = myOpen($outfile);
while (<$in_hdl>) 
{
	chomp;
	my($tag, $len, $id, $seq) = split;
	output() if($pre_tag!=$tag && $base>$max_base);
	push(@array, "$tag\t$id\t$seq");

	$pre_tag = $tag;
	$base += length($seq);
}
close $in_hdl;
# |
# ==============================================================================
`rm $outfile`;


# 
# ==============================================================================
# |
sub output
{
	return if(@array==0);
	my $num_dir = 50;

	$index++;
	my $dir_idx = int($index/$num_dir);
	my $sub_dir = "sub_$dir_idx";
	mkdir($sub_dir) if(!-e $sub_dir);
	
	my %id_hash;
	my $tag0 = "";
	my $out_hdl1 = myOpen(">$sub_dir/split_$index.fasta");
	my $out_hdl2 = myOpen(">$sub_dir/split_$index.subread.fasta");
	foreach(@array)
	{
		my($tag, $id, $seq) = split;
		next if(exists $id_hash{$id});
		$id_hash{$id} = 1;

		print $out_hdl1 ">$id\n$seq\n" if($tag ne $tag0);
		print $out_hdl2 ">$id\n$seq\n";
		$tag0 = $tag;
	}
	close $out_hdl1;
	close $out_hdl2;

	$base  = 0;
	@array = ();
}
# |
# ==============================================================================


# 
# ==============================================================================
# |
sub record{
	return if(@array==0);
	
	$index++;
	foreach(@array)
	{
		$pb_tag_hash{$_} .= "$index ";
	}	
	@array = ();
}
# |
# ==============================================================================

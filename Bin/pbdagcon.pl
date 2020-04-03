use strict;
use warnings;
use FindBin qw($Bin);
use File::Basename;
use Cwd qw(abs_path getcwd);
BEGIN {
    push (@INC,"$Bin");
}
use SeqIO;

die "perl $0 <backbone.fasta> <subread.fasta> <thread:4> <min_occ:1>\n" if(@ARGV==0);

# 读入参数
# ==============================================================================
# |
my($infile, $subread, $thread, $min_occ) = @ARGV;
$infile = abs_path($infile);
$subread = abs_path($subread);

$thread ||= 4;
$min_occ ||= 1;

# 如果已存在被纠错，则退出
exit(1) if(-e "$infile.consensus");

my $dir_name = basename($infile).".dir";
mkdir "$dir_name";
chdir($dir_name);
# |
# ==============================================================================




# 进行纠错
# ==============================================================================
# |
cutFiles("subread.fasta", $subread, $infile);

my $input = $infile;
my $outfile;
for(my $i=0; $i<1; $i++)
{
	my $cmd = "source $Bin/source.sh\n";
	$cmd .= "$Bin/blasr subread.fasta $input -out mapped.m5 ";
	$cmd .= "-fastSDP -minMatch 14 -bestn 1 -m 5  -nproc $thread\n";
	$cmd .= "sort -k 6,6 -k 8,8n mapped.m5 -o mapped.m5 -S 1G -T ./ \n";
	$cmd .= "pbdagcon -c $min_occ -j $thread mapped.m5 > consensus.$i.fasta 2>log\n";

	system($cmd);
	$outfile = "consensus.$i.fasta";
	$input  = $outfile;
}

# |
# ==============================================================================



# 整理结果
# ==============================================================================
# |
my $cmd  = "cd ..\n";
$cmd .= "mv $dir_name/$outfile $infile.consensus\n";
$cmd .= "rm -r $dir_name\n";

system($cmd);
# |
# ==============================================================================



# 对数据进行剪切
# ==============================================================================
# |
sub cutFiles{
	my $MAX_LEN = 50000;
    my $OVL = 1000;

	my ($outfile, @files) = @_;
	my $ou_hdl = myOpen(">$outfile");
	foreach my $file(@files){
		my($id, $seq);
		my $in_hdl = myOpen($file);
		while (getSeq($in_hdl, \$id, \$seq)!=-1)
		{
			my $len = length($seq);
			for(my $i=0; $i<$len; $i+=$MAX_LEN-$OVL)
			{
				my $sub_seq = substr($seq, $i, $MAX_LEN);
				print $ou_hdl ">$id|$i\n$sub_seq\n";
			}
		}
	}
	close $ou_hdl;
}
# |
# ==============================================================================

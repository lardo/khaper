use strict;
use warnings;
use FindBin qw($Bin);
use File::Basename;
use Cwd qw(abs_path getcwd);
BEGIN {
    push (@INC,"$Bin");
}
use SeqIO;

die "perl $0 <backbone.fasta> <contig.fasta> <thread:4> \n" if(@ARGV==0);

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
my $cmd = "source $Bin/source.sh\n";
$cmd .= "proovread -l $infile -u $subread --overwrite -t $thread"; #-m dazz-utg-noccs";
system($cmd);
# |
# ==============================================================================



# 整理结果
# ==============================================================================
# |
my $outfile = "correct.fa";
$cmd  = "perl $Bin/fq2fa.pl proovread/proovread.untrimmed.fq $outfile\n";
$cmd .= "cd ..\n";
$cmd .= "mv $dir_name/$outfile $infile.consensus\n";
$cmd .= "rm -r $dir_name\n";

system($cmd);
# |
# ==============================================================================

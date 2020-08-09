use strict;
use warnings;
use FindBin qw($Bin);
use File::Basename;
use Getopt::Long;
use Cwd qw(abs_path getcwd);

BEGIN {
    push (@INC,"$Bin");
}
use SeqIO;
use Qsub;

=head1 Usage

 perl remDup.pl <genome.fa> <outdir> <cutoff:0.7> 

 Options:
        --ref   <str> The ref genome to build kbit
        --kbit  <str> The unique kmer file
        --kmer  <int> the kmer size [15]
        --sort  <int> sort seq by length [1]

=head1 Description

 This script is to remove dupplcation seq


=cut

die `pod2text $0` if (@ARGV ==0);


# Get parameters
# ==============================================================================
#
my($kbit, $kmer_size, $ref, $sortseq) = ("-", 15, "-", 1);

GetOptions(
    "ref:s" => \$ref,
    "kbit:s" => \$kbit,
    "kmer:s" => \$kmer_size,
    "sort:i" => \$sortseq,
);
$kbit = abs_path($kbit);
$ref = abs_path($ref) if(-e $ref);

my($genome, $outdir, $cutoff) = @ARGV;

$genome = abs_path($genome);
$outdir = abs_path($outdir);
$cutoff||= 0.7;

mkdir($outdir);
chdir($outdir);
#
# ==============================================================================





#  
# ==============================================================================
# 
my $cmd = "";
my $BIN = "$Bin/";

$cmd .= "perl $BIN/sortSeq.pl $genome sort.fa 5 \n" if($sortseq==1);
$cmd .= "cp $genome sort.fa \n" if($sortseq!=1);

# 构建kmer表
if(!-e $kbit && !-e $ref)
{
    outShell2("sort.fa", "kmer.lst");
    
    $cmd .= "perl $BIN/Graph.pl pipe -i kmer.lst ";
    $cmd .= " -m 1 -k 15 -s 1,5 -d Kmer \n";
    $kbit = "Kmer/03.All_bit/kmer_15.bit";
    $kmer_size = 15;
}

# 构建kmer表
if(-e $ref)
{
    outShell2("ref.fa", "kmer.lst");

    $cmd .= "perl $BIN/sortSeq.pl $ref ref.fa 5 \n";
    $cmd .= "perl $BIN/Graph.pl pipe -i kmer.lst ";
    $cmd .= " -m 1 -k 15 -s 1,5 -d Kmer \n";
    $kbit = "Kmer/03.All_bit/kmer_15.bit";
    $kmer_size = 15;
}

# 压缩数据
outShell2("trinity sort.fa", "file.lst");

$cmd .= "perl $BIN/Compress.pl compress -i file.lst -m 1 ";
$cmd .=" -g $kbit -t $cutoff -n 1 -v 1 -k $kmer_size\n";

outShell2($cmd, "$outdir/removeDup.sh");
# 
# ==============================================================================

`sh removeDup.sh`;

#!/usr/bin/perl
use warnings;
use strict;
use FindBin qw($Bin);
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path getcwd);

BEGIN {
    push (@INC,"$Bin");
}
use Qsub;
use SeqIO;

=head1 Name

    Align_s.pl  --The Alignment Tool for illumina reads

=head1 Usage

    perl Align_s.pl [arguments] <reference.fa> <read.lst>

    Argument List:
                  -g <FILE>       the unique kmer(.h5 or .bit)
                  -k <INT>        the kmer size of the unique graph[17]
                  -u <INT>        the min unique kmer[3]

                  -s <INT>        split the reference file, in unit of M[100].
                  -n <INT>        max align number for query [2]
                  -t <INT>        thread number[4]

                  -d <DIR>        the output directory [Align]

                  -b <STR>        the pro_name [align]
                  -q <STR>        the queue of sge [dna.q,rna.q,reseq.q]
                  -p <STR>        the project of sge [og]

    <read.lst>:
                  pe450 read1
                  pe450 read2 

=head1 Example

perl Align -g k17.bit -k 17 -u 3 -n 10 ref.fa fq.lst

=cut

my(
$opt_g, $opt_k, $opt_u,   # kmer
$opt_s, $opt_s2,          # split size
$opt_j, $opt_n, $opt_x,   # scope
$opt_t, $opt_m,           # speed
$opt_c, $opt_d,           # output
$opt_b, $opt_q, $opt_p,   # qsub
$help
);

# Get parameters
# ==============================================================================
# |
GetOptions(
"g:s"     => \$opt_g,
"k:i"     => \$opt_k,
"u:i"     => \$opt_u,

"s:i"     => \$opt_s,
"s2:i"    => \$opt_s2,

"x:i"     => \$opt_x,
"j:i"     => \$opt_j,
"n:i"     => \$opt_n,

"t:i"     => \$opt_t,
"m:i"     => \$opt_m,

"c:s"     => \$opt_c,
"d:s"     => \$opt_d,

"b:s"     => \$opt_b,
"q:s"     => \$opt_q,
"p:s"     => \$opt_p,

"help|h"  => \$help,
);

my($ref, $query) = @ARGV;

checkParam();

$opt_g ||= -1;
$opt_k ||= 17;

$opt_g = abs_path($opt_g);

$opt_u ||= 3;
$opt_x ||= -1;

$opt_s ||= 100;
$opt_s *= 10**6;

$opt_s2 ||= 1000;
$opt_s2 *= 10**6;

$opt_j ||= 1;
$opt_n ||= 2;
$opt_t ||= 4;
$opt_m ||= 1;

$opt_d ||= "Align";

$opt_b ||= "Align";
$opt_q ||= "all.q";
$opt_p ||= "wrbio";

$ref   = abs_path($ref);
$query = abs_path($query);
$opt_d = abs_path($opt_d);
# |
# ==============================================================================


# Software path
# ==============================================================================
# |
my $SPLIT = "perl $Bin/splitFa.pl";
my $MERGE = "perl $Bin/mergeReads.pl";
my $BUILD = "$Bin/2bwt-builder";
my $ALIGN = "$Bin/align_bwt_s";
# |
# ==============================================================================


# Read files
# ==============================================================================
# |
my %read_hash;

my $lst_hdl = myOpen($query);
while (<$lst_hdl>) {
  chomp;
  my($lib, $file, $others) = split;
  die "no input" if(!-e $file);
  $file = abs_path($file);
  $read_hash{$lib} .= "$file\t";
}
close $lst_hdl;
# |
# ==============================================================================


# Direcotory
# ==============================================================================
#
my $shell_dir = "$opt_d/Shell";
my $result_dir = "$opt_d/Split_A";

mkdir "$opt_d";
mkdir($shell_dir);
mkdir($result_dir);

chdir($opt_d);
my $pwd = getcwd();
# |
# ==============================================================================



# build index
# ==============================================================================
#
my @refs = splitFile($ref, "Split_D", $opt_s, 20);

my $cmd_file = "$shell_dir/BUILD.sh";
my $ou_hdl = myOpen(">$cmd_file");
for(my $i=0; $i<scalar(@refs); $i++)
{
    my $ref = abs_path($refs[$i]);
    my $cmd = "$BUILD $ref";

    my $shell = "build_$i.sh";
    print $ou_hdl "$shell\t$cmd\n";
}
close $ou_hdl;

qsub($cmd_file, $shell_dir, "1G", 1,
            $opt_q, $opt_p, $opt_b, 50);
#
# ==============================================================================



# align file
# ==============================================================================
#
$cmd_file = "$shell_dir/ALIGN.SH";
$ou_hdl = myOpen(">$cmd_file");

my $merge_cmd = "$shell_dir/MERGE.SH";
my $cmd_hdl = myOpen(">$merge_cmd");

while (my($lib, $reads) = each %read_hash)
{
    my($read1, $read2) = split(/\s+/, $reads);
    for(my $j=0; $j<=$#refs; $j++)
    {
        # output dir
        my $out_dir = "$result_dir/$lib";
        mkdir($out_dir) unless (-e $out_dir);

        # command
        my $ref = abs_path($refs[$j]);

        my $cmd = "";
       
        # source
        $cmd .= "source $Bin/source.sh && ";
       
        # directory
        $cmd .= "cd $out_dir && ";
       
        # merge reads
        $cmd .= "$MERGE $read1 $read2 /dev/stdout | ";

        # do alignment
        $cmd .= "$ALIGN ";
        $cmd .= "-g $opt_g " if(-e $opt_g);
        $cmd .= "-k $opt_k -u $opt_u ";
        $cmd .= "-x $opt_x -j $opt_j -n $opt_n ";
        $cmd .= " -t $opt_t -m $opt_m ";
        $cmd .= "/dev/stdin $ref $out_dir/${lib}_$j.tab";

        # shell file
        my $align_sh = "align_${lib}_$j.sh";
        print $ou_hdl "$align_sh\t$cmd\n";
    }

    my $cmd = "cat $result_dir/$lib/*.tab > $pwd/$lib.tab && ";
    $cmd .= "sort -k 1,1 -k 2,2 -o $pwd/$lib.tab -S 1g $pwd/$lib.tab && ";
    $cmd .= "rm -r $result_dir/$lib\n";

    my $shll = "merge_$lib.sh";
    print $cmd_hdl "$shll\t$cmd\n";
}
close $ou_hdl;

# alignment
qsub($cmd_file, $shell_dir, "4G", $opt_t,
                $opt_q, $opt_p, $opt_b, 50);

# merge
qsub($merge_cmd, $shell_dir, "1G", 1,
                $opt_q, $opt_p, $opt_b, 50);

#
# ==============================================================================


# check parameters
# ==============================================================================
#
sub checkParam
{
	if ($help || @ARGV != 2) {
        die `pod2text $0`;
    }
}
#
# ==============================================================================

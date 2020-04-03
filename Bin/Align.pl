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

    Align.pl  --The Alignment Tool

=head1 Usage

    perl Align.pl [arguments] <reference.fa> <query.fa>

    Argument List:
                  -g <FILE>       the unique kmer(.h5 or .bit)
                  -k <INT>        the kmer size of the unique graph[17]
                  -u <INT>        the min unique kmer[3]

                  -s <INT>        split the reference file, in unit of M[100].
                  -s2<INT>        split the query file, in unit of M[1000].

                  -x <INT>        the scope to align[-1]
                  -j <INT>        the jump length to get kmer[1]
                  -n <INT>        max align number for query [20]

                  -t <INT>        thread number[4]
                  -m <INT>        align mode[1]
                                  1. align with LCS,for uncrrected reads
                                  2. align with kmer, for corrected reads

                  -d <DIR>        the output directory [Align]


                  -b <STR>        the pro_name [align]
                  -q <STR>        the queue of sge [dna.q,rna.q,reseq.q]
                  -p <STR>        the project of sge [og]


=head1 Example

perl Align -g k17.bit -k 17 -u 3 -n 10 ref.fa query.fa

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


$opt_n ||= 20;
$opt_t ||= 4;
$opt_m ||= 1;
$opt_j ||= 1 if($opt_m==1);
$opt_j ||= 50 if($opt_m==2);

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
my $BUILD = "$Bin/2bwt-builder";
my $ALIGN = "$Bin/align_bwt";
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


# prepare input data
# ==============================================================================
#
# split reference
my @refs = splitFile($ref, "Split_D", $opt_s, 20);

# split the query
my @query;
if($opt_s2>0)
{
    @query=splitFile($query, "Split_I", $opt_s2, 100);
}else{
    push(@query, $query);
}
# |
# ==============================================================================



# build index
# ==============================================================================
#
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
            $opt_q, $opt_p, $opt_b, 20);
#
# ==============================================================================



# align file
# ==============================================================================
#
$cmd_file = "$shell_dir/ALIGN.SH";
$ou_hdl = myOpen(">$cmd_file");

my $index = 0;
for(my $i=0; $i<=$#query; $i++)
{
    my $query = abs_path($query[$i]);
    for(my $j=0; $j<=$#refs; $j++)
    {
        # output dir
        $index++;
        my $idx = int($index/500);
        my $out_dir = "$result_dir/sub_$idx";
        mkdir($out_dir) unless (-e $out_dir);

        # output file
        my $result = basename("${query}_$j.al");
        my $ref = abs_path($refs[$j]);

        # command
        my $cmd = "source $Bin/source.sh && ";
        $cmd .= "cd $out_dir && ";
        $cmd .= "$ALIGN ";
        $cmd .= "-g $opt_g " if(-e $opt_g);
        $cmd .= "-k $opt_k -u $opt_u ";
        $cmd .= "-x $opt_x -j $opt_j -n $opt_n ";
        $cmd .= " -t $opt_t -m $opt_m ";
        $cmd .= "$query $ref $out_dir/$result";

        my $align_sh = "align_${i}_${j}.sh";
        print $ou_hdl "$align_sh\t$cmd\n";
    }
}
close $ou_hdl;

qsub($cmd_file, $shell_dir, "10G", $opt_t,
                $opt_q, $opt_p, $opt_b, 20);

#
# ==============================================================================


# merge file
# ==============================================================================
#
my $cmd = "cat $pwd/Split_A/sub_*/*.al > $pwd/align.al ";
$cmd .= " && rm -r $pwd/Split_* ";

my $merge_sh = "merge.sh";
$cmd_file = "$shell_dir/MERGE.SH";
outShell2("$merge_sh\t$cmd", $cmd_file);

qsub($cmd_file, $shell_dir, "1G", 1,
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

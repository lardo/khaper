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

die "perl $0 <dir_consensus> <queue:dna.q,rna.q,reseq.q> <Project:og> <pro_name:class_$$>\n" if(@ARGV==0);

my($dir_corr, $queue, $Project, $pro_name) = @ARGV;
$dir_corr = abs_path($dir_corr);

$queue ||= "all.q";
$Project ||= "wrbio";
$pro_name ||= "class_$$";

# Software path
# ==============================================================================
# |
my $PBDAGCON  = "perl $Bin/pbdagcon.pl";
# |
# ==============================================================================



# Do consensus
# ==============================================================================
#
mkdir("Shell");
chdir("Shell");

my $list = `ls $dir_corr/*/*.subread.fasta`;
my @files = split(/\n/, $list);
my @shell_cors;

# generate correct shells #
my $cmd = "";
my $index = 0;
foreach my $read (@files)
{
    my $sub_ref = $read;
    $sub_ref =~ s/\.subread//;

    my $dir = dirname($read);
    $cmd .= "cd $dir\n";
    $cmd .= "$PBDAGCON  $sub_ref $read 4 3\n";

    $index++;
    outCor()  if($index%20==0);
}
outCor();

# generate config file #
open FO, ">Correct.cfg" or die $!;
foreach my $shell(@shell_cors){
    print FO "$shell:5G\n";
}
close FO;

# qsub jobs
qsub2("Correct.cfg", 4, $queue, $Project, $pro_name, 50);

# cat files
$cmd  = "cd $dir_corr\n";
$cmd .= "cat $dir_corr/*/*.consensus > correct.fa";
system($cmd);

sub outCor
{
    return if($cmd eq "");
    my $index2 = int($index/20);
    my $shell = abs_path("pbdagcon_$index2.sh");

    $cmd .= "echo done";
    outShell2($cmd, $shell);
    push(@shell_cors, $shell);

    $cmd = "";
}
#|
#
# ==============================================================================


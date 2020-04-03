#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Cwd qw(abs_path getcwd);
use FindBin qw($Bin);
BEGIN {
    push (@INC,"$Bin");
}
use SeqIO;

=head1 Name

    pbsplit.pl  --split backbone and related reads

=head1 Usage

    perl pbsplit.pl [arguments]

    Arguments:
        -a <FILE>       the backbone file to be corrected
        -b <FILE>       the pacbio file used to correct backbone
        -l <FILE>       the relation list between backbone and pacbio reads
                        >bacbone1
                        pacb1
                        pacb2

        -d <STR>        the output dir

        -s <int>        the base to do consensens [2000000]


=head1 Example

	perl consensus.blasr.pl -a bacbone.fa -b pacbio.fa -l relation.lst -d Consensus

=cut

# Parameters
# ==============================================================================
#|
my($backbone, $pb, $relation, $outDir, $max_base);
my $help;

GetOptions(
    "a:s"     => \$backbone,
    "b:s"     => \$pb,
    "l:s"     => \$relation,
    "d:s"     => \$outDir,
    "s:i"     => \$max_base,
    "help|h"  => \$help,
);

checkParam();

$max_base ||= 2000000;
$outDir = abs_path($outDir);
#|
# ==============================================================================


# Directory
# ==============================================================================
#|
`rm -r $outDir` if(-e $outDir);
mkdir("$outDir");

# ==============================================================================


my($in_hdl, $ou_hdl, $ou2_hdl);
my %read_backbone;
# Record the related seqs with backbone
# ==============================================================================
#|
debug("record relationship");
$/ = ">";
$in_hdl = myOpen($relation);
while (<$in_hdl>)
{
    chomp;
    next if(length $_ == 0);

    my @array = split(/\n/, $_);
    my $backbone = $array[0];
    foreach my $id(@array[1..$#array])
    {
        $read_backbone{$id} .= "$backbone\t";
    }
}
close $in_hdl;
$/ = "\n";
#|
# ==============================================================================


my %backbone_hash;
my %len_hash;
# Record the backbone
# ==============================================================================
#|
my $bac_num = 0;
my($id, $seq);
$in_hdl = myOpen($backbone);
while (getSeq($in_hdl, \$id, \$seq)!=-1)
{
    $backbone_hash{$id} = $seq;

    my $len = length($seq);
    $len_hash{$id} = $len;
    $bac_num++;
}
close $in_hdl;
#|
# ==============================================================================


# Covert the format of pb reads
# ==============================================================================
#|
$in_hdl = myOpen($pb);
$ou_hdl = myOpen(">temp.tab");
while (getSeq($in_hdl, \$id, \$seq)!=-1)
{
    next if(!exists $read_backbone{$id});

    my $backbones = $read_backbone{$id};
    my @backbones = split(/\s+/, $backbones);
    for(my $i=0; $i<@backbones && $i<3; $i++)
    {
        my $backbone = $backbones[$i];
        print $ou_hdl "$backbone\t$id\t$seq\n";
    }
}
close $in_hdl;
close $ou_hdl;
`sort -k 1,1 -S 5g -o temp.tab temp.tab`;

my $temp = abs_path("temp.tab");
#|
# ==============================================================================


my @files;
# Split fasta file into a directory and with 20 files per subdirectory
# ==============================================================================
#|
debug("split files...");

chdir($outDir);

$bac_num = 0;

my $index = 0;
my $base  = 0;
my $pre_bac = "";

$in_hdl  = myOpen($temp);
$ou_hdl  = myOpen(">temp.fa");  # backbone
$ou2_hdl = myOpen(">temp2.fa"); # sub reads
while (<$in_hdl>)
{
    chomp;
    my($backbone, $id, $seq) = split;
    output() if($pre_bac ne $backbone && $base>$max_base);

    # output backbone
    if($pre_bac ne $backbone)
    {
        $base += $len_hash{$backbone};
        print $ou_hdl ">$backbone\n$backbone_hash{$backbone}\n";
        delete $backbone_hash{$backbone};
        $bac_num++;
    }
    # output subreads
    print $ou2_hdl ">$id\n$seq\n";
    $pre_bac = $backbone;
}
output();
close $in_hdl;

# 输出未覆盖的backbone #
while (my($id, $seq)=each %backbone_hash)
{
    $base += $len_hash{$id};
    print $ou_hdl ">$id\n$seq\n";
    print $ou2_hdl ">$id\n$seq\n";
    output() if($base>$max_base);
    $bac_num++;
}
output();
close $ou_hdl;
close $ou2_hdl;

`rm temp* $temp`;

sub output
{
    my $number = 20;

    return if($base==0);

    my $file1 = "split_$index.fasta";
    my $file2 = "split_$index.subreads.fasta";

    my $index2 = int($index/$number);
    my $dir = "dir_$index2";
    mkdir($dir) if(!-e $dir);

    # move file into directory
    close $ou_hdl;
    close $ou2_hdl;
    `mv temp.fa $dir/$file1`;
    `mv temp2.fa $dir/$file2`;

    push(@files, "$outDir/$dir/split_$index");

    # reset
    $ou_hdl  = myOpen(">temp.fa");
    $ou2_hdl = myOpen(">temp2.fa");

    $index++;
    $base = 0;
}
#|
# ==============================================================================


sub checkParam{
    if ($help || !defined $backbone)
    {
        die `pod2text $0`;
    }
}

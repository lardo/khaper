#!/usr/bin/perl
package SeqIO;
require Exporter;

use warnings;
use strict;
use File::Basename;
use Cwd qw(abs_path getcwd);
use List::Util qw(max min);

our @version = 1.0;
our @ISA     = qw(Exporter);
our @EXPORT  = qw(
                    splitFile
                    myOpen
                    getSeq
                    parseConfig
                    revCom
                    mergeConfigs
                    getFileType
                    debug
                    outShell2

                    align
);

# Parse config file
# ==============================================================================
#|
sub parseConfig
{
    my($file,$hash) = @_;
    my $in_hdl = myOpen($file);
    while (<$in_hdl>) {
        chomp;
        ${$hash}{$1}=$2 if(/(\S+)\s+(\S+)/);
        ${$hash}{$1}=$2 if(/(\S+)=(\S+)/);
    }
    close $in_hdl;
}
#|
# ==============================================================================

# Merge config files
# ==============================================================================
#|
sub mergeConfigs
{
    my($file_1, $file_2, $out_file) = @_;
    my %last_hash;
    my %first_hash;
    my %relation_hash;

    # get the last shells from file_1
    my $hdl_1 = myOpen($file_1);
    while (<$hdl_1>) {
        chomp;
        my($sh_1, $sh_2) = split;
        $relation_hash{$sh_1} = $sh_2 if(defined $sh_2);
    }
    close $hdl_1;

    $hdl_1 = myOpen($file_1);
    while (<$hdl_1>) {
        chomp;
        my($sh_1, $sh_2) = split;
        $last_hash{$sh_1} = 1 if(!exists $relation_hash{$sh_1});
        next if(!defined $sh_2);
        $last_hash{$sh_2} = 1 if(!exists $relation_hash{$sh_2});
    }
    close $hdl_1;
    %relation_hash = ();

    # get the first shells from file_2
    my $hdl_2 = myOpen($file_2);
    while (<$hdl_2>) {
        chomp;
        my($sh_1, $sh_2) = split;
        $relation_hash{$sh_2} = $sh_1 if(defined $sh_2);
    }
    close $hdl_2;

    $hdl_2 = myOpen($file_2);
    while (<$hdl_2>) {
        chomp;
        my($sh_1, $sh_2) = split;
        first_hash{$sh_1} if(!exists $relation_hash{$sh_1});
        next if(defined $sh_2);
    }
    close $hdl_2;

    # merge config file
    my $out_hdl = myOpen(">$out_file");
    $hdl_1 = myOpen($file_1);
    $hdl_2 = myOpen($file_2);
    while (<$hdl_1>) {
        print $out_hdl "$_";
    }
    while (<$hdl_2>) {
        print $out_hdl "$_";
    }
    foreach my $last_sh (keys %last_hash){
        foreach my $first_sh(keys %first_hash){
            print $out_hdl "$last_sh\t$first_sh\n";
        }
    }
    close $hdl_1;
    close $hdl_2;
    close $out_hdl;

}
#|
# ==============================================================================


# output msg
# ==============================================================================
#|
sub debug
{
    my $msg = shift;
    print STDERR "$msg\n";
}
#|
# ==============================================================================



# Reverse and Completement
# ==============================================================================
#|
sub revCom
{
    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/ATGCatgcNn/TACGtacgNn/;
    return $seq;
}
#|
# ==============================================================================



# Split fasta file into a directory and with 100 files per subdirectory
# ==============================================================================
#|
sub splitFile
{
	my @files;

    my($file, $outDir, $max_base, $file_num) = @_;
    $file_num ||= 200;

    # dir
    my $pwd = getcwd();
    mkdir $outDir;
    chdir($outDir);

    # split file
    my $index = 0;
    my $seq_base = 0;

    my ($id,$seq);
    my $in_hdl = myOpen($file);
    my $ou_hdl = myOpen(">temp.fa");
    while (getSeq($in_hdl,\$id,\$seq)!=-1)
    {
        my $seq_len = length($seq);
        print $ou_hdl ">$id\n$seq\n";
        $seq_base += length($seq);

        next if($seq_base<$max_base);
        $ou_hdl=output("temp.fa", $ou_hdl, \$seq_base, \$index, \@files, $file_num);
    }
    $ou_hdl=output("temp.fa", $ou_hdl, \$seq_base, \$index, \@files, $file_num);
    close $in_hdl;
    close $ou_hdl;

    `rm temp.fa` if(-e "temp.fa");
    chdir($pwd);
    return @files;
}

sub output
{
    my($infile, $ou_hdl, $base, $index, $files, $num) = @_;
    return if(${$base}==0);

    my $dir_idx = int(${$index}/$num);
    my $out_dir = abs_path("Split_$dir_idx");
    mkdir($out_dir) if(!-e $out_dir);
    close $ou_hdl;
    `mv $infile $out_dir/split.${$index}.fasta`;
    push(@{$files}, "$out_dir/split.${$index}.fasta");

    # reset
    ${$index}++;
    ${$base} = 0;
    $ou_hdl = myOpen(">temp.fa");
    return $ou_hdl;
}
#|
# ==============================================================================



# Open file in different format
# ==============================================================================
#|
sub myOpen
{
    my $file = shift;
    my $handle;
    if($file =~ /\.gz$/ && $file !~ /\>/){
        open $handle,"gzip -dc $file|" or die "can't open $file";
    }else{
         open $handle,"$file" or die "can't open $file";
    }
    return $handle;
}
#|
# ==============================================================================


# get the format of file
# 1:fasta 2:fastq 0:undefined
# ==============================================================================
#|
sub getFileType
{
    my $file = shift;
    my $handle = myOpen($file);
    my $line = <$handle>;
    close $handle;

    return 2 if($line =~ /^@/);
    return 1 if($line =~ /^>/);

    return 0;
}
#|
# ==============================================================================



# Get one sequence from fasta file each time
# ==============================================================================
#|
sub getSeq
{
    my ($hdl, $id, $seq) = @_;
    $/  = ">";
    while (<$hdl>) {
        chomp;
        next if(length $_ == 0);

        my @array = split(/\n/,$_);
        my $id1 = $array[0];
        $id1 =~ /(\S+)/;
        ${$id} = $1;
        ${$seq} = join("",@array[1..$#array]);

        $/ = "\n";

        return 1;
    }
    $/ = "\n";
    return -1;
}
#|
# ==============================================================================



# 输出SHELL文件
# ==============================================================================
# |
sub outShell2{
  my($cmd, $outfile, $app) = @_;
  $app ||= 0;

  my $out_hdl;
  $out_hdl = myOpen(">$outfile")  if($app==0);
  $out_hdl = myOpen(">>$outfile") if($app==1);

  print $out_hdl "$cmd\n";
  close $out_hdl;
}
# |
# ==============================================================================


# 比对
# ==============================================================================
# |

sub getComKmer
{
    my($seq1, $seq2, $k, $jump) = @_;
    $k ||= 17;
    $jump ||= $k;

    my %kmer_hash;
    for(my $i=0; $i<length($seq1)-$k; $i+=$jump)
    {
        my $kmer = substr($seq1, $i, $k);
        next if($kmer=~/N/);
        $kmer_hash{$kmer} = $i;
    }
    my $num = keys %kmer_hash;

    my @matrix;
    for(my $i=0; $i<length($seq2)-$k; $i++)
    {
        my $kmer = substr($seq2, $i, $k);
        next if(!exists $kmer_hash{$kmer});

        push(@matrix, "$kmer_hash{$kmer} $i");
    }
    return @matrix;
}

sub findAnchor{
    my $array = shift;
    my $anchor = "-1 -1";
    my $total = scalar @{$array};
    my $best_match = 0;

    my %match_hash;
    for(my $i=0; $i<@{$array}; $i++)
    {
        next if(exists $match_hash{$i});

        my $match = 1;
        my($x1, $y1) = split(/\s+/, ${$array}[$i]);

        $match_hash{$i} = 1;
        my @match_arr = ("$x1 $y1");
        for(my $j=$i+1; $j<@{$array}; $j++)
        {
            my($x1, $y1) = split(/\s+/, $match_arr[-1]);
            my($x2, $y2) = split(/\s+/, ${$array}[$j]);
            next if($x2==$x1);

            my $slope = ($y2-$y1)/($x2-$x1);
            next if($slope<0.9 or $slope>1.1);

            $match++;
            $match_hash{$j} = 1;
            push(@match_arr, "$x2 $y2");
        }

        if($match>$best_match){
            $best_match=$match;
            $anchor = "$x1 $y1";
            last if($best_match/$total>0.6);
        }
    }
    return $anchor;
}

sub simi{
    my($seq1, $seq2, $k) = @_;
    my $len = length($seq1);
    my $span = int($len/100)+1;
    my ($match, $total)= (0, 1);
    for(my $i=0; $i<$len-$k; $i+=$span)
    {
        $total++;
        my $kmer = substr($seq1, $i, $k);
        my $start = int(0.9*$i);

        my $len_ = int(0.2*$i)+2*$k;
        next if($start+$len_>$len);

        my $sub_2 = substr($seq2, $start, $len_);

        $match++ if(index($sub_2, $kmer)!=-1);
    }
    my $simi = $match/$total;
    return $simi;
}

sub align{
   my($seq1, $seq2, $k, $jump, $sort) = @_;

   $k    ||= 14;
   $jump ||= $k;
   $sort ||= 1;

   my @matrix = getComKmer($seq1, $seq2, $k, $jump);
   @matrix = sort by_pos1_inc(@matrix);
   @matrix = reverse(@matrix) if($sort == -1);

   my $anchor = findAnchor(\@matrix);
   return (0, $anchor) if($anchor eq "-1 -1");

   my($x, $y) = split(/\s+/, $anchor);

   my $len1 = length($seq1);
   my $len2 = length($seq2);
   my $tail1 = $len1-$x;
   my $tail2 = $len2-$y;

   my $head = min($x,$y);
   my $tail = min($tail1, $tail2);

   my $ovl = $head+$tail;
   my $sub_1 = substr($seq1,$x-$head,$ovl);
   my $sub_2 = substr($seq2,$y-$head,$ovl);

   my $similarity = simi($sub_1, $sub_2, 7);
   return ($similarity, $anchor);
}

sub by_pos1_inc{
    my($pos1, $pos2) = split(/\s+/,$a);
    my($pos3, $pos4) = split(/\s+/,$b);
    return $pos1<=>$pos3;
}
# |
# ==============================================================================

1;

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
use SeqIO;

die "perl $0 <genome.fa> <out> <len:10000> <cov:20>\n" if(@ARGV==0);

my($genome, $out, $cut, $cov) = @ARGV;
$cut ||= 10000;
$cov ||= 20;

my($id, $seq);
my $in_hdl = myOpen($genome);
my $ou_hdl = myOpen(">$out");
while (getSeq($in_hdl, \$id, \$seq)!=-1)
{
    my $len = length($seq);
    next if($len<$cut);

    my $total = 0;
    while ($total/$len<$cov) {
        my $start = int(rand($len));
        my $sub = substr($seq, $start, $cut);
        next if(length($sub)<500);
        $total += length($sub);
        my $new_id = "${id}_$start";
        print $ou_hdl ">$new_id\n$sub\n";
    }
}
close $in_hdl;
close $ou_hdl;

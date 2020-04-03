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


die "perl $0 <seq.fa> <out> <memory(G):1>\n" if(@ARGV==0);

my($file, $out, $memory) = @ARGV;
$memory ||= 1;

# covert format
my $in_hdl = myOpen($file);
open FO,">temp.$$" or die $!;
my($id, $seq);
while(getSeq($in_hdl,\$id,\$seq)!=-1)
{
	my $len = length($seq);
	print FO "$id\t$len\t$seq\n";		
}
close $in_hdl;
close FO;

# sort by length
my $cmd = "sort -k 2,2nr -T ./ -S ${memory}G -o temp.$$ temp.$$";
system($cmd);

# convert formate
$in_hdl = myOpen("temp.$$");
open FO,">$out" or die $!;
while(<$in_hdl>){
    chomp;
    my($id,$len,$seq) = split;
    print FO ">$id\n$seq\n";
}
close $in_hdl;
close FO;

$cmd = "rm temp.$$";
system($cmd);

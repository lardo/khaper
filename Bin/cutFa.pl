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

die "perl $0 <in.fasta> <out.fasta> <cut_len:20000>\n" if(@ARGV==0);

my($infile, $outfile, $cut_len) = @ARGV;
$cut_len ||= 20000;

my $in_hdl  = myOpen($infile);
my $out_hdl = myOpen(">$outfile");

my ($id,$seq);
while (getSeq($in_hdl,\$id,\$seq)!=-1) 
{
	my $len = length($seq);
	if($len<=$cut_len || $cut_len<=0)
	{
		print $out_hdl ">$id\n$seq\n";
		next;
	}
	for(my $i=0;$i<$len;$i+=$cut_len)
	{
		my $sub_seq = substr($seq,$i,$cut_len);
		print $out_hdl ">${id}_$i\n$sub_seq\n";
	}
	
}
close $in_hdl;
close $out_hdl;
$ou_hdl;

sub formatSeq{
    my $seq = shift;
    my $len = length($seq);
    my $sub_len = 80;
    for(my $i=0; $i<$len; $i+=$sub_len)
    {
        my $sub_seq = substr($seq, $i, $sub_len);
        print $ou_hdl "$sub_seq\n";
    }
}

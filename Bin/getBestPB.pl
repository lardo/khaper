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

die "perl $0 <in.fasta> <out.fasta> \n" if(@ARGV==0);

my($infile, $outfile) = @ARGV;

my $in_hdl  = myOpen($infile);
my $out_hdl = myOpen(">$outfile");

my ($id,$seq);
my $pre_id = "";
my @array;
while (getSeq($in_hdl,\$id,\$seq)!=-1) 
{
    my $cell = $id;
    $cell =~ s/\d+_\d+$//g;
    output() if($cell ne $pre_id);
    push(@array, "$id $seq");
    $pre_id = $cell;
}
output();
close $in_hdl;
close $out_hdl;

sub by_len{
    my($id1, $seq1) = split(/\s+/, $a);
    my($id2, $seq2) = split(/\s+/, $b);
    my $len1 = length($seq1);
    my $len2 = length($seq2);
    return $len2<=>$len1;
}

sub output{
    return if(@array==0);
    @array = sort by_len @array;
    my($id, $seq) = split(/\s+/, $array[0]);
    print $out_hdl ">$id\n$seq\n";
    @array = ();
}

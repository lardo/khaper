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


die "perl $0 <file.lst> <min_len> <out>\n" if(@ARGV==0);

my($list, $min_len, $out) = @ARGV;

open LIST,"$list" or die $!;
open FO,">$out" or die $!;
open LOG,">$out.log" or die $!;
my $index = 0;
while(<LIST>){
	chomp;
	my $fa = $_;
	my $in_hdl = myOpen($fa);
	my($id, $seq);
	while(getSeq($in_hdl,\$id,\$seq)!=-1)
	{
		my $len = length($seq);
		next if($len<$min_len);
        $index++;
        my $id2 = "C_$index";
        $seq =~ s/\r//g;
        $seq =~ s/[^ATGCN]//g;
		print FO ">$id2\n$seq\n";
        print LOG "$id2\t$id\n";		
	}
	close $in_hdl;
}
close LIST;
close FO;
close LOG;

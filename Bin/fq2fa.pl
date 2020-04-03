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

die "perl $0 <file.fq> <outfile>\n" if(@ARGV==0);

my($fq, $outfile, $min_len) = @ARGV;

my $in_hdl = myOpen($fq);
my $ou_hdl = myOpen(">$outfile");
while(my $id=<$in_hdl>){
    my $seq = <$in_hdl>;
    <$in_hdl>;
    <$in_hdl>;
    chomp $id;
    chomp $seq;

    $id =~ s/^@/>/;
    print $ou_hdl "$id\n$seq\n";
}
close $in_hdl;
close $ou_hdl;

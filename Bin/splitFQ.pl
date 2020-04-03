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

die "perl $0 <file.lst> <outfile> <min_len:5000>\n" if(@ARGV==0);

my($list, $outfile, $min_len) = @ARGV;
$min_len ||= 5000;

# Only output one read per cell
# ==============================================================================
# |
my @array = ();
my $pre_cell = "-1";

my $list_hdl = myOpen($list);
my $out_hdl = myOpen(">$outfile");

while (my $file=<$list_hdl>) 
{
    chomp $file;
    my $file_hdl = myOpen($file);
    while (my $id=<$file_hdl>) 
    {
        my $seq = <$file_hdl>;
        if($id=~/^@/){
            <$file_hdl>;
            <$file_hdl>;
        }

        chomp $id;
        chomp $seq;
        my $len = length($seq);
        next if($len<$min_len);

        $id =~ s/^@/>/;
        $id = (split(/\s+/, $id))[0];
        print $out_hdl "$id\n$seq\n"; 
    }
    close $file_hdl;
}
close $list_hdl;
close $out_hdl;
# |
# ==============================================================================


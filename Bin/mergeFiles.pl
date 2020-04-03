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

die "perl $0 <file.lst> <filt_len> <outfile>\n" if(@ARGV==0);

my($list, $filt_len, $outfile) = @ARGV;

my $list_hdl = myOpen($list);
my $out_hdl = myOpen(">$outfile");
while (my $file=<$list_hdl>) 
{
    chomp $file;
    my $type = getType($file);
    if($type eq "non")
    {
    	die "please check the format of $file\n";
    }

    my $file_hdl = myOpen($file);
    outFQ($file_hdl, $out_hdl) if($type eq "fastq");
    outFA($file_hdl, $out_hdl) if($type eq "fasta");
    close $file_hdl;
}
close $list_hdl;
close $out_hdl;

sub getType{
	my $file = shift;
	my $in_hdl = myOpen($file);
	my $line = <$in_hdl>;
	close $in_hdl;

	return "fasta" if($line =~ /^>/);
	return "fastq" if($line =~ /^@/);
	return "non";
}

sub outFQ{
	my($file_hdl,$out_hdl) = @_;
	while (my $id=<$file_hdl>) 
	{
        my $seq = <$file_hdl>;
        chomp $id;
        chomp $seq;
        my $len = length($seq);
        next if($len<$filt_len);
        $id =~ s/^[@|>]//;
        $id = (split(/\s+/,$id))[0];
        print $out_hdl ">$id\n$seq\n";

        <$file_hdl>;
        <$file_hdl>;
    }
}

sub outFA{
	my($file_hdl,$out_hdl) = @_;
	my($id,$seq);
	while (getSeq($file_hdl,\$id,\$seq)!=-1) 
	{
        my $len = length($seq);
        next if($len<$filt_len);
        print $out_hdl ">$id\n$seq\n";
    }
}

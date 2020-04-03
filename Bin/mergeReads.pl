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

die "perl $0 <read1.fasta(q)> <read2.fasta(q)> <outfile>\n" if(@ARGV==0);

my($read1, $read2, $outfile) = @ARGV;

my $index = 0;
my @reads = ($read1, $read2);
my $out_hdl = myOpen(">$outfile");
foreach my $file(@reads)
{ 
    my $type = getType($file);
    if($type eq "non")
    {
    	die "please check the format of $file\n";
    }

    $index++;
    my $file_hdl = myOpen($file);
    outFQ($file_hdl, $out_hdl, $index) if($type eq "fastq");
    outFA($file_hdl, $out_hdl, $index) if($type eq "fasta");
    close $file_hdl;
}
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
	my($file_hdl,$out_hdl,$prefix) = @_;
    my $index = 0;
	while (my $id=<$file_hdl>) 
	{
        my $seq = <$file_hdl>;
        
        chomp $seq;
        $index++;
        $id = "$index $prefix";
        print $out_hdl ">$id\n$seq\n";

        <$file_hdl>;
        <$file_hdl>;
    }
}

sub outFA{
	my($file_hdl,$out_hdl,$prefix) = @_;
	my $index = 0;
    my($id,$seq);
	while (getSeq($file_hdl,\$id,\$seq)!=-1) 
	{
        $index++;
        $id = "$index $prefix";
        print $out_hdl ">$id\n$seq\n";
    }
}

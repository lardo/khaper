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
my $out1_hdl = myOpen(">$outfile.1.fa");
my $out2_hdl = myOpen(">$outfile.2.fa");
while (my $file=<$list_hdl>)
{
    chomp $file;
    my $type = getType($file);
    if($type eq "non")
    {
    	die "please check the format of $file\n";
    }

    my $file_hdl = myOpen($file);
    outFQ($file_hdl, $out1_hdl, $out2_hdl) if($type eq "fastq");
    outFA($file_hdl, $out1_hdl, $out2_hdl) if($type eq "fasta");
    close $file_hdl;
}
close $list_hdl;
close $out1_hdl;
close $out2_hdl;

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
	my($file_hdl, $out1_hdl, $out2_hdl) = @_;
	while (my $id=<$file_hdl>)
	{
        my $seq = <$file_hdl>;
        chomp $id;
        chomp $seq;
        my $len = length($seq);
        next if($len<$filt_len);
        $id =~ s/^[@|>]//;
        $id = (split(/\s+/,$id))[0];

        # simulate solexa reads
        for(my $i=0; $i<$len; $i+=500){
            my $end = $i+$len;
            last if($end>=$len);

            my $read1 = substr($seq, $i, 150);
            my $read2 = substr($seq, $i+$filt_len-150, 150);
            $read2 = revCom($read2);

            print $out1_hdl ">$id:$i:$filt_len 1\n$read1\n";
            print $out2_hdl ">$id:$i:$filt_len 2\n$read2\n";
        }

        <$file_hdl>;
        <$file_hdl>;
    }
}

sub outFA{
	my($file_hdl, $out1_hdl, $out2_hdl) = @_;
	my($id,$seq);
	while (getSeq($file_hdl,\$id,\$seq)!=-1)
	{
        my $len = length($seq);
        next if($len<$filt_len);

        # simulate solexa reads
        for(my $i=0; $i<$len; $i+=500){
            my $end = $i+$filt_len;
            last if($end>=$len);

            my $read1 = substr($seq, $i, 150);
            my $read2 = substr($seq, $i+$filt_len-150, 150);
            $read2 = revCom($read2);

            print $out1_hdl ">$id:$i:$filt_len 1\n$read1\n";
            print $out2_hdl ">$id:$i:$filt_len 2\n$read2\n";
        }
    }
}

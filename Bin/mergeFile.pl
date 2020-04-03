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


die "perl $0 <infile.cfg>\n" if(@ARGV==0);

my($infile) = @ARGV;

my $in_hdl = myOpen($infile);
while (<$in_hdl>)
{
    chomp;
    debug($_);
    my($prefix, $read1, $read2 ) = split;

    # paired-end reads
    if(-e $read1 && -e $read2)
    {
        my $read1_hdl = myOpen($read1);
        my $read2_hdl = myOpen($read2);

        while (my $id1 = <$read1_hdl>) 
        {
            my $id2 = <$read2_hdl>;
            my $seq1 = <$read1_hdl>;
            my $seq2 = <$read2_hdl>;
            if($id1 =~ /^@/)
            {
                <$read1_hdl>;
                <$read1_hdl>;
                <$read2_hdl>;
                <$read2_hdl>;
            }

            chomp $id1;
            chomp $seq1;
            chomp $seq2;
        
            # convert seq
            my $len1 = length($seq1);
            my $len2 = length($seq2);
            my $seq = "$seq1$seq2";
            
            # convert id
            $id1 =~ /^[>|@](\S+)/;
            my $id = $1;
            $id .= "\t$prefix:$len1:$len2";

            # output sequence
            print ">$id\n$seq\n";
        }
        close $read1_hdl;
        close $read2_hdl;
    }
    
    # single end
    if(!defined $read2 || !-e $read2)
    {
        my $seq = 'N'x10000;
        print ">temp_1\t$prefix:0\n$seq\n";
        print ">temp_2\t$prefix:0\n$seq\n"; 
    
        my $read_hdl = myOpen($read1);
        while (my $line = <$read_hdl>) 
        {
            chomp $line;
            if($line =~ /^@\w+/ || $line =~ /^>\w+/)
            {
        	   $line .= "\t$prefix:0";
            }
            print "$line\n";
        }
        close $read_hdl;
    } 
}
close $in_hdl;

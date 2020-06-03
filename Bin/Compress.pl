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

=head1 Name

    Compress.pl  --The compress module for fastq(a) files

=head1 Usage

    perl Compress.pl <command> [arguments]

    Command should be one of the following command. Arguments depend on specific command.

    Command List:
                compress  Compress reads
                          Arguments:
                          -i <FILE>       the file list with format like:                                          
                                          Out_prefix read_1 read_2 
                                          Out_prefix read_n

                          -g <FILE>       the k-mer file from reads(.h5 or .bit)
                          -k <INT>        the k-mer size
                          -m <INT>        the min kmer num in read[3] 
                    
                          -t <FLOAT>      trim the reads if (occupied_kmer/unique_kmer)>=this_value[0.7]
                          -n <INT>        thread number [8]

                          -v <INT>        output the log file[0]
                                          0: don't record the log
                                          1: record the log

                pipe  Compress reads with different insert size
                          -i <FILE>       the file list with paired-end or mate reads                                         
                                          read_1  500
                                          read_2  500
                                          read_3  1000
                                          read_4  1000

                          -g <FILE>       the k-mer file from reads(.h5 or .bit)
                          -k <INT>        the k-mer size
                          -m <INT>        the min kmer num in read[3] 
                    
                          -t <FLOAT>      trim the reads if (occupied_kmer/unique_kmer)>=this_value[0.7]
                          -n <INT>        thread number [8]

                          -d <INT>        need to compress[1]
                                          1: need to compress
                                          0: out put all reads 

=head1 Example

 Compress data:
    perl Compress.pl compress -i file.lst -g kmer_17.h5 -k 17 -m 3 -t 0.7 -n 16

=cut


my($command) = @ARGV;

my( 
    $opt_i, $opt_w, $opt_j,
    $opt_g, $opt_k, $opt_m,
    $opt_t, $opt_n, $opt_d,
    $opt_v, $help 
);	

# Get parameters
# ==================================================================================================
# |
$opt_t = 0.7;
$opt_m = 3;
$opt_n = 8;
$opt_w = 500;
$opt_j = -1;
$opt_d = 1;
$opt_v = 0;

GetOptions(
    "i:s"     => \$opt_i,
    "w:i"     => \$opt_w,
    "j:i"     => \$opt_j,
    "d:i"     => \$opt_d,
    "g:s"     => \$opt_g,
    "k:i"     => \$opt_k,
    "m:i"     => \$opt_m,
    "t:f"     => \$opt_t,
    "n:i"     => \$opt_n,
    "v:i"     => \$opt_v,
    "help|h"  => \$help,
);

checkParam();
$opt_g = abs_path($opt_g);
# |
# ==================================================================================================


# Software path
# ==================================================================================================
# |
my $MERGE = "perl $Bin/mergeFile.pl ";
my $SPLIT = "perl $Bin/splitFA.pl ";
my $CMP = "$Bin/compress";
# |
# ==================================================================================================

# my $cmd = "source $Bin/source.sh\n";
my $cmd = "";

# compress data
# ==================================================================================================
# |
if($command eq "compress")
{
    $cmd .= "$MERGE $opt_i ";
    $cmd .= "| $CMP -g $opt_g -k $opt_k  -i /dev/stdin -o /dev/stdout ";
    $cmd .= "-m $opt_m -t $opt_n -c $opt_t -v $opt_v -d $opt_d ";
    $cmd .= "| $SPLIT /dev/stdin\n";

    debug($cmd);
    system($cmd);
}
# |
# ==================================================================================================

# compress data with each insert size
# ==================================================================================================
# |
if($command eq "pipe")
{
	# record insert size with related reads #
	my %insert_hash;
  my $in_hdl = myOpen($opt_i);
  while (my $read1 = <$in_hdl>) 
  {
    my $read2 = <$in_hdl>;
    chomp $read1;
    chomp $read2;
    my($fq1, $insert1) = split(/\s+/, $read1);
    my($fq2, $insert2) = split(/\s+/, $read2);

    $insert_hash{$insert1} .= "$fq1 $fq2\n";
  }
  close $in_hdl;

  my $qsub_sh = "qsub.sh";
  my $lib_lst = "sample.info";
  outShell2("", $qsub_sh);
  outShell2("", $lib_lst);

  # generate script for each insert size #
  while (my($insert, $reads) = each %insert_hash) 
  {
    my $dir = "lib_$insert";
    mkdir $dir;
    # file list #
    my @array = split(/\n/, $reads);
    open FO,">$dir/file.lst";
    foreach(@array){
    		my($fq1, $fq2) = split;
    		$fq1 = abs_path($fq1);
    		$fq2 = abs_path($fq2);
    		print FO "$dir\t$fq1\t$fq2\n";
    }
    close FO;
 		
   	# compress shell #
  	$dir = abs_path($dir);

  	my $cmd = "cd $dir && ";
    $cmd .= "source $Bin/source.sh && ";
    $cmd .= "perl $0 compress -i file.lst -g $opt_g -k $opt_k -m $opt_m -t $opt_t -n $opt_n -d $opt_d ";
    $cmd .= "&& echo done";
    outShell2($cmd, "$dir/lib_$insert.sh");

    # qsub shell #
    my $qsub = "qsub -cwd -l vf=10G,p=$opt_n -q all.q -S /bin/bash $dir/lib_$insert.sh";
    outShell2($qsub, $qsub_sh, 1);

    # output file #
    outShell2("lib_$insert $dir/lib_$insert.1.fasta.gz $dir/lib_$insert.2.fasta.gz $insert", $lib_lst, 1);
  }
}
# |
# ==================================================================================================


sub checkParam{
	if ($help || @ARGV != 1) {
    die `pod2text $0`;
  }
}


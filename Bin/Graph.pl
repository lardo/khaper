#!/usr/bin/perl
use warnings;
use strict;
use FindBin qw($Bin $RealScript);
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path getcwd);

BEGIN {
    push (@INC,"$Bin");
}
use SeqIO;

=head1 Name

    Graph.pl  --The De novo tool to build k-mer graph

=head1 Usage

    Graph.pl <command> [arguments]

    Command should be one of the following command. Arguments depend on specific command.

    Command List:
                count   Record k-mer and related occurrence by jellyfish
                        Arguments:
                          -i <FILE>       the file list to count kmer
                          -k <INT>        the k-mer size to store occupied k-mer [17]
                          -m <INT>        the mininum occurrence of kmer [3]

                graph   Build graph by using GATB tookit(for Lordec)
                          Arguments:
                          -i <FILE>       the kmer table file with format "ATGC 1"
                          -m <INT>        the mininum occurrence of kmer [3]
                          -x <INT>        the maxinum occurrence of kmer [-1]
                          -k <INT>        the k-mer size [17]

                bit     Record k-mer into bitset, this method is for k<=17.
                          -i <FILE>       the kmer table file with format "ATGC 1"
                          -m <INT>        the mininum occurrence of kmer [3]
                          -x <INT>        the maxinum occurrence of kmer [-1]
                          -k <INT>        the k-mer size to store occupied k-mer [17]

                pipe    Combine step(s) above
                          -i <FILE>       the file list
                          -m <INT>        the mininum occurrence of kmer [3]
                          -k <INT>        the k-mer size to store occupied k-mer [17]
                          -s <INTs>       the step to do
                                          1: count k-mer by jellyfish

                                          2: record unique k-mer into .h5 file
                                          3: record unique k-mer into .bit file

                                          4: record all k-mer into .h5 file
                                          5: record all k-mer into .bit file

                                          6: record all kmer into .bit with -m is 0.5 the peak

                                          7: get the genome size, repeate rate and hete rate

                -d      The output directory

=head1 Example

For k=17, we recommend:

  perl Graph.pl pipe -i fq.lst -m 2 -k 17 -s 1,3,5 -d Kmer_17

For k>17, we recommend:

  perl Graph.pl pipe -i fq.lst -m 2 -k 23 -s 1,2,4 -d Kmer_23

=cut


my($command) = @ARGV;

my(
$opt_i, $opt_m, $opt_x, $opt_k,
$opt_d, $opt_s, $help
);
# Get parameters
# ==============================================================================
# |
GetOptions(
"i=s"     => \$opt_i,
"m:i"     => \$opt_m,
"x:i"     => \$opt_x,
"k=i"     => \$opt_k,
"s:s"     => \$opt_s,
"d=s"     => \$opt_d,
"help|h"  => \$help,
);

$opt_m ||= 3;
$opt_x ||= -1;
$opt_k ||= 17;

checkParam();
# |
# ==============================================================================


# Software path
# ==============================================================================
# |
# jellyfish
my $KMER_COUNT = abs_path("$Bin/kmer_counter.pl");
# convert jellyfish format to fasta
my $KMER2FA = "perl $Bin/dump2fa.pl ";
# record kmer by GATB
my $GRAPH = "$Bin/buildGraph ";
# record kmer bit bitset
my $BUILD = "$Bin/buildIndex";

# estimate the genome size, repeate rate and
my $FORMAT = "perl $Bin/format.pl";
my $ANALYSE = "perl $Bin/GenomeAnalysis.pl";
# |
# ==============================================================================

my $cmd = "source $Bin/source.sh\n";


# Count Kmers By Jellyfish
# ==============================================================================
# |
if($command eq "count")
{
  # dir
  $opt_d ||= "Jellyfish";
  mkdir("$opt_d") if(!-e $opt_d);

  # kmer list
  my $list = abs_path("$opt_d/kmer.lst");
  my $out_hdl = myOpen(">$list");
  my $in_hdl = myOpen($opt_i);
  while (<$in_hdl>)
  {
      chomp;
      next if(length($_)==0);
      my $file = abs_path($_);
      print $out_hdl "$file\n";
  }
  close $in_hdl;
  close $out_hdl;

  $cmd .= "cd $opt_d\n";
  $cmd .= "perl $KMER_COUNT --k $opt_k --memory 10 --Prefix k$opt_k $list\n";
  $cmd .= "sh k$opt_k.sh\n";

  # save disk
  $cmd .= "rm k$opt_k.jf\n";
  $cmd .= "awk '\$2>1' k$opt_k.dump > filt.dump && rm k$opt_k.dump\n" if($opt_m>1);
  $cmd .= "mv k$opt_k.dump filt.dump\n" if($opt_m==1);
  $cmd .= "gzip  -f filt.dump\n";

  debug("$cmd");
  system($cmd);
}
# |
# ==============================================================================




# build graph
# ==============================================================================
# |
if($command eq "graph")
{
    # dir
    $opt_d ||= "GATB";
    mkdir("$opt_d") if(!-e $opt_d);
    my $outfile = abs_path("$opt_d/kmer_$opt_k");

    # input
    $opt_i = abs_path($opt_i);
    #convert to fasta
    $cmd .= "$KMER2FA $opt_i $opt_m $opt_x ";
    #build graph
    $cmd .= " | $GRAPH -k $opt_k -m 1 -i /dev/stdin -o $outfile -T 1";

    debug("$cmd");
    system($cmd);
}
# |
# ==============================================================================



# record unique k-mer into .bit file
# ==============================================================================
# |
if($command eq "bit")
{
    # the slope of kmer
    die "the kmer size must <=18" if($opt_k>18);

    # dir
    $opt_d ||= "Bit";
    mkdir("$opt_d") if(!-e $opt_d);

    # dir
    mkdir("$opt_d") if(!-e $opt_d);
    my $outfile = abs_path("$opt_d/kmer_$opt_k.bit");

    # input
    $opt_i = abs_path($opt_i);

    # build graph
    $cmd .= "$BUILD -d $opt_i -k $opt_k -m $opt_m -x $opt_x -b $outfile\n";

    debug("$cmd");
    system($cmd);
}
# |
# ==============================================================================



# pipeline
# ==============================================================================
# |
if($command eq "pipe")
{
    # dir
    $opt_d ||= "Kmer_$opt_k";
    mkdir("$opt_d") if(!-e $opt_d);

    my $count_dir = abs_path("$opt_d/01.KmerCount");
    my $dir_uni_h5  = abs_path "$opt_d/02.Uinque_h5";
    my $dir_uni_bit = abs_path "$opt_d/02.Uinque_bit";
    my $dir_all_h5  = abs_path "$opt_d/03.All_h5";
    my $dir_all_bit = abs_path "$opt_d/03.All_bit";
    my $dir_peak_bit = abs_path "$opt_d/04.peak_bit";
    my $dir_analyse = abs_path "$opt_d/05.Analyse";

    # input
    $opt_i = abs_path($opt_i);

    # steps
    my %step_hash;
    my @steps = split(/,/,$opt_s);
    $step_hash{$_}=1 foreach(@steps);

    # count kmer by jellyfish
    if(exists $step_hash{1})
    {
        my $cmd = "perl $0 count -k $opt_k -i $opt_i -m $opt_m -d $count_dir\n";
        debug($cmd);
        system($cmd);
    }

    # store unique kmer by GATB
    if(exists $step_hash{2})
    {
        my $input = "$count_dir/filt.dump.gz";
        my $histo = "$count_dir/k$opt_k.histo";

        my $peak = getPeak($histo);
        my $max = int(1.5*$peak)+1 if($peak>0);

        debug("peak: $peak");
        die "can't find peak" if($peak<0);

        my $cmd = "perl $0 graph -k $opt_k -i $input ";
        $cmd .= "-m $opt_m -x $max -d $dir_uni_h5\n";
        system($cmd);
    }

    # store the unique kmer by bitset
    if(exists $step_hash{3})
    {
        my $input = "$count_dir/filt.dump.gz";
        my $histo = "$count_dir/k$opt_k.histo";

        my $peak = getPeak($histo);
        my $max = int(1.5*$peak)+1 if($peak>0);

        debug("peak: $peak");
        die "can't find peak" if($peak<0);

        my $cmd = "perl $0 bit -k $opt_k -i $input ";
        $cmd .= "-m $opt_m -x $max -d $dir_uni_bit\n";
        system($cmd);
    }

    if(exists $step_hash{4})
    {
        my $input = "$count_dir/filt.dump.gz";
        my $histo = "$count_dir/k$opt_k.histo";

        my $cmd = "perl $0 graph -k $opt_k -i $input ";
        $cmd .= "-m $opt_m -x $opt_x -d $dir_all_h5\n";
        system($cmd);
    }

    if(exists $step_hash{5})
    {
        my $input = "$count_dir/filt.dump.gz";
        my $histo = "$count_dir/k$opt_k.histo";

        my $cmd = "perl $0 bit -k $opt_k -i $input ";
        $cmd .= "-m $opt_m -x $opt_x -d $dir_all_bit\n";
        system($cmd);
    }

    if(exists $step_hash{6})
    {
        my $input = "$count_dir/filt.dump.gz";
        my $histo = "$count_dir/k$opt_k.histo";

        my $peak = getPeak($histo);
        die "can't find the peak\n" if($peak<0);

        my $min_occ = $opt_m;
        $min_occ = int(0.5*$peak) if($peak>0);
        $min_occ = $opt_m if($min_occ<$opt_m);

        my $cmd = "perl $0 bit -k $opt_k -i $input ";
        $cmd .= "-m $min_occ -x $opt_x -d $dir_peak_bit\n";
        system($cmd);
    }

    if(exists $step_hash{7})
    {
        mkdir($dir_analyse);
        chdir($dir_analyse);

        my($node, $total) = getInfo("$count_dir/k$opt_k.stats");

        my $cmd = "cd $dir_analyse \n";
        $cmd .= "$FORMAT $count_dir/k$opt_k.histo kmer.xls $node\n";
        $cmd .= "$ANALYSE kmer.xls $opt_k $total > result.txt \n";

        system($cmd);
    }
}
# |
# ==============================================================================



# Get the node number and total kmer
# ==============================================================================
# |
sub getInfo{
    my $file = shift;
    my $in_hdl = myOpen($file);

    my($node, $total) = (0, 0);
    while (<$in_hdl>)
    {
        chomp;
        my($key, $value) = split;
        $node  = $value if($key eq "Distinct:");
        $total = $value if($key eq "Total:");
    }
    close $in_hdl;

    return($node, $total);
}
# |
# ==============================================================================


# Get The Peak from product（乘积曲线）
# ==============================================================================
# |
sub getPeak
{
    my($peak, $max_time) = (-1, -1);

    my ($file, $min, $max) = @_;
    $min ||= 3;
    $max ||= 1000;

    # get the peak
    my @array;
    my $in_hdl = myOpen($file);
    while (<$in_hdl>)
    {
        chomp;
        my($occ, $num) = split;
        next if($occ<$min);
        last if($occ>$max);
        my $time = $occ*$num;
        shift @array if(@array==3);
        push(@array, "$occ $time");
        if(@array==3)
        {
            my($occ0, $time0) = split(/\s+/,$array[0]);
            my($occ1, $time1) = split(/\s+/,$array[1]);
            my($occ2, $time2) = split(/\s+/,$array[2]);
            if($time1>$time0 && $time1>$time2 && $time1>$max_time)
            {
                $peak = $occ1;
                $max_time = $time1;
            }
        }
    }
    close $in_hdl;

    return $peak;
}
# |
# ==============================================================================


# Check parameters
# ==============================================================================
# |
sub checkParam
{
    if ($help || @ARGV != 1)
    {
        die `pod2text $0`;
    }
}
# |
# ==============================================================================

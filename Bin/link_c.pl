use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path getcwd);
use List::Util qw(max min);

BEGIN {
    push (@INC,"$Bin");
}
use SeqIO;

die "perl $0 <contig.fa> <pb.fa> <convt.best> <outfile> <simi:0.001>\n" if(@ARGV==0);


my($contig, $pb, $align, $outfile, $SIMI) = @ARGV;
$SIMI ||= 0.001;

my($in_hdl, $ou_hdl,$lg_hdl);




# get the best pb for each gap
# ==============================================================================
# |
my %pb_hash;
my %link_hash;
$in_hdl = myOpen($align);
while (<$in_hdl>) {
    chomp;
    my($scf_id, $index1, $index2,
        $q_id, $t_id, $t_id2, $strand,
        $gap_scf, $gap_pb, $anchor1, $anchor2, 
        $q_len, $t_len, $t_len2) = split;

    $pb_hash{$q_id} = 1;

    $link_hash{"$scf_id $index1"}{$index2} = "$q_id\t$strand\t$anchor1\t$anchor2";
}
close $in_hdl;
# |
# ==============================================================================




# store the pb reads 
# ==============================================================================
# |
debug("loading pabio reads...");

my($id, $seq);
$in_hdl = myOpen($pb);
while (getSeq($in_hdl,\$id,\$seq)!=-1)
{
    next if(!exists $pb_hash{$id});
    $pb_hash{$id} = $seq;
}
close $in_hdl;
# |
# ==============================================================================



# load and link contigs
# ==============================================================================
# |
debug("close gap...");

my @array = ("");
my $pre_scf = "";

$in_hdl = myOpen($contig);
$ou_hdl = myOpen(">$outfile");
$lg_hdl = myOpen(">$outfile.log");

while (getSeq($in_hdl,\$id,\$seq)!=-1)
{
    my($scf, $index, $start, $len) = parseID($id);

    linkCtg() if($scf ne $pre_scf);

    push(@array, "$id\t$seq");
    $pre_scf = $scf;
}
linkCtg();
close $in_hdl;
close $ou_hdl;
close $lg_hdl;
# |
# ==============================================================================






################################################################################
####################    subroutine   ###########################################
################################################################################

# get the id info
# ==============================================================================
# |
sub parseID{
    my $id = shift;
    $id =~ /(\S+)_(\d+)\|start=(\d+)\|length=(\d+)/;
    return($1, $2, $3, $4);
}
# |
# ==============================================================================



# get the gap size
# ==============================================================================
# |
sub getGap{
    my($id1, $id2) = @_;
    my($scf1, $index1, $start1, $len1) = parseID($id1);
    my($scf2, $index2, $start2, $len2) = parseID($id2);
    my $gap = $start2-$start1-$len1;
    return $gap;
}
# |
# ==============================================================================



# close gap
# ==============================================================================
# |
sub linkCtg
{
    return if(@array<=1);

    my($id, $seq) = split(/\s+/, $array[1]);
    my $scf = $seq;
    
    for(my $i=1; $i<@array; $i++)
    {
        my($id1, $seq1)  = split(/\s+/, $array[$i]);
        my $len1 = length($seq1);

        # close gap
        my $is_close = 0;
        my $key = "$pre_scf $i";
        if(exists $link_hash{$key})
        {
            # get the next sequence
            my @idxs = keys %{$link_hash{$key}};
            @idxs = sort {return $a<=>$b} @idxs;
            my $next_idx = $idxs[0];
            my ($id2, $next_seq) = split(/\s+/, $array[$next_idx]);

            my $len2 = length($next_seq);
            my $gap = getGap($id1, $id2);
            
            # the pb read
            my ($pb_id, $strand, $anchor1, $anchor2) =
                                    split(/\s+/, $link_hash{$key}{$next_idx});
            my $pb_seq = $pb_hash{$pb_id};
            if($strand eq "-"){
                $pb_seq = revCom($pb_seq);
            }

            # connect sequences
            my($pb_part, $merge);
            my $ovl = 400;
            $pb_part = cutSeq($pb_seq, $anchor1, $anchor2, $len1, $len2, $ovl);
            $merge = con_seq($scf, $pb_part, $next_seq, $ovl); 
           
            if($merge ne "")
            {
                my $gap2 = length($merge)-length($scf)-length($next_seq);
                print $lg_hdl "connect $id1 $id2 by $pb_id $gap2 $gap\n";

                $scf = $merge;
                $is_close = 1;
                $i = $next_idx-1;
            }else{
                last if($i==$#array);
                my ($id2, $next_seq) = split(/\s+/, $array[$i+1]);
                my $gap = getGap($id1, $id2);
                $scf .= "N"x$gap;
                $scf .= $next_seq;
                print $lg_hdl "cannot connect $id1 $id2 by $pb_id $gap\n";
            }
        }
    }

    print $ou_hdl ">$pre_scf\n$scf\n";
    @array=("");
}
# |
# ==============================================================================



# cut the sequence by overlap
# ==============================================================================
# |
sub cutSeq{
    my($pb_seq, $anchor1, $anchor2, $t1_len, $t2_len, $ovl) = @_;

    my($q_pos1, $t1_pos) = split(/,/, $anchor1);
    my($q_pos2, $t2_pos) = split(/,/, $anchor2);

    #                 t1_pos t1_len       t2_pos
    #   ---------------|------|      --|-------
    #             -----|------|--------|----
    #               q_pos1           q_pos2
    my $pb_len = length($pb_seq);

    $q_pos1 += $t1_len-$t1_pos-$ovl;
    $q_pos1  = 0 if($q_pos1<0);

    $q_pos2 -= $t2_pos-$ovl;
    $q_pos2  = $pb_len-1 if($q_pos2>=$pb_len);

    my $len  = $q_pos2-$q_pos1;
    my $sub_seq = substr($pb_seq, $q_pos1, $len);
    return $sub_seq;
}
# |
# ==============================================================================


# 连接3条序列
# ==============================================================================
# |
sub con_seq{
    my($left, $mid, $right, $ovl) = @_;
    my($seq1, $simi1) = &link($left, $mid, $ovl);

    return "" if($seq1 eq "");

    my($seq2, $simi2) = &link($seq1, $right, $ovl);
    return $seq2;
}
# |
# ==============================================================================


# 连接2条相邻的序列
# ==============================================================================
# |
sub link{
	my($seq1, $seq2, $ovl) = @_;
	my $len1 = length($seq1);
	my $len2 = length($seq2);
    $ovl = min($len1, $len2, $ovl);

	my $part1 = substr($seq1, -$ovl, $ovl);
    my $part2 = substr($seq2, 0, $ovl);
	my($simi, $anchor) = align($part1, $part2, 5);

    # die "$simi\n" if($simi<$SIMI);

	return ("", $simi) if($simi<$SIMI);
	return ("", $simi) if($anchor eq "-1 -1");

	my($pos1, $pos2) = split(/\s+/, $anchor);

	my $merge = substr($seq1,0,$len1-$ovl+$pos1);
	$merge .= substr($seq2, $pos2);

	return ($merge, $simi);
}
# |
# ==============================================================================






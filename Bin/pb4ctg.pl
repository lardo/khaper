use strict;
use warnings;
use FindBin qw($Bin);
use List::Util qw(max min);

BEGIN {
    push (@INC,"$Bin");
}
use SeqIO;

die "perl $0 <filt.al>  <out.pos>\n" if(@ARGV==0);

my($align, $outfile) = @ARGV;

my($in_hdl, $ou_hdl);

# Record the short pacbio reads
# ==============================================================================
# |
debug("record short reads");

my %len_hash;
my %short_hash;
my @align_arr = ();

$in_hdl = myOpen($align);
while (<$in_hdl>)
{
	chomp;
	my ($contig, $pacbio, $con_len, $pb_len, $strand,
				$overlap, $con_pos,  $pb_pos, $score) = split;

	next if(exists $short_hash{$pacbio});

	if($pb_len==$overlap)
	{
		$len_hash{$contig} = $con_len;
		$short_hash{$pacbio}=2;
		next;
	}
}
close $in_hdl;
# |
# ==============================================================================



# Record pacbio related contig and center position
# ==============================================================================
# |
debug("parsing alignment");

my $pre_pb = "";
my @contig_arr = ();
$in_hdl = myOpen($align);
$ou_hdl = myOpen(">$outfile");
while (<$in_hdl>) {
	chomp;
	my ($contig, $pacbio, $con_len, $pb_len, $strand,
				$overlap, $con_pos,  $pb_pos, $score) = split;

	next if(exists $short_hash{$pacbio});
	recordCtg() if($pre_pb ne $pacbio);

	my $cent_pos=$pb_pos+int($con_len/2-$con_pos);
	push(@contig_arr, "$contig,$strand,$cent_pos");

	# 输出的文件中包含这条contig
	delete $len_hash{$contig};

	$pre_pb = $pacbio;
}
recordCtg();
outSuperCtg();
close $in_hdl;
close $ou_hdl;
# |
# ==============================================================================


# output the contig position
# ==============================================================================
# |
sub recordCtg{
	return if($pre_pb eq "");
	return if(@contig_arr==0);
	@contig_arr = sort by_pos(@contig_arr);

	print $ou_hdl ">$pre_pb\n@contig_arr\n";
	@contig_arr = ();
}

sub by_pos{
	my($ctg1,$strand1,$pos1) = split(/,/, $a);
	my($ctg2,$strand2,$pos2) = split(/,/, $b);
	return $pos1<=>$pos2;
}

sub outSuperCtg{
	foreach my $contig(keys %len_hash){
		print $ou_hdl ">$contig\n$contig,+,0\n";
	}
}
# |
# ==============================================================================

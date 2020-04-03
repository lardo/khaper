use strict;
use warnings;
use FindBin qw($Bin);
BEGIN {
    push (@INC,"$Bin");
}
use SeqIO;

die "perl $0 <filt.al> <out.al>\n" if(@ARGV==0);

my($infile, $outfile) = @ARGV;

my @array;
my $pre_pb = "";

open FI, $infile or die $!;
open FO, ">$outfile" or die $!;
while (<FI>) 
{
	my($pb, $contig, $pb_len, $contig_len,
		$strand, $ovl, $pb_pos, $contig_pos, $score) = split;

    next if($pb_len==$ovl);

	outConnect() if($pb ne $pre_pb);
	if($strand eq "-")
	{
		$pb_pos = $pb_len-$pb_pos;
		$contig_pos = $contig_len-$contig_pos;
	}

	my $pb_cent = $pb_pos+($contig_len/2-$contig_pos);
	push(@array, "$pb_cent $contig $contig_len $strand");
	
	$pre_pb = $pb;
}
outConnect();
close FI;
close FO;

sub by_pos{
	my $pos1 = (split(/\s+/,$a))[0];
	my $pos2 = (split(/\s+/,$b))[0];
	return $pos1<=>$pos2;
}

sub outConnect
{
	if(@array>1){
		@array = sort by_pos @array;
		for(my $i=0; $i<$#array; $i++)
		{
			my($pb_cent1, $contig1, $contig1_len, $strand1) = split(/\s+/, $array[$i]);
			for(my $j=$i+1; $j<@array; $j++)
			{
				my($pb_cent2, $contig2, $contig2_len, $strand2) = split(/\s+/, $array[$j]);
				if($contig1 eq $contig2) # meet the barcode
				{	
					$i=$j;
					last;
					next;
				}
				my $dist = ($pb_cent2-$pb_cent1)-(0.5*$contig1_len+0.5*$contig2_len);

				print FO "$pre_pb $contig1,$strand1,$contig1_len $contig2,$strand2,$contig2_len $dist\n";
			}
		}
	}
	@array = ();
}

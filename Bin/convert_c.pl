use strict;
use warnings;

die "perl $0 <filt.al> <out.al>\n" if(@ARGV==0);

my($infile, $outfile) = @ARGV;

# Convert format
# ==============================================================================
# |
my @array = ();
open FI, $infile or die $!;
open FO,">$outfile" or die $!;
my($pre_pb, $pre_scf) = ("","");
while (<FI>)
{
	chomp;
	my($q_id, $t_id, $q_len, $t_len, $strand, 
					$ovl, $q_pos, $t_pos, $block, $scf_id, $idx) = split;

	if($q_id ne $pre_pb or $scf_id ne $pre_scf)
	{
		parseLink();
	}

	push(@array, $_);

	$pre_pb = $q_id;
	$pre_scf = $scf_id;
}
parseLink();
close FI;
close FO;
`sort -k 1,1 -k 2,2n -k 3,3n -k 9,9n -S 3G -o $outfile $outfile`;
# |
# ==============================================================================


# Get the best linkage
# ==============================================================================
# |
$pre_scf = "";
my($pre_idx1, $pre_idx2) = (-1, -1);
open FI, "$outfile" or die $!;
open FO, ">$outfile.best" or die $!;
while (<FI>)
{
	chomp;
	my($scf_id, $index1, $index2, $others) = split;

	if($scf_id ne $pre_scf or 
		$index1 ne $pre_idx1 or 
		$index2 ne $pre_idx2)
	{
		getBest();
	}

	push(@array, $_);

	$pre_scf = $scf_id;
	$pre_idx1 = $index1;
	$pre_idx2 = $index2;
}
getBest();
close FI;
close FO;
# |
# ==============================================================================



# Convert format
# ==============================================================================
# |
sub parseLink
{
	if(@array<2)
	{
		@array = ();
		return;
	}

	for(my $i=0; $i<@array; $i++)
	{
		my($q_id, $t_id, $q_len, $t_len, $strand, $ovl,
			$q_pos, $t_pos, $block, $scf_id, $index) = split(/\s+/, $array[$i]);

		$t_id =~ /(\S+)_(\d+)\|start=(\d+)\|length=(\d+)/;
		my($start, $len) = ($3, $4);
		
		# 进行连接
		for(my $j=$i+1; $j<@array; $j++)
		{
			my($q_id2, $t_id2, $q_len2, 
				$t_len2, $strand2, $ovl2,
				$q_pos2, $t_pos2, $block2, 
				$scf_id2, $index2) = split(/\s+/, $array[$j]);;

			next if($t_id eq $t_id2);
			next if($strand ne $strand2);

			# gap observed by scaffold
			$t_id2 =~ /(\S+)_(\d+)\|start=(\d+)\|length=(\d+)/;
			my($start2, $len2) = ($3, $4);
			my $gap_scf = $start2-$start-$len;

			# gap observed by pb
			#        t_pos1          t_pos2
			# --------|-----        --|----------
			#      ---|---------------|---
			#        q_pos1          q_pos2
			my $tail1 = $t_len-$t_pos;
			my $head1 = $t_pos2;
			my $gap_pb = $q_pos2-$q_pos-$head1-$tail1;
			next if($gap_pb<-1000);

			die "$array[$i]\n$array[$j]\n" if($scf_id ne $scf_id2);

			print FO "$scf_id\t$index\t$index2\t";
			print FO "$q_id\t$t_id\t$t_id2\t$strand\t$gap_scf\t$gap_pb\t";
			print FO "$q_pos,$t_pos\t$q_pos2,$t_pos2\t$q_len\t$t_len\t$t_len2\n";
		}
	}
	@array=();
}
# |
# ==============================================================================


sub getBest{
	return if(@array==0);
	my $mid = int(@array/2);
	print FO "$array[$mid]\n";
	@array = ();
}

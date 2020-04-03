use strict;
use warnings;

die "perl $0 <align.al> <outfile> <min_ovl:300> <score:0.7>\n" if(@ARGV==0);

my($in_file, $out_file, $min_ovl, $score) = @ARGV;
$min_ovl ||= 300;
$score ||= 0.7;


open FI, $in_file or die $!;
open FO,">$out_file" or die $!;
while (<FI>) {
	chomp;
	my($q_id, $t_id, $q_len, $t_len,
		$strand, $ovl, $match, $block,
		$q_pos, $t_pos, $uni
		) = split;

    next if($q_id eq $t_id);
	next if($ovl<$min_ovl);
	next if($block<$score);
	next if($match<$score);

	print FO "$q_id\t$t_id\t$q_len\t$t_len\t$strand\t$ovl\t$q_pos\t$t_pos\t$block\n";
}
close FI;
close FO;

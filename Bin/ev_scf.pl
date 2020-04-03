use strict;
use warnings;

die "perl $0 <test.al>\n" if(@ARGV==0);

my($align) = @ARGV;


open FI, $align or die $!;
open FO, ">temp.al";
while(<FI>){
    chomp;
    my($ctg, $chr, $ctg_len, $chr_len,
        $strand, $ovl, $match, $block,
        $q_pos, $t_pos, $uni_kmer) = split;

    next if($match<0.9);

    my($id, $start, $len) = split(/\|/, $ctg);
    $id =~ /(\w+)_(\d+)$/;
    my($scf, $idx) = ($1, $2);

    print FO "$scf\t$idx\t$chr\t$t_pos\n";
}
close FI;
close FO;

`sort -k 1,1 -k 2,2n -o temp.al temp.al`;

open FI, "temp.al";
open FO, ">temp2.al";
my($pre_scf, $pre_idx, $pre_chr, $pre_pos);
while (<FI>) {
    chomp;
    my $stat = 1;
    my($scf, $idx, $chr, $t_pos) = split;
    if(defined $pre_scf)
    {
        if($scf eq $pre_scf)
        {
            my $dist1 = $idx-$pre_idx;
            my $dist2 = $t_pos-$pre_pos;

            if($pre_chr ne $chr){
                $stat = -1;
            }
        }
    }
    print FO "$_\t$stat\n";

    ($pre_scf, $pre_idx, $pre_chr, $pre_pos) = ($scf, $idx, $chr, $t_pos);
}
close FI;
close FO;

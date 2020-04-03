use strict;
use warnings;

die "perl $0 <test.al>\n" if(@ARGV==0);

my($align) = @ARGV;

my $total = 0;
my $match = 0;
open FI, $align or die $!;
while(<FI>){
    chomp;
    my($id1, $id2) = (split)[0,1];
    
    my($s1, $idx1) = split(/_/,$id1);
    my($s2, $idx2) = split(/_/,$id2);
    
    $total++;
    $match++ if($s1 eq $s2);
}
close FI;

my $perc = $match/$total;

print "total: $total\nmatch: $match\npercent: $perc\n";

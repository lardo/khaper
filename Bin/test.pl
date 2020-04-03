use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path getcwd);

BEGIN {
    push (@INC,"$Bin");
}
use SeqIO;

use Inline C => Config => INC => '-I$Bin/../Source/DNA2';
use Inline C => Config => LIBS => '-LBin/../Source/DNA2 -lSeq';
use Inline 'C';

my $seq1 = "CCAACCAGTAGTTCCACCTCTGGTTCTTCTGAGAGCAAAACGAGTTCGGCTAGTTCTTCCTCTTCTTCCTCTTCTATCTCTTCTGAATCACCAAA";
my $seq2 = "TTCTGAGAGCAAAACGAGTTCGGCTAGTTCTTCCTCTTCTTCCTCTTCTATCTCTTCTGAATCACCAAA";
my $result = sayhello($seq1, $seq2);

print "$result\n";

__END__
__C__

#include "SeqTool.h"

int sayhello(char* seq1, char* seq2){
    return 1000;
}

int align(char* seq1, char* seq2){
    SeqTool tool;
    cerr<<seq1<<endl;
    Conn_info info = tool.alignment(seq1, seq2, 1);
    return 99;
}

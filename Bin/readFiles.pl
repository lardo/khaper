use strict;
use warnings;

die "perl $0 <read.lst> " if(@ARGV==0);

my($read_lst) = @ARGV;

my $lst_hdl = myOpen($read_lst);
while (my $file = <$lst_hdl>) {
	chomp $file;
	my $read_hdl = myOpen($file);
	while (<$read_hdl>) {
	    print $_;
    }
	close $read_hdl;
}
close $lst_hdl;



#########################################################################
sub myOpen{
    my $file = shift;
    my $handle;
    if($file !~ /\.gz$/){
        open $handle,"$file" or die "can't open $file";
    }else{
        open $handle,"gzip -dc $file|" or die "can't open $file";
    }
    return $handle;
}

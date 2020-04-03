#!/usr/bin/perl
package Qsub;
require Exporter;

use warnings;
use strict;
use File::Basename;
use Cwd qw(abs_path getcwd);

our @version = 1.0;
our @ISA     = qw(Exporter);
our @EXPORT  = qw(
    				qsub
                    qsub2
);

# software path
# ==============================================================================
#|
# forgive me for using abspath of LIB that the findBin don't work for perl module
my $LIB = "/nfs1/public1/Software/pymonitor/";

my $STAT = "$LIB/monitor stat ";
my $QSUB = "$LIB/monitor taskmonitor ";
my $UPDATE = "$LIB/monitor updateproject -p ";
my $RMPRO = "$LIB/monitor removeproject -m 1 -p ";
my $CLEAN = "$LIB/monitor removeproject -d";
#|
# ==============================================================================

my $pro_name;

# handle terminal signal
# ==============================================================================
#|
$SIG{INT}  = \&signal_handler;
$SIG{TERM} = \&signal_handler;
$SIG{KILL} = \&signal_handler;

sub signal_handler 
{
    my $pid = $$;
    my @group = getPGID();
    
    if(defined $pro_name)
    {
    	`$RMPRO $pro_name`;
        debug("$RMPRO $pro_name\n");
    }

    foreach(@group)
    {
    	next if($_==$pid);
    	`kill $_` if(isAlive($_));
    }
    `kill $pid`;
    exit(0);
}

sub isAlive
{
	my $id = shift;
	my $cmd = "ps -jf|awk '\$2==$id'";
	my $info = `$cmd`;
	chomp $info;
	return 1 if(length($info)>5);
	return 0;
}

sub getPGID
{
	my @pid_arr;
	
	# get the pgid
	my $cmd = "ps -jf |awk '\$2==$$'";
	my $info = `$cmd`;
	chomp $info;
	my($UID, $PID, $PPID, $PGID, $SID) = split(/\s+/, $info);
	
	# get the pids 
	$cmd = "ps -jf |awk '\$4==$PGID'";
	$info = `$cmd`;

	my @array = split(/\n/, $info);
	foreach(@array)
	{
		next if(length($_)<3);
		my $pid = (split)[1];
		push(@pid_arr, $pid);
	}

	return @pid_arr;
}
#|
# ==============================================================================


# output msg
# ==============================================================================
#|
sub debug
{
    my $msg = shift;
    print STDERR "$msg\n";
}
#|
# ==============================================================================


# Open file in different format
# ==============================================================================
#|
sub myOpen
{
    my $file = shift;
    my $handle;
    if($file =~ /\.gz$/ && $file !~ /\>/){
        open $handle,"gzip -dc $file|" or die "can't open $file";
    }else{
         open $handle,"$file" or die "can't open $file";
    }
    return $handle;
}
#|
# ==============================================================================


# 输出SHELL文件
# ==============================================================================
# |
sub outShell{
  my($cmd, $outfile, $app) = @_;
  $app ||= 0;

  my $out_hdl;
  $out_hdl = myOpen(">$outfile")  if($app==0);
  $out_hdl = myOpen(">>$outfile") if($app==1);

  print $out_hdl "$cmd\n";
  close $out_hdl;
}
# |
# ==============================================================================


# qsub
# file: script $cmd
# ==============================================================================
# |
sub qsub
{
	my($file, $outdir, $memory, $thread,
				$queue, $project, $pro_name_,  $max_job) = @_;

    $pro_name = "${pro_name_}_$$";

	mkdir($outdir) unless (-e $outdir);
    my %script_hash;
	open FO, ">$outdir/config$$.txt";
	my $in_hdl = myOpen($file);
	while (<$in_hdl>) 
	{
		chomp;
		my($script, @other) = split;
        next if(!defined $script);

        if(exists $script_hash{$script}){
            die "the script name occured before: $script\n";
        }

        $script_hash{$script} = 1;
        
		my $cmd = "@other";
		$cmd .= " && date && echo done";
		
		$script = abs_path("$outdir/$script") if($script !~ /^\//);
		outShell($cmd, $script);
		print FO "$script:$memory\n";
	}
	close $in_hdl;
	close FO;

    my $cmd = "$QSUB -i $outdir/config$$.txt  ";
    $cmd .= "-f 3 -s done -n $max_job -t $thread ";
    $cmd .= "-p $pro_name -q $queue -P $project ";
    
    debug("$cmd\n");    

    system($cmd);
    
    while (checkJobs($pro_name)==0)
    {
    	sleep(100);
	}
}
# |
# ==============================================================================


# qsub
# file: script $cmd
# ==============================================================================
# |
sub qsub2{
    my($cfg_file, $thread,
        $queue, $project, $pro_name_,  $max_job) = @_;
    $pro_name = "${pro_name_}_$$";

    my $cmd = "$QSUB -i $cfg_file  ";
    $cmd .= "-f 3 -s done -n $max_job -t $thread ";
    $cmd .= "-p $pro_name -q $queue -P $project";

    my $info = system($cmd);

    while (checkJobs($pro_name)==0)
    {
        sleep(100);
    }
}
# |
# ==============================================================================


sub checkJobs{
    my($pro_name) = @_;

    `$UPDATE $pro_name`;
    my $lines = `$STAT -p $pro_name`;
    my $info = (split(/\n/,$lines))[-1];

    my @array = split(/\s+/,$info);

    return 0 if(@array<7);
    my($name,$done,$total) = @array[0,7,9];

    return 0 if($done !~ /^\d+/);

    if($done==$total)
    {
        `$RMPRO $pro_name`;
        return 1 ;
    }
    return 0;
}

1;

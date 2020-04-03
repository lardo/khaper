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
use Qsub;

# Get parameters
# ==============================================================================
# |
die "perl $0 <input.cfg>\n" if(@ARGV==0);
my($cfg_file) = @ARGV;

#contig
my $contig;

# pb file
my($genome_size, $pb_lst);
my($filt_len);
# kmer
my($kmer_file, $kmer_size);
# sge
my($queue, $project, $pro_name, $job_num);
# thread
my $thread;
# recursive
my $recursive = 1;

parseCfg($cfg_file);
# |
# ==============================================================================


# Global Variables
# ==============================================================================
# |
my $SLEEP = 100;
# |
# ==============================================================================


# Software Paths
# ==============================================================================
# |
my $MERGE = "perl $Bin/splitFQ.pl ";
my $CTG   = "perl $Bin/get_scaftig.pl ";
my $CUT = "perl $Bin/cutFA.pl";

my $ALIGN = "perl $Bin/Align.pl";
my $FILTER = "perl $Bin/filtAl_m.pl";
my $CONVT = "perl $Bin/pb4ctg.pl";

my $GRAPH = "perl $Bin/buildGraph_m.pl";
my $GR2PB = "perl $Bin/graph2pb.pl";
my $LINK =  "perl $Bin/link_m.pl";

my $CORRECT = "perl $Bin/Correct.pl";
# |
# ==============================================================================

# 目录结构
# ==============================================================================
# |
mkdir("Data");
mkdir("Shell");
mkdir("Assemble");
mkdir("Consensus");

my $dir_input = abs_path("Data");
my $dir_shell = abs_path("Shell");
my $dir_align = abs_path("Align");
my $dir_assem = abs_path("Assemble");
my $dir_cor = abs_path("Consensus");
# |
# ==============================================================================


# prepare data
# ==============================================================================
# |
debug("prepare data");

my $pb_data = "$dir_input/pb.fasta";
my $contig_ = "$dir_input/contig.fasta";

my $cmd = "";
$cmd .= "$MERGE $pb_lst $pb_data $filt_len\n";
$cmd .= "$CTG $contig | $CUT /dev/stdin $contig_ -1\n";
system($cmd);

$contig = $contig_;
# |
# ==============================================================================

for(my $i=1; $i<=$recursive; $i++)
{
    # 将pacbio文件比对回contig中
    # ==========================================================================
    # |
    debug("do alignment");

    # 清除目录
    `rm -r Align` if(-e "Align");

    # 切分大小
    my $divide_size = int($genome_size/10**6)+1;
    $divide_size=400 if($divide_size>400);

    # 比对脚本
    $cmd  = "$ALIGN -k $kmer_size -g $kmer_file -u 3 ";
    $cmd .= "-s $divide_size -s2 2000 -n 1  -t $thread ";
    $cmd .= "-p $project -q $queue -b $pro_name ";
    $cmd .= "-d $dir_align  $contig $pb_data";
    debug($cmd);
    system($cmd);
    # |
    # ==========================================================================


    # 数据聚类
    # ==========================================================================
    # |
    debug("do assembly");

    # 进行组装
    $cmd  = "cd $dir_assem && ";
    # filter low quality alignment
    $cmd .= "$FILTER $dir_align/align.al filter.al 300 0.7 && ";
    # sort by pacbio
    $cmd .= "sort -k 2,2 -S 1g -o filter.al filter.al && ";
    # convert the pacbio format
    $cmd .= "$CONVT filter.al contig.pos && ";
    # build graph
    $cmd .= "$GRAPH contig.pos graph.info && ";
    $cmd .= "$GR2PB graph.info contig.pos graph_pb && ";
    # link sequence
    $cmd .= "$LINK graph_pb.lk backbone.fasta $contig $pb_data ";
    # |
    # ==========================================================================


    # 提交任务
    # ==========================================================================
    # |
    my $cmd_file = "$dir_shell/ASSEMBLE.SH";
    debug($cmd);
    outShell2("assemble.$i.sh\t$cmd", $cmd_file);
    qsub($cmd_file, $dir_shell, "5G", 1,
                    $queue, $project, $pro_name, 50);
    # |
    # ==========================================================================


    # 整理结果
    # ==========================================================================
    # |
    `cp $dir_assem/backbone.fasta $dir_input/assemble.$i.fa`;
    $contig = "$dir_input/assemble.$i.fa";
    # |
    # ==========================================================================
}


# 进行纠错
# ==============================================================================
# |
$cmd  = "$CORRECT -c 1 -g $kmer_file -k $kmer_size -u 3 -s 400 -s2 2000 ";
$cmd .= "-t $thread -m 1 -d $dir_cor -b $pro_name -q $queue -p $project ";
$cmd .= "-x -1 -j 1 -n 1 ";
$cmd .= "$contig $pb_data";
debug($cmd);
system($cmd);
# |
# ==============================================================================


# 从配置参数中取得参数
# ==============================================================================
# |
sub parseCfg{
    my $cfg_file = shift;
    my $cfg_hdl = myOpen($cfg_file);

    while (<$cfg_hdl>) {
        chomp;
        next if(/^#/);
        next if(!/^\[/);
        my($key, @values) = split;

        # contig
        $contig = $values[0] if($key eq "[contig]");
        # pb file
        $pb_lst = $values[0] if($key eq "[pb_lst]");
        $filt_len = $values[0] if($key eq "[filt_len]");
        $genome_size = $values[0] if($key eq "[genome_size]");
        # kmer
        ($kmer_file,$kmer_size) = @values if($key eq "[unique_kmer]");
        # qsub job
        $queue = $values[0] if($key eq "[queue]");
        $project = $values[0] if($key eq "[Project]");
        $pro_name = $values[0] if($key eq "[pro_name]");
        $job_num = $values[0] if($key eq "[max_job]");
        # thread
        $thread = $values[0] if($key eq "[thread]");
        # recursive
        $recursive = $values[0] if($key eq "[recursive]");
    }
    $genome_size *= 10**6;
    $kmer_file = abs_path($kmer_file);
    $thread ||= 8;
    $job_num ||= 50;
    $contig = abs_path($contig);

    close $cfg_hdl;
}
# |
# ==============================================================================

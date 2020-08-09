This is a software for remove dupplication for genome. 

**Note:**

* This software Need > 100G memory 

* This version on GitHub is only support reference <4000M

* Please contact us: tanger.zhang@gmail.com

  

**Usage:**

We firstly use >40x Illumina reads to build the kmer frequency table. Then use this kmer table to compress the reference.

1. install jellyfish

  ```
  conda install -c bioconda jellyfish
  ```

2. Prepare input files:

  ```
  Prepare:
  
  assemble.fasta 	# genemone assembly with dupplcated sequences.
  PE300_1.fq.gz		# read1
  PE300_2.fq.gz		# read2
  ```

3. Build the kmer frequency table:

  ```
  ls *.gz > fq.lst
  perl Bin/Graph.pl pipe -i fq.lst -m 2 -k 15 -s 1,3 -d Kmer_15
  
  #result:
  kmer bit file: Kmer_15/02.Uinque_bit/kmer_15.bit
  ```

  **Note:** 

  ​	a. k=15 is suitable for genome with size <100M.

  ​	b. k=17 is suitable for genome with size <10G.

  ​	c. This version is only support k<=17.

4. Compress the assembly file

  ````
  # compress the genome
  
  # Usage:
   perl remDup.pl <genome.fa> <outdir> <cutoff:0.7>
  
       Options:
              --ref   <str> The ref genome to build kbit
            --kbit  <str> The unique kmer file
              --kmer  <int> the kmer size [15]
            --sort  <int> sort seq by length [1]
  
Description
       This script is to remove dupplcation seq
  
  # Demo
  perl Bin/remDup.pl  --kbit Kmer_15/02.Uinque_bit/kmer_15.bit --kmer 15 assemble.fasta Compress 0.3
  
  # result:
  compress file: Compress/trinity.single.fasta.gz
  ````
  
  **Note:**
  
  ​	a. If the compress file is larger than estimated genome size, turn down the **"<cutoff>"** value
  
  ​	b. If the compress file is small than estimated genome size, turn up the **"<cutoff>"** value


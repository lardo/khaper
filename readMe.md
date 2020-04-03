This is a software for remove dupplication for genome. 

**Note:**

* This software Need > 100G memory 

* This version on GitHub is only support reference <100M

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
  ```

4. Compress the assembly file

   ````
   # compress the genome
   echo "clean assemble.fasta" > file.lst
perl Bin/Compress.pl compress -i file.lst -g Kmer_15/02.Uinque_bit/kmer_15.bit -k 15 -m 3 -t 0.3 -n 1
   ````
   


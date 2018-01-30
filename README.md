# Asite_LP_method
The following text describes the process to run the Linear Programming method to identify the A-site positions within ribosome protected fragments. 

The analysis is carried out in two steps. 

Step 1: An alignment file is processed to create read count files according to fragment length. These files will be used in Step 2 to run the Linear Programming algorithm. Genes are selected based on filtering criteria specified by the user.

Step 2: Using the read count files created in Step 1, Linear Programming algorithm is run and the output is a A-site offset table specific to fragment size and frame. 

### Step 1:  
There are two options to run Step 1 based on how the alignments have been made.

### Option A:  Alignment to genome
  This option is preferable for smaller genomes (e.g. E.coli, S. cerevisiae)  with simple gene structures with fewer introns and isoforms.

The input files for running Step 1 under option A are
-	SAM alignment file. A BAM file can also be given as input but it requires samtools to be installed for it to be converted to a SAM file.
-	Genome fasta file.
-	Annotation file in the following format
 
|Gene name	|Chromosome|	Strand|	CDS Length|	Number of CDS regions|	CDS Start 1| CDS End1  |CDS Start 2	|CDS End 2|
|:------:|:---------:|:------:|:---------:|:------:|:---------:|:------:|:---------:|:------:|
|	YAL008W |	chrI |	+ |	597	| 1 |	136914 | 137510 | | |
|	YBR255W |chrII |	+	| 2085 | 1 |	724456 | 726540 | | |
|	YJL142C |	chrX |	-	| 393	| 1 |	147819 | 148211 | | |
|	YFR045W |	chrVI |	+	| 930	| 2 |	241998 | 242009 | 242082 | 242999 |
|	YBL087C |	chrII|	- |	414 |	2 |	59822 |	60193 |	60698 |	60739 |

- Alternatively, a GFF file can also be given as input for annotations but the script is optimized for processing only GFF files for E. coli and S. cerevisiae. Ensure that the GFF file of your organism is in the same format as GFF files for the above organisms.

An example for running Step 1 under option A is shown below.

```
python src/asite_lp_preprocess.py -m 20 -x 35 -g data_files/sacCer3/cds_info.tab -e data_files/sacCer3/sacCer3_R64-2-1_genome.fa -i  input.sam
```

Additional input parameters are shown below
-	-a: Alignment mode: 1) genome (Default) 2) transcriptome 
-	-m: Minimum fragment length (Deafult: 20)
-	-x: Maximum fragment length (Default: 35)
-	-3:  If Yes, the script will quantify the reads from the 3’ end rather than 5’ end (Default).
-	-v: Number of nucleotides beyond the CDS region of a gene which are to be avoided to overlap with another gene. (Default:0 nt)
-	-o: output directory to write the output files. (Default: Current working directory)

The following output files will be written in the new directory called ‘output’
-	Read count files for each fragment size between minimum and maximum fragment length
-	Multiple mapped read counts. This file has to be used as input in Step 2
-	CDS table. This file will be used for creating an A-site table and run downstream analyses.

### Option B: Alignment to transcriptome
This option is preferable for complex genomes like mouse and human. The alignment file contains reads directly aligned to gene transcripts. The following input files are required to run Step 1 in the transcriptome mode.
- SAM alignment file with reads mapping to gene transcripts. A BAM file can also be given as input but it requires samtools to be installed for it to be converted to a SAM file.
- Fasta file containing the sequences of the transcripts.
- Annotation file in the following format

|Gene name	|Start index	|CDS Length|
|:------:|:---------:|:------:|
|uc008jxs.1	|79	|	1155|
|uc009fyo.2	|1048| 4023|	
|uc009ktx.2	|187|	1182|	
|uc007ehm.2	|332|	4173|	
|uc008qjj.2	|157|	1143|
 
 The start index is the index in the transcript where the CDS region starts upto a length specified in CDS Length column.

An example for running Step 1 under this option is shown below.

```
python src/asite_lp_preprocess.py -g data_files/Mus_musculus/mm10_start_index.tab -e data_files/Mus_musculus/mm10_transcript.fa -a transcriptome -m 20 -x 35 -i input.sam
```

Currently writing...

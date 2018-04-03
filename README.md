# Asite_LP_method
The following text describes the process to run the Linear Programming method to identify the A-site positions within ribosome protected fragments. 

The analysis is carried out in two steps. 

**Step 1:** An alignment file is processed to create read count files according to fragment length. These files will be used in Step 2 to run the Linear Programming algorithm. Genes are selected based on filtering criteria specified by the user.

**Step 2:** Using the read count files created in Step 1, Linear Programming algorithm is run and the output is a A-site offset table specific to fragment size and frame. 

### Step 1:  
There are two options to run Step 1 based on how the alignments have been made.

### Option A:  Alignment to genome
  This option is preferable for smaller genomes (e.g. *E.coli*, *S. cerevisiae*)  with simple gene structures with fewer introns and isoforms.

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
python src/asite_lp_preprocess.py -m 20 -x 35 -g data_files/sacCer3/CDS_info.tab -e data_files/sacCer3/sacCer3_R64-2-1_genome.fa -i  input.sam
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

There is option for generating the A-site read density profiles by providing a file containing offset table generated from result of LP algorithm. This is provided for those users who are using this script to parse their alignment file and directly generate the A-site profiles. The input offset file should contain offsets in rows of different fragment sizes and columns in frame 0,1 and 2. The first column should be fragment size.

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
python src/asite_lp_preprocess.py -g data_files/mm10/mm10_start_index.tab -e data_files/mm10/mm10_transcript.fa -a transcriptome -m 20 -x 35 -i input.sam
```

### Step 2:

The Read count files created under the directory `output/ ` will be used as input for Step2 where we will run the Linear Programming algorithm to determine fragment size and frame specific offsets for A-site within the ribosome protected fragments. This step is general to any organism and depends on the output created from Step 1.  
An example for running this step is shown below

```
time python src/asite_lp_run.py -i output/ -j output/Multiple_mapped_gene_read_counts_20_35.tab -m 20 -x 35
```

The parameters which can be given as an input to this step are
-	-i: Input directory containing the read count files (Usually output/ directory from Step1)
-	-j: Path to Multiple mapped read count files (Usually within the output/ directory from Step 1)
-	-m: Minimum fragment length
-	-x: Maximum fragment length
-	-o: output directory
-	-n: Threshold for average reads per codon for filtering genes(Default: 1 read per codon)
-	-t: Threshold for assigning a unique offset (Default: 70%)
-	-s: Threshold for secondary selection criteria (Deafult: 5). 5*R(1) < Average(R(2), R(3), R(4))
-	-3:  If given ‘Yes’, reads are mapped according to 3’ end rather than 5’ end (Default)

The output files created from Step 2 are

-	**Results_LP_algorithm.tab:**  This file contains the optimum offset table for A-site positions within ribosome protected fragments according to fragment size and frame. This file also contains  percentage of genes for unique offsets and top two percentages for combinations having ambiguous offsets.  The file also contains number of genes and number of reads mapped to each fragment size and frame.
-	 **Percentage_of_genes.tab:** This file contains the percentage of genes mapping to each offset between 0 and S which are multiple of 3. This information is written for each fragment length as well as for each frame within each fragment length.


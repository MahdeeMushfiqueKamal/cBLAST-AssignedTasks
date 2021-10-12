README

sgRNAcas9: a software package for designing CRISPR sgRNA and evaluating potential off-target cleavage sites

##########################################################################################

sgRNAcas9 is a software package, which contains eight scripts as follows: 
(1)sgRNAcas9.pl(main script), 
(2)combine_genome.pl,  #add at version 2.0.9
(3)format_genome.pl, 
(4)ot2gtf.pl,          #modified at version 2.0.8
(5)pot2gtf.pl
(6)check_sgRNA_seq.pl, 
(7)sgRPrimer.pl,
(8)extract_targetSeq.pl,

Copyright (C) 2014-2015, version 3.0.5 (all in one)
Software written by: Shengsong Xie
http://www.biootools.com

College of Veterinary Medicine, Huazhong Agricultural University
1st Shizishan Street, Hongshan, Wuhan 430070 China
http://www.hzau.edu.cn/en/home/

Institute of Biochemistry and Cell Biology, 
Shanghai Institutes for Biological Sciences, 
Chinese Academy of Sciences.
http://www.sibs.cas.cn

##########################################################################################
##########################################################################################

Anyone can use the source codes, documents or the excutable file of sgRNAcas9 free of charge 
for non-commercial use. For commercial use, please contact the author.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

If you use this program in your research, please cite:
Xie et al., sgRNAcas9: a software package for designing CRISPR sgRNA and evaluating 
potential off-target cleavage sites. 2014. PloS one.

Please send bug reports to: ssxieinfo@gmail.com or ssxie@mail.hzau.edu.cn
This README file covers the following topics:

1. Prerequest for sgRNAcas9
2. Prepare sgRNAcas9 input (fasta/genome) files
3. How to run sgRNAcas9
4. sgRNAcas9.pl output files explanation

##########################################################################################

1. Prerequest for sgRNAcas9

Perl (http://www.perl.org/, version: >= 5.8.8)
Perl 5 is a highly capable, feature-rich programming language with over 25 years of development. 
Perl 5 runs on over 100 platforms from portables to mainframes and is suitable for both rapid 
prototyping and large scale development projects. 

SeqMap(http://www-personal.umich.edu/~jianghui/SeqMap/)
SeqMap is a tool for mapping large amount of oligonucleotide to the genome. It is designed for 
finding all the places in a genome where an oligonucleotide could potentially come from. SeqMap 
can efficiently map as many as dozens of millions of short sequences to a genome of several 
billions of nucleotides. A typical mapping can be done in a few hours on an ordinary PC.

Latest update:
sgRNAcas9_3.0.5-All-in-One

For the sake of convenience, seqmap program has already been included in the sgRNAcas9 software.
(seqmap-1.0.12-linux, seqmap-1.0.12-linux-64, seqmap-1.0.12-mac, seqmap-1.0.12-mac-64,
seqmap-1.0.12-windows.exe)

When first using sgRNAcas9 in linux systerm,you need use chmod +x on the seqmap so linux will let 
you run it,

for example: chmod +x seqmap-1.0.12-linux-64;
             chmod +x seqmap-1.0.12-linux;

##########################################################################################

2. Prepare sgRNAcas9 input (fasta/genome) files

Example of input fasta file (hEMX1_example.txt):
>hEMX1_exon1
AGGTGAGCGGCGGCCAATGGGCGAGCGCGGGGCAGGTGCCCGCTAACTCGCGCCTCGCAGCGCTGGGCGGCCGGGGCTGGGCAGGGCAGTGCGGGGACAC
CGGGGGCTGGGGTCGGTCCCAGCGGGACTCCGAAAGGAGGGAGACGAGCTCAACCCTCGGGCCTTACTGGCAGCTCGCAGCCTAGCACGGAGCCCGCGCC
TGTGCGGGCGCCTGGAGCTGCCCGCTCCGCCGCAGCAGCCGCCGCGCCTGGCCGTACGCTGTGGCCGGACCCCGCGGTCGCTCGCTCACACACCCCTCGC
CGCTCCGCGCCTGGCTCGCCCGCGGGGGCCGAGCGCGAGCGGGCGGGCGGGGGAGGTGAGGGGTGCGGGCGGGTGTGCATGTGCCTGGCTGGGTGCACAC
CCCGCAAGGCGGCGGCGCCAGGACGCGGAGCGCTCCCCAGAGCCCGGCTGCCTCGCACAGCTCCCGCGGCTGCGACCATGTTCCAGCCCGCGGCCAAGCG
CGGCTTTACCATAGAGTCCTTGGTGGCCAAGGACGGCGGCACCGGCGGGGGCACTGGCGGCGGGGGCGCGGGCTCCCATCTCCTGGCGGCGGCCGCCTCC
GAGGAACCGCTCCGGCCCACGGCGCTCAACTACCCTCACCCCAGCGCGGCCGAGGCGGCCTTCGTGAGTGGCTTCCCTGCCGCGGCCGCCGCGGGCGCGG
GCCGCTCGCTCTACGGTGGGCCCGAGCTCGTGTTCCCCGAGGCCATGAACCACCCCGCGCTGACCGTGCATCCGGCGCACCAGCTGGGCGCCTCCCCGCT
GCAGCCCCCGCACTCCTTCTTCGGCGCCCAGCACCGGGACCCTCTCCATTTCTACCCCTGGGTCCTGCGGAACCGCTTCTTCGGCCACCGCTTCCAGG
>hEMX1_exon2
CCAGCGACGTGCCCCAGGACGGGCTGCTTCTGCACGGCCCCTTCGCACGCAAGCCCAAGCGGATCCGCACGGCCTTCTCGCCCTCGCAGCTGCTGCGGCT
GGAGCGCGCCTTCGAGAAGAACCACTACGTGGTGGGCGCCGAGCGGAAGCAGCTGGCCGGCAGTCTCAGCCTCTCCGAGACGCAG
>hEMX1_exon3
GTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAAGAAGGGCTCCCATCACA
TCAACCGGTGGCGCATTGCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAACCACAAACCCACGAGGGCAGAGT
GCTGCTTGCTGCTGGCCAGGCCCCTGCGTGGGCCCAAGCTGGACTCTGGCCACTCCCTGGCCAGGCTTTGGGGAGGCCTGGAGTCATGGCCCCACAGGGC
TTGAAGCCCGGGGCCGCCATTGACAGAGGGACAAGCAATGGGCTGGCTGAGGCCTGGGACCACTTGGCCTTCTCCTCGGAGAGCCTGCCTGCCTGGGCGG
GCCCGCCCGCCACCGCAGCCTCCCAGCTGCTCTCCGTGTCTCCAATCTCCCTTTTGTTTTGATGCATTTCTGTTTTAATTTATTTTCCAGGCACCACTGT
AGTTTAGTGATCCCCAGTGTCCCCCTTCCCTATGGGAATAATAAAAGTCTCTCTCTTAATGACACGGGCATCCAGCTCCAGCCCCAGAGCCTGGGGTGGT
AGATTCCGGCTCTGAGGGCCAGTGGGGGCTGGTAGAGCAAACGCGTTCAGGGCCTGGGAGCCTGGGGTGGGGTACTGGTGGAGGGGGTCAAGGGTAATTC
ATTAACTCCTCTCTTTTGTTGGGGGACCCTGGTCTCTACCTCCAGCTCCACAGCAGGAGAAACAGGCTAGACATAGGGAAGGGCCATCCTGTATCTTGAG
GGAGGACAGGCCCAGGTCTTTCTTAACGTATTGAGAGGTGGGAATCAGGCCCAGGTAGTTCAATGGGAGAGGGAGAGTGCTTCCCTCTGCCTAGAGACTC
TGGTGGCTTCTCCAGTTGAGGAGAAACCAGAGGAAAGGGGAGGATTGGGGTCTGGGGGAGGGAACACCATTCACAAAGGCTGACGGTTCCAGTCCGAAGT
CGTGGGCCCACCAGGATGCTCACCTGTCCTTGGAGAACCGCTGGGCAGGTTGAGACTGCAGAGACAGGGCTTAAGGCTGAGCCTGCAACCAGTCCCCAGT
GACTC


Genome files can be downloaded from ensembl ftp site (http://www.ensembl.org/info/data/ftp/index.html) 
or NCBI website (ftp://ftp.ncbi.nlm.nih.gov/genomes/) or other source.

The exon annotation file in GTF format of several speices can be downloaded from Ensembl.

Example:
Homo_sapiens.GRCh37.75.dna.chromosome 1-22, X, Y
ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/

Combine multiple fasta files into one big fasta file (Linux system)

gunzip Homo_sapiens.GRCh37.75.dna.chromosome* 

cat Homo_sapiens.GRCh37.75.dna.chromosome* >human_genome.fa

Homo_sapiens.GRCh37.75.gtf.gz (GTF file)
ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens 

Before perform sgRNAcas9.pl analysis, genome fasta header(description) need be treated by using perl 
script, please see below description



##########################################################################################

3. How to run sgRNAcas9

name:
combine_genome.pl

description: 
This script use to combine multiple fasta sequence files which under the same file folder into one big file.

Input:
suffix of Chromosome/genome sequences, such as "fa","fasta" etc.

Options:
-s <str>			file attributes/suffix, for example, "fa","txt"
-o <str>                        output file

Usage: perl combine_genome.pl -s <.fa> -o <filename>

Example:
Homo_sapiens.GRCh38.dna.chromosome.1.fa 
Homo_sapiens.GRCh38.dna.chromosome.2.fa
Homo_sapiens.GRCh38.dna.chromosome.3.fa
Homo_sapiens.GRCh38.dna.chromosome.4.fa
Homo_sapiens.GRCh38.dna.chromosome.5.fa
Homo_sapiens.GRCh38.dna.chromosome.6.fa
........

perl combine_genome.pl -s fa -o combine.genome.fa

Notes: when run this script, please pay attention that the non-genome files with the 
same suffix(.fa) to genome files should not including in the file folder.

------------------------------------------------------------------------------------------
name: 
format_genome.pl

description: 
This script use to format genome sequence before run sgRNAcas9.pl script.

Input:
a fasta file of the genome.

Options:
-i <str>			Input file

Usage: perl format_genome.pl -i <combine.genome.fa>

Output:
This command will generate a format genome file with removes whitespaces from id line in a fasta file,
only chromosome number will be retained.

Example:
before treated:
>1 dna:chromosome chromosome:GRCh37:1:1:249250621:1 REF

after perl script treated:
>1


------------------------------------------------------------------------------------------

name: 
sgRNAcas9.pl

description: 
This is a main script, which use to identify CRISPR/Cas9 target sequences within a given 
input sequence(s). For instance, if user choose paired gRNA searching mode, sgRNAcas9 will find potential 
target sequences that are/are not present as pairs that can be used as double nickases. Then call the SeqMap 
(Short Sequence Mapping Tool) program to search genome-wide for potential off-target sequences with 
maximum number of mismatches up to 5(default). Finally, the best candidate target for design your custom 
CRISPR guide RNAs will be given.

Input:
A fasta file with target sequences, a fasta file of the corresponding genome

Options ([], represents the default value):
-i <str>		 Input file
-x <int>                 Length of sgRNA[20]
-l <int>	         The minimum value of GC content [20]
-m <int>	         The maximum value of GC content [80]
-g <str>		 The reference genome sequence
-o <str>		 Searching CRISPR target sites using DNA strands based option(s/a/b)
                         [s, sense strand searching mode]
                         [a, anti-sense strand searching mode]
                         [b, both strand searching mode]
-t <str>		 Type of sgRNA searching mode(s/p)
		           [s, single-gRNA searching mode]
		           [p, paired-gRNA searching mode]
-v <str>		 Operation system(w/l/u/m/a)
                            [w, for windows-32, 64]
			    [l, for linux-64]
			    [u, for linux-32]
			    [m, for MacOSX-64]
			    [a, for MacOSX-32]
-n <int>		 Maximum number of mismatches [5]
-s <int>		 The minimum value of sgRNA offset [-2]
-e <int>		 The maximum value of sgRNA offset [32]
-p <str>                 Output path 

Note: for design truncated sgRNA, only permission size range from 17 to 20 nt (option -x).

Also a note on the confusion between "sense" and "antisense" strands: The strand names actually depend on which 
direction you are writing the sequence that contains the information for proteins (the "sense" information), 
not on which strand is on the top or bottom (that is arbitrary). The only real biological information that 
is important for labeling strands is the location of the 5' phosphate group and the 3' hydroxyl group 
because these ends determine the direction of transcription and translation. A sequence 5' CGCTAT 3' is 
equivalent to a sequence written 3' TATCGC 5' as long as the 5' and 3' ends are noted. (source from 
website: http://en.wikipedia.org/wiki/Sense_(molecular_biology)).

If not describe, all the sequence is written in the 5' to 3' direction. The strand you given as input file 
in the 5' to 3' direction will be called the ※sense strand§ in this program. To design paired gRNA, only 
one strand need be given. Make sure that input file without repeat sequences.

For testing purpose, we have provided a test dataset (hEMX1_example.txt) which contains a list of three 
exons from human EMX1 gene in fasta format under the directory "example". We also prepared a sample genome 
file in this directory (genome_example.fa).


Usage: perl sgRNAcas9_3.0.5.pl -i <fasta_input> -x <length_sgRNA> -l <GC%_min> -m <GC%_max> -g <genome.fa> 
-o <s/a/b> -t <s/p> -v <w/l/u/m/a> -s <offset_min> -e <offset_max> -n <number_mismatches> -p <output_path>

 
Example use 1:
The user wishes to design single CRISPR guide RNA to human EMX1 gene, using sense strand and single-gRNA 
searching mode.

perl sgRNAcas9_3.0.5.pl -i hEMX1_example.txt -x 20 -l 40 -m 80 -g genome_example.fa -o s -t s -v w -n 5 -p 2>log.txt

You can also using default options (size of sgRNA: 20 nt, GC content: 20% to 80%, number of mismatches: 5):

perl sgRNAcas9_3.0.5.pl -i hEMX1_example.txt -g genome_example.fa -o s -t s -v w -p

This command will generate a directory (sgRNAcas9.hEMX1_example.txt.report.s) with several intermediate 
directories and files. The important results are:※report_protospacer_single.txt§ showing the CRISPR 
target sequences, location and GC content, and ※final_count_POT.OT_candidate.protospacer.txt§ showing 
the candidate CRISPR target locus for design your custom CRISPR guide RNAs. For more detail output 
information, please see following output explanation.



Example use 2:

The user wishes to design single CRISPR guide RNA to human EMX1 gene, using anti-sense strand and 
single-gRNA searching mode.

perl sgRNAcas9_3.0.5.pl -i hEMX1_example.txt -x 20 -l 40 -m 80 -g genome_example.fa -o a -t s -v w -n 5 -p 2>log.txt

You can also using default options (size of sgRNA: 20 nt, GC content: 20% to 80%, number of mismatches: 5):

perl sgRNAcas9_3.0.5.pl -i hEMX1_example.txt -g genome_example.fa -o a -t s -v w -p

This command will generate the same type of files as example use 1 above.



Example use 3:

The user wishes to design single CRISPR guide RNA to human EMX1 gene, using both strand and 
single-gRNA searching mode.

perl sgRNAcas9_3.0.5.pl -i hEMX1_example.txt -x 20 -l 40 -m 80 -g genome_example.fa -o b -t s -v w -n 5 -p 2>log.txt

You can also using default options (size of sgRNA: 20 nt, GC content: 20% to 80%, number of mismatches: 5):

perl sgRNAcas9_3.0.5.pl -i hEMX1_example.txt -g genome_example.fa -o b -t s -v w -p

This command will generate the same type of files as example use 1 above.



Example use 4:

The user wishes to design paired CRISPR guide RNAs to human EMX1 gene, using both strand and 
paired-gRNA searching mode.

perl sgRNAcas9_3.0.5.pl -i hEMX1_example.txt -x 20 -l 40 -m 80 -g genome_example.fa -o b -t p -v w -n 5 -s 5 -e 35 -p 2>log.txt

You can also using default options (size of sgRNA: 20 nt, GC content: 20% to 80%, number of mismatches: 5, 
sgRNA offset: -2 to 32):

perl sgRNAcas9_3.0.5.pl -i hEMX1_example.txt -g genome_example.fa -o b -t p -v w -p

This command will generate a directory (sgRNAcas9.hEMX1_example.txt.report.b) with file 
※report_protospacer_pairs.txt§showing the CRISPR target sequences, location and GC content, 
and ※final_count_POT.OT_candidate.protospacer.txt§showing the candidate paired CRISPR target 
locus for design your custom paired CRISPR guide RNAs. For more detail output information, 
please see following output explanation.

------------------------------------------------------------------------------------------
name:
ot2gtf.pl

description: 
This script use to determine whether the off-target sites are located in exons.

Input:
A file must contain genome location information, which can be found in the directory of ※Sort_OT_byID§, 
GTF file can be downlaoded from ensemble ftp website.

for example:
hEMX1_exon1_S_109	POT1229	- g - - - - - t - - - - - - - a - a t - - - -	5M	10	7195112	+
hEMX1_exon1_S_109	POT3853	- - - c - - - c - a - - - - - - g g - - - - -	5M	10	30722547	-
hEMX1_exon1_S_109	POT6267	- a - - - t - - - - - - - c - - - - t c - - -	5M	10	50592942	-
hEMX1_exon1_S_109	POT7222	c - - - c a - - - - - - - - - - - - c t - - -	5M	10	61122316	+
hEMX1_exon1_S_109	POT10361	- - - - t - - t - - t - g - - - - g - - - - -	5M	10	82040102	-
hEMX1_exon1_S_109	POT16824	- t a - - g - - - - - - - - - - t - - c - - -	5M	10	126077810	-
hEMX1_exon1_S_109	POT19871	- - - - - - - t - - - - g - g g t - - - - - -	5M	11	993421	+
hEMX1_exon1_S_109	POT25310	- - - - - - a - - - - - g - - c g - - - c - -	4M	11	35640275	-
hEMX1_exon1_S_109	POT29655	c - - - - - a - - g - t - - - - - g - - - - -	5M	11	65343279	+
hEMX1_exon1_S_109	POT37906	- - - a - - - c - - - t - - t - - t - - - - -	5M	11	119784984	+
hEMX1_exon1_S_109	POT54231	- - a t - - - t - - - t - - - - - g - - - - -	5M	12	116824722	-
hEMX1_exon1_S_109	POT70899	- - - - - - - - - - - t - c t - - - a - a - -	4M	14	50359697	+
hEMX1_exon1_S_109	POT74456	- - - - - - - t - - - a - c - - g - a - - - -	5M	14	77328092	+
hEMX1_exon1_S_109	POT76614	- g - - - - t - - - - - - - g - - g - - a - -	4M	14	94566140	+
hEMX1_exon1_S_109	POT77918	- g - - - - - - - a g - - - - - g a - - - - -	5M	14	101001566	+
hEMX1_exon1_S_109	POT78409	- - - - t - a - - - t - - g a - - - - - - - -	5M	14	102510334	+
...............................................................................................................

Options:
-i                     Input file 
-g                     GTF file
-o                     Output file

Usage: perl ot2gtf.pl -i <input> -g <input2> -o <output> 

Example of output:
hEMX1_exon1_S_109	POT78409	- - - - t - a - - - t - - g a - - - - - - - -	5M	14	102510334	+	DYNC1H1
hEMX1_exon1_S_109	POT110542	- - - - - - - a - - - - g g - - g - t - - - -	5M	17	7623121	+	DNAH2
hEMX1_exon1_S_109	POT123622	- - - - - - - - - c a - - - - - g t - - - - -	4M	17	77767383	-	CBX8
hEMX1_exon1_S_109	POT147452	t - - - - - - a - - t - - - - - g t - - - - -	5M	19	51470476	-	KLK6
hEMX1_exon1_S_109	POT156608	- - - - a - - - a c - - - - - c - a - - - - -	5M	1	24105185	+	PITHD1
hEMX1_exon1_S_109	POT201077	- - - - - a - - g - g - - - - - g - a - - - -	5M	22	18775361	-	GGT3P
hEMX1_exon1_S_109	POT202118	- - - - - a - - g - g - - - - - g - a - - - -	5M	22	21576472	-	GGT2
hEMX1_exon1_S_109	POT203238	- - - - - a - - g - g - - - - - g - - - - - -	4M	22	25010752	+	GGT1
hEMX1_exon1_S_109	POT207247	- t - - - - - c - g - c - - - - - - - - g - -	4M	22	42084943	+	C22orf46
hEMX1_exon1_S_109	POT220544	- - - - - - - - - - - - - - - - - - - - - - -	0M	2	73145300	+	EMX1 (#Note: check on target)
hEMX1_exon1_S_109	POT274581	- - - c - g - - g c - - - - - - - t - - - - -	5M	4	85504393	+	CDS1
hEMX1_exon1_S_109	POT325007	a a - - - - - - - - - a - - g c - - - - - - -	5M	6	168394013	-	KIF25-AS1
hEMX1_exon1_S_109	POT338962	c - c - t - - - - - t - t - - - - - - - - - -	5M	7	100609348	+	MUC3A
hEMX1_exon1_S_109	POT344709	- - - - - - g - - - - - - g - - c a - - c - -	4M	7	150675206	-	KCNH2
hEMX1_exon1_S_109	POT349784	c - - t - - - t - - - a - - t - - - - - - - -	5M	8	21655136	+	OR6R2P
hEMX1_exon1_S_109	POT376535	- g - - - - - - - g - - - - - - - a - c c - -	4M	9	116918119	+	COL27A1


------------------------------------------------------------------------------------------
name: 
pot2gtf.pl

description: 
This script use to determine whether the potential off-target sites are located in exons.

Input:
A file must contain genome location information, which can be found in the directory of ※Sort_POT_byID§, 
GTF file can be downloaded from ensemble ftp website.

for example:
hEMX1_exon1_S_109.POT.txt
- - - - - a - - - - - - - - - c - - a - - a -	3M	hEMX1_exon1_S_109_POT105864	16	81490676	-	random_0_3M	
- t - - - - - c - g - c - - - - - - - - g - -	4M	hEMX1_exon1_S_109_POT207247	22	42084943	+	regionI_ident	
- - - - - - - - - - - - - - - - - - - - - - -	0M	hEMX1_exon1_S_109_POT220544	2	73145300	+	seed_ident	
c - - - c - a t - - g - - - - - - - - - - - -	5M	hEMX1_exon1_S_109_POT300654	5	149995015	+	regionI_ident	
c - c - t - - - - - t - t - - - - - - - - - -	5M	hEMX1_exon1_S_109_POT338962	7	100609348	+	regionI_ident	


Options:
-i                     Input file 
-g                     GTF file
-o                     Output file

description: 


Usage: perl pot2gtf.pl -i <input> -g <input2> -o <output> 


Example of output:
- t - - - - - c - g - c - - - - - - - - g - -	4M	hEMX1_exon1_S_109_POT207247	22	42084943	+	regionI_ident		C22orf46
- - - - - - - - - - - - - - - - - - - - - - -	0M	hEMX1_exon1_S_109_POT220544	2	73145300	+	seed_ident		EMX1(#Note: check on target)
c - c - t - - - - - t - t - - - - - - - - - -	5M	hEMX1_exon1_S_109_POT338962	7	100609348	+	regionI_ident		MUC3A


------------------------------------------------------------------------------------------
name: 
check_sgRNA_seq.pl

description: 
This script use to check whether each CRISPR target sequence has a certain unwanted 
nucleotides before design primers for construct sgRNA expression vector.

Input:
a fasta file with candidate CRISPR target sequences

Options:
-i <str>			Input file
-r <str>			restriction enzyme cutting sites

Usage: perl check_sgRNA_seq.pl -i <CRISPR target sequences> -r <restriction enzyme cutting sites>

BsaI
5'...GGTCTC(N)1...3'
3'...CCAGAG(N)5...5'
                                                                (5'-3')
Example: perl check_sgRNA_seq.pl 每i CRISPR.targets_pairs.fa -r "(GGTCTC)|(GAGACC)§

Output:
This command will generate an annotation file ※check.sgR_seq.txt§ with showing whether a CRISPR 
target sequence has a certain unwanted nucleotides, such as restriction enzyme cutting sites, more 
than four continuous T nucleotides (4-6 nucleotide poly(T) tract acts as a termination signal for 
RNA pol III), or other repeat nucleotides (more than 5 A or C or G, more than 6 two nucleotides 
or three nucleotides repeat).



------------------------------------------------------------------------------------------

name: 
sgRPrimer.pl

description: 
This script use to design primer pairs for construct sgRNA expression vector.

Input:
a fasta file of candidate CRISPR target sequences and a file of IDs for select to design CRISPR guide RNAs

Options:
-i <str>			Input file
-s <str>			a file of IDs
-l <str>                        length of sgRNA [20]
-f <str>			restriction enzyme cutting site for forward primer [accg]
-r <str>			restriction enzyme cutting site for reverse primer [aaac]

#design primer pairs for pGL3-U6-gRNA-Puromycin, Huang＊lab
#default, Bsa I: accg, aaac, length:20
#option:length:18, 17

Usage: perl sgRPrimer.pl 每i CRISPR.targets_pairs.fa 每s Select_ID.txt -l 20 -f accg -r aaac

Output:
This command will generate a file ※sgR.Primers.txt§ with given primer pairs for construct sgRNA expression vector.

Example of input:
CRISPR.targets_pairs.fa 
>hEMX1_exon1_S_4
CCAATGGGCGAGCGCGGGGCAGG

Select_ID.txt
hEMX1_exon1_S_4

Example of output:
#Primer pairs
hEMX1_exon1_S_4-20-FP	accgCCAATGGGCGAGCGCGGGGC	
hEMX1_exon1_S_4-20-RP	aaacGCCCCGCGCTCGCCCATTGG	



------------------------------------------------------------------------------------------

name: 
extract_targetSeq.pl

description: 
This script use to extract certain length of flank sequences from genome to down-stream validate on/off target effect.
 
Input:
A file must contain genome location information, which can be found in the directory of ※Sort_POT_byID§.

for example:
- - - - - a - - - - - - - - - c - - a - - a -	3M	hEMX1_exon1_S_109_POT105864	16	81490676	-	random_0_3M	
- t - - - - - c - g - c - - - - - - - - g - -	4M	hEMX1_exon1_S_109_POT207247	22	42084943	+	regionI_ident	
- - - - - - - - - - - - - - - - - - - - - - -	0M	hEMX1_exon1_S_109_POT220544	2	73145300	+	seed_ident	
c - - - c - a t - - g - - - - - - - - - - - -	5M	hEMX1_exon1_S_109_POT300654	5	149995015	+	regionI_ident	
c - c - t - - - - - t - t - - - - - - - - - -	5M	hEMX1_exon1_S_109_POT338962	7	100609348	+	regionI_ident	

Options:
-i <str>			Input file
-g <str>			The reference genome sequence
-l <int>			length of flank sequences [1000]

Usage: perl extract_targetSeq.pl 每i hEMX1_exon1_S_4.POT.txt -g <genome> -l 500

Output:
This command will generate a file with flank sequence from on/off CRISPR target sites for design PCR primer.

Example of output:
>hEMX1_exon1_S_109_POT105864 mismatch=3M random_0_3M 5'-3' >16 Chr_len=90354753, location=81490676, extract: 81490426-81490926 (501 bp) strand(-)
AAGAGAGGGGAGGTAAATGTGGAGTGGAGAAGGCAAGAGGCAAAGGCGACCTGGCTGCCCAGCTGGGGTGGTGCCGCGCCACAGGGGAGCAAGCCCTGGCTCCAACTGGGGAAGGTAGCCCTGGAAAACCACTGACTGGCGGGGG
GGTCGGAGGGGTTCCAAGCAGAAAGAAGAGGCGTGTGCCAAAGGTCAGAGGCATGGGGAAGCTGAGGCCAGAAGCCAGACAGGGCGGGACGCTCGCTCCACAGTAGACTTTGAACTTCGTGCTCACTGTGATCCAGTCCTCACCA
CCGGCTGTGACAGAGCGCTCTGCAGGGTCAGGTGGTCTCAGGCGAGGAGGCTGAGGTGACCACCCCAGGTTCAACGAGGGCAGGTGAGAAGCGATCCTCTGCCACGGCTCCTCCGGAGATGGGTCTGCAGGAAGCTGGCCCAGGA
AGCTCCCGGAAGGGGCCCCAGTGCAGTCAGGGTAGACCCTTCTGCTCTGTCCACACGCTCTCCCTT
>hEMX1_exon1_S_109_POT207247 mismatch=4M regionI_ident 5'-3' >22 Chr_len=51304566, location=42084943, extract: 42084693-42085193 (501 bp) strand(+)
GCGGATGAGACGGCGCAGGCTCTAAGCCTCCATCTAACAGGCTAGGGAGTGAACGCAACCCTGAGTCTCCCGCTTCCCGGTCCCTCGGCCGGCGCACGCACTCACCATGGCTGCGGTTCCGCGGGCTCAGCACTCCTAGGGGAGC
GCAGCTGACGTTTCAGAAGCACTCGCGTGCACCGGAAAAACTCACAGAAGCAGCAGCGGAAATGGCCCCGCGCGGCAGGAAGCGAAGGTGACGCTACGCGAGCGAGTGGGCCCCGCCCTCTACGGGGGCACTGGCGCGGAAACTG
GCCCTGTGTCGAAGAAGGAACGTACTTTGGCGTTCTTATGAGCTGCCTAGTACAAATTATTGGCAGAAACAATGAATTAAACATTCAAAAAGTAACCCACCGGCCGGCCGCAGCGGCTTATGGCTGTAATCCCGGCACTTTGGGA
GGCCGAGGCGGGCGCATCACCTGAGGACAGGAGTTCGAGACCAGCCTGGCCAACATGGCGAAACCT
>hEMX1_exon1_S_109_POT220544 mismatch=0M seed_ident 5'-3' >2 Chr_len=243199373, location=73145300, extract: 73145050-73145550 (501 bp) strand(+)
GCTGCCTCGCACAGCTCCCGCGGCTGCGACCATGTTCCAGCCCGCGGCCAAGCGCGGCTTTACCATAGAGTCCTTGGTGGCCAAGGACGGCGGCACCGGCGGGGGCACTGGCGGCGGGGGCGCGGGCTCCCATCTCCTGGCGGCG
GCCGCCTCCGAGGAACCGCTCCGGCCCACGGCGCTCAACTACCCTCACCCCAGCGCGGCCGAGGCGGCCTTCGTGAGTGGCTTCCCTGCCGCGGCCGCCGCGGGCGCGGGCCGCTCGCTCTACGGTGGGCCCGAGCTCGTGTTCC
CCGAGGCCATGAACCACCCCGCGCTGACCGTGCATCCGGCGCACCAGCTGGGCGCCTCCCCGCTGCAGCCCCCGCACTCCTTCTTCGGCGCCCAGCACCGGGACCCTCTCCATTTCTACCCCTGGGTCCTGCGGAACCGCTTCTT
CGGCCACCGCTTCCAGGGTGAGTGTCCACGCTGTGCCCGCCGAGGCGGCCGGCCGGCGCCCGTGCT
>hEMX1_exon1_S_109_POT300654 mismatch=5M regionI_ident 5'-3' >5 Chr_len=180915260, location=149995015, extract: 149994765-149995265 (501 bp) strand(+)
CCTTGGGAGGGAACCTTGAATTGGAGTGCTCCTGCCAACTCCTCCCTCCTGGCTTGGCAGCGTGAAGACAACGCCGGCCAAAAATATTTGTGTGGGTGGCGGAAAGTTAACTTTTCCGCCTCTCTCTTCTTTCCTGGAAACCGTGG
CCACCGTCAAATACATCAAAACAGAAGCAGGGATAGAGTGAGGCGAGGCCCCTGCGGGGAGTCAGAGGAGCCCTGCCCCCTCCCCACGCCCACTCCTTCTGCCTCCGGCCATCTGGCTCTACGGTGGGGGGATGACGGAGGTGGGT
ATCCAAGTGCTTAGAATCAGGGGATCCTGCAATGCCACAGTTGGACGGATCCTCCAGCATCATCTGGCTCAGTCTCACCATGGTGTCCTCCGGGAATCTGACTACAGAGTGCGAAAAGCACCTGCCCCAGTTCACAGTCCTGTCTC
TTCTGAGTGCTGCTGGCCTTTGTGAGATTCTAAGATTCTTTGAAACTCAAATGGCAGAAATCT
>hEMX1_exon1_S_109_POT338962 mismatch=5M regionI_ident 5'-3' >7 Chr_len=159138663, location=100609348, extract: 100609098-100609598 (501 bp) strand(+)
TCCTTGATGGGGTTGAGTCCAATCCCCTGGTTCTGGGATAGACCCCGCCCACTCATTCTAGGGTGGGGCCCCGCCCCTTCGTTCTAGGGCTGAACCTTGCCCCCTTCTTCTGGGGTGGAGCCCCGCCCCCTTGTTCTAGGGTGGATC
CCCGCCCCCTCCTTTTAGGGTGAAGCCCTGCCCACTTGATCTAAAGTGGAATCCCGCCCCCTCACCTAGGGTAGAGCCCCGCCCCCTCGTTCTAGGGTGGAGACCCGTCCGCTTGTTCTACGGTGGATTCCGGCCGCTTGTCTAGGG
TGGAACCCCCCAGCTTGCCCTAGGGTGGAACCCCCCCGCTGCCCTAGGCTGGAGCCCCGCCCCCTCACCCGCCCCCGCGGGGCCCAGGTGCACGCGTGGACCCCGAGCCCGGAGGTGAAGAGGGTCTGACCCTGCGATCTCCCGCAG
CTGCTACTCCACCGACACGCACTGGTTCTCTGGCCCGCGCTGCGAGGTGGCCGTCCACTG


##########################################################################################

4. sgRNAcas9.pl output files explanation

0utput of sgRNAcas9.pl
A.Sort_OT_byID
B.Sort_POT_byID
POT_seed.identity.xls
sgRNAcas9.Log.txt
sgRNAcas9_report.xls
report_protospacer_pairs.xls
TargetSeq

------------------------------------------------------------------------------------------
represent example files:

file: sgRNAcas9_report.xls

sgRID	Start	End	CRISPR_target_sequence(5'-3')	Length(nt)	GC%	OT	0M(on-/off-)	1M	2M	3M	4M	5M	No.OT	POT	0M(on-/off-)	1M	2M	3M	4M	5M	No.POT	Risk_evaluation
pig_NR1H4_cds_S_1	49	71	TGTCACAGCCACACGTCGCCTGG	23	65.20%	sgR_OT	1	0	0	3	103	245	351	sgR_POT	1	0	0	0	0	0	0	Best
pig_NR1H4_cds_S_2	50	72	GTCACAGCCACACGTCGCCTGGG	23	69.60%	sgR_OT	1	0	2	8	174	350	534	sgR_POT	1	0	0	0	2	0	2	moderate_risk
pig_NR1H4_cds_S_3	60	82	CACGTCGCCTGGGTTTACCATGG	23	60.90%	sgR_OT	1	0	0	5	35	111	151	sgR_POT	1	0	0	0	2	0	2	low_risk
pig_NR1H4_cds_S_4	88	110	ATGAGCATGAAGCCTGCCAAAGG	23	52.20%	sgR_OT	1	0	1	11	172	370	554	sgR_POT	1	0	1	1	2	0	4	moderate_risk
pig_NR1H4_cds_S_5	112	134	GTACTAACAGAACAAGCAGCAGG	23	47.80%	sgR_OT	1	0	1	11	137	380	529	sgR_POT	1	0	1	1	9	9	20	moderate_risk
pig_NR1H4_cds_S_6	120	142	AGAACAAGCAGCAGGTCCTCTGG	23	56.50%	sgR_OT	1	0	2	18	278	594	892	sgR_POT	1	0	0	0	9	4	13	moderate_risk

There are 23 columns in this output file. Their meaning are:
	
field                                                      meaning
sgRID                     ID of the CRISRP target site
Start                     start sites  	
End                       End sites
CRISPR_target_sequence    CRISRP target sequence(5'-3')     
GC%                       GC content  
OT                        annotation of off-targets with No.of mismatches <=5
0M(on-/off-)              on target site, if 0M >1, Repeat gene? or off-target sites with 0 mismatch
1M                        off-target sites with 1 mismatch
2M                        off-target sites with 2 mismatches
3M                        off-target sites with 3 mismatches
4M                        off-target sites with 4 mismatches
5M                        off-target sites with 5 mismatches
No.OT                     total number of off-target sites
POT                       annotation of off-targets with perfect match to the 12-nt seed region of the sgRNA (No.of mismatches <=5)
0M(on-/off-)              perfect match to the 12-nt seed region of the sgRNA, outer segments with 0 mismatch
1M                        perfect match to the 12-nt seed region of the sgRNA, outer segments with 1 mismatch
2M                        perfect match to the 12-nt seed region of the sgRNA, outer segments with 2 mismatches 
3M                        perfect match to the 12-nt seed region of the sgRNA, outer segments with 3 mismatches
4M                        perfect match to the 12-nt seed region of the sgRNA, outer segments with 4 mismatches
5M                        perfect match to the 12-nt seed region of the sgRNA, outer segments with 5 mismatches
No.POT                    total number of off-targets sites with perfect match to the 12-nt seed region of the sgRNA    
Risk_evaluation           evaluation of risk 
           

Discard > High_risk > moderate_risk > low_risk > repeat_sites_or_bad ? > Best  


Choose sgRNA:
The candidate CRISPR sgRNA with minimized off-target effect based on evaluation of risk. 
------------------------------------------------------------------------------------------

file: report_protospacer_pairs.xls:
		Pairs gRNA
sgRID_S	target_seq_S	Start_S	End_S	GC%_S		sgRID_A	target_seq_A	Start_A	End_A	GC%_A	sgRNA_offset(bp)
pig_NR1H4_cds_A_5	ATGACATAAGCATCTCTGCGTGG	1380	1402	47.80%	<->	pig_NR1H4_cds_S_74	CCCCGCTTCTCTGTGAAGTCTGG	1427	1449	60.90%	25
pig_NR1H4_cds_A_5	ATGACATAAGCATCTCTGCGTGG	1380	1402	47.80%	<->	pig_NR1H4_cds_S_75	CCCGCTTCTCTGTGAAGTCTGGG	1428	1450	60.90%	26
pig_NR1H4_cds_A_6	TAAGCATCTCTGCGTGGTGATGG	1374	1396	52.20%	<->	pig_NR1H4_cds_S_74	CCCCGCTTCTCTGTGAAGTCTGG	1427	1449	60.90%	31
pig_NR1H4_cds_A_6	TAAGCATCTCTGCGTGGTGATGG	1374	1396	52.20%	<->	pig_NR1H4_cds_S_75	CCCGCTTCTCTGTGAAGTCTGGG	1428	1450	60.90%	32
pig_NR1H4_cds_A_7	TGAATGTCCGCAACTCAGTCAGG	1350	1372	52.20%	<->	pig_NR1H4_cds_S_73	ACGCAGAGATGCTTATGTCATGG	1382	1404	47.80%	10

For detail information, please enter our website: BiooTools.com

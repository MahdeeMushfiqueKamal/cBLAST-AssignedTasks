# [cBLAST](https://bmb.du.ac.bd/cblast/) Fall'21 internship recruitment problem set and solution. 

## Task 1

RAP-DB is a public repository for rice genome annotation data. The database was last updated on the 10th of May, 2021. Use the data available from this database to complete this task.

Proteins are made up of 20 amino acids. These amino acids occur in varying degrees of abundance. Find the rarest amino acid in rice. Build a simple linear regression model to predict the content of the most abundant amino acid using protein length. Use the model to find the outlier protein that has the largest discrepancy between the prediction and the actual number.

BONUS: Find the amino acid that yields the most robust linear model.

## Task 2

Almost all prokaryotic genes have continuous coding sequences which directly translate to protein products. Eukaryotic genes are more obscure as they have intervening sequences called introns between their coding regions.  

Use the GenBank ID U00096.3 to retrieve the complete genome of Escherichia coli strain K-12 substrain MG1655. Find the discontinuous protein coding genes in this genome. 

Tip: All the data needed to solve this problem are present within the GenBank file and you will not need to look anywhere else.

## Task 3

CRISPR-Cas9 is a powerful and efficient tool for genome editing. To specifically downregulate a gene by knockout using a CRISPR-Cas9 system, a guide RNA sequence complementary to the target is designed. There are some special criteria that guide RNA sequences must fulfill to ensure maximum on-target and minimum off-target effects in the host.

Os03t0752300-01 is the RAP-DB ID of a two-pore potassium channel protein of special interest. Design an appropriate guide RNA for CRISPR knockout of this gene in rice.

Tip: There are open-source software developed for estimating on-target and genome-wide off-target scores of CRISPR gRNA sequences.

BONUS: Design a plasmid construct for delivery of the Cas9 system along with your designed gRNA using Agrobacteria-mediated transformation.

## Task 4 

Collect RNA-seq raw read count data from supplementary table 2 of [Formentin et al. 2018.](https://doi.org/10.3389/fpls.2018.00204) Find differentially expressed genes under their two experimental conditions at a significance level of 0.01. Create a heatmap to show the fold change in differentially expressed genes. 

BONUS: Test all these differentially expressed genes for homology with the two-pore potassium channel and find genes that share at least 60% overall homology.

---

**The python solutions, required data, figures and reports are in this repository.**


# Remarks: 
task1: In order to find the amino acid that yields the model I should have used Coefficient of determination ( R<sup>2</sup> ) instead of average absolute residuals. 

task4 Bonus: Should have used BLAST instead of global alignment using pairwise2 function for finding homology match.  
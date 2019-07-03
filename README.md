# RNA-Seq analysis of Breast Cancer Data.
Repository contians a detailed workflow of RNA-Seq analysis on matched tumour-normal breast cancer samples. Code included in repository for reproducability and learning excercises. 

# Raw Data Analysis
The first step in the workflow is to perform a standard Quality Control analysis of the raw reads provided by the researcher. This analysis is run on a high performance cluster (NUIG lugh) using FASTQC. FASTQC outputs HTML files for each file, of particular interest when performing RNA-Seq are information about the quality of reads, sequence length distribution, overrepresented sequences, and adapter contamination. Once HTML files have been output by FASTQC, the HTMLS can be concatenated into one comprehensive HTML report using MultiQC. 

Code used to run FASTQC can be found in: 
```
FASTQC
```

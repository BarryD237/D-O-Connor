# RNA-Seq analysis of Breast Cancer Data.
Repository contians a detailed workflow of RNA-Seq analysis on matched tumour-normal breast cancer samples. Code included in repository for reproducability and learning excercises. 
***
# Raw Data Analysis
The first step in the workflow is to perform a standard Quality Control analysis of the raw reads provided by the researcher. This analysis is run on a high performance cluster (NUIG lugh) using FASTQC. FASTQC outputs HTML files for each file, of particular interest when performing RNA-Seq are information about the quality of reads, sequence length distribution, overrepresented sequences, and adapter contamination. Once HTML files have been output by FASTQC, the HTMLS can be concatenated into one comprehensive HTML report using MultiQC. 

Code used to run FASTQC can be found in: 
```
Code/FASTQC.sh
```
### FASTQC Sequence Quality: 
![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/fastqc_per_base_sequence_quality_plot_before.png)
### FASTQC Adapter Content:
![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/fastqc_adapter_content_plot_before.png)

We can see that there is ~20% adapter contamination present in the reads. Adapters will effect the alignment of the reads to the human cDNA, as adapters are foreign sequences and will not align anywhere in the human cDNA genome. Ultimately, this results in a loss of data. 
***
# Adapter Trimming
To remove the aforementioned adapters from the reads, we must use a trimming tool. FASTQC reports the suspected source of adapter contamination in the reads (Illumina Universal Adapter). The full adapter sequences can be found at this [link](https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html). 
```
TruSeq LT and TruSeq HT-based kits:

  Read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
  Read 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```
Code used to perform trimming can be found in:
```
Code/Trimming.sh
```
In this script, the file names are provided through a tab seperated text file **file_list.txt** which can be found in the Code directory. It is not neccessary to provide the whole Adapter sequence for trimming, the first 12bp is sufficient as cutadapt will trim the remaining read once it encounters the 12bp sequence. Further information on the flags used for cutadapt are:
* -a Read 1 Adapter Sequence
* -A Read 2 Adapter Sequence
* -o Output for Read 1
* -p Output for Read 2
* -m Miminum Read length

The minimum read length is crucial for downstream quantification with Kallisto. Kallisto by default uses a kmer size of 31, thus if we have reads present that are shorter than this, Kallisto performs poorly. 
***
# FASTQC on Trimmed Data
After removing the adapters, run FASTQC on the trimmed data to make sure the trimming was performed correctly:
### FASTQC Sequence Distribution:
![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/fastqc_sequence_length_distribution_plot.png)
### FASTQC Overrepresented Sequences/Adapter Contamination:
![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/after_trimming.png)

The raw data has succesfully been cleaned and is now ready to be used for downstream transcript quantification. 
***
# Transcript Quantification
To perform transcriptome mapping of the reads, the human cDNA must be downloaded. It is important to use the cDNA, and not the human genome (DNA) to perform this step. Ensembl reference genomes and feature annotation files can be downloaded at the following [link](https://www.ensembl.org/Homo_sapiens/Info/Index). Download the file:
```
Homo_sapiens.GRCh38.cdna.all.fa.gz
```
### Indexing the reference genome
Before mapping reads to the human cDNA genome, an index of the reference file must first be created. This allows for fast random access to the bases in the genome when mapping. Use the mapping script found at:
```
Code/Index.sh
```

### Quantification
Once the reference genome is indexed, we can now map the transcripts in the trimmed files to the reference cDNA. This will output a directory containing:
1. abundance.h5
2. abundance.tsv
3. run_info.json

The two abundance files contain transcript quantification information. For downstream analysis we can use either the .h5 file or the .tsv file. The difference between the .h5 and .tsv file is the .h5 file contains bootstrapping information if the option was specified during **kallisto quant**. We did not perform bootstrapping for this analysis. The code used to run this analysis can be found at:
```
Code/Quantification.sh
```
***
This concludes the analysis performed on a high computing cluster. Moving forward, download the sample directories containing the abundance files via SFTP to your local machine. Further analysis will be conducted within the statistical package R. 
***
# Reading files into R
To read the abundance files into R, follow the code below:
```R
dir <- ("/Users/barrydigby/Desktop/D_OConnor/Quantification")

samples <- read.table(file.path(dir, "Exp_Design.csv"),sep=",", header=T, row.names = "samples")
 
files <- file.path(dir, rownames(samples), "abundance.h5")
names(files) <- paste0(rownames(samples))
```
Exp-Design.csv is a .csv file that contains the metadata for the experiment (Tumour status, grade, size etc..) which can be found at 
```
Code/Exp_Design.csv
```


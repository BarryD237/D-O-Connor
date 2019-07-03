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
It is not neccessary to provide the whole Adapter sequence for trimming, the first 12bp is sufficient as cutadapt will trim the remaining read once it encounters the 12bp sequence. Further information on the flags used for cutadapt are:
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

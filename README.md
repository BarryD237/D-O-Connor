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

#Correction added (R cannot use + and - in factor levels, treats both as X. resulting in duplicates) 
samples$ER <- revalue(samples$ER, c('ER-'='ER'))
samples$PR <- revalue(samples$PR, c('PR-'='PR'))
samples$LVI <- revalue(samples$LVI, c('LVI-'='LVI'))

files <- file.path(dir, rownames(samples), "abundance.h5")
names(files) <- paste0(rownames(samples))
```
Exp-Design.csv is a .csv file that contains the metadata for the experiment (Tumour status, grade, size etc..) which can be found at 
```
Code/Exp_Design.csv
```
# Matching Transcript ID's to Gene ID's
To match the mapped transcripts from the kallisto step to gene id's, we can use biomaRt. As the quantification was performed using ENSEMBL reference genome, the output for kallisto will be using ENSEMBL annotation. Convert as follows:
```R
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

results <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id"), mart = mart)
```
# Convert transcript quantification to gene-level abundances
For this step use the package **tximport**. Per the [bioconductor page](https://bioconductor.org/packages/release/bioc/html/tximport.html): "Imports transcript-level abundance, estimated counts and transcript lengths, and summarizes into matrices for use with downstream gene-level analysis packages. Average transcript length, weighted by sample-specific transcript abundance estimates, is provided as a matrix which can be used as an offset for different expression of gene-level counts."

```R
tx2gene <- results[, 1:2]

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
```
***
# Running DESeq2
Run the follwing R code to provide the read counts to the DDS object:
```R
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ Patient + Condition)
dds$type <- relevel(dds$Condition, ref = "Normal") 
dds<- DESeq(dds)
```
The design is specified as:
```
~ Patient + Condition
```
This design specifies pairwise comparisons for DESeq2, controlling for patient specific factors, resulting in Tumour vs. Normal comparisons. 
***
# PCA
PCA was used to identify further sources of variation that must be accounted for in the DESeq2 model. To perform PCA, first obtain the regularized logarithm of the counts matrix. Then run the following code:
```R
x <- rlog(counts(dds), blind=TRUE)

p <- pca(x, metadata = samples, removeVar = 0.1)

screeplot(p,
  components = getComponents(p, 1:24),
  hline = 80, vline = 10) +
  geom_text(aes(20, 80, label = '80% explained variation', vjust = -1))
```
![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/Scree_plot.png)

```R
biplot(p,
  colby = 'Condition', colkey = c('Tumour'='royalblue', 'Normal'='red1'),
  hline = 0, vline = 0,
  legendPosition = 'right', legendLabSize = 12, legendIconSize = 8.0,
  drawConnectors = TRUE,
  title = 'PCA bi-plot',
  subtitle = 'PC1 versus PC2')
```
![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/bi-plot-PCA.png)

```R
eigencorplot(p,
    components = getComponents(p, 1:10),
    metavars = c('Patient','Condition','Age','Size','Grade','Histology','ER','PR', 'LVI'),
    col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
    cexCorval = 0.7,
    colCorval = 'white',
    fontCorval = 2,
    posLab = 'bottomleft',
    rotLabX = 45,
    posColKey = 'top',
    cexLabColKey = 1.5,
    scale = TRUE,
    main = 'PC1-10 clinical correlations',
    colFrame = 'white',
    plotRsquared = FALSE)
```
![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/eigencor.png)

```R
  eigencorplot(p,
    components = getComponents(p, 1:10),
    metavars = c('Patient', 'Condition','Age','Size','Grade','Histology','ER','PR', 'LVI'),
    col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
    cexCorval = 1.2,
    fontCorval = 2,
    posLab = 'all',
    rotLabX = 45,
    scale = TRUE,
    main = bquote(Principal ~ component ~ Pearson ~ r^2 ~ clinical ~ correlates),
    plotRsquared = TRUE,
    corFUN = 'pearson',
    corUSE = 'pairwise.complete.obs',
    signifSymbols = c('****', '***', '**', '*', ''),
    signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))
```
![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/pearson_eigcor.png)
***
Given the results of the eigen correlation plots, we can see that **patient** and **Condition** are the source of main variation on PC1. These factors are already being accounted for in the DESeq2 model. Interestingly, **Age** has sufficient weight in PC2 to be factored into the DESeq2 model. **PR** also carries weight on PC2, however I have decided to omit this factor from the DESeq2 model. **PR** is an indicator of pathological cancer status, however the model already corrects for cancer status (Tumour/Normal) using **Condition**. Unless the researcher wishes to investigate the data further by comparing **ER / PR / LVI** status, it is recommended to use a basic model at first. 

### Corrected DESeq2 model:
```
design= ~Patient + Age + Condition
```
***
# Inspecting the Data
```R
hist(res$log2FoldChange, breaks=50, col="seagreen", xlab="(Fold change) Tumour vs Normal", main="Distribution of differential expression values")
abline(v=c(-2,2), col="black", lwd=2, lty=2)
legend("topright", "Fold change <-2 and >2", lwd=2, lty=2)
```
![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/dist_exp_val.png)

```R
hist(res$pvalue, breaks=50, col="seagreen", xlab="P-Value (Fold change) Tumour vs. Normal", main="Distribution of P-Values") 
abline(v=c(0.05),col="black",lwd=2,lty=2)
legend("topright", "P-Value <0.05",lwd=2,lty=2)
```
![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/dist_pval.png)

```R
hist(res$padj, breaks=50, col="seagreen", xlab="P-Adj (Fold change) Tumour vs. Normal", main="Distribution of AdjP-Values") 
abline(v=c(0.05),col="black",lwd=2,lty=2)
legend("topright", "P-Adj <0.05",lwd=2,lty=2)
```
![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/dist_adj-pval.png)

```R
resLFC <- lfcShrink(dds, coef="Condition_Tumour_vs_Normal", type="apeglm")

plotMA(resLFC, alpha = 0.05, ylim=c(-4,4))
```
![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/maplot.png)

```R
resOrdered <- res[order(res$padj),]

plot(resOrdered$log2FoldChange, resOrdered$padj, col="black", pch=1, cex=1.1,  main="Tumour vs. Normal Volcano Plot")
fsig <-rownames(subset(resOrdered, log2FoldChange > 1 | log2FoldChange < -1,))
fsig_resOrdered_plot=resOrdered[fsig,]
points(fsig_resOrdered_plot$log2FoldChange, fsig_resOrdered_plot$padj, col="green", pch=19, cex=1.1)
psig <-rownames(subset(resOrdered, pvalue < 0.05))
psig_resOrdered_plot=resOrdered[psig,]
points(psig_resOrdered_plot$log2FoldChange, psig_resOrdered_plot$padj, col="blue", pch=18, cex=0.5)
qsig <-rownames(subset(resOrdered, padj < 0.05))
qsig_resOrdered_plot=resOrdered[qsig,]
points(qsig_resOrdered_plot$log2FoldChange, qsig_resOrdered_plot$padj, col="red", pch=8, cex=0.03)
abline(v=c(1,-1), col="black", lwd=1)
abline(h=0.05, col="red",lwd=1)
legend_text = c("<-1 or >1 Log2-Fold-Change", "p-value<0.05", "Adj-P-Value<0.05")
legend_text2 = c("P-Adj <0.05", "Log2-Fold-Change -1 to 1")
legend("topleft", legend_text,bty="n",pch = c(19,18,8), col=c("green","blue","red"))
legend("topright", legend_text2, bty = "n", lty=c(1, 1), col=c("red", "black"))
```
![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/volcano.png)
***
# Differential Gene Expression Analysis
Analysis was conducted according to the following code:
```R
res <- results(dds, filterFun=ihw, alpha=0.05, name="Condition_Tumour_vs_Normal")
summary(res)
```
```
out of 29357 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 600, 2%
LFC < 0 (down)     : 608, 2.1%
outliers [1]       : 0, 0%
[1] see 'cooksCutoff' argument of ?results
see metadata(res)$ihwResult on hypothesis weighting
```
# Write DE genes to file
The list of up-regulated and down-regulated genes can be found in:
```
Results/DE/
```
These genes have been filtered according to:
* log fold change of 0.05
* Adjusted Pvalue of 0.05 

**this captures all genes returned by summary(res)**

The results are then ranked according to log fold change and statistical significance
```R
tmp=res
up = intersect(rownames(tmp)[which(tmp$log2FoldChange>=0.05)],rownames(tmp)[which(tmp$padj<=0.05)])
down = intersect(rownames(tmp)[which(tmp$log2FoldChange<=-0.05)],rownames(tmp)[which(tmp$padj<=0.05)])

Uptmp = as.data.frame((tmp)[which(rownames(tmp) %in% up),])
Up = Uptmp[order(Uptmp$log2FoldChange),]

Downtmp = as.data.frame((tmp)[which(rownames(tmp) %in% down),])
Down = Downtmp[order(Downtmp$log2FoldChange),]

write.csv(as.data.frame(Up),file = "Up_regulated_genes.csv")

results_csv <- "Up_regulated_genes.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "Up_regulated_genes.txt"

a <- read.table(results_txt, head=TRUE)

b <- getBM(attributes=c("ensembl_gene_id",
                        "external_gene_name",
                        "chromosome_name",
                        "start_position",
                        "end_position",
                        "description"
                        ),
           filters = c("ensembl_gene_id"),
           values = a$X,
           mart = mart)

biomart_results <- "Biomart.txt"

write.table(b,file=biomart_results)

m <- merge(b, a, by.x="ensembl_gene_id", by.y="X")

#change - for upreg, remove for downreg

final = m[order(-m$log2FoldChange,-m$padj),]

write.csv(final, file='Up_Regulated.csv')

# downreg

write.csv(as.data.frame(Down),file = "Up_regulated_genes.csv")


results_csv <- "Up_regulated_genes.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "Up_regulated_genes.txt"

a <- read.table(results_txt, head=TRUE)

b <- getBM(attributes=c("ensembl_gene_id",
                        "external_gene_name",
                        "chromosome_name",
                        "start_position",
                        "end_position",
                        "description"
                        ),
           filters = c("ensembl_gene_id"),
           values = a$X,
           mart = mart)

biomart_results <- "Biomart.txt"

write.table(b,file=biomart_results)

m <- merge(b, a, by.x="ensembl_gene_id", by.y="X")

#change - for upreg, remove for downreg

final = m[order(m$log2FoldChange,m$padj),]

write.csv(final, file='Down_Regulated.csv')
```
# Gene Set Enrichment Analysis (GSEA)
### Formatting files for GSEA input
The required files for GSEA analysis are the following:

1. Expression dataset in res, gct, pcl or txt format
2. Phenotype labels in cls format
3. Gene sets in gmx or gmt formt
4. Chip annotations

#### 1. Expression Dataset
The analysis we are conducting is using RNA-Seq data, thus we can use the normalized counts from DESeq2: 

```R
write.table(counts(dds), file="/Users/barrydigby/Desktop/D_Connor.counts.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
```















GSEA requires expression data, this can be obtained from the dds object in DESeq2 (R):
```R
write.table(counts(dds), file="/Users/barrydigby/Desktop/D_Connor.counts.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
```
Find the output file on your local machine, and edit it by adding "NAME" in row 1, column 1, and hit TAB. This can be done in a text editor. GSEA also requires a "Description" column (can be set to na), insert by  using the following commands in BASH:
```console
awk '{$1 = $1 OFS (NR==1?"Description":"na")}1' D_Connor.counts.txt > test.txt 

sed 's/ /\t/g' test.txt > D_Connor_counts.txt
```
GSEA also requires a .cls file, indicating the phentoype levels. The files used for GSEA input can be found in:
```
Code/D_Connor.cls
Code/D_Connor_counts.txt
```
***

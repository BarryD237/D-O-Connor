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
Before mapping reads to the human cDNA genome, an index of the reference file must first be created. This allows for fast random access to the bases in the genome when mapping. Use the indexing script found at:
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
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
keep <- rowSums(nc >= 10) >= 4
dds <- dds[keep,]
dds <- DESeq(dds)
```
Pre-filtering the dds object was necessary to rectify the error message given below. Two methods can be used to solve this; removing genes with low counts, or increasing the maximum number of iterations performed. I opted to remove genes with low counts, as this yielded a slightly higher number of differentially expressed genes. [link](https://support.bioconductor.org/p/65091/)

> 2 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest

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
Given the results of the eigen correlation plots, we can see that **patient** and **Condition** are the source of main variation on PC1. These factors are already being accounted for in the DESeq2 model. Interestingly, **PR** has sufficient weight in PC2 to warrant further investigation. If PCA shows that age is clustering into groups, then it is worth grouping age into meaningful groups (age bins) and accounting for them in our DESeq2 model. 

#### Investigating Age in PCA
![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/PCA-age.png) 

By splitting the age groups into bins, we can see that no clear clustering is formed. Age is this omitted from the model. Per Michael Love: "I don't like to add "age" alone to the design as it implies that log gene expression increases linearly with years of life.".


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
LFC > 0 (up)       : 576, 2%
LFC < 0 (down)     : 544, 1.9%
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
norm_counts <- counts(dds, normalized=TRUE)

write.table(norm_counts, file="/Users/barrydigby/Desktop/D_Connor.counts.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
```
Inspecting the file, we can see that the columns are shifted to the left:
![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/GSEA_raw.png)

This can be rectified using the following code:
```console
cat D_Connor.counts.txt \
| awk -F"\t" '{if(NR==1) $1="NAME"FS$1}1' OFS="\t" \
| awk '{$1 = $1 OFS (NR==1?"Description":"na")}1' \
| sed 's/ /\t/g' \
| awk 'BEGIN{print "#1.2""\n"40275"\t"24}1' > D.Connor.gct
```

The above code is summarized by the following steps:

```console
awk -F"\t" '{if(NR==1) $1="NAME"FS$1}1' OFS="\t" D_Connor.counts.txt > tmp1.txt
awk '{$1 = $1 OFS (NR==1?"Description":"na")}1' tmp1.txt > tmp2.txt
sed 's/ /\t/g' tmp2.txt > tmp3.txt
echo "$(cat tmp3.txt | wc -l)-1" | bc
cat tmp3.txt | awk '{print NF;exit}'
```

#### 2. Phenotype labels
```console
24 2 1
# Tumour Normal
Tumour Normal Tumour Normal Tumour Normal Tumour Normal Tumour Normal Tumour Normal Tumour Normal Tumour Normal Tumour Normal Tumour Normal Tumour Normal Tumour Normal
```

#### 3. Gene Sets
This can be selected in GSEA 

#### 4. Chip Annotations
The annotation format we are using is Ensembl Gene ID's. In GSEA select ENSEMBL_human_gene.chip
***

**The results of GSEA were not satisfactory, statistically significant pathways had pvalues of 0.00. I investigated other packages in R to carry out the same analysis**

# fgsea
Similar results and output plots can be generated using [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html):

## Native Workflow
In R, we can use the results of DESeq2 to capture information on the statistics of all genes in the differential expression analysis:

```R
res <- results(dds, tidy=TRUE)
write.csv(res, file="/Users/barrydigby/Desktop/deseq_results_tidy.csv")

library(tidyverse)
res <- read_csv("/Users/barrydigby/Desktop/deseq_results_tidy.csv")
```

```R
library(org.Hs.eg.db)
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=res$row, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
ens2symbol <- as_tibble(ens2symbol)
ens2symbol
```

```R
res <- inner_join(res, ens2symbol, by=c("row"="ENSEMBL"))
res
```

```R
res2 <- res %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))
res2

ranks <- deframe(res2)
```

The **fgsea** workflow can now fork to the gene sets of interest

### Hallmarks
"Hallmark gene sets summarize and represent specific well-defined biological states or processes." To use the hallmark gene sets on our dataset, the gene set can be downloaded from the following [link](http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.2/h.all.v6.2.symbols.gmt). 

This step assumed you have completed the fgsea native workflow up to ranks <- deframe(res2). 

```R
pathways.hallmark <- gmtPathways("/Users/barrydigby/Desktop/h.all.v6.2.symbols.gmt")
```

```R
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=10000) %>% 
  as_tibble() %>% 
  arrange(padj)
```

```R
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(pval) %>% 
  DT::datatable()
```

```R
outfile="/Users/barrydigby/Desktop/fgsea/hallmarks/Hallmark_Pathways.pdf"
pdf(file=outfile, width = 8)
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()
```

![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/hallmarks.png)

Using a pvalue cutoff of 0.05, Enrichment plots for each pathway was generated using the following code: 

```R
filtered_pathway <- subset(fgseaResTidy, pval < 0.05)

filt_p <- as.vector(filtered_pathway$pathway)

for (i in filt_p){
    pdf(paste0(i,".pdf"),height=5,width=7.5)
    plt <- plotEnrichment(pathway = pathways.hallmark[[i]], 
    gseaParam = 1, ticksSize = 0.3, stats= ranks) + 
    labs(title=i) + theme(plot.title = element_text(hjust = 0.5, face="bold"))
    print(plt)
    dev.off()
}
```

The results of the 'Hallmark' pathways can be found in:
> Results/GSEA/Hallmarks

Below is an example of an enrichment plot: 

![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/enrichment_plt.png)

This plot uses Normalized Enrichment Score (NES) to plot genes. Breifly, if a peak occurs to the left of the graph, the genes involved in this pathway have positive NES scores, and can be considered upregulated in the pathway. If a pathway is downregulated in Tumour cells, then the plot will have a negative peak on the right hand side, indicating genes involved have a negative NES and are thus downregulated in the pathway. 

The example above shows genes are upregulated in the epithelial to mesenchymal pathway, consistent with tumour cell activity. 

#### Using Enrichment Plots in Paper
To use the enrichment plot in the results of your experiement, I would suggest citing the NES and pvalue of the pathway. 

For each pathway, the NES, pvalue and genes involved in the pathway are included in a statistics report. The report was generated as follows....

```R
#leading edge is a list, need to collapse the list to vector
filtered_pathway$leadingEdge <- vapply(filtered_pathway$leadingEdge, paste, collapse = ",", character(1L))

stats <- data.frame(pathway=filtered_pathway$pathway) %>%
         data.frame(pvalue=filtered_pathway$pval) %>%
         data.frame(NES=filtered_pathway$NES) %>%
         data.frame(genes=filtered_pathway$leadingEdge)

write.csv(stats, file = "/Users/barrydigby/Desktop/fgsea/hallmarks/hallmark_pathway_stats.csv", quote = TRUE, col.names = TRUE, row.names = FALSE) #quote=T to avoid splitting by ','
```

#### Quick view Up/Down Regulated pathways
Two plots of up/down regulated pathways are generated using the log10(pvalue). An example of upregulated pathways is shown below:

```R
filtered_pathway <- subset(fgseaResTidy, pval < 0.05)

filt_up <- subset(filtered_pathway, NES > 0.0)
filt_up$log10 <- -log10(filt_up$pval)
filt_up$pathway <- str_replace_all(filt_up$pathway, c("HALLMARK"="", "_"=" "))

outfile="/Users/barrydigby/Desktop/fgsea/hallmarks/Hallmark_UPREG.pdf"
pdf(file=outfile, width = 8)
ggplot(data=filt_up, aes(x=reorder(pathway,log10), y=log10)) + 
      geom_bar(stat="identity", fill="steelblue", width = 0.7) +
      coord_flip() +
      labs(x="Pathways" y="pvalue (-log10)") +
      theme_minimal()
dev.off()
```

![alt text](https://github.com/BarryD237/D-O-Connor/blob/master/Images/upreg_hallmark.png)

# Williams Anaphylaxis Bulk RNA-seq

## Downstream analysis
[R code](code/20240329_WilliamsBulkRNASeq.preproc.Rmd), [Markdown](code/20240329_WilliamsBulkRNASeq.preproc.md)

<ins>Input</ins>  
raw read counts [TSV](input/williams.genecounts.tsv)  
read alignment statistics [TSV](input/ReadStats.txt)  
sample annotation file [XLSX](input/Williams01_02.15.2024.xlsx)  
gene annotation file *see data release*  
  
<ins>Output</ins>  
SeqExpressionSet with raw counts [RDA](output/williams.esetRaw.RData)  
SeqExpressionSet with normalized counts [RDA](output/williams.eset.RData)  
List of DGEGLM [RDA](output/williams.fits.RData)  

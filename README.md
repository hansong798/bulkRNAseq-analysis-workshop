
# bulk RNAseq-analysis-workshop using R
#### Compiled by Hansong Lee
This is for bulkRNAseq analysis workshop using R.

### Table of Content  
  * [Preparation](#preparation)
  * [Analysis](#analysis)
    * [Step 1. Preprocessing](#step-1-Preprocessing)
    * [Step 2. Filtering](#step-2-Filtering)
    * [Step 3. Run DESeq](#step-3-Run-DESeq)
    * [Step 4. Annotation](#step-4-Annotation)

## Preparation
### Install required packages 
```R
install.packages(c("ggplot2","dplyr","data.table","tidyverse","BiocManager"))
BiocManager::install("biomaRt")   
BiocManager::install("DESeq2")
BiocManager::install("org.Hs.eg.db") 
```
### Load required packages 
```R
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyverse)
library(BiocManager)
library(biomaRt)
library(DESeq2)
library(org.Hs.eg.db)
```

### Download and read GSE dataset
```R
getwd()
raw_count <- read.table('/content/GSE152418_p20047_Study1_RawCounts.txt', header = T, row.names = 1)
raw_count   # Ensembl Gene ID
dim(raw_count)
class(raw_count)

raw_count_matrix <- as.matrix.data.frame(raw_count)
str(raw_count_matrix)
```

## Analysis
### step 1. Preprocessing
#### Create group information
```R
colnames(raw_count)
col <- colnames(raw_count)
group <- col
group[grepl('COV',group)|grepl('CoV',group)] <- 'COVID19'
group[!grepl('COV',group) & !grepl('CoV',group)] <- 'Healthy'
group

group <- factor(group, levels = c('Healthy','COVID19'))  # 첫번째가 reference level이 됨
group
str(group)
table(group)
```

#### Create DESeq object
```R
colData = data.frame(sample = col, group = group)
dds <- DESeqDataSetFromMatrix(raw_count_matrix, colData = colData, design = ~group)
dds
```

### step 2. Filtering
```R
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds
```

### step 3. Run DESeq
```R
dds <- DESeq(dds)
res <- results(dds)
res

#baseMean: 모든 샘플의 normalized average
#log2FoldChange: case vs control log2FoldChange (+)
#lfeSE: log2FC standard error
#stat: Wald test statistics
#pvalue
#padj: corrected pvalue (Benjamini-Hochberg (BH))
```

#### Extract results
```R
summary(res)
res_filt <- results(dds, alpha = 0.01, lfcThreshold = 0.5)
summary(res_filt)

upregulated <- subset(res_filt, padj < 0.01 & log2FoldChange > 2)
downregulated <- subset(res_filt, padj < 0.01 & log2FoldChange < -2)

upregulated
downregulated
```

#### save results as txt files
```R
write.csv(rownames(downregulated), file = 'down_ensembl.txt', sep = '\t', row.names = F, col.names = F)
write.csv(rownames(upregulated), file = 'up_ensembl.txt', sep = '\t', row.names = F, col.names = F)
```

### step 4. Annotation
```R
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
searchDatasets(mart = ensembl, pattern = "Human")
ensembl <- useEnsembl(biomart = 'genes',
                      dataset = 'hsapiens_gene_ensembl')
ensembl_attributes <- listAttributes(ensembl)
head(ensembl_attributes, 20)
```

#### convert ensembl into gene symbol for each up and down DEGs
```R
up_annot <- getBM(attributes= c("ensembl_gene_id", "external_gene_name"),
                  filters = "ensembl_gene_id",
                  values = rownames(upregulated),
                  mart = ensembl)
head(up_annot)

down_annot <- getBM(attributes= c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = rownames(downregulated),
                    mart = ensembl)
head(down_annot)
```

#### save results
```R
write.csv(down_annot$external_gene_name, file = 'down_genename.txt', sep = '\t', row.names = F, col.names = F)
write.csv(up_annot$external_gene_name, file = 'up_genename.txt', sep = '\t', row.names = F, col.names = F)

save(upregulated, downregulated, up_annot, down_annot, file = 'your_path_to_save.RData')
```







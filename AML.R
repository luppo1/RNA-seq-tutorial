# Load all required packages
library("DESeq2")
library("edgeR")
library("limma")

# Set working directory
setwd("C:/Users/Albert Doughan/Desktop/Bioinformatics review/Comprehensive")
dir = "C:/Users/Albert Doughan/Desktop/Bioinformatics review/Comprehensive"

# Import metadata
metadatah = read.csv(file= "metadata_AML.csv", header=T, sep = ",")
head(metadatah)

# Reading counts data from featureCounts
counts <- read.csv("AML.csv", header=TRUE, sep = ",")
head(counts)

# Selecting columns of interest and removing irrelevant ones
countdata = counts[, c(7:28)]
countdata = countdata[, -c(2)]
head(countdata)

# Remove the Gene ID column
countdata <- countdata[, -c(1)]
head(countdata)

# Making "Geneid" column the rownames
rownames(countdata) <- counts[,1]
head(countdata)
dim(countdata)

# Check if the metadata and samples have the same names
table(colnames(countdata)==metadatah$SampleID)


##########################################
# Running differential expression analysis with DESeq2
##########################################

# Create the DESeqDataSet object from Matrix of counts and metadata
dds <- DESeqDataSetFromMatrix(countData = round(countdata), 
                              colData = metadatah, 
                              design = ~Condition)
nrow(dds) 

# Remove Genes with low counts
dds1 <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds1)

# Run DESeq function on the data to perform differential gene expression analysis
dds1 <- DESeq(dds1)
head(assay(dds1))

# Building out results table
res_table <- results(dds1)
summary(res_table)

# We can order our results table by the smallest p value:
order_results <- res_table[order(res_table$pvalue),]
head(order_results)

# How many adjusted p-values were less than 0.1?
sum(order_results$padj < 0.1, na.rm=TRUE)

# Working with alpha 0.05
res2 <- results(dds1, alpha=0.05)
summary(res2)

# How many adjusted p-values were less than 0.05?
sum(res2$padj < 0.05, na.rm=TRUE)

# We order our results table by the smallest p value:
res_small_p <- res2[order(res2$pvalue),]

# Select genes with p less than 0.05
res_sig <- subset(res_small_p, padj < 0.05)
dim(res_sig)

#Write final list to file
write.csv(as.data.frame(res_sig),"deseq2_AML.csv")



##########################################
# Running differential expression analysis with edgeR
##########################################

# Creating the edgeR gene list
d <- DGEList(counts=countdata,group=factor(metadatah$Condition))
dim(d)

# Total gene counts per sample
apply(d$counts, 2, sum) 

# Remove low count genes
keep <- rowSums(cpm(d)>10) >= 2
d <- d[keep,]
dim(d)

# After filtering, it is a good idea to reset the library sizes:
d$samples$lib.size <- colSums(d$counts)
d$samples

# Normalizing the data
d <- calcNormFactors(d)
d

# Estimating the Dispersion
d1 <- estimateCommonDisp(d, verbose=T) # estimating common dispersions
names(d1)
d1 <- estimateTagwiseDisp(d1) # estimating tagwise dispersions
names(d1)

# Testing for DE genes
# Classical Exact test approach
et <- exactTest(d1)
topTags(et, n=10)

#The total number of differentially expressed genes at FDR< 0:05 is:
de1 <- decideTestsDGE(et, adjust.method="BH", p.value=0.05)
summary(de1)

# Order by P-value  
res = et$table
order_res <- res[order(res$PValue),]
dim(order_res)

# Select only genes with P-value less that 0.05
sig_pvalue <- subset(order_res, PValue < 0.05)
dim(sig_pvalue)

# Write final list to file
write.csv(as.data.frame(sig_pvalue), "edgeR_AML.csv")


##########################################
# Running differential expression analysis with limma+voom
##########################################

#Creating the gene list through edgeR
dge <- DGEList(countdata)
dim(dge)

#Create a design matrix
design <- cbind("1"=1,"1vs2"=rep(c(1,2), each = nrow(metadatah)/2))

# Removing genes that are lowly expressed
# Number of genes with 0 count in all samples 
table(rowSums(dge$counts==0)==20)

# The filterByExpr function in the edgeR package provides an automatic way to 
# filter genes, while keeping as many genes as possible with worthwhile counts.
keep.exprs <- filterByExpr(dge, group=metadatah$Condition)
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge)

# Apply scale normalization to counts, TMM normalization method performs 
# well for comparative studies.
dge <- calcNormFactors(dge, method = "TMM")
dge$samples$norm.factors

# Data transformation with CPM
cpm <- cpm(dge)
lcpm <- cpm(dge, log=TRUE)

# Running the limma voom function
v <- voom(dge, design, plot=TRUE, normalize="quantile")

# After this, the usual limma pipelines for differential expression is be applied.
fit <- lmFit(v, design)
fit <- eBayes(fit)
res <- topTable(fit, coef=ncol(design),number=Inf)
res_pvalue <- as.data.frame(subset(res, adj.P.Val < 0.05))
dim(res_pvalue)

# We order our results table by the smallest p value:
order_res <- res_pvalue[order(res_pvalue$adj.P.Val),]
dim(order_res)

#Write final list to file
write.csv(order_res,"limma_AML.csv")


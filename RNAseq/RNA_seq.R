#load libraries
library(RUVSeq)
library(RColorBrewer)
library(Biobase)
library(tidyverse)
library(lme4)
library(reshape2)
library(ggplot2)
library(edgeR)
library(DESeq2)
library(ggfortify)
#Data load and exploration
#load the data, change the long sample name, and select only the columns that I will need.
data_hisat <- read.table("all_hista_counts_mod.txt", head = TRUE, sep="\t", row.names = 1)

data_hisat <- data_hisat %>% 
  dplyr::rename(
    S1_1_5 = X.scratch.project_2000350.genomics.rna.seq.bam.S1_1_5.bam,
    S2_1_5 = X.scratch.project_2000350.genomics.rna.seq.bam.S2_1_5.bam,
    S4_0_015 = X.scratch.project_2000350.genomics.rna.seq.bam.S4_0_015.bam,
    S5_0_015 = X.scratch.project_2000350.genomics.rna.seq.bam.S5_0_015.bam,
    S5_1_5 = X.scratch.project_2000350.genomics.rna.seq.bam.S5_1_5.bam,
    S6_0_015 = X.scratch.project_2000350.genomics.rna.seq.bam.S6_0_015.bam,
    S7_0_015 = X.scratch.project_2000350.genomics.rna.seq.bam.S7_0_015.bam,
    S7_1_5 = X.scratch.project_2000350.genomics.rna.seq.bam.S7_1_5.bam,
    S8_0_015 = X.scratch.project_2000350.genomics.rna.seq.bam.S8_0_015.bam,
    S8_1_5 = X.scratch.project_2000350.genomics.rna.seq.bam.S8_1_5.bam,
    S9_0_015 = X.scratch.project_2000350.genomics.rna.seq.bam.S9_0_015.bam,
    S9_1_5 = X.scratch.project_2000350.genomics.rna.seq.bam.S9_1_5.bam
  )

data1 <- data_hisat %>% select(-(Chr:Length))

#filter for the genes that have at least 5 reads in two samples
filter <- apply(data1, 1, function(x) length(x[x>5])>=2)
filtered <- data1[filter,]
#make a list of the name of the genes and the controls (spikes)
genes <- rownames(filtered)[grep("^NCU", rownames(filtered))] 
spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

#create a list of the conditions.
x2 <- as.factor(c("1.5", "1.5", "0.015", "0.015", "1.5", "0.015", "0.015", "1.5", "0.015", "1.5", "0.015", "1.5"))

#create a EDAseq object
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x2, row.names=colnames(filtered)))

colors <- brewer.pal(3, "Set2")

plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x2])
plotPCA(set, col=colors[x2], cex=1.2)

#Load the data for ERCC spikes
spikedata <- spikes
Y <- counts(set)
sY <- Y[8848:8925,] #Data from spikes (original data)
spikedata <- arrange(spikedata, ERCC_ID) #Now spikes and sY are in the same order
Beta <- spikedata$log2 #This is the fold changes that occurred in ERCC controls
sX <- matrix(rep(0, dim(sY)[1]*dim(sY)[2]), ncol = dim(sY)[2]) #Design matrix of which samples have mix 1
sX[,c(1,2,5,8,10,12)] <- rep(1, dim(sY)[1])

#Take logarithm of Y, and then replace infinite values with 0
logY <- log2(sY)
check <- is.infinite(logY)
logY[check] <- 0
#Calculate spike log normalized counts
spikenorm <- logY - sX*Beta

#For RUVg the data need to be as pseudocounts, so they need to be transformed back
#Store the corrected pseudocounts
counts(set)[8848:8925,] <- round(2^spikenorm)

#There are some spikes that have very low expression
#Drop these from the normalisation
goodspikes <- !apply(check, 1, any)

set <- betweenLaneNormalization(set, which="median") #normalize thedata using upper-quartile (UQ) normalization
set.SN <- RUVg(set, spikes[goodspikes], k=1) #remove unwanted variation
pData(set.SN) #see phenodata


#Create the PCA plot
db.1 <- counts(set.SN)
db.1<-as.data.frame(t(db.1))

db1 <- counts(set.SN)
db1<-as.data.frame(t(db1))
db1$treatment <- x2

pca.res <- prcomp(db.1, scale. = TRUE)
pca.plot <- autoplot(pca.res, data = db1, colour = 'treatment',
                     frame = TRUE, frame.type = 'norm')+
  scale_fill_manual(values=c("#56B4E9", "#FC4E07"))+
  scale_color_manual(values = c("#56B4E9", "#FC4E07"))+ 
  theme_classic()+
  theme(text = element_text(size = 9),
        legend.position = "none")


#Create RLE graph
colnames(filtered)
my_rle_data <- filtered[, c(1, 2, 5, 8, 10, 12, 3, 4, 6, 7, 9, 11)]
x3 <- as.factor(c("1.5", "1.5", "1.5", "1.5", "1.5", "1.5", "0.015", "0.015", "0.015", "0.015", "0.015", "0.015"))
set_graph <- newSeqExpressionSet(as.matrix(my_rle_data),
                                 phenoData = data.frame(x3, row.names=colnames(my_rle_data)))
plotRLE(set_graph, outline=FALSE, ylim=c(-3, 3), col=colors[x3], las=2)


#Differential Expression analysis edgeR
design <- model.matrix(~x2 + W_1, data=pData(set.SN))
y <- DGEList(counts=counts(set.SN), group=x2)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)
test <- lrt$table[which(lrt$table$PValue < 0.05),]

list_edge <- lrt$table
#Differential analysis with DESeq2
dds1 <- DESeqDataSetFromMatrix(countData = counts(set.SN),colData = pData(set.SN),
                               design = ~ W_1 + x2)
dds1 <- DESeq(dds1)
res <- results(dds1)
res
test1 <- res[which(res$padj < 0.05),]
hist(res$pvalue)

list_deseq <- as.data.frame(res)

#Volcano plot
# add a column of NAs
list_deseq$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
list_deseq$diffexpressed[list_deseq$log2FoldChange > 0.6 & list_deseq$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
list_deseq$diffexpressed[list_deseq$log2FoldChange < -0.6 & list_deseq$pvalue < 0.05] <- "DOWN"


volcano.plot <- ggplot(data=list_deseq, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) +
  geom_point(size= 0.2) + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype = "longdash", size = 0.7) +
  geom_hline(yintercept=-log10(0.001), col="darkorange", linetype = "solid", size = 1)+
  geom_hline(yintercept=-log10(0.05), col="black", linetype = "longdash", size = 0.7)+
  theme(text = element_text(size = 16),
        legend.position = "none",
        legend.title = element_blank())

#ANNOTATION FORGE
library(clusterProfiler)
library(rvcheck)
library(AnnotationHub)
library(AnnotationForge)
library(KEGGREST)
library(EnrichmentBrowser)
library(biomaRt)

hub <- AnnotationHub()
query(hub, c("Neurospora crassa","orgdb"))
orgdb <- hub[["AH86823"]]
#KEGG analysis
search_kegg_organism('Neurospora crassa', by='scientific_name') 
#open GeneList
head(list_deseq)
head(list_edge)
geneList <- list_deseq$X
geneList_edge <- list_edge$X

#histograph plot
hist.plot <- ggplot(list_deseq, aes(x=pvalue)) + 
  geom_histogram(color="black")+
  theme_classic()+
  theme(text = element_text(size = 10))+
  xlab("p-value")+
  ylab("Frequency")

#Prepare the data to RUN ORA
# we want the log2 fold change 
original_gene_list <- list_deseq$log2FoldChange
names(original_gene_list) <- list_deseq$X # name the vector
gene_list <-na.omit(original_gene_list) # omit any NA values
gene_list = sort(gene_list, decreasing = TRUE) # sort the list in decreasing order (required for clusterProfiler)

# Extract significant results (padj < 0.05)
sig_genes_df = subset(list_deseq, padj < 0.01)
genes <- sig_genes_df$log2FoldChange # From significant results, we want to filter on log2fold change
names(genes) <- sig_genes_df$X # Name the vector
genes <- na.omit(genes) # omit NA values
genes <- names(genes)[abs(genes) > 2] # filter on min log2fold change (log2FoldChange > 2)

####################Enrichment analysis data created from deseq##########
#The enrichment analysis with edge is the same as deseq.
#KEGG pathway over-representation analysis

kk <- enrichKEGG(gene         = genes,
                 universe = names(gene_list),
                 organism     = 'ncr',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05)
head(kk)
kk[,c(1,2,3,4,5,6,7,8)]
svg("ORA.dotplot.svg")
ora.plot <- dotplot(kk, showCategory = 20, title = NULL, font.size = 8)+
  theme(legend.position = "none")

#kegg gene set enrichment analysis
kk2 <- gseKEGG(geneList     = gene_list,
               organism     = 'ncr',
               minGSSize    = 20,
               pvalueCutoff = 0.05,
               eps = 0,
               pAdjustMethod = "BH",
               by = "fgsea",
               verbose      = FALSE)
head(kk2)
kk2[,c(1,2,3,4,5,6,7,8)]

gsea.plot <- dotplot(kk2, title = NULL , split=".sign", font.size = 8,
                     label_format = "NULL") + 
  facet_grid(.~.sign)+
  theme(text = element_text(size = 8))


gseaplot(kk2, by = "all", title = kk2$Description[14], geneSetID = 14)
browseKEGG(kk2, 'ncr00500')

########################################################################################################################
#Behaviour of the controls
#mean difference graphs
y1 <- DGEList(counts=counts(set), group=x2)
status <- rep("gene", 8925)
status[8847:8925] <- "spikes" 

#Do the data base to do the lowess regresion
logCPM <- cpm(y1, log=TRUE)
AveOfOthers <- rowMeans(logCPM[,-11,drop=FALSE],na.rm=TRUE)
Diff <- logCPM[,11]-AveOfOthers
Mean <- (logCPM[,11]+AveOfOthers)/2

Mean_spikes <- Mean[spikes]
Diff_spikes <- Diff[spikes]

#MD graphs (one per sample)
plotMD(y1,column = 11, status=status, main="S9 0.015%")
abline(h=0, col="blue", lty=2, lwd=2)
lines(lowess(Mean, Diff), col="green", lwd = 3)
lines(lowess(Mean_spikes, Diff_spikes), col="red", lwd=3)

#library size graph
count_lib <- filtered[spikes,]

set_spikes <- newSeqExpressionSet(as.matrix(count_lib),
                                  phenoData = data.frame(x2, row.names=colnames(count_lib)))

y2 <- DGEList(counts=counts(set_spikes), group=x2)
names <- y2$samples[order(y2$samples[,1],decreasing=TRUE),]


barplot(names$lib.size/1e06,names=rownames(names),las=2,col=x1)

#check the behaviour of the spike in controls
#extract the spike counts
data1.1 <- data_hisat %>% select(-(Chr:Strand))
data2 <- data1.1[spikes,]
dim(data2)
write.csv(data2,"~/path/data2.csv")
#data 3 is the data with the concentration modified(multiplied by 0.06)
data3 <- read.csv("~/path/ERCC_Controls_Analysis_mod.csv")
#I needed to load it again becasue I couldn't name the first column
data4 <- read.csv("~/path/data2.csv")
colnames(data4)[1] <- "ERCC_ID" 

data_spikes <- full_join(data3, data4, by="ERCC_ID") 
data_spikes <- na.omit(data_spikes)
write.csv(data_spikes,"~/path/data_spikes.csv")
dim(data_spikes)

#data set modified to be able to do the poisson model
spikes_mod <- read.csv("~/path/spikes_mod.csv",head = TRUE, sep = ",")

M1 <- spikes_mod[which(spikes_mod$Mix == "Mix1"),]
M2 <- spikes_mod[which(spikes_mod$Mix == "Mix2"),]

ggplot(M1, aes(x=log(concentration), y=log(counts+1), color = Sample)) + geom_point()
ggplot(M2, aes(x=log(concentration), y=log(counts+1), color = Sample)) + geom_point()

spikes_mod$Sample <- c(rep("S1 1.5%", 78), rep("S2 1.5%", 78), rep("S4 0.015%", 78),
                       rep("S5 0.015%", 78), rep("S5 1.5%", 78), rep("S6 0.015%", 78),
                       rep("S7 0.015%", 78), rep("S7 1.5%", 78), rep("S8 0.015%", 78),
                       rep("S8 1.5%", 78), rep("S9 0.015%", 78), rep("S9 1.5%", 78))
regression.plot <- ggplot(spikes_mod, aes(x=log(concentration), y=log(counts+1))) + 
  geom_jitter(aes(color = Sample), alpha = 0.6, width = 1) +
  geom_smooth(method = 'lm', method.args = list(family = 'poisson'), 
              col = 'black')+
  theme_classic()+
  theme(text = element_text(size = 10))

spikes_mod$counts_mod <- log(spikes_mod$counts+1)
spikes_mod$concentration_mod <- log(spikes_mod$concentration)

regresion <- lm(counts_mod~ concentration_mod, data= spikes_mod)
summary(regresion)
plot(effect('concentration_mod', regresion))

#Do a poisson model per sample
poisson <- glm(S1_1_5 ~ concentration_M1, data=data_spikes, family='poisson')
summary(poisson)$coefficient
plot(effect('concentration_M1', poisson))

#graph of the poisson regression coeficients
a <- as.numeric(exp(c(0.002370821, 0.002422699, 0.002370624, 0.002406579, 0.002388411, 0.00237606, 0.002520031,
                      0.002513442, 0.002516869, 0.002517968, 0.002524198, 0.002519361)))
b <- unique(spikes_mod$Sample)
c <- cbind(data.frame(b), data.frame(a))
c$d <- as.factor(c("1.5", "1.5", "1.5", "1.5", "1.5", "1.5", "0.015", "0.015", "0.015", "0.015", "0.015", "0.015"))

#Poisson graph
ggplot(c, aes(x = a, y = b, col=d))+
  geom_point(size=5)+
  geom_vline(xintercept = 1.0)+
  xlim(1,1.02)


#Number of genes per GSEA pathways
paths <- read.csv("genes_path.csv")
head(paths)

paths <- transform(paths, histone = factor(histone,
                                           levels = c("H3K9", "H3K27", "H3K36")))

paths <- transform(paths, path = factor(path,
                                        levels = c("Glutathione metabolism",
                                                   "Cysteine and methionine metabolism",
                                                   "Biosynthesis of amino acids",
                                                   " Valine, leucine and isoleucine degradation",
                                                   "Galactose metabolism",
                                                   "Starch and sucrose metabolism",
                                                   "Fructose and mannose metabolism",
                                                   "Glycerolipid metabolism",
                                                   "Ribosome biogenesis in eukaryotes",
                                                   "Fatty acid degradation",
                                                   "Peroxisome",
                                                   "RNA polymerase",
                                                   "Ribosome", 
                                                   "Proteasome")))


mycols <- c('#00CCFF', '#3168FF','#98CBF8')

path.plot <- ggplot(paths, aes(y=path, x=total_gene, fill= histone)) + 
  geom_bar( stat="identity")+
  theme_minimal()+
  theme(text = element_text(size = 12),
        legend.position = "top",
        legend.title = element_blank())+
  xlab("Number of genes")+
  ylab(NULL)+
  scale_fill_manual(values = mycols)


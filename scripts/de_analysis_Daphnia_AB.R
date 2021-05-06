
#creates a series of plots relating to differential expression analysis
# uses the bioconductor package 'limma'

## loading required packages (requires prior installation of each of the following)
require("limma") 
require("edgeR")
require("Rsubread")
require("Biobase")
require("gplots")

### set your path to the scripts directory
WD="/scratch/scwalls/T502_RNAseq/scripts"
setwd(WD)

#obtaining directory paths for both groups
bamDir <- "/scratch/scwalls/DaphniaDevel/STAR_RNA/fastqs"
dp_Annot <- "/scratch/scwalls/T502_RNAseq/annotation/PA42_4.0_onlygeneIDs.saf"
#dp_Annot <- "/scratch/scwalls/T502_RNAseq/annotation/PA42.4.0.test2.gtf"
#dp_Annot <- "/scratch/scwalls/STRIPEseq_pipelines/GoSTRIPES_sing/STRIPES/DpGENOME/PA42.4.0.test2.gtf"


#obtaining list of file names for both groups
## change variable namess to A, B, ect
## change pattern to grab .bam files for the daphnia
files_A <- list.files(bamDir, pattern='*A.*.fastq.gz.STAR.Aligned.sortedByCoord.out.bam$', full.names=TRUE)
files_B <- list.files(bamDir, pattern='*B.*.fastq.gz.STAR.Aligned.sortedByCoord.out.bam$', full.names=TRUE)
files_C <- list.files(bamDir, pattern='*C.*.fastq.gz.STAR.Aligned.sortedByCoord.out.bam$', full.names=TRUE)
files_D <- list.files(bamDir, pattern='*D.*.fastq.gz.STAR.Aligned.sortedByCoord.out.bam$', full.names=TRUE)
files_E <- list.files(bamDir, pattern='*E.*.fastq.gz.STAR.Aligned.sortedByCoord.out.bam$', full.names=TRUE)
files_F <- list.files(bamDir, pattern='*F.*.fastq.gz.STAR.Aligned.sortedByCoord.out.bam$', full.names=TRUE)

daphnia_files <- c(files_A, files_B, files_C, files_D, files_E, files_F)

#creating a count table
daphnia_fc <- featureCounts(daphnia_files, annot.ext=dp_Annot, useMetaFeatures=TRUE, strandSpecific=1, isPairedEnd=TRUE, nthreads=16, isGTFAnnotationFile=FALSE, GTF.featureType="gene", GTF.attrType="gene_id", primaryOnly=TRUE)

save(daphnia_fc, file="daphnia_DE.RData") #saving our featureCounts data to an R binary

## end of read counts section ##

#load("daphnia_DE.RData") #starting from an R binary containing the featureCounts list created using the commands above. To run the above commands simply uncomment them (remove the leading '#' from each individual command), and commend out this line.

dge <- DGEList(counts = daphnia_fc$counts,
               group = c(rep("files_A",4),rep("files_B",4),rep("files_C",2),rep("files_D",5),rep("files_E",3),rep("files_F",1)),
               genes = daphnia_fc$annotation$GeneID)

### Now we will apply TMM normalization
               
dge <- calcNormFactors(dge)

### Let's take a look at the normalization factors

dge$samples

### making a plot of the library counts data
barplot(dge$samples$lib.size, names=c("A1","A2","A3",
                                      "A4", "B1","B2", "B3",
                                      "B4", "C1", "C2", "D1", "D2", "D3", "D4",
                                      "D5", "E1", "E2", "E3", "F1"), las=2, ylim=c(0,1200000))
ggsave(file="library_counts_barplot.png", path="/scratch/scwalls/T502_RNAseq/plots")

### making a plot of the counts value
logcounts <- cpm(dge,log=TRUE)
boxplot(logcounts, xlab="", ylab="(log2) counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
ggsave(file="boxplot_count_data.png", path="/scratch/scwalls/T502_RNAseq/plots")

# Filtering out genes that have a low number of counts (i.e. are lowly-expressed)

keep <- rowSums(cpm(dge)>1) >= 2
dge <- dge[keep, keep.lib.sizes=FALSE]

# Creating a design matrix to model our experiment

#design <- model.matrix(~0 + dge$samples$group)
design <- model.matrix(~dge$samples$group)

colnames(design) <- c("A", "B", "C", "D", "E", "F")

design #what does this object look like?
save(design, file = "design_matrix.RData")

#estimate the dispersion

dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design) 

# evaluate the common dispersion

sqrt(dge$common.disp)

# Now we plot the tagwise dispersions against the log2-scaled counts-per million (CPM) values

plotBCV(dge)
ggsave(file="DP_tagwise_dispersions_nointercept.png", path="/scratch/scwalls/T502_RNAseq/plots")

# Now we performed the differential expression calculation

# fit <- glmFit(dge, design)
fit <- glmQLFit(dge, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)

summary(de <- decideTestsDGE(lrt, p=0.01, adjust="BH"))
de_tags <- rownames(decideTestsDGE(lrt, p=0.01, adjust="BH"))
de_tags <- rownames(dge)[as.logical(de)]

#mkaing a smear (i.e. a mean-difference) plot of our data

plotSmear(lrt, de.tags=de_tags)
abline(h=c(-2,2), col="blue")

# We can also make the (classic) volcano plot from our data
volcanoData <- cbind(lrt$table$logFC, -log10(lrt$table$PValue))
plot(volcanoData, pch=19)
abline(v=c(-2,2), col="red")

save(dge, file="daphniaDGE_2.RData") #saving the updated dge object to our working directory

daphnia_top_tags <- topTags(lrt, n=Inf, adjust.method="BH", sort.by="PValue", p.value=0.01)
head(daphnia_top_tags[[1]]) #shows the top results on the screen
write.csv(daphnia_top_tags[[1]], file="daphnia_top_tags_2.csv", row.names=FALSE) #writes a csv file to your working directory
save(daphnia_top_tags, file= "daphnia_top_tags_2.RData") #saves the daphnia_top_tags file as a p-value

##############
#Making a heatmap with the differentially-expressed genes
library(Biobase) #load this required package if you haven't already done so
library(gplots)

de_data <- dge$counts
colnames(de_data) <- c("A1", "A2","A3",
                       "A4", "B1","B2", "B3",
                       "B4", "C1", "C2", "D1", "D2", "D3", "D4",
                       "D5", "E1", "E2", "E3", "F1" )
head(de_data)

top_tags <- topTags(lrt, n= 18146, sort.by="none")

#differential analysis results
de_data <- cbind(de_data, top_tags[[1]])
head(de_data)

diff.genes = rownames(de_data[de_data$FDR<0.01, ])
head(diff.genes)
length(diff.genes)

dge.subset = dge[diff.genes, ]
colnames(dge.subset$counts) <- c("A1", "A2","A3",
                                 "A4", "B1","B2", "B3",
                                 "B4", "C1", "C2", "D1", "D2", "D3", "D4",
                                 "D5", "E1", "E2", "E3", "F1")
rownames(dge.subset$counts) <- NULL

# plotting the heatmap
heatmap.2(dge.subset$counts,symm=FALSE,symkey=FALSE, scale="row", 
          density.info="none",trace="none", key=TRUE,margins=c(3,3))

dev.off()

# plotting and saving the heatmap to a file
pdf("daphnia_dge_heatmap_2.pdf")
heatmap.2(dge.subset$counts,symm=FALSE,symkey=FALSE, scale="row", density.info="none",trace="none",
          key=TRUE,margins=c(3,3))
dev.off()

#### Done! ######

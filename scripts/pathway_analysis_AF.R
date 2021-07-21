### RNAseq Visualization Automation
#Install RVA from GitHub

devtools::install_github("THERMOSTATS/RVA")

#Load package for use

library(RVA)

#generating summary statistics table
df <- read.csv(file = "DSG_GO_fisher_CC_AB.txt", header = FALSE)



plot_pathway(
  data = df,
  comp.names = NULL,
  gene.id.type = "GO",
  FC.cutoff = 1.3,
  FDR.cutoff = 0.05,
  FCflag = "logFC",
  FDRflag = "adj.P.Val",
  Fisher.cutoff = 0.1,
  Fisher.up.cutoff = 0.1,
  Fisher.down.cutoff = 0.1,
  plot.save.to = NULL,
  pathway.db = "rWikiPathways"
)


# stitching together 
library(dplyr)
TF_file <- read.csv(file = "TF_List_FULL.csv", header = TRUE)
full_join(TF_file,daphnia_fc,copy=TRUE)

# visualization and heatmap
###You can plot a heatmap from raw data rather than a summary statistics table. plot_heatmap.expr has the ability to calculate average expression values and change from baseline. Importantly, these calculations do not calculate statistical signifance or correct for confounding factors - they should not be used as statistical analyses but as data overviews.  

###For this, you need a count table and annotation
###table. The count table should have the geneid as
###row names and the samples as column names. These 
###column names must match the sample.id column in 
###your annotation file:  
  
count <- RVA::count_table[,1:50]

count[1:6,1:5]

#saving counts table to an object

annot <- RVA::sample_annotation[1:50,]

#plotting the heatap of GO terms for gene expression

hm.expr <- plot_heatmap.expr(data = count, 
                             annot = annot,
                             sample.id = "sample_id",
                             annot.flags = c("day", "Treatment"),
                             ct.table.id.type = "ENSEMBL",
                             gene.id.type = "SYMBOL",
                             gene.names = NULL,
                             gene.count = 10,
                             title = "RVA Heatmap",
                             fill = "CPM",
                             baseline.flag = "day",
                             baseline.val = "0",
                             plot.save.to = NULL,
                             input.type = "count")

#saving heatmap and counts table to file

library(ComplexHeatmap)
png("heatmap_plots2cp.png", width = 800, height = 800)
draw(hm.expr$gp)
dev.off()

#plotting gene expression
gene.result <- plot_gene(ct, 
                         anno,
                         gene.names = c("AAAS", "A2ML1", "AADACL3", "AARS"),
                         ct.table.id.type = "ENSEMBL",
                         gene.id.type = "SYMBOL",
                         treatment = "Treatment",
                         sample.id = "sample_id",
                         time = "day",
                         log.option = TRUE,
                         plot.save.to = NULL,
                         input.type = "cpm")

#saving gene expression plot
library(ggplot2)
ggsave(gene.result, "gene_plots1_4.png", device = "png", width = 100, height = 100, dpi = 200, limitsize = FALSE)
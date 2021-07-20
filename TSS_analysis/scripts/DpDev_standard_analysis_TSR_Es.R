### TSRexploreR standard processing module

setwd("/scratch/scwalls/STRIPEseq_pipelines/GoSTRIPES_sing/STRIPES/tsr_Es") 

library(TSRchitect)
library(TSRexploreR)
library(viridis)
library(ggrastr)
library(ggplot2)
library(GenomicRanges)

#load("PdSTRIPE_testObj_tss.RData")

#update file names
tss.1 <- read.table(file="TSSset-1.txt", header=TRUE)
tss.2 <- read.table(file="TSSset-2.txt", header=TRUE)
tss.3 <- read.table(file="TSSset-3.txt", header=TRUE)

tsr.1 <- read.table(file="TSRset-1.tab", header=TRUE)
tsr.2 <- read.table(file="TSRset-2.tab", header=TRUE)
tsr.3 <- read.table(file="TSRset-3.tab", header=TRUE)

tss.1 <- tss.1[,1:4]
tss.2 <- tss.2[,1:4]
tss.3 <- tss.3[,1:4]

tss.1_1 <- tss.1[tss.1$seq=="scaffold_1",]
tss.1_2 <- tss.1[tss.1$seq=="scaffold_2",]
tss.1_3 <- tss.1[tss.1$seq=="scaffold_3",]
tss.1_4 <- tss.1[tss.1$seq=="scaffold_4",]
tss.1_5 <- tss.1[tss.1$seq=="scaffold_5",]
tss.1_6 <- tss.1[tss.1$seq=="scaffold_6",]
tss.1_8 <- tss.1[tss.1$seq=="scaffold_8",]
tss.1_9 <- tss.1[tss.1$seq=="scaffold_9",]
tss.1_10 <- tss.1[tss.1$seq=="scaffold_10",]
tss_r1_1_10 <- rbind(tss.1_1, tss.1_2, tss.1_3, tss.1_4, tss.1_5, tss.1_6, tss.1_8, tss.1_9, tss.1_10)

tss.2_1 <- tss.2[tss.2$seq=="scaffold_1",]
tss.2_2 <- tss.2[tss.2$seq=="scaffold_2",]
tss.2_3 <- tss.2[tss.2$seq=="scaffold_3",]
tss.2_4 <- tss.2[tss.2$seq=="scaffold_4",]
tss.2_5 <- tss.2[tss.2$seq=="scaffold_5",]
tss.2_6 <- tss.2[tss.2$seq=="scaffold_6",]
tss.2_8 <- tss.2[tss.2$seq=="scaffold_8",]
tss.2_9 <- tss.2[tss.2$seq=="scaffold_9",]
tss.2_10 <- tss.2[tss.2$seq=="scaffold_10",]
tss_r2_1_10 <- rbind(tss.2_1, tss.2_2, tss.2_3, tss.2_4, tss.2_5, tss.2_6, tss.2_8, tss.2_9, tss.2_10)

tss.3_1 <- tss.3[tss.3$seq=="scaffold_1",]
tss.3_2 <- tss.3[tss.3$seq=="scaffold_2",]
tss.3_3 <- tss.3[tss.3$seq=="scaffold_3",]
tss.3_4 <- tss.3[tss.3$seq=="scaffold_4",]
tss.3_5 <- tss.3[tss.3$seq=="scaffold_5",]
tss.3_6 <- tss.3[tss.3$seq=="scaffold_6",]
tss.3_8 <- tss.3[tss.3$seq=="scaffold_8",]
tss.3_9 <- tss.3[tss.3$seq=="scaffold_9",]
tss.3_10 <- tss.3[tss.3$seq=="scaffold_10",]
tss_r3_1_10 <- rbind(tss.3_1, tss.3_2, tss.3_3, tss.3_4, tss.3_5, tss.3_6, tss.3_8, tss.3_9, tss.3_10)

#update these ie path to files (DpGENOME)
Dp.annot <- "/scratch/scwalls/STRIPEseq_pipelines/GoSTRIPES_sing/STRIPES/DpGENOME/PA42.4.0_sans_5utr.gff"
Dp.assembly <- "/scratch/scwalls/STRIPEseq_pipelines/GoSTRIPES_sing/STRIPES/DpGENOME/PA42.4.1.fasta"
Dp.assembly.short <- "/scratch/scwalls/STRIPEseq_pipelines/GoSTRIPES_sing/STRIPES/DpGENOME/PA42.4.1_scaf1_10.fasta"

# Generate sample sheet
sample_sheet <- data.frame(
  sample_name=stringr::str_glue("DpDev_{seq_len(3)}"),
  file_1=NA, file_2=NA,
  condition=rep("V", 3)
)


#writing the tss files to the workspace
#tss.1 <- PdSTRIPE@tssCountData[[1]]
#tss.2 <- PdSTRIPE@tssCountData[[2]]
#tss.3 <- PdSTRIPE@tssCountData[[3]]

colnames(tss.1) <- c("seq","TSS", "strand", "score")
colnames(tss.2) <- c("seq","TSS", "strand", "score")
colnames(tss.3) <- c("seq","TSS", "strand", "score")

tsr.1 <- tsr.1[,1:6]
tsr.2 <- tsr.2[,1:6]
tsr.3 <- tsr.3[,1:6]

colnames(tsr.1) <- c("seq", "start", "end", "strand", "nTSS", "nTAGs")
colnames(tsr.2) <- c("seq", "start", "end", "strand", "nTSS", "nTAGs")
colnames(tsr.3) <- c("seq", "start", "end", "strand", "nTSS", "nTAGs")

#making granges files from tss data frames
tss.1.gr <- makeGRangesFromDataFrame(tss.1,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.2.gr <- makeGRangesFromDataFrame(tss.2,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.3.gr <- makeGRangesFromDataFrame(tss.3,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.1.short.gr <- makeGRangesFromDataFrame(tss_r1_1_10,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.2.short.gr <- makeGRangesFromDataFrame(tss_r2_1_10,
                                           keep.extra.columns = TRUE,
                                           seqnames.field="seq",
                                           start.field="TSS",
                                           end.field="TSS",
                                           strand.field = "strand"
)

tss.3.short.gr <- makeGRangesFromDataFrame(tss_r3_1_10,
                                           keep.extra.columns = TRUE,
                                           seqnames.field="seq",
                                           start.field="TSS",
                                           end.field="TSS",
                                           strand.field = "strand"
)

Dp.tss <- list(tss.1.gr, tss.2.gr, tss.3.gr)
names(Dp.tss) <- c("DpE_1", "DpE_2", "DpE_3")

Dp.tss.short <- list(tss.1.short.gr, tss.2.short.gr, tss.3.short.gr)
names(Dp.tss.short) <- c("DpE_1","DpE_2","DpE_3")

#making granges files from tsr data frames
tsr.1.gr <- makeGRangesFromDataFrame(tsr.1,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="start",
                                     end.field="end",
                                     strand.field = "strand"
)

tsr.2.gr <- makeGRangesFromDataFrame(tsr.2,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="start",
                                     end.field="end",
                                     strand.field = "strand"
)

tsr.3.gr <- makeGRangesFromDataFrame(tsr.3,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="start",
                                     end.field="end",
                                     strand.field = "strand"
)

Dp.tsr <- list(tsr.1.gr, tsr.2.gr, tsr.3.gr)
names(Dp.tsr) <- c("DpE_1", "DpE_2", "DpE_3")


#Creating the TSR explorer object
exp <- tsr_explorer(TSSs=Dp.tss, TSRs=Dp.tsr, 
                    genome_annotation=Dp.annot, genome_assembly=Dp.assembly,
                    sample_sheet = sample_sheet
)

#Initial TSS processing

exp <- format_counts(exp, data_type="tss")

#Normalize TSSs
exp <- normalize_counts(exp, data_type = "tss", method = "DESeq2")

## TSS annotation
exp <- annotate_features(exp, data_type = "tss", feature_type="transcript")

### plot correlation

plot_correlation(
  exp, data_type="tss",
  use_normalized=TRUE, font_size=12,
  heatmap_colors=viridis::viridis(100)
)

ggsave(file="DpE_correlation_matrix.png") #saving the plot

### plot genomic distribution of TSSs

plot_genomic_distribution(exp, data_type="tss", samples=c("DpE_1","DpE_2","DpE_3")) +
  scale_fill_viridis_d(direction=-1, name="Annotation") 

ggsave(file="DpE_genomic_distribution.png") #saving the plot

### plot density

plot_density(exp, data_type="tss", samples=c("DpE_1", "DpE_2", "DpE_3"))

#Current error
#Error in `[.data.table`(sample_data, , `:=`(samples, factor(samples, levels = samples))) : 
#Supplied 3 items to be assigned to 384133 items of column 'samples'. If you wish to 'recycle' the RHS please use rep() to make this intent clear to readers of your code.

ggsave(file="DpE_density.png") #saving the plot

### plot TSS heatmap

plot_heatmap(
  exp, data_type="tss", samples=c("DpE_1","DpE_2","DpE_3"),
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="DpE_tss_heatmap.png") #saving the plot

### formatting/normalizing tsrs
exp <- format_counts(exp, data_type="tsr")

## TSR annotation
exp <- annotate_features(exp, data_type = "tsr", feature_type="gene")

### plot TSR heatmap
## this doesn't work yet

plot_heatmap(
  exp, data_type="tsr", samples=c("DpE_1","DpE_2","DpE_3"),
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="DpE_tsr_heatmap.png") #saving the plot

##### Sequence analysis
## creating a truncated object for sequence analysis
## some intervals are too close to the edges of short scaffolds, so this was my workaround

exp_short <- tsr_explorer(TSSs=Dp.tss.short, 
                    genome_annotation=Dp.annot, genome_assembly=Dp.assembly.short,
                    sample_sheet = sample_sheet
)
exp_short <- format_counts(exp_short, data_type="tss")
exp_short <- normalize_counts(exp_short, data_type = "tss", method = "DESeq2")

#plotting sequence logo for all three replicates
plot_sequence_logo(exp_short, samples="DpE_1")
ggsave(file="DpulE_r1_sequenceLogo.png")
##### threshold exploration

plot_threshold_exploration(exp, samples="DpE_1", point_size=1) + scale_colour_viridis_c()
ggsave(file="DpulE_Count_Threshold.png")
exp <- apply_threshold(exp, threshold=3, n_samples=1)

plot_density(exp, data_type = "tss", samples = "DpE_1")
ggsave(file="DpulE_Densityplot.png")

plot_dinucleotide_frequencies(exp, samples = "DpE_1") + scale_fill_viridis_c()
ggsave(file="Dinucleotide_Frequencies.png")


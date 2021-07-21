library(topGO)


## load the gene annotation file
geneID2GO <- readMappings(file = "/scratch/scwalls/T502_RNAseq/scripts/blast2go_go_propagation_20190313_1159_col12.txt")

## read in the count table - all Eu and St DSG genes
#### REPLACE read.csv file between each comparison
ppa_dsg <- read.csv("/scratch/scwalls/T502_RNAseq/scripts/daphnia_top_tags_AF.csv", header = F, stringsAsFactors=FALSE,sep="\t",skip=1)

dsg_gene <- unique(ppa_dsg$Gene.ID)

geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% dsg_gene))
names(geneList) <- geneNames
str(geneList)

length(dsg_gene)

ppaGOdata <- new("topGOdata",
                 description = "Genes differentially spliced", 
                 ontology = "CC",
                 allGenes = geneList,
                 annot = annFUN.gene2GO,
                 nodeSize = 3,
                 gene2GO = geneID2GO)

resultFisher <- runTest(ppaGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(ppaGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(ppaGOdata, algorithm = "elim", statistic = "ks")
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(ppaGOdata, test.stat)
resultFisher

allRes <- GenTable(ppaGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 50)

printGraph(ppaGOdata, resultKS, firstSigNodes = 5, useInfo = "all", pdfSW = TRUE,fn.prefix="DSG_CC_AF")
write.table(allRes,file="DSG_GO_fisher_CC_AF.txt",col.names=TRUE,row.names=FALSE,quote=TRUE,sep="\t")

## Eu genes shared with eud1(OE) and nhr-40 mutant + seud1-mutant

ppa_eu <- read.csv("/scratch/scwalls/T502_RNAseq/scripts/description_per_gene.txt", header = F, sep = "", skip = 1)

eu_gene <- unique(ppa_eu$geneID)
geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% eu_gene))
names(geneList) <- geneNames
str(geneList)

length(eu_gene)

ppaGOdata <- new("topGOdata",
                 description = "Eu shared genes", 
                 ontology = "MF",allGenes = geneList,
                 annot = annFUN.gene2GO,nodeSize = 3,gene2GO = geneID2GO)

resultFisher <- runTest(ppaGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(ppaGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(ppaGOdata, algorithm = "elim", statistic = "ks")
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(ppaGOdata, test.stat)
resultFisher

allRes <- GenTable(ppaGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 50)

printGraph(ppaGOdata, resultKS, firstSigNodes = 5, useInfo = "all", pdfSW = TRUE,fn.prefix="Eu_shared_MF")
write.table(allRes,file="Eu_shared_GO_fisher_MF.txt",col.names=TRUE,row.names=FALSE,quote=TRUE,sep="\t")

## St genes shared with eud1(OE) and nhr-40 mutant + seud1-mutant

ppa_st <- read.csv("~/Desktop/Transcriptomic_related/St_shared.sorted.csv", header = T)

st_gene <- unique(ppa_st$geneID)
geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% st_gene))
names(geneList) <- geneNames
str(geneList)

length(st_gene)

ppaGOdata <- new("topGOdata",description = "St shared genes", ontology = "BP",allGenes = geneList,
                 annot = annFUN.gene2GO,nodeSize = 3,gene2GO = geneID2GO)

resultFisher <- runTest(ppaGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(ppaGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(ppaGOdata, algorithm = "elim", statistic = "ks")
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(ppaGOdata, test.stat)
resultFisher

allRes <- GenTable(ppaGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 50)

printGraph(ppaGOdata, resultKS, firstSigNodes = 5, useInfo = "all", pdfSW = TRUE,fn.prefix="St_shared_BP")
write.table(allRes,file="St_shared_GO_fisher_BP.txt",col.names=TRUE,row.names=FALSE,quote=TRUE,sep="\t")



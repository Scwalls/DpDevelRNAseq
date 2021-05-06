#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')

#BiocManager::install('PCAtools')

library(PCAtools)
library(Biobase)
library(GEOquery)

# load series and platform data from GEO
gset <- getGEO('GSE2990', GSEMatrix = TRUE, getGPL = FALSE)
mat <- exprs(gset[[1]])

# remove Affymetrix control probes
mat <- mat[-grep('^AFFX', rownames(mat)),]

# extract information of interest from the phenotype data (pdata)
idx <- which(colnames(pData(gset[[1]])) %in%
               c('relation', 'age:ch1', 'distant rfs:ch1', 'er:ch1',
                 'ggi:ch1', 'grade:ch1', 'size:ch1',
                 'time rfs:ch1'))
metadata <- data.frame(pData(gset[[1]])[,idx],
                       row.names = rownames(pData(gset[[1]])))

# tidy column names
colnames(metadata) <- c('Study', 'Age', 'Distant.RFS', 'ER', 'GGI', 'Grade',
                        'Size', 'Time.RFS')

# prepare certain phenotypes of interest
metadata$Study <- gsub('Reanalyzed by: ', '', as.character(metadata$Study))
metadata$Age <- as.numeric(gsub('^KJ', NA, as.character(metadata$Age)))
metadata$Distant.RFS <- factor(metadata$Distant.RFS,
                               levels = c(0,1))
metadata$ER <- factor(gsub('\\?', NA, as.character(metadata$ER)),
                      levels = c(0,1))
metadata$ER <- factor(ifelse(metadata$ER == 1, 'ER+', 'ER-'),
                      levels = c('ER-', 'ER+'))
metadata$GGI <- as.numeric(as.character(metadata$GGI))
metadata$Grade <- factor(gsub('\\?', NA, as.character(metadata$Grade)),
                         levels = c(1,2,3))
metadata$Grade <- gsub(1, 'Grade 1', gsub(2, 'Grade 2', gsub(3, 'Grade 3', metadata$Grade)))
metadata$Grade <- factor(metadata$Grade, levels = c('Grade 1', 'Grade 2', 'Grade 3'))
metadata$Size <- as.numeric(as.character(metadata$Size))
metadata$Time.RFS <- as.numeric(gsub('^KJX|^KJ', NA, metadata$Time.RFS))

# remove samples from the pdata that have any NA value
discard <- apply(metadata, 1, function(x) any(is.na(x)))
metadata <- metadata[!discard,]

# filter the expression data to match the samples in our pdata
mat <- mat[,which(colnames(mat) %in% rownames(metadata))]

# check that sample names match exactly between pdata and expression data 
all(colnames(mat) == rownames(metadata))

## Performing principal component analysis
p <- pca(mat, metadata = metadata, removeVar = 0.1)

save(metadata, file="metadata.RData")
save(mat, file="mat.RData")

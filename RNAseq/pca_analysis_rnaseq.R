library(PCAtools)
library(Biobase)

load("dp_matrix.RData")
load("dp_metadata.RData")

#check that sample names match exactly between pdata and expression data (should equal TRUE)
all(colnames(dp.mat2) == rownames(dp.metadata))
    
#Conduct principal component analysis (PCA):
p <- pca(dp.mat2, metadata = dp.metadata, removeVar = 0.1)
    
# a scree plot
screeplot(p, axisLabSize = 18, titleLabSize = 22)

# creating a biplot
biplot(p)

#creating a pairs plot
pairsplot(p)

library(clusterpval)
library(fastcluster)
library(DropletUtils)
library(scater)
library(ggfortify)
library(patchwork)

source("util.R")

##### Pre-processing data #### 
cd4.t <- read10xCounts("./raw/filtered_matrices_mex_memory/hg19")
cd19.b <- read10xCounts( "./raw/filtered_matrices_mex_bcell/hg19")
mono <-  read10xCounts( "./raw/filtered_matrices_mex_mono/hg19")

set.seed(1)
subset.t <- sample(1:ncol(cd4.t), 800)
subset.b <- sample(1:ncol(cd19.b), 200)
subset.m <- sample(1:ncol(mono), 200)

pure.t <- cd4.t[, subset.t[1:600]]
mix3 <- cbind(cd4.t[, subset.t[601:800]], cd19.b[,subset.b], mono[, subset.m])
colData(mix3)$CellType <- c(rep("Memory T cell", 200), rep("B cell", 200), rep("Monocyte", 200))

# Define mitochondrial, ribosomal protein genes
rowData(pure.t)$Mito <- grepl("^MT", rowData(pure.t)$Symbol)
rowData(pure.t)$ProteinRibo <- grepl("^RP", rowData(pure.t)$Symbol)

# Add count of unique genes
colData(pure.t)$total_features <- apply(counts(pure.t), 2, function(x) sum(x != 0))

# Calculate quality statistics
pure.t <- calculateQCMetrics(pure.t, feature_controls = list(mt=rowData(pure.t)$Mito , rbp=rowData(pure.t)$ProteinRibo))

# Delete bad cells (low library size, low gene count, high mitochondrial read proportion)
libsize.drop <- isOutlier(pure.t$total_counts, nmads = 3, type = "lower", log = TRUE)
feature.drop <- isOutlier(pure.t$total_features, nmads = 3, type = "both", log = TRUE)
mito.drop <- isOutlier(pure.t$pct_counts_mt, nmads = 3, type = "higher", log = TRUE)
pure.t <- pure.t[,!(libsize.drop | feature.drop |mito.drop)]

# Normalize library sizes, then log2 transformation with pseudo-count of 1
pure.t <- logNormCounts(pure.t)

# Subset to top 500 average count genes
X1 <- t(logcounts(pure.t))
X1 <- as.matrix(X1)
X1 <- X1[, order(rowData(pure.t)$mean_counts, decreasing=T)[1:500]]

# Define mitochondrial, ribosomal protein genes
rowData(mix3)$Mito <- grepl("^MT", rowData(mix3)$Symbol)
rowData(mix3)$ProteinRibo <- grepl("^RP", rowData(mix3)$Symbol)

# Add count of unique genes
colData(mix3)$total_features <- apply(counts(mix3), 2, function(x) sum(x != 0))

# Calculate quality statistics
mix3 <- calculateQCMetrics(mix3, feature_controls = list(mt=rowData(mix3)$Mito , rbp=rowData(mix3)$ProteinRibo))

# Delete bad cells (low library size, low gene count, high mitochondrial read proportion)
libsize.drop <- isOutlier(mix3$total_counts, nmads = 3, type = "lower", log = TRUE)
feature.drop <- isOutlier(mix3$total_features, nmads = 3, type = "both", log = TRUE)
mito.drop <- isOutlier(mix3$pct_counts_mt, nmads = 3, type = "higher", log = TRUE)
mix3 <- mix3[,!(libsize.drop | feature.drop |mito.drop)]

# Normalize library sizes, then log2 transformation with pseudo-count of 1
mix3 <- logNormCounts(mix3)

# Subset to top 500 average count genes
X2 <- t(logcounts(mix3))
X2 <- as.matrix(X2)
X2 <- X2[, order(rowData(mix3)$mean_counts, decreasing=T)[1:500]]

##### Data analysis #### 
# Negative control
hcl <- hclust(dist(X1)^2, method="ward.D")
plot(hcl)
K <- 3

# Visualize data and estimated clusters
pca_data <- prcomp(X1, rank=2)
ggplot(pca_data$x, aes(x=PC1, y=PC2, colour=as.factor(cutree(hcl, K)))) + geom_point() +
  coord_fixed(ratio=1.38/3.33) + ggtitle("Negative control") + ylab("PC2 (1.38%)") +
  xlab("PC1 (3.33%)") + scale_colour_discrete(name="Estimated clusters")  +  theme_bw(base_size=22)  +
  theme(legend.position="bottom")

# Use POET with 5 factors to estimate inverse of covariance matrix
SampleCovX1 <- cov(X1)
EigenSampleCovX1 <- eigen(SampleCovX1)
EigenSampleCovX1$values[1:10] # look for large spiked values to pick number of factors
# The POET paper says that when in doubt, more is better - robust to over-specifying factors, 
# not robust to under-specifying factors
# We pick five factors

p <- ncol(X1)
SigX <- POET(t(scale(X1, scale=F)), K=5, C=0.5, thres="soft", "vad")
Su <- SigX$SigmaU
B <- SigX$loadings
SuInv <- solve(Su)

SigX1 <- B%*%t(B) + Su
SigXInv1 <- SuInv - SuInv%*%B%*%solve(diag(5) + t(B)%*%SuInv%*%B)%*%t(B)%*%SuInv

# Wald tests
wald_test(X1, link="ward.D", K=3, k1=1, k2=2, iso=FALSE, SigInv=SigXInv1)
wald_test(X1, link="ward.D", K=3, k1=1, k2=3, iso=FALSE, SigInv=SigXInv1)
wald_test(X1, link="ward.D", K=3, k1=2, k2=3, iso=FALSE, SigInv=SigXInv1)

# Our tests
test_hier_clusters_exact(X1, link="ward.D", K=3, k1=1, k2=2, iso=FALSE, SigInv=SigXInv1, hcl=hcl)
test_hier_clusters_exact(X1, link="ward.D", K=3, k1=1, k2=3, iso=FALSE, SigInv=SigXInv1, hcl=hcl)
test_hier_clusters_exact(X1, link="ward.D", K=3, k1=2, k2=3, iso=FALSE, SigInv=SigXInv1, hcl=hcl)

# Positive control
# Visualize data 
pca_data2 <- prcomp(X2, rank=2)
pca_data2$x <- data.frame(PC1=pca_data2$x[, 1], PC2=pca_data2$x[, 2], CellType=colData(mix3)$CellType)
ggplot(pca_data2$x, aes(x=PC1, y=PC2, colour=CellType)) + geom_point() +
  coord_fixed(ratio=6.84/18.23) + ggtitle("Positive control") + xlab("PC2 (6.84%)") +
  ylab("PC1 (18.23%)") + scale_color_brewer(palette="Dark2", name="Cell types") + theme_bw(base_size=22)  +
  theme(legend.position="bottom") + guides(colour = guide_legend(nrow =2))

hcl2 <- hclust(dist(X2)^2, method="ward.D")
K <- 3
table(colData(mix3)$CellType, cutree(hcl2, K)) # estimated clusters line up well w/ cell types 

# Use POET with 5 factors to estimate inverse of covariance matrix
SampleCovX2 <- cov(X2)
EigenSampleCovX2 <- eigen(SampleCovX2)
EigenSampleCovX2$values # look for large spiked values to pick number of factors

p <- ncol(X2)
SigX <- POET(t(scale(X2, scale=F)), K=5, C=0.5, thres="soft", "vad")
Su <- SigX$SigmaU
B <- SigX$loadings
SuInv <- solve(Su)

SigX2 <- B%*%t(B) + Su
SigXInv2 <- SuInv - SuInv%*%B%*%solve(diag(5) + t(B)%*%SuInv%*%B)%*%t(B)%*%SuInv

# Wald tests
wald_test(X2, link="ward.D", K=3, k1=1, k2=2, iso=FALSE, SigInv=SigXInv2)
wald_test(X2, link="ward.D", K=3, k1=1, k2=3, iso=FALSE, SigInv=SigXInv2)
wald_test(X2, link="ward.D", K=3, k1=2, k2=3, iso=FALSE, SigInv=SigXInv2)

# Our tests
test_hier_clusters_exact(X2, link="ward.D", K=3, k1=1, k2=2, iso=FALSE, SigInv=SigXInv2, hcl=hcl2)
test_hier_clusters_exact(X2, link="ward.D", K=3, k1=1, k2=3, iso=FALSE, SigInv=SigXInv2, hcl=hcl2)
test_hier_clusters_exact(X2, link="ward.D", K=3, k1=2, k2=3, iso=FALSE, SigInv=SigXInv2, hcl=hcl2)


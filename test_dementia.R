library(pracma)
library(DESeq2)
library(plotly)
library(ggfortify)

# generating random integer matrix representing read count mapped to a given gene
count_matrix<-randi(imax=10000, n = 6, m = 10000)

# creating metadata table for target groups
A<-matrix(c('A','male','A','male','A','female','B','male',"B",'male','B','female'),nrow = 6,ncol = 2,byrow = TRUE)
rownames(A)<-c(1,2,3,4,5,6)
colnames(A)<-c("group", "sex")

# normalizing counts for PCA
f<-apply(count_matrix,2,norm<-function(x){return (x/sum(x))})
# performing PCA
pca_res <- prcomp(f, scale. = TRUE)
# plotting PCA to get outlook on data before performing DE
p<-autoplot(pca_res, data = A, colour = 'group')
ggplotly(p)

# creating DESeq object for differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = t(count_matrix),colData = A,design = ~ group)
# performing DE analysis
DESeq2obj <- DESeq(dds,fitType='local') 
# getting results from analysis
DESeq2res <- results(DESeq2obj, contrast = c("group", "A", "B"), pAdjustMethod = "BH", 
                     alpha = 0.05,
                     lfcThreshold = 0.58)
# printing results
print(summary(DESeq2res, alpha = 0.05))

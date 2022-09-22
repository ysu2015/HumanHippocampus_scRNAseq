library("Seurat")
library("ggplot2")
install.packages('randomForest')
source("./sourceCode/plotConfusionMatrix.R")
source("./sourceCode/plot_sample_dis.R")
source("./sourceCode/RF_train.R")
library(randomForest)

## using compare to vistal cortex on adult stage as example ##
## start to compare ## 
Hipp.Adult<-readRDS("./CurrentStudy_Hippocampus/Adult_HIPP.rds")
Hipp.Adult
table(Hipp.Adult$Sample) # 7 donors
DefaultAssay(Hipp.Adult)<-"integrated"
Idents(Hipp.Adult)<-"MajorCellTypes"
table(Hipp.Adult@meta.data$MajorCellTypes)

Adult.vis<-readRDS("./referenceData/Adult.VisCor.rds")
head(Adult.vis[[]])
table(Adult.vis[["sampleID"]])
Adult.vis_clusters<-Adult.vis@meta.data$Big_lake_clusters
Adult.vis<-NormalizeData(Adult.vis, normalization.method="LogNormalize", scale.factor=10000)
Adult.vis <- FindVariableFeatures(Adult.vis, selection.method = "vst", nfeatures = 2000)
var.genes<-intersect(VariableFeatures(Hipp.Adult),VariableFeatures(Adult.vis))
rf_model<-RF_train(Hipp.Adult,var.genes,do.scale = FALSE)

pred.labels<-predict(rf_model, t(t(scale(t(Adult.vis[["RNA"]]@data[var.genes,])))))
names(pred.labels)<-colnames(Adult.vis@assays$RNA@data)
plotConfusionMatrix(table(Adult.vis_clusters,pred.labels),order="Row")

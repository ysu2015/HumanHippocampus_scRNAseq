
library("randomForest")
library(Matrix)
library(reshape2)
RF_train <- function(train_object0, var.genes, do.scale=FALSE){
 predictor_Data = train_object0[["RNA"]]@counts[var.genes,]
# predictor_Data = as.matrix(train_object0[["SCT"]]@data[var.genes,])
  if (do.scale) predictor_Data = t(scale(t(predictor_Data)))
  max.cells.per.ident = 700; train.frac = 0.6
  training.set = c(); validation.set=c()
  training.label = c(); validation.label=c();
  print(paste0("Using mininum of ", 0.6*100, " percent cells or ", max.cells.per.ident, " cells per cluster for training"))
  for (i in as.character(levels(Idents(train_object0)))){
    cells.in.clust = WhichCells(train_object0,idents = i);
    n = min(max.cells.per.ident, round(length(cells.in.clust)*train.frac))
    train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
    validation.temp = setdiff(cells.in.clust, train.temp)
    training.set = c(training.set,train.temp); 
    validation.set=c(validation.set,validation.temp)
    training.label = c(training.label, rep(i,length(train.temp))); 
    validation.label = c(validation.label, rep(i, length(validation.temp)));
  }
  
  print(length(training.set))
  tmp = as.vector(table(training.label))
  sampsizes = rep(min(tmp),length(tmp))
  input.matrix<-as.matrix(t(predictor_Data[,training.set]))
  rf <- randomForest(x=input.matrix, 
                    y=factor(training.label), importance = TRUE, 
                    ntree = 301, proximity=TRUE, sampsize=sampsizes, 
                    keep.inbag=TRUE, replace=FALSE) 
  validation_pred_rf = predict(rf,t(predictor_Data[,validation.set]))
 plotConfusionMatrix(table(validation.label, validation_pred_rf))
  return(rf)
}

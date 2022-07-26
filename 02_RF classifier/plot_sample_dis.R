library(reshape2)
plot_sample_dist <- function(object, id.use = id.use, max.val.perc=100,
                             max.size=10,row.scale=FALSE,top.mar=1, left.mar = 4.5,
                             right.mar=1.5, bottom.mar = -1,do.transpose=FALSE,
                             ident.order=NULL,x.lab.rot=0,to.return=TRUE,alpha.use=1) 
{
  
  if (is.null(ident.order)){
    ident.order = levels(Idents(object))
  }
  
  
  # Get sample distribution
  ExpMat = table(object@meta.data[,id.use],Idents(object))
  ExpMat = ExpMat[,ident.order]
  if(row.scale){
    ExpMat = t(scale(t(ExpMat), center=FALSE, scale=rowSums(ExpMat)))*100
  } else {
    ExpMat = scale(ExpMat, center=FALSE, scale=colSums(ExpMat))*100
  }
  ExpVal = melt(ExpMat)
  head(ExpVal)
  colnames(ExpVal) = c("sample","cluster","perc")
  # ExpVal$sample = factor(ExpVal$sample, levels=rev(levels(object@meta.data[,id.use])))
  ExpVal$cluster = factor(ExpVal$cluster, levels= ident.order)
  p=ggplot(ExpVal, aes(y = factor(sample),  x = factor(cluster))) + geom_point(aes(colour = perc,  size =perc),alpha=alpha.use) + 
    scale_color_gradient(low ="darkgreen",   high = "red", limits=c(0, 100))+scale_size(range = c(1, 10), limits = c(0,100))+
    theme_bw() +nogrid
  p = p + xlab("Cluster") + ylab("Sample") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
    theme(axis.text.y=element_text(size=12, face="italic")) + theme(axis.text.x = element_text(angle = x.lab.rot))
  
  print(p)
  
  if (to.return){
    return(p)
  }
  
}
nogrid=theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sort.column=function(x, col) {
  return(x[order(x[,col]),])
}


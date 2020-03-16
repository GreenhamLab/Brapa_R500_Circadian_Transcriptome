

# This function takes in a numeric vector or a factor and a vector of colors
# It returns the color-binned version of the numeric vector
# If Ledge=T, a color ledgend vector is returned
ReturnColorMap<-function(NumVec,ColVec=c("darkblue","skyblue","yellow","red"),nbins=25,Ledge=F)
{ 
  
  if(class(NumVec)=="numeric" | class(NumVec)=="integer")
  {
    # bin the NumVec
    numCut<-cut(NumVec,breaks = nbins)
  }
  else
  {
    numCut<-as.factor(NumVec)
  }
  
  colab<-levels(numCut)
  # dont want more bins than levels
  if(nbins>length(colab))
    nbins<-length(colab)
  
  colFun<-colorRampPalette(colors = ColVec)
  colMap<-colFun(nbins)
  names(colMap)<-colab
  if(Ledge)
  {
    if(class(NumVec)=="numeric")
      nms<-round(seq(min(NumVec),max(NumVec),length.out = length(colMap)),0)
    
    else
      nms<-colab
    
    names(colMap)<-nms
    return(colMap)
  }
  ret<-colMap[numCut]
  return(ret)
}

PhyperOverlap<-function(set1,set2,backset,nlog10=T,FoldEnrich=F,lowerTail=F)
{
  
  if(FoldEnrich)
  {
    # This assumes set1 is the test set and set2 is the reference set
    bRatio<-length(intersect(set2,backset))/length(backset)
    tRatio<-length(intersect(set1,set2))/length(set1)
    #print(c(tRatio,bRatio))
    return(tRatio/bRatio)
  }
  int<-length(intersect(set1,set2))-1
  #print(paste("q:",int,"  m:",length(set1),"  n=",(length(backset)-length(set1)),"  k=",length(set2)))
  ph<-phyper(q=int, m=length(set1), n=(length(backset)-length(set1)), k=length(set2),lower.tail=lowerTail)
  #print(ph)
  if(nlog10)
    return(-log10(ph))
  else
    return(ph)
}


#This function creates a category -> gene mapping based a list of gene annotations
# It takes as input a dataframe with columns:
#     Accession - the gene acession
#     Annotation - the GO ID designating the corresponding annotation
# It outputs a list of 2
# 1- a category -> gene mapping
# 2- a gene -> category mapping
BuildEnrichMaps<-function(inDF)
{
  catMap<-tapply(inDF$Accession,INDEX = inDF$Annotation, function(x) x)
  geneMap<-tapply(inDF$Annotation,INDEX = inDF$Accession, function(x) x)
  
  return(list(catMap=catMap,geneMap=geneMap))
}

#This function takes an annotated ontology list (from BuildEnrichMaps) along with an input vector of accession #s (testvec),
# a vector of code descriptions with bincodes as names (codedesc) and a vector of protein Accessions to use as a reference set (refvec) (all genes by default)
# It returns a list of "EnrichObjs" where each object in the list represents a category that enrichment was calculated for
CalcEnrich<-function(MapBuild,testvec,codedesc,refvec=NULL)
{
  EnrichSet<-setClass(Class="EnrichSet",representation=list(GenesInSet="character",UnannotGenes="character",CatList="matrix"))
  
  if(is.null(refvec))
    refvec<-names(MapBuild[[2]])
  
  if(length(testvec) != length(unique(testvec)))
  {
    print("*** WARNING: Test set contains duplicates which will be removed for the analysis ***")
    testvec<-unique(testvec)
  }
  if(length(setdiff(testvec,refvec))>0)
  {
    print(setdiff(testvec,refvec))
    print("*** Warning: Accessions Found in The Test Set That Are Not Present In the Reference Set  ***")
  }
  
  refvec<-intersect(refvec,names(MapBuild[[2]]))
  
  Unann<-setdiff(testvec,names(MapBuild[[2]]))
  if(is.null(Unann))
    Unann<-""
  
  tmp<-testvec
  testvec<-intersect(testvec,names(MapBuild[[2]]))
  print(paste("Length Of TestVec:",length(testvec)))
  ret<-new("EnrichSet",GenesInSet=testvec,UnannotGenes=Unann)
  if(length(testvec)==0)
  {
    print("*** ERROR: None of the Test Set Accessions Were Found in the Annotation Build ***")
    #print(tmp)
    return(NA)
  }
 
  allcats<-unique(unlist(lapply(testvec,function(x) MapBuild[[2]][[x]])))
  
  print(paste("Number Of Categories in TestVec:",length(allcats)))
  
  retCats<-t(sapply(allcats,CaEn2,MapBuild=MapBuild,testvec=testvec,refvec=refvec,EnrichCat=EnrichCat))
  
  lnames<-sapply(allcats,function(x) paste(x,codedesc[x],sep=" - "))
  row.names(retCats)<-lnames
  
  AdjPval<-p.adjust(as.numeric(retCats[,1]),method="fdr")
  FoldEn<-apply(retCats,1,function(x) (as.numeric(x[5])/as.numeric(x[4]))/(as.numeric(x[3])/as.numeric(x[2])))
  OverUnder<-sapply(FoldEn,function(x) if(x>1) return("Over") else return("Under"))
  retCats<-cbind(retCats[,1,drop=F],AdjPval,FoldEn,OverUnder,retCats[,2:6,drop=F])
  
  if(ncol(retCats)==9)
    colnames(retCats)<-c("P-Value","Adj_P-Value","Fold Enrichment","Enrichment Direction","TotalInRef","AnnInRef","TotalInTest","AnnInTest","GenesWithAnn")
  
  if(nrow(retCats)>1)
    retCats<-retCats[order(as.numeric(retCats[,1])),]
  
  ret@CatList<-retCats
  return(ret)
}
CaEn2<-function(curcat,MapBuild,testvec,refvec=refvec,EnrichCat)
{
  #number of genes in the reference set
  Univ<-length(refvec)
  #number of genes that have curcat annotation in the reference set
  Ucat<-length(intersect(MapBuild[[1]][[curcat]],refvec))
  #number of genes in the test set
  Tset<-length(testvec)
  #number of genes that have curcat annotation in the test set
  Tcat<-length(intersect(MapBuild[[1]][[curcat]],testvec))
  #list of genes with curr ann in test set
  Tgenes<-intersect(MapBuild[[1]][[curcat]],testvec)
  Tgenes<-paste(Tgenes,collapse=" | ")
  #print(paste(Tcat,Ucat,Univ,Tset))
  hypg1<-phyper(q=Tcat-1,m=Ucat,n=(Univ-Ucat),k=Tset,lower.tail=F)
  hypg2<-phyper(q=Tcat,m=Ucat,n=(Univ-Ucat),k=Tset,lower.tail=T)
  hypg<-min(c(hypg1,hypg2))
  
  return(c(hypg,Univ,Ucat,Tset,Tcat,Tgenes))
}

# This function plots a ggplot2-style boxplot
# dList - a list where each element of the list is a numeric vector representing a population to be plotted as one density function
# xAxsNames - T/F print labels for each box on x-axis?
# jitter - numeric, if 0 or greater, actual data points are superimposed on top of the boxes. The jitter number controls the x-axis distribution of the points
# paired - T/F, if T, the corresponding points between two box plots are connected with lines
#   **** paired argument is only implemented for the first two distributions in the dList, and they must be the SAME LENGTH ***
ggboxplot<-function(dList, cols=rep("darkblue",length(dList)), notch=T, notchLines=F, alpha=0.5, main="", ylab="", xlab="", ledge=T, ylim=NULL, xAxsNames=T, xAxsNamesAngle=90, jitter=(-1), paired=F)
{
  theme_set(theme_bw())
  boxLst<-lapply(dList,function(x) boxplot(x,plot=F))
  
  # add whisker horizontle bars
  whiskDF<-data.frame(x1=as.numeric(sapply(1:length(boxLst),function(x) rep(x-0.25,2))),
                      x2=as.numeric(sapply(1:length(boxLst),function(x) rep(x+0.25,2))),
                      y1=as.numeric(sapply(boxLst,function(x) x$stats[c(1,5),1])),
                      y2=as.numeric(sapply(boxLst,function(x) x$stats[c(1,5),1])),
                      nt=as.numeric(sapply(boxLst,function(x) x$conf[,1])),
                      Category=rep(names(boxLst),each=2))
  #print(whiskDF)
  ggdf<-data.frame(vals=unlist(dList, use.names = F),Category = factor(rep(names(dList),times=sapply(dList,length)),levels = names(dList)))
  
  if(is.null(ylim))
  {
    ymin<-min(c(whiskDF$y1,whiskDF$y2))
    ymax<-max(c(whiskDF$y1,whiskDF$y2))
    ylim = c((ymin-ymin*0.05),(ymax+ymax*0.05))
  }
  # if notches are drawn, draw horizontal bars accross in order to compare notches
  if(notchLines)
    ggp<-ggplot(ggdf,aes(x = Category,vals)) + geom_boxplot(outlier.shape = NA, mapping = aes(color=Category, fill=Category), notch=notch) + geom_hline(aes(yintercept=nt, color=Category) , lty=2, size=0.5, data = whiskDF) + coord_cartesian(ylim) + labs(y=ylab, x=xlab) + scale_color_manual(values = cols) + scale_fill_manual(values=adjustcolor(col = cols, alpha.f = alpha))
  else
    ggp<-ggplot(ggdf,aes(x = Category,vals)) + geom_boxplot(outlier.shape = NA, mapping = aes(color=Category, fill=Category), notch=notch) + coord_cartesian(ylim = ylim) + labs(y=ylab, x=xlab) + scale_color_manual(values = cols) + scale_fill_manual(values=adjustcolor(col = cols, alpha.f = alpha)) 
  
  if(!ledge)
    ggp<-ggp+theme(legend.position="none")
  
  if(!(xAxsNames))
    ggp<-ggp+theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
  else
    ggp<-ggp+theme(axis.text.x=element_text(angle = xAxsNamesAngle, hjust = 1))
  
  ggp<-ggp + ggtitle(main) + geom_segment(data=whiskDF, mapping=aes(x=x1,y=y1,xend=x2,yend=y2,color=Category))
  
  if(jitter>=0)
    ggp<-ggp+geom_jitter(mapping = aes(color=Category),position = position_jitter(jitter))
  
  if(paired)
  {
    pairDf<-data.frame(X=rep(1,length(dList[[1]])), Xend=rep(2,length(dList[[2]])), Y=dList[[1]], Yend=dList[[2]])
    pDiff<-pairDf$Y-pairDf$Yend
    pCol<-rep(rgb(1,0,0,0.1),length(pDiff))
    pCol[which(pDiff<0)]<-rgb(0,0,1,0.1)
    pairDf<-cbind(pairDf,pCol)
    ggp<-ggp+geom_segment(mapping = aes(x=X,y=Y,xend=Xend,yend=Yend), color=pCol, data = pairDf, size=1.5)
    ggp<-ggp+guides()
  }
  
  return(ggp)
}

# This function is a slight modification of the GENIE3 algorithm that allows you to input completely distinct target and input matricies
# target.matrix - expression data of all potential target genes in matrix format (samples X genes) so genes as columns (not Rows as usual)
# input.matrix - expression data of all potential TF regulator genes in matrix format (samples X genes) so genes as columns (not Rows as usual)
# K - how to calculate the RandomForest "mtry" parameter
# nb.trees - number of RandomForest trees
# importance.measure - what method to use to calculate feature importance
# seed - random seed to use
# trace - T/F verbose output?
# nThr - number of threads
# ... - additional arguments to RandomForest
RS.Get.Weigth.Matrix<- function(target.matrix, input.matrix, K="sqrt", nb.trees=1000, importance.measure="%IncMSE", seed=NULL, trace=TRUE, nThr=1, ...)  
{
  require(randomForest)
  require(parallel)
  # set random number generator seed if seed is given
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # to be nice, report when parameter importance.measure is not correctly spelled
  if (importance.measure != "IncNodePurity" && importance.measure != "%IncMSE") {
    stop("Parameter importance.measure must be \"IncNodePurity\" or \"%IncMSE\"")
  }
  
  # normalize expression matrix
  target.matrix <- apply(target.matrix, 2, function(x) { (x - mean(x,na.rm = T)) / sd(x,na.rm = T) } )
  input.matrix <- apply(input.matrix, 2, function(x) { (x - mean(x,na.rm =T)) / sd(x,na.rm=T) } )
  # setup weight matrix
  num.samples <- dim(target.matrix)[1]
  num.targets <- dim(target.matrix)[2]
  num.inputs <- dim(input.matrix)[2]
  target.names <- colnames(target.matrix)
  input.names <- colnames(input.matrix)
  #print(input.names)
  
  weight.matrix <- matrix(0.0, nrow=num.targets, ncol=num.inputs)
  rownames(weight.matrix) <- target.names
  colnames(weight.matrix) <- input.names
  
  # set mtry
  if (class(K) == "numeric") {
    mtry <- K
  } else if (K == "sqrt") {
    mtry <- round(sqrt(num.inputs))
  } else if (K == "all") {
    mtry <- num.inputs-1
  } else {
    stop("Parameter K must be \"sqrt\", or \"all\", or an integer")
  }
  if (trace) {
    cat(paste("Starting RF computations with ", nb.trees,
              " trees/target gene,\nand ", mtry,
              " candidate input genes/tree node\n",
              sep=""))
    flush.console()
  }
  
  # compute importances for every target gene
  names(target.names)<-target.names
  
  if(nThr>1)
  {
    clst<-makeCluster(getOption("cl.cores",nThr),outfile="")
    clusterExport(cl = clst, varlist = c("importance","randomForest","RSGWM2","num.targets","target.names","input.matrix","target.matrix","trace","mtry","nb.trees","importance.measure",...),envir = environment())
    imList<-parLapplyLB(cl=clst, X=target.names, function(x) RSGWM2(x,num.targets,target.names,input.matrix,target.matrix,trace,mtry,nb.trees,importance.measure,...))
    stopCluster(cl=clst)
  }
  else
  {
    imList<-lapply(target.names,function(x) RSGWM2(x,num.targets,target.names,input.matrix,target.matrix,trace,mtry,nb.trees,importance.measure,...))
  }
  for(nm in names(imList))
  {
    tcols<-names(imList[[nm]])
    weight.matrix[nm,tcols] <- imList[[nm]]
  }
  # weight.matrix<-sapply(imList,function(x) x)
  
  return(weight.matrix / num.samples)
  #    return(list(Weight=weight.matrix,Model=model.matrix,PedictionCorrelations=cor.vec))
}       
# this is a helper function for RS.Get.Weight.Matrix()
RSGWM2<-function(target.gene.name,num.targets,target.names,input.matrix,target.matrix,trace,mtry,nb.trees,importance.measure,...)
{
  target.gene.idx<-which(target.names==target.gene.name)
  if (trace) 
  {
    cat(paste("Computing gene ", target.gene.idx, "/", num.targets, "\n", sep=""))
    flush.console()
  }
  #target.gene.name <- target.names[target.gene.idx]
  # remove target gene from input genes
  #these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
  temp.input.matrix<-input.matrix
  if(target.gene.name %in% colnames(input.matrix))
  {
    #rmind<-which(colnames(input.matrix)==target.gene.name)
    rmind<-grep(target.gene.name,colnames(input.matrix),fixed=T)
    temp.input.matrix<-input.matrix[,-rmind]
    print(paste("Removing:",target.gene.name,"fom Input| Index Number:",rmind," |  New Dimensions:",dim(temp.input.matrix)[1],"X",dim(temp.input.matrix)[2]))
  }
  x <- temp.input.matrix
  y <- target.matrix[,target.gene.name]
  incSamp<-names(y)[which(!is.na(y))]
  x<-x[incSamp,]
  y<-y[incSamp]
  
  #print(dim(x))
  
  rf <- randomForest(x = x, y = y, mtry=mtry, ntree=nb.trees, keep.forest=F, importance=TRUE,...)
  im <- importance(rf)[,importance.measure]   
  #im.names <- names(im)
  return(im)
}

# this function is used to identify significant edges in a GRN matrix. It takes the output of RS.Get.Weight.Matrix.
# It returns a matrix where the non-significnat edges are set to zero.
#tMat is the actual GRN matrix
#pMat is a permuted version of this
#cutP is a pvalue for the confidence of each edge
CutGRNMat<-function(tMat,pMat,cutPv)
{
  tNum<-as.numeric(tMat)
  pNum<-as.numeric(pMat)
  # calculate left and right-tail cutoffs
  rTcut<-findCut(tVec = tNum, pVec = pNum, type = "RT", fdrVal = cutPv)
  #lTcut<-findCut(tVec = tNum, pVec = pNum, type = "LT", fdrVal = cutPv)
  
  pInd<-which(tMat>rTcut)
  #nInd<-which(tMat<lTcut)
  #inds<-c(pInd,nInd)
  ret<-tMat
  ret[1:length(ret)]<-0
  ret[pInd]<-tMat[pInd]
  return(ret)
}

# This function takes in a set of test values and a set of permuted values
# Each set is analyze assuming the most extreme example is index [1]
# This function scans along the sets until it reaches the specified FDR
# type can be "RT" (right tail) or "LT" (left tail)
findCut<-function(tVec, pVec, fdrVal = 0.05, nCuts=100, type="RT")
{
  #require(parallel)
  # determine the value cutoffs to test the FDR on
  fMin<-min(c(min(tVec,na.rm = T),min(pVec,na.rm = T)))
  fMax<-max(c(max(tVec,na.rm = T),max(pVec,na.rm = T)))
  if(type=="RT")
    fCuts<-seq(fMax,fMin,length.out = nCuts)
  else
    fCuts<-seq(fMin,fMax,length.out = nCuts)
  
  names(fCuts)<-fCuts
  #print(fCuts)
  # apply over fCuts and calculate FDRs
  #clst<-makeCluster(getOption("cl.cores",4), outfile="")
  #clusterExport(cl = clst, varlist = c("pVec","tVec"),envir = environment())
  #tst<-parLapply(cl = clst, X = Br2sTFsTargList[1:8],function(x) TargetOverlapStats(TargList = x,OrthoMat = BrOrtho2_Rep))
  
  if(type=='RT')
  {
    #fdrS<-parSapplyLB(cl = clst, X = fCuts,function(x) (length(which(pVec>x))/length(which(tVec>x)))*(length(pvec)/length(tvec)))
    fdrS<-sapply(X = fCuts,function(x) (length(which(pVec>x))/length(which(tVec>x)))*(length(pVec)/length(tVec)))
  }
  else
  {
    #fdrS<-parSapplyLB(cl = clst, X = fCuts,function(x) (length(which(pVec<x))/length(which(tVec<x)))*(length(pvec)/length(tvec)))
    fdrS<-sapply(X = fCuts,function(x) (length(which(pVec<x))/length(which(tVec<x)))*(length(pVec)/length(tVec)))
  }
  #stopCluster(cl=clst)
  
  #plot(as.numeric(names(fdrS)),fdrS,xlab="Cutoff",ylab="FDR")
  
  if(any(fdrS<fdrVal, na.rm = T))
  {
    retInd<-max(which(fdrS<fdrVal))
    ret<-as.numeric(names(fdrS)[retInd])
    #abline(v=ret,lwd=3,col="red")
  }
  else
  {
    print("No Cutoff could be found for the desired FDR")
    return(NA)
  }
  return(ret)
}

# This function is used to compare one At target group to two Br target groups.
# This function converts the At targets to the corresponding expanded set of Br paralogs
# followed by a hypergeometric test for over enrichment
# This function takes in TargGroups: a list of 3 with names (Ath,Br1 and Br2) each representing a group of target genes
# along with OrthoMat: a matrix of two columns with ortholog mapping from Br (Column 1) to At (Column 2)
# If an OrthoMat is not needed (i.e. no mapping between target groups in required), leave this argument blank and the target groups will be compared directly
TargetOverlapStats<-function(TargGroups,OrthoMat=NA)
{
  # Map the At to Br
  if(!is.na(OrthoMat))
    AllBrGroups<-list(Ath=OrthoMat[which(OrthoMat[,2]%in%TargGroups$Ath),1],Br1=TargGroups$Br1,Br2=TargGroups$Br2)
  else
    AllBrGroups<-TargGroups
      
  # prepare the phyper parameters
  wht<-length(unique(AllBrGroups$Ath))
  #  blk<-length(unique(OrthoMat[,1]))-wht
  draw1<-length(AllBrGroups$Br1)
  draw2<-length(AllBrGroups$Br2)
  whtDraw1<-length(intersect(AllBrGroups$Ath,AllBrGroups$Br1))
  whtDraw2<-length(intersect(AllBrGroups$Ath,AllBrGroups$Br2))
  
  #  phyp1<-phyper(q = whtDraw1-1, m = wht, n = blk, k = draw1, lower.tail = F)
  #  phyp2<-phyper(q = whtDraw2-1, m = wht, n = blk, k = draw2, lower.tail = F)
  
  res<-c(Test1=1,Test2=1, AtLen=wht, Br1Len=draw1, Br2Len=draw2, AtBr1Len=whtDraw1, AtBr2Len=whtDraw2)
  return(res)
}

# this function takes in a target group list (list of 3 with Ath, Br1 and Br2 target groups)
# It creates a permuted group where the Br targets are randomly sampled without replacement
# to create two Br groups of the same size as the original
TargetGroupPerm<-function(targGrp)
{
  AllBrs<-unlist(targGrp[c("Br1","Br2")])
  AllBrs<-sample(AllBrs,length(AllBrs),replace = F)
  ret<-list(Ath=targGrp$Ath, Br1=AllBrs[1:length(targGrp$Br1)], Br2=AllBrs[(length(targGrp$Br1)+1):length(AllBrs)])
  return(ret)
}


ZhangSC<-function(A,M,Vm,clm,q,beta = 1){
  # Ak: Adjacency matrix
  # M: putative number of clusters
  # Vm: eigenvectors that correspond to the M smallest eigenvalues
  # clm: clustering method
  N=nrow(A[[1]])
  K=length(A)
  res=list()
  f.EG=paste("EG_q",q,"b",beta,".Rdata",sep = "")
  print(f.EG)
  print(paste("We have ",N,"genes and ",K," conditions."))
  #### Construct L matrices
  if(missing(Vm)){
    if(file.exists("L.Rdata")){
      print("Loading D,L matrices")
      load("D.Rdata")
      load("L.Rdata")
    }else{
      print("Constructing L matrices")
      D=list()
      L=list()
      for (k in names(A)){
        D[[k]]=diag(rowSums(A[[k]], na.rm = FALSE, dims = 1)
                    ,nrow = N,ncol = N)
        L[[k]]=D[[k]]-A[[k]]
      }
      save(D,file = "D.Rdata")
      save(L,file = "L.Rdata")
    }
    #### Construct C matrix
    print("Constructing C matrix")
    C=matrix(nrow = (K*N),ncol = (K*N),0)
    I=matrix(nrow = (K*N),ncol = (K*N),0)
    r=1;  c=1;
    for (k in names(A)){
      C[r:(r+N-1),r:(r+N-1)]=L[[k]]
      
      # were create the identity matrix (diagonal has ones)
      # then we reverse rows so that we have the anti-diagonal
      # we do this because we have to take the anti diagonal of I but retain the identity matrices at they original form
      anti.In = antidiag(diag(N))
      I[r:(r+N-1),r:(r+N-1)]=anti.In

      r=r+N
      c=c+N
    }
    #Free some space
    rm(L); rm(D);gc();
    print("Free some space")
    print(gc())
    print(f.EG)
    
    # we take the anti-diagonal of I so that the sub-identity matrices are not at the diagonal but on the secondary diagonal 
    I=antidiag(I)
    C=C-(beta*I)
    rm(I);gc();
    if(file.exists(f.EG)){
      print("Loading EigenVectors and EigenValues")
      load(f.EG)
    }else{
      print("Calculating eigenvectors")
      EG<-eigen(C)
      print("Exporting EigenVectors.")
      print(nrow(EG$vectors))
      row.names(EG$vectors)=rep(row.names(A[[1]]),2)
      
      save(EG,file = f.EG)
      # Open a pdf file and save the eigenvalue plot
      #pdf("EigenValues.pdf") 
      plot(EG$values[(N*2):1])
      # Close the pdf file
      #dev.off()
    }
    # select the eigenvectors of the M smallest eigenvalues
    Vm=EG$vectors[,(N-M+1):N]
    res[["eigenvalues"]]=EG$values
    res[["eigenvectorsM"]]=Vm
  }# end of if
  plot(EG$values[(N*2):1])
  if(clm=="kmeans"){
    #### K-means  ####
    print("kmeans clustering.")
    km<-kmeans(Vm,M,iter.max = 1000)
    # kmeans ended returning results
    print("kmeans ended returning results")
    res[["cluster"]] =km$cluster
  }
  if(clm=="hierarchical"){
    print("hierarchical clustering.")
    # Zhang's algorithm says spearman distance
    dis=dist(Vm,method = "euclidean")
    hc1 <- hclust(dis^2, method = "complete")
    print(dim(Vm))
    mycl <- cutree(hc1, h=max(hc1$height/1.5))
    res[["cluster"]] =mycl
  }
  
  return(res)
}

antidiag <- function(M) {
  ##M: matrix
  if(!is.matrix(M)){return(errorCondition("M parameter should be a matrix"))}
  else{
    s = dim(M)
    if(s[1]==0 || s[2]==0){
      return(errorCondition("M should have dimensions."))
    }
    if(s[1]!=s[2]){return(errorCondition("M should be square"))}
    else{N = s[1] }
  }
  for (r in 1:N) {
    M[r,1:N] = M[r,N:1]
  }
  return(M)
}


MergedCoexpGraph <- function(Cs1,Cs2,s1,s2,short.nm,logFC,f.prefx,module.num) {
  library("igraph")
  print(head(short.nm))
  f.txt=paste(f.prefx,"_NODES.txt",sep = "")
  SYMBOLS=rownames(Cs1)
  #Unconnected nodes in both networks are excluded
  strong_couples=which( (Cs1>0 | Cs2>0) & lower.tri(Cs1), arr.ind = TRUE)
  ind.strg.cpls=(Cs1>0 | Cs2>0) & lower.tri(Cs1)
  #n.n=nrow(strong_couples)
  print(paste("<MergedCoexpGraph> #Nodes: ",nrow(Cs1)))
  print(paste("<MergedCoexpGraph> Strong couples (#edges): ",nrow(strong_couples)))
  stages=paste(c("",s1)[(Cs1>0)+1],
               c("",s2)[(Cs2>0)+1],
               sep = ";",collapse = NULL
  )
  if(nrow(strong_couples)>0){
    wlogFC=log2( (Cs2+1)/(Cs1+1) )
    ##### Create graph from dataframe ####
    status=c("yes","no")
    temp=(wlogFC[ind.strg.cpls]<0)+1
    cor_decr=status[temp]
    GR.DF<-data.frame(source=SYMBOLS[strong_couples[,1]],
                      target=SYMBOLS[strong_couples[,2]],
                      stage=stages[ind.strg.cpls],
                      logw=wlogFC[ind.strg.cpls],
                      cor_decr=cor_decr)
    glist=unique(c(SYMBOLS[strong_couples[,1]],SYMBOLS[strong_couples[,2]]))
    write.table(paste(glist),file = f.txt,sep = "\t")
    #,    cor_decr=status[ (wlogFC[ind.strg.cpls]<0)*1 ] 
    GR<-graph.data.frame(GR.DF, directed=FALSE, vertices=NULL)
    
    #### Save some Vector features ####
    V<-as.matrix(V(GR))
    indx.SYMBOLS<-simplify2array(lapply(row.names(V), 
                                        function(x) which(SYMBOLS==x)))
    for (i in colnames(logFC)) {
      GR<-set.vertex.attribute(GR, name=i,index = V(GR),
                               value=logFC[indx.SYMBOLS,i])
    }

    GR<-set.vertex.attribute(GR, name="coal.name",index = V(GR),
                             value=as.matrix(short.nm[indx.SYMBOLS]) )
    GR<-set.vertex.attribute(GR, name="type",index = V(GR),
                             value=as.matrix(rep("protein",nrow(V))) )

    v.names = igraph::get.vertex.attribute(GR,"coalname")

    GR<-set.vertex.attribute(GR, name="module",index = V(GR),
                             value = as.matrix(module.num[v.names]) )#as.matrix(rep(module.num,nrow(V))))
    
    write_graph(GR, paste(f.prefx,".gml",sep = ""), format = "gml")
  }
  
 
  #library("varbvs")
}

Modules <- function(M,mthd,limma,q,beta = 1) {
  # CX: list of NxN matrices of coexpression values
  # CX.A: list of NxN adjacency matrices of coexpression values
  # number of modules
  get('CX', envir=.GlobalEnv)
  get('CX.A', envir=.GlobalEnv)

  conditions=names(CX.A)
  N=nrow(CX.A[[1]]);print(paste("Number of proteins:",N))
  
  SYMBOLS=row.names(CX.A[[1]])
  ext=paste(conditions,collapse = "to")
  
  merged.nodes=matrix(nrow = N,ncol = 1,FALSE)
  node.module=matrix(nrow = N,ncol = 1,-1) # Saves the number of module of the node
  row.names(node.module) = SYMBOLS
  ## File names
  f1 = paste("Eigenvalues_",ext,"_Ms",M,"_q",q,"_beta",beta,".txt",sep = "")
  f0 = paste(ext,"_ZhangSC_Ms",M,"_",mthd,"_q",q,"_beta",beta,".Rdata",sep = "")
  f2 = paste(ext,"_clusters_Ms",M,"_",mthd,"_q",q,"_beta",beta,".txt",sep = "")
  f3 = paste(ext,"_modules_Ms",M,"_",mthd,"_q",q,"_beta",beta,".txt",sep = "")
  
  ### Spectral Clustering / Save Results
  cls<-ZhangSC(CX.A,as.integer(M),clm=mthd,q=q, beta = beta)
  save(cls,file=f0)
  #### Save Eigenvalues in txt file
  write.table(cls$eigenvalues[(2*N):1],file = f1,sep = "\t")
  #### RSH file contains node clustering and logFC values from limma
  RSH.cls=matrix(cls$cluster,nrow = N, ncol = 2 )
  row.names(RSH.cls)=SYMBOLS
  colnames(RSH.cls)=conditions
    
  RSH=cbind(RSH.cls,as.matrix(limma$table))
  write.table(RSH,file = f2,sep = "\t")
  #########################
  n.modules=list()
  ## Iterate across clusters
  for (m in 1:M) {
    
    # based on Zhang et al 2018; a cluster and its counterpart may include different nodes. 
    # We take the union of them
    index=(RSH[,1]==m | RSH[,2]==m)

    n.nodes=sum(index)     
    n.modules[[as.character(m)]]=n.nodes
    
    if(n.nodes>8 & n.nodes<800){
      print(paste("###Module: ",m," Nodes: ",n.nodes," <Modules>"))

      merged.nodes[index] = TRUE
      node.module[index] = m
      
      MergedCoexpGraph(CX[[1]][index,index],
                       CX[[2]][index,index],
                       conditions[1],conditions[2],
                       RSH[index,grep("combo.name",colnames(RSH))],RSH[index,3:ncol(RSH)],
                       paste("M",m,"_",ext,"_q",q,"_b",beta,sep = ""),node.module)
      
      }
   }# module loop
  write.table(as.matrix(n.modules),file = f3,sep = "\t")
  print(node.module);
  ## Merged nodes of all Modules
  union.indx=which(merged.nodes)
  MergedCoexpGraph(CX[[1]][union.indx,union.indx],
                   CX[[2]][union.indx,union.indx],
                   conditions[1],conditions[2],
                   RSH[merged.nodes,grep("combo.name",colnames(RSH))],limma$logFC[merged.nodes,],
                   paste(ext,"_M",M,"_q",q,"_b",beta,sep = ""),node.module)
}

Map2limma <- function(limma.DEGs,SYMBOLS,pv.col) {
  # limma.DEGs: results of limma comparative analysis
  # SYMBOLS:    Adjacency matrix gene names
  # pv.col:     p-value and logFC index column  

  sel.genes = SYMBOLS[SYMBOLS %in% row.names(limma.DEGs)]
  sel.cols = grep(pv.col,colnames(limma.DEGs))
  df= limma.DEGs[sel.genes,sel.cols]
  
  print(sel.cols)
  #print(sel.genes)
  #print(dim(df))

  df = data.frame(df)
  df$combo.name=sel.genes

  degs=rownames(df)[ df[,grep("adj.P.val",colnames(limma.DEGs)[sel.cols])]<=0.05 ]
  print(grep("adj.P.val",colnames(limma.DEGs)[sel.cols]))

  print( paste("We have ",length(degs)," DEGS in ",pv.col) )
  
  res=list(degs=degs,table=df)
  return(res)
}

eigenvalue.step <- function(EGV) {
  
  k=diff(EGV)
  h = k[1:length(k)-1]/k[2:length(k)]
  
  return( which( h>(2*mean(h)) ) )
}
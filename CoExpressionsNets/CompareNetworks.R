Common.Modules <- function(conditions,DEGS, beta = 1, q = 0.95, n.modules = 10) {
  ##  Loads gene co-expression adjacency matrices saves in .Rdata files
  ##  The co-expression .Rdata files are named as "COX_6U133_{condition}.Rdata"
  ##  where condition suffix is contained in the condition input variable
  ##
  ##  conditions: a vector of two strings, the names of the conditions, e.g. c("MCL","Bcell")
  ##  q  = quantalile percentage value, based on wchich we flatten the coexpression network, 
  ##        i.e. set coexpression relations to zero
  ##  n.modules: set the number of co-expression modules to be searched
  ##  beta: the beta value as described in Shuqin Zhang (2015) (doi:10.1109/TCBB.2015.2396073)
  
  CX=list()
  print(conditions)
  for(m in conditions){
    load(paste("COX_6U133_",m,".Rdata",sep=""))
    CX[[m]]=cur.COX
  }
  transcripts = nrow(cur.COX)
  print(paste("Transcipts:",transcripts))
  rm(cur.COX)
  #-- Flatten weak coexpression relations
  for (k in names(CX)){
    qv<-quantile(CX[[k]],q)
    CX[[k]][CX[[k]]<qv]=0
  }
  #---------------------  adjacency Matrix
  CX.A=list();
  for (k in names(CX)){
    CX.A[[k]] = CX[[k]]
    CX.A[[k]][CX.A[[k]]>0] = 1 
  }
  
  SYMBOLS=rownames(CX.A[[1]])
  #-- save CX (coexpression matrices) to global ENviroment (instead of using it as function input)
  assign("CX", CX, envir = .GlobalEnv)
  assign("CX.A", CX.A, envir = .GlobalEnv)
  
  parent.dir =gsub("/CoExpressionsNets","",getwd())
  cmp= c(conditions,n.modules,paste(conditions[1],"_",conditions[2],sep = ""))
  print(cmp)
  
  #--- create output directory in case it doesn't exist
  subDir=file.path("Modules",paste(cmp[1],"_",cmp[2],"_",cmp[3],"_q",q,sep = ""))
    
  if (file.exists(file.path(parent.dir,"CoExpressionsNets", subDir))){
    setwd(file.path(parent.dir,"CoExpressionsNets", subDir))
  }
  else{
    cur.dir = file.path(parent.dir,"CoExpressionsNets", subDir)
    dir.create(cur.dir)
    setwd(cur.dir)
    
    print(paste("Creating directory",cur.dir))
  }
  
  
  #--- calculate C table, Shuqin Zhang (2015) (doi:10.1109/TCBB.2015.2396073)
  source("Zhang_spectral_clustering.R")
  
  limma.map=Map2limma(DEGS,SYMBOLS,cmp[4])
    
  Modules(cmp[3],
          mthd="hierarchical",
          limma.map,q,beta)
    setwd(file.path(parent.dir,"CoExpressionsNets"))
}

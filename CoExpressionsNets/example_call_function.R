########### FOR LOOP MODULES #############

source("CompareNetworks.R"); 
load("DEGS.Rdata")
Common.Modules(conditions = c("MCL","Bcell"),limma.DEGs,beta = 1,q=0.95, n.modules = 55)
Common.Modules(conditions = c("MCLe","Bcell"),limma.DEGs,beta = 2,q=0.95, n.modules = 55)

############### EIGEN GAP ################
#--------- # Modules, decide based on eigengap
load("Modules/ALKminus_Naive_100_q0.95/EG_q0.95b1.Rdata")
source("Zhang_spectral_clustering.R");gap.pos = eigenvalue.step(EG$values);
print(gap.pos)
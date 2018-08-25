############################################################################################
############################################################################################
# R functions corresponding to the manuscript entitled "Array Testing with Multiplex Assays"                                                                          
############################################################################################
############################################################################################

############################################################################################
# Note: some functions are written in standalone C++ file for faster calculation.  
#       Latest update of R (3.4.2), Rtools (Rtools34), R packages "Rcpp" and "RcppArmadillo"  
#       are required.
############################################################################################

library(Rcpp)
library(RcppArmadillo)
sourceCpp("Hierarchical.cpp")
sourceCpp("Array.cpp")

############################################################################################
#       1.  Efficiency and Accuracy calculation for hierarchical and Array algorithms
############################################################################################

############################################################################################
# Function 1: EFF_HIER(p,S,SE,SP,ns)
# Purpose: This function calculates the expected number of tests per-individual described in 
#          Section 3.1 of the manuscript "Hierarchical group testing for multiple infections".
# 
# Input:
# p = Vector of cell probabilities: c(p00,p10,p01,p11)
# S = Number of stages 
# SE = Sensitivity matrix; Dimension: S*2 (note: row=stage, column=infection)
# SP = Specificity matrix; Dimension: S*2
# ns = Vector of pool sizes at stage 1,2,...,S: c(n1,n2, ... ,nS)
#
# Output: Expected number of tests per individual E(T)/n1

############################################################################################
#              Example: 
Prevalence=c(0.95,0.03,0.01,0.01)
stagenum=3
sensitivity=matrix(0.95,stagenum,2)
specificity=matrix(0.99,stagenum,2)
groupsize=c(9,3,1)
EFF_HIER(p=Prevalence,S=stagenum,SE=sensitivity,SP=specificity,ns=groupsize)
############################################################################################




###############################################################################################
# Function 2: ACCU_HIER(p,S,SE,SP,ns)
# Purpose: This function calculates the classification accuracy measures described in 
#           Section 3.2 of the manuscript "Hierarchical group testing for multiple infections".
#
# Input:
# p = Vector of cell probabilities: c(p00,p10,p01,p11)
# S = Number of stages
# SE = Sensitivity matrix; Dimension: S*2 (note: row=stage, column=infection)
# SP = Specificity matrix; Dimension: S*2
# ns = Vector of pool sizes at stage 1,2,...,S: c(n1,n2, ... ,nS)
#
# Output: PSE1/PSE2: Pooling sensitivities for first/second infection
#         PSP1/PSP2: Pooling specificities for first/second infection
#         PPV1/PPV2: Pooling positive predictive values for first/second infection
#         NPV1/NPV2: Pooling negative predictive values for first/second infection
#         ACCU: Expected number of correct classification per individual: E(C)/n1

###############################################################################################
#              Example: 
Prevalence=c(0.999,0.0004,0.0004,0.0002)
stagenum=6
sensitivity=matrix(0.95,stagenum,2)
specificity=matrix(0.99,stagenum,2)
groupsize=c(96,48,24,12,4,1)
ACCU_HIER(p=Prevalence,S=stagenum,SE=sensitivity,SP=specificity,ns=groupsize)
###############################################################################################



###########################################################################################################
# Function 3: ARRAY(p,SE,SP,n)
# Purpose: This function calculates the expected number of tests per-individual and accuracies
#          for algorithms AT and ATM described in Section 3 of the manuscript "Array Testing 
#          with Multiplex Assays".
# 
# Input:
# p = Vector of cell probabilities: c(p00,p10,p01,p11)
# SE = Sensitivity vector; c(se1,se2)
# SP = Specificity vector; c(sp1,sp2)
# n = row/column size
#
# Output: ET_AT/ET_ATM: Expected number of tests per-individual for algorithm AT/ATM
#         PSE1_AT/PSE1_ATM: Pooling sensitivities of the first infection for algorithm AT/ATM 
#         PSE2_AT/PSE2_ATM: Pooling sensitivities of the second infection for algorithm AT/ATM 
#         PSP1_AT/PSP1_ATM: Pooling specificities of the first infection for algorithm AT/ATM 
#         PSP2_AT/PSP2_ATM: Pooling specificities of the second infection for algorithm AT/ATM 
#         PPV1_AT/PPV1_ATM: Pooling positive predictive values of the first infection for algorithm AT/ATM 
#         PPV2_AT/PPV2_ATM: Pooling positive predictive values of the second infection for algorithm AT/ATM
#         NPV1_AT/NPV1_ATM: Pooling negative predictive values of the first infection for algorithm AT/ATM 
#         NPV2_AT/NPV2_ATM: Pooling negative predictive values of the second infection for algorithm AT/ATM 
###########################################################################################################

###############################################################################################
#              Example: 
Prevalence=c(0.999,0.0004,0.0004,0.0002)
sensitivity=c(0.95, 0.95)
specificity=c(0.99, 0.99)
row_col=5
ARRAY(p=Prevalence,SE=sensitivity,SP=specificity,n=row_col)
###############################################################################################




##################################################################################################
#       2.  Optimal configuration and corresponding characteristics for algorithm H2, H3, AT, ATM
##################################################################################################

#################################################################################################################
# Function 4: OPT.H_VS_A(p,SE,SP,MAMPS,n_at,n_atm)
# Purpose: This function calculates the optimal configuration and the corresponding operating characteristics
#          base on minimizing expected number of tests per-individual for each of H2, H3, AT, and ATM algorithm.
#
# Input:
# p = Vector of cell probabilities: c(p00,p10,p01,p11)
# SE = Sensitivity vector; c(se1,se2)
# SP = Specificity vector; c(sp1,sp2)
# MAMPS = Maximal allowable master pool size for hierarchical algorithms
# n_at = Maximal row/column size for algorithm AT
# n_atm = Maximal row/column size for algorithm ATM

OPT.H_VS_A<-function(p,SE,SP,MAMPS,n_at,n_atm){
  SE_h = matrix(SE,20,2,byrow=T)
  SP_h = matrix(SP,20,2,byrow=T)
  H <- OPT(p,SE_h,SP_h,MAMPS,"minimize","optimal")
  H_size = H$PoolSize[1:2,1:3]; H_size [1,3] = NA; #Hierarchical optimal group sizes for H2 and H3
  colnames(H_size)=c("Stage 1","Stage 2", "Stage 3")
  rownames(H_size)=c("H2 pool sizes", "H3 pool sizes")
  OprChar <- H$OperatingCharacteristics[1:2,1:9]   #Hierarchical accuracy measures under the optimal designs for H2 and H3
  colnames(OprChar)=c("E(T)/n1","PSe1", "PSe2", "PSp1", "PSp2", "PPV1","PPV2","NPV1","NPV2")
  rownames(OprChar)=c("H2", "H3")
  array_res<-sapply(seq(2,max(n_at, n_atm),1),ARRAY,p=p,SE=SE,SP=SP) 
  AT_size<-which(unlist(array_res[1,1:(n_at-1)])==min(unlist(array_res[1,1:(n_at-1)])))+1
  AT_OprChar<-rbind(c(unlist(array_res[1:9,AT_size-1])))
  
  ATM_size<-which(unlist(array_res[10,1:(n_atm-1)])==min(unlist(array_res[10,1:(n_atm-1)])))+1
  ATM_OprChar<-rbind(c(unlist(array_res[10:18,ATM_size-1])))
  
  array_size = rbind(AT_size, ATM_size)
  rownames(array_size) <- c("AT", "ATM")
  colnames(array_size) <- "Row/Column size"
  
  array_OprChar = rbind(AT_OprChar, ATM_OprChar)
  colnames(array_OprChar)=c("E(T)/n^2","PSe1", "PSe2", "PSp1", "PSp2", "PPV1","PPV2","NPV1","NPV2")
  rownames(array_OprChar)=c("AT", "ATM")
  
  output<-list("Hierarchical optimal group sizes" = H_size, "Hierarchical accuracy measures under optimal algorithms" = OprChar, 
               "Array optimal row/column size" = array_size, "Array accuracy measures under optimal design" = array_OprChar)
  output
}

###############################################################################################################
#              Example: 
Prevalence=c(0.95,0.02,0.02,0.01)
sensitivity=c(0.95,0.95)                                                             
specificity=c(0.99,0.99)     
maximal=100
n_at=20
n_atm=10
OPT.H_VS_A(p=Prevalence,SE=sensitivity,SP=specificity,MAMPS=maximal, n_at=n_at, n_atm=n_atm)

##############################################################################################################



##########################################################################################################
#                                       3.  Simulation
##########################################################################################################

#############################################################################################################
# Function 5: HIER_SIM(size, p, S, SE, SP, ns, iter)
# Purpose: Conduct Monte Carlo simulations for given infection prevalences, assay sensitivities/specificities
#          and hierarchical algorithm design.  
#
#
# Input:
# size = Number of individuals
# p = Vector of cell probabilities: c(p00,p10,p01,p11)
# S = Number of stage
# SE = Sensitivity vector: c(se1,se2)
# SP = Specificity vector: c(sp1, sp2)
# ns = Vector of pool sizes at stage 1,2,...,S: c(n1,n2, ... ,nS)
# iter = Number of MC replications
#
# Output: Operating characteristics averaged over the number of replications specified in iter

###################################################################################################
#                       Example: 
size = 10000
Prevalence = c(0.90,0.04,0.04,0.02)
S=2
sensitivity = c(0.95,0.95)
specificity = c(0.99,0.99)
ns = c(3,1)
B = 1000
res<-HIER_SIM(size=size, p=Prevalence, S=S, SE=sensitivity, SP=specificity, ns=ns, iter=B)
mean(res$'Number of tests')/size; 
mean(res$'PSE1')
mean(res$'PSE2')
mean(res$'PSP1')
mean(res$'PSP2')
mean(res$'PPV1')
mean(res$'PPV2')
mean(res$'NPV1')
mean(res$'NPV2')
boxplot(res$'Number of tests'/size,range=0,xaxes=F,cex.axis=1.5,col="grey",xlab = paste("H",S,sep=""), 
        ylab = "Number of tests per-individual", cex.lab=1.5)

#############################################################################################################
# Function 6: AT_SIM(size,p,SE,SP,row,col,iter)
# Purpose: Conduct Monte Carlo simulations for given infection prevalences, assay sensitivities/specificities
#          and row and column size for array testing algorithm AT.  
#
#
# Input:
# size = Number of individuals
# p = Vector of cell probabilities: c(p00,p10,p01,p11)
# SE = Sensitivity vector: c(se1,se2)
# SP = Specificity vector: c(sp1, sp2)
# row = Row size
# col = Column size
# iter = Number of MC replications
#
# Output: Operating characteristics averaged over the number of replications specified in iter

###################################################################################################
#                       Example: 
size = 10000
Prevalence = c(0.90,0.04,0.04,0.02)
sensitivity = c(0.95,0.95)
specificity = c(0.99,0.99)
row = 4
col = 5
B = 1000
res<-AT_SIM(size=size,p=Prevalence,SE=sensitivity,SP=specificity,row=row,col=col,iter=B)
mean(res$'Number of tests')/size; 
mean(res$'PSE1')
mean(res$'PSE2')
mean(res$'PSP1')
mean(res$'PSP2')
mean(res$'PPV1')
mean(res$'PPV2')
mean(res$'NPV1')
mean(res$'NPV2')
boxplot(res$'Number of tests'/size,range=0,xaxes=F,cex.axis=1.5,col="grey", xlab = "AT", 
        ylab = "Number of tests per-individual", cex.lab=1.5)

#############################################################################################################
# Function 7: ATM_SIM(size,p,SE,SP,row,col,iter)
# Purpose: Conduct Monte Carlo simulations for given infection prevalences, assay sensitivities/specificities
#          and row and column size for array testing algorithm ATM.  
#
#
# Input:
# size = Number of individuals
# p = Vector of cell probabilities: c(p00,p10,p01,p11)
# SE = Sensitivity vector: c(se1,se2)
# SP = Specificity vector: c(sp1, sp2)
# row = Row size
# col = Column size
# iter = Number of MC replications
#
# Output: Operating characteristics averaged over the number of replications specified in iter

###################################################################################################
#                       Example: 
size = 10000
Prevalence = c(0.90,0.04,0.04,0.02)
sensitivity = c(0.95,0.95)
specificity = c(0.99,0.99)
row = 4
col = 5
B = 1000
res<-ATM_SIM(size=size,p=Prevalence,SE=sensitivity,SP=specificity,row=row,col=col,iter=B)
mean(res$'Number of tests')/size; 
mean(res$'PSE1')
mean(res$'PSE2')
mean(res$'PSP1')
mean(res$'PSP2')
mean(res$'PPV1')
mean(res$'PPV2')
mean(res$'NPV1')
mean(res$'NPV2')
boxplot(res$'Number of tests'/size,range=0,xaxes=F,cex.axis=1.5,col="grey", xlab = "ATM", 
        ylab = "Number of tests per-individual", cex.lab=1.5)


###################################################################################################################################
# Function 8: AT_3_SIM(size,p,SE,SP,row,col,iter)
# Purpose: For three infections, conduct Monte Carlo simulations for given infection prevalences, assay sensitivities/specificities
#          and row and column size for array testing algorithm AT.  
#
#
# Input:
# size = Number of individuals
# p = Vector of cell probabilities: c(p000,p100,p010,p001,p110,p101,p011,p111)
# SE = Sensitivity vector: c(se1,se2,se3)
# SP = Specificity vector: c(sp1,sp2,sp3)
# row = Row size
# col = Column size
# iter = Number of MC replications
#
# Output: Operating characteristics averaged over the number of replications specified in iter

##################################################################################################################################
#                       Example: 
size = 10000
Prevalence = c(0.900,0.02,0.02,0.02,0.01,0.01,0.01,0.01)
sensitivity = c(0.95,0.95,0.95)
specificity = c(0.99,0.99,0.99)
row = 10
col = 10
B = 1000
res<-AT_3_SIM(size=size,p=Prevalence,SE=sensitivity,SP=specificity,row=row,col=col,iter=B)
mean(res$'Number of tests')/size; 
mean(res$'PSE1')
mean(res$'PSE2')
mean(res$'PSP1')
mean(res$'PSP2')
mean(res$'PPV1')
mean(res$'PPV2')
mean(res$'NPV1')
mean(res$'NPV2')
boxplot(res$'Number of tests'/size,range=0,xaxes=F,cex.axis=1.5,col="grey", xlab = "AT", 
        ylab = "Number of tests per-individual", cex.lab=1.5)


###################################################################################################################################
# Function 9: ATM_3_SIM(size,p,SE,SP,row,col,iter)
# Purpose: For three infections, conduct Monte Carlo simulations for given infection prevalences, assay sensitivities/specificities
#          and row and column size for array testing algorithm ATM.  
#
#
# Input:
# size = Number of individuals
# p = Vector of cell probabilities: c(p000,p100,p010,p001,p110,p101,p011,p111)
# SE = Sensitivity vector: c(se1,se2,se3)
# SP = Specificity vector: c(sp1,sp2,sp3)
# row = Row size
# col = Column size
# iter = Number of MC replications
#
# Output: Operating characteristics averaged over the number of replications specified in iter

##################################################################################################################################
#                       Example: 
size = 10000
Prevalence = c(0.900,0.02,0.02,0.02,0.01,0.01,0.01,0.01)
sensitivity = c(0.95,0.95,0.95)
specificity = c(0.99,0.99,0.99)
row = 10
col = 10
B = 1000
res<-ATM_3_SIM(size=size,p=Prevalence,SE=sensitivity,SP=specificity,row=row,col=col,iter=B)
mean(res$'Number of tests')/size; 
mean(res$'PSE1')
mean(res$'PSE2')
mean(res$'PSP1')
mean(res$'PSP2')
mean(res$'PPV1')
mean(res$'PPV2')
mean(res$'NPV1')
mean(res$'NPV2')
boxplot(res$'Number of tests'/size,range=0,xaxes=F,cex.axis=1.5,col="grey", xlab = "AT", 
        ylab = "Number of tests per-individual", cex.lab=1.5)

############################################################################################
############################################################################################
# This R program reproduces the results (Tables 2 and 3) in Section 5 of the manuscript 
#                 "Array testing with multiplex assays."                                                                          
############################################################################################
############################################################################################

############################################################################################
# Please download all the R and cpp files from: https://github.com/harrindy/multiplex.
# Please make sure that R can access the Rcpp files: Array.cpp and Hierarchical.cpp, 
#                                and the data file: IowaSHLdata.csv.
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

#################################################################################################################
# Function: OPT.H_VS_A(p,SE,SP,MAMPS,n_at,n_atm)
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

############################################################################################
# Real data analysis
############################################################################################
# Read the data.
library(utils)
SHL.data=read.csv("IowaSHLdata.csv",col.names =c("Stratum","CT_N,NG_N","CT_P,NG_N","CT_N,NG_P","CT_P,NG_P","Total",
                                                 "Se1","Sp1","Se2","Sp2"))
############################################################################################
# Reproduce Table 2 in the paper.
Table2=SHL.data[,1:6]
############################################################################################
# Find the optimal design of each of H2, H3, AT, and ATM algorithm for the stratum "Female Swab."
ind=1
Prevalence = round(unlist(SHL.data[ind,3:5]/SHL.data[ind,6]),3)
Prevalence=c(1-sum(Prevalence),Prevalence)
sensitivity = unlist(SHL.data[ind,c(7,9)])
specificity = unlist(SHL.data[ind,c(8,10)])
maximal=10
n_at=20
n_atm=10
opt.FS=OPT.H_VS_A(p=Prevalence,SE=sensitivity,SP=specificity,MAMPS=maximal, n_at=n_at, n_atm=n_atm)

# Save results to Table2.
Table2$hat_p00[ind]=Prevalence[1]
Table2$hat_p10[ind]=Prevalence[2]
Table2$hat_p01[ind]=Prevalence[3]
Table2$hat_p11[ind]=Prevalence[4]
Table2=cbind(Table2,SHL.data[,7:10])
Table2$Optimal_H2[ind]=paste("H2(",opt.FS[[1]][1,1],":",opt.FS[[1]][1,2],")",sep="")
Table2$Optimal_H3[ind]=paste("H3(",opt.FS[[1]][2,1],":",opt.FS[[1]][2,2],":",opt.FS[[1]][2,3],")",sep="")
Table2$Optimal_AT[ind]=paste("AT(",opt.FS[[3]][1],"*",opt.FS[[3]][1],")",sep="")

############################################################################################
# Find the optimal design of each of H2, H3, AT, and ATM algorithm for the stratum "Female Urine."
ind=2
Prevalence = round(unlist(SHL.data[ind,3:5]/SHL.data[ind,6]),3)
Prevalence=c(1-sum(Prevalence),Prevalence)
sensitivity = unlist(SHL.data[ind,c(7,9)])
specificity = unlist(SHL.data[ind,c(8,10)])
maximal=10
n_at=20
n_atm=10
opt.FU=OPT.H_VS_A(p=Prevalence,SE=sensitivity,SP=specificity,MAMPS=maximal, n_at=n_at, n_atm=n_atm)

# Save results to Table2.
Table2$Optimal_H2[ind]=paste("H2(",opt.FU[[1]][1,1],":",opt.FU[[1]][1,2],")",sep="")
Table2$Optimal_H3[ind]=paste("H3(",opt.FU[[1]][2,1],":",opt.FU[[1]][2,2],":",opt.FU[[1]][2,3],")",sep="")
Table2$Optimal_AT[ind]=paste("AT(",opt.FU[[3]][1],"*",opt.FU[[3]][1],")",sep="")
Table2$hat_p00[ind]=Prevalence[1]
Table2$hat_p10[ind]=Prevalence[2]
Table2$hat_p01[ind]=Prevalence[3]
Table2$hat_p11[ind]=Prevalence[4]
############################################################################################
# Find the optimal design of each of H2, H3, AT, and ATM algorithm for the stratum "Male Swab."
ind=3
Prevalence = round(unlist(SHL.data[ind,3:5]/SHL.data[ind,6]),3)
Prevalence=c(1-sum(Prevalence),Prevalence)
sensitivity = unlist(SHL.data[ind,c(7,9)])
specificity = unlist(SHL.data[ind,c(8,10)])
maximal=10
n_at=20
n_atm=10
opt.MS=OPT.H_VS_A(p=Prevalence,SE=sensitivity,SP=specificity,MAMPS=maximal, n_at=n_at, n_atm=n_atm)

# Save results to Table2.
Table2$Optimal_H2[ind]=paste("H2(",opt.MS[[1]][1,1],":",opt.MS[[1]][1,2],")",sep="")
Table2$Optimal_H3[ind]=paste("H3(",opt.MS[[1]][2,1],":",opt.MS[[1]][2,2],":",opt.MS[[1]][2,3],")",sep="")
Table2$Optimal_AT[ind]=paste("AT(",opt.MS[[3]][1],"*",opt.MS[[3]][1],")",sep="")
Table2$hat_p00[ind]=Prevalence[1]
Table2$hat_p10[ind]=Prevalence[2]
Table2$hat_p01[ind]=Prevalence[3]
Table2$hat_p11[ind]=Prevalence[4]
############################################################################################
# Find the optimal design of each of H2, H3, AT, and ATM algorithm for the stratum "Male Urine."
ind=4
Prevalence = round(unlist(SHL.data[ind,3:5]/SHL.data[ind,6]),3)
Prevalence=c(1-sum(Prevalence),Prevalence)
sensitivity = unlist(SHL.data[ind,c(7,9)])
specificity = unlist(SHL.data[ind,c(8,10)])
maximal=10
n_at=20
n_atm=10
opt.MU=OPT.H_VS_A(p=Prevalence,SE=sensitivity,SP=specificity,MAMPS=maximal, n_at=n_at, n_atm=n_atm)

# Save results to Table2.
Table2$Optimal_H2[ind]=paste("H2(",opt.MU[[1]][1,1],":",opt.MU[[1]][1,2],")",sep="")
Table2$Optimal_H3[ind]=paste("H3(",opt.MU[[1]][2,1],":",opt.MU[[1]][2,2],":",opt.MU[[1]][2,3],")",sep="")
Table2$Optimal_AT[ind]=paste("AT(",opt.MU[[3]][1],"*",opt.MU[[3]][1],")",sep="")
Table2$hat_p00[ind]=Prevalence[1]
Table2$hat_p10[ind]=Prevalence[2]
Table2$hat_p01[ind]=Prevalence[3]
Table2$hat_p11[ind]=Prevalence[4]
##################################################
# Reproducing Table 2 in the paper is finished.
Table2
##################################################






############################################################################################
# Reproduce Table 3 in the paper and Figure E.1 in the supplementary.
Table3=SHL.data[rep(1:4,rep(3,4)),c(1,6)]
row.names(Table3)=1:12
############################################################################################
# Reproduce results related to the stratum "Female Swab."
ind=1
Prevalence = round(unlist(SHL.data[ind,3:5]/SHL.data[ind,6]),3)
Prevalence=c(1-sum(Prevalence),Prevalence)
sensitivity = unlist(SHL.data[ind,c(7,9)])
specificity = unlist(SHL.data[ind,c(8,10)])
size = unlist(SHL.data[ind,6])

opt=opt.FS

Table3$Algorithm[(ind-1)*3+1]=paste("H2(",opt[[1]][1,1],":",opt[[1]][1,2],")",sep="")
Table3$Algorithm[(ind-1)*3+2]=paste("AT(",opt[[3]][1],"*",opt[[3]][1],")",sep="")
Table3$Algorithm[(ind-1)*3+3]=paste("H3(",opt[[1]][2,1],":",opt[[1]][2,2],":",opt[[1]][2,3],")",sep="")

ns = c(opt[[1]][1,1],opt[[1]][1,2])
set.seed(10)
res<-HIER_SIM(size=size, p=Prevalence, S=length(ns), SE=sensitivity, SP=specificity, ns=ns, iter=5000)

FS.H2=c(round(mean(res$'Number of tests')),round(sd(res$'Number of tests')),
        formatC(mean(res$'Number of tests')/size,digit=3,format="f"),
        formatC(mean(res$'PSE1'),digit=3,format="f"),
        formatC(mean(res$'PSE2'),digit=3,format="f"),
        formatC(mean(res$'PSP1'),digit=3,format="f"),
        formatC(mean(res$'PSP2'),digit=3,format="f"),
        formatC(mean(res$'PPV1'),digit=3,format="f"),
        formatC(mean(res$'PPV2'),digit=3,format="f"),
        formatC(mean(res$'NPV1'),digit=3,format="f"),
        formatC(mean(res$'NPV2'),digit=3,format="f"))
FS.H2
res.H2=res

ns = c(opt[[1]][2,1],opt[[1]][2,2],opt[[1]][2,3])
set.seed(10)
res<-HIER_SIM(size=size, p=Prevalence, S=length(ns), SE=sensitivity, SP=specificity, ns=ns, iter=5000)

FS.H3=c(round(mean(res$'Number of tests')),round(sd(res$'Number of tests')),
        formatC(mean(res$'Number of tests')/size,digit=3,format="f"),
        formatC(mean(res$'PSE1'),digit=3,format="f"),
        formatC(mean(res$'PSE2'),digit=3,format="f"),
        formatC(mean(res$'PSP1'),digit=3,format="f"),
        formatC(mean(res$'PSP2'),digit=3,format="f"),
        formatC(mean(res$'PPV1'),digit=3,format="f"),
        formatC(mean(res$'PPV2'),digit=3,format="f"),
        formatC(mean(res$'NPV1'),digit=3,format="f"),
        formatC(mean(res$'NPV2'),digit=3,format="f"))
FS.H3
res.H3=res

set.seed(10)
res<-AT_SIM(size=size,p=Prevalence,SE=sensitivity,SP=specificity,row=opt[[3]][1],col=opt[[3]][1],iter=5000)

FS.AT=c(round(mean(res$'Number of tests')),round(sd(res$'Number of tests')),
        formatC(mean(res$'Number of tests')/size,digit=3,format="f"),
        formatC(mean(res$'PSE1'),digit=3,format="f"),
        formatC(mean(res$'PSE2'),digit=3,format="f"),
        formatC(mean(res$'PSP1'),digit=3,format="f"),
        formatC(mean(res$'PSP2'),digit=3,format="f"),
        formatC(mean(res$'PPV1'),digit=3,format="f"),
        formatC(mean(res$'PPV2'),digit=3,format="f"),
        formatC(mean(res$'NPV1'),digit=3,format="f"),
        formatC(mean(res$'NPV2'),digit=3,format="f"))
FS.AT
res.AT=res

par(mar=c(6.1,6.1,4.1,2.1))

res.box=cbind(res.H2$'Number of tests'/size,
              res.AT$'Number of tests'/size,
              res.H3$'Number of tests'/size)
boxplot(res.box[,1],res.box[,2],res.box[,3],range=0,cex.axis=1.5,col="grey",xaxes=F)#, xlab = "Algorithm", names=c("H2", "AT", "H3"),
#ylab = "", cex.lab=1.8)

mtext("Algorithm",side=1,line=4.5,cex=1.5)
mtext("Number of tests per individual",side=2,line=4.5,cex=1.5)
mtext("Female Swab",side=3,cex=1.5,line=1.5)
box()
axis(1,at=c(1,2,3),labels=c("H2","AT","H3"),tck=-0.01,cex.axis=1.4)
###################################################################
# Reproduce results related to the stratum "Female Urine."

ind=2
Prevalence = round(unlist(SHL.data[ind,3:5]/SHL.data[ind,6]),3)
Prevalence=c(1-sum(Prevalence),Prevalence)
sensitivity = unlist(SHL.data[ind,c(7,9)])
specificity = unlist(SHL.data[ind,c(8,10)])
size = unlist(SHL.data[ind,6])

opt=opt.FU

Table3$Algorithm[(ind-1)*3+1]=paste("H2(",opt[[1]][1,1],":",opt[[1]][1,2],")",sep="")
Table3$Algorithm[(ind-1)*3+2]=paste("AT(",opt[[3]][1],"*",opt[[3]][1],")",sep="")
Table3$Algorithm[(ind-1)*3+3]=paste("H3(",opt[[1]][2,1],":",opt[[1]][2,2],":",opt[[1]][2,3],")",sep="")

ns = c(opt[[1]][1,1],opt[[1]][1,2])
set.seed(10)
res<-HIER_SIM(size=size, p=Prevalence, S=length(ns), SE=sensitivity, SP=specificity, ns=ns, iter=5000)

FU.H2=c(round(mean(res$'Number of tests')),round(sd(res$'Number of tests')),
        formatC(mean(res$'Number of tests')/size,digit=3,format="f"),
        formatC(mean(res$'PSE1'),digit=3,format="f"),
        formatC(mean(res$'PSE2'),digit=3,format="f"),
        formatC(mean(res$'PSP1'),digit=3,format="f"),
        formatC(mean(res$'PSP2'),digit=3,format="f"),
        formatC(mean(res$'PPV1'),digit=3,format="f"),
        formatC(mean(res$'PPV2'),digit=3,format="f"),
        formatC(mean(res$'NPV1'),digit=3,format="f"),
        formatC(mean(res$'NPV2'),digit=3,format="f"))
FU.H2
res.H2=res

ns = c(opt[[1]][2,1],opt[[1]][2,2],opt[[1]][2,3])
set.seed(10)
res<-HIER_SIM(size=size, p=Prevalence, S=length(ns), SE=sensitivity, SP=specificity, ns=ns, iter=5000)

FU.H3=c(round(mean(res$'Number of tests')),round(sd(res$'Number of tests')),
        formatC(mean(res$'Number of tests')/size,digit=3,format="f"),
        formatC(mean(res$'PSE1'),digit=3,format="f"),
        formatC(mean(res$'PSE2'),digit=3,format="f"),
        formatC(mean(res$'PSP1'),digit=3,format="f"),
        formatC(mean(res$'PSP2'),digit=3,format="f"),
        formatC(mean(res$'PPV1'),digit=3,format="f"),
        formatC(mean(res$'PPV2'),digit=3,format="f"),
        formatC(mean(res$'NPV1'),digit=3,format="f"),
        formatC(mean(res$'NPV2'),digit=3,format="f"))
FU.H3
res.H3=res

set.seed(10)
res<-AT_SIM(size=size,p=Prevalence,SE=sensitivity,SP=specificity,row=opt[[3]][1],col=opt[[3]][1],iter=5000)

FU.AT=c(round(mean(res$'Number of tests')),round(sd(res$'Number of tests')),
        formatC(mean(res$'Number of tests')/size,digit=3,format="f"),
        formatC(mean(res$'PSE1'),digit=3,format="f"),
        formatC(mean(res$'PSE2'),digit=3,format="f"),
        formatC(mean(res$'PSP1'),digit=3,format="f"),
        formatC(mean(res$'PSP2'),digit=3,format="f"),
        formatC(mean(res$'PPV1'),digit=3,format="f"),
        formatC(mean(res$'PPV2'),digit=3,format="f"),
        formatC(mean(res$'NPV1'),digit=3,format="f"),
        formatC(mean(res$'NPV2'),digit=3,format="f"))
FU.AT
res.AT=res

par(mar=c(6.1,6.1,4.1,2.1))

res.box=cbind(res.H2$'Number of tests'/size,
              res.AT$'Number of tests'/size,
              res.H3$'Number of tests'/size)
boxplot(res.box[,1],res.box[,2],res.box[,3],range=0,cex.axis=1.5,col="grey",xaxes=F)#, xlab = "Algorithm", names=c("H2", "AT", "H3"),
#ylab = "", cex.lab=1.8)

mtext("Algorithm",side=1,line=4.5,cex=1.5)
mtext("Number of tests per individual",side=2,line=4.5,cex=1.5)
mtext("Female Urine",side=3,cex=1.5,line=1.5)
box()
axis(1,at=c(1,2,3),labels=c("H2","AT","H3"),tck=-0.01,cex.axis=1.4)
###################################################################
# Reproduce results related to the stratum "Male Swab."

ind=3
Prevalence = round(unlist(SHL.data[ind,3:5]/SHL.data[ind,6]),3)
Prevalence=c(1-sum(Prevalence),Prevalence)
sensitivity = unlist(SHL.data[ind,c(7,9)])
specificity = unlist(SHL.data[ind,c(8,10)])
size = unlist(SHL.data[ind,6])

opt=opt.MS

Table3$Algorithm[(ind-1)*3+1]=paste("H2(",opt[[1]][1,1],":",opt[[1]][1,2],")",sep="")
Table3$Algorithm[(ind-1)*3+2]=paste("AT(",opt[[3]][1],"*",opt[[3]][1],")",sep="")
Table3$Algorithm[(ind-1)*3+3]=paste("H3(",opt[[1]][2,1],":",opt[[1]][2,2],":",opt[[1]][2,3],")",sep="")

ns = c(opt[[1]][1,1],opt[[1]][1,2])
set.seed(10)
res<-HIER_SIM(size=size, p=Prevalence, S=length(ns), SE=sensitivity, SP=specificity, ns=ns, iter=5000)

MS.H2=c(round(mean(res$'Number of tests')),round(sd(res$'Number of tests')),
        formatC(mean(res$'Number of tests')/size,digit=3,format="f"),
        formatC(mean(res$'PSE1'),digit=3,format="f"),
        formatC(mean(res$'PSE2'),digit=3,format="f"),
        formatC(mean(res$'PSP1'),digit=3,format="f"),
        formatC(mean(res$'PSP2'),digit=3,format="f"),
        formatC(mean(res$'PPV1'),digit=3,format="f"),
        formatC(mean(res$'PPV2'),digit=3,format="f"),
        formatC(mean(res$'NPV1'),digit=3,format="f"),
        formatC(mean(res$'NPV2'),digit=3,format="f"))
MS.H2
res.H2=res

ns = c(opt[[1]][2,1],opt[[1]][2,2],opt[[1]][2,3])
set.seed(10)
res<-HIER_SIM(size=size, p=Prevalence, S=length(ns), SE=sensitivity, SP=specificity, ns=ns, iter=5000)

MS.H3=c(round(mean(res$'Number of tests')),round(sd(res$'Number of tests')),
        formatC(mean(res$'Number of tests')/size,digit=3,format="f"),
        formatC(mean(res$'PSE1'),digit=3,format="f"),
        formatC(mean(res$'PSE2'),digit=3,format="f"),
        formatC(mean(res$'PSP1'),digit=3,format="f"),
        formatC(mean(res$'PSP2'),digit=3,format="f"),
        formatC(mean(res$'PPV1'),digit=3,format="f"),
        formatC(mean(res$'PPV2'),digit=3,format="f"),
        formatC(mean(res$'NPV1'),digit=3,format="f"),
        formatC(mean(res$'NPV2'),digit=3,format="f"))
MS.H3
res.H3=res

set.seed(10)
res<-AT_SIM(size=size,p=Prevalence,SE=sensitivity,SP=specificity,row=opt[[3]][1],col=opt[[3]][1],iter=5000)

MS.AT=c(round(mean(res$'Number of tests')),round(sd(res$'Number of tests')),
        formatC(mean(res$'Number of tests')/size,digit=3,format="f"),
        formatC(mean(res$'PSE1'),digit=3,format="f"),
        formatC(mean(res$'PSE2'),digit=3,format="f"),
        formatC(mean(res$'PSP1'),digit=3,format="f"),
        formatC(mean(res$'PSP2'),digit=3,format="f"),
        formatC(mean(res$'PPV1'),digit=3,format="f"),
        formatC(mean(res$'PPV2'),digit=3,format="f"),
        formatC(mean(res$'NPV1'),digit=3,format="f"),
        formatC(mean(res$'NPV2'),digit=3,format="f"))
MS.AT
res.AT=res

par(mar=c(6.1,6.1,4.1,2.1))

res.box=cbind(res.H2$'Number of tests'/size,
              res.AT$'Number of tests'/size,
              res.H3$'Number of tests'/size)
boxplot(res.box[,1],res.box[,2],res.box[,3],range=0,cex.axis=1.5,col="grey",xaxes=F)#, xlab = "Algorithm", names=c("H2", "AT", "H3"),
#ylab = "", cex.lab=1.8)

mtext("Algorithm",side=1,line=4.5,cex=1.5)
mtext("Number of tests per individual",side=2,line=4.5,cex=1.5)
mtext("Male Swab",side=3,cex=1.5,line=1.5)
box()
axis(1,at=c(1,2,3),labels=c("H2","AT","H3"),tck=-0.01,cex.axis=1.4)
###################################################################
# Reproduce results related to the stratum "Male Urine."

ind=4
Prevalence = round(unlist(SHL.data[ind,3:5]/SHL.data[ind,6]),3)
Prevalence=c(1-sum(Prevalence),Prevalence)
sensitivity = unlist(SHL.data[ind,c(7,9)])
specificity = unlist(SHL.data[ind,c(8,10)])
size = unlist(SHL.data[ind,6])

opt=opt.MU

Table3$Algorithm[(ind-1)*3+1]=paste("H2(",opt[[1]][1,1],":",opt[[1]][1,2],")",sep="")
Table3$Algorithm[(ind-1)*3+2]=paste("AT(",opt[[3]][1],"*",opt[[3]][1],")",sep="")
Table3$Algorithm[(ind-1)*3+3]=paste("H3(",opt[[1]][2,1],":",opt[[1]][2,2],":",opt[[1]][2,3],")",sep="")

ns = c(opt[[1]][1,1],opt[[1]][1,2])
set.seed(10)
res<-HIER_SIM(size=size, p=Prevalence, S=length(ns), SE=sensitivity, SP=specificity, ns=ns, iter=5000)

MU.H2=c(round(mean(res$'Number of tests')),round(sd(res$'Number of tests')),
        formatC(mean(res$'Number of tests')/size,digit=3,format="f"),
        formatC(mean(res$'PSE1'),digit=3,format="f"),
        formatC(mean(res$'PSE2'),digit=3,format="f"),
        formatC(mean(res$'PSP1'),digit=3,format="f"),
        formatC(mean(res$'PSP2'),digit=3,format="f"),
        formatC(mean(res$'PPV1'),digit=3,format="f"),
        formatC(mean(res$'PPV2'),digit=3,format="f"),
        formatC(mean(res$'NPV1'),digit=3,format="f"),
        formatC(mean(res$'NPV2'),digit=3,format="f"))
MU.H2
res.H2=res

ns = c(opt[[1]][2,1],opt[[1]][2,2],opt[[1]][2,3])
set.seed(10)
res<-HIER_SIM(size=size, p=Prevalence, S=length(ns), SE=sensitivity, SP=specificity, ns=ns, iter=5000)

MU.H3=c(round(mean(res$'Number of tests')),round(sd(res$'Number of tests')),
        formatC(mean(res$'Number of tests')/size,digit=3,format="f"),
        formatC(mean(res$'PSE1'),digit=3,format="f"),
        formatC(mean(res$'PSE2'),digit=3,format="f"),
        formatC(mean(res$'PSP1'),digit=3,format="f"),
        formatC(mean(res$'PSP2'),digit=3,format="f"),
        formatC(mean(res$'PPV1'),digit=3,format="f"),
        formatC(mean(res$'PPV2'),digit=3,format="f"),
        formatC(mean(res$'NPV1'),digit=3,format="f"),
        formatC(mean(res$'NPV2'),digit=3,format="f"))
MU.H3
res.H3=res

set.seed(10)
res<-AT_SIM(size=size,p=Prevalence,SE=sensitivity,SP=specificity,row=opt[[3]][1],col=opt[[3]][1],iter=5000)

MU.AT=c(round(mean(res$'Number of tests')),round(sd(res$'Number of tests')),
        formatC(mean(res$'Number of tests')/size,digit=3,format="f"),
        formatC(mean(res$'PSE1'),digit=3,format="f"),
        formatC(mean(res$'PSE2'),digit=3,format="f"),
        formatC(mean(res$'PSP1'),digit=3,format="f"),
        formatC(mean(res$'PSP2'),digit=3,format="f"),
        formatC(mean(res$'PPV1'),digit=3,format="f"),
        formatC(mean(res$'PPV2'),digit=3,format="f"),
        formatC(mean(res$'NPV1'),digit=3,format="f"),
        formatC(mean(res$'NPV2'),digit=3,format="f"))
MU.AT
res.AT=res

par(mar=c(6.1,6.1,4.1,2.1))

res.box=cbind(res.H2$'Number of tests'/size,
              res.AT$'Number of tests'/size,
              res.H3$'Number of tests'/size)
boxplot(res.box[,1],res.box[,2],res.box[,3],range=0,cex.axis=1.5,col="grey",xaxes=F)#, xlab = "Algorithm", names=c("H2", "AT", "H3"),
        #ylab = "", cex.lab=1.8)

mtext("Algorithm",side=1,line=4.5,cex=1.5)
mtext("Number of tests per individual",side=2,line=4.5,cex=1.5)
mtext("Male Urine",side=3,cex=1.5,line=1.5)
box()
axis(1,at=c(1,2,3),labels=c("H2","AT","H3"),tck=-0.01,cex.axis=1.4)


##################################################
# Save results to Table 3.
Table3=cbind(Table3,rbind(FS.H2,FS.AT,FS.H3,FU.H2,FU.AT,FU.H3,MS.H2,MS.AT,MS.H3,MU.H2,MU.AT,MU.H3))
colnames(Table3)=c("Stratum","Total","Algorithm","Mean","SD","EFF","PSE1","PSE2","PSP1","PSP2","PPV1","PPV2","NPV1","NPV2")
##################################################
# Reproducing Table 3 in the paper is finished.
Table3
##################################################

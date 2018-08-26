############################################################################################
############################################################################################
# This R program will reproduce the results of the real data analysis in 
#                 "Array testing with Multiplex Assays."                                                                          
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

############################################################################################
# Read data
############################################################################################
library(utils)
SHL.data=read.csv("IowaSHLdata.csv",col.names =c("Stratum","CT_N,NG_N","CT_P,NG_N","CT_N,NG_P","CT_P,NG_P","Total",
                                                 "Se1","Sp1","Se2","Sp2"))
############################################################################################
# Find the optimal design of each of H2, H3, AT, and ATM algorithm for the stratum "Female Swab."
ind=1
Prevalence = round(unlist(SHL.data[ind,3:5]/SHL.data[ind,6]),3)
Prevalence=c(1-sum(Prevalence),Prevalence)
sensitivity = unlist(SHL.data[ind,c(7,9)])
specificity = unlist(SHL.data[ind,c(8,10)])
maximal=100
n_at=20
n_atm=10
OPT.H_VS_A(p=Prevalence,SE=sensitivity,SP=specificity,MAMPS=maximal, n_at=n_at, n_atm=n_atm)

############################################################################################
# Find the optimal design of each of H2, H3, AT, and ATM algorithm for the stratum "Female Urine."
ind=2
Prevalence = round(unlist(SHL.data[ind,3:5]/SHL.data[ind,6]),3)
Prevalence=c(1-sum(Prevalence),Prevalence)
sensitivity = unlist(SHL.data[ind,c(7,9)])
specificity = unlist(SHL.data[ind,c(8,10)])
maximal=100
n_at=20
n_atm=10
OPT.H_VS_A(p=Prevalence,SE=sensitivity,SP=specificity,MAMPS=maximal, n_at=n_at, n_atm=n_atm)

############################################################################################
# Find the optimal design of each of H2, H3, AT, and ATM algorithm for the stratum "Male Swab."
ind=3
Prevalence = round(unlist(SHL.data[ind,3:5]/SHL.data[ind,6]),3)
Prevalence=c(1-sum(Prevalence),Prevalence)
sensitivity = unlist(SHL.data[ind,c(7,9)])
specificity = unlist(SHL.data[ind,c(8,10)])
maximal=100
n_at=20
n_atm=10
OPT.H_VS_A(p=Prevalence,SE=sensitivity,SP=specificity,MAMPS=maximal, n_at=n_at, n_atm=n_atm)

############################################################################################
# Find the optimal design of each of H2, H3, AT, and ATM algorithm for the stratum "Male Urine."
ind=4
Prevalence = round(unlist(SHL.data[ind,3:5]/SHL.data[ind,6]),3)
Prevalence=c(1-sum(Prevalence),Prevalence)
sensitivity = unlist(SHL.data[ind,c(7,9)])
specificity = unlist(SHL.data[ind,c(8,10)])
maximal=100
n_at=20
n_atm=10
OPT.H_VS_A(p=Prevalence,SE=sensitivity,SP=specificity,MAMPS=maximal, n_at=n_at, n_atm=n_atm)

#install.packages('gplots')
library('mvmeta')
library(xtable)
library('gplots')


set.seed(1234)

alp<-0.05
################### Determining X:p\times n and U:pn\times pn
source("Data_CCAUV_V_K1.R")
PTB_BB_n<-ncol(PTB_BB)
PTB_BB_mean<-PTB_BB%*%rep(1,PTB_BB_n)/PTB_BB_n
PTB_BB_S<-PTB_BB%*%(diag(PTB_BB_n)-matrix(1,PTB_BB_n,PTB_BB_n)/PTB_BB_n)%*%t(PTB_BB)/(PTB_BB_n-1)

PTB_SE_n<-ncol(PTB_SE)
PTB_SE_mean<-PTB_SE%*%rep(1,PTB_SE_n)/PTB_SE_n
PTB_SE_S<-PTB_SE%*%(diag(PTB_SE_n)-matrix(1,PTB_SE_n,PTB_SE_n)/PTB_SE_n)%*%t(PTB_SE)/(PTB_SE_n-1)

###
BNM_CESTA_BB_n<-ncol(BNM_CESTA_BB)
BNM_CESTA_BB_mean<-BNM_CESTA_BB%*%rep(1,BNM_CESTA_BB_n)/BNM_CESTA_BB_n
BNM_CESTA_BB_S<-BNM_CESTA_BB%*%(diag(BNM_CESTA_BB_n)-matrix(1,BNM_CESTA_BB_n,BNM_CESTA_BB_n)/BNM_CESTA_BB_n)%*%t(BNM_CESTA_BB)/(BNM_CESTA_BB_n-1)

BNM_CESTA_SE_n<-ncol(BNM_CESTA_SE)
BNM_CESTA_SE_mean<-BNM_CESTA_SE%*%rep(1,BNM_CESTA_SE_n)/BNM_CESTA_SE_n
BNM_CESTA_SE_S<-BNM_CESTA_SE%*%(diag(BNM_CESTA_SE_n)-matrix(1,BNM_CESTA_SE_n,BNM_CESTA_SE_n)/BNM_CESTA_SE_n)%*%t(BNM_CESTA_SE)/(BNM_CESTA_SE_n-1)

###
CSIRO_NML_BB_n<-ncol(CSIRO_NML_BB)
CSIRO_NML_BB_mean<-CSIRO_NML_BB%*%rep(1,CSIRO_NML_BB_n)/CSIRO_NML_BB_n
CSIRO_NML_BB_S<-CSIRO_NML_BB%*%(diag(CSIRO_NML_BB_n)-matrix(1,CSIRO_NML_BB_n,CSIRO_NML_BB_n)/CSIRO_NML_BB_n)%*%t(CSIRO_NML_BB)/(CSIRO_NML_BB_n-1)

CSIRO_NML_SE_n<-ncol(CSIRO_NML_SE)
CSIRO_NML_SE_mean<-CSIRO_NML_SE%*%rep(1,CSIRO_NML_SE_n)/CSIRO_NML_SE_n
CSIRO_NML_SE_S<-CSIRO_NML_SE%*%(diag(CSIRO_NML_SE_n)-matrix(1,CSIRO_NML_SE_n,CSIRO_NML_SE_n)/CSIRO_NML_SE_n)%*%t(CSIRO_NML_SE)/(CSIRO_NML_SE_n-1)

###
CMI_BB_n<-ncol(CMI_BB)
CMI_BB_mean<-CMI_BB%*%rep(1,CMI_BB_n)/CMI_BB_n
CMI_BB_S<-CMI_BB%*%(diag(CMI_BB_n)-matrix(1,CMI_BB_n,CMI_BB_n)/CMI_BB_n)%*%t(CMI_BB)/(CMI_BB_n-1)

CMI_SE_n<-ncol(CMI_SE)
CMI_SE_mean<-CMI_SE%*%rep(1,CMI_SE_n)/CMI_SE_n
CMI_SE_S<-CMI_SE%*%(diag(CMI_SE_n)-matrix(1,CMI_SE_n,CMI_SE_n)/CMI_SE_n)%*%t(CMI_SE)/(CMI_SE_n-1)

###
CSIR_NML_BB_n<-ncol(CSIR_NML_BB)
CSIR_NML_BB_mean<-CSIR_NML_BB%*%rep(1,CSIR_NML_BB_n)/CSIR_NML_BB_n
CSIR_NML_BB_S<-CSIR_NML_BB%*%(diag(CSIR_NML_BB_n)-matrix(1,CSIR_NML_BB_n,CSIR_NML_BB_n)/CSIR_NML_BB_n)%*%t(CSIR_NML_BB)/(CSIR_NML_BB_n-1)

CSIR_NML_SE_n<-ncol(CSIR_NML_SE)
CSIR_NML_SE_mean<-CSIR_NML_SE%*%rep(1,CSIR_NML_SE_n)/CSIR_NML_SE_n
CSIR_NML_SE_S<-CSIR_NML_SE%*%(diag(CSIR_NML_SE_n)-matrix(1,CSIR_NML_SE_n,CSIR_NML_SE_n)/CSIR_NML_SE_n)%*%t(CSIR_NML_SE)/(CSIR_NML_SE_n-1)

###
CENAM_BB_n<-ncol(CENAM_BB)
CENAM_BB_mean<-CENAM_BB%*%rep(1,CENAM_BB_n)/CENAM_BB_n
CENAM_BB_S<-CENAM_BB%*%(diag(CENAM_BB_n)-matrix(1,CENAM_BB_n,CENAM_BB_n)/CENAM_BB_n)%*%t(CENAM_BB)/(CENAM_BB_n-1)

CENAM_SE_n<-ncol(CENAM_SE)
CENAM_SE_mean<-CENAM_SE%*%rep(1,CENAM_SE_n)/CENAM_SE_n
CENAM_SE_S<-CENAM_SE%*%(diag(CENAM_SE_n)-matrix(1,CENAM_SE_n,CENAM_SE_n)/CENAM_SE_n)%*%t(CENAM_SE)/(CENAM_SE_n-1)

###
NRC_BB_n<-ncol(NRC_BB)
NRC_BB_mean<-NRC_BB%*%rep(1,NRC_BB_n)/NRC_BB_n
NRC_BB_S<-NRC_BB%*%(diag(NRC_BB_n)-matrix(1,NRC_BB_n,NRC_BB_n)/NRC_BB_n)%*%t(NRC_BB)/(NRC_BB_n-1)

NRC_SE_n<-ncol(NRC_SE)
NRC_SE_mean<-NRC_SE%*%rep(1,NRC_SE_n)/NRC_SE_n
NRC_SE_S<-NRC_SE%*%(diag(NRC_SE_n)-matrix(1,NRC_SE_n,NRC_SE_n)/NRC_SE_n)%*%t(NRC_SE)/(NRC_SE_n-1)

###
KRISS_BB_n<-ncol(KRISS_BB)
KRISS_BB_mean<-KRISS_BB%*%rep(1,KRISS_BB_n)/KRISS_BB_n
KRISS_BB_S<-KRISS_BB%*%(diag(KRISS_BB_n)-matrix(1,KRISS_BB_n,KRISS_BB_n)/KRISS_BB_n)%*%t(KRISS_BB)/(KRISS_BB_n-1)

KRISS_SE_n<-ncol(KRISS_SE)
KRISS_SE_mean<-KRISS_SE%*%rep(1,KRISS_SE_n)/KRISS_SE_n
KRISS_SE_S<-KRISS_SE%*%(diag(KRISS_SE_n)-matrix(1,KRISS_SE_n,KRISS_SE_n)/KRISS_SE_n)%*%t(KRISS_SE)/(KRISS_SE_n-1)

###
NMIJ_BB_n<-ncol(NMIJ_BB)
NMIJ_BB_mean<-NMIJ_BB%*%rep(1,NMIJ_BB_n)/NMIJ_BB_n
NMIJ_BB_S<-NMIJ_BB%*%(diag(NMIJ_BB_n)-matrix(1,NMIJ_BB_n,NMIJ_BB_n)/NMIJ_BB_n)%*%t(NMIJ_BB)/(NMIJ_BB_n-1)

NMIJ_SE_n<-ncol(NMIJ_SE)
NMIJ_SE_mean<-NMIJ_SE%*%rep(1,NMIJ_SE_n)/NMIJ_SE_n
NMIJ_SE_S<-NMIJ_SE%*%(diag(NMIJ_SE_n)-matrix(1,NMIJ_SE_n,NMIJ_SE_n)/NMIJ_SE_n)%*%t(NMIJ_SE)/(NMIJ_SE_n-1)

###
VNIIM_BB_n<-ncol(VNIIM_BB)
VNIIM_BB_mean<-VNIIM_BB%*%rep(1,VNIIM_BB_n)/VNIIM_BB_n
VNIIM_BB_S<-VNIIM_BB%*%(diag(VNIIM_BB_n)-matrix(1,VNIIM_BB_n,VNIIM_BB_n)/VNIIM_BB_n)%*%t(VNIIM_BB)/(VNIIM_BB_n-1)

VNIIM_SE_n<-ncol(VNIIM_SE)
VNIIM_SE_mean<-VNIIM_SE%*%rep(1,VNIIM_SE_n)/VNIIM_SE_n
VNIIM_SE_S<-VNIIM_SE%*%(diag(VNIIM_SE_n)-matrix(1,VNIIM_SE_n,VNIIM_SE_n)/VNIIM_SE_n)%*%t(VNIIM_SE)/(VNIIM_SE_n-1)

###
NIST_BB_n<-ncol(NIST_BB)
NIST_BB_mean<-NIST_BB%*%rep(1,NIST_BB_n)/NIST_BB_n
NIST_BB_S<-NIST_BB%*%(diag(NIST_BB_n)-matrix(1,NIST_BB_n,NIST_BB_n)/NIST_BB_n)%*%t(NIST_BB)/(NIST_BB_n-1)

NIST_SE_n<-ncol(NIST_SE)
NIST_SE_mean<-NIST_SE%*%rep(1,NIST_SE_n)/NIST_SE_n
NIST_SE_S<-NIST_SE%*%(diag(NIST_SE_n)-matrix(1,NIST_SE_n,NIST_SE_n)/NIST_SE_n)%*%t(NIST_SE)/(NIST_SE_n-1)

###
NMI_VSL_BB_n<-ncol(NMI_VSL_BB)
NMI_VSL_BB_mean<-NMI_VSL_BB%*%rep(1,NMI_VSL_BB_n)/NMI_VSL_BB_n
NMI_VSL_BB_S<-NMI_VSL_BB%*%(diag(NMI_VSL_BB_n)-matrix(1,NMI_VSL_BB_n,NMI_VSL_BB_n)/NMI_VSL_BB_n)%*%t(NMI_VSL_BB)/(NMI_VSL_BB_n-1)

NMI_VSL_SE_n<-ncol(NMI_VSL_SE)
NMI_VSL_SE_mean<-NMI_VSL_SE%*%rep(1,NMI_VSL_SE_n)/NMI_VSL_SE_n
NMI_VSL_SE_S<-NMI_VSL_SE%*%(diag(NMI_VSL_SE_n)-matrix(1,NMI_VSL_SE_n,NMI_VSL_SE_n)/NMI_VSL_SE_n)%*%t(NMI_VSL_SE)/(NMI_VSL_SE_n-1)


X<-t(cbind(PTB_BB_mean[c(12,13)], BNM_CESTA_BB_mean[c(12,13)], CSIRO_NML_BB_mean[c(12,13)], CMI_BB_mean[c(12,13)], CSIR_NML_BB_mean[c(12,13)], CENAM_BB_mean[c(12,13)], NRC_BB_mean[c(12,13)], KRISS_BB_mean[c(12,13)], NMIJ_BB_mean[c(12,13)], VNIIM_BB_mean[c(12,13)], NIST_BB_mean[c(12,13)], NMI_VSL_BB_mean[c(12,13)]))
S<-t(cbind(as.vector(PTB_BB_S[c(12,13),c(12,13)]), as.vector(BNM_CESTA_BB_S[c(12,13),c(12,13)]), as.vector(CSIRO_NML_BB_S[c(12,13),c(12,13)]), as.vector(CMI_BB_S[c(12,13),c(12,13)]), as.vector(CSIR_NML_BB_S[c(12,13),c(12,13)]), as.vector(CENAM_BB_S[c(12,13),c(12,13)]), as.vector(NRC_BB_S[c(12,13),c(12,13)]), as.vector(KRISS_BB_S[c(12,13),c(12,13)]), as.vector(NMIJ_BB_S[c(12,13),c(12,13)]), as.vector(VNIIM_BB_S[c(12,13),c(12,13)]), as.vector(NIST_BB_S[c(12,13),c(12,13)]), as.vector(NMI_VSL_BB_S[c(12,13),c(12,13)])))
S<-S[,-3]

res_ML<-summary(mvmeta(X,S=S,method="ml"))
res_REML<-summary(mvmeta(X,S=S,method="reml"))
res_MM<-summary(mvmeta(X,S=S,method="mm"))
print(res_ML, digits=12)
print(res_REML, digits=12)
print(res_MM, digits=12)

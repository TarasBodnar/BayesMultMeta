library(xtable)

set.seed(1234)
source("MH_sample_post.R")

alp<-0.05
Np<-10^5 # number of observations drawn from the proposal distribution
B<-10^4 # burn in 
d<-7 ## degrees of freedom in t multivariate random effects model
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

###
X_BB=cbind(PTB_BB_mean, BNM_CESTA_BB_mean, CSIRO_NML_BB_mean, CMI_BB_mean, CSIR_NML_BB_mean, CENAM_BB_mean, NRC_BB_mean, KRISS_BB_mean, NMIJ_BB_mean, VNIIM_BB_mean, NIST_BB_mean, NMI_VSL_BB_mean)
X_SE=cbind(PTB_SE_mean, BNM_CESTA_SE_mean, CSIRO_NML_SE_mean, CMI_SE_mean, CSIR_NML_SE_mean, CENAM_SE_mean, NRC_SE_mean, KRISS_SE_mean, NMIJ_SE_mean, VNIIM_SE_mean, NIST_SE_mean, NMI_VSL_SE_mean)
###
S_BB=sqrt(cbind(diag(PTB_BB_S), diag(BNM_CESTA_BB_S), diag(CSIRO_NML_BB_S), diag(CMI_BB_S), diag(CSIR_NML_BB_S), diag(CENAM_BB_S), diag(NRC_BB_S), diag(KRISS_BB_S), diag(NMIJ_BB_S), diag(VNIIM_BB_S), diag(NIST_BB_S), diag(NMI_VSL_BB_S)))
S_SE=sqrt(cbind(diag(PTB_SE_S), diag(BNM_CESTA_SE_S), diag(CSIRO_NML_SE_S), diag(CMI_SE_S), diag(CSIR_NML_SE_S), diag(CENAM_SE_S), diag(NRC_SE_S), diag(KRISS_SE_S), diag(NMIJ_SE_S), diag(VNIIM_SE_S), diag(NIST_SE_S), diag(NMI_VSL_SE_S)))

#dataREM<-hyp
################### Matrix X
#X<-t(cbind(dataREM$sbp,dataREM$dbp))
index<-c(11,12,13,14)
X<-cbind(PTB_BB_mean[index], BNM_CESTA_BB_mean[index], CSIRO_NML_BB_mean[index], CMI_BB_mean[index], CSIR_NML_BB_mean[index], CENAM_BB_mean[index], NRC_BB_mean[index], KRISS_BB_mean[index], NMIJ_BB_mean[index], NIST_BB_mean[index], NMI_VSL_BB_mean[index])
p<-nrow(X)  # model dimension
n<-ncol(X)  # sample size

################### Matrix U
U<-matrix(0,n*p,n*p)
U[1:4,1:4]<-PTB_BB_S[index,index]/PTB_BB_n 
U[5:8,5:8]<-BNM_CESTA_BB_S[index,index]/BNM_CESTA_BB_n
U[9:12,9:12]<-CSIRO_NML_BB_S[index,index]/CSIRO_NML_BB_n
U[13:16,13:16]<-CMI_BB_S[index,index]/CMI_BB_n
U[17:20,17:20]<-CSIR_NML_BB_S[index,index]/CSIR_NML_BB_n
U[21:24,21:24]<-CENAM_BB_S[index,index]/CENAM_BB_n
U[25:28,25:28]<-NRC_BB_S[index,index]/NRC_BB_n
U[29:32,29:32]<-KRISS_BB_S[index,index]/KRISS_BB_n
U[33:36,33:36]<-NMIJ_BB_S[index,index]/NMIJ_BB_n
U[37:40,37:40]<-NIST_BB_S[index,index]/NIST_BB_n
U[41:44,41:44]<-NMI_VSL_BB_S[index,index]/NMI_VSL_BB_n
#U[45:48,45:48]<-VNIIM_BB_S[index,index]/VNIIM_BB_n


############## addtional definitons  
bi_n<-rep(1,n)
tbi_n<-t(bi_n)
In<-diag(bi_n)
Jn<-matrix(1,n,n)
Ip<-diag(rep(1,p))

####### duplication matrix
Gp<-Dp(p)
Lp<-Gp%*%solve(t(Gp)%*%Gp)

################################### ranks plots #################################################################################
M<-4
Np<-10^5 # number of observations drawn from the proposal distribution
B<-10^4 # burn in 

chains_mu<-NULL
chains_Psi<-NULL
i_m<-0
while (i_m < M)
{chain<-sample_post_nor_jef_marg_mu(X,U,Np+B)
print(sum(chain[[3]]))
if (sum(chain[[3]])>3000)
{chains_mu<-cbind(chains_mu,chain[[1]][,(B+1):(Np+B)])
chains_Psi<-cbind(chains_Psi,chain[[2]][,(B+1):(Np+B)])
i_m<-i_m+1
}
}

ranks_mu1<-rank(chains_mu[1,], ties.method = "average")
ranks_mu2<-rank(chains_mu[2,], ties.method = "average")
ranks_mu3<-rank(chains_mu[3,], ties.method = "average")
ranks_mu4<-rank(chains_mu[4,], ties.method = "average")
ranks_Psi11<-rank(chains_Psi[1,], ties.method = "average")
ranks_Psi21<-rank(chains_Psi[2,], ties.method = "average")
ranks_Psi31<-rank(chains_Psi[3,], ties.method = "average")
ranks_Psi41<-rank(chains_Psi[4,], ties.method = "average")
ranks_Psi22<-rank(chains_Psi[6,], ties.method = "average")
ranks_Psi32<-rank(chains_Psi[7,], ties.method = "average")
ranks_Psi42<-rank(chains_Psi[8,], ties.method = "average")
ranks_Psi33<-rank(chains_Psi[11,], ties.method = "average")
ranks_Psi43<-rank(chains_Psi[12,], ties.method = "average")
ranks_Psi44<-rank(chains_Psi[16,], ties.method = "average")

ranks_nor_jef_mu_mu1<-matrix(ranks_mu1,Np,M)
ranks_nor_jef_mu_mu2<-matrix(ranks_mu2,Np,M)
ranks_nor_jef_mu_mu3<-matrix(ranks_mu3,Np,M)
ranks_nor_jef_mu_mu4<-matrix(ranks_mu4,Np,M)
ranks_nor_jef_mu_Psi11<-matrix(ranks_Psi11,Np,M)
ranks_nor_jef_mu_Psi21<-matrix(ranks_Psi21,Np,M)
ranks_nor_jef_mu_Psi31<-matrix(ranks_Psi31,Np,M)
ranks_nor_jef_mu_Psi41<-matrix(ranks_Psi41,Np,M)
ranks_nor_jef_mu_Psi22<-matrix(ranks_Psi22,Np,M)
ranks_nor_jef_mu_Psi32<-matrix(ranks_Psi32,Np,M)
ranks_nor_jef_mu_Psi42<-matrix(ranks_Psi42,Np,M)
ranks_nor_jef_mu_Psi33<-matrix(ranks_Psi33,Np,M)
ranks_nor_jef_mu_Psi43<-matrix(ranks_Psi43,Np,M)
ranks_nor_jef_mu_Psi44<-matrix(ranks_Psi44,Np,M)

pdf(file="ranks_nor_jef_mu_mu1_1.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_mu1[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_mu1_2.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_mu1[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_mu1_3.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_mu1[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_mu1_4.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_mu1[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_mu2_1.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_mu2[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_mu2_2.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_mu2[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_mu2_3.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_mu2[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_mu2_4.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_mu2[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()


pdf(file="ranks_nor_jef_mu_mu3_1.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_mu3[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_mu3_2.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_mu3[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_mu3_3.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_mu3[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_mu3_4.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_mu3[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_mu4_1.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_mu4[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_mu4_2.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_mu4[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_mu4_3.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_mu4[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_mu4_4.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_mu4[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi11_1.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi11[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi11_2.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi11[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi11_3.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi11[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi11_4.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi11[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi21_1.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi21[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi21_2.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi21[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi21_3.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi21[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi21_4.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi21[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi31_1.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi31[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi31_2.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi31[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi31_3.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi31[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[311]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi31_4.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi31[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi41_1.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi41[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi41_2.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi41[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi41_3.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi41[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi41_4.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi41[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi22_1.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi22[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi22_2.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi22[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi22_3.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi22[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi22_4.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi22[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi32_1.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi32[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi32_2.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi32[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi32_3.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi32[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi32_4.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi32[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi42_1.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi42[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi42_2.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi42[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi42_3.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi42[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi42_4.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi42[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi33_1.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi33[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi33_2.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi33[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi33_3.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi33[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi33_4.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi33[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi43_1.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi43[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi43_2.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi43[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi43_3.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi43[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi43_4.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi43[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi44_1.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi44[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi44_2.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi44[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi44_3.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi44[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_mu_Psi44_4.pdf", width=9, height=7)
hist(ranks_nor_jef_mu_Psi44[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()
###########

chains_mu<-NULL
chains_Psi<-NULL
i_m<-0
while (i_m < M)
{chain<-sample_post_nor_jef_marg_Psi(X,U,Np+B)
print(sum(chain[[3]]))
if (sum(chain[[3]])>3000)
{chains_mu<-cbind(chains_mu,chain[[1]][,(B+1):(Np+B)])
chains_Psi<-cbind(chains_Psi,chain[[2]][,(B+1):(Np+B)])
i_m<-i_m+1
}
}

ranks_mu1<-rank(chains_mu[1,], ties.method = "average")
ranks_mu2<-rank(chains_mu[2,], ties.method = "average")
ranks_mu3<-rank(chains_mu[3,], ties.method = "average")
ranks_mu4<-rank(chains_mu[4,], ties.method = "average")
ranks_Psi11<-rank(chains_Psi[1,], ties.method = "average")
ranks_Psi21<-rank(chains_Psi[2,], ties.method = "average")
ranks_Psi31<-rank(chains_Psi[3,], ties.method = "average")
ranks_Psi41<-rank(chains_Psi[4,], ties.method = "average")
ranks_Psi22<-rank(chains_Psi[6,], ties.method = "average")
ranks_Psi32<-rank(chains_Psi[7,], ties.method = "average")
ranks_Psi42<-rank(chains_Psi[8,], ties.method = "average")
ranks_Psi33<-rank(chains_Psi[11,], ties.method = "average")
ranks_Psi43<-rank(chains_Psi[12,], ties.method = "average")
ranks_Psi44<-rank(chains_Psi[16,], ties.method = "average")

ranks_nor_jef_Psi_mu1<-matrix(ranks_mu1,Np,M)
ranks_nor_jef_Psi_mu2<-matrix(ranks_mu2,Np,M)
ranks_nor_jef_Psi_mu3<-matrix(ranks_mu3,Np,M)
ranks_nor_jef_Psi_mu4<-matrix(ranks_mu4,Np,M)
ranks_nor_jef_Psi_Psi11<-matrix(ranks_Psi11,Np,M)
ranks_nor_jef_Psi_Psi21<-matrix(ranks_Psi21,Np,M)
ranks_nor_jef_Psi_Psi31<-matrix(ranks_Psi31,Np,M)
ranks_nor_jef_Psi_Psi41<-matrix(ranks_Psi41,Np,M)
ranks_nor_jef_Psi_Psi22<-matrix(ranks_Psi22,Np,M)
ranks_nor_jef_Psi_Psi32<-matrix(ranks_Psi32,Np,M)
ranks_nor_jef_Psi_Psi42<-matrix(ranks_Psi42,Np,M)
ranks_nor_jef_Psi_Psi33<-matrix(ranks_Psi33,Np,M)
ranks_nor_jef_Psi_Psi43<-matrix(ranks_Psi43,Np,M)
ranks_nor_jef_Psi_Psi44<-matrix(ranks_Psi44,Np,M)

pdf(file="ranks_nor_jef_Psi_mu1_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_mu1[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_mu1_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_mu1[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_mu1_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_mu1[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_mu1_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_mu1[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_mu2_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_mu2[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_mu2_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_mu2[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_mu2_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_mu2[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_mu2_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_mu2[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()


pdf(file="ranks_nor_jef_Psi_mu3_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_mu3[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_mu3_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_mu3[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_mu3_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_mu3[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_mu3_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_mu3[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_mu4_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_mu4[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_mu4_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_mu4[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_mu4_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_mu4[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_mu4_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_mu4[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi11_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi11[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi11_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi11[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi11_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi11[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi11_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi11[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi21_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi21[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi21_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi21[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi21_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi21[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi21_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi21[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi31_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi31[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi31_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi31[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi31_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi31[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[311]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi31_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi31[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi41_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi41[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi41_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi41[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi41_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi41[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi41_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi41[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi22_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi22[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi22_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi22[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi22_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi22[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi22_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi22[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi32_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi32[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi32_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi32[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi32_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi32[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi32_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi32[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi42_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi42[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi42_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi42[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi42_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi42[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi42_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi42[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi33_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi33[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi33_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi33[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi33_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi33[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi33_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi33[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi43_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi43[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi43_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi43[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi43_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi43[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi43_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi43[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi44_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi44[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi44_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi44[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi44_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi44[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Psi_Psi44_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Psi_Psi44[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()
###########################################################################

###########################################################################

chains_mu<-NULL
chains_Psi<-NULL
i_m<-0
while (i_m < M)
{chain<-sample_post_nor_ref_marg_mu(X,U,Np+B)
print(sum(chain[[3]]))
if (sum(chain[[3]])>3000)
{chains_mu<-cbind(chains_mu,chain[[1]][,(B+1):(Np+B)])
chains_Psi<-cbind(chains_Psi,chain[[2]][,(B+1):(Np+B)])
i_m<-i_m+1
}
}

ranks_mu1<-rank(chains_mu[1,], ties.method = "average")
ranks_mu2<-rank(chains_mu[2,], ties.method = "average")
ranks_mu3<-rank(chains_mu[3,], ties.method = "average")
ranks_mu4<-rank(chains_mu[4,], ties.method = "average")
ranks_Psi11<-rank(chains_Psi[1,], ties.method = "average")
ranks_Psi21<-rank(chains_Psi[2,], ties.method = "average")
ranks_Psi31<-rank(chains_Psi[3,], ties.method = "average")
ranks_Psi41<-rank(chains_Psi[4,], ties.method = "average")
ranks_Psi22<-rank(chains_Psi[6,], ties.method = "average")
ranks_Psi32<-rank(chains_Psi[7,], ties.method = "average")
ranks_Psi42<-rank(chains_Psi[8,], ties.method = "average")
ranks_Psi33<-rank(chains_Psi[11,], ties.method = "average")
ranks_Psi43<-rank(chains_Psi[12,], ties.method = "average")
ranks_Psi44<-rank(chains_Psi[16,], ties.method = "average")

ranks_nor_ref_mu_mu1<-matrix(ranks_mu1,Np,M)
ranks_nor_ref_mu_mu2<-matrix(ranks_mu2,Np,M)
ranks_nor_ref_mu_mu3<-matrix(ranks_mu3,Np,M)
ranks_nor_ref_mu_mu4<-matrix(ranks_mu4,Np,M)
ranks_nor_ref_mu_Psi11<-matrix(ranks_Psi11,Np,M)
ranks_nor_ref_mu_Psi21<-matrix(ranks_Psi21,Np,M)
ranks_nor_ref_mu_Psi31<-matrix(ranks_Psi31,Np,M)
ranks_nor_ref_mu_Psi41<-matrix(ranks_Psi41,Np,M)
ranks_nor_ref_mu_Psi22<-matrix(ranks_Psi22,Np,M)
ranks_nor_ref_mu_Psi32<-matrix(ranks_Psi32,Np,M)
ranks_nor_ref_mu_Psi42<-matrix(ranks_Psi42,Np,M)
ranks_nor_ref_mu_Psi33<-matrix(ranks_Psi33,Np,M)
ranks_nor_ref_mu_Psi43<-matrix(ranks_Psi43,Np,M)
ranks_nor_ref_mu_Psi44<-matrix(ranks_Psi44,Np,M)

pdf(file="ranks_nor_ref_mu_mu1_1.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_mu1[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_mu1_2.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_mu1[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_mu1_3.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_mu1[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_mu1_4.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_mu1[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_mu2_1.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_mu2[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_mu2_2.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_mu2[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_mu2_3.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_mu2[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_mu2_4.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_mu2[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()


pdf(file="ranks_nor_ref_mu_mu3_1.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_mu3[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_mu3_2.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_mu3[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_mu3_3.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_mu3[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_mu3_4.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_mu3[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_mu4_1.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_mu4[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_mu4_2.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_mu4[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_mu4_3.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_mu4[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_mu4_4.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_mu4[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi11_1.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi11[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi11_2.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi11[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi11_3.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi11[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi11_4.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi11[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi21_1.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi21[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi21_2.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi21[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi21_3.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi21[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi21_4.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi21[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi31_1.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi31[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi31_2.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi31[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi31_3.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi31[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[311]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi31_4.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi31[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi41_1.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi41[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi41_2.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi41[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi41_3.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi41[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi41_4.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi41[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi22_1.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi22[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi22_2.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi22[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi22_3.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi22[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi22_4.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi22[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi32_1.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi32[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi32_2.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi32[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi32_3.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi32[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi32_4.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi32[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi42_1.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi42[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi42_2.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi42[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi42_3.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi42[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi42_4.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi42[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi33_1.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi33[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi33_2.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi33[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi33_3.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi33[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi33_4.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi33[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi43_1.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi43[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi43_2.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi43[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi43_3.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi43[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi43_4.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi43[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi44_1.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi44[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi44_2.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi44[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi44_3.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi44[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_mu_Psi44_4.pdf", width=9, height=7)
hist(ranks_nor_ref_mu_Psi44[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()
###########

chains_mu<-NULL
chains_Psi<-NULL
i_m<-0
while (i_m < M)
{chain<-sample_post_nor_ref_marg_Psi(X,U,Np+B)
print(sum(chain[[3]]))
if (sum(chain[[3]])>3000)
{chains_mu<-cbind(chains_mu,chain[[1]][,(B+1):(Np+B)])
chains_Psi<-cbind(chains_Psi,chain[[2]][,(B+1):(Np+B)])
i_m<-i_m+1
}
}

ranks_mu1<-rank(chains_mu[1,], ties.method = "average")
ranks_mu2<-rank(chains_mu[2,], ties.method = "average")
ranks_mu3<-rank(chains_mu[3,], ties.method = "average")
ranks_mu4<-rank(chains_mu[4,], ties.method = "average")
ranks_Psi11<-rank(chains_Psi[1,], ties.method = "average")
ranks_Psi21<-rank(chains_Psi[2,], ties.method = "average")
ranks_Psi31<-rank(chains_Psi[3,], ties.method = "average")
ranks_Psi41<-rank(chains_Psi[4,], ties.method = "average")
ranks_Psi22<-rank(chains_Psi[6,], ties.method = "average")
ranks_Psi32<-rank(chains_Psi[7,], ties.method = "average")
ranks_Psi42<-rank(chains_Psi[8,], ties.method = "average")
ranks_Psi33<-rank(chains_Psi[11,], ties.method = "average")
ranks_Psi43<-rank(chains_Psi[12,], ties.method = "average")
ranks_Psi44<-rank(chains_Psi[16,], ties.method = "average")

ranks_nor_ref_Psi_mu1<-matrix(ranks_mu1,Np,M)
ranks_nor_ref_Psi_mu2<-matrix(ranks_mu2,Np,M)
ranks_nor_ref_Psi_mu3<-matrix(ranks_mu3,Np,M)
ranks_nor_ref_Psi_mu4<-matrix(ranks_mu4,Np,M)
ranks_nor_ref_Psi_Psi11<-matrix(ranks_Psi11,Np,M)
ranks_nor_ref_Psi_Psi21<-matrix(ranks_Psi21,Np,M)
ranks_nor_ref_Psi_Psi31<-matrix(ranks_Psi31,Np,M)
ranks_nor_ref_Psi_Psi41<-matrix(ranks_Psi41,Np,M)
ranks_nor_ref_Psi_Psi22<-matrix(ranks_Psi22,Np,M)
ranks_nor_ref_Psi_Psi32<-matrix(ranks_Psi32,Np,M)
ranks_nor_ref_Psi_Psi42<-matrix(ranks_Psi42,Np,M)
ranks_nor_ref_Psi_Psi33<-matrix(ranks_Psi33,Np,M)
ranks_nor_ref_Psi_Psi43<-matrix(ranks_Psi43,Np,M)
ranks_nor_ref_Psi_Psi44<-matrix(ranks_Psi44,Np,M)

pdf(file="ranks_nor_ref_Psi_mu1_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_mu1[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_mu1_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_mu1[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_mu1_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_mu1[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_mu1_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_mu1[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_mu2_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_mu2[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_mu2_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_mu2[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_mu2_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_mu2[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_mu2_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_mu2[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()


pdf(file="ranks_nor_ref_Psi_mu3_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_mu3[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_mu3_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_mu3[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_mu3_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_mu3[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_mu3_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_mu3[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_mu4_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_mu4[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_mu4_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_mu4[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_mu4_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_mu4[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_mu4_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_mu4[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi11_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi11[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi11_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi11[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi11_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi11[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi11_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi11[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi21_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi21[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi21_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi21[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi21_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi21[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi21_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi21[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi31_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi31[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi31_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi31[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi31_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi31[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[311]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi31_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi31[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi41_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi41[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi41_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi41[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi41_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi41[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi41_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi41[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi22_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi22[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi22_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi22[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi22_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi22[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi22_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi22[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi32_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi32[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi32_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi32[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi32_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi32[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi32_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi32[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi42_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi42[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi42_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi42[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi42_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi42[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi42_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi42[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi33_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi33[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi33_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi33[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi33_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi33[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi33_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi33[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi43_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi43[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi43_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi43[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi43_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi43[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi43_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi43[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi44_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi44[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi44_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi44[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi44_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi44[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Psi_Psi44_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Psi_Psi44[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

####new

chains_mu<-NULL
chains_Psi<-NULL
i_m<-0
while (i_m < M)
{chain<-sample_post_nor_ref_Gibbs(X,U,Np+B)
print(sum(chain[[3]]))
if (sum(chain[[3]])>3000)
{chains_mu<-cbind(chains_mu,chain[[1]][,(B+1):(Np+B)])
chains_Psi<-cbind(chains_Psi,chain[[2]][,(B+1):(Np+B)])
i_m<-i_m+1
}
}

ranks_mu1<-rank(chains_mu[1,], ties.method = "average")
ranks_mu2<-rank(chains_mu[2,], ties.method = "average")
ranks_mu3<-rank(chains_mu[3,], ties.method = "average")
ranks_mu4<-rank(chains_mu[4,], ties.method = "average")
ranks_Psi11<-rank(chains_Psi[1,], ties.method = "average")
ranks_Psi21<-rank(chains_Psi[2,], ties.method = "average")
ranks_Psi31<-rank(chains_Psi[3,], ties.method = "average")
ranks_Psi41<-rank(chains_Psi[4,], ties.method = "average")
ranks_Psi22<-rank(chains_Psi[6,], ties.method = "average")
ranks_Psi32<-rank(chains_Psi[7,], ties.method = "average")
ranks_Psi42<-rank(chains_Psi[8,], ties.method = "average")
ranks_Psi33<-rank(chains_Psi[11,], ties.method = "average")
ranks_Psi43<-rank(chains_Psi[12,], ties.method = "average")
ranks_Psi44<-rank(chains_Psi[16,], ties.method = "average")

ranks_nor_ref_Gibbs_mu1<-matrix(ranks_mu1,Np,M)
ranks_nor_ref_Gibbs_mu2<-matrix(ranks_mu2,Np,M)
ranks_nor_ref_Gibbs_mu3<-matrix(ranks_mu3,Np,M)
ranks_nor_ref_Gibbs_mu4<-matrix(ranks_mu4,Np,M)
ranks_nor_ref_Gibbs_Psi11<-matrix(ranks_Psi11,Np,M)
ranks_nor_ref_Gibbs_Psi21<-matrix(ranks_Psi21,Np,M)
ranks_nor_ref_Gibbs_Psi31<-matrix(ranks_Psi31,Np,M)
ranks_nor_ref_Gibbs_Psi41<-matrix(ranks_Psi41,Np,M)
ranks_nor_ref_Gibbs_Psi22<-matrix(ranks_Psi22,Np,M)
ranks_nor_ref_Gibbs_Psi32<-matrix(ranks_Psi32,Np,M)
ranks_nor_ref_Gibbs_Psi42<-matrix(ranks_Psi42,Np,M)
ranks_nor_ref_Gibbs_Psi33<-matrix(ranks_Psi33,Np,M)
ranks_nor_ref_Gibbs_Psi43<-matrix(ranks_Psi43,Np,M)
ranks_nor_ref_Gibbs_Psi44<-matrix(ranks_Psi44,Np,M)

pdf(file="ranks_nor_ref_Gibbs_mu1_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_mu1[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_mu1_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_mu1[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_mu1_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_mu1[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_mu1_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_mu1[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_mu2_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_mu2[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_mu2_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_mu2[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_mu2_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_mu2[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_mu2_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_mu2[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()


pdf(file="ranks_nor_ref_Gibbs_mu3_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_mu3[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_mu3_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_mu3[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_mu3_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_mu3[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_mu3_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_mu3[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_mu4_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_mu4[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_mu4_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_mu4[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_mu4_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_mu4[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_mu4_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_mu4[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi11_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi11[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi11_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi11[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi11_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi11[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi11_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi11[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi21_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi21[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi21_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi21[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi21_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi21[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi21_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi21[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi31_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi31[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi31_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi31[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi31_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi31[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[311]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi31_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi31[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi41_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi41[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi41_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi41[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi41_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi41[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi41_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi41[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi22_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi22[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi22_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi22[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi22_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi22[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi22_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi22[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi32_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi32[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi32_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi32[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi32_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi32[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi32_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi32[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi42_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi42[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi42_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi42[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi42_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi42[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi42_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi42[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi33_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi33[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi33_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi33[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi33_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi33[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi33_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi33[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi43_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi43[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi43_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi43[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi43_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi43[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi43_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi43[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi44_1.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi44[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi44_2.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi44[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi44_3.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi44[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_ref_Gibbs_Psi44_4.pdf", width=9, height=7)
hist(ranks_nor_ref_Gibbs_Psi44[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

####jef, Gibbs
chains_mu<-NULL
chains_Psi<-NULL
i_m<-0
while (i_m < M)
{chain<-sample_post_nor_jef_Gibbs(X,U,Np+B)
print(sum(chain[[3]]))
if (sum(chain[[3]])>3000)
{chains_mu<-cbind(chains_mu,chain[[1]][,(B+1):(Np+B)])
chains_Psi<-cbind(chains_Psi,chain[[2]][,(B+1):(Np+B)])
i_m<-i_m+1
}
}

ranks_mu1<-rank(chains_mu[1,], ties.method = "average")
ranks_mu2<-rank(chains_mu[2,], ties.method = "average")
ranks_mu3<-rank(chains_mu[3,], ties.method = "average")
ranks_mu4<-rank(chains_mu[4,], ties.method = "average")
ranks_Psi11<-rank(chains_Psi[1,], ties.method = "average")
ranks_Psi21<-rank(chains_Psi[2,], ties.method = "average")
ranks_Psi31<-rank(chains_Psi[3,], ties.method = "average")
ranks_Psi41<-rank(chains_Psi[4,], ties.method = "average")
ranks_Psi22<-rank(chains_Psi[6,], ties.method = "average")
ranks_Psi32<-rank(chains_Psi[7,], ties.method = "average")
ranks_Psi42<-rank(chains_Psi[8,], ties.method = "average")
ranks_Psi33<-rank(chains_Psi[11,], ties.method = "average")
ranks_Psi43<-rank(chains_Psi[12,], ties.method = "average")
ranks_Psi44<-rank(chains_Psi[16,], ties.method = "average")

ranks_nor_jef_Gibbs_mu1<-matrix(ranks_mu1,Np,M)
ranks_nor_jef_Gibbs_mu2<-matrix(ranks_mu2,Np,M)
ranks_nor_jef_Gibbs_mu3<-matrix(ranks_mu3,Np,M)
ranks_nor_jef_Gibbs_mu4<-matrix(ranks_mu4,Np,M)
ranks_nor_jef_Gibbs_Psi11<-matrix(ranks_Psi11,Np,M)
ranks_nor_jef_Gibbs_Psi21<-matrix(ranks_Psi21,Np,M)
ranks_nor_jef_Gibbs_Psi31<-matrix(ranks_Psi31,Np,M)
ranks_nor_jef_Gibbs_Psi41<-matrix(ranks_Psi41,Np,M)
ranks_nor_jef_Gibbs_Psi22<-matrix(ranks_Psi22,Np,M)
ranks_nor_jef_Gibbs_Psi32<-matrix(ranks_Psi32,Np,M)
ranks_nor_jef_Gibbs_Psi42<-matrix(ranks_Psi42,Np,M)
ranks_nor_jef_Gibbs_Psi33<-matrix(ranks_Psi33,Np,M)
ranks_nor_jef_Gibbs_Psi43<-matrix(ranks_Psi43,Np,M)
ranks_nor_jef_Gibbs_Psi44<-matrix(ranks_Psi44,Np,M)

pdf(file="ranks_nor_jef_Gibbs_mu1_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_mu1[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_mu1_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_mu1[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_mu1_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_mu1[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_mu1_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_mu1[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_mu2_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_mu2[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_mu2_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_mu2[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_mu2_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_mu2[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_mu2_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_mu2[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()


pdf(file="ranks_nor_jef_Gibbs_mu3_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_mu3[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_mu3_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_mu3[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_mu3_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_mu3[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_mu3_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_mu3[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_mu4_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_mu4[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_mu4_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_mu4[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_mu4_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_mu4[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_mu4_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_mu4[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi11_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi11[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi11_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi11[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi11_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi11[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi11_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi11[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi21_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi21[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi21_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi21[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi21_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi21[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi21_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi21[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi31_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi31[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi31_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi31[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi31_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi31[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[311]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi31_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi31[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi41_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi41[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi41_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi41[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi41_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi41[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi41_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi41[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi22_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi22[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi22_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi22[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi22_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi22[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi22_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi22[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi32_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi32[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi32_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi32[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi32_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi32[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi32_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi32[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi42_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi42[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi42_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi42[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi42_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi42[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi42_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi42[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi33_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi33[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi33_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi33[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi33_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi33[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi33_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi33[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi43_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi43[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi43_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi43[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi43_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi43[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi43_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi43[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi44_1.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi44[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi44_2.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi44[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi44_3.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi44[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_nor_jef_Gibbs_Psi44_4.pdf", width=9, height=7)
hist(ranks_nor_jef_Gibbs_Psi44[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

###########################################################################
Ut<-U*(d-2)/d

chains_mu<-NULL
chains_Psi<-NULL
i_m<-0
while (i_m < M)
{chain<-sample_post_t_jef_marg_mu(X,Ut,d,Np+B)
print(sum(chain[[3]]))
if (sum(chain[[3]])>3000)
{chains_mu<-cbind(chains_mu,chain[[1]][,(B+1):(Np+B)])
chains_Psi<-cbind(chains_Psi,chain[[2]][,(B+1):(Np+B)])
i_m<-i_m+1
}
}

ranks_mu1<-rank(chains_mu[1,], ties.method = "average")
ranks_mu2<-rank(chains_mu[2,], ties.method = "average")
ranks_mu3<-rank(chains_mu[3,], ties.method = "average")
ranks_mu4<-rank(chains_mu[4,], ties.method = "average")
ranks_Psi11<-rank(chains_Psi[1,], ties.method = "average")
ranks_Psi21<-rank(chains_Psi[2,], ties.method = "average")
ranks_Psi31<-rank(chains_Psi[3,], ties.method = "average")
ranks_Psi41<-rank(chains_Psi[4,], ties.method = "average")
ranks_Psi22<-rank(chains_Psi[6,], ties.method = "average")
ranks_Psi32<-rank(chains_Psi[7,], ties.method = "average")
ranks_Psi42<-rank(chains_Psi[8,], ties.method = "average")
ranks_Psi33<-rank(chains_Psi[11,], ties.method = "average")
ranks_Psi43<-rank(chains_Psi[12,], ties.method = "average")
ranks_Psi44<-rank(chains_Psi[16,], ties.method = "average")

ranks_t_jef_mu_mu1<-matrix(ranks_mu1,Np,M)
ranks_t_jef_mu_mu2<-matrix(ranks_mu2,Np,M)
ranks_t_jef_mu_mu3<-matrix(ranks_mu3,Np,M)
ranks_t_jef_mu_mu4<-matrix(ranks_mu4,Np,M)
ranks_t_jef_mu_Psi11<-matrix(ranks_Psi11,Np,M)
ranks_t_jef_mu_Psi21<-matrix(ranks_Psi21,Np,M)
ranks_t_jef_mu_Psi31<-matrix(ranks_Psi31,Np,M)
ranks_t_jef_mu_Psi41<-matrix(ranks_Psi41,Np,M)
ranks_t_jef_mu_Psi22<-matrix(ranks_Psi22,Np,M)
ranks_t_jef_mu_Psi32<-matrix(ranks_Psi32,Np,M)
ranks_t_jef_mu_Psi42<-matrix(ranks_Psi42,Np,M)
ranks_t_jef_mu_Psi33<-matrix(ranks_Psi33,Np,M)
ranks_t_jef_mu_Psi43<-matrix(ranks_Psi43,Np,M)
ranks_t_jef_mu_Psi44<-matrix(ranks_Psi44,Np,M)

pdf(file="ranks_t_jef_mu_mu1_1.pdf", width=9, height=7)
hist(ranks_t_jef_mu_mu1[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_mu1_2.pdf", width=9, height=7)
hist(ranks_t_jef_mu_mu1[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_mu1_3.pdf", width=9, height=7)
hist(ranks_t_jef_mu_mu1[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_mu1_4.pdf", width=9, height=7)
hist(ranks_t_jef_mu_mu1[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_mu2_1.pdf", width=9, height=7)
hist(ranks_t_jef_mu_mu2[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_mu2_2.pdf", width=9, height=7)
hist(ranks_t_jef_mu_mu2[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_mu2_3.pdf", width=9, height=7)
hist(ranks_t_jef_mu_mu2[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_mu2_4.pdf", width=9, height=7)
hist(ranks_t_jef_mu_mu2[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()


pdf(file="ranks_t_jef_mu_mu3_1.pdf", width=9, height=7)
hist(ranks_t_jef_mu_mu3[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_mu3_2.pdf", width=9, height=7)
hist(ranks_t_jef_mu_mu3[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_mu3_3.pdf", width=9, height=7)
hist(ranks_t_jef_mu_mu3[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_mu3_4.pdf", width=9, height=7)
hist(ranks_t_jef_mu_mu3[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_mu4_1.pdf", width=9, height=7)
hist(ranks_t_jef_mu_mu4[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_mu4_2.pdf", width=9, height=7)
hist(ranks_t_jef_mu_mu4[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_mu4_3.pdf", width=9, height=7)
hist(ranks_t_jef_mu_mu4[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_mu4_4.pdf", width=9, height=7)
hist(ranks_t_jef_mu_mu4[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi11_1.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi11[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi11_2.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi11[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi11_3.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi11[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi11_4.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi11[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi21_1.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi21[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi21_2.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi21[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi21_3.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi21[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi21_4.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi21[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi31_1.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi31[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi31_2.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi31[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi31_3.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi31[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[311]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi31_4.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi31[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi41_1.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi41[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi41_2.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi41[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi41_3.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi41[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi41_4.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi41[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi22_1.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi22[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi22_2.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi22[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi22_3.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi22[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi22_4.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi22[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi32_1.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi32[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi32_2.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi32[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi32_3.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi32[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi32_4.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi32[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi42_1.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi42[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi42_2.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi42[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi42_3.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi42[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi42_4.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi42[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi33_1.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi33[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi33_2.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi33[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi33_3.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi33[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi33_4.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi33[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi43_1.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi43[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi43_2.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi43[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi43_3.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi43[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi43_4.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi43[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi44_1.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi44[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi44_2.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi44[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi44_3.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi44[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_mu_Psi44_4.pdf", width=9, height=7)
hist(ranks_t_jef_mu_Psi44[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()
###########

chains_mu<-NULL
chains_Psi<-NULL
i_m<-0
while (i_m < M)
{chain<-sample_post_t_jef_marg_Psi(X,Ut,d,Np+B)
print(sum(chain[[3]]))
if (sum(chain[[3]])>3000)
{chains_mu<-cbind(chains_mu,chain[[1]][,(B+1):(Np+B)])
chains_Psi<-cbind(chains_Psi,chain[[2]][,(B+1):(Np+B)])
i_m<-i_m+1
}
}

ranks_mu1<-rank(chains_mu[1,], ties.method = "average")
ranks_mu2<-rank(chains_mu[2,], ties.method = "average")
ranks_mu3<-rank(chains_mu[3,], ties.method = "average")
ranks_mu4<-rank(chains_mu[4,], ties.method = "average")
ranks_Psi11<-rank(chains_Psi[1,], ties.method = "average")
ranks_Psi21<-rank(chains_Psi[2,], ties.method = "average")
ranks_Psi31<-rank(chains_Psi[3,], ties.method = "average")
ranks_Psi41<-rank(chains_Psi[4,], ties.method = "average")
ranks_Psi22<-rank(chains_Psi[6,], ties.method = "average")
ranks_Psi32<-rank(chains_Psi[7,], ties.method = "average")
ranks_Psi42<-rank(chains_Psi[8,], ties.method = "average")
ranks_Psi33<-rank(chains_Psi[11,], ties.method = "average")
ranks_Psi43<-rank(chains_Psi[12,], ties.method = "average")
ranks_Psi44<-rank(chains_Psi[16,], ties.method = "average")

ranks_t_jef_Psi_mu1<-matrix(ranks_mu1,Np,M)
ranks_t_jef_Psi_mu2<-matrix(ranks_mu2,Np,M)
ranks_t_jef_Psi_mu3<-matrix(ranks_mu3,Np,M)
ranks_t_jef_Psi_mu4<-matrix(ranks_mu4,Np,M)
ranks_t_jef_Psi_Psi11<-matrix(ranks_Psi11,Np,M)
ranks_t_jef_Psi_Psi21<-matrix(ranks_Psi21,Np,M)
ranks_t_jef_Psi_Psi31<-matrix(ranks_Psi31,Np,M)
ranks_t_jef_Psi_Psi41<-matrix(ranks_Psi41,Np,M)
ranks_t_jef_Psi_Psi22<-matrix(ranks_Psi22,Np,M)
ranks_t_jef_Psi_Psi32<-matrix(ranks_Psi32,Np,M)
ranks_t_jef_Psi_Psi42<-matrix(ranks_Psi42,Np,M)
ranks_t_jef_Psi_Psi33<-matrix(ranks_Psi33,Np,M)
ranks_t_jef_Psi_Psi43<-matrix(ranks_Psi43,Np,M)
ranks_t_jef_Psi_Psi44<-matrix(ranks_Psi44,Np,M)

pdf(file="ranks_t_jef_Psi_mu1_1.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_mu1[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_mu1_2.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_mu1[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_mu1_3.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_mu1[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_mu1_4.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_mu1[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_mu2_1.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_mu2[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_mu2_2.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_mu2[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_mu2_3.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_mu2[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_mu2_4.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_mu2[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()


pdf(file="ranks_t_jef_Psi_mu3_1.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_mu3[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_mu3_2.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_mu3[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_mu3_3.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_mu3[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_mu3_4.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_mu3[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_mu4_1.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_mu4[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_mu4_2.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_mu4[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_mu4_3.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_mu4[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_mu4_4.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_mu4[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi11_1.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi11[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi11_2.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi11[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi11_3.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi11[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi11_4.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi11[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi21_1.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi21[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi21_2.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi21[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi21_3.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi21[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi21_4.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi21[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi31_1.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi31[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi31_2.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi31[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi31_3.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi31[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[311]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi31_4.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi31[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi41_1.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi41[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi41_2.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi41[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi41_3.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi41[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi41_4.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi41[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi22_1.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi22[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi22_2.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi22[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi22_3.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi22[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi22_4.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi22[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi32_1.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi32[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi32_2.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi32[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi32_3.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi32[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi32_4.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi32[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi42_1.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi42[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi42_2.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi42[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi42_3.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi42[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi42_4.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi42[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi33_1.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi33[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi33_2.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi33[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi33_3.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi33[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi33_4.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi33[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi43_1.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi43[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi43_2.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi43[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi43_3.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi43[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi43_4.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi43[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi44_1.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi44[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi44_2.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi44[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi44_3.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi44[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Psi_Psi44_4.pdf", width=9, height=7)
hist(ranks_t_jef_Psi_Psi44[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()
###########################################################################

###########################################################################

chains_mu<-NULL
chains_Psi<-NULL
i_m<-0
while (i_m < M)
{chain<-sample_post_t_ref_marg_mu(X,Ut,d,Np+B)
print(sum(chain[[3]]))
if (sum(chain[[3]])>3000)
{chains_mu<-cbind(chains_mu,chain[[1]][,(B+1):(Np+B)])
chains_Psi<-cbind(chains_Psi,chain[[2]][,(B+1):(Np+B)])
i_m<-i_m+1
}
}

ranks_mu1<-rank(chains_mu[1,], ties.method = "average")
ranks_mu2<-rank(chains_mu[2,], ties.method = "average")
ranks_mu3<-rank(chains_mu[3,], ties.method = "average")
ranks_mu4<-rank(chains_mu[4,], ties.method = "average")
ranks_Psi11<-rank(chains_Psi[1,], ties.method = "average")
ranks_Psi21<-rank(chains_Psi[2,], ties.method = "average")
ranks_Psi31<-rank(chains_Psi[3,], ties.method = "average")
ranks_Psi41<-rank(chains_Psi[4,], ties.method = "average")
ranks_Psi22<-rank(chains_Psi[6,], ties.method = "average")
ranks_Psi32<-rank(chains_Psi[7,], ties.method = "average")
ranks_Psi42<-rank(chains_Psi[8,], ties.method = "average")
ranks_Psi33<-rank(chains_Psi[11,], ties.method = "average")
ranks_Psi43<-rank(chains_Psi[12,], ties.method = "average")
ranks_Psi44<-rank(chains_Psi[16,], ties.method = "average")

ranks_t_ref_mu_mu1<-matrix(ranks_mu1,Np,M)
ranks_t_ref_mu_mu2<-matrix(ranks_mu2,Np,M)
ranks_t_ref_mu_mu3<-matrix(ranks_mu3,Np,M)
ranks_t_ref_mu_mu4<-matrix(ranks_mu4,Np,M)
ranks_t_ref_mu_Psi11<-matrix(ranks_Psi11,Np,M)
ranks_t_ref_mu_Psi21<-matrix(ranks_Psi21,Np,M)
ranks_t_ref_mu_Psi31<-matrix(ranks_Psi31,Np,M)
ranks_t_ref_mu_Psi41<-matrix(ranks_Psi41,Np,M)
ranks_t_ref_mu_Psi22<-matrix(ranks_Psi22,Np,M)
ranks_t_ref_mu_Psi32<-matrix(ranks_Psi32,Np,M)
ranks_t_ref_mu_Psi42<-matrix(ranks_Psi42,Np,M)
ranks_t_ref_mu_Psi33<-matrix(ranks_Psi33,Np,M)
ranks_t_ref_mu_Psi43<-matrix(ranks_Psi43,Np,M)
ranks_t_ref_mu_Psi44<-matrix(ranks_Psi44,Np,M)

pdf(file="ranks_t_ref_mu_mu1_1.pdf", width=9, height=7)
hist(ranks_t_ref_mu_mu1[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_mu1_2.pdf", width=9, height=7)
hist(ranks_t_ref_mu_mu1[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_mu1_3.pdf", width=9, height=7)
hist(ranks_t_ref_mu_mu1[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_mu1_4.pdf", width=9, height=7)
hist(ranks_t_ref_mu_mu1[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_mu2_1.pdf", width=9, height=7)
hist(ranks_t_ref_mu_mu2[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_mu2_2.pdf", width=9, height=7)
hist(ranks_t_ref_mu_mu2[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_mu2_3.pdf", width=9, height=7)
hist(ranks_t_ref_mu_mu2[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_mu2_4.pdf", width=9, height=7)
hist(ranks_t_ref_mu_mu2[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()


pdf(file="ranks_t_ref_mu_mu3_1.pdf", width=9, height=7)
hist(ranks_t_ref_mu_mu3[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_mu3_2.pdf", width=9, height=7)
hist(ranks_t_ref_mu_mu3[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_mu3_3.pdf", width=9, height=7)
hist(ranks_t_ref_mu_mu3[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_mu3_4.pdf", width=9, height=7)
hist(ranks_t_ref_mu_mu3[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_mu4_1.pdf", width=9, height=7)
hist(ranks_t_ref_mu_mu4[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_mu4_2.pdf", width=9, height=7)
hist(ranks_t_ref_mu_mu4[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_mu4_3.pdf", width=9, height=7)
hist(ranks_t_ref_mu_mu4[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_mu4_4.pdf", width=9, height=7)
hist(ranks_t_ref_mu_mu4[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi11_1.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi11[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi11_2.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi11[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi11_3.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi11[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi11_4.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi11[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi21_1.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi21[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi21_2.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi21[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi21_3.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi21[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi21_4.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi21[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi31_1.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi31[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi31_2.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi31[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi31_3.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi31[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[311]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi31_4.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi31[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi41_1.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi41[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi41_2.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi41[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi41_3.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi41[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi41_4.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi41[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi22_1.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi22[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi22_2.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi22[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi22_3.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi22[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi22_4.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi22[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi32_1.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi32[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi32_2.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi32[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi32_3.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi32[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi32_4.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi32[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi42_1.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi42[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi42_2.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi42[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi42_3.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi42[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi42_4.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi42[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi33_1.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi33[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi33_2.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi33[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi33_3.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi33[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi33_4.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi33[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi43_1.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi43[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi43_2.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi43[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi43_3.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi43[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi43_4.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi43[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi44_1.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi44[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi44_2.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi44[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi44_3.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi44[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_mu_Psi44_4.pdf", width=9, height=7)
hist(ranks_t_ref_mu_Psi44[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()
###########

chains_mu<-NULL
chains_Psi<-NULL
i_m<-0
while (i_m < M)
{chain<-sample_post_t_ref_marg_Psi(X,Ut,d,Np+B)
print(sum(chain[[3]]))
if (sum(chain[[3]])>3000)
{chains_mu<-cbind(chains_mu,chain[[1]][,(B+1):(Np+B)])
chains_Psi<-cbind(chains_Psi,chain[[2]][,(B+1):(Np+B)])
i_m<-i_m+1
}
}

ranks_mu1<-rank(chains_mu[1,], ties.method = "average")
ranks_mu2<-rank(chains_mu[2,], ties.method = "average")
ranks_mu3<-rank(chains_mu[3,], ties.method = "average")
ranks_mu4<-rank(chains_mu[4,], ties.method = "average")
ranks_Psi11<-rank(chains_Psi[1,], ties.method = "average")
ranks_Psi21<-rank(chains_Psi[2,], ties.method = "average")
ranks_Psi31<-rank(chains_Psi[3,], ties.method = "average")
ranks_Psi41<-rank(chains_Psi[4,], ties.method = "average")
ranks_Psi22<-rank(chains_Psi[6,], ties.method = "average")
ranks_Psi32<-rank(chains_Psi[7,], ties.method = "average")
ranks_Psi42<-rank(chains_Psi[8,], ties.method = "average")
ranks_Psi33<-rank(chains_Psi[11,], ties.method = "average")
ranks_Psi43<-rank(chains_Psi[12,], ties.method = "average")
ranks_Psi44<-rank(chains_Psi[16,], ties.method = "average")

ranks_t_ref_Psi_mu1<-matrix(ranks_mu1,Np,M)
ranks_t_ref_Psi_mu2<-matrix(ranks_mu2,Np,M)
ranks_t_ref_Psi_mu3<-matrix(ranks_mu3,Np,M)
ranks_t_ref_Psi_mu4<-matrix(ranks_mu4,Np,M)
ranks_t_ref_Psi_Psi11<-matrix(ranks_Psi11,Np,M)
ranks_t_ref_Psi_Psi21<-matrix(ranks_Psi21,Np,M)
ranks_t_ref_Psi_Psi31<-matrix(ranks_Psi31,Np,M)
ranks_t_ref_Psi_Psi41<-matrix(ranks_Psi41,Np,M)
ranks_t_ref_Psi_Psi22<-matrix(ranks_Psi22,Np,M)
ranks_t_ref_Psi_Psi32<-matrix(ranks_Psi32,Np,M)
ranks_t_ref_Psi_Psi42<-matrix(ranks_Psi42,Np,M)
ranks_t_ref_Psi_Psi33<-matrix(ranks_Psi33,Np,M)
ranks_t_ref_Psi_Psi43<-matrix(ranks_Psi43,Np,M)
ranks_t_ref_Psi_Psi44<-matrix(ranks_Psi44,Np,M)

pdf(file="ranks_t_ref_Psi_mu1_1.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_mu1[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_mu1_2.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_mu1[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_mu1_3.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_mu1[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_mu1_4.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_mu1[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_mu2_1.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_mu2[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_mu2_2.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_mu2[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_mu2_3.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_mu2[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_mu2_4.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_mu2[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()


pdf(file="ranks_t_ref_Psi_mu3_1.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_mu3[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_mu3_2.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_mu3[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_mu3_3.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_mu3[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_mu3_4.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_mu3[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_mu4_1.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_mu4[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_mu4_2.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_mu4[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_mu4_3.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_mu4[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_mu4_4.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_mu4[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi11_1.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi11[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi11_2.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi11[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi11_3.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi11[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi11_4.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi11[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi21_1.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi21[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi21_2.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi21[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi21_3.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi21[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi21_4.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi21[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi31_1.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi31[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi31_2.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi31[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi31_3.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi31[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[311]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi31_4.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi31[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi41_1.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi41[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi41_2.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi41[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi41_3.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi41[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi41_4.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi41[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi22_1.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi22[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi22_2.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi22[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi22_3.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi22[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi22_4.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi22[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi32_1.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi32[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi32_2.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi32[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi32_3.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi32[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi32_4.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi32[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi42_1.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi42[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi42_2.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi42[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi42_3.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi42[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi42_4.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi42[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi33_1.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi33[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi33_2.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi33[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi33_3.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi33[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi33_4.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi33[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi43_1.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi43[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi43_2.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi43[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi43_3.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi43[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi43_4.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi43[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi44_1.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi44[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi44_2.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi44[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi44_3.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi44[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Psi_Psi44_4.pdf", width=9, height=7)
hist(ranks_t_ref_Psi_Psi44[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

####new

chains_mu<-NULL
chains_Psi<-NULL
i_m<-0
while (i_m < M)
{chain<-sample_post_t_ref_Gibbs(X,Ut,d,Np+B)
print(sum(chain[[3]]))
if (sum(chain[[3]])>3000)
{chains_mu<-cbind(chains_mu,chain[[1]][,(B+1):(Np+B)])
chains_Psi<-cbind(chains_Psi,chain[[2]][,(B+1):(Np+B)])
i_m<-i_m+1
}
}

ranks_mu1<-rank(chains_mu[1,], ties.method = "average")
ranks_mu2<-rank(chains_mu[2,], ties.method = "average")
ranks_mu3<-rank(chains_mu[3,], ties.method = "average")
ranks_mu4<-rank(chains_mu[4,], ties.method = "average")
ranks_Psi11<-rank(chains_Psi[1,], ties.method = "average")
ranks_Psi21<-rank(chains_Psi[2,], ties.method = "average")
ranks_Psi31<-rank(chains_Psi[3,], ties.method = "average")
ranks_Psi41<-rank(chains_Psi[4,], ties.method = "average")
ranks_Psi22<-rank(chains_Psi[6,], ties.method = "average")
ranks_Psi32<-rank(chains_Psi[7,], ties.method = "average")
ranks_Psi42<-rank(chains_Psi[8,], ties.method = "average")
ranks_Psi33<-rank(chains_Psi[11,], ties.method = "average")
ranks_Psi43<-rank(chains_Psi[12,], ties.method = "average")
ranks_Psi44<-rank(chains_Psi[16,], ties.method = "average")

ranks_t_ref_Gibbs_mu1<-matrix(ranks_mu1,Np,M)
ranks_t_ref_Gibbs_mu2<-matrix(ranks_mu2,Np,M)
ranks_t_ref_Gibbs_mu3<-matrix(ranks_mu3,Np,M)
ranks_t_ref_Gibbs_mu4<-matrix(ranks_mu4,Np,M)
ranks_t_ref_Gibbs_Psi11<-matrix(ranks_Psi11,Np,M)
ranks_t_ref_Gibbs_Psi21<-matrix(ranks_Psi21,Np,M)
ranks_t_ref_Gibbs_Psi31<-matrix(ranks_Psi31,Np,M)
ranks_t_ref_Gibbs_Psi41<-matrix(ranks_Psi41,Np,M)
ranks_t_ref_Gibbs_Psi22<-matrix(ranks_Psi22,Np,M)
ranks_t_ref_Gibbs_Psi32<-matrix(ranks_Psi32,Np,M)
ranks_t_ref_Gibbs_Psi42<-matrix(ranks_Psi42,Np,M)
ranks_t_ref_Gibbs_Psi33<-matrix(ranks_Psi33,Np,M)
ranks_t_ref_Gibbs_Psi43<-matrix(ranks_Psi43,Np,M)
ranks_t_ref_Gibbs_Psi44<-matrix(ranks_Psi44,Np,M)

pdf(file="ranks_t_ref_Gibbs_mu1_1.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_mu1[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_mu1_2.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_mu1[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_mu1_3.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_mu1[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_mu1_4.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_mu1[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_mu2_1.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_mu2[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_mu2_2.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_mu2[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_mu2_3.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_mu2[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_mu2_4.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_mu2[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()


pdf(file="ranks_t_ref_Gibbs_mu3_1.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_mu3[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_mu3_2.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_mu3[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_mu3_3.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_mu3[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_mu3_4.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_mu3[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_mu4_1.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_mu4[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_mu4_2.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_mu4[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_mu4_3.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_mu4[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_mu4_4.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_mu4[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi11_1.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi11[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi11_2.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi11[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi11_3.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi11[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi11_4.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi11[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi21_1.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi21[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi21_2.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi21[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi21_3.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi21[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi21_4.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi21[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi31_1.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi31[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi31_2.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi31[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi31_3.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi31[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[311]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi31_4.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi31[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi41_1.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi41[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi41_2.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi41[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi41_3.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi41[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi41_4.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi41[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi22_1.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi22[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi22_2.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi22[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi22_3.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi22[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi22_4.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi22[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi32_1.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi32[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi32_2.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi32[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi32_3.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi32[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi32_4.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi32[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi42_1.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi42[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi42_2.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi42[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi42_3.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi42[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi42_4.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi42[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi33_1.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi33[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi33_2.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi33[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi33_3.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi33[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi33_4.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi33[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi43_1.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi43[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi43_2.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi43[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi43_3.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi43[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi43_4.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi43[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi44_1.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi44[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi44_2.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi44[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi44_3.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi44[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_ref_Gibbs_Psi44_4.pdf", width=9, height=7)
hist(ranks_t_ref_Gibbs_Psi44[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

####jef, Gibbs
chains_mu<-NULL
chains_Psi<-NULL
i_m<-0
while (i_m < M)
{chain<-sample_post_t_jef_Gibbs(X,Ut,d,Np+B)
print(sum(chain[[3]]))
if (sum(chain[[3]])>3000)
{chains_mu<-cbind(chains_mu,chain[[1]][,(B+1):(Np+B)])
chains_Psi<-cbind(chains_Psi,chain[[2]][,(B+1):(Np+B)])
i_m<-i_m+1
}
}

ranks_mu1<-rank(chains_mu[1,], ties.method = "average")
ranks_mu2<-rank(chains_mu[2,], ties.method = "average")
ranks_mu3<-rank(chains_mu[3,], ties.method = "average")
ranks_mu4<-rank(chains_mu[4,], ties.method = "average")
ranks_Psi11<-rank(chains_Psi[1,], ties.method = "average")
ranks_Psi21<-rank(chains_Psi[2,], ties.method = "average")
ranks_Psi31<-rank(chains_Psi[3,], ties.method = "average")
ranks_Psi41<-rank(chains_Psi[4,], ties.method = "average")
ranks_Psi22<-rank(chains_Psi[6,], ties.method = "average")
ranks_Psi32<-rank(chains_Psi[7,], ties.method = "average")
ranks_Psi42<-rank(chains_Psi[8,], ties.method = "average")
ranks_Psi33<-rank(chains_Psi[11,], ties.method = "average")
ranks_Psi43<-rank(chains_Psi[12,], ties.method = "average")
ranks_Psi44<-rank(chains_Psi[16,], ties.method = "average")

ranks_t_jef_Gibbs_mu1<-matrix(ranks_mu1,Np,M)
ranks_t_jef_Gibbs_mu2<-matrix(ranks_mu2,Np,M)
ranks_t_jef_Gibbs_mu3<-matrix(ranks_mu3,Np,M)
ranks_t_jef_Gibbs_mu4<-matrix(ranks_mu4,Np,M)
ranks_t_jef_Gibbs_Psi11<-matrix(ranks_Psi11,Np,M)
ranks_t_jef_Gibbs_Psi21<-matrix(ranks_Psi21,Np,M)
ranks_t_jef_Gibbs_Psi31<-matrix(ranks_Psi31,Np,M)
ranks_t_jef_Gibbs_Psi41<-matrix(ranks_Psi41,Np,M)
ranks_t_jef_Gibbs_Psi22<-matrix(ranks_Psi22,Np,M)
ranks_t_jef_Gibbs_Psi32<-matrix(ranks_Psi32,Np,M)
ranks_t_jef_Gibbs_Psi42<-matrix(ranks_Psi42,Np,M)
ranks_t_jef_Gibbs_Psi33<-matrix(ranks_Psi33,Np,M)
ranks_t_jef_Gibbs_Psi43<-matrix(ranks_Psi43,Np,M)
ranks_t_jef_Gibbs_Psi44<-matrix(ranks_Psi44,Np,M)

pdf(file="ranks_t_jef_Gibbs_mu1_1.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_mu1[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_mu1_2.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_mu1[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_mu1_3.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_mu1[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_mu1_4.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_mu1[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[1]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_mu2_1.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_mu2[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_mu2_2.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_mu2[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_mu2_3.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_mu2[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_mu2_4.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_mu2[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[2]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()


pdf(file="ranks_t_jef_Gibbs_mu3_1.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_mu3[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_mu3_2.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_mu3[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_mu3_3.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_mu3[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_mu3_4.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_mu3[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[3]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_mu4_1.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_mu4[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_mu4_2.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_mu4[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_mu4_3.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_mu4[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_mu4_4.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_mu4[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~mu[4]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi11_1.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi11[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi11_2.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi11[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi11_3.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi11[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi11_4.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi11[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[11]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi21_1.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi21[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi21_2.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi21[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi21_3.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi21[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi21_4.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi21[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[21]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi31_1.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi31[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi31_2.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi31[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi31_3.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi31[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[311]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi31_4.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi31[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[31]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi41_1.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi41[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi41_2.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi41[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi41_3.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi41[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi41_4.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi41[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[41]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi22_1.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi22[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi22_2.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi22[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi22_3.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi22[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi22_4.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi22[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[22]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi32_1.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi32[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi32_2.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi32[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi32_3.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi32[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi32_4.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi32[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[32]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi42_1.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi42[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi42_2.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi42[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi42_3.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi42[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi42_4.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi42[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[42]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi33_1.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi33[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi33_2.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi33[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi33_3.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi33[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi33_4.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi33[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[33]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi43_1.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi43[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi43_2.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi43[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi43_3.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi43[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi43_4.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi43[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[43]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi44_1.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi44[,1],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 1,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi44_2.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi44[,2],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 2,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi44_3.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi44[,3],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 3,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

pdf(file="ranks_t_jef_Gibbs_Psi44_4.pdf", width=9, height=7)
hist(ranks_t_jef_Gibbs_Psi44[,4],breaks=25,prob=TRUE, labels = F, border = "dark blue", col = "light blue", main = expression("Chain 4,"~Psi[44]), xlab = expression(), ylab = expression(),cex.axis=2.2,cex.main=3.0,font=2)
dev.off()

###########################################################################

############### hAtR ######################################################
hatR<-matrix(0,12,14)

##############################################
x_ranks_full<-qnorm((ranks_nor_jef_mu_mu1-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[1,1]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_mu_mu2-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[1,2]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_mu_mu3-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[1,3]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_mu_mu4-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[1,4]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_mu_Psi11-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[1,5]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_mu_Psi21-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[1,6]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_mu_Psi31-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[1,7]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_mu_Psi41-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[1,8]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_mu_Psi22-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[1,9]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_mu_Psi32-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[1,10]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_mu_Psi42-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[1,11]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_mu_Psi33-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[1,12]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_mu_Psi43-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[1,13]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_mu_Psi44-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[1,14]<-sqrt(1-2/Np+(2/Np)*BR/WR)
##############################################
x_ranks_full<-qnorm((ranks_nor_jef_Psi_mu1-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[2,1]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Psi_mu2-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[2,2]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Psi_mu3-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[2,3]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Psi_mu4-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[2,4]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Psi_Psi11-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[2,5]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Psi_Psi21-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[2,6]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Psi_Psi31-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[2,7]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Psi_Psi41-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[2,8]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Psi_Psi22-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[2,9]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Psi_Psi32-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[2,10]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Psi_Psi42-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[2,11]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Psi_Psi33-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[2,12]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Psi_Psi43-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[2,13]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Psi_Psi44-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[2,14]<-sqrt(1-2/Np+(2/Np)*BR/WR)
##############################################
x_ranks_full<-qnorm((ranks_nor_jef_Gibbs_mu1-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[3,1]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Gibbs_mu2-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[3,2]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Gibbs_mu3-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[3,3]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Gibbs_mu4-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[3,4]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Gibbs_Psi11-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[3,5]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Gibbs_Psi21-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[3,6]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Gibbs_Psi31-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[3,7]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Gibbs_Psi41-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[3,8]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Gibbs_Psi22-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[3,9]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Gibbs_Psi32-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[3,10]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Gibbs_Psi42-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[3,11]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Gibbs_Psi33-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[3,12]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Gibbs_Psi43-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[3,13]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_jef_Gibbs_Psi44-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[3,14]<-sqrt(1-2/Np+(2/Np)*BR/WR)
##############################################
x_ranks_full<-qnorm((ranks_nor_ref_mu_mu1-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[4,1]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_mu_mu2-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[4,2]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_mu_mu3-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[4,3]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_mu_mu4-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[4,4]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_mu_Psi11-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[4,5]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_mu_Psi21-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[4,6]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_mu_Psi31-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[4,7]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_mu_Psi41-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[4,8]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_mu_Psi22-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[4,9]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_mu_Psi32-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[4,10]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_mu_Psi42-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[4,11]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_mu_Psi33-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[4,12]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_mu_Psi43-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[4,13]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_mu_Psi44-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[4,14]<-sqrt(1-2/Np+(2/Np)*BR/WR)
##############################################
x_ranks_full<-qnorm((ranks_nor_ref_Psi_mu1-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[5,1]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Psi_mu2-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[5,2]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Psi_mu3-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[5,3]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Psi_mu4-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[5,4]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Psi_Psi11-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[5,5]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Psi_Psi21-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[5,6]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Psi_Psi31-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[5,7]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Psi_Psi41-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[5,8]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Psi_Psi22-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[5,9]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Psi_Psi32-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[5,10]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Psi_Psi42-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[5,11]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Psi_Psi33-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[5,12]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Psi_Psi43-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[5,13]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Psi_Psi44-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[5,14]<-sqrt(1-2/Np+(2/Np)*BR/WR)
##############################################
##############################################
x_ranks_full<-qnorm((ranks_nor_ref_Gibbs_mu1-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[6,1]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Gibbs_mu2-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[6,2]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Gibbs_mu3-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[6,3]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Gibbs_mu4-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[6,4]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Gibbs_Psi11-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[6,5]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Gibbs_Psi21-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[6,6]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Gibbs_Psi31-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[6,7]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Gibbs_Psi41-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[6,8]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Gibbs_Psi22-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[6,9]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Gibbs_Psi32-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[6,10]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Gibbs_Psi42-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[6,11]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Gibbs_Psi33-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[6,12]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Gibbs_Psi43-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[6,13]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_nor_ref_Gibbs_Psi44-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[6,14]<-sqrt(1-2/Np+(2/Np)*BR/WR)
##############################################

##############################################
x_ranks_full<-qnorm((ranks_t_jef_mu_mu1-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[7,1]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_mu_mu2-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[7,2]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_mu_mu3-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[7,3]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_mu_mu4-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[7,4]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_mu_Psi11-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[7,5]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_mu_Psi21-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[7,6]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_mu_Psi31-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[7,7]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_mu_Psi41-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[7,8]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_mu_Psi22-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[7,9]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_mu_Psi32-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[7,10]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_mu_Psi42-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[7,11]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_mu_Psi33-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[7,12]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_mu_Psi43-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[7,13]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_mu_Psi44-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[7,14]<-sqrt(1-2/Np+(2/Np)*BR/WR)
##############################################
x_ranks_full<-qnorm((ranks_t_jef_Psi_mu1-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[8,1]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Psi_mu2-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[8,2]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Psi_mu3-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[8,3]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Psi_mu4-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[8,4]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Psi_Psi11-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[8,5]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Psi_Psi21-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[8,6]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Psi_Psi31-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[8,7]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Psi_Psi41-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[8,8]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Psi_Psi22-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[8,9]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Psi_Psi32-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[8,10]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Psi_Psi42-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[8,11]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Psi_Psi33-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[8,12]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Psi_Psi43-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[8,13]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Psi_Psi44-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[8,14]<-sqrt(1-2/Np+(2/Np)*BR/WR)
##############################################
x_ranks_full<-qnorm((ranks_t_jef_Gibbs_mu1-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[9,1]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Gibbs_mu2-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[9,2]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Gibbs_mu3-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[9,3]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Gibbs_mu4-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[9,4]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Gibbs_Psi11-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[9,5]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Gibbs_Psi21-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[9,6]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Gibbs_Psi31-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[9,7]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Gibbs_Psi41-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[9,8]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Gibbs_Psi22-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[9,9]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Gibbs_Psi32-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[9,10]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Gibbs_Psi42-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[9,11]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Gibbs_Psi33-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[9,12]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Gibbs_Psi43-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[9,13]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_jef_Gibbs_Psi44-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[9,14]<-sqrt(1-2/Np+(2/Np)*BR/WR)
##############################################
x_ranks_full<-qnorm((ranks_t_ref_mu_mu1-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[10,1]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_mu_mu2-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[10,2]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_mu_mu3-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[10,3]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_mu_mu4-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[10,4]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_mu_Psi11-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[10,5]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_mu_Psi21-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[10,6]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_mu_Psi31-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[10,7]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_mu_Psi41-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[10,8]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_mu_Psi22-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[10,9]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_mu_Psi32-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[10,10]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_mu_Psi42-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[10,11]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_mu_Psi33-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[10,12]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_mu_Psi43-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[10,13]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_mu_Psi44-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[10,14]<-sqrt(1-2/Np+(2/Np)*BR/WR)
##############################################
x_ranks_full<-qnorm((ranks_t_ref_Psi_mu1-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[11,1]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Psi_mu2-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[11,2]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Psi_mu3-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[11,3]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Psi_mu4-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[11,4]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Psi_Psi11-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[11,5]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Psi_Psi21-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[11,6]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Psi_Psi31-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[11,7]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Psi_Psi41-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[11,8]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Psi_Psi22-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[11,9]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Psi_Psi32-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[11,10]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Psi_Psi42-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[11,11]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Psi_Psi33-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[11,12]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Psi_Psi43-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[11,13]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Psi_Psi44-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[11,14]<-sqrt(1-2/Np+(2/Np)*BR/WR)
##############################################
##############################################
x_ranks_full<-qnorm((ranks_t_ref_Gibbs_mu1-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[12,1]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Gibbs_mu2-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[12,2]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Gibbs_mu3-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[12,3]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Gibbs_mu4-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[12,4]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Gibbs_Psi11-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[12,5]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Gibbs_Psi21-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[12,6]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Gibbs_Psi31-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[12,7]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Gibbs_Psi41-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[12,8]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Gibbs_Psi22-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[12,9]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Gibbs_Psi32-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[12,10]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Gibbs_Psi42-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[12,11]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Gibbs_Psi33-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[12,12]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Gibbs_Psi43-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[12,13]<-sqrt(1-2/Np+(2/Np)*BR/WR)
######
x_ranks_full<-qnorm((ranks_t_ref_Gibbs_Psi44-3/8)/(M*Np+1/4))
x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
means<-apply(x_ranks,2,mean)
BR<-Np*var(means)/2
vars<-apply(x_ranks,2,var)
WR<-sum(vars)/(2*M)
hatR[12,14]<-sqrt(1-2/Np+(2/Np)*BR/WR)
##############################################

##############################################

write.table(hatR,file="hatR.csv", sep=",", dec=".", row.names = F,col.names = F)
print(xtable(hatR, type = "latex",digits =4), file = "hatR.tex")
print(round(hatR, digits = 4))

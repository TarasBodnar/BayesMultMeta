library('xtable')
library('vioplot')

set.seed(123)
source("MH_sample_post.R")

alp<-0.05
bet<-0.0001
d<-7 ## degrees of freedom in t multivariate random effects model
################### Determining X:p\times n and U:pn\times pn

set.seed(1234)
source("MH_sample_post.R")

alp<-0.05
Np<-10^5 # number of observations drawn from the proposal distribution
B<-10^4 # burn in 
d<-3 ## degrees of freedom in t multivariate random effects model
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

####################################Markov chains
Np<-10^5 # number of observations drawn from the proposal distribution
B<-Np/10 # burn in 

i_m<-0
while (i_m<0.17*Np)
{chain<-sample_post_nor_jef_marg_mu(X,U,Np+B)
print(sum(chain[[3]]))
i_m<-sum(chain[[3]])
}
chains_mu_nor_jef_Alg_A<-chain[[1]][,(B+1):(Np+B)]
chains_Psi_nor_jef_Alg_A<-chain[[2]][,(B+1):(Np+B)]

i_m<-0
while (i_m<0.17*Np)
{chain<-sample_post_nor_jef_marg_Psi(X,U,Np+B)
print(sum(chain[[3]]))
i_m<-sum(chain[[3]])
}
chains_mu_nor_jef_Alg_B<-chain[[1]][,(B+1):(Np+B)]
chains_Psi_nor_jef_Alg_B<-chain[[2]][,(B+1):(Np+B)]

i_m<-0
while (i_m<0.17*Np)
{chain<-sample_post_nor_jef_Gibbs(X,U,Np+B)
print(sum(chain[[3]]))
i_m<-sum(chain[[3]])
}
chains_mu_nor_jef_Alg_C<-chain[[1]][,(B+1):(Np+B)]
chains_Psi_nor_jef_Alg_C<-chain[[2]][,(B+1):(Np+B)]

i_m<-0
while (i_m<0.17*Np)
{chain<-sample_post_nor_ref_marg_mu(X,U,Np+B)
print(sum(chain[[3]]))
i_m<-sum(chain[[3]])
}
chains_mu_nor_ref_Alg_A<-chain[[1]][,(B+1):(Np+B)]
chains_Psi_nor_ref_Alg_A<-chain[[2]][,(B+1):(Np+B)]

i_m<-0
while (i_m<0.17*Np)
{chain<-sample_post_nor_ref_marg_Psi(X,U,Np+B)
print(sum(chain[[3]]))
i_m<-sum(chain[[3]])
}
chains_mu_nor_ref_Alg_B<-chain[[1]][,(B+1):(Np+B)]
chains_Psi_nor_ref_Alg_B<-chain[[2]][,(B+1):(Np+B)]

i_m<-0
while (i_m<0.17*Np)
{chain<-sample_post_nor_ref_Gibbs(X,U,Np+B)
print(sum(chain[[3]]))
i_m<-sum(chain[[3]])
}
chains_mu_nor_ref_Alg_C<-chain[[1]][,(B+1):(Np+B)]
chains_Psi_nor_ref_Alg_C<-chain[[2]][,(B+1):(Np+B)]

####### t distribution
Ut<-U*(d-2)/d
set.seed(321)
i_m<-0
while (i_m<0.17*Np)
{chain<-sample_post_t_jef_marg_mu(X,Ut,d,Np+B)
print(sum(chain[[3]]))
i_m<-sum(chain[[3]])
}
chains_mu_t_jef_Alg_A<-chain[[1]][,(B+1):(Np+B)]
chains_Psi_t_jef_Alg_A<-chain[[2]][,(B+1):(Np+B)]

i_m<-0
while (i_m<0.17*Np)
{chain<-sample_post_t_jef_marg_Psi(X,Ut,d,Np+B)
print(sum(chain[[3]]))
i_m<-sum(chain[[3]])
}
chains_mu_t_jef_Alg_B<-chain[[1]][,(B+1):(Np+B)]
chains_Psi_t_jef_Alg_B<-chain[[2]][,(B+1):(Np+B)]

i_m<-0
while (i_m<0.17*Np)
{chain<-sample_post_t_jef_Gibbs(X,Ut,d,Np+B)
print(sum(chain[[3]]))
i_m<-sum(chain[[3]])
}
chains_mu_t_jef_Alg_C<-chain[[1]][,(B+1):(Np+B)]
chains_Psi_t_jef_Alg_C<-chain[[2]][,(B+1):(Np+B)]

i_m<-0
while (i_m<0.17*Np)
{chain<-sample_post_t_ref_marg_mu(X,Ut,d,Np+B)
print(sum(chain[[3]]))
i_m<-sum(chain[[3]])
}
chains_mu_t_ref_Alg_A<-chain[[1]][,(B+1):(Np+B)]
chains_Psi_t_ref_Alg_A<-chain[[2]][,(B+1):(Np+B)]

i_m<-0
while (i_m<0.17*Np)
{chain<-sample_post_t_ref_marg_Psi(X,Ut,d,Np+B)
print(sum(chain[[3]]))
i_m<-sum(chain[[3]])
}
chains_mu_t_ref_Alg_B<-chain[[1]][,(B+1):(Np+B)]
chains_Psi_t_ref_Alg_B<-chain[[2]][,(B+1):(Np+B)]

i_m<-0
while (i_m<0.17*Np)
{chain<-sample_post_t_ref_Gibbs(X,Ut,d,Np+B)
print(sum(chain[[3]]))
i_m<-sum(chain[[3]])
}
chains_mu_t_ref_Alg_C<-chain[[1]][,(B+1):(Np+B)]
chains_Psi_t_ref_Alg_C<-chain[[2]][,(B+1):(Np+B)]

################################ plots with kernel densities
ind<- seq(50,Np,50)

pdf("mu1_violin_nor.pdf", height=5.5, width=7)
vioplot(chains_mu_nor_jef_Alg_A[1,ind],chains_mu_nor_jef_Alg_B[1,ind],chains_mu_nor_jef_Alg_C[1,ind],chains_mu_nor_ref_Alg_A[1,ind],chains_mu_nor_ref_Alg_B[1,ind],chains_mu_nor_ref_Alg_C[1,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("mu1_violin_t.pdf", height=5.5, width=7)
vioplot(chains_mu_t_jef_Alg_A[1,ind],chains_mu_t_jef_Alg_B[1,ind],chains_mu_t_jef_Alg_C[1,ind],chains_mu_t_ref_Alg_A[1,ind],chains_mu_t_ref_Alg_B[1,ind],chains_mu_t_ref_Alg_C[1,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("mu2_violin_nor.pdf", height=5.5, width=7)
vioplot(chains_mu_nor_jef_Alg_A[2,ind],chains_mu_nor_jef_Alg_B[2,ind],chains_mu_nor_jef_Alg_C[2,ind],chains_mu_nor_ref_Alg_A[2,ind],chains_mu_nor_ref_Alg_B[2,ind],chains_mu_nor_ref_Alg_C[2,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("mu2_violin_t.pdf", height=5.5, width=7)
vioplot(chains_mu_t_jef_Alg_A[2,ind],chains_mu_t_jef_Alg_B[2,ind],chains_mu_t_jef_Alg_C[2,ind],chains_mu_t_ref_Alg_A[2,ind],chains_mu_t_ref_Alg_B[2,ind],chains_mu_t_ref_Alg_C[2,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("mu3_violin_nor.pdf", height=5.5, width=7)
vioplot(chains_mu_nor_jef_Alg_A[3,ind],chains_mu_nor_jef_Alg_B[3,ind],chains_mu_nor_jef_Alg_C[3,ind],chains_mu_nor_ref_Alg_A[3,ind],chains_mu_nor_ref_Alg_B[3,ind],chains_mu_nor_ref_Alg_C[3,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("mu3_violin_t.pdf", height=5.5, width=7)
vioplot(chains_mu_t_jef_Alg_A[3,ind],chains_mu_t_jef_Alg_B[3,ind],chains_mu_t_jef_Alg_C[3,ind],chains_mu_t_ref_Alg_A[3,ind],chains_mu_t_ref_Alg_B[3,ind],chains_mu_t_ref_Alg_C[3,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("mu4_violin_nor.pdf", height=5.5, width=7)
vioplot(chains_mu_nor_jef_Alg_A[4,ind],chains_mu_nor_jef_Alg_B[4,ind],chains_mu_nor_jef_Alg_C[4,ind],chains_mu_nor_ref_Alg_A[4,ind],chains_mu_nor_ref_Alg_B[4,ind],chains_mu_nor_ref_Alg_C[4,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("mu4_violin_t.pdf", height=5.5, width=7)
vioplot(chains_mu_t_jef_Alg_A[4,ind],chains_mu_t_jef_Alg_B[4,ind],chains_mu_t_jef_Alg_C[4,ind],chains_mu_t_ref_Alg_A[4,ind],chains_mu_t_ref_Alg_B[4,ind],chains_mu_t_ref_Alg_C[4,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 



pdf("Psi11_violin_nor.pdf", height=5.5, width=7)
vioplot(chains_Psi_nor_jef_Alg_A[1,ind],chains_Psi_nor_jef_Alg_B[1,ind],chains_Psi_nor_jef_Alg_C[1,ind],chains_Psi_nor_ref_Alg_A[1,ind],chains_Psi_nor_ref_Alg_B[1,ind],chains_Psi_nor_ref_Alg_C[1,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("Psi11_violin_t.pdf", height=5.5, width=7)
vioplot(chains_Psi_t_jef_Alg_A[1,ind],chains_Psi_t_jef_Alg_B[1,ind],chains_Psi_t_jef_Alg_C[1,ind],chains_Psi_t_ref_Alg_A[1,ind],chains_Psi_t_ref_Alg_B[1,ind],chains_Psi_t_ref_Alg_C[1,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 


pdf("Psi21_violin_nor.pdf", height=5.5, width=7)
vioplot(chains_Psi_nor_jef_Alg_A[2,ind],chains_Psi_nor_jef_Alg_B[2,ind],chains_Psi_nor_jef_Alg_C[2,ind],chains_Psi_nor_ref_Alg_A[2,ind],chains_Psi_nor_ref_Alg_B[2,ind],chains_Psi_nor_ref_Alg_C[2,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("Psi21_violin_t.pdf", height=5.5, width=7)
vioplot(chains_Psi_t_jef_Alg_A[2,ind],chains_Psi_t_jef_Alg_B[2,ind],chains_Psi_t_jef_Alg_C[2,ind],chains_Psi_t_ref_Alg_A[2,ind],chains_Psi_t_ref_Alg_B[2,ind],chains_Psi_t_ref_Alg_C[2,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("Psi31_violin_nor.pdf", height=5.5, width=7)
vioplot(chains_Psi_nor_jef_Alg_A[3,ind],chains_Psi_nor_jef_Alg_B[3,ind],chains_Psi_nor_jef_Alg_C[3,ind],chains_Psi_nor_ref_Alg_A[3,ind],chains_Psi_nor_ref_Alg_B[3,ind],chains_Psi_nor_ref_Alg_C[3,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("Psi31_violin_t.pdf", height=5.5, width=7)
vioplot(chains_Psi_t_jef_Alg_A[3,ind],chains_Psi_t_jef_Alg_B[3,ind],chains_Psi_t_jef_Alg_C[3,ind],chains_Psi_t_ref_Alg_A[3,ind],chains_Psi_t_ref_Alg_B[3,ind],chains_Psi_t_ref_Alg_C[3,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("Psi41_violin_nor.pdf", height=5.5, width=7)
vioplot(chains_Psi_nor_jef_Alg_A[4,ind],chains_Psi_nor_jef_Alg_B[4,ind],chains_Psi_nor_jef_Alg_C[4,ind],chains_Psi_nor_ref_Alg_A[4,ind],chains_Psi_nor_ref_Alg_B[4,ind],chains_Psi_nor_ref_Alg_C[4,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("Psi41_violin_t.pdf", height=5.5, width=7)
vioplot(chains_Psi_t_jef_Alg_A[4,ind],chains_Psi_t_jef_Alg_B[4,ind],chains_Psi_t_jef_Alg_C[4,ind],chains_Psi_t_ref_Alg_A[4,ind],chains_Psi_t_ref_Alg_B[4,ind],chains_Psi_t_ref_Alg_C[4,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off()

pdf("Psi22_violin_nor.pdf", height=5.5, width=7)
vioplot(chains_Psi_nor_jef_Alg_A[6,ind],chains_Psi_nor_jef_Alg_B[6,ind],chains_Psi_nor_jef_Alg_C[6,ind],chains_Psi_nor_ref_Alg_A[6,ind],chains_Psi_nor_ref_Alg_B[6,ind],chains_Psi_nor_ref_Alg_C[6,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("Psi22_violin_t.pdf", height=5.5, width=7)
vioplot(chains_Psi_t_jef_Alg_A[6,ind],chains_Psi_t_jef_Alg_B[6,ind],chains_Psi_t_jef_Alg_C[6,ind],chains_Psi_t_ref_Alg_A[6,ind],chains_Psi_t_ref_Alg_B[6,ind],chains_Psi_t_ref_Alg_C[6,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("Psi32_violin_nor.pdf", height=5.5, width=7)
vioplot(chains_Psi_nor_jef_Alg_A[7,ind],chains_Psi_nor_jef_Alg_B[7,ind],chains_Psi_nor_jef_Alg_C[7,ind],chains_Psi_nor_ref_Alg_A[7,ind],chains_Psi_nor_ref_Alg_B[7,ind],chains_Psi_nor_ref_Alg_C[7,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("Psi32_violin_t.pdf", height=5.5, width=7)
vioplot(chains_Psi_t_jef_Alg_A[7,ind],chains_Psi_t_jef_Alg_B[7,ind],chains_Psi_t_jef_Alg_C[7,ind],chains_Psi_t_ref_Alg_A[7,ind],chains_Psi_t_ref_Alg_B[7,ind],chains_Psi_t_ref_Alg_C[7,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("Psi42_violin_nor.pdf", height=5.5, width=7)
vioplot(chains_Psi_nor_jef_Alg_A[8,ind],chains_Psi_nor_jef_Alg_B[8,ind],chains_Psi_nor_jef_Alg_C[8,ind],chains_Psi_nor_ref_Alg_A[8,ind],chains_Psi_nor_ref_Alg_B[8,ind],chains_Psi_nor_ref_Alg_C[8,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("Psi42_violin_t.pdf", height=5.5, width=7)
vioplot(chains_Psi_t_jef_Alg_A[8,ind],chains_Psi_t_jef_Alg_B[8,ind],chains_Psi_t_jef_Alg_C[8,ind],chains_Psi_t_ref_Alg_A[8,ind],chains_Psi_t_ref_Alg_B[8,ind],chains_Psi_t_ref_Alg_C[8,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off()

pdf("Psi33_violin_nor.pdf", height=5.5, width=7)
vioplot(chains_Psi_nor_jef_Alg_A[11,ind],chains_Psi_nor_jef_Alg_B[11,ind],chains_Psi_nor_jef_Alg_C[11,ind],chains_Psi_nor_ref_Alg_A[11,ind],chains_Psi_nor_ref_Alg_B[11,ind],chains_Psi_nor_ref_Alg_C[11,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("Psi33_violin_t.pdf", height=5.5, width=7)
vioplot(chains_Psi_t_jef_Alg_A[11,ind],chains_Psi_t_jef_Alg_B[11,ind],chains_Psi_t_jef_Alg_C[11,ind],chains_Psi_t_ref_Alg_A[11,ind],chains_Psi_t_ref_Alg_B[11,ind],chains_Psi_t_ref_Alg_C[11,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("Psi43_violin_nor.pdf", height=5.5, width=7)
vioplot(chains_Psi_nor_jef_Alg_A[12,ind],chains_Psi_nor_jef_Alg_B[12,ind],chains_Psi_nor_jef_Alg_C[12,ind],chains_Psi_nor_ref_Alg_A[12,ind],chains_Psi_nor_ref_Alg_B[12,ind],chains_Psi_nor_ref_Alg_C[12,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("Psi43_violin_t.pdf", height=5.5, width=7)
vioplot(chains_Psi_t_jef_Alg_A[12,ind],chains_Psi_t_jef_Alg_B[12,ind],chains_Psi_t_jef_Alg_C[12,ind],chains_Psi_t_ref_Alg_A[12,ind],chains_Psi_t_ref_Alg_B[12,ind],chains_Psi_t_ref_Alg_C[12,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("Psi44_violin_nor.pdf", height=5.5, width=7)
vioplot(chains_Psi_nor_jef_Alg_A[16,ind],chains_Psi_nor_jef_Alg_B[16,ind],chains_Psi_nor_jef_Alg_C[16,ind],chains_Psi_nor_ref_Alg_A[16,ind],chains_Psi_nor_ref_Alg_B[16,ind],chains_Psi_nor_ref_Alg_C[16,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

pdf("Psi44_violin_t.pdf", height=5.5, width=7)
vioplot(chains_Psi_t_jef_Alg_A[16,ind],chains_Psi_t_jef_Alg_B[16,ind],chains_Psi_t_jef_Alg_C[16,ind],chains_Psi_t_ref_Alg_A[16,ind],chains_Psi_t_ref_Alg_B[16,ind],chains_Psi_t_ref_Alg_C[16,ind],col=c("black","red","green","darkgrey","orange","blue"),xlab="",ylab = "", main="",xaxt = "n")
legend("bottom",fill=c("black","red","green","darkgrey","orange","blue"), legend=c("MH1, Jeffreys","MH2, Jeffreys","Gibbs, Jeffreys","MH1, reference","MH2, reference","Gibbs, reference"),ncol=3, inset=c(0, -.2), xpd=TRUE)
dev.off() 

############# Bayes inference
alp<-0.05
bet<-0.0001
#### mu
res_mu_nor_jef_Alg_A<-Bayes_inference(chains_mu_nor_jef_Alg_A[,ind],alp)
res_mu_nor_jef_Alg_B<-Bayes_inference(chains_mu_nor_jef_Alg_B[,ind],alp)
res_mu_nor_jef_Alg_C<-Bayes_inference(chains_mu_nor_jef_Alg_C[,ind],alp)

res_mu_nor_ref_Alg_A<-Bayes_inference(chains_mu_nor_ref_Alg_A[,ind],alp)
res_mu_nor_ref_Alg_B<-Bayes_inference(chains_mu_nor_ref_Alg_B[,ind],alp)
res_mu_nor_ref_Alg_C<-Bayes_inference(chains_mu_nor_ref_Alg_C[,ind],alp)

res_mu_t_jef_Alg_A<-Bayes_inference(chains_mu_t_jef_Alg_A[,ind],alp)
res_mu_t_jef_Alg_B<-Bayes_inference(chains_mu_t_jef_Alg_B[,ind],alp)
res_mu_t_jef_Alg_C<-Bayes_inference(chains_mu_t_jef_Alg_C[,ind],alp)

res_mu_t_ref_Alg_A<-Bayes_inference(chains_mu_t_ref_Alg_A[,ind],alp)
res_mu_t_ref_Alg_B<-Bayes_inference(chains_mu_t_ref_Alg_B[,ind],alp)
res_mu_t_ref_Alg_C<-Bayes_inference(chains_mu_t_ref_Alg_C[,ind],alp)

#### Psi
res_Psi_nor_jef_Alg_A<-Bayes_inference_Psi(chains_Psi_nor_jef_Alg_A[,ind],alp,bet)
res_Psi_nor_jef_Alg_B<-Bayes_inference_Psi(chains_Psi_nor_jef_Alg_B[,ind],alp,bet)
res_Psi_nor_jef_Alg_C<-Bayes_inference_Psi(chains_Psi_nor_jef_Alg_C[,ind],alp,bet)

res_Psi_nor_ref_Alg_A<-Bayes_inference_Psi(chains_Psi_nor_ref_Alg_A[,ind],alp,bet)
res_Psi_nor_ref_Alg_B<-Bayes_inference_Psi(chains_Psi_nor_ref_Alg_B[,ind],alp,bet)
res_Psi_nor_ref_Alg_C<-Bayes_inference_Psi(chains_Psi_nor_ref_Alg_C[,ind],alp,bet)

res_Psi_t_jef_Alg_A<-Bayes_inference_Psi(chains_Psi_t_jef_Alg_A[,ind],alp,bet)
res_Psi_t_jef_Alg_B<-Bayes_inference_Psi(chains_Psi_t_jef_Alg_B[,ind],alp,bet)
res_Psi_t_jef_Alg_C<-Bayes_inference_Psi(chains_Psi_t_jef_Alg_C[,ind],alp,bet)

res_Psi_t_ref_Alg_A<-Bayes_inference_Psi(chains_Psi_t_ref_Alg_A[,ind],alp,bet)
res_Psi_t_ref_Alg_B<-Bayes_inference_Psi(chains_Psi_t_ref_Alg_B[,ind],alp,bet)
res_Psi_t_ref_Alg_C<-Bayes_inference_Psi(chains_Psi_t_ref_Alg_C[,ind],alp,bet)

library('xtable')

#######################################################################################################
Bayes_result<-round(rbind(cbind(res_mu_nor_jef_Alg_A,res_Psi_nor_jef_Alg_A),
                          cbind(res_mu_nor_jef_Alg_B,res_Psi_nor_jef_Alg_B),
                          cbind(res_mu_nor_jef_Alg_C,res_Psi_nor_jef_Alg_C),
                          cbind(res_mu_nor_ref_Alg_A,res_Psi_nor_ref_Alg_A),
                          cbind(res_mu_nor_ref_Alg_B,res_Psi_nor_ref_Alg_B),
                          cbind(res_mu_nor_ref_Alg_C,res_Psi_nor_ref_Alg_C),
                          cbind(res_mu_t_jef_Alg_A,res_Psi_t_jef_Alg_A),
                          cbind(res_mu_t_jef_Alg_B,res_Psi_t_jef_Alg_B),
                          cbind(res_mu_t_jef_Alg_C,res_Psi_t_jef_Alg_C),
                          cbind(res_mu_t_ref_Alg_A,res_Psi_t_ref_Alg_A),
                          cbind(res_mu_t_ref_Alg_B,res_Psi_t_ref_Alg_B),
                          cbind(res_mu_t_ref_Alg_C,res_Psi_t_ref_Alg_C)),digits = 2)

Bayes_result1<-rbind(cbind(res_mu_nor_jef_Alg_A,res_Psi_nor_jef_Alg_A),
                           cbind(res_mu_nor_jef_Alg_B,res_Psi_nor_jef_Alg_B),
                           cbind(res_mu_nor_jef_Alg_C,res_Psi_nor_jef_Alg_C),
                           cbind(res_mu_nor_ref_Alg_A,res_Psi_nor_ref_Alg_A),
                           cbind(res_mu_nor_ref_Alg_B,res_Psi_nor_ref_Alg_B),
                           cbind(res_mu_nor_ref_Alg_C,res_Psi_nor_ref_Alg_C),
                           cbind(res_mu_t_jef_Alg_A,res_Psi_t_jef_Alg_A),
                           cbind(res_mu_t_jef_Alg_B,res_Psi_t_jef_Alg_B),
                           cbind(res_mu_t_jef_Alg_C,res_Psi_t_jef_Alg_C),
                           cbind(res_mu_t_ref_Alg_A,res_Psi_t_ref_Alg_A),
                           cbind(res_mu_t_ref_Alg_B,res_Psi_t_ref_Alg_B),
                           cbind(res_mu_t_ref_Alg_C,res_Psi_t_ref_Alg_C))
print(xtable(Bayes_result, type = "latex",digits =2), file = "Bayes_result.tex")

write.table(Bayes_result1,file="Bayes_result1.csv", sep=",", dec=".", row.names = F,col.names = F)

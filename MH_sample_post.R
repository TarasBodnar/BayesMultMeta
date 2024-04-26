#########################################################################################################
############## Duplication matrix #######################################################################
#########################################################################################################
Dp <- function(x){
  mat <- diag(x)
  index <- seq(x*(x+1)/2)
  mat[ lower.tri( mat , TRUE ) ] <- index
  mat[ upper.tri( mat ) ] <- t( mat )[ upper.tri( mat ) ]
  outer(c(mat), index , function( x , y ) ifelse(x==y, 1, 0 ) )
}

#########################################################################################################
############## Baesian inference #######################################################################
#########################################################################################################
Bayes_inference <- function(x,alp){
  x_mean<-apply(x,1,mean)
  x_med<-apply(x,1,quantile,probs=0.5)
  x_ql<-apply(x,1,quantile,probs=alp/2)
  x_qu<-apply(x,1,quantile,probs=1-alp/2)
  x_sd<-sqrt(apply(x,1,var))
  rbind(x_mean,x_med,x_sd,x_ql,x_qu)  
}

#########################################################################################################
### algorithmA, Jeffreys normal, X:p\times n data matrix, U: pn \times pn matrix with uncertainties #####
#########################################################################################################
sample_post_nor_jef_marg_mu<-function(X,U,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size
  
  ############## addtional definitons  
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-Dp(p)
  tGp<-t(Gp)
  
  trans_m<-NULL
  mu_m<-NULL
  Psi_m<-NULL
  
  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  tcholS<-t(chol(S))
  
  ### generating an initial value for mu and Psi new draw from proposal
  mu0<-bar_X+sqrt((n-1)/n/(n-p+1))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p+1)/(n-p+1))
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)
  cCov_p<-chol(Cov_p)
  Z<-matrix(rnorm(p*n+p),p,n+1)
  Psi0<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p
  
  ### initial value of the proposal
  q0<-det(Psi0)^(-0.5*(n+p+2))*exp(-0.5*sum(diag(solve(Psi0)%*%Cov_p)))
  
  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  mat_numJ<-matrix(0,p,p)
  det_num<-1
  exp_num<-1
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  mat_numJ<-mat_numJ+iPsiUi
  det_num<-det_num*sqrt(det(iPsiUi))
  exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
  }
  p0<-sqrt(det(mat_numJ/n))*sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num*exp_num
  
  for (j in 1:Np)
  {
    ### generating a new draw Psi from the proposal 
    mu_p<-bar_X+sqrt((n-1)/n/(n-p+1))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p+1)/(n-p+1))
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)
    cCov_p<-chol(Cov_p)
    Z<-matrix(rnorm(p*n+p),p,n+1)
    Psi_p<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p
    
    ### value of the proposal at new draw
    q1 <-det(Psi_p)^(-0.5*(n+p+2))*exp(-0.5*sum(diag(solve(Psi_p)%*%Cov_p)))
    
    ### value of the posterior
    mat_num<-matrix(0,p^2,p^2)
    mat_numJ<-matrix(0,p,p)
    det_num<-1
    exp_num<-1
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    mat_numJ<-mat_numJ+iPsiUi
    det_num<-det_num*sqrt(det(iPsiUi))
    exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
    }
    p1<-sqrt(det(mat_numJ/n))*sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num*exp_num
    
    ### MH ratio
    ratio_MH<-p1*q0/p0/q1
    
    trans_p<-0
    if (runif(1) <= ratio_MH)
    {trans_p<-1
    mu0<-mu_p
    Psi0<-Psi_p
    q0<-q1
    p0<-p1
    }
    
    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu0)
    trans_m<-cbind(trans_m,trans_p)  
  }
  output<-list(mu_m,Psi_m,trans_m)  
}

#########################################################################################################
### algorithmB, Jeffreys normal, X:p\times n data matrix, U: pn \times pn matrix with uncertainties #####
#########################################################################################################
sample_post_nor_jef_marg_Psi<-function(X,U,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size
  
  ############## addtional definitons  
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-Dp(p)
  tGp<-t(Gp)
  
  trans_m<-NULL
  mu_m<-NULL
  Psi_m<-NULL
  
  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  cS<-sqrt(n-1)*chol(S)
  
  ### generating an initial value for mu and Psi new draw from proposal
  Z<-matrix(rnorm(p*n),p,n)
  Psi0<-t(cS)%*%solve(Z%*%t(Z))%*%cS
  mu0<-bar_X+t(chol(Psi0))%*%rnorm(p)/sqrt(n)
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)
  
  ### initial value of the proposal
  exp_prop<-exp(-0.5*sum(diag(solve(Psi0)%*%Cov_p)))
  q0<-det(Psi0)^(-0.5*(n+p+2))*exp_prop
  
  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  mat_numJ<-matrix(0,p,p)
  det_num<-1
  exp_num<-1
  #exp_prop<-1
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  mat_numJ<-mat_numJ+iPsiUi
  det_num<-det_num*sqrt(det(iPsiUi))
  exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
  }
  p0<-sqrt(det(mat_numJ/n))*sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num*exp_num
  
  for (j in 1:Np)
  {
    ### generating a new draw Psi from the proposal 
    Z<-matrix(rnorm(p*n),p,n)
    Psi_p<-t(cS)%*%solve(Z%*%t(Z))%*%cS
    mu_p<-bar_X+t(chol(Psi_p))%*%rnorm(p)/sqrt(n)
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)
    
    ### value of the proposal at new draw
    exp_prop<-exp(-0.5*sum(diag(solve(Psi_p)%*%Cov_p)))
    q1<-det(Psi_p)^(-0.5*(n+p+2))*exp_prop
    
    ### value of the posterior
    mat_num<-matrix(0,p^2,p^2)
    mat_numJ<-matrix(0,p,p)
    det_num<-1
    exp_num<-1
    #exp_prop<-1
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    mat_numJ<-mat_numJ+iPsiUi
    det_num<-det_num*sqrt(det(iPsiUi))
    exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
    }
    p1<-sqrt(det(mat_numJ/n))*sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num*exp_num
    
    ### MH ratio
    ratio_MH<-p1*q0/p0/q1
    
    trans_p<-0
    if (runif(1) <= ratio_MH)
    {trans_p<-1
    mu0<-mu_p
    Psi0<-Psi_p
    q0<-q1
    p0<-p1
    }
    
    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu0)
    trans_m<-cbind(trans_m,trans_p)  
  }
  output<-list(mu_m,Psi_m,trans_m)  
}

#########################################################################################################
### algorithmA, Jeffreys t, X:p\times n data matrix, U: pn \times pn matrix with uncertainties ##########
#########################################################################################################
sample_post_t_jef_marg_mu<-function(X,U,d,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size
  
  ############## addtional definitons  
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-Dp(p)
  tGp<-t(Gp)
  
  trans_m<-NULL
  mu_m<-NULL
  Psi_m<-NULL
  
  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  tcholS<-t(chol(S))
  
  ### generating an initial value for mu and Psi new draw from proposal
  mu0<-bar_X+sqrt((n-1)/n/(n-p+1))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p+1)/(n-p+1))
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)
  cCov_p<-chol(Cov_p)
  Z<-matrix(rnorm(p*n+p),p,n+1)
  Psi0<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p*sqrt(rchisq(1,d)/(d))
  
  ### initial value of the proposal
  q0<-det(Psi0)^(-0.5*(n+p+2))*(1+sum(diag(solve(Psi0)%*%Cov_p))/d)^(-0.5*(p*n+d))
  
  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  mat_numJ<-matrix(0,p,p)
  det_num<-1
  exp_num<-1
  t_num<-0
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  mat_numJ<-mat_numJ+iPsiUi
  det_num<-det_num*sqrt(det(iPsiUi))
  t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
  }
  p0<-sqrt(det(mat_numJ/n))*sqrt(det(tGp%*%(mat_num-as.vector(mat_numJ)%*%t(as.vector(mat_numJ))/(n*p+d))%*%Gp/n))*det_num*(1+t_num)^(-0.5*(p*n+d))
  
  for (j in 1:Np)
  {
    ### generating a new draw Psi from the proposal 
    mu_p<-bar_X+sqrt((n-1)/n/(n-p+1))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p+1)/(n-p+1))
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)
    cCov_p<-chol(Cov_p)
    Z<-matrix(rnorm(p*n+p),p,n+1)
    Psi_p<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p*sqrt(rchisq(1,d)/(d))
    
    ### value of the proposal at new draw
    q1 <-det(Psi_p)^(-0.5*(n+p+2))*(1+sum(diag(solve(Psi_p)%*%Cov_p))/d)^(-0.5*(p*n+d))
    
    ### value of the posterior
    mat_num<-matrix(0,p^2,p^2)
    mat_numJ<-matrix(0,p,p)
    det_num<-1
    exp_num<-1
    t_num<-0
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    mat_numJ<-mat_numJ+iPsiUi
    det_num<-det_num*sqrt(det(iPsiUi))
    t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
    }
    p1<-sqrt(det(mat_numJ/n))*sqrt(det(tGp%*%(mat_num-as.vector(mat_numJ)%*%t(as.vector(mat_numJ))/(n*p+d))%*%Gp/n))*det_num*(1+t_num)^(-0.5*(p*n+d))
    
    ### MH ratio
    ratio_MH<-p1*q0/p0/q1
    
    trans_p<-0
    if (runif(1) <= ratio_MH)
    {trans_p<-1
    mu0<-mu_p
    Psi0<-Psi_p
    q0<-q1
    p0<-p1
    }
    
    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu0)
    trans_m<-cbind(trans_m,trans_p)  
  }
  output<-list(mu_m,Psi_m,trans_m)  
}

#########################################################################################################
### algorithmB, Jeffreys t, X:p\times n data matrix, U: pn \times pn matrix with uncertainties ##########
#########################################################################################################
sample_post_t_jef_marg_Psi<-function(X,U,d,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size
  
  ############## addtional definitons  
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-Dp(p)
  tGp<-t(Gp)
  
  trans_m<-NULL
  mu_m<-NULL
  Psi_m<-NULL
  
  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  cS<-sqrt(n-1)*chol(S)
  
  ### generating an initial value for mu and Psi new draw from proposal
  Z<-matrix(rnorm(p*n),p,n)
  Psi0<-t(cS)%*%solve(Z%*%t(Z))%*%cS*sqrt(rchisq(1,d)/(d))
  mu0<-bar_X+sqrt((d+(n-1)*sum(diag(solve(Psi0)%*%S)))/n/(d+p*n-p))*t(chol(Psi0))%*%rnorm(p)/sqrt(rchisq(1,d+p*n-p)/(d+p*n-p))
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)
  
  ### initial value of the proposal
  q0<-det(Psi0)^(-0.5*(n+p+2))*(1+sum(diag(solve(Psi0)%*%Cov_p))/d)^(-0.5*(p*n+d))
  
  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  mat_numJ<-matrix(0,p,p)
  det_num<-1
  exp_num<-1
  t_num<-0
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  mat_numJ<-mat_numJ+iPsiUi
  det_num<-det_num*sqrt(det(iPsiUi))
  t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
  }
  p0<-sqrt(det(mat_numJ/n))*sqrt(det(tGp%*%(mat_num-as.vector(mat_numJ)%*%t(as.vector(mat_numJ))/(n*p+d))%*%Gp/n))*det_num*(1+t_num)^(-0.5*(p*n+d))
  
  for (j in 1:Np)
  {
    ### generating a new draw Psi from the proposal 
    Z<-matrix(rnorm(p*n),p,n)
    Psi_p<-t(cS)%*%solve(Z%*%t(Z))%*%cS*sqrt(rchisq(1,d)/(d))
    mu_p<-bar_X+sqrt((d+(n-1)*sum(diag(solve(Psi_p)%*%S)))/n/(d+p*n-p))*t(chol(Psi_p))%*%rnorm(p)/sqrt(rchisq(1,d+p*n-p)/(d+p*n-p))
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)
    
    ### value of the proposal at new draw
    q1<-det(Psi_p)^(-0.5*(n+p+2))*(1+sum(diag(solve(Psi_p)%*%Cov_p))/d)^(-0.5*(p*n+d))
    
    ### value of the posterior
    mat_num<-matrix(0,p^2,p^2)
    mat_numJ<-matrix(0,p,p)
    det_num<-1
    exp_num<-1
    t_num<-0
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    mat_numJ<-mat_numJ+iPsiUi
    det_num<-det_num*sqrt(det(iPsiUi))
    t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
    }
    p1<-sqrt(det(mat_numJ/n))*sqrt(det(tGp%*%(mat_num-as.vector(mat_numJ)%*%t(as.vector(mat_numJ))/(n*p+d))%*%Gp/n))*det_num*(1+t_num)^(-0.5*(p*n+d))
    
    ### MH ratio
    ratio_MH<-p1*q0/p0/q1
    
    trans_p<-0
    if (runif(1) <= ratio_MH)
    {trans_p<-1
    mu0<-mu_p
    Psi0<-Psi_p
    q0<-q1
    p0<-p1
    }
    
    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu0)
    trans_m<-cbind(trans_m,trans_p)  
  }
  output<-list(mu_m,Psi_m,trans_m)  
}
#########################################################################################################
### algorithmA, reference normal, X:p\times n data matrix, U: pn \times pn matrix with uncertainties ####
#########################################################################################################
sample_post_nor_ref_marg_mu<-function(X,U,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size
  
  ############## addtional definitons  
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-Dp(p)
  tGp<-t(Gp)
  
  trans_m<-NULL
  mu_m<-NULL
  Psi_m<-NULL
  
  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  tcholS<-t(chol(S))
  
  ### generating an initial value for mu and Psi new draw from proposal
  mu0<-bar_X+sqrt((n-1)/n/(n-p))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p)/(n-p))
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)
  cCov_p<-chol(Cov_p)
  Z<-matrix(rnorm(p*n),p,n)
  Psi0<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p
  
  ### initial value of the proposal
  q0<-det(Psi0)^(-0.5*(n+p+1))*exp(-0.5*sum(diag(solve(Psi0)%*%Cov_p)))
  
  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  mat_numJ<-matrix(0,p,p)
  det_num<-1
  exp_num<-1
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  det_num<-det_num*sqrt(det(iPsiUi))
  exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
  }
  p0<-sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num*exp_num
  
  for (j in 1:Np)
  {
    ### generating a new draw Psi from the proposal 
    mu_p<-bar_X+sqrt((n-1)/n/(n-p))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p)/(n-p))
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)
    cCov_p<-chol(Cov_p)
    Z<-matrix(rnorm(p*n),p,n)
    Psi_p<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p
    
    ### value of the proposal at new draw
    q1 <-det(Psi_p)^(-0.5*(n+p+1))*exp(-0.5*sum(diag(solve(Psi_p)%*%Cov_p)))
    
    ### value of the posterior
    mat_num<-matrix(0,p^2,p^2)
    mat_numJ<-matrix(0,p,p)
    det_num<-1
    exp_num<-1
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    det_num<-det_num*sqrt(det(iPsiUi))
    exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
    }
    p1<-sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num*exp_num
    
    ### MH ratio
    ratio_MH<-p1*q0/p0/q1
    
    trans_p<-0
    if (runif(1) <= ratio_MH)
    {trans_p<-1
    mu0<-mu_p
    Psi0<-Psi_p
    q0<-q1
    p0<-p1
    }
    
    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu0)
    trans_m<-cbind(trans_m,trans_p)  
  }
  output<-list(mu_m,Psi_m,trans_m)  
}


#########################################################################################################
### algorithmB, reference normal, X:p\times n data matrix, U: pn \times pn matrix with uncertainties ####
#########################################################################################################
sample_post_nor_ref_marg_Psi<-function(X,U,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size
  
  ############## addtional definitons  
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-Dp(p)
  tGp<-t(Gp)
  
  trans_m<-NULL
  mu_m<-NULL
  Psi_m<-NULL
  
  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  cS<-sqrt(n-1)*chol(S)
  
  ### generating an initial value for mu and Psi new draw from proposal
  Z<-matrix(rnorm(p*n-p),p,n-1)
  Psi0<-t(cS)%*%solve(Z%*%t(Z))%*%cS
  mu0<-bar_X+t(chol(Psi0))%*%rnorm(p)/sqrt(n)
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)
  
  ### initial value of the proposal
  exp_prop<-exp(-0.5*sum(diag(solve(Psi0)%*%Cov_p)))
  q0<-det(Psi0)^(-0.5*(n+p+1))*exp_prop
  
  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  mat_numJ<-matrix(0,p,p)
  det_num<-1
  exp_num<-1
  #exp_prop<-1
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  det_num<-det_num*sqrt(det(iPsiUi))
  exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
  }
  p0<-sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num*exp_num
  
  for (j in 1:Np)
  {
    ### generating a new draw Psi from the proposal 
    Z<-matrix(rnorm(p*n-p),p,n-1)
    Psi_p<-t(cS)%*%solve(Z%*%t(Z))%*%cS
    mu_p<-bar_X+t(chol(Psi_p))%*%rnorm(p)/sqrt(n)
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)
    
    ### value of the proposal at new draw
    exp_prop<-exp(-0.5*sum(diag(solve(Psi_p)%*%Cov_p)))
    q1<-det(Psi_p)^(-0.5*(n+p+1))*exp_prop
    
    ### value of the posterior
    mat_num<-matrix(0,p^2,p^2)
    mat_numJ<-matrix(0,p,p)
    det_num<-1
    exp_num<-1
    #exp_prop<-1
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    det_num<-det_num*sqrt(det(iPsiUi))
    exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
    }
    p1<-sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num*exp_num
    
    ### MH ratio
    ratio_MH<-p1*q0/p0/q1
    
    trans_p<-0
    if (runif(1) <= ratio_MH)
    {trans_p<-1
    mu0<-mu_p
    Psi0<-Psi_p
    q0<-q1
    p0<-p1
    }
    
    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu0)
    trans_m<-cbind(trans_m,trans_p)  
  }
  output<-list(mu_m,Psi_m,trans_m)  
}

#########################################################################################################
### algorithmA, reference t, X:p\times n data matrix, U: pn \times pn matrix with uncertainties #########
#########################################################################################################
sample_post_t_ref_marg_mu<-function(X,U,d,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size
  
  ############## addtional definitons  
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-Dp(p)
  tGp<-t(Gp)
  
  trans_m<-NULL
  mu_m<-NULL
  Psi_m<-NULL
  
  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  tcholS<-t(chol(S))
  
  ### generating an initial value for mu and Psi new draw from proposal
  mu0<-bar_X+sqrt((n-1)/n/(n-p))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p)/(n-p))
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)
  cCov_p<-chol(Cov_p)
  Z<-matrix(rnorm(p*n),p,n)
  Psi0<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p*sqrt(rchisq(1,d)/(d))
  
  ### initial value of the proposal
  q0<-det(Psi0)^(-0.5*(n+p+1))*(1+sum(diag(solve(Psi0)%*%Cov_p))/d)^(-0.5*(p*n+d))
  
  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  mat_numJ<-matrix(0,p,p)
  det_num<-1
  exp_num<-1
  t_num<-0
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  mat_numJ<-mat_numJ+iPsiUi
  det_num<-det_num*sqrt(det(iPsiUi))
  t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
  }
  p0<-sqrt(det(tGp%*%(mat_num-as.vector(mat_numJ)%*%t(as.vector(mat_numJ))/(n*p+d))%*%Gp/n))*det_num*(1+t_num)^(-0.5*(p*n+d))
  
  for (j in 1:Np)
  {
    ### generating a new draw Psi from the proposal 
    mu_p<-bar_X+sqrt((n-1)/n/(n-p))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p)/(n-p))
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)
    cCov_p<-chol(Cov_p)
    Z<-matrix(rnorm(p*n),p,n)
    Psi_p<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p*sqrt(rchisq(1,d)/(d))
    
    ### value of the proposal at new draw
    q1 <-det(Psi_p)^(-0.5*(n+p+1))*(1+sum(diag(solve(Psi_p)%*%Cov_p))/d)^(-0.5*(p*n+d))
    
    ### value of the posterior
    mat_num<-matrix(0,p^2,p^2)
    mat_numJ<-matrix(0,p,p)
    det_num<-1
    exp_num<-1
    t_num<-0
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    mat_numJ<-mat_numJ+iPsiUi
    det_num<-det_num*sqrt(det(iPsiUi))
    t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
    }
    p1<-sqrt(det(tGp%*%(mat_num-as.vector(mat_numJ)%*%t(as.vector(mat_numJ))/(n*p+d))%*%Gp/n))*det_num*(1+t_num)^(-0.5*(p*n+d))
    
    ### MH ratio
    ratio_MH<-p1*q0/p0/q1
    
    trans_p<-0
    if (runif(1) <= ratio_MH)
    {trans_p<-1
    mu0<-mu_p
    Psi0<-Psi_p
    q0<-q1
    p0<-p1
    }
    
    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu0)
    trans_m<-cbind(trans_m,trans_p)  
  }
  output<-list(mu_m,Psi_m,trans_m)  
}

#########################################################################################################
### algorithmB, reference t, X:p\times n data matrix, U: pn \times pn matrix with uncertainties #########
#########################################################################################################
sample_post_t_ref_marg_Psi<-function(X,U,d,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size
  
  ############## addtional definitons  
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-Dp(p)
  tGp<-t(Gp)
  
  trans_m<-NULL
  mu_m<-NULL
  Psi_m<-NULL
  
  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  cS<-sqrt(n-1)*chol(S)
  
  ### generating an initial value for mu and Psi new draw from proposal
  Z<-matrix(rnorm(p*(n-1)),p,n-1)
  Psi0<-t(cS)%*%solve(Z%*%t(Z))%*%cS*sqrt(rchisq(1,d)/(d))
  mu0<-bar_X+sqrt((d+(n-1)*sum(diag(solve(Psi0)%*%S)))/n/(d+p*n-p))*t(chol(Psi0))%*%rnorm(p)/sqrt(rchisq(1,d+p*n-p)/(d+p*n-p))
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)
  
  ### initial value of the proposal
  q0<-det(Psi0)^(-0.5*(n+p+1))*(1+sum(diag(solve(Psi0)%*%Cov_p))/d)^(-0.5*(p*n+d))
  
  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  mat_numJ<-matrix(0,p,p)
  det_num<-1
  exp_num<-1
  t_num<-0
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  mat_numJ<-mat_numJ+iPsiUi
  det_num<-det_num*sqrt(det(iPsiUi))
  t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
  }
  p0<-sqrt(det(tGp%*%(mat_num-as.vector(mat_numJ)%*%t(as.vector(mat_numJ))/(n*p+d))%*%Gp/n))*det_num*(1+t_num)^(-0.5*(p*n+d))
  
  for (j in 1:Np)
  {
    ### generating a new draw Psi from the proposal 
    Z<-matrix(rnorm(p*(n-1)),p,n-1)
    Psi_p<-t(cS)%*%solve(Z%*%t(Z))%*%cS*sqrt(rchisq(1,d)/(d))
    mu_p<-bar_X+sqrt((d+(n-1)*sum(diag(solve(Psi_p)%*%S)))/n/(d+p*n-p))*t(chol(Psi_p))%*%rnorm(p)/sqrt(rchisq(1,d+p*n-p)/(d+p*n-p))
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)
    
    ### value of the proposal at new draw
    q1<-det(Psi_p)^(-0.5*(n+p+1))*(1+sum(diag(solve(Psi_p)%*%Cov_p))/d)^(-0.5*(p*n+d))
    
    ### value of the posterior
    mat_num<-matrix(0,p^2,p^2)
    mat_numJ<-matrix(0,p,p)
    det_num<-1
    exp_num<-1
    t_num<-0
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    mat_numJ<-mat_numJ+iPsiUi
    det_num<-det_num*sqrt(det(iPsiUi))
    t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
    }
    p1<-sqrt(det(tGp%*%(mat_num-as.vector(mat_numJ)%*%t(as.vector(mat_numJ))/(n*p+d))%*%Gp/n))*det_num*(1+t_num)^(-0.5*(p*n+d))
    
    ### MH ratio
    ratio_MH<-p1*q0/p0/q1
    
    trans_p<-0
    if (runif(1) <= ratio_MH)
    {trans_p<-1
    mu0<-mu_p
    Psi0<-Psi_p
    q0<-q1
    p0<-p1
    }
    
    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu0)
    trans_m<-cbind(trans_m,trans_p)  
  }
  output<-list(mu_m,Psi_m,trans_m)  
}


#########################################################################################################
### algorithmC, Jeffreys normal, X:p\times n data matrix, U: pn \times pn matrix with uncertainties ####
#########################################################################################################
sample_post_nor_jef_Gibbs<-function(X,U,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size
  
  ############## addtional definitons  
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-Dp(p)
  tGp<-t(Gp)
  
  trans_m<-NULL
  mu_m<-NULL
  Psi_m<-NULL
  
  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  tcholS<-t(chol(S))
  
  ### generating an initial value for mu and Psi new draw from proposal
  mu0<-bar_X+sqrt((n-1)/n/(n-p))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p)/(n-p))
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)
  cCov_p<-chol(Cov_p)
  Z<-matrix(rnorm(p*n+p),p,n+1)
  Psi0<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p
  
  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  det_num<-1
  bar_PsiX<-matrix(0,p,1)
  V_PsiX<-matrix(0,p,p)
  
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  bar_PsiX<-bar_PsiX+iPsiUi%*%X[,i_n]
  V_PsiX<-V_PsiX+iPsiUi
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  det_num<-det_num*sqrt(det(iPsiUi))
  }
  p0Psi<-sqrt(det(tGp%*%mat_num%*%Gp/n))*sqrt(det(V_PsiX))*det_num
  iV_PsiX0<-solve(V_PsiX)
  tx_PsiX0<-iV_PsiX0%*%bar_PsiX
  
  for (j in 1:Np)
  {
    ### generating a new draw mu and Psi from the proposal 
    mu_p<-tx_PsiX0+t(chol(iV_PsiX0))%*%rnorm(p)
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)
    cCov_p<-chol(Cov_p)
    Z<-matrix(rnorm(p*n+p),p,n+1)
    Psi_p<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p
    
    ### value of the proposal at new draw
    q0 <-det(Psi0)^(-0.5*(n+p+2))*exp(-0.5*sum(diag(solve(Psi0)%*%Cov_p)))
    q1 <-det(Psi_p)^(-0.5*(n+p+2))*exp(-0.5*sum(diag(solve(Psi_p)%*%Cov_p)))
    
    ### value of the posterior
    exp_num0<-1
    mat_num<-matrix(0,p^2,p^2)
    det_num<-1
    exp_num<-1
    bar_PsiX<-matrix(0,p,1)
    V_PsiX<-matrix(0,p,p)
    
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    bar_PsiX<-bar_PsiX+iPsiUi%*%X[,i_n]
    V_PsiX<-V_PsiX+iPsiUi
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    det_num<-det_num*sqrt(det(iPsiUi))
    exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
    iPsiUi_0<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    exp_num0<-exp_num0*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi_0%*%Xn[,i_n]))
    }
    
    p1Psi<-sqrt(det(tGp%*%mat_num%*%Gp/n))*sqrt(det(V_PsiX))*det_num
    p0<-p0Psi*exp_num0
    p1<-p1Psi*exp_num
    iV_PsiX<-solve(V_PsiX)
    tx_PsiX<-iV_PsiX%*%bar_PsiX
    
    ### MH ratio
    ratio_MH<-p1*q0/p0/q1
    
    trans_p<-0
    if (runif(1) <= ratio_MH)
    {trans_p<-1
    Psi0<-Psi_p
    iV_PsiX0<-iV_PsiX
    tx_PsiX0<-tx_PsiX
    p0Psi<-p1Psi
    }
    
    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu_p)
    trans_m<-cbind(trans_m,trans_p)  
  }
  output<-list(mu_m,Psi_m,trans_m)  
}

#########################################################################################################
### algorithmC, reference normal, X:p\times n data matrix, U: pn \times pn matrix with uncertainties ####
#########################################################################################################
sample_post_nor_ref_Gibbs<-function(X,U,Np){
    p<-nrow(X)  # model dimension
    n<-ncol(X)  # sample size
    
    ############## addtional definitons  
    bi_n<-rep(1,n)
    tbi_n<-t(bi_n)
    In<-diag(bi_n)
    Jn<-matrix(1,n,n)
    Ip<-diag(rep(1,p))
    Gp<-Dp(p)
    tGp<-t(Gp)
    
    trans_m<-NULL
    mu_m<-NULL
    Psi_m<-NULL
    
    bar_X<-X%*%bi_n/n
    S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
    tcholS<-t(chol(S))
    
    ### generating an initial value for mu and Psi new draw from proposal
    mu0<-bar_X+sqrt((n-1)/n/(n-p))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p)/(n-p))
    Xn<-X-mu0%*%tbi_n
    Cov_p<-Xn%*%t(Xn)
    cCov_p<-chol(Cov_p)
    Z<-matrix(rnorm(p*n),p,n)
    Psi0<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p
    
    ### initial value of the proposal
    q0<-det(Psi0)^(-0.5*(n+p+1))*exp(-0.5*sum(diag(solve(Psi0)%*%Cov_p)))
    
    ### initial value of the target (posterior)
    mat_num<-matrix(0,p^2,p^2)
    det_num<-1
    bar_PsiX<-matrix(0,p,1)
    V_PsiX<-matrix(0,p,p)
    
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    bar_PsiX<-bar_PsiX+iPsiUi%*%X[,i_n]
    V_PsiX<-V_PsiX+iPsiUi
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    det_num<-det_num*sqrt(det(iPsiUi))
    }
    p0Psi<-sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num
    iV_PsiX0<-solve(V_PsiX)
    tx_PsiX0<-iV_PsiX0%*%bar_PsiX
    
    for (j in 1:Np)
    {
      ### generating a new draw mu and Psi from the proposal 
      mu_p<-tx_PsiX0+t(chol(iV_PsiX0))%*%rnorm(p)
      Xn<-X-mu_p%*%tbi_n
      Cov_p<-Xn%*%t(Xn)
      cCov_p<-chol(Cov_p)
      Z<-matrix(rnorm(p*n),p,n)
      Psi_p<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p
      
      ### value of the proposal at new draw
      q0 <-det(Psi0)^(-0.5*(n+p+1))*exp(-0.5*sum(diag(solve(Psi0)%*%Cov_p)))
      q1 <-det(Psi_p)^(-0.5*(n+p+1))*exp(-0.5*sum(diag(solve(Psi_p)%*%Cov_p)))
      
      ### value of the posterior
      exp_num0<-1
      mat_num<-matrix(0,p^2,p^2)
      det_num<-1
      exp_num<-1
      bar_PsiX<-matrix(0,p,1)
      V_PsiX<-matrix(0,p,p)
      
      for (i_n in 1:n)
      {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
      bar_PsiX<-bar_PsiX+iPsiUi%*%X[,i_n]
      V_PsiX<-V_PsiX+iPsiUi
      mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
      det_num<-det_num*sqrt(det(iPsiUi))
      exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
      iPsiUi_0<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
      exp_num0<-exp_num0*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi_0%*%Xn[,i_n]))
      }
      p1Psi<-sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num
      p0<-p0Psi*exp_num0
      p1<-p1Psi*exp_num
      iV_PsiX<-solve(V_PsiX)
      tx_PsiX<-iV_PsiX%*%bar_PsiX
      
      ### MH ratio
      ratio_MH<-p1*q0/p0/q1
      
      trans_p<-0
      if (runif(1) <= ratio_MH)
      {trans_p<-1
      Psi0<-Psi_p
      iV_PsiX0<-iV_PsiX
      tx_PsiX0<-tx_PsiX
      p0Psi<-p1Psi
      }
      
      Psi_m<-cbind(Psi_m,as.vector(Psi0))
      mu_m<-cbind(mu_m,mu_p)
      trans_m<-cbind(trans_m,trans_p)  
    }
    output<-list(mu_m,Psi_m,trans_m)  
  }

#########################################################################################################
### algorithmC, Jeffreys t, X:p\times n data matrix, U: pn \times pn matrix with uncertainties ##########
#########################################################################################################
sample_post_t_jef_Gibbs<-function(X,U,d,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size
  
  ############## addtional definitons  
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-Dp(p)
  tGp<-t(Gp)
  
  trans_m<-NULL
  mu_m<-NULL
  Psi_m<-NULL
  
  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  tcholS<-t(chol(S))
  
  ### generating an initial value for mu and Psi new draw from proposal
  mu0<-bar_X+sqrt((n-1)/n/(n-p+1))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p+1)/(n-p+1))
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)
  cCov_p<-chol(Cov_p)
  Z<-matrix(rnorm(p*n+p),p,n+1)
  Psi0<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p*sqrt(rchisq(1,d)/(d))
  
  ### initial value of the proposal
  q0<-det(Psi0)^(-0.5*(n+p+2))*(1+sum(diag(solve(Psi0)%*%Cov_p))/d)^(-0.5*(p*n+d))
  
  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  det_num<-1
  bar_PsiX<-matrix(0,p,1)
  V_PsiX<-matrix(0,p,p)
  
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  bar_PsiX<-bar_PsiX+iPsiUi%*%X[,i_n]
  V_PsiX<-V_PsiX+iPsiUi
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  det_num<-det_num*sqrt(det(iPsiUi))
  }
  p0Psi<-sqrt(det(V_PsiX))*sqrt(det(tGp%*%(mat_num-as.vector(V_PsiX)%*%t(as.vector(V_PsiX))/(n*p+d))%*%Gp/n))*det_num
  iV_PsiX0<-solve(V_PsiX)
  tx_PsiX0<-iV_PsiX0%*%bar_PsiX
  
  for (j in 1:Np)
  {
    ### generating a new draw mu and Psi from the proposal 
    t_mu<-0
    Xnt<-X-tx_PsiX0%*%tbi_n
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    t_mu<-t_mu+sum(t(Xnt[,i_n])%*%iPsiUi%*%Xnt[,i_n])
    }
    
    mu_p<-tx_PsiX0+sqrt((d+t_mu)/(p*n+d-p))*t(chol(iV_PsiX0))%*%rnorm(p)/sqrt(rchisq(1,n*p+d-p)/(n*p+d-p))
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)
    cCov_p<-chol(Cov_p)
    Z<-matrix(rnorm(p*n+p),p,n+1)
    Psi_p<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p*sqrt(rchisq(1,d)/(d))
    
    ### value of the proposal at new draw
    q0 <-det(Psi0)^(-0.5*(n+p+2))*(1+sum(diag(solve(Psi0)%*%Cov_p))/d)^(-0.5*(p*n+d))
    q1 <-det(Psi_p)^(-0.5*(n+p+2))*(1+sum(diag(solve(Psi_p)%*%Cov_p))/d)^(-0.5*(p*n+d))
    
    ### value of the posterior
    t_num0<-0
    mat_num<-matrix(0,p^2,p^2)
    det_num<-1
    t_num<-0
    bar_PsiX<-matrix(0,p,1)
    V_PsiX<-matrix(0,p,p)
    
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    bar_PsiX<-bar_PsiX+iPsiUi%*%X[,i_n]
    V_PsiX<-V_PsiX+iPsiUi
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    det_num<-det_num*sqrt(det(iPsiUi))
    t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
    iPsiUi_0<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    t_num0<-t_num0+sum(t(Xn[,i_n])%*%iPsiUi_0%*%Xn[,i_n])/d
    }
    p1Psi<-sqrt(det(V_PsiX))*sqrt(det(tGp%*%(mat_num-as.vector(V_PsiX)%*%t(as.vector(V_PsiX))/(n*p+d))%*%Gp/n))*det_num
    p0<-p0Psi*(1+t_num0)^(-0.5*(p*n+d))
    p1<-p1Psi*(1+t_num)^(-0.5*(p*n+d))
    iV_PsiX<-solve(V_PsiX)
    tx_PsiX<-iV_PsiX%*%bar_PsiX
    
    ### MH ratio
    ratio_MH<-p1*q0/p0/q1
    
    trans_p<-0
    if (runif(1) <= ratio_MH)
    {trans_p<-1
    #mu0<-mu_p
    Psi0<-Psi_p
    iV_PsiX0<-iV_PsiX
    tx_PsiX0<-tx_PsiX
    p0Psi<-p1Psi
    }
    
    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu_p)
    trans_m<-cbind(trans_m,trans_p)  
  }
  output<-list(mu_m,Psi_m,trans_m)  
}



#########################################################################################################
### algorithmC, reference t, X:p\times n data matrix, U: pn \times pn matrix with uncertainties #########
#########################################################################################################
sample_post_t_ref_Gibbs<-function(X,U,d,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size
  
  ############## addtional definitons  
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-Dp(p)
  tGp<-t(Gp)
  
  trans_m<-NULL
  mu_m<-NULL
  Psi_m<-NULL
  
  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  tcholS<-t(chol(S))
  
  ### generating an initial value for mu and Psi new draw from proposal
  mu0<-bar_X+sqrt((n-1)/n/(n-p))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p)/(n-p))
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)
  cCov_p<-chol(Cov_p)
  Z<-matrix(rnorm(p*n),p,n)
  Psi0<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p*sqrt(rchisq(1,d)/(d))
  
  ### initial value of the proposal
  q0<-det(Psi0)^(-0.5*(n+p+1))*(1+sum(diag(solve(Psi0)%*%Cov_p))/d)^(-0.5*(p*n+d))
  
  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  det_num<-1
  bar_PsiX<-matrix(0,p,1)
  V_PsiX<-matrix(0,p,p)
  
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  bar_PsiX<-bar_PsiX+iPsiUi%*%X[,i_n]
  V_PsiX<-V_PsiX+iPsiUi
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  det_num<-det_num*sqrt(det(iPsiUi))
  }
  p0Psi<-sqrt(det(tGp%*%(mat_num-as.vector(V_PsiX)%*%t(as.vector(V_PsiX))/(n*p+d))%*%Gp/n))*det_num
  iV_PsiX0<-solve(V_PsiX)
  tx_PsiX0<-iV_PsiX0%*%bar_PsiX
  
  for (j in 1:Np)
  {
    ### generating a new draw mu and Psi from the proposal 
    t_mu<-0
    Xnt<-X-tx_PsiX0%*%tbi_n
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    t_mu<-t_mu+sum(t(Xnt[,i_n])%*%iPsiUi%*%Xnt[,i_n])
    }
    
    mu_p<-tx_PsiX0+sqrt((d+t_mu)/(p*n+d-p))*t(chol(iV_PsiX0))%*%rnorm(p)/sqrt(rchisq(1,n*p+d-p)/(n*p+d-p))
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)
    cCov_p<-chol(Cov_p)
    Z<-matrix(rnorm(p*n),p,n)
    Psi_p<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p*sqrt(rchisq(1,d)/(d))
    
    ### value of the proposal at new draw
    q0 <-det(Psi0)^(-0.5*(n+p+1))*(1+sum(diag(solve(Psi0)%*%Cov_p))/d)^(-0.5*(p*n+d))
    q1 <-det(Psi_p)^(-0.5*(n+p+1))*(1+sum(diag(solve(Psi_p)%*%Cov_p))/d)^(-0.5*(p*n+d))
    
    ### value of the posterior
    t_num0<-0
    mat_num<-matrix(0,p^2,p^2)
    det_num<-1
    t_num<-0
    bar_PsiX<-matrix(0,p,1)
    V_PsiX<-matrix(0,p,p)
    
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    bar_PsiX<-bar_PsiX+iPsiUi%*%X[,i_n]
    V_PsiX<-V_PsiX+iPsiUi
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    det_num<-det_num*sqrt(det(iPsiUi))
    t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
    iPsiUi_0<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    t_num0<-t_num0+sum(t(Xn[,i_n])%*%iPsiUi_0%*%Xn[,i_n])/d
    }
    p1Psi<-sqrt(det(tGp%*%(mat_num-as.vector(V_PsiX)%*%t(as.vector(V_PsiX))/(n*p+d))%*%Gp/n))*det_num
    p0<-p0Psi*(1+t_num0)^(-0.5*(p*n+d))
    p1<-p1Psi*(1+t_num)^(-0.5*(p*n+d))
    iV_PsiX<-solve(V_PsiX)
    tx_PsiX<-iV_PsiX%*%bar_PsiX
    
    ### MH ratio
    ratio_MH<-p1*q0/p0/q1
    
    trans_p<-0
    if (runif(1) <= ratio_MH)
    {trans_p<-1
    #mu0<-mu_p
    Psi0<-Psi_p
    iV_PsiX0<-iV_PsiX
    tx_PsiX0<-tx_PsiX
    p0Psi<-p1Psi
    }
    
    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu_p)
    trans_m<-cbind(trans_m,trans_p)  
  }
  output<-list(mu_m,Psi_m,trans_m)  
}


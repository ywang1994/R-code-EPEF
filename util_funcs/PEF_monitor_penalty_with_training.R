

thehat_training <- function(x,initval,M2, M1,p,q,momente,lambda,a,tau=0,tw = NA) {
  #x <- sim$y
  #initval <- est
  #tau=0
  #c <- 2.632
  n <- length(x)
  pdim <- 1+p+q   # number of parameters to be estimated
  #estimate <- matrix(0, ncol = pdim, nrow = iteration)
  err <- rep(1,pdim)
  mue <- momente[1]
  vare <- momente[2]
  skewe <- momente[3]
  kurte <- momente[4]
  
  ########################################################################    
  ## INITIALIZATION OF ARRAYS
  ########################################################################
  ## Identity matrix
  iden = diag(pdim)
  
  ## psi hat
  psih <- rep(1,n)
  
  ## k matrix (variance-covariance) and k inverse (observed information)
  kmat = array(NA, dim = c(pdim, pdim, n)) 
  kinv = array(NA, dim = c(pdim, pdim, n))
  kmat2 = array(NA, dim = c(pdim, pdim, n)) 
  kinv2 = array(NA, dim = c(pdim, pdim, n))
  termk1 = array(NA, dim=c(pdim,pdim,n))
  termk2 = array(NA, dim=c(pdim,pdim,n))
  
  ## Parameter estimates for each iteration
  thehat = array(NA, dim = c(pdim, 1, n))
  
  ## Derivative of psi and second derivative of psi
  derpsi<-matrix(rep(0),pdim,1)
  der2psi<-matrix(rep(0),pdim,pdim)
  
  ## Derivative of mu, sigsq, gamma, kappa; second derivates of mu, sigsq				
  dermu<-matrix(rep(0),pdim,1)				
  dersigsq<-matrix(rep(0),pdim,1)
  dergamma<-matrix(rep(0),pdim,1)
  derkappa<-matrix(rep(0),pdim,1)				
  der2mu<-matrix(rep(0),pdim,pdim)				
  der2sigsq<-matrix(rep(0),pdim,pdim)
  
  ## Derivative of m(t) and M(t) 
  derm<-matrix(rep(0),pdim,1)
  derqm<-matrix(rep(0),pdim,1)
  
  ## Derivative of quadratic variation, quadratic covariation, eta, rho
  dervm<-matrix(rep(0),pdim,1)
  dervqm<-matrix(rep(0),pdim,1)
  dervmqm<-matrix(rep(0),pdim,1)
  dereta<-matrix(rep(0),pdim,1)
  derrho<-matrix(rep(0),pdim,1)
  
  
  ## Optimal astr and bstr
  astr<-matrix(rep(0),pdim,1)
  bstr<-matrix(rep(0),pdim,1)
  
  ## Derivative of astr and bstr
  derastr<-matrix(rep(0),pdim,pdim)
  derbstr<-matrix(rep(0),pdim,pdim)
  
  
  ########################################################################
  ## INITIAL VALUES: for omega, alpha, beta, passed from main program
  ########################################################################
  
  #initial<- as.numeric(c(initval[1], initval[2], initval[3]))
  initial <- initval
  
  ########################################################################
  ## RECURSIVE FORMULAS FROM ESTIMATING EQUATIONS
  ########################################################################
  ## Put initial values into initial positions of arrays
  
  thehat[,,1] <- initial
  init.mat <- InitMat(x[1:M2],p,q)
  gg <- matrix(0, nrow=pdim, ncol=pdim)
  gg_training <- c() ## g functions matrix for bootstrap
  monitor <- matrix(0, nrow=pdim, ncol= n)
  CK <-  matrix(0, nrow=pdim, ncol= n) ## Only one occurrence in the code
  bound <- rep(0, n)
  
  
  sum_termk <- 0 ## Only one occurrence in the code
  mxpq <- max(p,q)
  for (t in 1:mxpq){
    psih[t] = thehat[1,1,1]    # omega 
    thehat[,1,t] <- thehat[,1,1]
    kinv[,,t] <- init.mat     ###### KEY PART #############
    #kinv2[,,t] <- init.mat
    kmat[,,t] <- solve(kinv[,,t])
    #kmat2[,,t] <- solve(kinv2[,,t])
  }
  
  pp <- NULL
  
  ## t=(mxpq+1):n					
  for (t in (mxpq+1):M2)
  {	
    ## Make sure the theta only updated to M2
    if(t>0) {
      omega <- as.matrix(thehat[1,1,t-1],ncol=1,nrow=1)
      alpha <- as.matrix(thehat[2:(1+p),1,t-1],ncol=1,nrow=p)
      if(q != 0){
        beta <- as.matrix(thehat[(p+2):pdim,1,t-1],ncol=1,nrow=q)
      }
      # print(omega[1,1])
      # print(alpha[,1])
    }
    else{
      omega <- as.matrix(thehat[1,1,M2],ncol=1,nrow=1)
      alpha <- as.matrix(thehat[2:(1+p),1,M2],ncol=1,nrow=p)
      if(q != 0){
        beta <- as.matrix(thehat[(p+2):pdim,1,M2],ncol=1,nrow=q)
      }
    }
    
    pos1 <- t-1
    pos2 <- t-p
    pos3 <- t-q
    
    psih[t]=psit(x,psih,omega, alpha, beta, pos1, pos2, pos3, q)
    if(abs(psih[t]) > log(.Machine$double.xmax/max(momente))/4) {
      psih[t] <- sign(psih[t])*log(.Machine$double.xmax/max(momente))/8
    }
    
    ## Define derivatives of psih(t) wrt theta: pdim*1 vector							
    ## First derivative									
    derpsi=derpsit(x, psih, pos1, pos2, pos3, pdim, q)
    
    
    
    ## Second derivative: pdim*pdim matrix									
    der2psi=der2psit(pdim)   
    
    ## mu(t), sigsq(t), gamma(t), kappa(t)							
    mu=exp(psih[t])*mue 							
    sigsq=vare*exp(2*psih[t])			
    gamma=skewe*exp(3*psih[t])	 # recall this is third central moment							
    kappa=kurte*exp(4*psih[t])	 # recall this is fourth central moment
    
    ## First Derivatives of mu(t),sigsq(t), gamma(t) and kappa(t)												
    dermu=exp(psih[t])*derpsi*mue							
    dersigsq=2*vare*exp(2*psih[t])*derpsi
    dergamma=3*skewe*exp(3*psih[t])*derpsi	
    derkappa=4*kurte*exp(4*psih[t])*derpsi
    
    ## Second Derivatives of mu(t) and sigsq(t)							
    der2mu <- mu*((derpsi)%*%t(derpsi) + der2psi)
    der2sigsq <- 2*sigsq*(der2psi + 2*(derpsi)%*%t(derpsi))
    
    ## Compute m(t) and M(t)							
    m = x[t]-mu       							
    qm = m**2-sigsq
    
    ## Compute Quadratic variations of m(t) and M(t) and 
    ## covariance of (m(t), M(t))							
    vm = sigsq
    vqm = kappa-vm**2  
    vmqm = gamma
    
    ## Define rho^2(t) 															
    termr = 1-((vmqm**2)/(vm*vqm))								
    rho = 1/termr
    if(is.na(rho) | is.infinite(rho))
      rho <- 1
    
    ## Define eta 								
    eta = vmqm/(vm*vqm)	
    if(is.na(eta) | is.infinite((eta)))
      eta = 0
    
    ## Define vectors astr and bstr						
    astr = rho*(-dermu/vm + dersigsq*eta)							
    bstr = rho*(dermu*eta - dersigsq/vqm)
    
    ## Define Derivatives of m(t) and M(t) 							
    derm = -dermu
    derqm = 2*m*derm - dersigsq	
    
    ## Derivatives of variations <m>(t), <M>(t) and <m,M>(t)							
    dervm = dersigsq						
    dervqm = derkappa - 2*sigsq*dersigsq
    dervmqm = dergamma
    
    #Note: derrho=0 - and it is!
    ru=vm*vqm
    rv=vm*vqm-vmqm**2
    rdu=vm*dervqm+vqm*dervm
    rdv=vm*dervqm+vqm*dervm -2*vmqm*dervmqm
    derrho=(rv*rdu-ru*rdv)/(rv**2)
    if(sum(is.na(derrho)) >0)
      derrho <- matrix(rep(0),pdim,1)
    ## Derivative of eta(t)							
    num1=vm*vqm*dervmqm
    num2=vmqm*(vm*dervqm + vqm*dervm)
    den=(vm**2)*(vqm**2)
    dereta=(num1-num2)/den
    if(sum(is.infinite(dereta)) >0 | sum(is.na(dereta))>0)
      dereta <- matrix(rep(0),pdim,1)
    ## Derivatives of astr and bstr							
    ## For astr						
    terma1 = der2mu/vm - dermu%*%t(dervm)/(vm**2)
    terma2 = der2sigsq*eta + dersigsq%*%t(dereta)
    terma3 = -dermu/vm + dersigsq*eta							
    derastr = -rho*terma1 + rho*terma2 + derrho%*%t(terma3)
    #derastr = -rho*terma1 + rho*terma2 
    
    ## For bstr						
    termb1 = der2mu*eta + dermu%*%t(dereta)													
    termb21 = der2sigsq/vqm 
    termb22 = (dersigsq%*%t(dervqm))/(vqm**2)
    if(sum(is.infinite(termb22)>0 | sum(is.na(termb22))>0))
      termb22 <- matrix(0, ncol=pdim, nrow = pdim)
    termb2 = termb21 - termb22
    termb3 = dermu*eta - dersigsq/vqm
    derbstr = rho*termb1 - rho*termb2 + derrho%*%t(termb3)
    
    
    ### Recursive Formulas						
    ## Compute Kinv(t): Information 							
    termk1[,,t] = astr%*%t(derm) + m * derastr 							
    termk2[,,t] = bstr%*%t(derqm) + qm * derbstr 
    
    ## add derivative of penalty here
    if(t>0){
      penalty <- p.lam.prime(thehat[,,t-1], thehat[,,1] ,lambda, a, tau)
    }
    else{
      penalty <- p.lam.prime(thehat[,,M2], thehat[,,1] ,lambda, a, tau)
    }
    
    
    pp <- rbind(pp, penalty)
    #penalty <- p.lam.prime(thehat[,,t-1],lambda, a)
    penalty <- matrix(penalty, ncol = 1, nrow = pdim)
    
    ## Need Numerical FIX when termk1 and termk2 are NA
    if ( sum(as.numeric(is.na(termk1[,,t]))) > 0){
      foundFlag1 = 0
      for (w in 1:(t-1)) {
        if ( (foundFlag1 == 0) & sum(as.numeric(is.na(termk1[,,t-w]))) == 0 ) {
          termk1[,,t] = termk1[,,(t-w)]
          foundFlag1 = 1
        }
        else
          termk1[,,t]= matrix(rep(.Machine$double.xmax/3,pdim*pdim),pdim,pdim)
      }
    }
    
    if ( sum(as.numeric(is.na(termk2[,,t]))) > 0){
      foundFlag2 = 0
      for (w in 1:(t-1)) {
        if ( (foundFlag2 == 0) & sum(as.numeric(is.na(termk2[,,t-w]))) == 0 ) {
          termk2[,,t] = termk2[,,(t-w)]
          foundFlag2 = 1
        }
        else
          termk2[,,t]= matrix(rep(.Machine$double.xmax/3,pdim*pdim),pdim,pdim)
      }
    }
    
    
    p_lambda <- p.lam.2prime(thehat[,,1], lambda, a, tau)
    p_lambda <- diag(p_lambda)
    if(t >= 0*n){
      kinv[,,t] = kinv[,,t-1] - (termk1[,,t]+termk2[,,t] - p_lambda)
      # kinv[,,t] = kinv[,,t-1] - (termk1[,,t]+termk2[,,t])
      
    }else{
      kinv[,,t] = kinv[,,t-1] - (termk1[,,t]+termk2[,,t])
      
    }
    
    
    s <- svd(kinv[,,t])
    
    inf.index <- which(is.infinite(s$d))
    D <- diag(1/s$d)
    for(i in 1:length(s$d)){
      if(is.infinite(D[i,i])) 
        D[i,i] <- .Machine$double.xmin*100
    }
    kmat[,,t]<- s$v%*%D%*%t(s$u)
    
    termt = astr*m + bstr*qm
    if(sum(as.numeric(is.na(termt))) >0 ) {
      break
    }
    thehat[,,t] = thehat[,,t-1] + kmat[,,t]%*%(termt - penalty)
    # thehat[,,t] = thehat[,,t-1] + kmat[,,t]%*%(termt)
    
  }
  return(thehat[,,M2])
}



cal_G_func_AR =function(dat,thehat,momente,M1,M2,
                        est,lambda,a,tau=0,p.lam.prime){
  #compute psih AR once and for all!!!
  
  ################ Compute psih
  n=length(dat)
  pdim = length(thehat)
  psih = rep(NA,n)
  psih[1:(pdim-1)] = thehat[1]
  dat_embed=embed(dat,dimension = pdim) ## embedded dat for psih matrix calc
  dat_embed[,1]=1
  psih[pdim:n]=c(dat_embed%*%thehat)
  ################
  
  
  ##################Create the n X 1 vectors for g func
  mue <- momente[1]
  vare <- momente[2]
  skewe <- momente[3]
  kurte <- momente[4]
  
  
  ### n X 1 vecs
  mu=exp(psih)*mue 							
  sigsq=vare*exp(2*psih)			
  gamma=skewe*exp(3*psih)	 # recall this is third central moment							
  kappa=kurte*exp(4*psih)	 # recall this is fourth central moment
  
  m = dat-mu       							
  qm = m**2-sigsq
  
  ## Compute Quadratic variations of m(t) and M(t) and 
  ## covariance of (m(t), M(t))							
  vm = sigsq
  vqm = kappa-vm**2  
  vmqm = gamma
  
  ## Define rho^2(t) (n X 1)															
  termr = 1-((vmqm**2)/(vm*vqm))								
  rho = 1/termr
  rho[is.na(rho)|is.infinite(rho)]=1
  
  ## Define eta (n X 1)								
  eta = vmqm/(vm*vqm)	
  eta[is.na(eta) | is.infinite((eta))]=0
  
  #### derivative of psi (n X pdim) vectors
  derpsi = rbind(matrix(rep(c(1,rep(0,pdim-1)),times=pdim-1),nrow=pdim-1,byrow = TRUE),dat_embed)
  
  dermu=matrix(rep(exp(psih),times=pdim),ncol=pdim)*derpsi*mue							
  dersigsq=2*vare*matrix(rep(exp(2*psih),times=pdim),ncol=pdim)*derpsi
  
  #### Define the astr*m + bstr*qm - penalty
  ### (n X pdim)
  astr = matrix(rep(rho,times=pdim),ncol=pdim)*(-dermu/
                                                  matrix(rep(vm,times=pdim),ncol=pdim)
                                                + dersigsq*
                                                  matrix(rep(eta,times=pdim),ncol=pdim))							
  bstr = matrix(rep(rho,times=pdim),ncol=pdim)*(dermu*
                                                  matrix(rep(eta,times=pdim),ncol=pdim)
                                                - dersigsq/
                                                  matrix(rep(vqm,times=pdim),ncol=pdim))
  ### (n X pdim)
  penalty <- matrix(rep(p.lam.prime(thehat, est,lambda, a, tau=0),times=n),ncol=pdim,byrow = TRUE)
  
  
  score = astr*matrix(rep(m,times=pdim),ncol=pdim)+
    bstr*matrix(rep(qm,times=pdim),ncol=pdim)-penalty
  
  Dm_diag = sqrt(diag(t(score[(M1+1):M2,])%*%(score[(M1+1):M2,])/(M2-M1))) ## Scaling factor on the diagnonal
  
  G_func = (apply(score[(M2+1):n,],2,cumsum))/
    matrix(rep(Dm_diag,times=(n-M2)),ncol=pdim,byrow = TRUE) ### Detector
  return(matrix(G_func,ncol=1))
}

cal_G_func_AR_corrected=function(dat,thehat,momente,M1,M2,
                                 est,lambda,a,tau=0,p.lam.prime){
  n=length(dat)
  pdim = length(thehat)
  psih = rep(NA,n)
  psih[1:(pdim-1)] = thehat[1]
  dat_embed=embed(dat,dimension = pdim) ## embedded dat for psih matrix calc
  dat_embed[,1]=1
  psih[pdim:n]=c(dat_embed%*%thehat)
  ################
  
  
  ##################Create the n X 1 vectors for g func
  mue <- momente[1]
  vare <- momente[2]
  skewe <- momente[3]
  kurte <- momente[4]
  
  
  ### n X 1 vecs
  mu=exp(psih)*mue 							
  sigsq=vare*exp(2*psih)			
  gamma=skewe*exp(3*psih)	 # recall this is third central moment							
  kappa=kurte*exp(4*psih)	 # recall this is fourth central moment
  
  m = dat-mu       							
  qm = m**2-sigsq
  
  ## Compute Quadratic variations of m(t) and M(t) and 
  ## covariance of (m(t), M(t))							
  vm = sigsq
  vqm = kappa-vm**2  
  vmqm = gamma
  
  ## Define rho^2(t) (n X 1)															
  termr = 1-((vmqm**2)/(vm*vqm))								
  rho = 1/termr
  rho[is.na(rho)|is.infinite(rho)]=1
  
  ## Define eta (n X 1)								
  eta = vmqm/(vm*vqm)	
  eta[is.na(eta) | is.infinite((eta))]=0
  
  #### derivative of psi (n X pdim) vectors
  derpsi = rbind(matrix(rep(c(1,rep(0,pdim-1)),times=pdim-1),nrow=pdim-1,byrow = TRUE),dat_embed)
  
  dermu=matrix(rep(exp(psih),times=pdim),ncol=pdim)*derpsi*mue							
  dersigsq=2*vare*matrix(rep(exp(2*psih),times=pdim),ncol=pdim)*derpsi
  
  #### Define the astr*m + bstr*qm - penalty
  ### (n X pdim)
  astr = matrix(rep(rho,times=pdim),ncol=pdim)*(-dermu/
                                                  matrix(rep(vm,times=pdim),ncol=pdim)
                                                + dersigsq*
                                                  matrix(rep(eta,times=pdim),ncol=pdim))							
  bstr = matrix(rep(rho,times=pdim),ncol=pdim)*(dermu*
                                                  matrix(rep(eta,times=pdim),ncol=pdim)
                                                - dersigsq/
                                                  matrix(rep(vqm,times=pdim),ncol=pdim))
  ### (n X pdim)
  penalty <- matrix(rep(p.lam.prime(thehat, est,lambda, a, tau=0),times=n),ncol=pdim,byrow = TRUE)
  
  
  score = astr*matrix(rep(m,times=pdim),ncol=pdim)+
    bstr*matrix(rep(qm,times=pdim),ncol=pdim)-penalty
  
  ### Calculating the Dm^(-1/2)
  Dm= t(score[(M1+1):M2,])%*%(score[(M1+1):M2,])/(M2-M1) ## no burn-in
  eigen_decomp =  eigen(Dm)
  eigen_decomp_val = eigen_decomp$values
  eigen_decomp_vecs = eigen_decomp$vectors
  Dm_neg_sqrt = eigen_decomp_vecs%*%diag(eigen_decomp_val**-.5)%*%t(eigen_decomp_vecs)
  
  G_func = apply(abs(apply(score[(M2+1):n,],2,cumsum)%*%Dm_neg_sqrt),1,max)
  return((G_func))
}

cal_g_func_AR =function(dat,thehat,momente,M1,M2,
                        est,lambda,a,tau=0){
  #compute psih AR once and for all!!!
  
  ################ Compute psih
  n=length(dat)
  pdim = length(thehat)
  psih = rep(NA,n)
  psih[1:(pdim-1)] = thehat[1]
  dat_embed=embed(dat,dimension = pdim) ## embedded dat for psih matrix calc
  dat_embed[,1]=1
  psih[pdim:n]=c(dat_embed%*%thehat)
  ################
  
  
  ##################Create the n X 1 vectors for g func
  mue <- momente[1]
  vare <- momente[2]
  skewe <- momente[3]
  kurte <- momente[4]
  
  
  ### n X 1 vecs
  mu=exp(psih)*mue 							
  sigsq=vare*exp(2*psih)			
  gamma=skewe*exp(3*psih)	 # recall this is third central moment							
  kappa=kurte*exp(4*psih)	 # recall this is fourth central moment
  
  m = dat-mu       							
  qm = m**2-sigsq
  
  ## Compute Quadratic variations of m(t) and M(t) and 
  ## covariance of (m(t), M(t))							
  vm = sigsq
  vqm = kappa-vm**2  
  vmqm = gamma
  
  ## Define rho^2(t) (n X 1)															
  termr = 1-((vmqm**2)/(vm*vqm))								
  rho = 1/termr
  rho[(rho>1)|is.na(rho)|is.infinite(rho)]=1
  
  ## Define eta (n X 1)								
  eta = vmqm/(vm*vqm)	
  eta[is.na(eta) | is.infinite((eta))]=0
  
  #### derivative of psi (n X pdim) vectors
  derpsi = rbind(matrix(rep(c(1,rep(0,pdim-1)),times=pdim-1),nrow=pdim-1,byrow = TRUE),dat_embed)
  
  dermu=matrix(rep(exp(psih),times=pdim),ncol=pdim)*derpsi*mue							
  dersigsq=2*vare*matrix(rep(exp(2*psih),times=pdim),ncol=pdim)*derpsi
  
  #### Define the astr*m + bstr*qm - penalty
  ### (n X pdim)
  astr = matrix(rep(rho,times=pdim),ncol=pdim)*(-dermu/
                                                  matrix(rep(vm,times=pdim),ncol=pdim)
                                                + dersigsq*
                                                  matrix(rep(eta,times=pdim),ncol=pdim))							
  bstr = matrix(rep(rho,times=pdim),ncol=pdim)*(dermu*
                                                  matrix(rep(eta,times=pdim),ncol=pdim)
                                                - dersigsq/
                                                  matrix(rep(vqm,times=pdim),ncol=pdim))
  ### (n X pdim)
  penalty <- matrix(rep(p.lam.prime(thehat, est,lambda, a, tau=0),times=n),ncol=pdim,byrow = TRUE)
  
  
  score = astr*matrix(rep(m,times=pdim),ncol=pdim)+
    bstr*matrix(rep(qm,times=pdim),ncol=pdim)-penalty
  
  
  g_func = score
  ### Detector
  return(c(g_func))
}

ts0_boot= function(
  boot_ind,
  csv_temp,
  len_ts_mon = 5000,
  mean_l = 20,
  r=2,## Ratio of resampling size
  block_resample_indexing
){
  ts0= csv_temp$V1 ## Original ts in training period
  # M2=length(ts0) ## Serving to keep only the monitoring period
  blocks=1+rgeom(r*length(ts0)/mean_l,p=1/mean_l)## create blocks of non-zero length
  block_id =1:(sum(cumsum(blocks)<length(ts0))+1) ## Vec
  blocks_to_sample = blocks[block_id]
  blocks_to_sample[length(blocks_to_sample)] = length(ts0)-sum(blocks[block_id[-length(block_id)]])
  starts = c(1,1+cumsum(blocks_to_sample)[-length(blocks_to_sample)])
  ends = c(cumsum(blocks_to_sample)[-length(blocks_to_sample)],length(ts0))
  blocks_to_sample_mat = cbind(blocks_to_sample,starts,ends)
  
  ## Begin bootstrapping of block id
  block_id_resampled_ALL = sample(length(blocks_to_sample),size = r*(len_ts_mon)/mean_l,replace = TRUE)
  block_id_resampled_needed_ind = 1:which(cumsum(blocks_to_sample_mat[block_id_resampled_ALL,1])>=(len_ts_mon))[1]
  block_id_resampled_needed =block_id_resampled_ALL[block_id_resampled_needed_ind]
  
  ts0_boot0 = Reduce(c,sapply(block_id_resampled_needed,block_resample_indexing,
                              blocks_to_sample_mat=blocks_to_sample_mat,ts0=ts0))
  return(ts0_boot0[1:len_ts_mon])
}


block_resample_indexing = function(block_id_resampled_needed_i,blocks_to_sample_mat,ts0){
  return(ts0[blocks_to_sample_mat[block_id_resampled_needed_i,2]:
               blocks_to_sample_mat[block_id_resampled_needed_i,3]])
} ## Help extract the resampled time series


ts0_boot_corrected= function(
  boot_ind,
  csv_temp,
  len_ts_mon = 5000,
  mean_l = 20,
  r=2## Ratio of resampling size
){
  ts0= csv_temp$V1 ## Original ts in training period
  ts0_long_vector_for_circles =rep(ts0,2)
  # M2=length(ts0) ## Serving to keep only the monitoring period
  blocks=1+rgeom(r*(len_ts_mon)/mean_l,p=1/mean_l)## create blocks of non-zero length
  block_id =1:(sum(cumsum(blocks)<len_ts_mon)+1) ## Vec
  blocks_to_sample = blocks[block_id]
  blocks_to_sample[length(blocks_to_sample)] = len_ts_mon-sum(blocks[block_id[-length(block_id)]])
  ts0_ind_start = sample(1:length(ts0),size=length(block_id),replace = TRUE)
  ts0_ind_end = (ts0_ind_start+blocks_to_sample-1)
  blocks_to_sample_mat = cbind(ts0_ind_start,ts0_ind_end)
  blocks_to_sample_mat_list=split(blocks_to_sample_mat,row(blocks_to_sample_mat))
  ts0_list=lapply(blocks_to_sample_mat_list,function(x) return(ts0_long_vector_for_circles[x[1]:x[2]]))
  
  return(Reduce(c,ts0_list))
}




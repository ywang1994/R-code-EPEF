# This is Yanzhao's user define functions to reduce the number of lines in main script
####################################################################################
# Step 1 fit the training data and select the optimal lambda
# I plan to use foreach function doparallel
lambda_tuning_indiv = function(lambda_temp,
                               duras_ts, # duration time series
                               model.order,
                               M1, # burning 
                               M2,# training
                               sparsity_num,sparsity_thres,a){
  dat <- duras_ts[1:M2]
  p <- model.order[1]
  q <- model.order[2]
  # p <- ps[1]
  # q <- qs[1]
  pdim <- p+q+1
  est <- finitval(dat,p,q)
  if (q==0) est <- est[-(pdim+1)]
  a=3.7
  psih <- psi_cal(dat,est,p,q)
  eps <- dat/exp(psih)
  ed <- eps_dist2(eps)
  momente <- c(ed$mue, ed$vare, ed$skewe, ed$kurte)
  difference <- log(momente[1]) - mean(log(eps))
  thehat = thehat_training(dat,est,M2, M1, p,q,momente,lambda=lambda_temp,a,tau=0)
  thehat.used_temp  = thehat
  thehat.used_temp[1] =thehat.used_temp[1] + difference
  MAD_temp=mean(abs(dat[-(1:p)]-cbind(1,embed(dat[-M2],dimension = p))%*%thehat.used_temp))
  Sparse_or_not = as.numeric(sum(abs(thehat.used_temp[-1])>=sparsity_thres)<=sparsity_num)
  return(c(lambda_temp,MAD_temp,Sparse_or_not))
}

####################################################################################
####################################################################################

# Step 2 Decide the best monitoring window length ## Do it in the Real application
# Now the best lambda has been selected




####################################################################################
# Create a big list for all scenarios included in the simulation section
create_all_scenarios = function(){
  Scenarios_list=list(
    case1 = list(
      para = list(seq1=c(0.2,c(0.1,0.2)),
                  seq2=c(0.2,c(0.1,0.5))),
      ps = c(2,2),    ## values of p at different period
      qs = c(0,0) ,   ## values of q at different period
      feps = rep('weibull', 2),   # distribution of epsilon at different period
      par1 = rep(.6, 2),         # values of par1 at different period
      par2 = rep(.7, 2)         # values of par2 at different period
    ),
    case2 = list(
      para = list(seq1=c(0.2,c(0.1)),
                  seq2=c(0.2,c(0.4))),
      ps = c(1,1),    ## values of p at different period
      qs = c(0,0) ,   ## values of q at different period
      feps = rep('gamma', 2),   # distribution of epsilon at different period
      par1 = rep(.5, 2),         # values of par1 at different period
      par2 = rep(.5, 2)         # values of par2 at different period
    ),
    case3 = list(
      para = list(seq1=c(0.2,c(0.1,0.2)),
                  seq2=c(0.2,c(0.1,0.5))),
      ps = c(1,2),    ## values of p at different period
      qs = c(1,0) ,   ## values of q at different period
      feps = rep('gamma', 2),   # distribution of epsilon at different period
      par1 = rep(.5, 2),         # values of par1 at different period
      par2 = rep(.5, 2)         # values of par2 at different period
    ),
    case4 = list(
      para = list(seq1=c(0.6,c(0.1,0.2)),
                  seq2=c(0.2,c(0.1,0.2))),
      ps = c(2,2),    ## values of p at different period
      qs = c(0,0) ,   ## values of q at different period
      feps = rep('weibull', 2),   # distribution of epsilon at different period
      par1 = rep(.8, 2),         # values of par1 at different period
      par2 = rep(.9, 2)         # values of par2 at different period
    )
  )
  return(Scenarios_list)
}

create_all_scenarios_for_case2 = function(){
  Scenarios_list=list(
    case2_v1 = list(
      para = list(seq1=c(0.2,c(0.1)),
                  seq2=c(0.2,c(0.2))),
      ps = c(1,1),    ## values of p at different period
      qs = c(0,0) ,   ## values of q at different period
      feps = rep('gamma', 2),   # distribution of epsilon at different period
      par1 = rep(.5, 2),         # values of par1 at different period
      par2 = rep(.5, 2)         # values of par2 at different period
    ),
    case2_v2 = list(
      para = list(seq1=c(0.2,c(0.1)),
                  seq2=c(0.2,c(0.4))),
      ps = c(1,1),    ## values of p at different period
      qs = c(0,0) ,   ## values of q at different period
      feps = rep('gamma', 2),   # distribution of epsilon at different period
      par1 = rep(.5, 2),         # values of par1 at different period
      par2 = rep(.5, 2)         # values of par2 at different period
    ),
    case2_v3 = list(
      para = list(seq1=c(0.2,c(0.1)),
                  seq2=c(0.2,c(0.6))),
      ps = c(1,1),    ## values of p at different period
      qs = c(0,0) ,   ## values of q at different period
      feps = rep('gamma', 2),   # distribution of epsilon at different period
      par1 = rep(.5, 2),         # values of par1 at different period
      par2 = rep(.5, 2)         # values of par2 at different period
    ),
    case2_v4 = list(
      para = list(seq1=c(0.2,c(0.1)),
                  seq2=c(0.2,c(0.8))),
      ps = c(1,1),    ## values of p at different period
      qs = c(0,0) ,   ## values of q at different period
      feps = rep('gamma', 2),   # distribution of epsilon at different period
      par1 = rep(.5, 2),         # values of par1 at different period
      par2 = rep(.5, 2)         # values of par2 at different period
    )
  )
  return(Scenarios_list)
}
######################################################################################################################
## Define the function count_delay()
## Input: single binary delta vector (Break indicator), t_cut (For TYPE I ERROR control), tau (True break time point)
## Output: Delay for this Single vector
count_delay =function(delta_vec,t_cut,tau){
  ## Need to find the first non-zero entry in delta after t_cut
  idx_eq_1= which(delta_vec==1) ## All the idx where delta ==1
  idx_delay = idx_eq_1[which(idx_eq_1>=t_cut)[1]]
  ## Delay
  Delay = idx_delay+M2-tau
  return(Delay)
}


######################################################################################################################


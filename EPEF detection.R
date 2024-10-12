library("ggplot2")
library(dplyr)
library(parallel)
library(doParallel)
source("./util_funcs/core_funs.R")
source("./util_funcs/user_define.R")
source("./util_funcs/core_funs_monitor.R")
source("./util_funcs/PEF_monitor_penalty_with_training.R")
source("./util_funcs/funcs_boot_testing.R")
source("./util_funcs/user_define_yz.R")
#########################################################################################
num_cores = detectCores()-1
cl=makeCluster(num_cores)
registerDoParallel(cl)
  dur.day1 = c(rgamma(2500,1,1),rgamma(1000,2,2)) ## Duration
  # ts.plot(dur.day1);ts.plot(Price_benchmark)
  
  #### Online detection algorithm
  tw = 1 ## Sliding window
  
  M1 <- 1500
  M2 <- 2500
  n_boot=2000
  num_drop = 10 ## Drop first 10 tpts when calculating the Maha dist
  p <- 20
  q <- 0
  pdim <- p+q + 1
  a <- 3.7
  model.order = c(p,q)
  pdim = sum(model.order)+1
  sparsity_num = 5; stopifnot(sparsity_num<=pdim)
  sparsity_thres = 0.05
  lambda_list=c(seq(0.1,1,by=0.1),seq(2,10,by=1),seq(20,100,by=10)) # lambda in scald penalty


  boot_thres = c(.90,.95,.99)[2]
  ## For Empirical method
  alpha_maha = (1-boot_thres)/(pdim+1);print(alpha_maha)
  alpha_B = (1-boot_thres)/(pdim+1);print(alpha_B)
  
  
  paras0 = factor(c("omega",paste("alpha",1:(pdim-1),sep='')),
                  levels=c("omega",paste("alpha",1:(pdim-1),sep='')))
  
  sim = dur.day1
  sim_n_total = length(sim) # number of observations in one day
  len_ts_mon = sim_n_total-M2
  
  
  
  sim_temp = sim
  
  lambda_candidates = foreach(lambda_temp = lambda_list,.combine = rbind,.packages = c('moments')) %dopar%
    lambda_tuning_indiv(lambda_temp,sim, # duration time series
                        model.order,
                        M1, # burning 
                        M2,# training
                        sparsity_num,sparsity_thres,a)
  
  
  lambda=lambda_candidates[(which.min(lambda_candidates[lambda_candidates[,3]==1,2])),1] ## Best lambda
  Perform_break_list=Perform_break(sim_temp,M2,model.order,
                                   lambda,a,M1,n_boot,len_ts_mon,alpha_B,cl)

  mon_orig = Perform_break_list$mon_orig
  mon_boot =Perform_break_list$mon_boot
  mon_boot_thres_Upper=Perform_break_list$mon_boot_thres_Upper
  mon_boot_thres_Lower=Perform_break_list$mon_boot_thres_Lower
  
  gen_pic_G_mon_list_BOOTSTRAP_only=gen_pic_G_mon_Boot_only(mon_orig,mon_boot_thres_Upper,
                                                            mon_boot_thres_Lower,
                                                            crit,M2,M1,pdim,
                                                            paras0,len_ts_mon)

  
  mon_BOOTSTRAP_only_logic=gen_pic_G_mon_list_BOOTSTRAP_only$mon_WB_logic ## Logic vector from BOOSTRAP only threshold
 
  
  ### Mahalanobis distance
  Ensemble_detection_fixed_boot_thres_temp_BOOTSTRAP_only = Ensemble_detection_fixed_boot_thres(mon_orig,len_ts_mon,pdim,
                                                                                                mon_boot,n_boot,num_drop,
                                                                                                cl,mahalanobis_for_Lapply,
                                                                                                alpha_maha,mon_BOOTSTRAP_only_logic
  )
  stopCluster(cl)
  
  Ensemble_logic=Ensemble_detection_fixed_boot_thres_temp_BOOTSTRAP_only$Ensemble_logic
  Ensemble_logic_smoothed = apply(embed(Ensemble_logic,tw),1,mean)
  
  break_tpt = ifelse(sum(Ensemble_logic_smoothed==1,na.rm = TRUE)==0,"No break",which(Ensemble_logic_smoothed==1)[1]+tw+M2)

break_tpt


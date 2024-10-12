# The code is for collecting EPEF break detection result from simulation study
# Author: Yanzhao Wang
# Date : 10/21/2022
####################################################################################
# Loading pakcages
# library(plyr)
library(dplyr)
library(parallel)
library(doParallel)
library(ggplot2)
# Sourcing use defined functions
wd = "~/Chapter1/code/"
setwd(wd)
num_cores = detectCores()-1
source("./EPEF_simulation/util_funcs/mix_logacd_sim.R")
source("./EPEF_simulation/util_funcs/core_funs.R")
source("./EPEF_simulation/util_funcs/user_define.R")
source("./EPEF_simulation/util_funcs/core_funs_monitor.R")
source("./EPEF_simulation/util_funcs/PEF_monitor_penalty_with_training.R")
source("./EPEF_simulation/util_funcs/funcs_boot_testing.R")
source("./EPEF_simulation/util_funcs/user_define_yz.R")
####################################################################################
# Parameters
case_idx = 4; # Which scenario you are looking at
tau_idx = 4; # Which break point you are looking at, (2700,3000,3500,4000)
Test_NON_stationary = TRUE

n_total_length = 7500 # total length of the time series
M1 <- 1500 # burning 
M2 <- 2500 # training
n_boot=2000
num_drop = 10 ## Drop first 10 tpts when calculating the Maha dist
non_stationary_coef_mat = rbind(c(0.05,0),
                                c(0.1,0),
                                c(0.05,-0.01),
                                c(0.1,-0.01)
)
non_stationary_coef=non_stationary_coef_mat[case_idx,]
tau_list = c(2700,3000,3500,4000)
tau = tau_list[tau_idx] # true break points

model.order = c(5,0); pdim =sum(model.order)+1
sparsity_num = 5; stopifnot(sparsity_num<=pdim)
sparsity_thres = 0.05
lambda_list=c(seq(0.1,1,by=0.1),seq(2,10,by=1),seq(20,100,by=10)) # lambda in scald penalty
a=3.7 # a in SCALD penalty
boot_thres = c(.90,.95,.99)[2] ## Type I error control

## For Empirical method
alpha_maha = (1-boot_thres)/(pdim+1);print(alpha_maha)
alpha_B = (1-boot_thres)/(pdim+1);print(alpha_B)

## For existing method only
alpha_W =(1-boot_thres)/pdim;print(alpha_W)


len_ts_mon = n_total_length-M2
paras0 = factor(c("omega",paste("alpha",1:(pdim-1),sep='')),
                levels=c("omega",paste("alpha",1:(pdim-1),sep='')))
### Scenario
Scenario_list= create_all_scenarios()
stopifnot(case_idx%in%1:4) ## Currently there are four cases
Scenario_temp = Scenario_list[[case_idx]]
Scenario_name = names(Scenario_list)[case_idx]
# cat(paste0("tau is ",tau),'\n',paste0("Working on ",Scenario_name))
### Generate the simulation right here!!!
para <- Scenario_temp$para
cpts <- tau/n_total_length    ## position of change point
ps <- Scenario_temp$ps    ## values of p at different period
qs <- Scenario_temp$qs    ## values of q at different period
feps <- Scenario_temp$feps   # distribution of epsilon at different period
par1 <- Scenario_temp$par1        # values of par1 at different period
par2 <- Scenario_temp$par2         # values of par2 at different period

set.seed(1125412)
B=500
Ensemble_fixed_thres_mat_BOOTSTRAP_only =c()
Ensemble_fixed_thres_mat_WIENER_only =c()
cl = makeCluster(num_cores)
registerDoParallel(cl)
for(i in 1:B){
  print(paste0("We are in case_",case_idx,"_tau_",tau,"_B_",i))
  if(Test_NON_stationary){
    sim = mix.logacd(n_total_length, 500, para, cpts, ps, qs, feps, par1, par2,
                     non_stationary_coef)$y
  }
  else{
    sim = mix.logacd(n_total_length, 500, para, cpts, ps, qs, feps, par1, par2,
                     non_stationary_coef=c(0,0))$y
  }
  lambda_candidates = foreach(lambda_temp = lambda_list,.combine = rbind,.packages = c('moments')) %dopar%
    lambda_tuning_indiv(lambda_temp,sim, # duration time series
                        model.order,
                        M1, # burning 
                        M2,# training
                        sparsity_num,sparsity_thres,a)
  
  
  lambda=lambda_candidates[(which.min(lambda_candidates[lambda_candidates[,3]==1,2])),1] ## Best lambda
  ### Break detection START
  Perform_break_list=Perform_break(sim,M2,model.order,
                                   lambda,a,M1,n_boot,len_ts_mon,alpha_B,cl)
  
  
  
  ### Break detection END
  
  ### Generate ggplot for G function with threshold START
  mon_orig = Perform_break_list$mon_orig
  mon_boot =Perform_break_list$mon_boot
  mon_boot_thres_Upper=Perform_break_list$mon_boot_thres_Upper
  mon_boot_thres_Lower=Perform_break_list$mon_boot_thres_Lower
  
  gen_pic_G_mon_list_BOOTSTRAP_only=gen_pic_G_mon_Boot_only(mon_orig,mon_boot_thres_Upper,
                                                            mon_boot_thres_Lower,
                                                            crit,M2,M1,pdim,
                                                            paras0,len_ts_mon)
  gen_pic_G_mon_list_WIENER_only=gen_pic_G_mon_Wiener_only(mon_orig,crit,alpha_W,M2,M1,pdim,
                                                           paras0,len_ts_mon)
  
  
  mon_BOOTSTRAP_only_logic=gen_pic_G_mon_list_BOOTSTRAP_only$mon_WB_logic ## Logic vector from BOOSTRAP only threshold
  mon_WIENER_only_logic=gen_pic_G_mon_list_WIENER_only$mon_WB_logic ## Logic vector from BOOSTRAP only threshold
  
  
  ### Mahalanobis distance
  Ensemble_detection_fixed_boot_thres_temp_BOOTSTRAP_only = Ensemble_detection_fixed_boot_thres(mon_orig,len_ts_mon,pdim,
                                                                                                mon_boot,n_boot,num_drop,
                                                                                                cl,mahalanobis_for_Lapply,
                                                                                                alpha_maha,mon_BOOTSTRAP_only_logic
  )
  
  Ensemble_detection_fixed_boot_thres_temp_WIENER_only = Ensemble_detection_fixed_boot_thres(mon_orig,len_ts_mon,pdim,
                                                                                             mon_boot,n_boot,num_drop,
                                                                                             cl,mahalanobis_for_Lapply,
                                                                                             alpha_maha,mon_WIENER_only_logic
  )
  
  Ensemble_fixed_thres_mat_BOOTSTRAP_only = cbind(Ensemble_fixed_thres_mat_BOOTSTRAP_only,Ensemble_detection_fixed_boot_thres_temp_BOOTSTRAP_only$Ensemble_logic)
  Ensemble_fixed_thres_mat_WIENER_only = cbind(Ensemble_fixed_thres_mat_WIENER_only,Ensemble_detection_fixed_boot_thres_temp_WIENER_only$Ensemble_logic)
}
stopCluster(cl)


if(Test_NON_stationary){
  file_name_prob_BOOTSTRAP_ONLY = paste0(Scenario_name,"_tau_",tau,"_alpha_",1-boot_thres,"_BOOTSTRAP_ONLY_NonStationary.csv")  
  file_name_prob_WIENER_ONLY = paste0(Scenario_name,"_tau_",tau,"_alpha_",1-boot_thres,"_WIENER_ONLY_NonStationary.csv")  
}
if(!Test_NON_stationary){
  file_name_prob_BOOTSTRAP_ONLY = paste0(Scenario_name,"_tau_",tau,"_alpha_",1-boot_thres,"_BOOTSTRAP_ONLY.csv")  
  file_name_prob_WIENER_ONLY = paste0(Scenario_name,"_tau_",tau,"_alpha_",1-boot_thres,"_WIENER_ONLY.csv")  
}
save_path = "./EPEF_simulation/Turing_saves"
write.csv(Ensemble_fixed_thres_mat_BOOTSTRAP_only,file = paste0(save_path,"/",file_name_prob_BOOTSTRAP_ONLY))
write.csv(Ensemble_fixed_thres_mat_WIENER_only,file = paste0(save_path,"/",file_name_prob_WIENER_ONLY))
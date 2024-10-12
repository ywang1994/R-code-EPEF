
gen_pic_raw = function(sim,M2,len_ts_mon,tau){
  pic_raw=tibble(
    time = 1:(M2+len_ts_mon),
    duration = c(sim)
  )%>%
    ggplot()+
    geom_line(aes(x=time,y=duration))+
    geom_vline(xintercept = tau,col="red",linetype=2)
  return(pic_raw)
}

Perform_break = function(sim,M2,model.order,
                         lambda,a,M1,n_boot,len_ts_mon,alpha_B,cl){
  dat <- sim[1:M2]
  p <- model.order[1]
  q <- model.order[2]
  # p <- ps[1]
  # q <- qs[1]
  pdim <- p+q+1
  est <- finitval(dat,p,q)
  if (q==0) est <- est[-(pdim+1)]
  
  psih <- psi_cal(dat,est,p,q)
  eps <- dat/exp(psih)
  ed <- eps_dist2(eps)
  momente <- c(ed$mue, ed$vare, ed$skewe, ed$kurte)
  difference <- log(momente[1]) - mean(log(eps))
  thehat = thehat_training(dat,est,M2, M1, p,q,momente,lambda,a,tau=0)
  thehat.used  = thehat
  thehat.used[1] =thehat.used[1] + difference

  
  lambda.used = lambda

  
  
  csv_temp = data.frame(V1=dat)
  boot_inds= 1:n_boot
  l_B_list = parLapply(cl,as.list(boot_inds),ts0_boot_corrected,
                       csv_temp=csv_temp,
                       len_ts_mon = len_ts_mon,
                       mean_l = 20,
                       r=2)
  
  l_B = do.call(cbind,l_B_list)
  mon_orig <- cal_G_func_AR(sim[1:(M2+len_ts_mon)],
                            thehat = thehat.used,
                            momente = momente,
                            M1=M1,M2=M2,
                            est = est,lambda=lambda.used,a,tau=0,p.lam.prime)
  
  
  sim_y_boot = rbind(matrix(rep(sim[1:M2],n_boot),ncol=n_boot),l_B)
  
  mon_boot_list = parLapply(cl,split(sim_y_boot,col(sim_y_boot)),cal_G_func_AR,
                            thehat = thehat.used,momente=momente,
                            M1=M1,M2=M2,
                            est = est,lambda=lambda.used,a,tau=0,p.lam.prime)
  
  mon_boot =do.call(cbind,mon_boot_list)
  mon_boot_threses = apply(simplify2array(mon_boot_list),
                           1:2,quantile,probs = c(alpha_B/2,1-alpha_B/2))
  
  
  return(list(mon_orig=mon_orig,
              mon_boot=mon_boot,
              mon_boot_thres_Upper=mon_boot_threses[2,,1],
              mon_boot_thres_Lower=mon_boot_threses[1,,1]))
}



gen_pic_G_mon_updated=function(mon_orig,
                               mon_boot_thres_Upper,
                               mon_boot_thres_Lower,
                               crit,alpha_W,M2,M1,pdim,
                               paras0,len_ts_mon){
  
  critical = uniroot(crit,c(0,5),alpha = alpha_W)
  mon_Wiener_thres = rep((M2-M1)**.5*(1+(1:len_ts_mon)/(M2-M1))*critical$root,pdim)
  mon_tbl = tibble(
    mon_orig0 = mon_orig,
    mon_boot_thres_U=mon_boot_thres_Upper,
    mon_boot_thres_L = mon_boot_thres_Lower,
    mon_Wiener_thres_U = mon_Wiener_thres,
    mon_Wiener_thres_L = -mon_Wiener_thres,
    mon_logic_B = (mon_orig<mon_boot_thres_L|mon_orig>mon_boot_thres_U)*1,
    mon_logic_W = (mon_orig<mon_Wiener_thres_L|mon_orig>mon_Wiener_thres_U)*1,
    mon_logic = mon_logic_B*mon_logic_W,
    # mon_p_val = (mon_p_val<0)*1,
    paras= rep(paras0,each = length(mon_orig)/length(paras0)),
    time = rep((M2+1):length(sim[1:(M2+len_ts_mon)]),length(paras0))
  )
  pic_G_thres_MAX = mon_tbl%>%
    ggplot()+
    geom_line(aes(x=time,y=mon_orig0,group=paras))+
    geom_ribbon(aes(x=time,ymin = mon_thres_L,ymax = mon_thres_U),fill='red',alpha=0.2)+
    facet_wrap(~paras)+
    ylab("Observed G")+
    theme(axis.text.x = element_text(angle = 90))
  pic_G_thres_Wiener = mon_tbl%>%
    ggplot()+
    geom_line(aes(x=time,y=mon_orig0,group=paras))+
    geom_ribbon(aes(x=time,ymin = mon_Wiener_thres_L,ymax = mon_Wiener_thres_U),fill='green',alpha=0.2)+
    facet_wrap(~paras)+
    ylab("Observed G")+
    theme(axis.text.x = element_text(angle = 90))
  pic_G_thres_boot = mon_tbl%>%
    ggplot()+
    geom_line(aes(x=time,y=mon_orig0,group=paras))+
    geom_ribbon(aes(x=time,ymin = mon_boot_thres_L,ymax = mon_boot_thres_U),fill='blue',alpha=0.2)+
    facet_wrap(~paras)+
    ylab("Observed G")+
    theme(axis.text.x = element_text(angle = 90))
  pic_G_logic_MAX = mon_tbl%>%
    ggplot()+
    geom_line(aes(x=time,y=mon_logic,group=paras))+
    facet_wrap(~paras)+
    ylab("Detection indicator")+
    theme(axis.text.x = element_text(angle = 90))
  return(list(pic_G_thres_MAX=pic_G_thres_MAX,
              pic_G_thres_boot = pic_G_thres_boot,
              pic_G_thres_Wiener=pic_G_thres_Wiener,
              pic_G_logic_MAX = pic_G_logic_MAX,
              mon_WB_logic = mon_tbl$mon_logic)
  )
}

gen_pic_G_mon_Boot_only=function(mon_orig,
                                 mon_boot_thres_Upper,
                                 mon_boot_thres_Lower,
                                 crit,M2,M1,pdim,
                                 paras0,len_ts_mon){
  mon_tbl = tibble(
    mon_orig0 = mon_orig,
    mon_boot_thres_U=mon_boot_thres_Upper,
    mon_boot_thres_L = mon_boot_thres_Lower,
    mon_logic = ((mon_orig>mon_boot_thres_U) |(mon_orig<mon_boot_thres_L) )*1,
    # mon_p_val = (mon_p_val<0)*1,
    paras= rep(paras0,each = length(mon_orig)/length(paras0)),
    time = rep((M2+1):length(sim[1:(M2+len_ts_mon)]),length(paras0))
  )
  pic_G_thres_MAX = mon_tbl%>%
    ggplot()+
    geom_line(aes(x=time,y=mon_orig0,group=paras))+
    geom_ribbon(aes(x=time,ymin = mon_thres_L,ymax = mon_thres_U),fill='red',alpha=0.2)+
    facet_wrap(~paras)+
    ylab("Observed G")+
    theme(axis.text.x = element_text(angle = 90))
  pic_G_thres_boot = mon_tbl%>%
    ggplot()+
    geom_line(aes(x=time,y=mon_orig0,group=paras))+
    geom_ribbon(aes(x=time,ymin = mon_boot_thres_L,ymax = mon_boot_thres_U),fill='blue',alpha=0.2)+
    facet_wrap(~paras)+
    ylab("Observed G")+
    theme(axis.text.x = element_text(angle = 90))
  pic_G_logic_MAX = mon_tbl%>%
    ggplot()+
    geom_line(aes(x=time,y=mon_logic,group=paras))+
    facet_wrap(~paras)+
    ylab("Detection indicator")+
    theme(axis.text.x = element_text(angle = 90))
  return(list(pic_G_thres_MAX=pic_G_thres_MAX,
              pic_G_thres_boot = pic_G_thres_boot,
              pic_G_logic_MAX = pic_G_logic_MAX,
              mon_WB_logic = mon_tbl$mon_logic)
  )
}

gen_pic_G_mon_Wiener_only=function(mon_orig,
                                   crit,alpha_W,M2,M1,pdim,
                                   paras0,len_ts_mon){
  
  critical = uniroot(crit,c(0,5),alpha = alpha_W)
  mon_Wiener_thres = rep((M2-M1)**.5*(1+(1:len_ts_mon)/(M2-M1))*critical$root,pdim)
  mon_tbl = tibble(
    mon_orig0 = mon_orig,
    mon_Wiener_thres_U = mon_Wiener_thres,
    mon_Wiener_thres_L = -mon_Wiener_thres,
    mon_logic = ((mon_orig>mon_Wiener_thres_U) |(mon_orig<mon_Wiener_thres_L) )*1,
    # mon_p_val = (mon_p_val<0)*1,
    paras= rep(paras0,each = length(mon_orig)/length(paras0)),
    time = rep((M2+1):length(sim[1:(M2+len_ts_mon)]),length(paras0))
  )
  pic_G_thres_Wiener = mon_tbl%>%
    ggplot()+
    geom_line(aes(x=time,y=mon_orig0,group=paras))+
    geom_ribbon(aes(x=time,ymin = mon_Wiener_thres_L,ymax = mon_Wiener_thres_U),fill='green',alpha=0.2)+
    facet_wrap(~paras)+
    ylab("Observed G")+
    theme(axis.text.x = element_text(angle = 90))
  pic_G_logic_MAX = mon_tbl%>%
    ggplot()+
    geom_line(aes(x=time,y=mon_logic,group=paras))+
    facet_wrap(~paras)+
    ylab("Detection indicator")+
    theme(axis.text.x = element_text(angle = 90))
  return(list(
    pic_G_thres_Wiener=pic_G_thres_Wiener,
    pic_G_logic_MAX = pic_G_logic_MAX,
    mon_WB_logic = mon_tbl$mon_logic)
  )
}


gen_G_hist_temp = function(tpt_temp,mon_boot,n_boot,paras0,mon_orig,M2){
  tl_boot = tibble(
    G_boot = c(t(mon_boot[(tpt_temp-M2)+(0:(pdim-1))*len_ts_mon,])),
    para= rep(paras0,each =n_boot))
  tl_orig = tibble(
    G_obs = mon_orig[(tpt_temp-M2)+(0:(pdim-1))*len_ts_mon],
    para = paras0
  )
  G_hist_temp = tl_boot%>%
    ggplot()+
    geom_histogram(aes(x=G_boot,group=para))+
    facet_wrap(~para)+
    geom_vline(data = tl_orig,aes(xintercept = G_obs,group=para),col='red')+
    facet_wrap(~para)
  return(G_hist_temp)
}




gen_pic_Maha_d2=function(mon_orig,mon_boot,len_ts_mon,num_drop,pdim
){
  
  Maha_d2 =c()
  for(tpt in (num_drop+1):len_ts_mon){
    multi_G_orig = (mon_orig[tpt+(0:(pdim-1))*len_ts_mon])
    multi_G_mat =  t(mon_boot[tpt+(0:(pdim-1))*len_ts_mon,])
    
    multi_G_mat_used = multi_G_mat
    multi_G_orig_used = multi_G_orig
    
    xbar = apply(multi_G_mat_used,2,mean)
    Sigma = cov(multi_G_mat_used)
    num_para = pdim
    
    d2=mahalanobis(multi_G_orig_used,center = xbar,cov = Sigma)
    # d2_dynamic = unname(quantile(mahalanobis(multi_G_mat_used,center = xbar,cov = Sigma),prob=1-alpha_maha))
    Maha_d2=c(Maha_d2,d2)
    # Dynamic_thres = c(Dynamic_thres,d2_dynamic)
  }
  Maha_d2 = c(rep(0,num_drop),Maha_d2)
  # Dynamic_thres = c(rep(0,num_drop),Dynamic_thres)
  
  
  d2_tbl = tibble(
    time = M2+1:len_ts_mon,
    # d2_dynamic = Dynamic_thres,
    d2=c(Maha_d2),
    # d2_logic = (Maha_d2>qchisq(1-alpha_maha,df=pdim)*1),
    # d2_logic_dynamic = Maha_d2 > d2_dynamic
  )
  
  pic_Maha_d2=d2_tbl%>%
    ggplot()+
    geom_line(aes(x=time,y=d2))+
    geom_hline(yintercept = qchisq(1-alpha_maha,df=pdim),col='red')
  
  # pic_Maha_d2_dynamic=d2_tbl%>%
  #   ggplot()+
  #   geom_line(aes(x=time,y=d2))+
  #   geom_line(aes(x=time,y=d2_dynamic),col='blue')
  # pic_Maha_d2_logic = d2_tbl%>%
  #   ggplot()+
  #   geom_line(aes(x=time,y=as.numeric(d2_logic)))
  # pic_Maha_d2_logic_dynamic = d2_tbl%>%
  #   ggplot()+
  #   geom_line(aes(x=time,y=as.numeric(d2_logic_dynamic)))
  return(list(
    # pic_Maha_d2=pic_Maha_d2,
    #           pic_Maha_d2_dynamic=pic_Maha_d2_dynamic,
    #           pic_Maha_d2_logic = pic_Maha_d2_logic,
    #           pic_Maha_d2_logic_dynamic=pic_Maha_d2_logic_dynamic,
    Maha_d2=Maha_d2
    # Dynamic_thres =Dynamic_thres,
    # mon_Maha_logic = Maha_d2>Dynamic_thres
  )
  )
}

Gen_scatter_plot = function(tpt_temp,mon_boot,mon_orig,pdim,len_ts_mon,n_boot,
                            paras0,M2){
  multi_G_orig = (mon_orig[tpt_temp-M2+(0:(pdim-1))*len_ts_mon])
  multi_G_mat =  t(mon_boot[tpt_temp-M2+(0:(pdim-1))*len_ts_mon,])
  G_tbl = data.frame(
    rbind(multi_G_mat,
          multi_G_orig)
  )
  colnames(G_tbl) = paste0("G_",paras0)
  cols = c(rep("black",n_boot),'red')
  pchs = c(rep(1,n_boot),19)
  pairs(G_tbl[,1:(pdim)],col=cols,pch=pchs)
}



mahalanobis_for_Lapply =function(x_list){
  ## compute mahalanobis distance once and for all
  x_orig = t(x_list[[1]])
  x_boot = t(x_list[[2]])
  d2=mahalanobis(x_orig,center=apply(x_boot,2,mean),cov=cov(x_boot))
  return(d2)
}




Ensemble_detection_fixed_boot_thres=function(mon_orig,len_ts_mon,pdim,
                                             mon_boot,n_boot,num_drop,
                                             cl,mahalanobis_for_Lapply,
                                             alpha_maha,
                                             mon_WB_logic
){
  
  mon_orig_array = array(mon_orig,dim=c(len_ts_mon,pdim)) 
  mon_boot_array = array(mon_boot,dim = c(len_ts_mon,pdim,n_boot))
  mon_boot_array_split = asplit(mon_boot_array,1)
  mon_boot_array_list =lapply((num_drop+1):len_ts_mon,function(x) list(
    mon_boot_array_split[[x]],mon_boot_array_split[[x]]))
  mon_orig_array_list = lapply((num_drop+1):len_ts_mon,function(x) list(
    (mon_orig_array[x,]),mon_boot_array_split[[x]]
  ))
  # cl=makeCluster(num_cores)
  d2_boot_list = parLapply(cl,mon_boot_array_list,mahalanobis_for_Lapply) ## n_boot X len_ts_mon
  d2_orig_list = parLapply(cl,mon_orig_array_list,mahalanobis_for_Lapply) ## 1 X len_ts_mon
  # stopCluster(cl)
  
  
  d2_threses = apply(simplify2array(d2_boot_list),2,quantile,prob = 1-alpha_maha)
  d2_orig_mat = unlist(d2_orig_list)
  
  # plot(x= M2+1:(len_ts_mon-num_drop),
  #      y=(d2_orig_mat),typ='l')
  # lines(x= M2+1:(len_ts_mon-num_drop),
  #       y=(d2_threses),col='red',lty=2)
  
  
  mon_WB_logic_mat= matrix(mon_WB_logic,ncol=pdim) ## len_ts_mon by pdim
  
  
  mon_WB_logic_mat_for_ensemble = (apply(mon_WB_logic_mat,1,sum)>=1)[-(1:num_drop)]
  # plot(x= M2+1:(len_ts_mon-num_drop),
  #      y=mon_WB_logic_mat_for_ensemble,typ='l')
  # 
  
  
  Ensemble_logic = c(matrix(rep(FALSE,num_drop),ncol = num_drop),
                     (mon_WB_logic_mat_for_ensemble)*
                       (d2_orig_mat>d2_threses))
  
  
  return(list(Ensemble_logic=Ensemble_logic,
              d2_orig_mat = d2_orig_mat,
              d2_threses =d2_threses))
}



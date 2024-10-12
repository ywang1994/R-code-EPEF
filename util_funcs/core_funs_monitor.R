

mov_rej <- function(y, critical, M2, N, Start){
  out <- which(y > critical)
  rej <- rep(NA,N)
  for(k in 1:N){
    if( sum(out < M2 + Start + k*50 ) == 0) {
      rej[k] <- 0
    }else{
      rej[k] <- 1
    }
  }
  return(rej)
}

crit <- function(c, alpha){
  k <- seq(0, 1000, 1)
  t1 <- 4/pi
  t2 <- -pi^2*(2*k+1)^2/(8*c^2)
  t3 <- sum((-1)^k/(2*k+1)*exp(t2))
  stat <- t1*t3
  return(1- stat - alpha)
}


delay_time <- function(y, critical, cp){
  out <- which(y > critical)
  if(length(out) ==0){
    return(-999)
  }else{
    return(out[1] - cp)
  }
}

cc <- seq(0,5,0.01)
aa <- rep(0, length(cc))
for(i in 1:length(cc)){
  aa[i] <- crit(cc[i], 0.05)
}
critical <- uniroot(crit, c(0,5), alpha=0.05)

##### hypothesis test
binary.hypo.test=function(Ymat, delta0){
  N=2*length(Ymat[,1]) # the number of subjects. 
  if(is.integer(delta0*N)==0){
    delta0=round(delta0*N)/N
  }
  ## Create the full potential outcome matrix = (Y_T(1), Y_T(0), Y_C(1), Y_C(0))
  Full.mat=rbind(c(0,0,0,0), c(0,0,1,0), c(0,1,0,0), c(0,1,1,0),
                 c(0,0,0,1), c(0,0,1,1), c(0,1,0,1), c(0,1,1,1),
                 c(1,0,0,0), c(1,0,1,0), c(1,1,0,0), c(1,1,1,0),
                 c(1,0,0,1), c(1,0,1,1), c(1,1,0,1), c(1,1,1,1))
  
  ### compute delta and Delta. 
  delta.mat=cbind(Full.mat[,1]-Full.mat[,2], Full.mat[,3]-Full.mat[,4])
  Delta.vec=delta.mat[,1]+delta.mat[,2]
  
  #### compute v_[sp]
  r_T1=Full.mat[,1]  # potential outcome for Subject 1 when he/she is treated.
  r_T2=Full.mat[,3]  # potential outcome for Subject 2 when he/she is treated.
  r_C1=Full.mat[,2]  # potential outcome for Subject 1 when he/she is control.
  r_C2=Full.mat[,4]  # potential outcome for Subject 2 when he/she is control.
  delta1=delta.mat[,1]
  delta2=delta.mat[,2]
  avg.r_T=(r_T1+r_T2)/2
  avg.r_C=(r_C1+r_C2)/2
  avg.delta=(delta1+delta2)/2
  
  # compute S_T[sp]^2
  S2_T=(r_T1-avg.r_T)^2+(r_T2-avg.r_T)^2
  
  # compute S_C[sp]^2
  S2_C=(r_C1-avg.r_C)^2+(r_C2-avg.r_C)^2
  
  # compute S_delta[sp]^2
  S2_delta=(delta1-avg.delta)^2+(delta2-avg.delta)^2
  
  # compute v_[sp]
  v=(2^2/N^2)*(S2_T+S2_C-S2_delta/2)
  
  ## There are four unique tables of C = (T00, T01, T10, T11) 
  # s=1 => C=(1,0,1,0)
  M1=sum(Ymat[,1]==0 & Ymat[,2]==0) ## compute M_s when s=1
  
  # s=2 => C=(0,1,1,0)
  M2=sum(Ymat[,1]==0 & Ymat[,2]==1)
  
  # s=3 => C=(1,0,0,1)
  M3=sum(Ymat[,1]==1 & Ymat[,2]==0)
  
  # s=4 => C=(0,1,0,1)
  M4=sum(Ymat[,1]==1 & Ymat[,2]==1)
  
  #########################################
  ##### solver based on "lpSolve"
  library("lpSolve")
  
  find.x=function(v, m, N, Delta, delta0){
    f.obj=v
    f.con=rbind(diag(rep(1,16)),
                c(rep(1,4), rep(0,12)),
                c(rep(0,4), rep(1,4), rep(0,8)),
                c(rep(0,8), rep(1,4), rep(0,4)),
                c(rep(0,12), rep(1,4)),
                Delta)
    f.dir=c(rep(">=", 16), rep("=", 5))
    f.rhs=c(rep(0,16), m, N*delta0)
    x.solve=lp ("max", f.obj, f.con, f.dir, f.rhs, int.vec=1:16)$solution
    
    return(x.solve)
  }
  
  
  
  #########################################
  ##### Hypothesis test of delta=delta0
  est.delta=sum(Ymat[,1]-Ymat[,2])/(N/2)
  
  #delta0=0
  max.x=find.x(v=v, m=c(M1,M2,M3,M4), N=N, Delta=Delta.vec, delta0=delta0)
  max.var=sum(v*max.x)
  
  deviate=(est.delta-delta0)/sqrt(max.var)
  pval=1-pnorm(abs(deviate))
  return(list(pvalue=pval, max.var=max.var, deviate=deviate))
}

##### Find the confidence interval
binary.conf.int=function(Ymat, alpha=0.05){
  N=2*length(Ymat[,1])
  est.delta=sum(Ymat[,1]-Ymat[,2])/(N/2)
  
  var.max=binary.hypo.test(Ymat=Ymat, delta0=est.delta)$max.var
  rough.int=c(est.delta-qnorm(1-alpha/2)*sqrt(var.max), est.delta+qnorm(1-alpha/2)*sqrt(var.max))
  
  ### find the lower end
  rough.int.lower.start=rough.int[1]-(rough.int[2]-rough.int[1])*0.1
  rough.int.lower.end=rough.int[1]+(rough.int[2]-rough.int[1])*0.1
  lower.start=round(rough.int.lower.start*N)
  lower.end=round(rough.int.lower.end*N)
  lower.seq=seq(lower.start, lower.end, by=1)
  
  pval.lower=rep(NA, length(lower.seq))
  for(i in 1:length(lower.seq)){
    delta0=lower.seq[i]/N
    pval.lower[i]=2*binary.hypo.test(Ymat=Ymat, delta0=delta0)$pvalue
  }
  conf.int.lower=min(lower.seq[pval.lower>=0.05])/N
  
  ### find the upper end
  rough.int.upper.start=rough.int[2]-(rough.int[2]-rough.int[1])*0.1
  rough.int.upper.end=rough.int[2]+(rough.int[2]-rough.int[1])*0.1
  upper.start=round(rough.int.upper.start*N)
  upper.end=round(rough.int.upper.end*N)
  upper.seq=seq(upper.start, upper.end, by=1)
  
  pval.upper=rep(NA, length(upper.seq))
  for(i in 1:length(upper.seq)){
    delta0=upper.seq[i]/N
    pval.upper[i]=2*binary.hypo.test(Ymat=Ymat, delta0=delta0)$pvalue
  }
  conf.int.upper=max(upper.seq[pval.upper>=0.05])/N
  return(list(exact=c(conf.int.lower, conf.int.upper), approx=rough.int))
}

# Ymat=cbind(test.matched.pairs$death5, test.matched.pairs$death5.1)
# #binary.conf.int(Ymat=Ymat, alpha=0.05)
# 
# N=2*length(Ymat[,1])
# est.delta=sum(Ymat[,1]-Ymat[,2])/(N/2)
# alpha=0.05
# var.max=binary.hypo.test(Ymat=Ymat, delta0=est.delta)$max.var
# rough.int=c(est.delta-qnorm(1-alpha/2)*sqrt(var.max), est.delta+qnorm(1-alpha/2)*sqrt(var.max))
# rough.int

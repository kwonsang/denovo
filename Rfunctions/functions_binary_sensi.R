
sort.new = function(x)
{
  temp = sort(unique(x))
  new = 1:length(temp)
  ret = rep(0,length(x))
  for(i in new)
  {
    ret[x == temp[i]] = i
  }
  ret
}	

##########################
# CRDbinary
################
CRDbinary = function(index, treatment, outcome, null = 0, DE = "both", alternative = "two.sided", conf.int = T, conf.level = .95, continuous.relax = F)
{
  
  require(Matrix)
  
  index = sort.new(index)
  outcome = 1*outcome
  if(any(outcome!=0 & outcome!= 1) == T)
  {
    stop("Outcomes must be binary")
  }
  
  if(alternative != "two.sided" & alternative != "greater" & alternative != "less")
  {
    stop("Alternative options are two.sided, greater, or less")
  }
  
  
  ns = table(index)
  ms = table(index[treatment*1==1]) 
  N.total = sum(ns)
  
  nostratum = length(unique(index))
  treatment = 1*treatment
  if(any(treatment != 1 & treatment != 0))
  {
    stop("Treatment Vector Must be Binary")
  }
  treatment = (1*treatment == 1)
  
  max.ATE = (sum(outcome[treatment]) + sum(!treatment) - sum(outcome[!treatment]))/length(outcome)
  min.ATE = (sum(outcome[treatment])  - sum(outcome[!treatment]) - sum(treatment))/length(outcome)
  if(null < min.ATE || null > max.ATE)
  {
    stop("Null is infeasible given observed data. No allocation of unobserved potential outcomes could satisfy it.")
  }
  
  gur = suppressWarnings(require(gurobi))
  
  
  wATE.per.strat = (ns)/sum(ns)*(tapply((treatment)*outcome, index, sum)/ms - tapply((1-treatment)*outcome, index, sum)/(ns-ms))
  ATE.est = sum(wATE.per.strat)
  
  ntnc = cbind(ms, ns-ms)
  mins = apply(ntnc, 1, min)
  if(!all(mins == 1))
  {
    stop("Each stratum must have either one treated and many controls or one control and many treateds")
  }
  
  if(abs(null) < 1)
  {
    null = round(null*N.total)
  }
  
  
  max.e = (N.total*ATE.est > null)
  
  vec.11 = tapply(outcome*treatment, index, sum)
  vec.01 = tapply((!treatment)*outcome, index, sum)
  vec.10 = tapply((treatment)*(!outcome), index, sum)
  vec.00 = tapply((!treatment)*(!outcome), index, sum)
  V = cbind(vec.11, vec.10, vec.00, vec.01)
  id = apply(V, 1, paste, collapse = "-")
  
  num.id =  sort.new(xtfrm(id))
  nosymm = length(unique(num.id))
  cc = table(num.id)
  bb = tapply(ns, num.id, mean)
  N.sr = sum(cc-1)
  m.01 = vec.01+1
  m.00 = vec.00+1
  m.10 = vec.10 + 1
  m.11 = vec.11 + 1
  
  mult.00 = tapply(m.00, num.id, mean)
  mult.01 = tapply(m.01, num.id, mean)
  mult.11 = tapply(m.11, num.id, mean)
  mult.10 = tapply(m.10, num.id, mean)
  mult.all = mult.11*mult.10*mult.00*mult.01
  if(DE == "nonpositive")
  {
    mult.all = mult.10*mult.01
  }
  if(DE=="nonnegative")
  {
    mult.all = mult.11*mult.00
  }
  ns.type = tapply(ns, num.id, mean)
  N.vars = sum(mult.all)
  index.symm = rep(1:nosymm, mult.all)
  n.types = (mult.all)
  Diff = rep(0, N.vars)
  #V.list = vector("list", N.vars)
  
  PV = rep(0, N.vars)
  n.per = cc
  row.ind = rep(0, 2*N.vars)
  col.ind = row.ind
  values = row.ind
  b = rep(0, nosymm + 1)
  for(kk in 1:nosymm)
  {
    row.ind[which(index.symm==kk)] = rep(kk, n.types[kk])
    col.ind[which(index.symm==kk)] = which(index.symm==kk)  
    values[which(index.symm==kk)] = rep(1, (n.types[kk]))
    b[kk] = n.per[kk]
  }
  row.ind[(N.vars+1):(2*N.vars)] = rep(nosymm + 1, N.vars)
  col.ind[(N.vars+1):(2*N.vars)] = 1:N.vars
  #row.ind[(2*N.vars+1):(3*N.vars)] = rep(nosymm + 2, N.vars+1)
  
  Gamma.sens = 1
  
  for(kk in 1:nosymm)
  {
    i = which(num.id==kk)[1]
    symmgroup = which(index.symm == kk)
    ind = which(index==i)
    treatstrat = treatment[ind]
    outstrat = outcome[ind]
    outsymm = c(sort(outstrat[treatstrat==F]), sort(outstrat[treatstrat==T]))
    PO = matrix(0, n.types[kk], ns[i])
    count = 1
    mt.01 = mult.01[kk]
    mt.00 = mult.00[kk]
    mt.10 = mult.10[kk]
    mt.11 = mult.11[kk]
    T.00 = matrix(1, mt.00-1, mt.00-1)
    T.00[lower.tri(T.00)] = 0
    T.00 = rbind(T.00, c(rep(0, mt.00-1)) )
    T.01 = matrix(1, mt.01-1, mt.01-1)
    T.01[lower.tri(T.01)] = 0
    T.01 = rbind(T.01, c(rep(0, mt.01-1)) )
    
    T.10 = matrix(1, mt.10-1, mt.10-1)
    T.10[lower.tri(T.10)] = 0
    T.10 = rbind(T.10, c(rep(0, mt.10-1)) )
    T.11 = matrix(1, mt.11-1, mt.11-1)
    T.11[lower.tri(T.11)] = 0
    T.11 = rbind(T.11, c(rep(0, mt.11-1)) )
    c.00 = mt.00
    c.01 = mt.01
    c.10 = mt.10
    c.11 = mt.11
    if(DE == "nonnegative")
    {
      T.10 = matrix(0, 1, mt.10-1)
      c.10 = 1
      T.01 = matrix(1, 1, mt.01-1)
      c.01 = 1
    }
    if(DE == "nonpositive")
    {
      T.11 = matrix(1, 1, mt.11-1)
      c.11 = 1
      T.00 = matrix(0, 1, mt.00-1)
      c.00 = 1
    }
    
    for(ll in 1:c.00)
    {
      for(mm in 1:c.01)
      {
        for(jj in 1:c.10)
        {
          for(oo in 1:c.11)
          {
            
            tempvec = c(T.00[ll,], T.01[mm,], T.10[jj,], T.11[oo,])
            PO[count,] = tempvec[!is.na(tempvec)]
            count = count+1
          }
        }
        
      }
    }
    for(jj in 1:n.types[kk])
    {
      ind.jj = (((jj-1)*(ns[i]-1))+1):(jj*(ns[i]-1))
      po.symm = PO[jj,]
      treatsymm = c(rep(F, ns[i]-ms[i]), rep(T, ms[i]))
      outcontrol = outsymm*(1-treatsymm) + po.symm*(treatsymm)
      outtreat = outsymm*(treatsymm) + po.symm*(1-treatsymm)           
      sum.cont = sum(outcontrol)/(ns[i]-1)
      Q = (outtreat + outcontrol/(ns[i]-1) - sum.cont)*ns[i]
      if(sum(treatstrat)>1)
      {
        sum.cont = sum(outtreat)/(ns[i]-1)
        Q = -(outtreat/(ns[i]-1) + outcontrol - sum.cont)*ns[i]  
      }
      
      qi = Q*max.e - Q*(!max.e)
      ord = order(qi)
      qi.sort = sort(qi)
      
      
      Gamma.sens = 1
      
      
      
      mu = mean(qi.sort)
      sigma2 = mean(qi.sort^2) - mu^2
      
      mu[abs(mu) < 1e-8] = 0
      sigma2[sigma2 < 1e-8] = 0
      
      PV[symmgroup[jj]] = (sigma2)
      Diff[symmgroup[jj]] = sum(outtreat - outcontrol)  
    }
  }
  
  values[(N.vars+1):(2*N.vars)] = Diff
  b[nosymm+1] = 0
  const.dir = c(rep("=", nosymm+1))
  A = sparseMatrix(row.ind, col.ind, x=values)
  
  res = ATEtest(ATE.est, N.total, b, A, const.dir,PV, Diff, null, alternative)
  pval = res$pval
  tstat = res$tstat
  SE = res$SE
  TE.wald = round(ATE.est*N.total)
  if(conf.int == T)
  {
    alpha.temp = (1-conf.level)/2
    
    SE.wald =  ATEtest(ATE.est, N.total, b,A, const.dir,PV, Diff, TE.wald, alternative)$SE
    ub = round(TE.wald - qnorm(alpha.temp)*SE.wald)
    lb = round(TE.wald + qnorm(alpha.temp)*SE.wald)
    if(DE == "nonpositive")
    {
      ub = min(0, ub)
    }
    if(DE == "nonnegative")
    {
      lb = max(0, lb)
    }
    
    upval = ATEtest(ATE.est, N.total, b, A, const.dir,PV, Diff, ub, "less")$pval
    RU = (upval < alpha.temp)
    diff = 10
    DU = (-1)*RU + 1*(!RU)
    ex = !(RU)
    if(DE == "nonpositive" & RU == F & ub == 0)
    {
      diff = 1
      ex = 0
    }
    while(diff != 1)
    {
      ub = ub + DU
      upval = ATEtest(ATE.est, N.total, b, A, const.dir,PV, Diff, ub, "less")$pval
      RU1 = (upval < alpha.temp)
      diff = abs(RU1 - RU)	
      RU = RU1  		
      if(DE == "nonpositive" & ub==0 & RU == F)
      {
        diff=1
        ex = 0
      }
    }
    ub = ub - ex
    
    
    loval = ATEtest(ATE.est, N.total, b, A, const.dir,PV, Diff, lb, "greater")$pval
    RU = (loval < alpha.temp)
    diff = 10
    DU = (1)*RU + -1*(!RU)
    ex = !(RU)
    if(DE == "nonnegative" & RU == F & lb == 0)
    {
      diff = 1
      ex = 0
    }
    while(diff != 1)
    {
      lb = lb + DU
      loval = ATEtest(ATE.est, N.total, b, A, const.dir, PV, Diff, lb, "greater")$pval
      RU1 = (loval < alpha.temp)
      diff = abs(RU1 - RU)
      if(DE == "nonnegative" & lb==0 & RU == F)
      {
        diff=1
        ex = 0
      }
      RU = RU1  		
    }
    lb = lb + ex
    
    CI = c(lb, ub) 	
    return(list(ATE.est = ATE.est, SE = SE/N.total, null = null/N.total, pval = pval, alternative = alternative, CI = CI/N.total, conf.level = conf.level))
  }
  else
  {
    
    return(list(ATE.est = ATE.est, SE = SE/N.total, null = null/N.total, pval = pval, alternative = alternative, CI = NULL, conf.level = NULL))
  }
} 



ATEtest = function(ATE.est, N.total, b, A, const.dir,PV, Diff, null, alternative, continuous.relax = F)
{
  nosymm = length(b)-1
  N.vars = length(PV)
  b[nosymm+1] = null
  
  model = list()
  
  model$A = A  	
  model$obj = c(PV)
  model$sense = const.dir
  model$rhs = b
  model$vtype = c(rep("I", N.vars))
  if(continuous.relax == T){model$vtype = c(rep("C", N.vars))}
  model$modelsense = "max"
  
  solm = gurobi(model, params = list(OutputFlag = 0))
  SE = sqrt(solm$objval)
  
  tstat = (N.total*ATE.est - null)/SE
  
  if(alternative == "two.sided")
  {
    pval = 2*pnorm(-abs(tstat))
  }
  if(alternative == "greater")
  {
    pval = 1 - pnorm((tstat))
  }
  if(alternative == "less")
  {
    pval = pnorm((tstat))
  }
  return(list(tstat = tstat, SE = SE, pval = pval))
}



########
#sensCRD
#########

sensCRD = function(index, treatment, outcome, null = 0, DE = "both", alternative = "two.sided", alpha = 0.05, Gamma.vec = 1, calculate.pval = T, continuous.relax = F)
{
  
  PVAL = calculate.pval
  ns = table(index)
  ms = table(index[treatment*1==1]) 
  N.total = sum(ns)
  nostratum = length(unique(index))
  treatment = 1*treatment
  
  treatment = (1*treatment == 1)
  
  wATE.per.strat = (ns)/sum(ns)*(tapply((treatment)*outcome, index, sum)/ms - tapply((1-treatment)*outcome, index, sum)/(ns-ms))
  ATE.est = sum(wATE.per.strat)
  max.e = (N.total*ATE.est > null)
  
  vec.11 = tapply(outcome*treatment, index, sum)
  vec.01 = tapply((!treatment)*outcome, index, sum)
  vec.10 = tapply((treatment)*(!outcome), index, sum)
  vec.00 = tapply((!treatment)*(!outcome), index, sum)
  V = cbind(vec.11, vec.10, vec.00, vec.01)
  id = apply(V, 1, paste, collapse = "-")
  
  num.id =  sort.new(xtfrm(id))
  nosymm = length(unique(num.id))
  cc = table(num.id)
  bb = tapply(ns, num.id, mean)
  N.sr = sum(cc-1)
  m.01 = vec.01+1
  m.00 = vec.00+1
  m.10 = vec.10 + 1
  m.11 = vec.11 + 1
  
  mult.00 = tapply(m.00, num.id, mean)
  mult.01 = tapply(m.01, num.id, mean)
  mult.11 = tapply(m.11, num.id, mean)
  mult.10 = tapply(m.10, num.id, mean)
  mult.all = mult.11*mult.10*mult.00*mult.01
  if(DE == "nonpositive")
  {
    mult.all = mult.10*mult.01
  }
  if(DE=="nonnegative")
  {
    mult.all = mult.11*mult.00
  }
  ns.type = tapply(ns, num.id, mean)
  N.vars = sum(mult.all*(ns.type-1))
  index.symm = rep(1:nosymm, mult.all*(ns.type-1))
  n.types = (mult.all)
  Diff = rep(0, N.vars)
  #V.list = vector("list", N.vars)
  
  PM = rep(0, N.vars)
  PV = rep(0, N.vars)
  n.per = cc
  row.ind = rep(0, 3*N.vars+1)
  col.ind = row.ind
  values = row.ind
  b = rep(0, nosymm + 2)
  for(kk in 1:nosymm)
  {
    row.ind[which(index.symm==kk)] = rep(kk, (ns.type[kk]-1)*n.types[kk])
    col.ind[which(index.symm==kk)] = which(index.symm==kk)  
    values[which(index.symm==kk)] = rep(1, (ns.type[kk]-1)*(n.types[kk]))
    b[kk] = n.per[kk]
  }
  row.ind[(N.vars+1):(2*N.vars)] = rep(nosymm + 1, N.vars)
  col.ind[(N.vars+1):(2*N.vars)] = 1:N.vars
  row.ind[(2*N.vars+1):(3*N.vars+1)] = rep(nosymm + 2, N.vars+1)
  col.ind[(2*N.vars+1):(3*N.vars+1)] = 1:(N.vars+1)
  zscore = rep(0, length(Gamma.vec))
  Rejectvec = zscore
  pvalvec = zscore
  kappavec = zscore
  
  max.var.vec=rep(0, length(Gamma.vec))
  mean.vec=rep(0, length(Gamma.vec))
  test.stat.vec=rep(0, length(Gamma.vec))
  for(ee in 1:length(Gamma.vec))
  {
    Gamma.sens = Gamma.vec[ee]
    
    for(kk in 1:nosymm)
    {
      i = which(num.id==kk)[1]
      symmgroup = which(index.symm == kk)
      ind = which(index==i)
      treatstrat = treatment[ind]
      outstrat = outcome[ind]
      outsymm = c(sort(outstrat[treatstrat==F]), sort(outstrat[treatstrat==T]))
      PO = matrix(0, n.types[kk], ns[i])
      count = 1
      mt.01 = mult.01[kk]
      mt.00 = mult.00[kk]
      mt.10 = mult.10[kk]
      mt.11 = mult.11[kk]
      T.00 = matrix(1, mt.00-1, mt.00-1)
      T.00[lower.tri(T.00)] = 0
      T.00 = rbind(T.00, c(rep(0, mt.00-1)) )
      T.01 = matrix(1, mt.01-1, mt.01-1)
      T.01[lower.tri(T.01)] = 0
      T.01 = rbind(T.01, c(rep(0, mt.01-1)) )
      
      T.10 = matrix(1, mt.10-1, mt.10-1)
      T.10[lower.tri(T.10)] = 0
      T.10 = rbind(T.10, c(rep(0, mt.10-1)) )
      T.11 = matrix(1, mt.11-1, mt.11-1)
      T.11[lower.tri(T.11)] = 0
      T.11 = rbind(T.11, c(rep(0, mt.11-1)) )
      c.00 = mt.00
      c.01 = mt.01
      c.10 = mt.10
      c.11 = mt.11
      if(DE == "nonnegative")
      {
        T.10 = matrix(0, 1, mt.10-1)
        c.10 = 1
        T.01 = matrix(1, 1, mt.01-1)
        c.01 = 1
      }
      if(DE == "nonpositive")
      {
        T.11 = matrix(1, 1, mt.11-1)
        c.11 = 1
        T.00 = matrix(0, 1, mt.00-1)
        c.00 = 1
      }
      
      
      for(ll in 1:c.00)
      {
        for(mm in 1:c.01)
        {
          for(jj in 1:c.10)
          {
            for(oo in 1:c.11)
            {
              
              tempvec = c(T.00[ll,], T.01[mm,], T.10[jj,], T.11[oo,])
              PO[count,] = tempvec[!is.na(tempvec)]
              count = count+1
            }
          }
          
        }
      }
      for(jj in 1:n.types[kk])
      {
        ind.jj = (((jj-1)*(ns[i]-1))+1):(jj*(ns[i]-1))
        po.symm = PO[jj,]
        treatsymm = c(rep(F, ns[i]-ms[i]), rep(T, ms[i]))
        outcontrol = outsymm*(1-treatsymm) + po.symm*(treatsymm)
        outtreat = outsymm*(treatsymm) + po.symm*(1-treatsymm) 
        sum.cont = sum(outcontrol)/(ns[i]-1)
        Q = (outtreat + outcontrol/(ns[i]-1) - sum.cont)*ns[i]
        if(sum(treatstrat)>1)
        {
          sum.cont = sum(outtreat)/(ns[i]-1)
          Q = -(outtreat/(ns[i]-1) + outcontrol - sum.cont)*ns[i]  
        }
        qi = Q*max.e - Q*(!max.e)
        ord = order(qi)
        qi.sort = sort(qi)
        
        
        mu = rep(0, length(ind)-1)
        sigma2 = rep(0, length(ind)-1)
        
        
        for(j in 1:(length(ind)-1))
        {
          mu[j] = (sum(qi.sort[1:(j)]) + Gamma.sens*sum(qi.sort[(j+1):(length(ind))]))/((j) + Gamma.sens*(ns[i]-(j)))
          sigma2[j] = (sum(qi.sort[1:(j)]^2) + Gamma.sens*sum(qi.sort[(j+1):(length(ind))]^2))/((j) + Gamma.sens*(ns[i]-(j))) - mu[j]^2
        }
        mu[abs(mu) < 1e-8] = 0
        sigma2[sigma2 < 1e-8] = 0
        
        
        PM[symmgroup[ind.jj]] = mu*(max.e) - mu*(!max.e)
        PV[symmgroup[ind.jj]] = (sigma2)
        Diff[symmgroup[ind.jj]] = sum(outtreat - outcontrol)
        
        
      }
    }
    
    values[(N.vars+1):(2*N.vars)] = Diff
    values[(2*N.vars+1):(3*N.vars+1)] = c(-PM, 1)
    b[nosymm+1] = null
    b[nosymm+2] = 0
    
    const.dir = c(rep("=", nosymm+2))
    model = list()
    alpha.opt = alpha
    if(alternative != "two.sided")
    {
      alpha.opt = 2*alpha
    }
    if(Gamma.sens==1)
    {
      model$A = sparseMatrix(row.ind[1:(2*N.vars)], col.ind[1:(2*N.vars)], x=values[1:(2*N.vars)])
      model$obj = c(PV)
      model$sense = const.dir[1:(nosymm+1)]
      model$rhs = b[1:(nosymm+1)]
      model$vtype = c(rep("I", N.vars))
      if(continuous.relax == T){model$vtype = c(rep("C", N.vars))}
      model$modelsense = "max"
      
      
      solm = gurobi(model, params = list(OutputFlag = 0))
      zed = ((N.total*ATE.est - null)/sqrt(solm$objval))
      tstat = zed
      kappa = (N.total*ATE.est - null)^2 - qchisq(1-alpha.opt, 1)*solm$objval
      kappavec[ee] = kappa
      max.var=solm$objval
      mean.val=null
      pval = 0
      if(alternative == "two.sided")
      {
        pval = 2*pnorm(-abs(tstat))
      }
      if(alternative == "greater")
      {
        pval = 1 - pnorm((tstat))
      }
      if(alternative == "less")
      {
        pval = pnorm((tstat))
      }
      Reject = (pval < alpha)
    }
    if(Gamma.sens != 1)
    {
      diff = 200
      kappa = qchisq(1-alpha.opt, 1)
      count=0
      sim.count=0
      while(diff > 1e-8 & sim.count <= 100)
      {
        Plin = -2*N.total*ATE.est*PM - kappa*PV 
        rowind.q =  1:(N.vars+1)
        colind.q = 1:(N.vars+1)
        values.q = c(rep(0, N.vars),1)
        Q = sparseMatrix(rowind.q, colind.q, x=values.q)
        model$A = sparseMatrix(row.ind, col.ind, x=values)
        model$obj = c(Plin,0)
        model$Q = Q
        model$sense = const.dir
        model$rhs = b
        model$vtype = c(rep("I", N.vars), "C")
        if(continuous.relax == T){model$vtype = c(rep("C", N.vars+1))}
        model$lb = c(rep(0, N.vars), -Inf)
        model$modelsense = "min"
        solm = gurobi(model, params = list(OutputFlag = 0))
        x = solm$x[1:N.vars]
        if(sum(PV*x)==0){
          kappa.new=0
        }else{
          kappa.new = (N.total*ATE.est - sum(PM*x))^2/sum(PV*x)
        }
        diff = abs(kappa.new - kappa)
        pval = 0
        if(PVAL == F)
        {
          diff = 0
          Reject = (kappa.new > kappa)
        }
        kappa = kappa.new
        kappavec[ee] = (N.total*ATE.est - sum(PM*x))^2 - qchisq(1-alpha.opt, 1)*sum(PV*x)
        sim.count=sim.count+1
      }
      
      zed = sqrt((N.total*ATE.est - sum(PM*x))^2/sum(PV*x))
      
      if(alternative == "less")
      {
        zed = -zed
      }
      zscore[ee] = zed
      tstat = zed
      max.var=sum(PV*x)
      mean.val=(sum(PM*x))
      
      if(PVAL == T)
      {
        if(alternative == "two.sided")
        {
          pval = 2*pnorm(-abs(tstat))
        }
        if(alternative == "greater")
        {
          pval = 1 - pnorm((tstat))
        }
        if(alternative == "less")
        {
          pval = pnorm((tstat))
        }
        Reject = (pval < alpha)
      }
      
      
      if(sign(-sum(PM*x) + N.total*ATE.est)!=sign(N.total*ATE.est - null))
      {
        Reject = F
        pval = 0.5
        if(alternative == "two.sided")
        {
          pval = 1
        }
      }
      
      if(alternative == "greater" & sum(PM*x) < null)
      {
        pval = .5
      }
      if(alternative == "less" & sum(PM*x) > null)
      {
        pval = .5
      }
      
      
      
    }
    pvalvec[ee] = pval
    Rejectvec[ee] =  Reject
    
    max.var.vec[ee]=max.var
    test.stat.vec[ee]=tstat
    mean.vec[ee]=mean.val
  }
  if(PVAL == F)
  {
    return(list(Gamma.vec = Gamma.vec, Reject = Rejectvec, kappa = kappavec))
  }
  if(PVAL == T)
  {
    return(list(Gamma.vec = Gamma.vec, pval = pvalvec, test.stat=test.stat.vec, maxvar=max.var.vec, mean=mean.vec))
  }
}  
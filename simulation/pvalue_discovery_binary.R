
#################################################################
## Compared to the Causal Tree method (Athey and Imbens, PNAS)
#################################################################
library(causalTree)
library(mvtnorm)
library(rpart.utils)
library(sensitivitymv)

source("../Rfunctions/basic_functions.R")
source("../Rfunctions/functions_binary.R")


########################################################

nsim=1

pval.matrix=matrix(NA, nrow=nsim, ncol=12)
check.matrix=matrix(NA, nrow=nsim, ncol=30)

time.before=Sys.time()
for(k in 1:nsim){
  # Data generation
  n=2000
  x1=rbinom(n, 1, prob=0.5)
  x2=rbinom(n, 1, prob=0.5)
  x3=rbinom(n, 1, prob=0.5)
  x4=rbinom(n, 1, prob=0.5)
  x5=rbinom(n, 1, prob=0.5)
  
  # for(i in 1:10){
  #   assign(paste0("x", i), rbinom(n, 1, prob=0.5))
  # }
  
  ### First case: no effect modification
  y_t=rep(NA, n)
  sd=sqrt(1/2)
  #y_t=rnorm(n, mean=0, sd=sd)
  y_t[x1==0]=rbinom(sum(x1==0), 1, prob=0.4)
  y_t[x1==1]=rbinom(sum(x1==1), 1, prob=0.6)
  
  y_c=rbinom(n, 1, prob=0.5)
  matched.data=cbind(y_t, y_c, x1, x2, x3, x4, x5)
  full.data=c(y_t, y_c)
  full.data=as.data.frame(full.data)
  full.data$z=c(rep(1,n), rep(0,n))
  full.data$x1=c(x1,x1)
  full.data$x2=c(x2,x2)
  full.data$x3=c(x3,x3)
  full.data$x4=c(x4,x4)
  full.data$x5=c(x5,x5)
  colnames(full.data)[1]="y_obs"
  
  training.data.index=sample(1:n, n/10, replace=F)
  train=c(training.data.index, training.data.index+n)
  ##################################################
  ##### 1.Select training sample & Est sample (10% vs. 90%)
  ##################################################
  tra.matched.data=matched.data[training.data.index, ]
  est.matched.data=matched.data[-training.data.index,]
  tra.matched.data=as.data.frame(tra.matched.data)
  est.matched.data=as.data.frame(est.matched.data)
  
  traData=full.data[train,] # training data
  testData=full.data[-train, ] # test data
  
  ######### 1.1. Discovering tree structures from causal tree
  chosen.check=rep(0,5)
  ### training sample -> Create tree using causaltree
  tree=causalTree(y_obs ~ x1+x2+x3+x4+x5, data=traData, treatment = traData$z, HonestSampleSize=length(testData[,1]),
                  split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, xval = 5, 
                  cp = 0, minsize = 20, propensity = 0.5)
  
  opcp=tree$cptable[,1][which.min(tree$cptable[,4])]
  opfit=prune(tree, opcp) ## obtained tree from a training sample
  
  var.names=opfit$frame$var[opfit$frame$var!="<leaf>"]
  var.numbers=as.numeric(substr(var.names, start=2, stop=4))
  chosen.check.matrix.causaltree=chosen.check
  chosen.check.matrix.causaltree[var.numbers]=1
  
  
  ######## using causaltree
  ### Test sample -> hypothesis test
  est.tree.stru=tree.stru(opfit, est.matched.data)
  subgroup.num=est.tree.stru$new.x$subgroup.num
  C=est.tree.stru$C
  unique.subg=as.numeric(colnames(est.tree.stru$C))
  num.of.subgroups=length(unique.subg)
  
  test.matched.pairs.ymat=est.matched.data
  test.matched.pairs.ymat$sub.num=subgroup.num
  
  ##########################################
  n.test=length(test.matched.pairs.ymat[,1])
  n.vec=table(test.matched.pairs.ymat$sub.num)
  
  # Since the population treatment effect size (mu) is not known, 
  ## get the 99% confidence interval of delta
  est.delta=sum(test.matched.pairs.ymat[,1]-test.matched.pairs.ymat[,2])/n.test
  alpha=0.01
  var.max=binary.hypo.test(Ymat=test.matched.pairs.ymat, delta0=est.delta)$max.var
  rough.int=c(est.delta-qnorm(1-alpha/2)*sqrt(var.max), est.delta+qnorm(1-alpha/2)*sqrt(var.max))
  #rough.int
  delta.end.points=round(rough.int*n*2)
  delta.seq=round(seq(delta.end.points[1], delta.end.points[2], length.out=20))
  
  ### subgroup estimates
  est.delta.vec=rep(NA, length(C[1,]))
  for(j in 1:length(C[1,])){
    est.delta.vec[j]=mean(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1]-test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],2])
  }
  #est.delta.vec
  
  library(mvtnorm)
  max.dev.vec=max.dev.pval.vec=pval.vec=rep(NA, length(delta.seq))
  for(i in 1:length(delta.seq)){
    temp.delta0=delta.seq[i]
    
    ### submax 
    temp.maxvar.vec=rep(NA, length(C[1,]))
    temp.delta0.vec=rep(NA, length(C[1,]))
    temp.pval.vec=rep(NA, length(C[1,]))
    for(j in 1:length(C[1,])){
      fraction.val=(temp.delta0/(2*n.test))*2*n.vec[j]
      closest.two.points=c(floor(fraction.val), ceiling(fraction.val))
      
      temp.hypo1=binary.hypo.test(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1:2], delta0=closest.two.points[1]/(2*n.vec[j]))
      temp.hypo2=binary.hypo.test(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1:2], delta0=closest.two.points[2]/(2*n.vec[j]))
      temp.maxvar1=temp.hypo1$max.var
      temp.maxvar2=temp.hypo2$max.var
      max.temp.maxvar=max(temp.maxvar1, temp.maxvar2)
      temp.maxvar.vec[j]=max.temp.maxvar
      temp.delta0.vec[j]=closest.two.points[which.max(c(temp.maxvar1, temp.maxvar2))]
      
      temp.pval.vec[j]=max(temp.hypo1$pvalue, temp.hypo2$pvalue)*2
    }
    
    T.vec=est.delta.vec*n.vec*2
    test= C %*% T.vec
    if(dim(C)[1]==1){
      var.mat=as.numeric(temp.maxvar.vec*(n.vec^2)*4)
      theta0= temp.delta0.vec
      om0= var.mat
      corr0= C
      deviate=(test-theta0)/sqrt(om0)
    }else{
      var.mat=diag(temp.maxvar.vec*(n.vec^2)*4)
      theta0= C %*% temp.delta0.vec
      om0=C %*% var.mat %*% t(C)
      deviate=(test-theta0)/sqrt(diag(om0))
      corr0=om0/sqrt(outer(diag(om0),diag(om0),"*"))
    }
    
    max.deviate=max(abs(deviate))
    
    max.dev.vec[i]=max.deviate
    if(dim(C)[1]==1){
      max.dev.pval.vec[i]=2*(1-pnorm(max.deviate))
    }else{
      max.dev.pval.vec[i]=2*(1-pmvnorm(lower=-Inf, upper=max.deviate, mean=rep(0, length(corr0[,1])), corr=corr0)[1])
    }
    pval.vec[i]=truncatedP(temp.pval.vec, trunc=0.1)
    #critical.constant.vec[i]=qmvnorm(0.975+alpha/2, corr=corr0)[1]
    if(i%%100==0)cat("..", i)
  }
  pval.vec[is.nan(pval.vec)]=0
  pval.trunc=max(pval.vec)+alpha
  if(pval.trunc > 1){
    pval.trunc=1
  }
  pval.max=max(max.dev.pval.vec)+alpha
  if(pval.max > 1){
    pval.max=1
  }
  
  pval.matrix[k,1:2]=c(pval.trunc, pval.max)
  
  ######### 1.2. from rpart
  chosen.check=rep(0,5)
  ### training sample -> Create tree using causaltree
  res=rpart((y_t-y_c) ~ x1+x2+x3+x4+x5, data=tra.matched.data, method="anova", control=rpart.control(cp=0))
  opt.cp=res$cptable[,1][which.min(res$cptable[,4])]
  pruned=prune(res, opt.cp) ## obtained tree from a training sample
  
  var.names=pruned$frame$var[pruned$frame$var!="<leaf>"]
  var.numbers=as.numeric(substr(var.names, start=2, stop=4))
  chosen.check.matrix.rpart=chosen.check
  chosen.check.matrix.rpart[var.numbers]=1
  
  ######## using rpart tree
  ### Test sample -> hypothesis test
  est.tree.stru=tree.stru(pruned, est.matched.data)
  subgroup.num=est.tree.stru$new.x$subgroup.num
  C=est.tree.stru$C
  unique.subg=as.numeric(colnames(est.tree.stru$C))
  num.of.subgroups=length(unique.subg)
  
  test.matched.pairs.ymat=est.matched.data
  test.matched.pairs.ymat$sub.num=subgroup.num
  
  ##########################################
  n.test=length(test.matched.pairs.ymat[,1])
  n.vec=table(test.matched.pairs.ymat$sub.num)
  
  # Since the population treatment effect size (mu) is not known, 
  ## get the 99% confidence interval of delta
  est.delta=sum(test.matched.pairs.ymat[,1]-test.matched.pairs.ymat[,2])/n.test
  alpha=0.01
  var.max=binary.hypo.test(Ymat=test.matched.pairs.ymat, delta0=est.delta)$max.var
  rough.int=c(est.delta-qnorm(1-alpha/2)*sqrt(var.max), est.delta+qnorm(1-alpha/2)*sqrt(var.max))
  #rough.int
  delta.end.points=round(rough.int*n*2)
  delta.seq=round(seq(delta.end.points[1], delta.end.points[2], length.out=20))
  
  ### subgroup estimates
  est.delta.vec=rep(NA, length(C[1,]))
  for(j in 1:length(C[1,])){
    est.delta.vec[j]=mean(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1]-test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],2])
  }
  #est.delta.vec
  
  library(mvtnorm)
  max.dev.vec=max.dev.pval.vec=pval.vec=rep(NA, length(delta.seq))
  for(i in 1:length(delta.seq)){
    temp.delta0=delta.seq[i]
    
    ### submax 
    temp.maxvar.vec=rep(NA, length(C[1,]))
    temp.delta0.vec=rep(NA, length(C[1,]))
    temp.pval.vec=rep(NA, length(C[1,]))
    for(j in 1:length(C[1,])){
      fraction.val=(temp.delta0/(2*n.test))*2*n.vec[j]
      closest.two.points=c(floor(fraction.val), ceiling(fraction.val))
      
      temp.hypo1=binary.hypo.test(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1:2], delta0=closest.two.points[1]/(2*n.vec[j]))
      temp.hypo2=binary.hypo.test(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1:2], delta0=closest.two.points[2]/(2*n.vec[j]))
      temp.maxvar1=temp.hypo1$max.var
      temp.maxvar2=temp.hypo2$max.var
      max.temp.maxvar=max(temp.maxvar1, temp.maxvar2)
      temp.maxvar.vec[j]=max.temp.maxvar
      temp.delta0.vec[j]=closest.two.points[which.max(c(temp.maxvar1, temp.maxvar2))]
      
      temp.pval.vec[j]=max(temp.hypo1$pvalue, temp.hypo2$pvalue)*2
    }
    
    T.vec=est.delta.vec*n.vec*2
    test= C %*% T.vec
    if(dim(C)[1]==1){
      var.mat=as.numeric(temp.maxvar.vec*(n.vec^2)*4)
      theta0= temp.delta0.vec
      om0= var.mat
      corr0= C
      deviate=(test-theta0)/sqrt(om0)
    }else{
      var.mat=diag(temp.maxvar.vec*(n.vec^2)*4)
      theta0= C %*% temp.delta0.vec
      om0=C %*% var.mat %*% t(C)
      deviate=(test-theta0)/sqrt(diag(om0))
      corr0=om0/sqrt(outer(diag(om0),diag(om0),"*"))
    }
    
    max.deviate=max(abs(deviate))
    
    max.dev.vec[i]=max.deviate
    if(dim(C)[1]==1){
      max.dev.pval.vec[i]=2*(1-pnorm(max.deviate))
    }else{
      max.dev.pval.vec[i]=2*(1-pmvnorm(lower=-Inf, upper=max.deviate, mean=rep(0, length(corr0[,1])), corr=corr0)[1])
    }
    pval.vec[i]=truncatedP(temp.pval.vec, trunc=0.1)
    #critical.constant.vec[i]=qmvnorm(0.975+alpha/2, corr=corr0)[1]
    if(i%%100==0)cat("..", i)
  }
  pval.vec[is.nan(pval.vec)]=0
  pval.trunc=max(pval.vec)+alpha
  if(pval.trunc > 1){
    pval.trunc=1
  }
  pval.max=max(max.dev.pval.vec)+alpha
  if(pval.max > 1){
    pval.max=1
  }
  
  pval.matrix[k,3:4]=c(pval.trunc, pval.max)
  check.matrix[k,1:10]=c(chosen.check.matrix.causaltree, chosen.check.matrix.rpart)
  
  ##################################################
  ## 2. Select training sample & Est sample (25% vs. 75%)
  ##################################################
  training.data.index=sample(1:n, n/4, replace=F)
  train=c(training.data.index, training.data.index+n)
  tra.matched.data=matched.data[training.data.index, ]
  est.matched.data=matched.data[-training.data.index,]
  tra.matched.data=as.data.frame(tra.matched.data)
  est.matched.data=as.data.frame(est.matched.data)
  
  traData=full.data[train,] # training data
  testData=full.data[-train, ] # test data
  
  ######### 2.1. Discovering tree structures from causal tree
  chosen.check=rep(0,5)
  ### training sample -> Create tree using causaltree
  tree=causalTree(y_obs ~ x1+x2+x3+x4+x5, data=traData, treatment = traData$z, HonestSampleSize=length(testData[,1]),
                  split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, xval = 5, 
                  cp = 0, minsize = 20, propensity = 0.5)
  
  opcp=tree$cptable[,1][which.min(tree$cptable[,4])]
  opfit=prune(tree, opcp) ## obtained tree from a training sample
  
  var.names=opfit$frame$var[opfit$frame$var!="<leaf>"]
  var.numbers=as.numeric(substr(var.names, start=2, stop=4))
  chosen.check.matrix.causaltree=chosen.check
  chosen.check.matrix.causaltree[var.numbers]=1
  
  
  ######## using causaltree
  ### Test sample -> hypothesis test
  est.tree.stru=tree.stru(opfit, est.matched.data)
  subgroup.num=est.tree.stru$new.x$subgroup.num
  C=est.tree.stru$C
  unique.subg=as.numeric(colnames(est.tree.stru$C))
  num.of.subgroups=length(unique.subg)
  
  test.matched.pairs.ymat=est.matched.data
  test.matched.pairs.ymat$sub.num=subgroup.num
  
  ##########################################
  n.test=length(test.matched.pairs.ymat[,1])
  n.vec=table(test.matched.pairs.ymat$sub.num)
  
  # Since the population treatment effect size (mu) is not known, 
  ## get the 99% confidence interval of delta
  est.delta=sum(test.matched.pairs.ymat[,1]-test.matched.pairs.ymat[,2])/n.test
  alpha=0.01
  var.max=binary.hypo.test(Ymat=test.matched.pairs.ymat, delta0=est.delta)$max.var
  rough.int=c(est.delta-qnorm(1-alpha/2)*sqrt(var.max), est.delta+qnorm(1-alpha/2)*sqrt(var.max))
  #rough.int
  delta.end.points=round(rough.int*n*2)
  delta.seq=round(seq(delta.end.points[1], delta.end.points[2], length.out=20))
  
  ### subgroup estimates
  est.delta.vec=rep(NA, length(C[1,]))
  for(j in 1:length(C[1,])){
    est.delta.vec[j]=mean(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1]-test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],2])
  }
  #est.delta.vec
  
  library(mvtnorm)
  max.dev.vec=max.dev.pval.vec=pval.vec=rep(NA, length(delta.seq))
  for(i in 1:length(delta.seq)){
    temp.delta0=delta.seq[i]
    
    ### submax 
    temp.maxvar.vec=rep(NA, length(C[1,]))
    temp.delta0.vec=rep(NA, length(C[1,]))
    temp.pval.vec=rep(NA, length(C[1,]))
    for(j in 1:length(C[1,])){
      fraction.val=(temp.delta0/(2*n.test))*2*n.vec[j]
      closest.two.points=c(floor(fraction.val), ceiling(fraction.val))
      
      temp.hypo1=binary.hypo.test(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1:2], delta0=closest.two.points[1]/(2*n.vec[j]))
      temp.hypo2=binary.hypo.test(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1:2], delta0=closest.two.points[2]/(2*n.vec[j]))
      temp.maxvar1=temp.hypo1$max.var
      temp.maxvar2=temp.hypo2$max.var
      max.temp.maxvar=max(temp.maxvar1, temp.maxvar2)
      temp.maxvar.vec[j]=max.temp.maxvar
      temp.delta0.vec[j]=closest.two.points[which.max(c(temp.maxvar1, temp.maxvar2))]
      
      temp.pval.vec[j]=max(temp.hypo1$pvalue, temp.hypo2$pvalue)*2
    }
    
    T.vec=est.delta.vec*n.vec*2
    test= C %*% T.vec
    if(dim(C)[1]==1){
      var.mat=as.numeric(temp.maxvar.vec*(n.vec^2)*4)
      theta0= temp.delta0.vec
      om0= var.mat
      corr0= C
      deviate=(test-theta0)/sqrt(om0)
    }else{
      var.mat=diag(temp.maxvar.vec*(n.vec^2)*4)
      theta0= C %*% temp.delta0.vec
      om0=C %*% var.mat %*% t(C)
      deviate=(test-theta0)/sqrt(diag(om0))
      corr0=om0/sqrt(outer(diag(om0),diag(om0),"*"))
    }
    
    max.deviate=max(abs(deviate))
    
    max.dev.vec[i]=max.deviate
    if(dim(C)[1]==1){
      max.dev.pval.vec[i]=2*(1-pnorm(max.deviate))
    }else{
      max.dev.pval.vec[i]=2*(1-pmvnorm(lower=-Inf, upper=max.deviate, mean=rep(0, length(corr0[,1])), corr=corr0)[1])
    }
    pval.vec[i]=truncatedP(temp.pval.vec, trunc=0.1)
    #critical.constant.vec[i]=qmvnorm(0.975+alpha/2, corr=corr0)[1]
    if(i%%100==0)cat("..", i)
  }
  pval.vec[is.nan(pval.vec)]=0
  pval.trunc=max(pval.vec)+alpha
  if(pval.trunc > 1){
    pval.trunc=1
  }
  pval.max=max(max.dev.pval.vec)+alpha
  if(pval.max > 1){
    pval.max=1
  }
  
  pval.matrix[k,5:6]=c(pval.trunc, pval.max)
  
  ######### 2.2. from rpart
  chosen.check=rep(0,5)
  ### training sample -> Create tree using causaltree
  res=rpart((y_t-y_c) ~ x1+x2+x3+x4+x5, data=tra.matched.data, method="anova", control=rpart.control(cp=0))
  opt.cp=res$cptable[,1][which.min(res$cptable[,4])]
  pruned=prune(res, opt.cp) ## obtained tree from a training sample
  
  var.names=pruned$frame$var[pruned$frame$var!="<leaf>"]
  var.numbers=as.numeric(substr(var.names, start=2, stop=4))
  chosen.check.matrix.rpart=chosen.check
  chosen.check.matrix.rpart[var.numbers]=1
  
  ######## using rpart tree
  ### Test sample -> hypothesis test
  est.tree.stru=tree.stru(pruned, est.matched.data)
  subgroup.num=est.tree.stru$new.x$subgroup.num
  C=est.tree.stru$C
  unique.subg=as.numeric(colnames(est.tree.stru$C))
  num.of.subgroups=length(unique.subg)
  
  test.matched.pairs.ymat=est.matched.data
  test.matched.pairs.ymat$sub.num=subgroup.num
  
  ##########################################
  n.test=length(test.matched.pairs.ymat[,1])
  n.vec=table(test.matched.pairs.ymat$sub.num)
  
  # Since the population treatment effect size (mu) is not known, 
  ## get the 99% confidence interval of delta
  est.delta=sum(test.matched.pairs.ymat[,1]-test.matched.pairs.ymat[,2])/n.test
  alpha=0.01
  var.max=binary.hypo.test(Ymat=test.matched.pairs.ymat, delta0=est.delta)$max.var
  rough.int=c(est.delta-qnorm(1-alpha/2)*sqrt(var.max), est.delta+qnorm(1-alpha/2)*sqrt(var.max))
  #rough.int
  delta.end.points=round(rough.int*n*2)
  delta.seq=round(seq(delta.end.points[1], delta.end.points[2], length.out=20))
  
  ### subgroup estimates
  est.delta.vec=rep(NA, length(C[1,]))
  for(j in 1:length(C[1,])){
    est.delta.vec[j]=mean(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1]-test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],2])
  }
  #est.delta.vec
  
  library(mvtnorm)
  max.dev.vec=max.dev.pval.vec=pval.vec=rep(NA, length(delta.seq))
  for(i in 1:length(delta.seq)){
    temp.delta0=delta.seq[i]
    
    ### submax 
    temp.maxvar.vec=rep(NA, length(C[1,]))
    temp.delta0.vec=rep(NA, length(C[1,]))
    temp.pval.vec=rep(NA, length(C[1,]))
    for(j in 1:length(C[1,])){
      fraction.val=(temp.delta0/(2*n.test))*2*n.vec[j]
      closest.two.points=c(floor(fraction.val), ceiling(fraction.val))
      
      temp.hypo1=binary.hypo.test(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1:2], delta0=closest.two.points[1]/(2*n.vec[j]))
      temp.hypo2=binary.hypo.test(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1:2], delta0=closest.two.points[2]/(2*n.vec[j]))
      temp.maxvar1=temp.hypo1$max.var
      temp.maxvar2=temp.hypo2$max.var
      max.temp.maxvar=max(temp.maxvar1, temp.maxvar2)
      temp.maxvar.vec[j]=max.temp.maxvar
      temp.delta0.vec[j]=closest.two.points[which.max(c(temp.maxvar1, temp.maxvar2))]
      
      temp.pval.vec[j]=max(temp.hypo1$pvalue, temp.hypo2$pvalue)*2
    }
    
    T.vec=est.delta.vec*n.vec*2
    test= C %*% T.vec
    if(dim(C)[1]==1){
      var.mat=as.numeric(temp.maxvar.vec*(n.vec^2)*4)
      theta0= temp.delta0.vec
      om0= var.mat
      corr0= C
      deviate=(test-theta0)/sqrt(om0)
    }else{
      var.mat=diag(temp.maxvar.vec*(n.vec^2)*4)
      theta0= C %*% temp.delta0.vec
      om0=C %*% var.mat %*% t(C)
      deviate=(test-theta0)/sqrt(diag(om0))
      corr0=om0/sqrt(outer(diag(om0),diag(om0),"*"))
    }
    
    max.deviate=max(abs(deviate))
    
    max.dev.vec[i]=max.deviate
    if(dim(C)[1]==1){
      max.dev.pval.vec[i]=2*(1-pnorm(max.deviate))
    }else{
      max.dev.pval.vec[i]=2*(1-pmvnorm(lower=-Inf, upper=max.deviate, mean=rep(0, length(corr0[,1])), corr=corr0)[1])
    }
    pval.vec[i]=truncatedP(temp.pval.vec, trunc=0.1)
    #critical.constant.vec[i]=qmvnorm(0.975+alpha/2, corr=corr0)[1]
    if(i%%100==0)cat("..", i)
  }
  pval.vec[is.nan(pval.vec)]=0
  pval.trunc=max(pval.vec)+alpha
  if(pval.trunc > 1){
    pval.trunc=1
  }
  pval.max=max(max.dev.pval.vec)+alpha
  if(pval.max > 1){
    pval.max=1
  }
  
  pval.matrix[k,7:8]=c(pval.trunc, pval.max)
  check.matrix[k,11:20]=c(chosen.check.matrix.causaltree, chosen.check.matrix.rpart)
  
  
  ##################################################
  ## 3. Select training sample & Est sample (50% vs. 50%)
  ##################################################
  training.data.index=sample(1:n, n/2, replace=F)
  train=c(training.data.index, training.data.index+n)
  tra.matched.data=matched.data[training.data.index, ]
  est.matched.data=matched.data[-training.data.index,]
  tra.matched.data=as.data.frame(tra.matched.data)
  est.matched.data=as.data.frame(est.matched.data)
  
  traData=full.data[train,] # training data
  testData=full.data[-train, ] # test data
  
  ######### 3.1. Discovering tree structures from causal tree
  chosen.check=rep(0,5)
  ### training sample -> Create tree using causaltree
  tree=causalTree(y_obs ~ x1+x2+x3+x4+x5, data=traData, treatment = traData$z, HonestSampleSize=length(testData[,1]),
                  split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, xval = 5, 
                  cp = 0, minsize = 20, propensity = 0.5)
  
  opcp=tree$cptable[,1][which.min(tree$cptable[,4])]
  opfit=prune(tree, opcp) ## obtained tree from a training sample
  
  var.names=opfit$frame$var[opfit$frame$var!="<leaf>"]
  var.numbers=as.numeric(substr(var.names, start=2, stop=4))
  chosen.check.matrix.causaltree=chosen.check
  chosen.check.matrix.causaltree[var.numbers]=1
  
  
  ######## using causaltree
  ### Test sample -> hypothesis test
  est.tree.stru=tree.stru(opfit, est.matched.data)
  subgroup.num=est.tree.stru$new.x$subgroup.num
  C=est.tree.stru$C
  unique.subg=as.numeric(colnames(est.tree.stru$C))
  num.of.subgroups=length(unique.subg)
  
  test.matched.pairs.ymat=est.matched.data
  test.matched.pairs.ymat$sub.num=subgroup.num
  
  ##########################################
  n.test=length(test.matched.pairs.ymat[,1])
  n.vec=table(test.matched.pairs.ymat$sub.num)
  
  # Since the population treatment effect size (mu) is not known, 
  ## get the 99% confidence interval of delta
  est.delta=sum(test.matched.pairs.ymat[,1]-test.matched.pairs.ymat[,2])/n.test
  alpha=0.01
  var.max=binary.hypo.test(Ymat=test.matched.pairs.ymat, delta0=est.delta)$max.var
  rough.int=c(est.delta-qnorm(1-alpha/2)*sqrt(var.max), est.delta+qnorm(1-alpha/2)*sqrt(var.max))
  #rough.int
  delta.end.points=round(rough.int*n*2)
  delta.seq=round(seq(delta.end.points[1], delta.end.points[2], length.out=20))
  
  ### subgroup estimates
  est.delta.vec=rep(NA, length(C[1,]))
  for(j in 1:length(C[1,])){
    est.delta.vec[j]=mean(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1]-test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],2])
  }
  #est.delta.vec
  
  library(mvtnorm)
  max.dev.vec=max.dev.pval.vec=pval.vec=rep(NA, length(delta.seq))
  for(i in 1:length(delta.seq)){
    temp.delta0=delta.seq[i]
    
    ### submax 
    temp.maxvar.vec=rep(NA, length(C[1,]))
    temp.delta0.vec=rep(NA, length(C[1,]))
    temp.pval.vec=rep(NA, length(C[1,]))
    for(j in 1:length(C[1,])){
      fraction.val=(temp.delta0/(2*n.test))*2*n.vec[j]
      closest.two.points=c(floor(fraction.val), ceiling(fraction.val))
      
      temp.hypo1=binary.hypo.test(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1:2], delta0=closest.two.points[1]/(2*n.vec[j]))
      temp.hypo2=binary.hypo.test(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1:2], delta0=closest.two.points[2]/(2*n.vec[j]))
      temp.maxvar1=temp.hypo1$max.var
      temp.maxvar2=temp.hypo2$max.var
      max.temp.maxvar=max(temp.maxvar1, temp.maxvar2)
      temp.maxvar.vec[j]=max.temp.maxvar
      temp.delta0.vec[j]=closest.two.points[which.max(c(temp.maxvar1, temp.maxvar2))]
      
      temp.pval.vec[j]=max(temp.hypo1$pvalue, temp.hypo2$pvalue)*2
    }
    
    T.vec=est.delta.vec*n.vec*2
    test= C %*% T.vec
    if(dim(C)[1]==1){
      var.mat=as.numeric(temp.maxvar.vec*(n.vec^2)*4)
      theta0= temp.delta0.vec
      om0= var.mat
      corr0= C
      deviate=(test-theta0)/sqrt(om0)
    }else{
      var.mat=diag(temp.maxvar.vec*(n.vec^2)*4)
      theta0= C %*% temp.delta0.vec
      om0=C %*% var.mat %*% t(C)
      deviate=(test-theta0)/sqrt(diag(om0))
      corr0=om0/sqrt(outer(diag(om0),diag(om0),"*"))
    }
    
    max.deviate=max(abs(deviate))
    
    max.dev.vec[i]=max.deviate
    if(dim(C)[1]==1){
      max.dev.pval.vec[i]=2*(1-pnorm(max.deviate))
    }else{
      max.dev.pval.vec[i]=2*(1-pmvnorm(lower=-Inf, upper=max.deviate, mean=rep(0, length(corr0[,1])), corr=corr0)[1])
    }
    pval.vec[i]=truncatedP(temp.pval.vec, trunc=0.1)
    #critical.constant.vec[i]=qmvnorm(0.975+alpha/2, corr=corr0)[1]
    if(i%%100==0)cat("..", i)
  }
  pval.vec[is.nan(pval.vec)]=0
  pval.trunc=max(pval.vec)+alpha
  if(pval.trunc > 1){
    pval.trunc=1
  }
  pval.max=max(max.dev.pval.vec)+alpha
  if(pval.max > 1){
    pval.max=1
  }
  
  pval.matrix[k,9:10]=c(pval.trunc, pval.max)
  
  ######### 3.2. from rpart
  chosen.check=rep(0,5)
  ### training sample -> Create tree using causaltree
  res=rpart((y_t-y_c) ~ x1+x2+x3+x4+x5, data=tra.matched.data, method="anova", control=rpart.control(cp=0))
  opt.cp=res$cptable[,1][which.min(res$cptable[,4])]
  pruned=prune(res, opt.cp) ## obtained tree from a training sample
  
  var.names=pruned$frame$var[pruned$frame$var!="<leaf>"]
  var.numbers=as.numeric(substr(var.names, start=2, stop=4))
  chosen.check.matrix.rpart=chosen.check
  chosen.check.matrix.rpart[var.numbers]=1
  
  ######## using rpart tree
  ### Test sample -> hypothesis test
  est.tree.stru=tree.stru(pruned, est.matched.data)
  subgroup.num=est.tree.stru$new.x$subgroup.num
  C=est.tree.stru$C
  unique.subg=as.numeric(colnames(est.tree.stru$C))
  num.of.subgroups=length(unique.subg)
  
  test.matched.pairs.ymat=est.matched.data
  test.matched.pairs.ymat$sub.num=subgroup.num
  
  ##########################################
  n.test=length(test.matched.pairs.ymat[,1])
  n.vec=table(test.matched.pairs.ymat$sub.num)
  
  # Since the population treatment effect size (mu) is not known, 
  ## get the 99% confidence interval of delta
  est.delta=sum(test.matched.pairs.ymat[,1]-test.matched.pairs.ymat[,2])/n.test
  alpha=0.01
  var.max=binary.hypo.test(Ymat=test.matched.pairs.ymat, delta0=est.delta)$max.var
  rough.int=c(est.delta-qnorm(1-alpha/2)*sqrt(var.max), est.delta+qnorm(1-alpha/2)*sqrt(var.max))
  #rough.int
  delta.end.points=round(rough.int*n*2)
  delta.seq=round(seq(delta.end.points[1], delta.end.points[2], length.out=20))
  
  ### subgroup estimates
  est.delta.vec=rep(NA, length(C[1,]))
  for(j in 1:length(C[1,])){
    est.delta.vec[j]=mean(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1]-test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],2])
  }
  #est.delta.vec
  
  library(mvtnorm)
  max.dev.vec=max.dev.pval.vec=pval.vec=rep(NA, length(delta.seq))
  for(i in 1:length(delta.seq)){
    temp.delta0=delta.seq[i]
    
    ### submax 
    temp.maxvar.vec=rep(NA, length(C[1,]))
    temp.delta0.vec=rep(NA, length(C[1,]))
    temp.pval.vec=rep(NA, length(C[1,]))
    for(j in 1:length(C[1,])){
      fraction.val=(temp.delta0/(2*n.test))*2*n.vec[j]
      closest.two.points=c(floor(fraction.val), ceiling(fraction.val))
      
      temp.hypo1=binary.hypo.test(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1:2], delta0=closest.two.points[1]/(2*n.vec[j]))
      temp.hypo2=binary.hypo.test(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1:2], delta0=closest.two.points[2]/(2*n.vec[j]))
      temp.maxvar1=temp.hypo1$max.var
      temp.maxvar2=temp.hypo2$max.var
      max.temp.maxvar=max(temp.maxvar1, temp.maxvar2)
      temp.maxvar.vec[j]=max.temp.maxvar
      temp.delta0.vec[j]=closest.two.points[which.max(c(temp.maxvar1, temp.maxvar2))]
      
      temp.pval.vec[j]=max(temp.hypo1$pvalue, temp.hypo2$pvalue)*2
    }
    
    T.vec=est.delta.vec*n.vec*2
    test= C %*% T.vec
    if(dim(C)[1]==1){
      var.mat=as.numeric(temp.maxvar.vec*(n.vec^2)*4)
      theta0= temp.delta0.vec
      om0= var.mat
      corr0= C
      deviate=(test-theta0)/sqrt(om0)
    }else{
      var.mat=diag(temp.maxvar.vec*(n.vec^2)*4)
      theta0= C %*% temp.delta0.vec
      om0=C %*% var.mat %*% t(C)
      deviate=(test-theta0)/sqrt(diag(om0))
      corr0=om0/sqrt(outer(diag(om0),diag(om0),"*"))
    }
    
    max.deviate=max(abs(deviate))
    
    max.dev.vec[i]=max.deviate
    if(dim(C)[1]==1){
      max.dev.pval.vec[i]=2*(1-pnorm(max.deviate))
    }else{
      max.dev.pval.vec[i]=2*(1-pmvnorm(lower=-Inf, upper=max.deviate, mean=rep(0, length(corr0[,1])), corr=corr0)[1])
    }
    pval.vec[i]=truncatedP(temp.pval.vec, trunc=0.1)
    #critical.constant.vec[i]=qmvnorm(0.975+alpha/2, corr=corr0)[1]
    if(i%%100==0)cat("..", i)
  }
  pval.vec[is.nan(pval.vec)]=0
  pval.trunc=max(pval.vec)+alpha
  if(pval.trunc > 1){
    pval.trunc=1
  }
  pval.max=max(max.dev.pval.vec)+alpha
  if(pval.max > 1){
    pval.max=1
  }
  
  pval.matrix[k,11:12]=c(pval.trunc, pval.max)
  check.matrix[k,21:30]=c(chosen.check.matrix.causaltree, chosen.check.matrix.rpart)
  
  print(k)
}
time.after=Sys.time()
time.after-time.before

pval.matrix
check.matrix

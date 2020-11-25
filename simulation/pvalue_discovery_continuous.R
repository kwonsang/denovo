
#################################################################
## Compared to the Causal Tree method (Athey and Imbens, PNAS)
#################################################################
library(causalTree)
library(mvtnorm)
library(rpart.utils)
library(sensitivitymv)

source("../Rfunctions/basic_functions.R")


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
  y_t[x1==0]=rnorm(sum(x1==0), mean=0.4, sd=sd)
  y_t[x1==1]=rnorm(sum(x1==1), mean=0.6, sd=sd)
  
  y_c=rnorm(n, mean=0, sd=sd)
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
  est.tree.stru=tree.stru(opfit, testData)
  subgroup.num=est.tree.stru$new.x$subgroup.num
  unique.subg=as.numeric(colnames(est.tree.stru$C))
  num.of.subgroups=length(unique.subg)
  
  # Since the population treatment effect size (mu) is not known, we estimate the 95\% CI of mu. 
  tau.vec=get.tau.vector(testData[,1], testData[,2])
  
  ##### Truncated Product method
  pval.tau.vec=rep(NA, length(tau.vec))
  for(j in 1:length(tau.vec)){
    temp.tau=tau.vec[j]
    
    ## compute the p-value for each subgroup
    temp.pval.vec=rep(NA, num.of.subgroups)
    for(i in 1:num.of.subgroups){
      subgroup=testData[subgroup.num==unique.subg[i],]
      treated.y=subgroup[1:(length(subgroup[,1])/2),1]
      control.y=subgroup[-(1:(length(subgroup[,1])/2)),1]
      temp.pval.vec[i]=wilcox.test(treated.y-temp.tau, control.y, paired=T, alternative="two.sided")$p.value
      
    }
    pval.tau.vec[j]=truncatedP(temp.pval.vec, trunc=0.1)
  }
  pval.trunc=max(pval.tau.vec)+0.001
  if(pval.trunc > 1){
    pval.trunc=1
  }
  
  
  ##### Submax
  C=est.tree.stru$C
  
  null.values=matrix(NA, nrow=num.of.subgroups, ncol=2)
  for(i in 1:num.of.subgroups){
    ## null values
    null.values[i,]=as.vector(unlist(wilcoxSenMoments(sum(subgroup.num==unique.subg[i])/2, gamma=1)))[c(1,3)]
  }
  mu=null.values[,1]; V=null.values[,2];
  if(dim(C)[1]==1){
    theta0= mu
    sigma0= V
    corr0= C
  }else{
    theta0= C %*% mu
    sigma0= C %*% diag(V) %*% t(C)
    corr0=sigma0/sqrt(outer(diag(sigma0),diag(sigma0),"*"))
  }
  
  max.dev.vec=max.dev.pval.vec=rep(NA, length(tau.vec))
  for(j in 1:length(tau.vec)){
    temp.tau=tau.vec[j]
    
    ## compute the test statistic for each subgroup
    temp.t.vec=rep(NA, num.of.subgroups)
    null.values=matrix(NA, nrow=num.of.subgroups, ncol=2)
    for(i in 1:num.of.subgroups){
      subgroup=testData[subgroup.num==unique.subg[i],]
      treated.y=subgroup[1:(length(subgroup[,1])/2),1]
      control.y=subgroup[-(1:(length(subgroup[,1])/2)),1]
      temp.t.vec[i]=wilcox.test(treated.y-temp.tau, control.y, paired=T, alternative="two.sided")$statistic
    }
    
    temp.test= C %*% temp.t.vec
    if(dim(C)[1]==1){
      temp.deviate=(temp.test-theta0)/sqrt(sigma0)
    }else{
      temp.deviate=(temp.test-theta0)/sqrt(diag(sigma0))
    }
    
    max.deviate=max(abs(temp.deviate))
    max.dev.vec[j]=max.deviate
    if(dim(C)[1]==1){
      max.dev.pval.vec[j]=2*(1-pnorm(max.deviate))
    }else{
      max.dev.pval.vec[j]=2*(1-pmvnorm(lower=-Inf, upper=max.deviate, mean=rep(0, length(corr0[,1])), corr=corr0)[1])
    }
    
  }
  pval.max=max(max.dev.pval.vec)+0.001
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
  unique.subg=as.numeric(colnames(est.tree.stru$C))
  num.of.subgroups=length(unique.subg)
  
  # Since the population treatment effect size (mu) is not known, we estimate the 95\% CI of mu. 
  tau.vec=get.tau.vector(testData[,1], testData[,2])
  ##### Truncated Product method
  pval.tau.vec=rep(NA, length(tau.vec))
  for(j in 1:length(tau.vec)){
    temp.tau=tau.vec[j]
    
    ## compute the p-value for each subgroup
    temp.pval.vec=rep(NA, num.of.subgroups)
    for(i in 1:num.of.subgroups){
      subgroup=testData[subgroup.num==unique.subg[i],]
      treated.y=subgroup[1:(length(subgroup[,1])/2),1]
      control.y=subgroup[-(1:(length(subgroup[,1])/2)),1]
      temp.pval.vec[i]=wilcox.test(treated.y-temp.tau, control.y, paired=T, alternative="two.sided")$p.value
      
    }
    pval.tau.vec[j]=truncatedP(temp.pval.vec, trunc=0.1)
  }
  pval.trunc=max(pval.tau.vec)+0.001
  if(pval.trunc > 1){
    pval.trunc=1
  }
  
  
  ##### Submax
  C=est.tree.stru$C
  
  null.values=matrix(NA, nrow=num.of.subgroups, ncol=2)
  for(i in 1:num.of.subgroups){
    ## null values
    null.values[i,]=as.vector(unlist(wilcoxSenMoments(sum(subgroup.num==unique.subg[i]), gamma=1)))[c(1,3)]
  }
  mu=null.values[,1]; V=null.values[,2];
  if(dim(C)[1]==1){
    theta0= mu
    sigma0= V
    corr0= C
  }else{
    theta0= C %*% mu
    sigma0= C %*% diag(V) %*% t(C)
    corr0=sigma0/sqrt(outer(diag(sigma0),diag(sigma0),"*"))
  }
  
  
  max.dev.vec=max.dev.pval.vec=rep(NA, length(tau.vec))
  for(j in 1:length(tau.vec)){
    temp.tau=tau.vec[j]
    
    ## compute the test statistic for each subgroup
    temp.t.vec=rep(NA, num.of.subgroups)
    null.values=matrix(NA, nrow=num.of.subgroups, ncol=2)
    for(i in 1:num.of.subgroups){
      subgroup=testData[subgroup.num==unique.subg[i],]
      treated.y=subgroup[1:(length(subgroup[,1])/2),1]
      control.y=subgroup[-(1:(length(subgroup[,1])/2)),1]
      temp.t.vec[i]=wilcox.test(treated.y-temp.tau, control.y, paired=T, alternative="two.sided")$statistic
    }
    
    temp.test= C %*% temp.t.vec
    if(dim(C)[1]==1){
      temp.deviate=(temp.test-theta0)/sqrt(sigma0)
    }else{
      temp.deviate=(temp.test-theta0)/sqrt(diag(sigma0))
    }
    
    max.deviate=max(abs(temp.deviate))
    max.dev.vec[j]=max.deviate
    if(dim(C)[1]==1){
      max.dev.pval.vec[j]=2*(1-pnorm(max.deviate))
    }else{
      max.dev.pval.vec[j]=2*(1-pmvnorm(lower=-Inf, upper=max.deviate, mean=rep(0, length(corr0[,1])), corr=corr0)[1])
    }
  }
  pval.max=max(max.dev.pval.vec)+0.001
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
  est.tree.stru=tree.stru(opfit, testData)
  subgroup.num=est.tree.stru$new.x$subgroup.num
  unique.subg=as.numeric(colnames(est.tree.stru$C))
  num.of.subgroups=length(unique.subg)
  
  # Since the population treatment effect size (mu) is not known, we estimate the 95\% CI of mu. 
  tau.vec=get.tau.vector(testData[,1], testData[,2])
  
  ##### Truncated Product method
  pval.tau.vec=rep(NA, length(tau.vec))
  for(j in 1:length(tau.vec)){
    temp.tau=tau.vec[j]
    
    ## compute the p-value for each subgroup
    temp.pval.vec=rep(NA, num.of.subgroups)
    for(i in 1:num.of.subgroups){
      subgroup=testData[subgroup.num==unique.subg[i],]
      treated.y=subgroup[1:(length(subgroup[,1])/2),1]
      control.y=subgroup[-(1:(length(subgroup[,1])/2)),1]
      temp.pval.vec[i]=wilcox.test(treated.y-temp.tau, control.y, paired=T, alternative="two.sided")$p.value
      
    }
    pval.tau.vec[j]=truncatedP(temp.pval.vec, trunc=0.1)
  }
  pval.trunc=max(pval.tau.vec)+0.001
  if(pval.trunc > 1){
    pval.trunc=1
  }
  
  
  ##### Submax
  C=est.tree.stru$C
  
  null.values=matrix(NA, nrow=num.of.subgroups, ncol=2)
  for(i in 1:num.of.subgroups){
    ## null values
    null.values[i,]=as.vector(unlist(wilcoxSenMoments(sum(subgroup.num==unique.subg[i])/2, gamma=1)))[c(1,3)]
  }
  mu=null.values[,1]; V=null.values[,2];
  if(dim(C)[1]==1){
    theta0= mu
    sigma0= V
    corr0= C
  }else{
    theta0= C %*% mu
    sigma0= C %*% diag(V) %*% t(C)
    corr0=sigma0/sqrt(outer(diag(sigma0),diag(sigma0),"*"))
  }
  
  max.dev.vec=max.dev.pval.vec=rep(NA, length(tau.vec))
  for(j in 1:length(tau.vec)){
    temp.tau=tau.vec[j]
    
    ## compute the test statistic for each subgroup
    temp.t.vec=rep(NA, num.of.subgroups)
    null.values=matrix(NA, nrow=num.of.subgroups, ncol=2)
    for(i in 1:num.of.subgroups){
      subgroup=testData[subgroup.num==unique.subg[i],]
      treated.y=subgroup[1:(length(subgroup[,1])/2),1]
      control.y=subgroup[-(1:(length(subgroup[,1])/2)),1]
      temp.t.vec[i]=wilcox.test(treated.y-temp.tau, control.y, paired=T, alternative="two.sided")$statistic
    }
    
    temp.test= C %*% temp.t.vec
    if(dim(C)[1]==1){
      temp.deviate=(temp.test-theta0)/sqrt(sigma0)
    }else{
      temp.deviate=(temp.test-theta0)/sqrt(diag(sigma0))
    }
    
    max.deviate=max(abs(temp.deviate))
    max.dev.vec[j]=max.deviate
    if(dim(C)[1]==1){
      max.dev.pval.vec[j]=2*(1-pnorm(max.deviate))
    }else{
      max.dev.pval.vec[j]=2*(1-pmvnorm(lower=-Inf, upper=max.deviate, mean=rep(0, length(corr0[,1])), corr=corr0)[1])
    }
  }
  pval.max=max(max.dev.pval.vec)+0.001
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
  unique.subg=as.numeric(colnames(est.tree.stru$C))
  num.of.subgroups=length(unique.subg)
  
  # Since the population treatment effect size (mu) is not known, we estimate the 95\% CI of mu. 
  tau.vec=get.tau.vector(testData[,1], testData[,2])
  
  ##### Truncated Product method
  pval.tau.vec=rep(NA, length(tau.vec))
  for(j in 1:length(tau.vec)){
    temp.tau=tau.vec[j]
    
    ## compute the p-value for each subgroup
    temp.pval.vec=rep(NA, num.of.subgroups)
    for(i in 1:num.of.subgroups){
      subgroup=testData[subgroup.num==unique.subg[i],]
      treated.y=subgroup[1:(length(subgroup[,1])/2),1]
      control.y=subgroup[-(1:(length(subgroup[,1])/2)),1]
      temp.pval.vec[i]=wilcox.test(treated.y-temp.tau, control.y, paired=T, alternative="two.sided")$p.value
      
    }
    pval.tau.vec[j]=truncatedP(temp.pval.vec, trunc=0.1)
  }
  pval.trunc=max(pval.tau.vec)+0.001
  if(pval.trunc > 1){
    pval.trunc=1
  }
  
  
  ##### Submax
  C=est.tree.stru$C
  
  null.values=matrix(NA, nrow=num.of.subgroups, ncol=2)
  for(i in 1:num.of.subgroups){
    ## null values
    null.values[i,]=as.vector(unlist(wilcoxSenMoments(sum(subgroup.num==unique.subg[i]), gamma=1)))[c(1,3)]
  }
  mu=null.values[,1]; V=null.values[,2];
  if(dim(C)[1]==1){
    theta0= mu
    sigma0= V
    corr0= C
  }else{
    theta0= C %*% mu
    sigma0= C %*% diag(V) %*% t(C)
    corr0=sigma0/sqrt(outer(diag(sigma0),diag(sigma0),"*"))
  }
  
  max.dev.vec=max.dev.pval.vec=rep(NA, length(tau.vec))
  for(j in 1:length(tau.vec)){
    temp.tau=tau.vec[j]
    
    ## compute the test statistic for each subgroup
    temp.t.vec=rep(NA, num.of.subgroups)
    null.values=matrix(NA, nrow=num.of.subgroups, ncol=2)
    for(i in 1:num.of.subgroups){
      subgroup=testData[subgroup.num==unique.subg[i],]
      treated.y=subgroup[1:(length(subgroup[,1])/2),1]
      control.y=subgroup[-(1:(length(subgroup[,1])/2)),1]
      temp.t.vec[i]=wilcox.test(treated.y-temp.tau, control.y, paired=T, alternative="two.sided")$statistic
    }
    
    temp.test= C %*% temp.t.vec
    if(dim(C)[1]==1){
      temp.deviate=(temp.test-theta0)/sqrt(sigma0)
    }else{
      temp.deviate=(temp.test-theta0)/sqrt(diag(sigma0))
    }
    
    max.deviate=max(abs(temp.deviate))
    max.dev.vec[j]=max.deviate
    if(dim(C)[1]==1){
      max.dev.pval.vec[j]=2*(1-pnorm(max.deviate))
    }else{
      max.dev.pval.vec[j]=2*(1-pmvnorm(lower=-Inf, upper=max.deviate, mean=rep(0, length(corr0[,1])), corr=corr0)[1])
    }
  }
  pval.max=max(max.dev.pval.vec)+0.001
  if(pval.max > 1){
    pval.max=1
  }
  
  pval.matrix[k,7:8]=c(pval.trunc, pval.max)
  check.matrix[k,11:20]=c(chosen.check.matrix.causaltree, chosen.check.matrix.rpart)
  
  ##################################################
  ## 3. Select training sample & Est sample (25% vs. 75%)
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
  est.tree.stru=tree.stru(opfit, testData)
  subgroup.num=est.tree.stru$new.x$subgroup.num
  unique.subg=as.numeric(colnames(est.tree.stru$C))
  num.of.subgroups=length(unique.subg)
  
  # Since the population treatment effect size (mu) is not known, we estimate the 95\% CI of mu. 
  tau.vec=get.tau.vector(testData[,1], testData[,2])
  
  ##### Truncated Product method
  pval.tau.vec=rep(NA, length(tau.vec))
  for(j in 1:length(tau.vec)){
    temp.tau=tau.vec[j]
    
    ## compute the p-value for each subgroup
    temp.pval.vec=rep(NA, num.of.subgroups)
    for(i in 1:num.of.subgroups){
      subgroup=testData[subgroup.num==unique.subg[i],]
      treated.y=subgroup[1:(length(subgroup[,1])/2),1]
      control.y=subgroup[-(1:(length(subgroup[,1])/2)),1]
      temp.pval.vec[i]=wilcox.test(treated.y-temp.tau, control.y, paired=T, alternative="two.sided")$p.value
      
    }
    pval.tau.vec[j]=truncatedP(temp.pval.vec, trunc=0.1)
  }
  pval.trunc=max(pval.tau.vec)+0.001
  if(pval.trunc > 1){
    pval.trunc=1
  }
  
  
  ##### Submax
  C=est.tree.stru$C
  
  null.values=matrix(NA, nrow=num.of.subgroups, ncol=2)
  for(i in 1:num.of.subgroups){
    ## null values
    null.values[i,]=as.vector(unlist(wilcoxSenMoments(sum(subgroup.num==unique.subg[i])/2, gamma=1)))[c(1,3)]
  }
  mu=null.values[,1]; V=null.values[,2];
  if(dim(C)[1]==1){
    theta0= mu
    sigma0= V
    corr0= C
  }else{
    theta0= C %*% mu
    sigma0= C %*% diag(V) %*% t(C)
    corr0=sigma0/sqrt(outer(diag(sigma0),diag(sigma0),"*"))
  }
  
  max.dev.vec=max.dev.pval.vec=rep(NA, length(tau.vec))
  for(j in 1:length(tau.vec)){
    temp.tau=tau.vec[j]
    
    ## compute the test statistic for each subgroup
    temp.t.vec=rep(NA, num.of.subgroups)
    null.values=matrix(NA, nrow=num.of.subgroups, ncol=2)
    for(i in 1:num.of.subgroups){
      subgroup=testData[subgroup.num==unique.subg[i],]
      treated.y=subgroup[1:(length(subgroup[,1])/2),1]
      control.y=subgroup[-(1:(length(subgroup[,1])/2)),1]
      temp.t.vec[i]=wilcox.test(treated.y-temp.tau, control.y, paired=T, alternative="two.sided")$statistic
    }
    
    temp.test= C %*% temp.t.vec
    if(dim(C)[1]==1){
      temp.deviate=(temp.test-theta0)/sqrt(sigma0)
    }else{
      temp.deviate=(temp.test-theta0)/sqrt(diag(sigma0))
    }
    
    max.deviate=max(abs(temp.deviate))
    max.dev.vec[j]=max.deviate
    if(dim(C)[1]==1){
      max.dev.pval.vec[j]=2*(1-pnorm(max.deviate))
    }else{
      max.dev.pval.vec[j]=2*(1-pmvnorm(lower=-Inf, upper=max.deviate, mean=rep(0, length(corr0[,1])), corr=corr0)[1])
    }
  }
  pval.max=max(max.dev.pval.vec)+0.001
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
  unique.subg=as.numeric(colnames(est.tree.stru$C))
  num.of.subgroups=length(unique.subg)
  
  # Since the population treatment effect size (mu) is not known, we estimate the 95\% CI of mu. 
  tau.vec=get.tau.vector(testData[,1], testData[,2])
  
  ##### Truncated Product method
  pval.tau.vec=rep(NA, length(tau.vec))
  for(j in 1:length(tau.vec)){
    temp.tau=tau.vec[j]
    
    ## compute the p-value for each subgroup
    temp.pval.vec=rep(NA, num.of.subgroups)
    for(i in 1:num.of.subgroups){
      subgroup=testData[subgroup.num==unique.subg[i],]
      treated.y=subgroup[1:(length(subgroup[,1])/2),1]
      control.y=subgroup[-(1:(length(subgroup[,1])/2)),1]
      temp.pval.vec[i]=wilcox.test(treated.y-temp.tau, control.y, paired=T, alternative="two.sided")$p.value
      
    }
    pval.tau.vec[j]=truncatedP(temp.pval.vec, trunc=0.1)
  }
  pval.trunc=max(pval.tau.vec)+0.001
  if(pval.trunc > 1){
    pval.trunc=1
  }
  
  
  ##### Submax
  C=est.tree.stru$C
  
  null.values=matrix(NA, nrow=num.of.subgroups, ncol=2)
  for(i in 1:num.of.subgroups){
    ## null values
    null.values[i,]=as.vector(unlist(wilcoxSenMoments(sum(subgroup.num==unique.subg[i]), gamma=1)))[c(1,3)]
  }
  mu=null.values[,1]; V=null.values[,2];
  if(dim(C)[1]==1){
    theta0= mu
    sigma0= V
    corr0= C
  }else{
    theta0= C %*% mu
    sigma0= C %*% diag(V) %*% t(C)
    corr0=sigma0/sqrt(outer(diag(sigma0),diag(sigma0),"*"))
  }
  
  max.dev.vec=max.dev.pval.vec=rep(NA, length(tau.vec))
  for(j in 1:length(tau.vec)){
    temp.tau=tau.vec[j]
    
    ## compute the test statistic for each subgroup
    temp.t.vec=rep(NA, num.of.subgroups)
    null.values=matrix(NA, nrow=num.of.subgroups, ncol=2)
    for(i in 1:num.of.subgroups){
      subgroup=testData[subgroup.num==unique.subg[i],]
      treated.y=subgroup[1:(length(subgroup[,1])/2),1]
      control.y=subgroup[-(1:(length(subgroup[,1])/2)),1]
      temp.t.vec[i]=wilcox.test(treated.y-temp.tau, control.y, paired=T, alternative="two.sided")$statistic
    }
    
    temp.test= C %*% temp.t.vec
    if(dim(C)[1]==1){
      temp.deviate=(temp.test-theta0)/sqrt(sigma0)
    }else{
      temp.deviate=(temp.test-theta0)/sqrt(diag(sigma0))
    }
    
    max.deviate=max(abs(temp.deviate))
    max.dev.vec[j]=max.deviate
    if(dim(C)[1]==1){
      max.dev.pval.vec[j]=2*(1-pnorm(max.deviate))
    }else{
      max.dev.pval.vec[j]=2*(1-pmvnorm(lower=-Inf, upper=max.deviate, mean=rep(0, length(corr0[,1])), corr=corr0)[1])
    }
  }
  pval.max=max(max.dev.pval.vec)+0.001
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
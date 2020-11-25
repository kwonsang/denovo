#' Find an effect modification structure from data.
#'
#' @param Data input matched pairs. The first column must contain outcomes for treated subjects and the second column must contain outcomes for control subjects. The remaining columns contain shared covariates.
#' @param split.ratio set a splitting ratio. 0.25 means using 25\% of data as a training sample.
#' @return A list with the elements:
#' \item{training.index}{an index set for indicating a training sample that is randomly chosen.}
#' \item{test.index}{an index set for indicating a test sample.}
#' \item{tree}{discovered effect modification tree structure.}
#' @import rpart
#' @export
denovo.training=function(Data, split.ratio=0.25){
  require(rpart)
  
  n=dim(Data)[1] # number of pairs
  m=dim(Data)[2]
  
  colnames(Data)[1:2]=c("Y_t", "Y_c")
  
  training.data.index=sample(1:n, round(n*split.ratio), replace=F)
  test.data.index=(1:n)[-training.data.index]
  
  tra.matched.data=Data[training.data.index, ]
  tra.matched.data=as.data.frame(tra.matched.data)
  
  ################################################################################################
  ######### Discovering tree structures in the first subsample (training sample)
  ################################################################################################
  
  ### training sample -> Create tree using CART
  res=rpart((Y_t-Y_c)~., data=tra.matched.data, method="anova", control=rpart.control(cp=0))
  opt.cp=res$cptable[,1][which.min(res$cptable[,4])]
  pruned=prune(res, opt.cp) ## obtained tree from a training sample
  return(list(training.index=training.data.index, test.index=test.data.index, tree=pruned))
}

#' Conduct a hypothesis test under no unmeasured confounder assumption.
#'
#' @param Data a total sample of matched pairs.
#' @param test.index an index set that indicates a test sample in the total sample.
#' @param tree an effect modification tree that needs to be tested.
#' @param total.significance a total significance level.
#' @param gamma a significance level that is used for estimating the confidence interval.
#' @return A result matrix with deviates is reported. The `Max` column represents the maximum value of deviates.
#' The `kappa` column represents the critical value.
#' @import mvtnorm
#' @export
denovo.test=function(Data, test.index, tree, total.significance=0.05, gamma=0.01){
  require(mvtnorm)
  test.matched.data=Data[test.index,]
  test.matched.data=as.data.frame(test.matched.data)
  
  test.tree.stru=tree.stru(tree=tree, x=test.matched.data)
  
  subgroup.num=test.tree.stru$new.x$subgroup.num
  unique.subg=as.numeric(colnames(test.tree.stru$C))
  num.of.subgroups=length(unique.subg)
  
  ########################################################################################################
  ##### Use the second sample (test sample)
  ##### test of effect modification when Gamma=1 (under no unmeasured confounder assumption)
  ########################################################################################################
  C=test.tree.stru$C
  
  # set significance levels
  alpha=total.significance-gamma
  
  # Since the population treatment effect size (tau) is not known, we estimate the 100(1-gamma)% CI of tau.
  tau.vec=get.tau.vector.wilcox(test.matched.data, gamma=gamma)
  
  null.values=matrix(NA, nrow=num.of.subgroups, ncol=3)
  for(i in 1:num.of.subgroups){
    ## null values
    null.values[i,]=as.vector(unlist(wilcoxSenMoments(sum(subgroup.num==unique.subg[i]), gamma=1)))
  }
  mu=null.values[,1]; V=null.values[,3];
  if(dim(C)[1]==1){
    theta0= mu
    sigma0= V
    corr0= C
  }else{
    theta0= C %*% mu
    sigma0= C %*% diag(V) %*% t(C)
    corr0=sigma0/sqrt(outer(diag(sigma0),diag(sigma0),"*"))
  }
  
  max.dev.vec=rep(NA, length(tau.vec))
  dev.mat=matrix(NA, nrow=length(tau.vec), ncol=dim(C)[1])
  for(j in 1:length(tau.vec)){
    temp.tau=tau.vec[j]
    
    ## compute the test statistic T for each subgroup
    temp.t.vec=rep(NA, num.of.subgroups)
    for(i in 1:num.of.subgroups){
      subgroup=test.matched.data[subgroup.num==unique.subg[i],]
      treated.y=subgroup[,1]
      control.y=subgroup[,2]
      temp.t.vec[i]=wilcox.test(treated.y-temp.tau, control.y, paired=T, alternative="two.sided")$statistic
    }
    
    ## compute the comparisons using S=CT
    temp.test= C %*% temp.t.vec
    if(dim(C)[1]==1){
      temp.deviate=(temp.test-theta0)/sqrt(sigma0)
    }else{
      temp.deviate=(temp.test-theta0)/sqrt(diag(sigma0))
    }
    dev.mat[j,]=temp.deviate
    max.deviate=max(abs(temp.deviate))
    max.dev.vec[j]=max.deviate
  }
  ## Critical value at alpha.
  if(dim(C)[1]==1){
    critical.val=qnorm(1-alpha/2)
  }else{
    critical.val=qmvnorm(1-alpha/2, mean=rep(0, length(corr0[,1])), corr=corr0)$quantile
  }
  
  dev.res.mat=cbind(dev.mat, max.dev.vec, rep(critical.val, length(max.dev.vec))) #
  colnames(dev.res.mat)=c(rownames(C), "Max", "kappa")
  dev.res.mat=as.data.frame(dev.res.mat)
  dev.res.mat$tau=tau.vec
  return(dev.res.mat)
}

#' Conduct a sensitivity analysis
#'
#' @param Data a total sample of matched pairs.
#' @param test.index an index set that indicates a test sample in the total sample.
#' @param tree an effect modification tree that needs to be tested.
#' @param Gamma.vec a vector of sensitivity parameter \eqn{\Gamma}.
#' @param total.significance a total significance level.
#' @param gamma a significance level that is used for estimating the confidence interval.
#' @return A result matrix with deviates is reported. The `Max` column represents the maximum value of deviates.
#' The `kappa` column represents the critical value. If `Max` is greater than 'kappa', we reject the null hypothesis.
#' @import mvtnorm
#' @export
denovo.test.sensi=function(Data, test.index, tree, Gamma.vec, total.significance=0.05, gamma=0.0001){
  require(mvtnorm)
  test.matched.data=Data[test.index,]
  test.matched.data=as.data.frame(test.matched.data)
  
  test.tree.stru=tree.stru(tree=tree, x=test.matched.data)
  
  subgroup.num=test.tree.stru$new.x$subgroup.num
  unique.subg=as.numeric(colnames(test.tree.stru$C))
  num.of.subgroups=length(unique.subg)
  C=test.tree.stru$C
  
  ########################################################################
  ########## Sensitivity analysis with various values of Gamma. ##########
  ########################################################################
  sensi.param.vec=Gamma.vec # need to specify this vector.
  
  alpha=total.significance-gamma
  # Instead of estimating the 100(1-gamma)% CI of tau, choose a wide enough interval for tau
  tau.vec=get.tau.vector.wilcox(test.matched.data, gamma=gamma)
  
  sensi.mat=matrix(NA, nrow=length(sensi.param.vec), ncol=dim(C)[1]+2)
  for(k in 1:length(sensi.param.vec)){
    Gamma=sensi.param.vec[k]
    
    null.values=matrix(NA, nrow=num.of.subgroups, ncol=3)
    for(i in 1:num.of.subgroups){
      ## null values
      null.values[i,]=as.vector(unlist(wilcoxSenMoments(sum(subgroup.num==unique.subg[i]), gamma=Gamma)))
    }
    mu.upper=null.values[,1] # lower bound of Exp.
    mu.lower=null.values[,2] # upper bound of Exp.
    V=null.values[,3]
    if(dim(C)[1]==1){
      theta0.upper= mu.upper
      theta0.lower= mu.lower
      sigma0= V
      corr0= C
    }else{
      theta0.upper= C %*% mu.upper
      theta0.lower= C %*% mu.lower
      sigma0= C %*% diag(V) %*% t(C)
      corr0=sigma0/sqrt(outer(diag(sigma0),diag(sigma0),"*"))
    }
    
    max.dev.vec=rep(NA, length(tau.vec))
    dev.mat=matrix(NA, nrow=length(tau.vec), ncol=dim(C)[1])
    for(j in 1:length(tau.vec)){
      temp.tau=tau.vec[j]
      
      ## compute the test statistic for each subgroup
      temp.t.vec=rep(NA, num.of.subgroups)
      null.values=matrix(NA, nrow=num.of.subgroups, ncol=2)
      for(i in 1:num.of.subgroups){
        subgroup=test.matched.data[subgroup.num==unique.subg[i],]
        treated.y=subgroup[,1]
        control.y=subgroup[,2]
        temp.t.vec[i]=wilcox.test(treated.y-temp.tau, control.y, paired=T, alternative="two.sided")$statistic
      }
      
      temp.test= C %*% temp.t.vec
      if(dim(C)[1]==1){
        temp.deviate.upper=(temp.test-theta0.upper)/sqrt(sigma0)
        temp.deviate.lower=(temp.test-theta0.lower)/sqrt(sigma0)
      }else{
        temp.deviate.upper=(temp.test-theta0.upper)/sqrt(diag(sigma0))
        temp.deviate.lower=(temp.test-theta0.lower)/sqrt(diag(sigma0))
      }
      # if two deviate bounds have the same sign, choose the minimum; otherwise, give 0.
      same.sign=(temp.deviate.upper*temp.deviate.lower >= 0) # check whether two deviate bounds have the same sign.
      temp.deviate=rep(0, length(temp.deviate.upper))
      temp.deviate[same.sign==1]=pmin(abs(temp.deviate.upper), abs(temp.deviate.lower))[same.sign==1]
      dev.mat[j,]=temp.deviate
      max.deviate=max(abs(temp.deviate))
      max.dev.vec[j]=max.deviate
    }
    critical.val=qmvnorm(1-alpha/2, mean=rep(0, length(corr0[,1])), corr=corr0)$quantile
    
    sensi.mat[k,]=c(dev.mat[which.min(max.dev.vec),], min(max.dev.vec), critical.val)
    
  }
  sensi.mat=as.data.frame(sensi.mat)
  colnames(sensi.mat)[1:(dim(C)[1])] <- rownames(C)
  colnames(sensi.mat)[(dim(sensi.mat)[2]-1)]="Max"
  colnames(sensi.mat)[dim(sensi.mat)[2]]="kappa"
  sensi.mat$Gamma=sensi.param.vec
  
  return(sensi.mat)
}

#' De novo discovery of effect modification under no unmeasured confounder assumption.
#'
#' @param Data a total sample of matched pairs
#' @param split.ratio set a splitting ratio. 0.25 means using 25\% of data as a training sample.
#' @param total.significance a total significance level.
#' @param gamma a significance level that is used for estimating the confidence interval.
#' @return A result matrix with deviates is reported. The `Max` column represents the maximum value of deviates.
#' The `kappa` column represents the critical value. If `Max` is greater than 'kappa', we reject the null hypothesis.
#' @export
denovo=function(Data, split.ratio=0.25, total.significance=0.05, gamma=0.01){
  # creating a training sample and obtain a tree structure
  learning.from.training=denovo.training(Data=Data, split.ratio=split.ratio)
  
  # analysis under no unmeasured confounder assumption
  analysis=denovo.test(Data=Data, test.index=learning.from.training$test.index, tree=learning.from.training$tree, total.significance=total.significance, gamma=gamma)
  minmax.dev.vec=analysis[which.min(analysis$Max),]
  return(list(tree=learning.from.training$tree, deviate.mat=analysis, Dminmax=minmax.dev.vec))
}

#' De novo discovery of effect modification with a sensitivity analysis.
#'
#' @param Data a total sample of matched pairs.
#' @param Gamma.vec a vector of sensitivity parameter \eqn{\Gamma}.
#' @param split.ratio set a splitting ratio. 0.25 means using 25\% of data as a training sample.
#' @param total.significance a total significance level.
#' @param gamma a significance level that is used for estimating the confidence interval.
#' @return A result matrix with deviates is reported. The `Max` column represents the maximum value of deviates.
#' The `kappa` column represents the critical value. If `Max` is greater than 'kappa', we reject the null hypothesis.
#' @export
denovo.sensi=function(Data, Gamma.vec, split.ratio=0.25, total.significance=0.05, gamma=0.0001){
  # creating a training sample and obtain a tree structure
  learning.from.training=denovo.training(Data=Data, split.ratio=split.ratio)
  
  # Sensitivity analysis
  sensi.analysis=denovo.test.sensi(Data=Data, Gamma.vec=Gamma.vec, test.index=learning.from.training$test.index, tree=learning.from.training$tree, total.significance=total.significance, gamma=gamma)
  return(list(tree=learning.from.training$tree, deviate.mat=sensi.analysis))
}


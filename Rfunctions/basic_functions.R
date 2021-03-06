#' Combine a tree structure with a test sample
#'
#' @param tree a tree object
#' @param x a test sample with the same data structure that is used for creating `tree`.
#' @return a list with the elements
#' \item{C}{the matrix representation of a tree.}
#' \item{new.x}{a new test sample by adding a vector of subgroup memberships.}
#' @description The fuction combines a tree object (obtained from a training sample) with a test sample.
#' @import rpart.utils
tree.stru=function(tree, x){
  require(rpart.utils)
  # rpart.rules.table(tree)
  # rpart.subrules.table(tree)
  
  tree.frame=tree$frame
  nodes=as.numeric(rownames(tree.frame))
  terminal.node.ind=(tree.frame$var=="<leaf>")
  terminal.node.num=sum(terminal.node.ind)
  terminal.nodes=nodes[terminal.node.ind]
  
  if(length(nodes)==1){
    matrix.C=matrix(1, nrow=1, ncol=1)
    rownames(matrix.C)=1
    colnames(matrix.C)=1
    x$subgroup.num=rep(1, length(x[,1]))
    return(list(C=matrix.C, new.x=x))
  }
  
  parent.node.list=vector("list", terminal.node.num)
  for(i in 1:terminal.node.num){
    temp.t.node.num=(nodes[terminal.node.ind==1])[i]
    temp.parent.node.num=c()
    temp.parent.node.num[1]=temp.t.node.num
    j=1
    temp.p.node.end=0
    while(temp.p.node.end != 1){
      j=j+1
      if(temp.t.node.num %% 2 ==0){
        temp.p.node.end=temp.t.node.num/2
        temp.parent.node.num[j]=temp.p.node.end
      }else{
        temp.p.node.end=(temp.t.node.num-1)/2
        temp.parent.node.num[j]=temp.p.node.end
      }
      temp.t.node.num=temp.p.node.end
    }
    parent.node.list[[i]]=temp.parent.node.num
  }
  ## create C matrix
  matrix.C=matrix(0, nrow=length(nodes), ncol=length(terminal.nodes))
  rownames(matrix.C)=nodes
  colnames(matrix.C)=terminal.nodes
  
  for(i in 1:terminal.node.num){
    matrix.C[rownames(matrix.C) %in% parent.node.list[[i]], i]=1
  }
  matrix.C=matrix.C[-1,] # delete the first row
  
  
  ## assign the terminal node name
  decision.rule=rpart.rules.table(tree)
  var.rule=rpart.subrules.table(tree)
  
  x$subgroup.num=rep(NA, length(x[,1]))
  for(i in 1:terminal.node.num){
    temp.t.node=terminal.nodes[i]
    temp.d.rule=decision.rule[decision.rule$Rule==temp.t.node, "Subrule"]
    
    temp.v.rule=var.rule[var.rule$Subrule %in% temp.d.rule, ]
    temp.ind.matrix=matrix(NA, nrow=length(temp.d.rule), ncol=length(x[,1]))
    for(j in 1:length(temp.d.rule)){
      temp.var.name=temp.v.rule$Variable[j]
      x.column.num=which(colnames(x)==temp.var.name)
      if(is.na(temp.v.rule$Less[j])==1){
        temp.ind.matrix[j,]=(x[,x.column.num] >= as.numeric(as.character(temp.v.rule$Greater[j])))
      }else{
        temp.ind.matrix[j,]=(x[,x.column.num] < as.numeric(as.character(temp.v.rule$Less[j])))
      }
    }
    temp.ind=apply(temp.ind.matrix, 2, prod)
    x[temp.ind==1, "subgroup.num"]=rep(temp.t.node, sum(temp.ind))
  }
  return(list(C=matrix.C, new.x=x))
}


#' The Wilcoxon's signed rank sum test in a sensitivity analysis.
#'
#' @param N the number of matched pairs
#' @param gamma the value of sensitivity parameter \eqn{\Gamma}
#' @return
#' \item{expect.upper}{Expectation of the upper bound.}
#' \item{expect.lower}{Expectation of the lower bound.}
#' \item{var}{Variance of both upper and lower bounds.}
wilcoxSenMoments=function(N,gamma){
  #Computes null expectation and variance of Wilcoxon's statistic
  #in a sensitivity analysis with parameter Gamma
  #Uses formula (4.11) in Rosenbaum (2002) Observational Studies
  pplus.u=gamma/(1+gamma)
  pplus.l=1/(1+gamma)
  expect.u=pplus.u*N*(N+1)/2
  expect.l=pplus.l*N*(N+1)/2
  vari=pplus.u*(1-pplus.u)*N*(N+1)*(2*N+1)/6
  list(expect.upper=expect.u,expect.lower=expect.l,var=vari)
}

#' Estimation of the confidence interval for the population mean by inverting the Wilcoxon signed rank sum test.
#'
#' @param data input a test sample.
#' @param gamma set a significance level for the confidence interval estimation.
#' @param grid.size the number of values within the confidence interval.
#' @import stats
#' @return The \eqn{100(1-\gamma)}% confidence interval for the population mean is estimated.
get.tau.vector.wilcox=function(data, gamma=0.001, grid.size=21){
  wilcox.res=wilcox.test(data[,1], data[,2], paired=T, alternative="two.sided", conf.int=T, conf.level=1-gamma)
  conf.int.vec=seq(wilcox.res$conf.int[1], wilcox.res$conf.int[2], by=(wilcox.res$conf.int[2]-wilcox.res$conf.int[1])/(grid.size-1))
  return(conf.int.vec)
}

get.tau.vector = function( Y, Z, X = NULL, gamma=0.001, grid.size=21, grid.gamma = 100*gamma ) {
  
  if ( is.null( X ) ) {
    
    te.hat <- mean(Y[Z == 1]) - mean(Y[Z == 0])
    te.se <- sqrt( var(Y[Z == 1])/sum(Z) + var(Y[Z == 0])/sum(1 - Z))
    
  } else {
    lm.tau <- lm( Y ~ Z + X )
    
    te.hat <- as.numeric(coef(lm.tau)["Z"])
    te.se <- as.numeric(sqrt( diag(vcov(lm.tau))["Z"]))
  }
  
  # oversample points near the estimated tau-hat
  te.MOE <- qnorm(1 - gamma/2)*te.se
  te.vec <- te.hat + (te.MOE/qnorm(1-grid.gamma)) * qnorm( seq( grid.gamma, 1-grid.gamma, length.out=grid.size ) )
  te.vec
  
  attr( te.vec, "te.hat" ) <- te.hat
  attr( te.vec, "te.se" ) <- te.se
  attr( te.vec, "te.MOE" ) <- te.MOE
  
  te.vec
}
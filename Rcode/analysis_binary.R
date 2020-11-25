##############################################################################
### Application of the denovo method with a simulated dataset 
##############################################################################
#### 
source("../Rfunctions/basic_functions.R")
source("../Rfunctions/denovo.R")
source("../Rfunctions/functions_binary.R")
source("../Rfunctions/functions_binary_sensi.R")

## Import the simulated dataset (binary outcomes)
matched.pairs=read.csv("simulated_dataset.csv")

## Use individual-level covariates only since exact matching is used for these covariates. 
matched.pairs.indiv = matched.pairs[,1:6]
##################################
## Binary outcomes require more intensive computation.
## So, application of the denovo method is slightly more complicated.
## We will look at the discovery step and the inference step separately.  
##################################
## Discovery Step
##################################
## Using CART
## (a) If you want to randomly select the discovery (training) subsample, 
##      then use denovo.training function with the full data
## (b) If you want to use a pre-selected index indicating the discovery subsample, 
##      then use the denovo.training function with the selecting subsample with split.ratio = 1
## Find a regression tree to illustrate the structure of heterogenous effects

## (a) use a randomly selected index vector 
denovo.discovery = denovo.training(Data = matched.pairs.indiv)
denovo.discovery$tree

discovery.index = denovo.discovery$training.index
inference.index = denovo.discovery$test.index

## (b) use a pre-selected index vector
discovery.index = read.csv("discovery_set_index.csv")$x
denovo.discovery = denovo.training(Data = matched.pairs.indiv[discovery.index,], split.ratio = 1)
inference.index = seq(1:dim(matched.pairs.indiv)[1])[-discovery.index]
rpart.plot(denovo.discovery$tree)
#################
## Causal Tree
library(causalTree) # to install this package, see the README file. 

# to use the Causal Tree method, need to transform the dataset
matched.set = matched.pairs.indiv
treated.matched.pair = matched.set[,c(1,3,4,5,6)]
control.matched.pair = matched.set[,c(2,3,4,5,6)]
colnames(treated.matched.pair)[1] = "death5" # 5-year mortality
colnames(control.matched.pair)[1] = "death5"
ct.matched.pair = rbind(treated.matched.pair, control.matched.pair)
ct.matched.pair$treatment = c(rep(1,dim(treated.matched.pair)[1]), rep(0,dim(control.matched.pair)[1]))

N = dim(matched.pairs.indiv)[1]
ct.matched.pair.discovery = ct.matched.pair[c(discovery.index, N + discovery.index),]
ct.matched.pair.inference = ct.matched.pair[c(inference.index, N + inference.index),]

tree=causalTree(death5 ~ age + male + white + medicaid,
                data=ct.matched.pair.discovery, treatment = ct.matched.pair.discovery$treatment,
                HonestSampleSize=2*length(inference.index),
                split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, xval = 5,
                cp = 0, minsize = 100, propensity = 0.5)

opcp=tree$cptable[,1][which.min(tree$cptable[,4])]
opfit=prune(tree, opcp) ## obtained tree from a training sample
rpart.plot(opfit)


###################
## Need to create the C matrix as in the manuscript. 
## At the same time, we add a new variable of the subgroup membership
## Adding this new variable to the inference subsample!

est.ct.tree = tree.stru(tree = opfit, x = ct.matched.pair.inference)
est.cart.tree = tree.stru(tree = denovo.discovery$tree, x = matched.pairs.indiv[inference.index,])

##################################
## Inference Step
##################################

### Test sample -> hypothesis test
subgroup.num=est.cart.tree$new.x$subgroup.num
C=est.cart.tree$C
unique.subg=as.numeric(colnames(est.cart.tree$C))
num.of.subgroups=length(unique.subg)

test.matched.pairs.ymat = matched.pairs.indiv[inference.index,]
test.matched.pairs.ymat$sub.num=subgroup.num

##########################################
n.test=length(test.matched.pairs.ymat[,1])
n.vec= rep(NA, num.of.subgroups)
for(i in 1:num.of.subgroups){
  n.vec[i] = sum(test.matched.pairs.ymat$sub.num==unique.subg[i])
}

# Since the population treatment effect size (mu) is not known, 
## get the 99% confidence interval of delta
est.delta=sum(test.matched.pairs.ymat[,1]-test.matched.pairs.ymat[,2])/n.test
eta=0.01
var.max=binary.hypo.test(Ymat=test.matched.pairs.ymat, delta0=est.delta)$max.var
rough.int=c(est.delta-qnorm(1-eta/2)*sqrt(var.max), est.delta+qnorm(1-eta/2)*sqrt(var.max))

delta.end.points=round(rough.int*n.test*2)
delta.seq=round(seq(delta.end.points[1], delta.end.points[2], length.out=20))

### subgroup estimates
est.delta.vec=rep(NA, length(C[1,]))
for(j in 1:length(C[1,])){
  est.delta.vec[j]=mean(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],1]-test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num==unique.subg[j],2])
}

### compute deviates when Gamma = 1 (no unmeasured confounding bias)

library(mvtnorm)
max.dev.vec=max.dev.pval.vec=rep(NA, length(delta.seq))
deviate.mat=matrix(NA, nrow=length(delta.seq), ncol=2*length(C[1,])-2)
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
  
  deviate.mat[i,]=deviate
  max.dev.vec[i]=max.deviate
  if(dim(C)[1]==1){
    max.dev.pval.vec[i]=2*(1-pnorm(max.deviate))
  }else{
    max.dev.pval.vec[i]=2*(1-pmvnorm(lower=-Inf, upper=max.deviate, mean=rep(0, length(corr0[,1])), corr=corr0)[1])
  }
}
pval.max=max(max.dev.pval.vec)+eta
if(pval.max > 1){
  pval.max=1
}



dev.matrix=cbind(round(delta.seq/(2*n.test), 4), round(deviate.mat, 2), round(max.dev.vec, 2))
colnames(dev.matrix) = c("Delta", rownames(C), "Dmax")

dev.matrix
2*(1-pmvnorm(lower=-Inf, upper=max(max.dev.vec), mean=rep(0, 2*length(C[1,])-2), corr=corr0)[1])+eta
qmvnorm(0.975+eta/2, corr=corr0)[1]

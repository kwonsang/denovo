# "gurobi" package needs to be installed. 
# https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation.html
# install.packages("c:/gurobi910/win64/R/gurobi_9.1-0.zip", repos = NULL)

########################################################################
## Sensitivity Analysis (need to run "analysis_binary.R" first)
#######################################################################
library(Matrix)
library(gurobi)
sub.delta.seq=seq(0, 6800, by=400)

gamma.vec=c(1, 1.05, 1.10)

total.array = array(NA, dim = c(length(sub.delta.seq), 2*length(C[1,]), length(gamma.vec)))
total.critical.mat=NULL
for(k in 1:length(gamma.vec)){
  gamma=gamma.vec[k]
  library(mvtnorm)
  max.dev.vec=rep(NA, length(sub.delta.seq))
  deviate.mat=matrix(NA, nrow=length(sub.delta.seq), ncol=2*length(C[1,])-2)
  critical.constant.vec=rep(NA, length(sub.delta.seq))
  for(i in 1:length(sub.delta.seq)){
    temp.delta0=sub.delta.seq[i]
    
    temp.maxvar.vec=rep(NA, length(C[1,]))
    temp.delta0.vec=rep(NA, length(C[1,]))
    for(j in 1:length(C[1,])){
      fraction.val=(temp.delta0/(2*n.test))*2*n.vec[j]
      closest.two.points=c(floor(fraction.val), ceiling(fraction.val))
      
      index.vec=rep(NA, 2*n.vec[j])
      index.vec[seq(1, 2*n.vec[j]-1, by=2)]=1:n.vec[j]
      index.vec[seq(2, 2*n.vec[j], by=2)]=1:n.vec[j]
      
      R.vec=as.vector(t(test.matched.pairs.ymat[test.matched.pairs.ymat$sub.num == unique.subg[j],1:2]))
      
      Z.vec=rep(0, 2*n.vec[j])
      Z.vec[seq(1, 2*n.vec[j]-1, by=2)]=1
      
      temp1=sensCRD(index=index.vec, outcome=R.vec, treatment=Z.vec, null=closest.two.points[1], Gamma.vec=gamma)
      temp2=sensCRD(index=index.vec, outcome=R.vec, treatment=Z.vec, null=closest.two.points[2], Gamma.vec=gamma)
      
      max.temp.maxvar=max(temp1$maxvar, temp2$maxvar)
      if(which.max(c(temp1$maxvar, temp2$maxvar))==1){
        temp.mean=temp1$mean
      }else{
        temp.mean=temp2$mean
      }
      temp.maxvar.vec[j]=max.temp.maxvar
      temp.delta0.vec[j]=temp.mean
      #print(c(i,j))
    }
    T.vec=est.delta.vec*n.vec*2
    var.mat=diag(temp.maxvar.vec)
    test= C %*% T.vec
    theta0= C %*% temp.delta0.vec
    om0=C %*% var.mat %*% t(C)
    deviate=(test-theta0)/sqrt(diag(om0))
    corr0=om0/sqrt(outer(diag(om0),diag(om0),"*"))
    max.deviate=max(abs(deviate))
    
    #max.pval[i]=2*(1-pmvnorm(lower=-Inf, upper=max.deviate, mean=rep(0, 11), corr=corr0)[1])
    deviate.mat[i,]=deviate
    max.dev.vec[i]=max.deviate
    critical.constant.vec[i]=qmvnorm(0.975, corr=corr0)[1]
    #if(i%%100==0)cat("..", i)
  }

  total.array[,,k] = cbind(round(sub.delta.seq/(2*n.test),4), round(deviate.mat, 2), round(max.dev.vec, 2))
  total.critical.mat=rbind(total.critical.mat, unlist(critical.constant.vec))
  print(k)
}


### Make a plot 
plot(total.array[,1,1], total.array[,2*length(C[1,]),1], type="n", xlim=c(0, 0.047), ylim=c(0, 14), xlab=expression(delta[0]), ylab=expression(D[paste(Gamma, "max")]))
for(i in c(1:length(gamma.vec))){
  lines(total.array[,1,i], total.array[,2*length(C[1,]),i], lwd=1.5)
}
abline(h=total.critical.mat[1,1], lty=2)


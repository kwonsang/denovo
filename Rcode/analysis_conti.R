##############################################################################
### Application of the denovo method with a simulated dataset 
##############################################################################
#### 
source("../Rfunctions/basic_functions.R")
source("../Rfunctions/denovo.R")
source("../Rfunctions/functions_binary.R")
source("../Rfunctions/functions_binary_sensi.R")

## Import the simulated dataset (continuous outcomes)
matched.pairs=read.csv("simulated_dataset_continuous.csv")

## Use individual-level covariates only since exact matching is used for these covariates. 
matched.pairs.indiv = matched.pairs[,1:6]
#################
library(rpart.plot)
## simple demonstration of the denovo method (without sensitivity analysis)
denovo.without.sensi = denovo(Data = matched.pairs.indiv)
denovo.without.sensi 
rpart.plot(denovo.without.sensi$tree)

## simple demonstration of the denovo method (with sensitivity analysis)
denovo.with.sensi = denovo.sensi(Data = matched.pairs.indiv, Gamma.vec = c(1, 1.02, 1.04, 1.06, 1.08, 1.1))
denovo.with.sensi
rpart.plot(denovo.with.sensi$tree)

#################
## denovo.without.sensi and denovo.with.sensi may provide different results. 
## This difference is because we didn't specify the index of the test dataset. 
## We can see the discovery & inference steps separately. 


## training (discovery step)
denovo.discovery = denovo.training(Data = matched.pairs.indiv)
denovo.discovery$tree

discovery.index = denovo.discovery$training.index
inference.index = denovo.discovery$test.index

## test (inference step)
denovo.inference = denovo.test(Data = matched.pairs.indiv,
                               test.index = inference.index,
                               tree = denovo.discovery$tree)
denovo.inference.sensi = denovo.test.sensi(Data = matched.pairs.indiv,
                                           test.index = inference.index,
                                           tree = denovo.discovery$tree,
                                           Gamma.vec = c(1,1.02,1.04,1.06,1.08,1.1))

denovo.inference
denovo.inference.sensi


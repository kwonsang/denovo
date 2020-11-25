

This zip file aims to achieve two main goals: (1) application of the denovo method with simulated datasets and (2) replication of simulation/application results.
----------------------------------------------------------------
Overview 

(1) Application
This part is aimed to replicating the results from the Medicare dataset discussed in Section 5. We are unable to share the original data, however, we include a simulated dataset, called "simulated_dataset.csv." We provide a script that uses the simulated dataset to conduct all the analysis used in the manuscript. Specifically, Figures 2,3 and Table 4 can be replicated in a similar fashion. 

(2) Simulation 
This part is aimed to replicating the simulation results such as Table 1 in the manuscript and Table 5 in the online supplementary materials.  


----------------------------------------------------------------
Description of Directories

There are three directories (i) /Rcode, (ii) /simulation, (iii) /Rfunctions
(i) Rcode - contains two simulated datasets and three R scripts. 
	(a) "analysis_conti.R" - application of the denovo method with continuous outcomes (with sensitivity analysis)
	(b) "analysis_binary.R" - application of the denovo method with binary outcomes (without sensitivity analysis)
	(c) "analysis_binary_sensi.R"- sensitivity analysis when outcomes are binary
	(d) "simulated_dataset.csv" - this dataset contains 110,091 matched pairs with 4 individual level covariates (age, sex, race, medicaid eligibility) and 11 zip code-level covariates. The Zip code-level covariates are the within-pair averages. The outcomes are simulated by the authors. 
	(e) "simulated_dataset_continuous.csv" - this dataset is similar to the previous dataset. However, the outcomes are generated as continuous outcomes instead. This dataset will be used with "analysis_conti.R."
	(f) "discovery_set_index.csv" - this includes an index vector indicating that which subjects were used as the discovery subsample. The same index vector was used in analyzing the actual Medicare data used in the manuscript.  

"analysis_conti.R" uses the simulated dataset "simulated_dataset_continuous.csv." Since handling continuous outcomes is much simpler than handling binary outcomes, we strongly encourage users start with this script first. In this script, we demonstrates how to implement the denovo method in order to discover a tree structure and conduct hypothesis tests. 

"analysis_binary.R" and "analysis_binary_sensi.R" can then be considered to deal with binary outcomes. These two R script demonstrate the denovo method discussed in Section 3.5 of the manuscript. 

--
(ii) simulation - contains two R scripts. 
	(a) "pvalue_discovery_continuous.R" 
	(b) "pvalue_discovery_binary.R"

The data generating process is discussed in Section 4. We considered 5 covariates, and two are effect modifiers among five. Also, we considered three different splitting ratios (10%, 90%), (25%, 75$), (50%, 50%). 

The outputs for both R scripts consist of two matrices: (1) pval.matrix and (2) check.matrix. pval.matrix includes the p-values that was used for power computation like in Table 1. check.matrix includes the discovery ratio of each covariate. This was used for Table 5 in the supplementary materials. 

First, the number of columns for pval.matrix is 12. Each splitting ratio takes 4 columns. Among these four columns, first two are from CART and the next two are from CT. The first two columns represent the p-values by the truncated product method and the denovo method. Therefore, for each splitting ratio, (1) p-value from truncated product with CART, (2) p-value from denovo with CART, (3) p-value from truncatedP with CT, and (4) p-value from denovo with CT. 

Second, the number of columns for check.matrix is 30. Each splitting ratio takes 10 columns. Among these 10 columns, first five are from CART and the next five are from CT. The first five columns represent whether each covariates x_i is used for split or not in CART. Similarly, the next five columns are defined for CT.  

--
(iii) Rfunctions - contains four R functions
	(a) "basic_fucntions.R" - containing all necessary functions to run other R scripts discussed above. 
	(b) "denovo.R" - containing functions for the denovo method with continuous outcomes 
	(c) "functions_binary.R" - containing functions for the denovo method with binary outcomes 
	(d) "functions_binary_sensi.R - containing functions for sensitivity analysis with binary outcomes 
 
Both "functions_binary.R" and "functions_binary_sensi.R" are modified from the R scripts provided by the supplementary materials of Fogarty et al. (2016) and Fogarty et al. (2017)


----------------------------------------------------------------
Computation Time 

(i) Rcode 
	(a) "analysis_conti.R" 
		(a1) the denovo function with 110,091 matched pairs takes 6.7 sec.
		(a2) the denovo.sensi function with 110,091 matched pairs & 6 values of Gamma takes 20.4 sec. 
	(b) "analysis_binary.R" takes 1.03 sec. 
	(c) "analysis_binary_sensi.R" takes 2.84 minutes (110,091 matched pairs & 3 values of Gamma)

(2) simulation 
	(a) "pvalue_discovery_continuous.R" takes 3.4 sec. for one simulation/repetition. The authors used 10,000 repetitions. 
	(b) "pvalue_discovery_binary.R" takes 4.6 sec. for one simulation. 

Computer specification
	Model Name: iMac Pro
	Processor Name: 8-Core Intel Xeon W
	Processor Speed: 3.2 GHz
	Number of Processors:	1
	Total Number of Cores: 8
	Memory: 64 GB
----------------------------------------------------------------
Computing Information:
R version 3.6.3 (2020-02-29)Platform: x86_64-apple-darwin15.6.0 (64-bit)Running under: macOS Catalina 10.15.7


To run the R scripts, several packages should be installed in advance. 

Attached packages:
sensitivitymv_1.4.3 
gurobi_9.1-0
slam_0.1-47
causalTree_0.0
data.table_1.13.2
Matrix_1.2-18
lpSolve_5.6.15
rpart.plot_3.0.9 
rpart.utils_0.5
mvtnorm_1.1-1
rpart_4.1-15     
----------------------------------------------------------------
There are two R packages (causalTree & gurobi) that cannot be installed from CRAN. Users need to install these packages manually.

How to install the "causalTree" package

	Use the following instruction: https://github.com/susanathey/causalTree

	install.packages("devtools")
	library(devtools) 
	install_github("susanathey/causalTree")

How to install the "gurobi" package

	Use the instruction provided by CRAN: https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation.html


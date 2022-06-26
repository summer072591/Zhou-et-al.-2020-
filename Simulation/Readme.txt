Each folder contains all R code needed for the simulation setting, where the true coefficient
vectors are sampled from (-1, 1, 0) ("discrete_dist_beta") or 
N(0, 1) + \delta_0 ("norm_dist_beta") with probability (0.005, 0.005, 0.99) ("XXX_sparse5")or 
(0.05, 0.05, 0.9) ("XXX_sparse50"). 

Folders were created for parallel computing on VSC. The R files named "WrapFull_wrapper_XXX.R"
and the batch files "WrapFull_XXX.pbs" are used for running jobs on VSC. 


To use implementations locally, please load the files "simu_test_compare_XXX.R", then use the function
"Wrap_Full(s)" where s refers to the seed number.
The bandwidth value depend on the loss function, signal-to-noise ratio; should be adjusted accordingly 
in order to achieve convergence

Same files are used in the real data section, but "simu_test_compare_XXX.R" should be replaced by 
"simu_test_compare.R" under the folder "real data section".

NOTICE!!! Original data used for producing the tables in Zhou et al. are stored in the folders named
"data" under the directories for each simulation setting.
Folders were created for parallel computing on VSC. 
Same files are used in simulation section, but "simu_test_compare_XXX.R" in the simulation 
section should be replaced by "simu_test_compare.R" under the folder "real data section".


To use implementations locally, please load the files "simu_test_compare.R", then use the function
"Wrap_Full(s)" where s refers to the seed number.
The bandwidth value depend on the loss function, signal-to-noise ratio; should be adjusted accordingly 
in order to achieve convergence.

For obtaining model-averaged estimators using Bates-Granger type weights, please replace the "main_body_of_the_simulation.R" file with the "main_body_of_the_simulation_BatesGranger.R" file.


NOTICE!!! Original data used for producing figures and Tables in the real data section 
          in Zhou et al.(2019) are stored in the folders named "Data".
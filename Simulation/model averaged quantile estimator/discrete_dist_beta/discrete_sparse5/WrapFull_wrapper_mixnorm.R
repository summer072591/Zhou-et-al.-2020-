args = commandArgs(TRUE)
s    = as.double(args[1])




library("mvtnorm")
library("lpSolve")
library("glmnet")
library("quadprog")
library("methods")

source("simu_test_compare_mixnorm.R")
WrapFull(s)
# R3.6
rm(list=ls())
condaENV <- "/home/chenzh/miniconda3/envs/R3.6"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

condaENV <- "/home/chenzh/miniconda3/envs/R3.6"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(pheatmap)
  library(ggplot2)
  library(foreach)
  library(doParallel)
  library(SCopeLoomR)
  
  
})
suppressMessages(library(loomR))

options(digits = 4)
options(future.globals.maxSize= 3001289600)
numCores <- 10
registerDoParallel(numCores)

#!/usr/bin/env Rscript

library(data.table)
library(logger)
library(ggplot2)
library(dplyr)
library(doMC)
registerDoMC(4)

setwd("/media/nathan/T7/path_integral/simulations/out/trees")
list.files()

start_freqs <- data.table()

count <- 0
for(cur_sim in list.files()){
  cur_seed <- as.integer(strsplit(cur_sim, split = "_")[[1]][2])
  count <- count + 1
  print(paste(count, ":", cur_sim))
  raw <- as.matrix(fread(paste(cur_sim,"/start_freqs.csv", sep = ""), skip = 2))
  df <- foreach(i = 1:nrow(raw), .combine = rbind) %dopar% {
    site_start <- raw[i,]
    data.frame(seed = cur_seed,
               deme = 1:100,
               start = site_start)
  }
  start_freqs <- dplyr::bind_rows(start_freqs, df)
}



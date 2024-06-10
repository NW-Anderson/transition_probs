library(data.table)
library(dplyr)
library(ggplot2)

setwd("/media/nathan/T7/path_integral/unlinkedSims")
list.files()
master <- data.frame()
for(file in list.files()){
  params <- strsplit(file, split = "r.")[[1]][2] 
  params <- strsplit(params, split = ".t")[[1]][1]
  print(params)
  print("\n")
  df <- fread(file)
  df <- df %>% filter(x0 < 0.11,
                      x0 > 0.09,
                      a > 0)
  start <- df$x0
  
  end_freqs <- c()
  deme <- c()
  for(i in 5:104){
    end_freqs <- c(end_freqs, df[,..i][[1]])
    deme <- c(deme, names(df)[i])
  }
  
  tmp <- data.frame(deme, 
                    start, 
                    end_freqs,
                    group_id  = params )
  master <- dplyr::bind_rows(master, tmp)
}

rm(tmp)

ggplot(master, aes(x = end_freqs)) + geom_density() + facet_wrap(vars(group_id))

fwrite(master, "unlinked_positive_eff_09_11.csv.gz")

#############################################

setwd("/media/nathan/T7/path_integral/unlinkedSims")
list.files()
master <- data.frame()
for(file in list.files()){
  params <- strsplit(file, split = "r.")[[1]][2] 
  params <- strsplit(params, split = ".t")[[1]][1]
  print(params)
  print("\n")
  df <- fread(file)
  repVG <- df %>% group_by(rep) %>% summarize(meanVG = unique(meanVG))
  tmp <- data.frame(VG = mean(repVG$meanVG),
                    group_id  = params )
  master <- dplyr::bind_rows(master, tmp)
}
  
library(data.table)
library(tidyr)
library(dplyr)

setwd("/home/nathan/Documents/GitHub/path_integral/simulations/")
params <- fread("params.txt")
colnames(params) <- c("seed", "U", "a")

setwd("/media/nathan/T7/path_integral/simulations/out/trees(neut)")

master <- data.frame()

count <- 0
for(file in list.files()){
  count <- count + 1
  print(paste(count, " : ", file))
  
  cur_seed <- as.integer(strsplit(file, split = "_")[[1]][2])
  cur_par <- params %>% filter(seed == cur_seed)
  cur_par <- paste("U=", cur_par$U, "_a=", cur_par$a, sep = "")
  
  end_freqs <- fread(paste(file, "/end_freqs.csv", sep = ""), skip = 2) %>% 
    mutate(site = 1:nrow(.)) %>% 
    melt(id.vars = "site")
  colnames(end_freqs) <- c("site", "deme", "end_freq")
  
  start_freqs <- fread(paste(file, "/start_freqs.csv", sep = ""), skip = 2) %>%
    mutate(site = 1:nrow(.)) %>% 
    melt(id.vars = "site")
  colnames(start_freqs) <- c("site", "deme", "start_freq")
  freqs <- merge(start_freqs, end_freqs) %>% 
    filter(start_freq < 0.11,
           start_freq > 0.09) %>%
    mutate(seed = cur_seed,
           par = cur_par)
  
  master <- dplyr::bind_rows(master, freqs)
}

fwrite(master, file = "../linked_neutral_09_11.csv.gz")

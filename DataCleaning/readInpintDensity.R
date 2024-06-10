library(data.table)
library(dplyr)
library(ggplot2)

setwd("/media/nathan/T7/path_integral/simulations/out/pintComparison")

df <- data.frame()
for(f in list.files()){
  if(!grepl("StatDist",f)){
    tmp <- strsplit(f, split="_")
    scen <- tmp[[1]][2]
    end_freq <- tmp[[1]][3]
    scen <- as.integer(gsub("Scenario","", scen))
    end_freq <- end_freq %>% gsub("end","", .) %>%
      gsub(".csv","", .) %>% as.numeric()
    dens <- fread(f)[[1,1]]
    df <- dplyr::bind_rows(df, data.frame(Scenario = scen,
                                          end_freq = end_freq,
                                          Dens = dens))
  }
}
# write.csv(df, file = "linked_sim_pint_densities.csv")

setwd("/media/nathan/T7/path_integral/simulations/out/KimuraComparison")

df <- data.frame()

for(f in list.files()){
  if(!grepl("StatDist",f)){
    tmp <- strsplit(f, split="_")
    scen <- tmp[[1]][2]
    end_freq <- tmp[[1]][3]
    end_freq <- end_freq %>% gsub("end","", .) %>%
      gsub(".csv","", .) %>% as.numeric()
    dens <- fread(f)[[1,1]]
    df <- dplyr::bind_rows(df, data.frame(Scenario = scen,
                                          end_freq = end_freq,
                                          Dens = dens))
  }
}

write.csv(df, file = "linked_sim_kimura_densities.csv")

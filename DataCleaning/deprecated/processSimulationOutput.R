library(ggraptR)
library(poolSeq)
library(patchwork)
library(dplyr)
library(ggridges)
library(tidyr)
library(viridis)
library(colorspace)
library(gganimate)
library(gifski)
library(av)
library(png)
library(magick)
library(scico)
library(plotrix)
library(ggdensity)
library(ggplot2) # needs to be version ≥ 2.1.0
library(scales)
library(devtools)
library(network)
library(sna)
library(GGally)
library(geomnet)
library(ggnetwork)
library(igraph)
library(ggraph)
library(tidygraph)
library(vcfR)
library(ggplotify)
library(pheatmap)
library(ggbreak)
library(ggimage)
library(rsvg)
library(ggrepel)
library(shadowtext)
dev.off()


rawData <- fread("test.csv", header = T)
master <- data.table()
for(site in 1:nrow(rawData)){
  obs <- rawData[site]
  for(rep in 1:ncol(obs)){
    tpl <- as.character(obs[,..rep])
    tpl <- tpl %>% gsub('\\[', '', .) %>%
      gsub('\\]', '', .) %>%
      strsplit(.,split = " ")
    start = as.numeric(tpl[[1]][1])
    end = as.numeric(tpl[[1]][2])
    master <- dplyr::bind_rows(master, 
                               data.frame(population = rep - 1,
                                          site = site,
                                          start = start,
                                          end = end))
  }
}
rawData <- fread("effect_sizes.csv")
effect_sizes <- data.table()
for(i in 1:ncol(rawData)){
  effect_sizes <- dplyr::bind_rows(effect_sizes,
                                   data.frame(site = i,
                                              a = rawData[[i]]))
}
master <- dplyr::left_join(master, effect_sizes, by = c("site"))

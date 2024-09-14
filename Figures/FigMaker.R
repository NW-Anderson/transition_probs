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
library(Cairo)
library(extrafont)
dev.off()


####################
#### Comparison ####
####################
rm(list=ls())

setwd("/media/nathan/T7/path_integral/mainComparison/x-0.2")

master <- data.frame()
for(file in list.files()){
  end <- file %>% strsplit(., split = "_")
  end <- end[[1]][2] %>% 
    gsub(".csv", "", .) %>% 
    as.numeric()
  
  tmp <- fread(file) %>% unlist()
  master <- dplyr::bind_rows(master, 
                             data.frame(dens = tmp,
                                        vg = factor(c(1e-4, 1e-3, 1e-2, 
                                                      1e-1, "Genic", "Neutral"),
                                                    levels = c("Genic", 1e-4, 1e-3, 1e-2, 
                                                               1e-1, "Neutral")),
                                        end = end))
}

textScale <- 1.5

ggplot(master, aes(x = end, y = dens)) + 
  geom_line(aes(color = vg),
            alpha = 0.75,
            size = 1.5) +
  theme_bw() +
  guides(color=guide_legend(title = bquote(italic(V[G])),
                            override.aes = list(alpha=1),
                            ncol = 2)) +
  scale_x_continuous("Ending frequency",
                     expand = c(0,0)) + 
  scale_y_continuous("Density",
                     expand = c(0,0),
                     limits = c(0,2.35)) + 
  scale_color_manual(values = turbo(10)[c(2,3,4,7,8,9)]) + 
  theme(text = element_text(family = "LM Roman 10"),
        axis.title = element_text(size= textScale * 12),
        title = element_text(size = textScale * 12),
        axis.text = element_text(size = textScale * 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = textScale * 10),
        legend.title = element_text(size = textScale * 12),
        legend.justification = c("right", "top"),
        legend.position = c(.98,.98),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid = element_blank()) 

         
####################################
####### pDetection Alpha VG #######
####################################

rm(list=ls())

setwd("~/Documents/GitHub/path_integral/results/improvedPDetectionAlphaVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "start", 
                   "time", "thresh", "pintAUC", 
                   "pintDetected", "numAUC", "numDetected",
                   "statDist")
master$popalpha = 1000 * master$alpha
rm(tmp, file)

master$statDist <- statDist <- master$statDist %>% 
  gsub('[{}]', '', .) %>% 
  gsub("\\*\\^", "e", .) %>% 
  as.numeric()

pintDf <- master %>% filter(statDist < 0.05)

numDf <- master %>% filter(VG == 1e-4,
                           popalpha >= 15)
numDf$yval <- c(numDf[1,]$pintDetected,
                numDf[2,]$numDetected)

p1 <- ggplot(data = pintDf, aes(x = popalpha, y = pintDetected)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  geom_hline(yintercept=0.01, 
             linetype="dashed", 
             size  = 0.75,
             alpha = 0.6) + 
  theme_bw() +
  scale_x_continuous(bquote(italic(2 * N[e] * alpha * Lambda * "/" *  W)), 
                     breaks = c(0, 1, 5, 10, 15, 20), 
                     labels = c(0, 1, 5, 10, 15, 20),
                     expand = c(0.05,0)) +
  scale_y_continuous("Probability detected\n(Q)",
                     # breaks = seq(0,0.3,by = 0.05),
                     expand = c(0,0), 
                     limits = c(0,max(numDf$yval)*1.05)) +
  geom_point(data = numDf, aes(x = popalpha, 
                               y = yval,
                               color = as.factor(VG)),
             alpha=c(0,1),
             size = 2.5,
             shape = 17,
             show.legend = F) +
  geom_line(data = numDf, aes(x = popalpha, 
                              y = yval,
                              color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5,
            linetype = "dashed")  + 
  guides(color=guide_legend(title = bquote(italic(V[G])),
                            override.aes = list(alpha=1,
                                                size = 1.25,
                                                linewidth = 0.75))) +
  scale_color_manual(values = turbo(10)[c(2,4,7,9)]) + 
  theme(text = element_text(family = "LM Roman 10"),
        axis.title = element_text(size=12),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification = c("left", "top"),
        legend.position = c(.01,.99),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) 

p1

########################
#### Error Alpha VG ####
########################

# rm(list=ls())

setwd("~/Documents/GitHub/path_integral/results/improvedPDetectionAlphaVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "start", 
                   "time", "thresh", "pintAUC", 
                   "pintDetected", "numAUC", "numDetected",
                   "statDist")
master$popalpha = 1000 * master$alpha
rm(tmp, file)

master$statDist <- statDist <- master$statDist %>% 
  gsub('[{}]', '', .) %>% 
  gsub("\\*\\^", "e", .) %>% 
  as.numeric()

summary(master$statDist)
# plot(density(master$statDist))

# master <- master %>% filter(VG != 1e-04)

e1 <- ggplot(data = master, aes(x = popalpha, y = statDist)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  theme_bw() +
  guides(color=guide_legend(title = bquote(V[G]),
                            override.aes = list(alpha=1,
                                                size = 1.25,
                                                linewidth = 0.75))) +
  scale_x_continuous(bquote(2 * N[e] * "\u03b1 \u039b / W"), 
                     breaks = c(0, 1, 5, 10, 15, 20), 
                     labels = c(0, 1, 5, 10, 15, 20)) +
  ylab("Statistical Distance") +
  scale_y_log10() + 
  scale_color_manual(values = turbo(10)[c(2,4,7,9)]) + 
  theme(axis.title = element_text(size=12),
      title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
                                 colour = c(rep("black",1),
                                            "red",
                                            rep("black",4))),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10),
      legend.justification = c("left", "top"),
      legend.position = c(.01,.99),
      legend.box.background = element_rect(colour = "black"),
      legend.spacing.y = unit(0, 'cm'),
      legend.key.size = unit(0.8, "line"),
      panel.grid.minor.x = element_blank())  + 
  annotation_logticks(sides = "l") 

e1

#########################
#### pDetection Time ####
#########################

# rm(list=ls())

setwd("~/Documents/GitHub/path_integral/results/improvedPDetectedTimeVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "start", 
                   "time", "thresh", "pintAUC", 
                   "pintDetected", "numAUC", "numDetected",
                   "statDist")
master$popalpha = 1000 * master$alpha
rm(tmp, file)

master$statDist <- statDist <- master$statDist %>% 
  gsub('[{}]', '', .) %>% 
  gsub("\\*\\^", "e", .) %>% 
  as.numeric()

pintDf <- master %>% filter(statDist < 0.05)

numDf <- master %>% filter(VG == 1e-4,
                           time >= 0.15)
numDf$yval <- c(numDf[1,]$pintDetected,
                numDf[2,]$numDetected)

# pintDf <- pintDf %>% filter(VG == 0.1)

p2 <- ggplot(data = pintDf, aes(x = time, y = pintDetected)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  geom_hline(yintercept=0.01, 
             linetype="dashed", 
             size  = 0.75,
             alpha = 0.6) + 
  theme_bw() +
  scale_x_continuous("Time\n(genomic units)",
                     expand = c(0.05,0)) +
  scale_y_continuous("P(detected)",
                     expand = c(0,0),
                     limits = c(0,max(numDf$yval) * 1.05)) +
  geom_point(data = numDf, aes(x = time,
                               y = yval,
                               color = as.factor(VG)),
             alpha=c(0,1),
             size = 2.5,
             shape = 17,
             show.legend = F) +
  geom_line(data = numDf, aes(x = time,
                              y = yval,
                              color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5,
            linetype = "dashed") +
  guides(color = F,
         shape = F) +
  scale_color_manual(values = turbo(10)[c(2,4,7,9)]) + 
  theme(text = element_text(family = "LM Roman 10"),
        axis.title = element_text(size=12),
        axis.title.y = element_blank(),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.justification = c("left", "top"),
        legend.position = c(.02,.98),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) 

p2

####################
#### Error Time ####
####################

# rm(list=ls())

setwd("~/Documents/GitHub/path_integral/results/improvedPDetectedTimeVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "start", 
                   "time", "thresh", "pintAUC", 
                   "pintDetected", "numAUC", "numDetected",
                   "statDist")
master$popalpha = 1000 * master$alpha
rm(tmp, file)

master$statDist <- statDist <- master$statDist %>% 
  gsub('[{}]', '', .) %>% 
  gsub("\\*\\^", "e", .) %>% 
  as.numeric()

summary(master$statDist)
# plot(density(master$statDist))

# master <- master %>% filter(VG != 1e-04)

e2 <- ggplot(data = master, aes(x = time, y = statDist)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  theme_bw() +
  guides(color = F) +
  xlab("Time\n(Genomic Units)") + 
  ylab("Statistical Distance") +
  scale_y_log10() + 
  scale_color_manual(values = turbo(10)[c(2,4,7,9)]) + 
  theme(axis.title = element_text(size=12),
        axis.title.y = element_blank(),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.justification = c("left", "top"),
        legend.position = c(.02,.98),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())  + 
  annotation_logticks(sides = "l") 


e2 

##########################
#### pDetection Start ####
##########################

# rm(list=ls())

setwd("~/Documents/GitHub/path_integral/results/improvedPDetectedStartVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "start", 
                   "time", "thresh", "pintAUC", 
                   "pintDetected", "numAUC", "numDetected",
                   "statDist")
master$popalpha = 1000 * master$alpha
rm(tmp, file)

master$statDist <- statDist <- master$statDist %>% 
  gsub('[{}]', '', .) %>% 
  gsub("\\*\\^", "e", .) %>% 
  as.numeric()

# pintDf <- master %>% filter(statDist < 0.05)
# 
# numDf <- master %>% filter(VG == 1e-4,
#                            time >= 0.15)
# numDf$yval <- c(numDf[1,]$pintDetected,
#                 numDf[2,]$numDetected)

p3 <- ggplot(data = master, aes(x = start, y = pintDetected)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  geom_hline(yintercept=0.01, 
             linetype="dashed",
             size  = 0.75) + 
  theme_bw() +
  scale_x_continuous("Starting frequency", 
                     breaks = c(0.025,0.05,0.1,0.15,0.20), 
                     labels = c(0.025,0.05,0.1,0.15,0.20),
                     expand = c(0.05,0)) +
  scale_y_continuous("P(detected)",
                     expand = c(0,0),
                     limits = c(0,max(master$pintDetected) * 1.05)) +
  guides(color = F) + 
  scale_color_manual(values = turbo(10)[c(2,4,7,9)]) +
  theme(text = element_text(family = "LM Roman 10"),
        axis.title = element_text(size=12),
        axis.title.y = element_blank(),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.justification = c("left", "top"),
        legend.position = c(.02,.98),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) 


p3

#####################
#### Error Start ####
#####################

# rm(list=ls())

setwd("~/Documents/GitHub/path_integral/results/improvedPDetectedStartVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "start", 
                   "time", "thresh", "pintAUC", 
                   "pintDetected", "numAUC", "numDetected",
                   "statDist")
master$popalpha = 1000 * master$alpha
rm(tmp, file)

master$statDist <- statDist <- master$statDist %>% 
  gsub('[{}]', '', .) %>% 
  gsub("\\*\\^", "e", .) %>% 
  as.numeric()

summary(master$statDist)
# plot(density(master$statDist))

# master <- master %>% filter(VG != 1e-04)

e3 <- ggplot(data = master, aes(x = start, y = statDist)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  theme_bw() +
  guides(color=F) +
  scale_x_continuous("Starting Frequency", 
                     breaks = c(0.025,0.05,0.1,0.15,0.20), 
                     labels = c(0.025,0.05,0.1,0.15,0.20)) +
  ylab("Statistical Distance") +
  scale_y_log10() + 
  scale_color_manual(values = turbo(10)[c(2,4,7,9)]) + 
  theme(axis.title = element_text(size=12),
        axis.title.y = element_blank(),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.justification = c("left", "top"),
        legend.position = c(.02,.98),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())  + 
  annotation_logticks(sides = "l") 

e3 

########################
#### pDetected Main ####
########################

p1 + p2 + p3 + 
  plot_annotation(tag_levels = 'A')  & 
  theme(plot.tag = element_text(size = 12))

#########################
#### pDetected Error ####
#########################

e1 + e2 + e3 + 
  plot_annotation(tag_levels = 'A')  & 
  theme(plot.tag = element_text(size = 12))

############################
#### Convergence 20 Gen ####
############################

rm(list = ls())

setwd("/media/nathan/T7/path_integral/comparison20gen")

master <- data.frame()
count <- 0
for(file in list.files()){
  count <- count + 1
  print(paste(count, ":", file))
  
  param <- strsplit(file, "_")[[1]]
  selcoef <- param[1] %>% gsub("2na", "", .) %>%
    as.numeric()
  geg <- param[2] %>% gsub("geg", "", .) %>%
    as.numeric()
  time <- param[3] %>% gsub("time", "", .) %>%
    as.numeric()
  end <- param[4] %>% gsub("end.csv", "", .) %>%
    as.numeric()
  
  tmp <- fread(file)
  
  if(nrow(tmp) == 0){
    tmp <- names(tmp)
  }
  
  tmp <- unlist(tmp) %>%
    gsub('\"', "", .) %>%
    gsub('[{}]', '', .) %>%
    gsub("\\*\\^", "e", .) %>%
    as.numeric()
  df <- data.frame(dens = tmp,
                   k = c(0:5, "num"),
                   time = time,
                   geg = paste("geg",geg,sep=""),
                   end = end, 
                   selcoef = selcoef)
  
  if(sum(is.na(df)) > 0) break()
  
  master <- dplyr::bind_rows(master, df)
}

setwd("/media/nathan/T7/path_integral/MC/data")
tmp <- fread("2na1_vg-3_gen20_1_transMats.csv")
tmp = tmp$V1
df <- data.frame(dens = tmp,
                 k = "MC",
                 time = 0.02,
                 geg = rep(paste("geg", c(10,50), sep = ""), 
                           each = length(tmp)),
                 end = 1:99/100,
                 selcoef = 1)
master <- dplyr::bind_rows(master, df)
tmp <- fread("2na10_vg-3_gen20-50-150_1_transMats.csv")
tmp = tmp$V1
df <- data.frame(dens = tmp,
                 k = "MC",
                 time = 0.02,
                 geg = rep(paste("geg", c(10,50), sep = ""), 
                           each = length(tmp)),
                 end = 1:99/100,
                 selcoef = 10)
master <- dplyr::bind_rows(master, df)

my_labeller = as_labeller(
  c("1" = "2 * N[e] * `\u03b1 \u039b / W = 1`",
    "10" = "2 * N[e] * `\u03b1 \u039b / W = 10`",
    "geg10" = "m[max]==10",
    "geg50" = "m[max]==50"),
  default = label_parsed
)

# tmp <- master %>% filter(k %in% c("1"))

ggplot(master, aes(x = end, y = dens)) +
  geom_hline(yintercept=0, 
             linetype="dashed", size  = 0.75) +
  geom_line(aes(color = k, 
                linewidth = k, 
                linetype = k),
            alpha = 0.6,
            # linewidth = 1.5
            ) + 
  facet_grid(rows=vars(geg),
             cols=vars(selcoef),
             labeller = my_labeller) +
  coord_cartesian(ylim=c(-2,10)) +
  theme_bw() + 
  scale_color_manual(values = turbo(11)[c(2,3,4,6,7,8,9,10)]) + 
  scale_linetype_manual(values = c(rep("solid",6), "11", "33")) +
  scale_linewidth_manual(values = c(rep(1.5, 6), rep(1,2))) + 
  guides(color=guide_legend(title = bquote(k[max]),
                            override.aes = list(linewidth = 0.75),
                            ncol = 2),
         linetype = "none",
         linewidth = "none") +
  xlab("Ending Frequency") + 
  ylab("Density") +
  labs(title = "20 Generations") +
  theme(axis.title = element_text(size=12),
      title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 11),
      legend.justification = c("right", "top"),
      legend.position = c(.48,.98),
      legend.box.background = element_rect(colour = "black"),
      legend.spacing.y = unit(0, 'cm'),
      legend.key.size = unit(0.8, "line"),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank()) 

##############################
#### Convergence Alpha VG ####
##############################

rm(list = ls())

setwd("/media/nathan/T7/path_integral/comparisonAlphaVG")

master <- data.frame()
count <- 0
for(file in list.files()){
  count <- count + 1
  print(paste(count, ":", file))
  
  param <- strsplit(file, "_")[[1]]
  selcoef <- param[1] %>% gsub("2na", "", .) %>%
    as.numeric()
  geg <- param[2] %>% gsub("geg", "", .) %>%
    as.numeric()
  vg <- param[3] %>% gsub("VG", "", .) %>%
    as.numeric()
  end <- param[4] %>% gsub("end.csv", "", .) %>%
    as.numeric()
  
  tmp <- fread(file)
  
  if(nrow(tmp) == 0){
    tmp <- names(tmp)
  }
  
  tmp <- unlist(tmp) %>%
    gsub('\"', "", .) %>%
    gsub('[{}]', '', .) %>%
    gsub("\\*\\^", "e", .) %>%
    as.numeric()
  df <- data.frame(dens = tmp,
                   k = c(0:5, "num"),
                   vg = vg,
                   geg = paste("geg",geg,sep=""),
                   end = end, 
                   selcoef = selcoef)
  
  if(sum(is.na(df)) > 0) break()
  
  master <- dplyr::bind_rows(master, df)
}

setwd("/media/nathan/T7/path_integral/MC/data")
tmp <- fread("2na10_vg-2_gen200_1_transMats.csv")
tmp = tmp$V1
df <- data.frame(dens = tmp,
                 k = "MC",
                 vg = 1e-2,
                 geg = "geg50",
                 end = 1:99/100,
                 selcoef = 10)
master <- dplyr::bind_rows(master, df)
tmp <- fread("2na10_vg-4_gen200_1_transMats.csv")
tmp = tmp$V1
df <- data.frame(dens = tmp,
                 k = "MC",
                 vg = 1e-4,
                 geg = "geg50",
                 end = 1:99/100,
                 selcoef = 10)
master <- dplyr::bind_rows(master, df)
tmp <- fread("2na1_vg-2_gen200_1_transMats.csv")
tmp = tmp$V1
df <- data.frame(dens = tmp,
                 k = "MC",
                 vg = 1e-2,
                 geg = "geg50",
                 end = 1:99/100,
                 selcoef = 1)
master <- dplyr::bind_rows(master, df)
tmp <- fread("2na1_vg-4_gen200_1_transMats.csv")
tmp = tmp$V1
df <- data.frame(dens = tmp,
                 k = "MC",
                 vg = 1e-4,
                 geg = "geg50",
                 end = 1:99/100,
                 selcoef = 1)
master <- dplyr::bind_rows(master, df)

my_labeller = as_labeller(
  c("1" = "2 * N[e] * `\u03b1 \u039b / W = 1`",
    "10" = "2 * N[e] * `\u03b1 \u039b / W = 10`",
    "1e-04" = "V[G] == 10^-4",
    "0.01" = "V[G] == 10^-2"),
  default = label_parsed
)

ggplot(master, aes(x = end, y = dens)) +
  geom_hline(yintercept=0, 
             linetype="dashed", size  = 0.75) +
  geom_line(aes(color = k,
                linetype = k,
                linewidth = k),
            alpha = 0.6,
            # linewidth = 1.5
            ) + 
  facet_grid(rows=vars(vg),
             cols=vars(selcoef),
             labeller = my_labeller) +
  coord_cartesian(ylim=c(-2,10)) +
  theme_bw() + 
  scale_color_manual(values = turbo(11)[c(2,3,4,6,7,8,9,10)]) + 
  scale_linetype_manual(values = c(rep("solid",6), "11", "33")) +
  scale_linewidth_manual(values = c(rep(1.5, 6), rep(1,2))) + 
  guides(color=guide_legend(title = bquote(k[max]),
                            override.aes = list(linewidth = 0.75),
                            ncol = 2),
         linewidth = "none",
         linetype = "none") +
  xlab("Ending Frequency") + 
  ylab("Density") +
  labs(title = "200 Generations") +
  theme(axis.title = element_text(size=12),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 11),
        legend.justification = c("right", "top"),
        legend.position = c(.48,.98),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) 

############################
#### Convergence 2Na 5 ####
############################

rm(list = ls())

setwd("/media/nathan/T7/path_integral/comparison2na5")

master <- data.frame()
count <- 0
for(file in list.files()){
  count <- count + 1
  print(paste(count, ":", file))
  
  param <- strsplit(file, "_")[[1]]
  selcoef <- param[1] %>% gsub("2na", "", .) %>%
    as.numeric()
  geg <- param[2] %>% gsub("geg", "", .) %>%
    as.numeric()
  if(geg == 5) geg <- 10
  time <- param[3] %>% gsub("time", "", .) %>%
    as.numeric()
  end <- param[4] %>% gsub("end.csv", "", .) %>%
    as.numeric()
  
  tmp <- fread(file)
  
  if(nrow(tmp) == 0){
    tmp <- names(tmp)
  }
  
  tmp <- unlist(tmp) %>%
    gsub('\"', "", .) %>%
    gsub('[{}]', '', .) %>%
    gsub("\\*\\^", "e", .) %>%
    as.numeric()
  df <- data.frame(dens = tmp,
                   k = c(0:5, "num"),
                   time = time,
                   geg = paste("geg",geg,sep=""),
                   end = end, 
                   selcoef = selcoef)
  
  if(sum(is.na(df)) > 0) break()
  
  master <- dplyr::bind_rows(master, df)
}

my_labeller = as_labeller(
  c("0.25" = "t==0.25",
    "0.75" = "t==0.75",
    "geg10" = "m[max]==10",
    "geg50" = "m[max]==50"),
  default = label_parsed
)

ggplot(master, aes(x = end, y = dens)) +
  geom_hline(yintercept=0, 
             linetype="dashed",
             size  = 0.75) +
  geom_line(aes(color = k,
                linetype = k,
                linewidth = k),
            alpha = 0.6,
            linewidth = 1.5) + 
  facet_grid(rows=vars(geg),
             cols=vars(time),
             labeller = my_labeller) +
  coord_cartesian(ylim=c(-0.25,2)) +
  theme_bw() + 
  scale_color_manual(values = turbo(11)[c(2,3,4,6,7,8,10)]) + 
  scale_linetype_manual(values = c(rep("solid",6), "33")) +
  scale_linewidth_manual(values = c(rep(1.5, 6), rep(1,1))) + 
  guides(color=guide_legend(title = bquote(k[max]),
                            override.aes = list(linewidth = 0.75),
                            ncol = 2),
         linetype = "none",
         linewidth = "none") +
  xlab("Ending Frequency") + 
  ylab("Density") +
  labs(title = bquote(2 * N[e] * "\u03b1 \u039b / W = 5")) +
  theme(axis.title = element_text(size=12),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 11),
        legend.justification = c("right", "top"),
        legend.position = c(.48,.98),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) 

############################
#### Convergence 2Na 10 ####
############################

rm(list = ls())

setwd("/media/nathan/T7/path_integral/comparison2na10")

list.files()

master <- data.frame()
count <- 0
for(file in list.files()){
  count <- count + 1
  print(paste(count, ":", file))
  
  param <- strsplit(file, "_")[[1]]
  geg <- param[2] %>% gsub("geg", "", .) %>%
    as.numeric()
  time <- param[3] %>% gsub("time", "", .) %>%
    as.numeric()
  end <- param[4] %>% gsub("end.csv", "", .) %>%
    as.numeric()
  
  tmp <- fread(file)
  
  if(nrow(tmp) == 0){
    tmp <- names(tmp)
  }
  
  tmp <- unlist(tmp) %>%
    gsub('\"', "", .) %>%
    gsub('[{}]', '', .) %>%
    gsub("\\*\\^", "e", .) %>%
    as.numeric()
  df <- data.frame(dens = tmp,
                   k = c(0:5, "num"),
                   time = time,
                   geg = geg,
                   end = end)
  
  if(sum(is.na(df)) > 0) break()
  
  master <- dplyr::bind_rows(master, df)
}

setwd("/media/nathan/T7/path_integral/MC/data")
tmp <- fread("2na10_vg-3_gen20-50-150_2_transMats.csv")
tmp = tmp$V1
df <- data.frame(dens = tmp,
                 k = "MC",
                 time = 0.05,
                 geg = rep(c(10,50), 
                           each = length(tmp)),
                 end = 1:99/100)
master <- dplyr::bind_rows(master, df)
tmp <- fread("2na10_vg-3_gen20-50-150_3_transMats.csv")
tmp = tmp$V1
df <- data.frame(dens = tmp,
                 k = "MC",
                 time = 0.15,
                 geg = rep(c(10,50), 
                           each = length(tmp)),
                 end = 1:99/100)
master <- dplyr::bind_rows(master, df)

my_labeller = as_labeller(
  c("0.05" = "t==0.05",
    "0.15" = "t==0.15",
    "10" = "m[max]==10",
    "50" = "m[max]==50"), 
  default = label_parsed
)

ggplot(master, aes(x = end, y = dens)) +
  geom_hline(yintercept=0, 
             linetype="dashed",
             size  = 0.75) +
  geom_line(aes(color = k,
                linetype = k,
                linewidth = k),
            alpha = 0.6,
            # linewidth = 1.5
            ) + 
  facet_grid(rows=vars(geg),
             cols=vars(time),
             labeller = my_labeller) + 
  coord_cartesian(ylim=c(-2,6)) +
  theme_bw() + 
  scale_color_manual(values = turbo(11)[c(2,3,4,6,7,8,9,10)]) + 
  scale_linetype_manual(values = c(rep("solid",6), "11", "33")) +
  scale_linewidth_manual(values = c(rep(1.5, 6), rep(1,2))) + 
  guides(color=guide_legend(title = bquote(k[max]),
                            override.aes = list(linewidth = 0.75),
                            ncol = 2),
         linetype = "none",
         linewidth = "none") +
  xlab("Ending Frequency") + 
  ylab("Density") +
  labs(title = bquote(2 * N[e] * "\u03b1 \u039b / W = 10")) +
  theme(axis.title = element_text(size=12),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 11),
        legend.justification = c("right", "top"),
        legend.position = c(.47,.48),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) 

#######################
#### pDetection Ne ####
#######################

rm(list=ls())

setwd("/media/nathan/T7/path_integral/improvedPDetectedAlphaNE")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("selCoef", "Ne", "time", "start", 
                   "thresh", "pintAUC", "pintdetected",
                   "numAUC", "numdetected", "statDist")
master$popalpha = 2 * master$Ne * master$selCoef
rm(tmp)

p1 <- ggplot(data = master, aes(x = Ne, y = pintdetected)) + 
  geom_line(aes(color = as.factor(selCoef)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(selCoef)),
             alpha=1,
             size = 2.5) +
  geom_hline(yintercept=0.01,
             linetype="dashed",
             color = turbo(11)[11], size  = 0.75) +
  theme_bw() +
  scale_x_continuous(bquote(N[e]),
                     breaks = c(100,200, 250, 500, 1000)) +
  scale_y_continuous("Probability Detected\n(Q)") +
  guides(color=guide_legend(title = "\u03b1",
                            override.aes = list(linewidth = 0.75))) +
  scale_color_manual(values = turbo(10)[c(2,4,6,7,9)]) +     
  theme(axis.title = element_text(size=12),
        axis.title.x = element_blank(),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 11),
        legend.justification = c("left", "top"),
        legend.position = c(0.02,0.98),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        # legend.key.size = unit(0.8, "line"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) 

p1
##################
#### Error Ne ####
##################

# rm(list=ls())

setwd("/media/nathan/T7/path_integral/improvedPDetectedAlphaNE")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("selCoef", "Ne", "time", "start", 
                   "thresh", "pintAUC", "pintdetected",
                   "numAUC", "numdetected", "statDist")
master$popalpha = 2 * master$Ne * master$selCoef
rm(tmp)

master$statDist <- statDist <- master$statDist %>% 
  gsub('[{}]', '', .) %>% 
  gsub("\\*\\^", "e", .) %>% 
  as.numeric()

p2 <- ggplot(data = master, aes(x = Ne, y = statDist)) + 
  geom_line(aes(color = as.factor(selCoef)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(selCoef)),
             alpha=1,
             size = 2.5) +
  theme_bw() +
  guides(color=F) +
  scale_x_continuous(bquote(N[e]),
                     breaks = c(100,200, 250, 500, 1000)) +
  scale_y_log10("Statistical Distance") +
  scale_color_manual(values = turbo(10)[c(2,4,6,7,9)]) +  
  theme(axis.title = element_text(size=12),
        axis.title.x = element_blank(),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 11),
        legend.justification = c("left", "top"),
        legend.position = c(0.02,0.98),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        # legend.key.size = unit(0.8, "line"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  annotation_logticks(sides = "l") 

#####################
#### Ne combined ####
#####################

p3 <- ggplot(data.frame(l = "x", 
                        x = 1, y = 1)) +
  geom_text(aes(x, y, label = deparse(bquote(N[e]))), 
            # angle = 90,
            size = 5 / 14 * 12,
            parse = T) + 
  theme_void() +
  coord_cartesian(clip = "off")

(p1 | p2) / p3 + 
  # plot_annotation(tag_levels = 'A')  & 
  # theme(plot.tag = element_text(size = 12)) & 
  plot_layout(heights = c(25,1))


#############
#### RFS ####
#############

rm(list=ls())

setwd("/media/nathan/T7/path_integral/improvedPDetectedAlphaNE")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("selCoef", "Ne", "time", "start", 
                   "thresh", "pintAUC", "pintdetected", 
                   "numAUC", "numdetected", "statDist")
master$popalpha = 2 * master$Ne * master$selCoef
rm(tmp)

# master <- master %>% filter(selCoef == 0.01)
master$reps = 5000 / master$Ne

df <- data.table()
for(i in 1:nrow(master)){
  tmp <- master[i,]
  for(j in 0:tmp$reps){
    foo <- data.table(Ne = tmp$Ne,
                      reps = paste(tmp$reps),
                      selCoef = tmp$selCoef,
                      reps_number = tmp$reps,
                      bin = j,
                      prob = choose(tmp$reps, j) * 
                        tmp$pintdetected^j * 
                        (1 - tmp$pintdetected)^(tmp$reps - j))
    df <- dplyr::bind_rows(df, foo)
  }
}
rm(foo, master, tmp, file, i, j)

df <- df %>% filter(Ne != 250) %>%
  mutate(reps = as.factor(reps),
         selCoef = as.factor(selCoef))
df$reps <- factor(df$reps,
                  levels = levels(df$reps)[c(3,1,2,4)])


for(i in 1:4){
  d <- df %>% filter(reps == levels(df$reps)[i])
  d <- d[order(d$selCoef),]
  
  p <- ggplot(d, aes(y=prob, x=bin, color = selCoef,
                     group = selCoef)) + 
    geom_line(alpha = 0.6,
              size = 1.5) + 
    geom_point(alpha=1,
               size = 2.5) +
    theme_bw() + 
    scale_x_continuous("Number of Replicates Detected",
                       expand = c(0.05,0)) +
    scale_y_continuous("Probability",
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     limits = c(0,1), 
                     expand = c(0,0)) + 
    scale_color_manual(values = turbo(10)[c(2,4,6,7,9)],
                       labels = levels(d$selCoef)) +  
    guides(color = F) + 
    theme(text = element_text(family = "LM Roman 10"),
          axis.title = element_blank(),
          title = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.justification = c("right", "top"),
          legend.position = c(.98,.98),
          legend.box.background = element_rect(colour = "black"),
          legend.spacing.y = unit(0, 'cm'),
          legend.key.size = unit(0.8, "line"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank()) 
  
  if(i == 1) p1 <- p + labs(title = bquote(5 ~ `replicates,` ~ italic(N[e]) == 1000))
  if(i == 2) p2 <- p + scale_x_continuous(breaks =  seq(0,10,by = 2)) + 
    labs(title = bquote(10 ~ "replicates," ~ italic(N[e]) == ~ "500")) + 
    guides(color=guide_legend(title = expression(italic(alpha)),
                              override.aes = list(alpha=1,
                                                  linewidth = 0.75,
                                                  size = 1.25))) 
  if(i == 3) p3 <- p + labs(title = bquote(25 ~ `replicates,` ~ italic(N[e]) == 200)) + 
      scale_x_continuous(breaks = seq(0,10,by = 2), limits = c(0,10)) 
  if(i == 4) p4 <- p + labs(title = bquote(50 ~ `replicates,` ~ italic(N[e]) == 100)) +
    scale_x_continuous(breaks = seq(0,10,by = 2), limits = c(0,10)) 
}

p5 <- ggplot(data.frame(l = "Number of replicates detected", 
                        x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), 
            # angle = 90,
            size = 5 / 14 * 12,
            family = "LM Roman 10") + 
  theme_void() +
  coord_cartesian(clip = "off")

p6 <- ggplot(data.frame(l = "Probability", 
                        x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), 
            angle = 90,
            size = 5 / 14 * 12,
            family = "LM Roman 10") + 
  theme_void() +
  coord_cartesian(clip = "off")

layout <- "
ABC
ADE
#FF
"

p6 + p1 + p2 + p3 + p4 + p5 + plot_layout(design = layout,
                                          widths = c(2,25,25),
                                          heights = c(25,25,1))


#############
#### RFS ####
#############

rm(list=ls())

setwd("/media/nathan/T7/path_integral/improvedPDetectedAlphaNE")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("selCoef", "Ne", "time", "start", 
                   "thresh", "pintAUC", "pintdetected", 
                   "numAUC", "numdetected", "statDist")
master$popalpha = 2 * master$Ne * master$selCoef
rm(tmp)

master <- master %>% filter(selCoef == 1e-2)

df <- data.table()
for(i in 1:nrow(master)){
  tmp <- master[i,]
  for(reps in c(5, 10, 15, 20)){
    for(j in 0:reps){
      foo <- data.table(reps = paste(reps, " Replicates"),
                        reps_number = reps,
                        bin = j,
                        selCoef = tmp$selCoef,
                        Ne = tmp$Ne,
                        prob = choose(reps, j) * 
                          tmp$pintdetected^j * 
                          (1 - tmp$pintdetected)^(reps - j))
      df <- dplyr::bind_rows(df, foo)
    }
  }
}
rm(foo, tmp, file, i, j)

df <- df %>% filter(Ne != 250) %>%
  mutate(reps = as.factor(reps),
         Ne = as.factor(Ne))
df$reps <- factor(df$reps,
                  levels = levels(df$reps)[c(4,1,2,3)])

 
for(i in 1:4){
  d <- df %>% filter(reps == levels(df$reps)[i])
  d <- d[order(d$Ne),]
  
  p <- ggplot(d, aes(y=prob, x=bin, color = Ne,
                               group = Ne)) + 
    geom_line(alpha = 0.6,
              size = 1.5) + 
    geom_point(alpha=1,
               size = 2.5) +
    theme_bw() + 
    scale_x_continuous("Number of Replicates Detected",
                       expand = c(0.05,0)) +
    scale_y_continuous("Probability",
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       limits = c(0,1), 
                       expand = c(0,0)) + 
    scale_color_manual(values = turbo(10)[c(2,4,7,9)],
                       labels = levels(d$Ne)) +  
    guides(color = F) + 
    theme(axis.title = element_blank(),
          title = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.justification = c("right", "top"),
          legend.position = c(.98,.98),
          legend.box.background = element_rect(colour = "black"),
          legend.spacing.y = unit(0, 'cm'),
          legend.key.size = unit(0.8, "line"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank()) 
  
  if(i == 1) p1 <- p + labs(title = unique(d$reps))
  if(i == 2) p2 <- p + scale_x_continuous(breaks =  seq(0,10,by = 2)) + 
    labs(title = unique(d$reps)) + 
    guides(color=guide_legend(title = bquote(italic(N[e])),
                              override.aes = list(alpha=1,
                                                  linewidth = 0.75,
                                                  size = 1.25))) 
  if(i == 3) p3 <- p + labs(title = unique(d$reps)) + 
    scale_x_continuous(breaks = seq(0,10,by = 2), limits = c(0,10)) 
  if(i == 4) p4 <- p + labs(title = unique(d$reps)) +
    scale_x_continuous(breaks = seq(0,10,by = 2), limits = c(0,10)) 
}

p5 <- ggplot(data.frame(l = "Number of Replicates Detected", 
                        x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), 
            # angle = 90,
            size = 5 / 14 * 12) + 
  theme_void() +
  coord_cartesian(clip = "off")

p6 <- ggplot(data.frame(l = "Probability", 
                        x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), 
            angle = 90,
            size = 5 / 14 * 12) + 
  theme_void() +
  coord_cartesian(clip = "off")

layout <- "
ABC
ADE
#FF
"

p6 + p1 + p2 + p3 + p4 + p5 + plot_layout(design = layout,
                                          widths = c(2,25,25),
                                          heights = c(25,25,1))

######################
#### model figure ####
######################

rm(list=ls())

x <- seq(from = -2, to = 3, length.out = 500)

phenodist <- dnorm(x, sd = sqrt(10^-2))
phenodist <- phenodist / max(phenodist)

olddist <- dnorm(x, sd = 1)
olddist <- olddist / max(olddist)

newdist <- dnorm(x, mean = 1, sd = 1)
newdist <- newdist / max(newdist)

df <- dplyr::bind_rows(data.table(x = x,
                                  y = phenodist,
                                  class = "pheno",
                                  cols = 1),
                       data.table(x = x,
                                  y = olddist,
                                  class = "old",
                                  cols = 2),
                       data.table(x = x,
                                  y = newdist,
                                  class = "new",
                                  cols = 2))

areadf <- data.table(x = x,
                     y = phenodist,
                     class = "pheno",
                     cols = 1)

the_colors <- c(rgb(0.831964, 0.810543, 0.372854),
                rgb(0.35082, 0.595178, 0.853742))

ggplot(df, aes(x = x, y = y)) + 
  geom_line(aes(color = class,
                linetype = class,
                group = class),
            linewidth = 2,
            alpha = 0.75) +
  scale_color_manual(values = c(rep(the_colors[1],2), the_colors[2]),
                     labels=c("Shifted fitness\n function",
                              "Initial fitness\n function",
                              "Initial trait\n distribution"),
                     name = "") +
  scale_linetype_manual(values = c("dotted", rep("solid",2)),
                        labels=c("Shifted fitness\n function",
                                 "Initial fitness\n function",
                                 "Initial trait\n distribution"),
                        name = "") +
  theme_bw() +
  theme(text = element_text(family = "LM Roman 10"),
        axis.title = element_text(size=12),
        axis.text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank()) + 
  scale_x_continuous(breaks = c(0,1),
                     labels = c("0" = "Old\noptimum", 
                                "1" = "New\noptimum"),
                     name = "Trait value") +
  guides(linetype = guide_legend(override.aes = list(linewidth = 2,
                                                     linetype = "solid"),
                                 nrow = 1)) +
  geom_area(data = areadf, 
            fill = the_colors[2],
            alpha = 0.5)

##################
#### Error VG ####
##################

setwd("~/Documents/GitHub/path_integral/results/numErrorAlpha")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "pintAUC", "numAUC", "statDist")
master$popalpha = 1000 * master$alpha
rm(tmp, file)

master$statDist <- statDist <- master$statDist %>% 
  gsub('[{}]', '', .) %>% 
  gsub("\\*\\^", "e", .) %>% 
  as.numeric()

p1 <- ggplot(data = master, aes(x = popalpha, y = statDist)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  theme_bw() +
  guides(color=F) +
  scale_x_continuous(bquote(2 * N[e] * "\u03b1 \u039b / W"), 
                     breaks = c(0, 1, 5, 10, 15, 20), 
                     labels = c(0, 1, 5, 10, 15, 20)) +
  ylab("Statistical Distance") +
  scale_y_log10() + 
  scale_color_manual(values = turbo(10)[c(2,4,7,9)]) + 
  labs(tag = "A") + 
  theme(axis.title = element_text(size=12),
        axis.title.y = element_blank(),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification = c("left", "top"),
        legend.position = c(.01,.99),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.tag = element_text(size = 12))  + 
  annotation_logticks(sides = "l") 

p1

####################
#### Error Time ####
####################

# rm(list = ls())

setwd("~/Documents/GitHub/path_integral/results/numErrorTime")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("time", "VG", "pintAUC", "numAUC", "statDist")
rm(tmp, file)

master$statDist <- statDist <- master$statDist %>% 
  gsub('[{}]', '', .) %>% 
  gsub("\\*\\^", "e", .) %>% 
  as.numeric()

p2 <- ggplot(data = master, aes(x = time, y = statDist)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  theme_bw() +
  guides(color=F) +
  xlab("Time (Genomic Units)") +
  ylab("Statistical Distance") +
  scale_y_log10() + 
  scale_color_manual(values = turbo(10)[c(2,4,7,9)]) + 
  labs(tag = "B") + 
  theme(axis.title = element_text(size=12),
        axis.title.y = element_blank(),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification = c("left", "top"),
        legend.position = c(.01,.99),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.tag = element_text(size = 12))  + 
  annotation_logticks(sides = "l") 

p2
#################
#### Error k ####
#################

# rm(list = ls())

setwd("~/Documents/GitHub/path_integral/results/numErrorK")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("k", "VG", "pintAUC", "numAUC", "statDist")
rm(tmp, file)

master$statDist <- statDist <- master$statDist %>% 
  gsub('[{}]', '', .) %>% 
  gsub("\\*\\^", "e", .) %>% 
  as.numeric()


p3 <- ggplot(data = master, aes(x = k, y = statDist)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  theme_bw() +
  guides(color=guide_legend(title = bquote(V[G]),
                            override.aes = list(alpha=1,
                                                size = 1.25,
                                                linewidth = 0.75),
                            nrow = 2)) +
  scale_x_continuous(bquote(k[max])) +
  scale_y_log10("Statistical Distance") + 
  scale_color_manual(values = turbo(10)[c(2,4,7,9)]) + 
  labs(tag = "C") + 
  theme(axis.title = element_text(size=12),
        axis.title.y = element_blank(),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification = c("left", "bottom"),
        legend.position = c(.01,.01),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.tag = element_text(size = 12))  + 
  annotation_logticks(sides = "l") 

p3

#####################
#### Error Start ####
#####################

# rm(list = ls())

setwd("~/Documents/GitHub/path_integral/results/numErrorStart")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("start", "VG", "pintAUC", "numAUC", "statDist")
rm(tmp, file)

master$statDist <- statDist <- master$statDist %>% 
  gsub('[{}]', '', .) %>% 
  gsub("\\*\\^", "e", .) %>% 
  as.numeric()


p4 <- ggplot(data = master, aes(x = start, y = statDist)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  theme_bw() +
  guides(color=F) +
  scale_x_continuous("Starting Frequency",
                     breaks = (1:9 / 10)) + 
  scale_y_log10("Statistical Distance") + 
  scale_color_manual(values = turbo(10)[c(2,4,7,9)]) + 
  labs(tag = "D") + 
  theme(axis.title = element_text(size=12),
        axis.title.y = element_blank(),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification = c("left", "top"),
        legend.position = c(.01,.99),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.tag = element_text(size = 12))  + 
  annotation_logticks(sides = "l") 

p4

########################
#### alt error main ####
########################

p5 <- ggplot(data.frame(l = "Statistical Distance", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90,
            size = 5 / 14 * 12) + 
  theme_void() +
  coord_cartesian(clip = "off")

p5 + ((p1 | p2) / (p3 | p4)) + plot_layout(widths = c(1, 25))


##################################
########### Deprecated ########### 
##################################
#######################
#### pDetection VG ####
#######################

rm(list=ls())

setwd("~/Documents/GitHub/path_integral/results/pDetectionAlphaVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "start", "thresh", "totalP", "Pdetected")
master$popalpha = 1000 * master$alpha
rm(tmp)

ggplot(data = master, aes(x = popalpha, y = Pdetected)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  geom_hline(yintercept=0.01, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) + 
  geom_vline(xintercept=1, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) + 
  theme_bw() +
  guides(color=guide_legend(title = "Genetic Variance",
                            override.aes = list(alpha=1))) +
  scale_x_continuous("Population Scaled Selection Coefficient", 
                     breaks = c(0,0.5,1,5,10), 
                     limits = c(0,10),
                     labels = c(0,0.5,1,5,10)) +
  scale_y_continuous("P(detected)",
                     breaks = c(0.01,0.02,0.04,0.06,0.08)) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
                                   colour = c(rep("black",2),"red",rep("black",2)))) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18))

#########################
#### pDetection Time ####
#########################

rm(list = ls())
setwd("~/Documents/GitHub/path_integral/results/pDetectionAlphaTime")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "time", "start", "thresh", "totalP", "Pdetected")
master$popalpha = 1000 * master$alpha
rm(tmp)

ggplot(data = master, aes(x = popalpha, y = Pdetected)) + 
  geom_line(aes(color = as.factor(time)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(time)),
             alpha=1,
             size = 2.5) +
  geom_hline(yintercept=0.01, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) + 
  geom_vline(xintercept=1, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) + 
  theme_bw() +
  guides(color=guide_legend(title = "Time (Genomic Units)",
                            override.aes = list(alpha=1))) +
  # scale_x_break(c(1,2), scales = 0.9) 
  scale_x_continuous("Population Scaled Selection Coefficient", 
                     breaks = c(0,0.5,1,5,10), 
                     limits = c(0,10),
                     labels = c(0,0.5,1,5,10)) +
  scale_y_continuous("P(detected)") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
                                   colour = c(rep("black",2),"red",rep("black",2)))) + 
  scale_color_viridis(discrete = T)

##########################
#### pDetection Start ####
##########################

rm(list = ls())

setwd("~/Documents/GitHub/path_integral/results/pDetectionAlphaStart")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "start", "thresh", "totalP", "Pdetected")
master$popalpha = 1000 * master$alpha
rm(tmp)

ggplot(data = master, aes(x = popalpha, y = Pdetected)) + 
  geom_line(aes(color = as.factor(start)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(start)),
             alpha=1,
             size = 2.5) +
  geom_hline(yintercept=0.01, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) + 
  geom_vline(xintercept=1, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) + 
  theme_bw() +
  guides(color=guide_legend(title = "Starting freq",
                            override.aes = list(alpha=1))) +
  # scale_x_break(c(1,2), scales = 0.9) 
  scale_x_continuous("Population Scaled Selection Coefficient", 
                     breaks = c(0,0.5,1,5,10), 
                     limits = c(0,10),
                     labels = c(0,0.5,1,5,10)) +
  scale_y_continuous("P(detected)") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
                                   colour = c(rep("black",2),"red",rep("black",2)))) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18))

##################
#### Error VG ####
##################

setwd("~/Documents/GitHub/path_integral/results/numErrorAlphaVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "pintAUC", "numAUC")
master$popalpha = 1000 * master$alpha
master$error = (master$pintAUC - master$numAUC) / master$numAUC
rm(tmp, file)

summary(master$error)
plot(density(master$error))

ggplot(data = master, aes(x = popalpha, y = error)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  theme_bw() +
  guides(color=guide_legend(title = "Genetic Variance",
                            override.aes = list(alpha=1))) +
  scale_x_continuous("Population Scaled Selection Coefficient") +
  scale_y_continuous("Error") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18)) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) 

# mster <- master %>% filter(VG<=0.001)
# ggplot(data = mster, aes(x = popalpha, y = error)) +
#   geom_line(aes(color = as.factor(VG)),
#             alpha = 0.6,
#             size = 1.5) +
#   geom_point(aes(color = as.factor(VG)),
#              alpha=1,
#              size = 2.5) +
#   theme_bw() +
#   guides(color=guide_legend(title = "Genetic Variance",
#                             override.aes = list(alpha=1))) +
#   scale_x_continuous("Population Scaled Selection Coefficient") +
#   scale_y_continuous("Absolute Error") +
#   theme(panel.grid.minor.y = element_blank(),
#         panel.grid.minor.x = element_blank()) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   scale_color_viridis(discrete = T)

####################
#### Error Time ####
####################

rm(list = ls())

setwd("~/Documents/GitHub/path_integral/results/numErrorTimeVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("time", "VG", "pintAUC", "numAUC")
master$error = (master$pintAUC - master$numAUC) / master$numAUC
rm(tmp, file)

summary(master$error)
plot(density(master$error))

ggplot(data = master, aes(x = time, y = error)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  theme_bw() +
  guides(color=guide_legend(title = "Genetic Variance",
                            override.aes = list(alpha=1))) +
  scale_x_continuous("Time (Genomic Units)") +
  scale_y_continuous("Error") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18)) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) 

# mster <- master %>% filter(VG<=0.001)
# ggplot(data = mster, aes(x = time, y = error)) +
#   geom_line(aes(color = as.factor(VG)),
#             alpha = 0.6,
#             size = 1.5) +
#   geom_point(aes(color = as.factor(VG)),
#              alpha=1,
#              size = 2.5) +
#   theme_bw() +
#   guides(color=guide_legend(title = "Genetic Variance",
#                             override.aes = list(alpha=1))) +
#   scale_x_continuous("Population Scaled Selection Coefficient") +
#   scale_y_continuous("Absolute Error") +
#   theme(panel.grid.minor.y = element_blank(),
#         panel.grid.minor.x = element_blank()) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   scale_color_viridis(discrete = T)
#################
#### Error k ####
#################

rm(list = ls())

setwd("~/Documents/GitHub/path_integral/results/numErrorKVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("k", "VG", "pintAUC", "numAUC")
master$error = (master$pintAUC - master$numAUC) / master$numAUC
rm(tmp, file)

summary(master$error)
plot(density(master$error))

ggplot(data = master, aes(x = k, y = error)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  theme_bw() +
  guides(color=guide_legend(title = "Genetic Variance",
                            override.aes = list(alpha=1))) +
  scale_x_continuous("k_max") +
  scale_y_continuous("Error") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18)) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) 

#####################
#### Error Start ####
#####################

rm(list = ls())

setwd("~/Documents/GitHub/path_integral/results/numErrorStartVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("start", "VG", "pintAUC", "numAUC")
master$error = (master$pintAUC - master$numAUC) / master$numAUC
rm(tmp, file)

summary(master$error)
plot(density(master$error))

ggplot(data = master, aes(x = start, y = error)) + 
  geom_line(aes(color = as.factor(VG)),
            alpha = 0.6,
            size = 1.5) + 
  geom_point(aes(color = as.factor(VG)),
             alpha=1,
             size = 2.5) +
  theme_bw() +
  guides(color=guide_legend(title = "Genetic Variance",
                            override.aes = list(alpha=1))) +
  scale_x_continuous("Starting Frequency",
                     breaks = (1:9 / 10)) +
  scale_y_continuous("Error") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_viridis(discrete = T) + 
  theme(axis.title = element_text(size=18)) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) 

############################
#### Convergence 20 Gen ####
############################

rm(list = ls())

setwd("~/Documents/GitHub/path_integral/results/convergence20gen")

df <- data.frame(sel=c(),
                 geg=c(),
                 fig=c())

for(sel in c(1, 5,10)){
  for(geg in c(10,30,50)){
    df <- dplyr::bind_rows(df,
                           data.frame(sel,
                                      geg,
                                      path.expand(paste("2na",sel,"_", geg,"geg","_0.02time.svg",sep=""))))
  }
}

df <- df %>% rename(fig = path.expand.paste..2na...sel..._...geg...geg...._0.02time.svg...,
                    m_max = geg,
                    PopScaledSelection = sel)

the_colors <- c(rgb(0.396811, 0.31014, 0.2041), 
                rgb(0.64274, 0.330577, 0.1540755),
                rgb(0.726732, 0.538136, 0.31593),
                rgb(0.817882, 0.7260905, 0.426991),
                rgb(0.831964, 0.810543, 0.372854),
                rgb(0.6419975, 0.7183185, 0.366907),
                rgb(0.35082, 0.595178, 0.853742))

pointdf <- data.table(x = 1,
                      y = 1,
                      cols = c(0:5, "NDSolve"),
                      lntyp = c(rep(1,6), 2)) 
ggplot(df) +
  geom_image(aes(x=1,
                 y=1,
                 image=fig),
             size=Inf) +
  facet_grid(rows=vars(m_max),
             cols=vars(PopScaledSelection),
             scales = "free",
             labeller = label_both) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position="bottom") + 
  xlab("Ending Frequency") +
  ylab("Density") + 
  scale_color_manual(limits = names(the_colors), 
                     values = the_colors) +
  geom_line(data = pointdf, aes(x = x, 
                                y = y,
                                color = cols,
                                linetype = cols),
            alpha = 0) + 
  guides(color = guide_legend(override.aes = list(alpha = 1),
                              nrow = 1)) +
  labs(color = "k_max",
       linetype = "k_max") + 
  scale_color_manual(values = the_colors, 
                     name = "k_max") +
  scale_linetype_manual(values=c(rep("solid",6), "dashed"),
                        name="k_max") + 
  theme(axis.title = element_text(size=18))



# ggplot(pointdf, aes(x = x, y = y, color = cols)) +
#   geom_line(aes(linetype = cols),
#              alpha = 0) + 
#   guides(color = guide_legend(override.aes = list(alpha = 1),
#                               nrow = 1)) +
#   theme_bw() + 
#   labs(color = "k_max",
#        linetype = "k_max") + 
#   scale_color_manual(values = the_colors, name = "k_max") +
#   scale_linetype_manual(values=c(rep("solid",6), "dashed"), name="k_max") + 
#   theme(legend.position="bottom") 
##############################
#### Convergence Alpha VG ####
##############################

rm(list = ls())
setwd("~/Documents/GitHub/path_integral/results/convergenceAlphaVG")
df <- data.frame(sel=c(),
                 genVar=c(),
                 fig=c())
for(sel in c(1, 5,10)){
  for(genVar in c("0.0001","0.001","0.01")){
    df <- dplyr::bind_rows(df,
                           data.frame(sel,
                                      genVar,
                                      path.expand(paste("2na",sel,"._50geg","_",genVar,"VG.svg",sep=""))))
  }
}
df <- df %>% rename(fig = path.expand.paste..2na...sel...._50geg...._...genVar...VG.svg...,
                    VG = genVar,
                    PopScaledSelection = sel)

the_colors <- c(rgb(0.396811, 0.31014, 0.2041), 
                rgb(0.64274, 0.330577, 0.1540755),
                rgb(0.726732, 0.538136, 0.31593),
                rgb(0.817882, 0.7260905, 0.426991),
                rgb(0.831964, 0.810543, 0.372854),
                rgb(0.6419975, 0.7183185, 0.366907),
                rgb(0.35082, 0.595178, 0.853742))

pointdf <- data.table(x = 1,
                      y = 1,
                      cols = c(0:5, "NDSolve"),
                      lntyp = c(rep(1,6), 2)) 
ggplot(df) +
  geom_image(aes(x=1,
                 y=1,
                 image=fig),
             size=Inf) +
  facet_grid(rows=vars(VG),
             cols=vars(PopScaledSelection),
             scales = "free",
             labeller = label_both) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position="bottom") +
  xlab("Ending Frequency") +
  ylab("Density") + 
  scale_color_manual(limits = names(the_colors), 
                     values = the_colors) +
  geom_line(data = pointdf, aes(x = x, 
                                y = y,
                                color = cols,
                                linetype = cols),
            alpha = 0) + 
  guides(color = guide_legend(override.aes = list(alpha = 1),
                              nrow = 1)) +
  labs(color = "k_max",
       linetype = "k_max") + 
  scale_color_manual(values = the_colors, 
                     name = "k_max") +
  scale_linetype_manual(values=c(rep("solid",6), "dashed"),
                        name="k_max") + 
  theme(axis.title = element_text(size=18))

############################
#### Convergence 2Na 5 ####
############################

rm(list = ls())
setwd("~/Documents/GitHub/path_integral/results/convergence2Na5")
list.files()
df <- data.frame(tme=c(),
                 geg=c(),
                 fig=c())
for(tme in c(0.5, 0.25, 0.75)){
  for(geg in c(10,30,50)){
    df <- dplyr::bind_rows(df,
                           data.frame(tme,
                                      geg,
                                      path.expand(paste("2na5_", geg,"geg","_",tme,"time.svg",sep=""))))
  }
}
df <- df %>% rename(fig = path.expand.paste..2na5_...geg...geg...._...tme...time.svg...,
                    m_max = geg,
                    Time = tme)

the_colors <- c(rgb(0.396811, 0.31014, 0.2041), 
                rgb(0.64274, 0.330577, 0.1540755),
                rgb(0.726732, 0.538136, 0.31593),
                rgb(0.817882, 0.7260905, 0.426991),
                rgb(0.831964, 0.810543, 0.372854),
                rgb(0.6419975, 0.7183185, 0.366907),
                rgb(0.35082, 0.595178, 0.853742))

pointdf <- data.table(x = 1,
                      y = 1,
                      cols = c(0:5, "NDSolve"),
                      lntyp = c(rep(1,6), 2)) 

ggplot(df) +
  geom_image(aes(x=1,
                 y=1,
                 image=fig),
             size=Inf) +
  facet_grid(rows=vars(m_max),
             cols=vars(Time),
             scales = "free",
             labeller = label_both) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position="bottom") +
  xlab("Ending Frequency") +
  ylab("Density") + 
  scale_color_manual(limits = names(the_colors), 
                     values = the_colors) +
  geom_line(data = pointdf, aes(x = x, 
                                y = y,
                                color = cols,
                                linetype = cols),
            alpha = 0) + 
  guides(color = guide_legend(override.aes = list(alpha = 1),
                              nrow = 1)) +
  labs(color = "k_max",
       linetype = "k_max") + 
  scale_color_manual(values = the_colors, 
                     name = "k_max") +
  scale_linetype_manual(values=c(rep("solid",6), "dashed"),
                        name="k_max") + 
  theme(axis.title = element_text(size=18))

############################
#### Convergence 2Na 10 ####
############################

rm(list = ls())
setwd("~/Documents/GitHub/path_integral/results/convergence2Na10")
list.files()
df <- data.frame(tme=c(),
                 geg=c(),
                 fig=c())
for(tme in c(0.05, 0.1, 0.15)){
  for(geg in c(10,30,50)){
    df <- dplyr::bind_rows(df,
                           data.frame(tme,
                                      geg,
                                      path.expand(paste("2na10_", geg,"geg","_",tme,"time.svg",sep=""))))
  }
}
df <- df %>% rename(fig = path.expand.paste..2na10_...geg...geg...._...tme...time.svg...,
                    m_max = geg,
                    Time = tme)

the_colors <- c(rgb(0.396811, 0.31014, 0.2041), 
                rgb(0.64274, 0.330577, 0.1540755),
                rgb(0.726732, 0.538136, 0.31593),
                rgb(0.817882, 0.7260905, 0.426991),
                rgb(0.831964, 0.810543, 0.372854),
                rgb(0.6419975, 0.7183185, 0.366907),
                rgb(0.35082, 0.595178, 0.853742))

pointdf <- data.table(x = 1,
                      y = 1,
                      cols = c(0:5, "NDSolve"),
                      lntyp = c(rep(1,6), 2)) 

ggplot(df) +
  geom_image(aes(x=1,
                 y=1,
                 image=fig),
             size=Inf) +
  facet_grid(rows=vars(m_max),
             cols=vars(Time),
             scales = "free",
             labeller = label_both) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position="bottom") +
  xlab("Ending Frequency") +
  ylab("Density") + 
  scale_color_manual(limits = names(the_colors), 
                     values = the_colors) +
  geom_line(data = pointdf, aes(x = x, 
                                y = y,
                                color = cols,
                                linetype = cols),
            alpha = 0) + 
  guides(color = guide_legend(override.aes = list(alpha = 1),
                              nrow = 1)) +
  labs(color = "k_max",
       linetype = "k_max") + 
  scale_color_manual(values = the_colors, 
                     name = "k_max") +
  scale_linetype_manual(values=c(rep("solid",6), "dashed"),
                        name="k_max") + 
  theme(axis.title = element_text(size=18))
#############
#### RFS ####
#############

rm(list=ls())

setwd("/media/nathan/T7/path_integral/improvedPDetectedAlphaNE")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("selCoef", "Ne", "time", "start", 
                   "thresh", "pintAUC", "pintdetected", 
                   "numAUC", "numdetected", "statDist")
master$popalpha = 2 * master$Ne * master$selCoef
rm(tmp)

master <- master %>% filter(selCoef == 0.01)
master$reps = 1000 / master$Ne

df <- data.table()
for(i in 1:nrow(master)){
  tmp <- master[i,]
  for(j in 0:tmp$reps){
    foo <- data.table(Ne = tmp$Ne,
                      reps = paste(tmp$reps, 
                                   " Replicates of Size ", 
                                   tmp$Ne),
                      reps_number = tmp$reps,
                      bin = j,
                      prob = choose(tmp$reps, j) * 
                        tmp$pintdetected^j * 
                        (1 - tmp$pintdetected)^(tmp$reps - j))
    df <- dplyr::bind_rows(df, foo)
  }
}
rm(foo, master, tmp, file, i, j)

dfdetected <- df %>% filter(bin > 0) %>% 
  group_by(reps) %>% 
  mutate(pintdetected = sum(prob)) %>% 
  mutate(prob = prob / pintdetected) %>% ungroup() 
dfnotdetected <- df %>% filter(bin == 0)
# dfnotdetected$reps_number <- c(1,3,2.5,1.5)

ggplot(dfdetected, aes(y=prob, x=as.factor(bin))) + 
  geom_bar(stat="identity",
           fill = viridis(4)[2],
           alpha = 0.5) + 
  facet_wrap(~ reps, scales="free_x") + 
  theme_bw() + 
  xlab("Number of Replicates Detected") + 
  ylab("Probability Given Detected at least Once") + 
  geom_text(data = dfnotdetected, aes(x = reps_number, 
                                      y = 1.1, 
                                      label = paste("P(detected > 0 times) : ",round(1- prob,3))),
            hjust = 0.5, vjust = 1) + 
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size=18),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12)) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))

########################
#### RFS all 1 size ####
########################

rm(list=ls())

setwd("~/Documents/GitHub/path_integral/results/improvedPDetectionAlphaVG")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("alpha", "VG", "start", 
                   "time", "thresh", "pintAUC", 
                   "pintDetected", "numAUC", "numDetected",
                   "statDist")
master$popalpha = 1000 * master$alpha
rm(tmp, file)

master$statDist <- statDist <- master$statDist %>% 
  gsub('[{}]', '', .) %>% 
  gsub("\\*\\^", "e", .) %>% 
  as.numeric()

master <- master %>% filter(VG  == 1e-3)

df <- data.table()
for(i in 1:nrow(master)){
  tmp <- master[i,]
  for(reps in c(5, 10, 15, 20)){
    for(j in 0:reps){
      foo <- data.table(reps = paste(reps, " Replicates"),
                        reps_number = reps,
                        bin = j,
                        popalpha = tmp$popalpha,
                        prob = choose(reps, j) * 
                          tmp$pintDetected^j * 
                          (1 - tmp$pintDetected)^(reps - j))
      df <- dplyr::bind_rows(df, foo)
    }
  }
}
rm(foo, tmp, file, i, j)

dfdetected <- df %>% filter(bin > 0) %>% 
  group_by(reps, popalpha) %>% 
  mutate(Pdetected = sum(prob)) %>% 
  mutate(prob = prob / Pdetected,
         # # bin = as.factor(bin),
         # Ne = as.factor(Ne),
         reps = as.factor(reps),
         popalpha = as.factor(popalpha)) %>% ungroup() 
dfdetected$reps <- factor(dfdetected$reps, 
                          levels = levels(dfdetected$reps)[c(4,1,2,3)])

dfnotdetected <- df %>% filter(bin == 0) %>% 
  mutate(reps = as.factor(reps),
         popalpha = as.factor(popalpha),
         pdetected = paste("P(detected > 0 times) : ",formatC(round(1-prob,3),3,format="f"))) %>%
  mutate(pdetected = as.factor(pdetected))
dfnotdetected$reps_number <- rep(c(3.5,
                                   6.75,
                                   9.75,
                                   13),
                                 6)
dfnotdetected$vjust = rep((1.5 * 1:6), each = 4)


breaksfun <- function(x){
  1:max(x)
}

ggplot(dfdetected, aes(y = prob,
                       x = bin,
                       color = popalpha)) + 
  geom_line(alpha = 0.6,
            size = 1.5) + 
  geom_point(alpha=1,
             size = 2.5) +
  facet_wrap(~ reps, scales="free") + 
  theme_bw() + 
  scale_x_continuous(breaks = breaksfun) + 
  scale_color_viridis(discrete = T) + 
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size=18),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12)) + 
  xlab("Number of Replicates Detected") + 
  ylab("Probability Given Detected at least Once") + 
  geom_text(data = dfnotdetected, aes(x = reps_number,
                                      y = 1,
                                      label = pdetected,
                                      vjust = vjust)) + 
  guides(color=guide_legend(title = "Selection\nCoefficient",
                            override.aes = list(alpha=1))) 

#############
#### RFS ####
#############

rm(list=ls())

setwd("/media/nathan/T7/path_integral/improvedPDetectedAlphaNE")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("selCoef", "Ne", "time", "start", 
                   "thresh", "pintAUC", "pintdetected", 
                   "numAUC", "numdetected", "statDist")
master$popalpha = 2 * master$Ne * master$selCoef
rm(tmp)

# master <- master %>% filter(selCoef == 0.01)
master$reps = 5000 / master$Ne

df <- data.table()
for(i in 1:nrow(master)){
  tmp <- master[i,]
  for(j in 0:tmp$reps){
    foo <- data.table(Ne = tmp$Ne,
                      reps = paste(tmp$reps),
                      selCoef = tmp$selCoef,
                      reps_number = tmp$reps,
                      bin = j,
                      prob = choose(tmp$reps, j) * 
                        tmp$pintdetected^j * 
                        (1 - tmp$pintdetected)^(tmp$reps - j))
    df <- dplyr::bind_rows(df, foo)
  }
}
rm(foo, master, tmp, file, i, j)

df <- df %>% filter(Ne != 250)

dfdetected <- df %>% filter(bin > 0) %>% 
  group_by(reps, selCoef) %>% 
  mutate(pintdetected = sum(prob)) %>% 
  mutate(prob = prob / pintdetected,
         selCoef = as.factor(selCoef),
         reps = as.factor(reps)) %>% ungroup() 
dfdetected$reps <- factor(dfdetected$reps,
                          levels = levels(dfdetected$reps)[c(3,1,2,4)])

dfnotdetected <- df %>% filter(bin == 0) %>% 
  mutate(reps = as.factor(reps),
         selCoef = as.factor(selCoef))

my_labeller = as_labeller(
  c("5" = "5 ~ `Replicates,` ~ N[e] == 1000",
    "10" = "10 ~ `Replicates,` ~ N[e] == 500",
    "20" = "20 ~ `Replicates,` ~ N[e] == 250",
    "25" = "25 ~ `Replicates,` ~ N[e] == 200",
    "50"  = "50 ~ `Replicates,` ~ N[e] == 100"),
  default = label_parsed
)


for(i in 1:4){
  tmpdetected <- dfdetected %>% filter(reps == levels(dfdetected$reps)[i])
  tmpnotdetected <- dfnotdetected %>% filter(reps == levels(dfdetected$reps)[i])
  
  p <- ggplot(tmpdetected, aes(y=prob, x=bin, color = selCoef,
                               group = selCoef)) + 
    geom_line(alpha = 0.6,
              size = 1.5) + 
    geom_point(alpha=1,
               size = 2.5) +
    facet_wrap(~ factor(reps,
                        levels = unique(dfdetected$reps)[c(2,4,3,1)]), 
               scales="free_x",
               labeller = my_labeller) + 
    theme_bw() + 
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
    scale_color_manual(values = turbo(10)[c(2,4,6,7,9)],
                       labels = paste(levels(dfnotdetected$selCoef),
                                      "              ",
                                      formatC(round(1-tmpnotdetected$prob[c(5,1:4)],3),
                                              3,format="f"), 
                                      "         ",
                                      sep ="")) +  
    xlab("Number of Replicates Detected") + 
    ylab("Probability Given Detected at least Once") + 
    guides(color=guide_legend(title = "   \u03b1    P(detected > 0 times)",
                              override.aes = list(alpha=1,
                                                  linewidth = 0.75,
                                                  size = 1.25),
                              label.position = "left")) +
    theme(axis.title = element_blank(),
          title = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.justification = c("right", "top"),
          legend.position = c(.98,.98),
          legend.box.background = element_rect(colour = "black"),
          legend.spacing.y = unit(0, 'cm'),
          legend.key.size = unit(0.8, "line"),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) 
  
  if(i == 1) p1 <- p
  if(i == 2) p2 <- p + scale_x_continuous(breaks = c(1, seq(2,10,by = 2))) 
  if(i == 3) p3 <- p + scale_x_continuous(breaks = c(1, seq(2,10,by = 2)), limits = c(1,10)) 
  if(i == 4) p4 <- p + scale_x_continuous(breaks = c(1, seq(2,10,by = 2)), limits = c(1,10)) 
}

p5 <- ggplot(data.frame(l = "Number of Replicates Detected", 
                        x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), 
            # angle = 90,
            size = 5 / 14 * 12) + 
  theme_void() +
  coord_cartesian(clip = "off")

p6 <- ggplot(data.frame(l = "Probability Given Detected at least Once", 
                        x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), 
            angle = 90,
            size = 5 / 14 * 12) + 
  theme_void() +
  coord_cartesian(clip = "off")

layout <- "
ABC
ADE
#FF
"

p6 + p1 + p2 + p3 + p4 + p5 + plot_layout(design = layout,
                                          widths = c(1,25,25),
                                          heights = c(25,25,1))


#############
#### RFS ####
#############

rm(list=ls())

setwd("/media/nathan/T7/path_integral/improvedPDetectedAlphaNE")
list.files()
master <- data.frame()
for(file in list.files()){
  tmp <- fread(file)
  master <- dplyr::bind_rows(master, tmp)
}
names(master) <- c("selCoef", "Ne", "time", "start", 
                   "thresh", "pintAUC", "pintdetected", 
                   "numAUC", "numdetected", "statDist")
master$popalpha = 2 * master$Ne * master$selCoef
rm(tmp)

master <- master %>% filter(selCoef == 1e-2)

df <- data.table()
for(i in 1:nrow(master)){
  tmp <- master[i,]
  for(reps in c(5, 10, 15, 20)){
    for(j in 0:reps){
      foo <- data.table(reps = paste(reps, " Replicates"),
                        reps_number = reps,
                        bin = j,
                        selCoef = tmp$selCoef,
                        Ne = tmp$Ne,
                        prob = choose(reps, j) * 
                          tmp$pintdetected^j * 
                          (1 - tmp$pintdetected)^(reps - j))
      df <- dplyr::bind_rows(df, foo)
    }
  }
}
rm(foo, tmp, file, i, j)

df <- df %>% filter(Ne != 250)

dfdetected <- df %>% filter(bin > 0) %>% 
  group_by(reps, Ne) %>% 
  mutate(pintdetected = sum(prob)) %>% 
  mutate(prob = prob / pintdetected,
         # # bin = as.factor(bin),
         Ne = as.factor(Ne),
         reps = as.factor(reps)) %>% ungroup() 
dfdetected$reps <- factor(dfdetected$reps, 
                          levels = levels(dfdetected$reps)[c(4,1,2,3)])

dfnotdetected <- df %>% filter(bin == 0) %>% 
  mutate(reps = as.factor(reps),
         Ne = as.factor(Ne))

for(i in 1:4){
  tmpdetected <- dfdetected %>% filter(reps == levels(dfdetected$reps)[i])
  tmpnotdetected <- dfnotdetected %>% filter(reps == levels(dfdetected$reps)[i])
  
  p <- ggplot(tmpdetected, aes(y=prob, x=bin, color = Ne,
                               group = Ne)) + 
    geom_line(alpha = 0.6,
              size = 1.5) + 
    geom_point(alpha=1,
               size = 2.5) +
    facet_wrap(~ reps, 
               scales="free_x") + 
    theme_bw() + 
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
    scale_color_manual(values = turbo(10)[c(2,4,7,9)],
                       labels = paste(levels(dfnotdetected$Ne),
                                      "              ",
                                      formatC(round(1-tmpnotdetected$prob[c(1,3,4,2)],3),
                                              3,format="f"), 
                                      "             ",
                                      sep ="")) +  
    xlab("Number of Replicates Detected") + 
    ylab("Probability Given Detected at least Once") + 
    guides(color=guide_legend(title = bquote("   "~N[e]~"    P(detected > 0 times)"),
                              override.aes = list(alpha=1,
                                                  linewidth = 0.75,
                                                  size = 1.25),
                              label.position = "left")) +
    theme(axis.title = element_blank(),
          title = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.justification = c("right", "top"),
          legend.position = c(.98,.98),
          legend.box.background = element_rect(colour = "black"),
          legend.spacing.y = unit(0, 'cm'),
          legend.key.size = unit(0.8, "line"),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) 
  
  if(i == 1) p1 <- p
  if(i == 2) p2 <- p + scale_x_continuous(breaks = c(1, seq(2,10,by = 2)))
  if(i == 3) p3 <- p + scale_x_continuous(breaks = seq(1,15,by = 2))
  if(i == 4) p4 <- p + scale_x_continuous(breaks = c(1, seq(2,20,by = 2)))
}

p5 <- ggplot(data.frame(l = "Number of Replicates Detected", 
                        x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), 
            # angle = 90,
            size = 5 / 14 * 12) + 
  theme_void() +
  coord_cartesian(clip = "off")

p6 <- ggplot(data.frame(l = "Probability Given Detected at least Once", 
                        x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), 
            angle = 90,
            size = 5 / 14 * 12) + 
  theme_void() +
  coord_cartesian(clip = "off")

layout <- "
ABC
ADE
#FF
"

p6 + p1 + p2 + p3 + p4 + p5 + plot_layout(design = layout,
                                          widths = c(1,25,25),
                                          heights = c(25,25,1))

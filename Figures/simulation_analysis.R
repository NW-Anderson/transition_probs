library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)
library(reshape2)

################
#### Linked ####
################
setwd("/media/nathan/T7/path_integral/simulations/out")

allele_freqs <- fread("linked_positive_eff_09_11.csv.gz")  %>% 
  mutate(clr = "1",
         bigU = case_when(group_id %in% c("U=0.0025_a=0.005", "U=0.0025_a=0.01") ~ "U = 0.0025",
                          group_id %in% c("U=0.025_a=0.005", "U=0.025_a=0.01") ~ "U = 0.025"),
         selCoef = case_when(group_id %in% c("U=0.025_a=0.01", "U=0.0025_a=0.01") ~ "\u03b1 = 0.01",
                             group_id %in% c("U=0.025_a=0.005", "U=0.0025_a=0.005") ~ "\u03b1 = 0.005"))

setwd("/media/nathan/T7/path_integral/simulations/out/pintComparison")

pints <- fread("linked_sim_pint_densities.csv") %>% 
  mutate(clr = "2",
         bigU = case_when(Scenario %in% c(1,2) ~ "U = 0.0025",
                          Scenario %in% c(3,4) ~ "U = 0.025"),
         selCoef = case_when(Scenario %in% c(1,3) ~ "\u03b1 = 0.005",
                             Scenario %in% c(2,4) ~ "\u03b1 = 0.01"))

setwd("/media/nathan/T7/path_integral/genicSimComparison/k=5")

genic <- data.frame()
for(file in list.files()){
  if(file != "k=5"){
    params <- strsplit(file, split = "_")[[1]] 
    end <- params[2] %>% 
      gsub("end", "", .) %>% 
      as.numeric()
    selCoef <- params[3] %>% 
      gsub("selCoef.csv", "", .) %>% 
      as.numeric()
    
    tmp <- fread(file) %>%
      unlist()
    genic <- dplyr::bind_rows(genic, data.frame(dens = tmp,
                                                end = end, 
                                                selCoef = selCoef))
  }
}

genic <- genic %>% 
  mutate(selCoef = case_when(selCoef == 5 ~ "\u03b1 = 0.005",
                             selCoef == 10 ~ "\u03b1 = 0.01"),
         clr = "3")

# pints <- pints %>% filter(end_freq <= 0.5, end_freq > 0)
# allele_freqs <- allele_freqs %>% filter(end <= 0.5, end > 0)
# genic <- genic %>% filter(end <= 0.5, end > 0)

pints <- pints %>% filter(end_freq > 0)
allele_freqs <- allele_freqs %>% filter(end > 0)
genic <- genic %>% filter(end > 0)

pint_means <- pints %>% group_by(bigU, selCoef) %>% 
  mutate(total_dens = sum(Dens)) %>% 
  mutate(norm_dens = Dens / total_dens) %>%
  summarize(expectation = sum(norm_dens * end_freq)) %>%
  ungroup()

group_means <- allele_freqs %>% group_by(bigU, selCoef) %>%
  summarize(mn = mean(end)) %>%
  ungroup()

genic_means <- genic %>% group_by(selCoef) %>%
  mutate(total_dens = sum(dens)) %>%
  mutate(norm_dens = dens / total_dens) %>%
  summarize(expectation = sum(norm_dens * end)) %>%
  ungroup()

ggplot(data = allele_freqs, aes(x = end, color = clr)) + 
  geom_vline(data = group_means, aes(xintercept = mn), 
             linetype = "dashed",
             color = turbo(10)[2]) +
  geom_vline(data = pint_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(10)[8]) + 
  geom_vline(data = genic_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(10)[4]) + 
  geom_density(data = allele_freqs, show.legend = F, size = 1.25,
               alpha = 1.75) + 
  geom_line(data = pints, aes(x = end_freq, y = Dens, color = clr),
            size = 1.25,
            alpha = 0.75) +
  geom_line(data = genic, aes(x = end, y = dens, color = clr),
            size = 1.25,
            alpha = 0.75) + 
  xlim(0,1) + 
  labs(
    title = "Linked: Ending Frequency for Alleles Starting Between 0.09 and 0.11") +
  scale_x_continuous("Ending Frequency",
                     breaks = seq(0,1,0.2)) + 
  scale_y_continuous("Density") + 
  facet_grid(rows = vars(bigU),
             cols = vars(selCoef)) + 
  theme_bw() +
  guides(color=guide_legend(title = element_blank(),
                            override.aes = list(alpha=1,
                                                size = 1.25,
                                                linewidth = 0.75))) +
  scale_color_manual(values = c(turbo(10)[2],
                                turbo(10)[8],  
                                turbo(10)[4]),
                     name = "",
                     labels = c("Linked\n Simulation",
                                "Polygenic\n Selection",
                                "Genic")) + 
  theme(axis.title = element_text(size=12),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification = c("right", "top"),
        legend.position = c(.99,.99),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 12))  

########
# main #
########

allele_freqs <- allele_freqs %>% 
  filter(group_id == "U=0.025_a=0.01",
         end <= 0.5)

pints <- pints %>% filter(Scenario == 4,
                          end_freq <= 0.5)

genic <- genic %>% filter(selCoef == "α = 0.01",
                          end <= 0.5)
genic_means <-  genic_means %>% filter(selCoef == "α = 0.01")
group_means <- group_means %>% filter(bigU == "U = 0.025",
                                      selCoef == "α = 0.01")
pint_means <- pint_means %>% filter(bigU == "U = 0.025",
                                    selCoef == "α = 0.01")

tmp <- as.data.frame(density(allele_freqs$end)[1:2]) %>% 
  mutate(clr = "1") %>%
  filter(x <= 0.5,
         x >= 0.01)
  

p1 <- ggplot(data = allele_freqs, aes(x = end, color = clr)) + 
  geom_vline(data = pint_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(10)[8]) + 
  geom_vline(data = genic_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(10)[4]) + 
  geom_vline(data = group_means, aes(xintercept = mn), 
             linetype = "dashed",
             color = turbo(10)[2]) +
  # geom_vline(xintercept = 0.1,
  #            linetype = "dashed",
  #            color = "gray80") +
  # geom_density(data = allele_freqs, show.legend = F, size = 1.25,
  #              alpha = 0.75) + 
  geom_line(data = tmp, aes(x = x, y = y, color = clr),
            size = 1.25,
            alpha = 0.75) +
  geom_line(data = pints, aes(x = end_freq, y = Dens, color = clr),
            size = 1.25,
            alpha = 0.75) +
  geom_line(data = genic, aes(x = end, y = dens, color = clr),
            size = 1.25,
            alpha = 0.75) +
  labs(
    title = "Linked simulation") +
  guides(color=guide_legend(title = element_blank(),
                            override.aes = list(alpha=1,
                                                size = 1.25,
                                                linewidth = 0.75)
                            # ,
                            # by_row =T,
                            # keyheight = 1.75
                            )) +
  scale_x_continuous("Ending frequency",
                     breaks = seq(0,1,0.1),
                     labels = c(0,expression(italic(x) == 0.1), seq(0.2,1,by = 0.1)),
                     limits = c(0,0.5), 
                     expand = c(0, 0)) + 
  scale_y_continuous("Density", 
                     limits = c(0,max(pints$Dens) * 1.05),
                     expand = (c(0,0))) + 
  theme_bw() +
  scale_color_manual(values = c(turbo(10)[2], 
                                turbo(10)[8],  
                                turbo(10)[4]),
                     name = "",
                     labels = c("Linked simulation",
                                "Polygenic selection",
                                "Genic")) + 
  theme(text = element_text(family = "LM Roman 10"),
        axis.title = element_text(size=10),
        title = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification = c("right", "top"),
        legend.position = c(.99,.99),
        # legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 10))  

p1

##################
#### Unlinked ####
##################

setwd("/media/nathan/T7/path_integral/simulations/out")

allele_freqs <- fread("unlinked_positive_eff_09_11.csv.gz")  %>% 
  mutate(clr = "1",
         bigU = case_when(group_id %in% c("mu_0.0025.a_0.005", "mu_0.0025.a_0.01") ~ "U = 0.0025",
                          group_id %in% c("mu_0.025.a_0.005", "mu_0.025.a_0.01") ~ "U = 0.025"),
         selCoef = case_when(group_id %in% c("mu_0.025.a_0.01", "mu_0.0025.a_0.01") ~ "\u03b1 = 0.01",
                             group_id %in% c("mu_0.025.a_0.005", "mu_0.0025.a_0.005") ~ "\u03b1 = 0.005"))

setwd("/media/nathan/T7/path_integral/simulations/out/pintComparison")

pints <- fread("unlinked_sim_densities.csv.gz") %>% 
  mutate(clr = "2",
         bigU = case_when(scenario %in% c(1,2) ~ "U = 0.0025",
                          scenario %in% c(3,4) ~ "U = 0.025"),
         selCoef = case_when(scenario %in% c(1,3) ~ "\u03b1 = 0.005",
                             scenario %in% c(2,4) ~ "\u03b1 = 0.01"))

setwd("/media/nathan/T7/path_integral/genicSimComparison/k=5")

genic <- data.frame()
for(file in list.files()){
  if(file != "k=5"){
    params <- strsplit(file, split = "_")[[1]] 
    end <- params[2] %>% 
      gsub("end", "", .) %>% 
      as.numeric()
    selCoef <- params[3] %>% 
      gsub("selCoef.csv", "", .) %>% 
      as.numeric()
    
    tmp <- fread(file) %>%
      unlist()
    genic <- dplyr::bind_rows(genic, data.frame(dens = tmp,
                                                end = end, 
                                                selCoef = selCoef))
  }
}

genic <- genic %>% 
  mutate(selCoef = case_when(selCoef == 5 ~ "\u03b1 = 0.005",
                             selCoef == 10 ~ "\u03b1 = 0.01"),
         clr = "3")

# pints <- pints %>% filter(end_freq <= 0.5, end_freq > 0)
# allele_freqs <- allele_freqs %>% filter(end_freqs <= 0.5, end_freqs > 0)
# genic <- genic %>% filter(end <= 0.5, end > 0)

pints <- pints %>% filter(end_freq > 0)
allele_freqs <- allele_freqs %>% filter(end_freqs > 0)
genic <- genic %>% filter(end > 0)

pint_means <- pints %>% group_by(bigU, selCoef) %>% 
  mutate(total_dens = sum(density)) %>% 
  mutate(norm_dens = density / total_dens) %>%
  summarize(expectation = sum(norm_dens * end_freq)) %>% 
  ungroup()

group_means <- allele_freqs %>% group_by(bigU, selCoef) %>%
  summarize(mn = mean(end_freqs))

genic_means <- genic %>% group_by(selCoef) %>%
  mutate(total_dens = sum(dens)) %>%
  mutate(norm_dens = dens / total_dens) %>%
  summarize(expectation = sum(norm_dens * end)) %>%
  ungroup()

ggplot(data = allele_freqs, aes(x = end_freqs, color = clr)) + 
  geom_vline(data = group_means, aes(xintercept = mn), 
             linetype = "dashed",
             color = turbo(10)[2]) + 
  geom_vline(data = pint_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(10)[8]) + 
  geom_vline(data = genic_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(10)[4]) +
  geom_density(data = allele_freqs, 
               show.legend = F, 
               size = 1.25,
               alpha = 1.75) + 
  geom_line(data = pints, aes(x = end_freq, y = density, color = clr),
            size = 1.25,
            alpha = 0.75) +
  geom_line(data = genic, aes(x = end, y = dens, color = clr),
            size = 1.25,
            alpha = 0.75) + 
  xlim(0,1) + 
  labs(
    title = "Unlinked: Ending Frequency for Alleles Starting Between 0.09 and 0.11") +
  scale_x_continuous("Ending Frequency",
                     breaks = seq(0,1,0.2)) + 
  scale_y_continuous("Density") + 
  facet_grid(cols = vars(selCoef),
             rows = vars(bigU)) + 
  theme_bw() +
  guides(color=guide_legend(title = element_blank(),
                            override.aes = list(alpha=1,
                                                size = 1.25,
                                                linewidth = 0.75))) +
  scale_color_manual(values = c(turbo(10)[2],  
                                turbo(10)[8], 
                                turbo(10)[4]),
                     name = "",
                     labels = c("Unlinked\n Simulation",
                                "Polygenic\n Selection",
                                "Genic")) + 
  theme(axis.title = element_text(size=12),
        title = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification = c("right", "top"),
        legend.position = c(.99,.99),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 12)) 

########
# main #
########


allele_freqs <- allele_freqs %>% 
  filter(group_id == "mu_0.025.a_0.01",
         end_freqs <= 0.5)

pints <- pints %>% filter(scenario == 4,
                          end_freq <= 0.5,
                          end_freq >= 0.01)

genic <- genic %>% filter(selCoef == "α = 0.01",
                          end <= 0.5)
genic_means <-  genic_means %>% filter(selCoef == "α = 0.01")
group_means <- group_means %>% filter(bigU == "U = 0.025",
                                      selCoef == "α = 0.01")
pint_means <- pint_means %>% filter(bigU == "U = 0.025",
                                    selCoef == "α = 0.01")

tmp <- as.data.frame(density(allele_freqs$end)[1:2]) %>% 
  mutate(clr = "1") %>%
  filter(x <= 0.5,
         x >= 0.01)

p2 <- ggplot(data = allele_freqs, aes(x = end_freqs, color = clr)) + 
  geom_vline(data = pint_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(10)[8]) + 
  geom_vline(data = genic_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(10)[4]) +
  geom_vline(data = group_means, aes(xintercept = mn), 
             linetype = "dashed",
             color = turbo(10)[2]) + 
  # geom_vline(xintercept = 0.1,
  #            linetype = "dashed",
  #            color = "gray80") +
  # geom_density(data = allele_freqs, show.legend = F, size = 1.25,
  #              alpha = 0.75) + 
  geom_line(data = tmp, aes(x = x, y = y, color = clr),
            size = 1.25,
            alpha = 0.75) +
  geom_line(data = pints, aes(x = end_freq, y = density, color = clr),
            size = 1.25,
            alpha = 0.75) +
  geom_line(data = genic, aes(x = end, y = dens, color = clr),
            size = 1.25,
            alpha = 0.75) + 
  labs(
    title = "Unlinked simulation") +
  scale_x_continuous("Ending frequency",
                     breaks = seq(0,1,0.1),
                     labels = c(0,expression(italic(x) == 0.1), seq(0.2,1,by = 0.1)),
                     limits = c(0,0.5), 
                     expand = c(0, 0)) + 
  scale_y_continuous("Density", 
                     limits = c(0,max(pints$density) * 1.05),
                     expand = (c(0,0))) +  
  theme_bw() +
  guides(color=guide_legend(title = element_blank(),
                            override.aes = list(alpha=1,
                                                size = 1.25,
                                                linewidth = 0.75)
                            # ,
                            # by_row =T,
                            # keyheight = 1.75
                            )) +
  scale_color_manual(values = c(turbo(10)[2],  
                                turbo(10)[8], 
                                turbo(10)[4]),
                     name = "",
                     labels = c("Unlinked simulation",
                                "Polygenic selection",
                                "Genic")) + 
  theme(text = element_text(family = "LM Roman 10"),
        axis.title = element_text(size=10),
        title = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification = c("right", "top"),
        legend.position = c(.99,.99),
        # legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 10)) 

p2

# p1 + p2 + plot_layout(guides = "collect") + 
#   plot_annotation(tag_levels = 'A')  & 
#   theme(plot.tag = element_text(size = 24),
#         legend.position = "bottom")

#############
## Neutral ##
#############

setwd("/media/nathan/T7/path_integral/simulations/out")

allele_freqs <- fread("linked_neutral_09_11.csv.gz")  %>% 
  mutate(clr = "1",
         bigU = case_when(par %in% c("U=0.0025_a=0.005", "U=0.0025_a=0.01") ~ "U = 0.0025",
                          par %in% c("U=0.025_a=0.005", "U=0.025_a=0.01") ~ "U = 0.025"),
         selCoef = case_when(par %in% c("U=0.025_a=0.01", "U=0.0025_a=0.01") ~ "\u03b1 = 0.01",
                             par %in% c("U=0.025_a=0.005", "U=0.0025_a=0.005") ~ "\u03b1 = 0.005"))

setwd("/media/nathan/T7/path_integral/simulations/out/KimuraComparison")

pints <- fread("linked_sim_kimura_densities.csv") %>% 
  mutate(clr = "2")

setwd("/media/nathan/T7/path_integral/trueNeutralSims")

true_neut <- fread("trueNeutralMaster.csv.gz") %>%
  mutate(clr = "3")


# pints <- pints %>% filter(end_freq <= 0.5, end_freq > 0)
# allele_freqs <- allele_freqs %>% filter(end_freq <= 0.5, end_freq > 0)
# true_neut <- true_neut %>% filter(freq <= 0.5, freq > 0)

pints <- pints %>% filter(end_freq > 0)
allele_freqs <- allele_freqs %>% filter(end_freq > 0)
true_neut <- true_neut %>% filter(freq > 0)

pint_means <- pints %>% group_by(Scenario) %>% 
  mutate(total_dens = sum(Dens)) %>% 
  mutate(norm_dens = Dens / total_dens) %>%
  summarize(expectation = sum(norm_dens * end_freq))

group_means <- allele_freqs %>% group_by(bigU, selCoef) %>%
  summarize(mn = mean(end_freq))

true_neut_means <- true_neut %>% summarize(mn = mean(freq))

ggplot(data = allele_freqs, aes(x = end_freq, color = clr)) + 
  geom_vline(data = pint_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(10)[4]) +
  geom_vline(data = group_means, aes(xintercept = mn),
             linetype = "dashed",
             color = turbo(10)[2]) +
  geom_vline(data = true_neut_means, aes(xintercept = mn),
             linetype = "dashed",
             color = turbo(10)[8]) +
  geom_density(data = true_neut, aes(x = freq, color =clr),
               size = 1.25,
               alpha = 0.75,
               show.legend = F) + 
  geom_line(data = pints, aes(x = end_freq, y = Dens, color = clr),
            size = 1.25,
            alpha = 0.75) +
  geom_density(show.legend = F, size = 1.25,
               alpha = 0.75) +
  xlim(0,1) + 
  labs(
    title = "Neutral: Ending Frequency for Alleles Starting Between 0.09 and 0.11") +
  scale_x_continuous("Ending Frequency",
                     breaks = seq(0,1,0.2)) + 
  scale_y_continuous("Density") + 
  facet_grid(cols = vars(selCoef),
             rows = vars(bigU)) + 
  theme_bw() +
  guides(color=guide_legend(title = element_blank(),
                            override.aes = list(alpha=1,
                                                size = 1.25,
                                                linewidth = 0.75))) +
  scale_color_manual(values = c(turbo(10)[8],  
                                turbo(10)[4], 
                                turbo(10)[2]),
                     name = "",
                     labels = c("Neutral + Linked\n Simulation",
                                "Neutral\n Simulation",
                                "Kimura's\n Solution")) + 
  theme(axis.title = element_text(size=12),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification = c("right", "top"),
        legend.position = c(.99,.99),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 12)) 

########
# main #
########

allele_freqs <- allele_freqs %>% filter(par == "U=0.025_a=0.01",
                                        end_freq <= 0.5)
group_means <- group_means %>% filter(bigU == "U = 0.025",
                                      selCoef == "α = 0.01")
pints <- pints %>% filter(end_freq <= 0.5)
true_neut <- true_neut %>% filter(freq <= 0.5)

tmp <- as.data.frame(density(allele_freqs$end)[1:2]) %>% 
  mutate(clr = "1") %>%
  filter(x <= 0.5,
         x >= 0.01)

trueNeuttmp <- as.data.frame(density(true_neut$freq)[1:2]) %>% 
  mutate(clr = "3") %>%
  filter(x <= 0.5,
         x >= 0.01)

p3 <- ggplot(data = allele_freqs, aes(x = end_freq, color = clr)) + 
  geom_vline(data = group_means, aes(xintercept = mn),
             linetype = "dashed",
             color = turbo(10)[8]) +
  geom_vline(data = true_neut_means, aes(xintercept = mn),
             linetype = "dashed",
             color = turbo(10)[2]) +
  geom_vline(data = pint_means, aes(xintercept = expectation),
             linetype = "dashed",
             color = turbo(10)[4]) +
  # geom_vline(xintercept = 0.1,
  #            linetype = "dashed",
  #            color = "gray80") +
  geom_line(data = tmp, aes(x = x, y = y, color = clr),
            size = 1.25,
            alpha = 0.75) +
  geom_line(data = trueNeuttmp, aes(x = x, y = y, color = clr),
            size = 1.25,
            alpha  = 0.75) + 
  # geom_density(data = allele_freqs, show.legend = F, size = 1.25,
  #              alpha = 0.75) +
  # geom_density(data = true_neut, aes(x = freq, color =clr),
  #              size = 1.25,
  #              alpha = 0.75,
  #              show.legend = F) +
  geom_line(data = pints, aes(x = end_freq, y = Dens, color = clr),
            size = 1.25,
            alpha = 0.75) +
  scale_x_continuous("Ending frequency",
                     breaks = seq(0,1,0.1),
                     labels = c(0,expression(italic(x) == 0.1), 
                                seq(0.2,1,by = 0.1)),
                     limits = c(0,0.5), 
                     expand = c(0, 0)) + 
  scale_y_continuous("Density", 
                     limits = c(0,max(pints$Dens) * 1.05),
                     expand = (c(0,0))) +
  labs(
    title = "Neutral simulation") +
  theme_bw() +
  guides(color=guide_legend(title = element_blank(),
                            override.aes = list(alpha=1,
                                                size = 1.25,
                                                linewidth = 0.75)
                            # ,
                            # by_row =T,
                            # keyheight = 1.75
                            )) +
  scale_color_manual(values = c(turbo(10)[8],
                                turbo(10)[4],
                                turbo(10)[2]),
                     name = "",
                     labels = c("Neutral + linked\n simulation",
                                "Kimura's solution",
                                "Neutral simulation")) +
  theme(text = element_text(family = "LM Roman 10"),
        axis.title = element_text(size=10),
        title = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification = c("right", "top"),
        legend.position = c(.99,.99),
        # legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 10)) 

p3

################
#### LD Fig ####
################

setwd("/home/nathan/Documents/GitHub/path_integral/simulations")

params  <- fread("params.txt")

setwd("/media/nathan/T7/path_integral/sigmaD-linlog")

master <- data.frame()
bincount_data <- data.frame()
count <- 0
for(file in list.files()){
  count <- count + 1
  print(paste(count, "XXX", file))
  
  cur_seed <- strsplit(file, split = "_")[[1]][2]
  
  par <- params %>% filter(V1 == cur_seed)
  raw_df <- fread(file)
  
  df <- raw_df %>% 
    select(-c("pos_bincount", 
              "neg_bincount",
              "posneg_bincount")) %>% 
    melt(id.vars = c("bin_start"),
         variable.name = "class") %>% 
    mutate(bigU = par$V2,
           selCoef = par$V3)
  master <- dplyr::bind_rows(master, df)
  
  df <- raw_df %>% 
    select(-c("pos", 
              "neg", 
              "posneg")) %>%
    melt(id.vars = "bin_start",
         variable.name = "class") %>%
    mutate(bigU = par$V2,
           selCoef = par$V3)
  bincount_data <- dplyr::bind_rows(bincount_data, df)
}

bins <- c(unique(master$bin_start), 1e8)
tmp <- unique(master$bin_start) + diff(bins) / 2
df <- data.frame(bin_start = unique(master$bin_start),
                 bin_mid = tmp)
# master$bin_start <- master$bin_start / 1e6 + 0.5
d <- master %>% 
  merge(., df, by = "bin_start") %>%
  mutate(bin_mid  = bin_mid * 1e-8  * 4e4) %>%
  group_by(bin_start,
           class, 
           bin_mid, 
           bigU, 
           selCoef) %>%
  summarise(mn = mean(value, na.rm = T)) %>%
  mutate(bigU = paste("U = ", bigU, sep = ""),
         selCoef = paste("\u03b1 = ", selCoef, sep = ""))

# master <- master %>% filter(bin_start <= 5e7)

ggplot(d, aes(x = bin_mid, y = mn, color = class)) + 
  geom_line(size = 1.25,
            alpha = 0.75) + 
  facet_grid(cols = vars(selCoef),
             rows = vars(bigU),
             scales = "free") + 
  theme_bw() + 
  xlab(bquote("Recombination Distance (4" * N * r * ")")) + 
  ylab(bquote(bar(sigma)[D]^1)) +
  labs(title = "Linkage Disequilibrium") + 
  scale_x_log10() +
  theme(legend.position = "bottom") +  
  theme(axis.title = element_text(size=12),
        title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification = c("right", "top"),
        legend.position = c(.99,.99),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 12)) +
  scale_color_manual(values = c(turbo(10)[8],  
                                turbo(10)[2], 
                                turbo(10)[4]),
                     name = "",
                     labels = c("+/+",
                                "-/-",
                                "+/-")) +
  guides(color=guide_legend(title = element_blank(),
                            override.aes = list(alpha=1,
                                                size = 1.25,
                                                linewidth = 0.75))) 

b <- bincount_data %>% 
  merge(., df, by = "bin_start") %>%
  # mutate(bin_mid  = bin_mid * 1e-6) %>%
  mutate(bin_mid  = bin_mid * 1e-8  * 4e4) %>%
  group_by(bin_start,
           class, 
           bin_mid, 
           bigU, 
           selCoef) %>%
  summarise(mn = mean(value, na.rm = T)) %>%
  mutate(bigU = paste("U = ", bigU, sep = ""),
         selCoef = paste("\u03b1 = ", selCoef, sep = ""))

ggplot(b, aes(x = bin_mid, y = mn, color = class)) + 
  geom_line(size = 1.25,
            alpha = 0.75) + 
  facet_grid(cols = vars(selCoef),
             rows = vars(bigU),
             scales = "free") + 
  theme_bw() + 
  xlab(bquote("Recombination Distance (4" * N * r * ")")) + 
  ylab("Mean Number of observations") +
  labs(title = "Linkage Disequilibrium") + 
  scale_x_log10(breaks = c(1e-4,
                           1e-2,
                           1,
                           1e2)) +
  scale_y_log10() + 
  theme(legend.position = "bottom") + 
  theme(axis.title = element_text(size=18),
        title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_color_manual(values = c(turbo(10)[9],  
                                turbo(10)[3], 
                                turbo(10)[6]),
                     name = "",
                     labels = c("+/+",
                                "-/-",
                                "+/-"))
########
# Main #
########

d <- d %>% filter(bigU == "U = 0.025",
                            selCoef == "α = 0.01",
                  bin_mid <= 500)

setwd("/media/nathan/T7/path_integral/sigmaD_moments/take2")
rhos <- fread("rhoval.csv") %>% unlist()
LD_diff <- fread("MomentsLDdiff.csv") %>% unlist()
LD_same <- fread("MomentsLDsame.csv") %>% unlist()

momentsDf<- dplyr::bind_rows(data.frame(rho = rhos,
                                        ld = LD_diff,
                                        class = "diff"),
                             data.frame(rho = rhos,
                                        ld = LD_same,
                                        class = "same")) %>%
  filter(rho >= 1.315,
         rho < 50) 

# %>%
#   mutate(ld = ld/2) # factor of two problem

p4 <- ggplot(d, aes(x = bin_mid, y = mn, color = class, linetype = class)) +
  geom_line(size = 1.25,
            alpha = 0.75) + 
  geom_line(data = momentsDf, aes(x = rho, y = ld, color = class, linetype = class),
            size = 3, alpha = 0.5) +
  theme_bw() + 
  xlab(bquote("Recombination distance (4" * italic(N) * italic(r) * ")")) + 
  ylab(bquote(bar(sigma)[italic(D)]^1)) +
  labs(title = "Linkage disequilibrium",
       linetype = "class",
       color = "class") +
  scale_x_log10() + 
  scale_color_manual(values = c(turbo(10)[8],
                                turbo(10)[2],
                                turbo(10)[4],
                                turbo(10)[7],
                                turbo(10)[3]),
                     name = "class",
                     labels = c("+/-",
                                "-/-",
                                "+/+",
                                "Pred. diff.",
                                "Pred. same"),
                     breaks = c("posneg",
                                "neg",
                                "pos",
                                "diff",
                                "same")
                     ) +
  scale_linetype_manual(name = "class",
                        values = c(rep("solid",3),
                                   rep("11",2)),
                        labels = c("+/-",
                                   "-/-",
                                   "+/+",
                                   "Pred. diff.",
                                   "Pred. same"),
                        breaks = c("posneg",
                                   "neg",
                                   "pos",
                                   "diff",
                                   "same")) +
  theme(text = element_text(family = "LM Roman 10"),
        axis.title = element_text(size=10),
        title = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification = c("right", "top"),
        legend.position = c(.99,.99),
        # legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0.8, "line"),
        panel.grid= element_blank(),
        strip.text = element_text(size = 10)) + 
  guides(color = guide_legend(title = element_blank(),
                              ncol = 2,
                              override.aes = list(alpha = 1,
                                                  linewidth = 2)),
         linetype = guide_legend(title = element_blank(),
                                 ncol = 2,
                                 override.aes = list(alpha = 1,
                                                     size = 1)))

p4

(p1 + p2) / (p4 + p3) + 
  plot_annotation(tag_levels = 'A')  &
  theme(plot.tag = element_text(size = 12),
        plot.margin = margin(0, 5.5, 0, 5.5, "pt")
        )

# setwd("/media/nathan/T7/path_integral")
# 
# eff_master <- fread("effect_sizes.csv")
# 
# setwd("/media/nathan/T7/path_integral/LD-matricies-linked/results")
# 
# breaks <- 0:100 * 10^6
# 
# completed_seeds <- c()
# master <- data.frame()
# count <- 0
# for(file in list.files()){
#   
#   #find out current seed
#   cur_seed <- strsplit(file, split = "_")[[1]][2]
#   
#   if(!cur_seed %in% completed_seeds){
#     count <- count + 1
#     completed_seeds <- c(completed_seeds, cur_seed)
#     
#     print(paste(count, "XXX", cur_seed))
#     
#     # read in site info 
#     sites <- fread(paste("result", cur_seed, "sites.csv.gz", sep = "_"))
#     sites$site_id <- sites$site_id + 1
#     
#     # filter effect sizes for current seed, and for the sites in site info and merge
#     eff <- eff_master %>% filter(seed == cur_seed,
#                                  site %in% sites$site_id) %>% 
#       merge(., sites, by.x = "site", by.y = "site_id")
#     
#     rm(sites)
#     
#     # read in ld matrix, reshape into data frame and filter for lower triangle
#     d <- fread(paste("result", 
#                      cur_seed, 
#                      "D.csv.gz", 
#                      sep = "_")) %>% 
#       as.matrix() %>% 
#       reshape2::melt(value.name = "d") %>%
#       mutate(Var2 = as.integer(Var2)) %>%
#       filter(Var1 > Var2)
#     
#     # read in pi2 matrix, reshape into data frame and filter for lower triangle
#     pi2 <- fread(paste("result", 
#                        cur_seed, 
#                        "pi2.csv.gz", 
#                        sep = "_")) %>% 
#       as.matrix() %>% 
#       reshape2::melt(value.name = "pi2") %>%
#       mutate(Var2 = as.integer(Var2)) %>%
#       filter(Var1 > Var2)
#     
#     # merge (col are in same order)
#     d$pi2 <- pi2$pi2
#     
#     # clean up memory
#     rm(pi2)
#     gc()
#     
#     # filter fixed sites
#     d <- d %>% filter(pi2 > 0)
#     
#     # merge with effect size info for each site
#     d <- d %>% mutate(site_i = eff$site[Var1],
#                       eff_i = eff$effect_size[Var1],
#                       pos_i = eff$site_position[Var1],
#                       site_j = eff$site[Var2],
#                       eff_j = eff$effect_size[Var2],
#                       pos_j = eff$site_position[Var2]) %>%
#       select(-c("Var1", "Var2")) %>% 
#       mutate(dist = pos_i - pos_j)
#     
#     d$bin <- 10^6 * (cut(d$dist, breaks = breaks, labels = F) - 1) + 5e5
#     
#     d <- d %>% mutate(group = case_when(eff_i > 0 & eff_j > 0 ~ "+/+",
#                                         eff_i > 0 & eff_j < 0 ~ "+/-",
#                                         eff_i < 0 & eff_j > 0 ~ "+/-",
#                                         eff_i < 0 & eff_j < 0 ~ "-/-")) %>%
#       group_by(bin, group) %>% 
#       summarize(sigma = sum(d)/ sum(pi2)) %>% 
#       mutate(seed = cur_seed)
#     
#     master <- dplyr::bind_rows(master, d)
#     
#     fwrite(d, file = paste("../results_",
#                            cur_seed,
#                            ".csv",
#                            sep = ""))
#   }else{
#     print("womp")
#   }
# }
# 
# # ggplot(d_sum, aes(x = bin, y = sigma, color = group)) + 
# #   geom_line()
# 
# 
# 
# fwrite(master, file = "../master.csv")
# table(master$bin)







####################
#### Linked RFS ####
####################
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

thresh <- 0.1 + unique(master$thresh)
rm(master)

setwd("/media/nathan/T7/path_integral/simulations/out")

allele_freqs <- fread("linked_positive_eff_09_11.csv.gz")  %>% 
  mutate(clr = "1",
         bigU = case_when(group_id %in% c("U=0.0025_a=0.005", "U=0.0025_a=0.01") ~ "U = 0.0025",
                          group_id %in% c("U=0.025_a=0.005", "U=0.025_a=0.01") ~ "U = 0.025"),
         selCoef = case_when(group_id %in% c("U=0.025_a=0.01", "U=0.0025_a=0.01") ~ "\u03b1 = 0.01",
                             group_id %in% c("U=0.025_a=0.005", "U=0.0025_a=0.005") ~ "\u03b1 = 0.005"))

cases <- unique(allele_freqs$group_id)
master <- data.frame()
for(i in cases){
  df <- allele_freqs %>% filter(group_id == i)
  for(j in 1:10000){
    tmp <- df %>% filter(seed == sample(seed,1)) %>% 
      filter(site == sample(site,1),
             deme %in% sample(1:100, 10)) %>%
      mutate(detected = end >= thresh)
    
    print(paste(i,j, sum(tmp$detected), sep = " : "))
    master <- dplyr::bind_rows(master,
                               data.frame(case = i,
                                          rep = j,
                                          num_detected = sum(tmp$detected)))
    
                         
  }
}

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
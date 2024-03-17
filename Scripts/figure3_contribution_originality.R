


# This is code is adapted from M. Paola Barajas Barbosa et al (2022) Assembly of functional diversity 
# in an oceanic island flora

# Figure 3 

library(tidyr) 
library(hypervolume) 
library(BAT) 
library(FactoMineR) 
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(scales)
library(rstatix)

# data
traits_orig <- read.csv("Desktop/Projects_Git/lianas_GEB/data_prepared/traits_orig.csv")

data_lianas_imp <- read.csv("Desktop/Projects_Git/lianas_GEB/data_prepared/Trait_imputed_phylo8.csv" ) 

data_lianas <- left_join (data_lianas, data_lianas_imp, by = "IDmaster")

data_lianas <- data_lianas %>%   # replace NAs with the imputed data
  mutate(sla.x = coalesce(sla.x, sla.y)) %>%
  mutate(la.x = coalesce(la.x, la.y)) %>%
  mutate(seed.x = coalesce(seed.x, seed.y)) %>%
  mutate(leafnitro.x = coalesce(leafnitro.x, leafnitro.y)) %>%
  mutate(wd.x = coalesce(wd.x, wd.y)) 

anyNA(data_lianas)
data_lianas$source <- 'lianas_data'

data_lianas <- data.frame(
  Species1      = data_lianas$accepted_taxon_name,
  
  `Leaf area`     = scale(log10(data_lianas$la.x)),
  SLA           = scale(log10(data_lianas$sla.x)),
  `Leaf N`        = scale(log10(data_lianas$leafnitro.x)),
  `Seed mass`     = scale(log10(data_lianas$seed.x)),
  `Stem density`  = scale(log10(data_lianas$wd.x)),
  
  Climbing_mechanism = data_lianas$habit,
  source       = data_lianas$source) 

# Null model ----
# PCA analysis
pca <- FactoMineR::PCA(data_lianas[, c("Leaf area", "SLA", "Leaf N", "Seed mass", 
                                       "Stem density")],
                       scale.unit = FALSE, graph = FALSE)
pca$eig
sp_coord <- as.data.frame(pca$ind$coord[, 1:3])
sp_coord <- cbind(sp_coord, data_lianas$Species1)
sp_coord$Species1 <- sp_coord$`data_lianas$Species1`
sp_coord[,4] <- NULL

data_lianas <- dplyr::left_join(data_lianas, sp_coord, by = "Species1")
rownames(data_lianas) <- data_lianas$Species1

# We calculate the three following indices: kernel.alpha, kernel.dispersion and kernel.evennes
# Indices based on hypervolume
all_bw <- estimate_bandwidth(data_lianas[, c("Dim.1", "Dim.2", "Dim.3")], # Estimate bandwidth for lianas
                             method = "cross-validation")

# ####
min(table(data_lianas$habit))
min_rar_bio <- 30
nbperm <- 999 # Number of permutations
####

act  <- data_lianas[which(data_lianas$Climbing_mechanism == "active"),  c("Dim.1", "Dim.2", "Dim.3")]
pas <- data_lianas[which(data_lianas$Climbing_mechanism == "passive"), c("Dim.1", "Dim.2", "Dim.3")]
all_sp <- data_lianas[, c("Dim.1", "Dim.2", "Dim.3")]

rar_bio_comp <- c()
for (i in 1:nbperm){
  # All species - null community
  all_sp_i <- all_sp[sample(1:nrow(all_sp), size = min_rar_bio),
                     c("Dim.1", "Dim.2", "Dim.3")]
  
  all_sp_i <-  hypervolume_gaussian(
    all_sp_i, kde.bandwidth = all_bw,
    quantile.requested = 0.95, quantile.requested.type = "probability",
    verbose = FALSE)
  
  all_sp_i <- data.frame(perm = i,
                         status = "All_species",
                         rich = kernel.alpha(all_sp_i),
                         eve = kernel.evenness(all_sp_i),
                         div = kernel.dispersion(all_sp_i))
  
  # Active species
  ACT_i <- act[sample(1:nrow(act), size = min_rar_bio),
               c("Dim.1", "Dim.2", "Dim.3")]
  
  ACT_i <-  hypervolume_gaussian(
    ACT_i, kde.bandwidth = all_bw,
    quantile.requested = 0.95, quantile.requested.type = "probability",
    verbose = FALSE)
  
  ACT_i <- data.frame(perm = i,
                      status = "active",
                      rich = kernel.alpha(ACT_i),
                      eve = kernel.evenness(ACT_i),
                      div = kernel.dispersion(ACT_i))
  
  # Passive species
  PAS_i <- pas[sample(1:nrow(pas), size = min_rar_bio),
               c("Dim.1", "Dim.2", "Dim.3")]
  
  PAS_i <-  hypervolume_gaussian(
    PAS_i, kde.bandwidth = all_bw,
    quantile.requested = 0.95, quantile.requested.type = "probability",
    verbose = FALSE)
  
  PAS_i <- data.frame(perm = i,
                      status = "passive",
                      rich = kernel.alpha(PAS_i),
                      eve = kernel.evenness(PAS_i),
                      div = kernel.dispersion(PAS_i))
  
  # Bind results
  rar_bio_comp <- rbind(rar_bio_comp, all_sp_i, ACT_i, PAS_i)
  
  cat(paste0(round(100*i/nbperm, 0), " %; "))
}


# # # 
save(rar_bio_comp, file = "Desktop/Projects_Git/lianas_GEB/data_prepared/Figure_3a.RData")

load("Desktop/Projects_Git/lianas_GEB/data_prepared/Figure_3a.RData")

rar_bio_tidy <- gather(rar_bio_comp, metric, val, c("rich", "eve", "div"))

FR <- filter(rar_bio_tidy, metric == "rich") 
FE <- filter(rar_bio_tidy, metric == "eve") 
FD <- filter(rar_bio_tidy, metric == "div") 

sum_fr <- FR %>%       
  group_by(status) %>%
  summarise( n=n(), mean_pao = mean(val), sd_pao = sd(val) ) %>%
  mutate( se = sd_pao/sqrt(n))  %>%
  mutate( ic = se * qt((1-0.05)/2 + .5, n-1))  %>%
  mutate( mean_pao, metric= "Functional Richness") 

sum_fe <- FE %>%
  group_by(status) %>%
  summarise( n=n(), mean_pao = mean(val), sd_pao = sd(val) ) %>%
  mutate( se = sd_pao/sqrt(n))  %>%
  mutate( ic = se * qt((1-0.05)/2 + .5, n-1)) %>%
  mutate( mean_pao, metric= "Functional Evenness") 

sum_fd <- FD %>%
  group_by(status) %>%
  summarise( n=n(), mean_pao = mean(val), sd_pao = sd(val) ) %>%
  mutate( se = sd_pao/sqrt(n))  %>%
  mutate( ic = se * qt((1-0.05)/2 + .5, n-1)) %>%
  mutate( mean_pao, metric= "Functional Dispersion") 
sum_metrics <-  rbind(sum_fr, sum_fe, sum_fd)

# Figure 3 Null model / Standardized values 
p_my_theme2 <-  theme(axis.text.x=element_text(colour=c("black"),face="bold",size=10),
                      axis.text.y=element_text(colour=c("black"),face="bold",size=10),
                      panel.background =element_rect(fill="transparent",colour="black"),
                      panel.grid.minor=element_blank(),
                      panel.border=element_rect(fill=NA,colour="grey"))

coco = c("dodgerblue3", "turquoise4", "gold2")

Figure_3a <- ggplot(sum_metrics) +
  geom_linerange(aes(factor(status, level = c ("active", "passive", "All_species")), 
                     ymin= mean_pao-ic, ymax= mean_pao+ic, color = status),  size = 0.3,  
                 show.legend = FALSE) +
  geom_point    (aes(x = status,  
                     y=mean_pao, color = status),  size = 3,  
                 show.legend = FALSE) +
  scale_colour_manual(values = coco) +
  facet_wrap(~ metric, scales = "free") +
  labs(title = "(a)", x = "", y = "") +
  p_my_theme2 +
  # Set the custom labels for the x-axis
  scale_x_discrete(labels = c("Active", "Passive", "Null"))

Figure_3a

# Figure 3 b

data_i        <- traits_orig[, c ("IDmaster","accepted_taxon_name","sla","la","seed",
                                  "leafnitro","wd","habit")]


data_i <- left_join (data_i, data_lianas_imp, by = "IDmaster")
data_i <- data_i %>%   # replace the few NAs with the imputed data
  mutate(sla.x = coalesce(sla.x, sla.y)) %>%
  mutate(la.x = coalesce(la.x, la.y)) %>%
  mutate(seed.x = coalesce(seed.x, seed.y)) %>%
  mutate(leafnitro.x = coalesce(leafnitro.x, leafnitro.y)) %>%
  mutate(wd.x = coalesce(wd.x, wd.y)) 

anyNA(data_i)

data_i_trans <- data.frame(
  IDmaster      = data_i$IDmaster,
  Species1      = data_i$accepted_taxon_name,
  Leaf_area     = scale(log10(data_i$la.x)),
  SLA           = scale(log10(data_i$sla.x)),
  Leaf_N        = scale(log10(data_i$leafnitro.x)),
  Seed_mass     = scale(log10(data_i$seed.x)),
  Stem_density  = scale(log10(data_i$wd.x)),
  clbimbing_status = data_i$habit)

rownames(data_i_trans) <- data_i_trans$Species1

# 
pca_2 <- FactoMineR::PCA( data_i_trans[, c("Leaf_area","SLA","Leaf_N","Seed_mass","Stem_density")]
                          , scale.unit = FALSE, graph = FALSE)

sp_coord_2 <- as.data.frame(pca_2$ind$coord[, 1:3])

sp_coord_2$Species1 <- rownames(sp_coord_2)
sp_tra <- dplyr::left_join(data_i_trans, sp_coord_2, by = "Species1")
rownames(sp_tra) <- sp_tra$Species1

PCA_2         <- prcomp(sp_tra[, c("Leaf_area","SLA","Leaf_N","Seed_mass","Stem_density")])
PCAloadings_2 <- data.frame(Variables = rownames(PCA_2$rotation), PCA_2$rotation)

# Hypervolume
all_bw <- estimate_bandwidth(sp_tra[, c("Dim.1", "Dim.2", "Dim.3")],
                             method = "cross-validation")

set.seed(8)

hv_lianas <- hypervolume_gaussian(
  sp_tra[, c("Dim.1", "Dim.2", "Dim.3")],
  name = "lianas_volume",
  kde.bandwidth = all_bw,
  quantile.requested = 0.95,
  quantile.requested.type = "probability", verbose = FALSE)

sp_contrib <- BAT::kernel.contribution(hv_lianas)
sp_orig    <- BAT::kernel.originality(hv_lianas)

sp_tra$contrib <- sp_contrib 
sp_tra$orig    <- sp_orig


save(sp_tra, file = "Desktop/Projects_Git/lianas_GEB/data_prepared/sp_tra_semi.Rdata")

load("Desktop/Projects_Git/lianas_GEB/data_prepared/sp_tra_semi.RData") 

# ggplot2 theme
tema <- theme(panel.background =element_rect(fill="transparent",colour="black"),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(fill=NA,colour="grey"),
              axis.text.y =element_text(colour=c("black"), size=10),
              axis.text.x =element_text(colour=c("black"), size=10))

# Is the functional contribution and originality different across groups? 

coco_2 = c("dodgerblue3", "gold2") 

con <- ggplot(sp_tra) +
  geom_boxplot(aes(x= factor(clbimbing_status, level = c ("active","passive")),
                   y=contrib, color = clbimbing_status), show.legend = FALSE) +
  scale_colour_manual(values = coco_2) +
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  ylab("Contribution") + xlab("") + tema + labs(title = "(b)") +
  scale_x_discrete(labels = stringr::str_wrap(c("Active (n = 545)", "Passive (n = 157)", width = 9))) 

ori <- ggplot(sp_tra) +
  geom_boxplot(aes(x= factor(clbimbing_status, level = c ("active","passive")),
                   y=orig, color = clbimbing_status), show.legend = FALSE) +
  scale_colour_manual(values = coco_2) +
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  ylab("Originality") + xlab("") + tema + labs(title = "(c)") +
  annotate(
    "text", label = "*",
    x = 1.5, y = 3, size = 10, colour = "black"
  ) +
  scale_x_discrete(labels = stringr::str_wrap(c("Active (n = 545)", "Passive (n = 157)", width = 9)))


# t-test contribution and orginality between groups
t_contrib <- sp_tra %>%
  t_test(., contrib ~ clbimbing_status) # P = 0.6

t_orig <- sp_tra %>% 
  t_test(., orig ~ clbimbing_status) # p = 0.01

####

fig3 <- ggarrange(Figure_3a,                                                 # First row with scatter plot
                  ggarrange(con, ori, ncol = 2), # Second row with box and dot plots
                  nrow = 2)                                     


ggsave(fig3, filename = "fig3_2023.png", bg = "white",
       path = "Desktop/Projects_Git/lianas_GEB/figures/",
       height = 9, width = 15, scale = .75)


# END
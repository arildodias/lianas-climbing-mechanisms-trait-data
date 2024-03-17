

# This is code is adapted from M. Paola Barajas Barbosa et al (2022) Assembly of functional diversity 
# in an oceanic island flora

# # Figure 1

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggExtra)
library(cowplot)


p_my_theme1 <-  theme( axis.title.x=element_text(colour="black",face="bold",size=10),
                       axis.title.y=element_text(colour="black",face="bold",size=10),
                       axis.text.x=element_text(colour=c("black"),face="bold",size=10),
                       axis.text.y=element_text(colour=c("black"),face="bold",size=10),
                       legend.position  = c(0.11,0.15), legend.direction="vertical",    # position
                       legend.key       = element_rect(fill="transparent"),
                       legend.key.size  = unit(.5,"line"),
                       legend.title     = element_blank(), 
                       legend.text      = element_text(size= 10, color="black"),
                       panel.background = element_rect(fill="transparent",colour="black"),
                       panel.grid.minor = element_blank(),
                       panel.border     = element_rect(fill=NA,colour="grey"))

p_my_theme2 <-  theme( axis.text.x=element_text(colour=c("white"),face="bold",size=10),
                       axis.text.y=element_text(colour=c("white"),face="bold",size=10),
                       legend.position  = c(0.11,0.15), legend.direction="vertical",    
                       legend.key       = element_rect(fill="transparent"),
                       legend.key.size  = unit(.3,"line"),
                       legend.title     = element_blank(), 
                       legend.text      = element_text(size= 9, color="black"),
                       panel.background = element_rect(fill="transparent",colour="black"),
                       panel.grid.minor = element_blank(),
                       panel.border     = element_rect(fill=NA,colour="grey"))


# data
data_lianas <- read.csv("Desktop/Projects_Git/lianas_GEB/data_prepared/traits_orig.csv")

data_lianas_imp <- read.csv("Desktop/Projects_Git/lianas_GEB/data_prepared/Trait_imputed_phylo8.csv" ) 

data_lianas <- left_join (data_lianas, data_lianas_imp, by = "IDmaster")

data_lianas <- data_lianas %>%   # replace NAs with the imputed data
  mutate(sla.x = coalesce(sla.x, sla.y)) %>%
  mutate(la.x = coalesce(la.x, la.y)) %>%
  mutate(seed.x = coalesce(seed.x, seed.y)) %>%
  mutate(leafnitro.x = coalesce(leafnitro.x, leafnitro.y)) %>%
  mutate(wd.x = coalesce(wd.x, wd.y)) 

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


data_lianas %>% group_by(Climbing_mechanism) %>% summarise(n = n())


# PCA 
PCA             <- prcomp(data_lianas[, c("Leaf area","SLA","Leaf N","Seed mass","Stem density")] )
summary(PCA)

PCAvalues       <- data.frame(Species = data_lianas$Species1, 
                              Climbing_mechanism = data_lianas$Climbing_mechanism, 
                              PCA$x)
PCAvalues$PC1   <- PCAvalues$PC1*-1 ; PCAvalues$PC2 <- PCAvalues$PC2*-1 # for visualization purposes

PCAloadings     <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)  
PCAloadings$PC1 <- PCAloadings$PC1*-1 ; PCAloadings$PC2 <- PCAloadings$PC2*-1


# 1. Plots Figure 1  
LIANAS <-
  PCAvalues %>% ggplot(aes(PC1,  PC2), size = 1)+ # Plot Lianas PCA 
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)),  
                  colour = "gray", bins = 34) +
  
  scale_fill_distiller(palette = "BuGn", direction = 1) +
  geom_jitter(alpha=0.6,  size = .7  , colour = "turquoise4") +     #  Display the points 
  geom_text(data = PCAloadings, aes(x = PC1*4.7, y = PC2*4.7, label = Variables), size = 4) +
  geom_segment(data = PCAloadings, size = 0.2,    # Plots the loadings, i.e., traits 
               aes(x = 0, xend = PC1*4.2, y = 0, yend = PC2*4.2),
               arrow = arrow(length = unit(0.1, "cm")),colour = "black")   + 
  xlab("PC1 (50%)") + ylab("PC2 (20%)")   +
  xlim(-5.75 , 5.75) + ylim(-5.75, 5.75) +
  labs(title = "(a)  Liana trait space") +
  #xlim(-3 , 3) + ylim(-3, 3) +
  p_my_theme1

LIANAS


# Plot each group separately
ACTIVE <- dplyr::filter(PCAvalues, Climbing_mechanism == "active")
ACTIVE <- ACTIVE  %>% ggplot(aes(PC1,  PC2))+ 
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "gray", bins = 10) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_jitter(alpha=0.5,  size = 0.7  , colour = "dodgerblue3") +    
  xlim(-5.75 , 5.75) + ylim(-5.75, 5.75) +
  xlab("PC1 (50%)") + ylab("PC2 (20%)") + 
  labs(title = "(b)  Active climbing species") +
  p_my_theme1

ACTIVE

PASSIVE <- dplyr::filter(PCAvalues, Climbing_mechanism == "passive")
PASSIVE <- PASSIVE %>% ggplot(aes(PC1,  PC2))+ 
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "gray", bins = 10) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_jitter(alpha=0.5,  size = 0.7  , colour = "gold2") +  
  xlim(-5.75 , 5.75) + ylim(-5.75, 5.75) +
  xlab("PC1 (50%)") + ylab("PC2 (20%)") +
  labs(title =  "(c)  Passive climbing species") +
  p_my_theme1

PASSIVE

Lianas <- ggExtra:: ggMarginal(LIANAS, type = "density", fill="transparent", size = 15)
Active <- ggExtra:: ggMarginal(ACTIVE, type = "density", fill="transparent", size = 15) 
Passive <- ggExtra:: ggMarginal(PASSIVE, type = "density", fill="transparent", size = 15)


fig1 <- plot_grid(Lianas, Active, Passive,
                  ncol = 3, nrow = 1, rel_widths = c(1,1,1,1))

ggsave(fig1, filename = "fig1_2023.png", bg = "white",
       path = "Desktop/Projects_Git/lianas_GEB/figures/",
       height = 9, width = 15, scale = .75)

# END
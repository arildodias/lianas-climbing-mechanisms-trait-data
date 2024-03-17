


# load libraries

library(terra)
library(sf)
library(tidyverse)
library(rWCVP)
library(readxl)
library(phyloregion)
library(ape)
library(jtools)
library(lme4)
library(MuMIn)
library(r2glmm)
library(DHARMa) 
library(cowplot)

# colors: passive gold2, active dodgerblue3

# ==============================================================================
# Analyses
# ==============================================================================


data_lianas <- read.csv("Desktop/Projects_Git/lianas_GEB/data_prepared/data_lianas.csv", header = T)


data_epi <- read_xlsx("Desktop/Projects_Git/lianas_GEB/data_raw/vascular_epiphytes/Supplementary_GEB.xlsx", sheet = 5)

# phylogenetic tree
tree_lianas <- read.tree("Desktop/Projects_Git/lianas_GEB/data_prepared/phylotree_lianas.txt")

# data for the sparse matrix
comm_bot_countries <- st_read("Desktop/Projects_Git/lianas_GEB/data_prepared/liana_countries.shp") %>% 
  st_make_valid()


data_glm <- data_epi %>% 
  dplyr::select(geo_entity, biome, longitude, latitude, lat_absolute, kingdom_WWF,
                area, TF_Miocene, TF_current, TF_LGM, prec, PrecSeas, temp_mean,
                ElevRange, LGM_ice) %>% 
  left_join(data_lianas, by = c("geo_entity" = "LEVEL3_N")) %>% 
  dplyr::select(!c(longitude.y, latitude.y)) %>% 
  rename(longitude = longitude.x, latitude = latitude.x) %>% 
  rowwise() %>% 
  mutate(rich = sum(active, passive, na.rm = TRUE)) %>% 
  ungroup() 


##############
# Phylogenetic 

com_bot_countries <- comm_bot_countries %>% 
  st_centroid() %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(country = comm_bot_countries$LEVEL3_N) %>% # Replace with the actual column name
  rename(longitude = X, latitude = Y) %>% 
  left_join(comm_bot_countries, by = c("country" = "LEVEL3_N")) %>% 
  dplyr::select(country, longitude, latitude, specs_d, habit) 

# we perform phylobetadiversity separately for active and passive species
com_active <- com_bot_countries %>% 
  dplyr::filter(habit == "active") %>% 
  dplyr::select(!habit) %>% 
  rename(species = specs_d) %>% 
  mutate(species = str_replace(species, " ", "_")) # to match taxon name with phylogenetic tree names

com_passive <- com_bot_countries %>% 
  dplyr::filter(habit == "passive") %>% 
  dplyr::select(!habit) %>% 
  rename(species = specs_d) %>% 
  mutate(species = str_replace(species, " ", "_")) # to match taxon name with phylogenetic tree names


rWCVPdata::wgsrpd3 %>% 
  st_join(st_set_crs(shp_mx_active, 4326)) %>% 
  rename(beta_ses_active = cur_sim) %>% 
  ggplot() +
  geom_sf(aes(fill = beta_ses_active))


# phylobetadiversity

sparse_com_active <- points2comm(dat = com_active, res = 0.5, lon = "longitude", lat = "latitude",
                                 species = "species") 

sparse_com_passive <- points2comm(dat = com_passive, res = 0.5, lon = "longitude", lat = "latitude",
                                  species = "species") # We select sparse_com_passive$comm_dat as the matrix distance to subsequent analyses

# estimate betaphylodiversity effect sizes
phylobeta_active <- phylobeta_ses(sparse_com_active$comm_dat, tree_lianas, model = "tipshuffle")

phylobeta_passive <- phylobeta_ses(sparse_com_passive$comm_dat, tree_lianas, model = "tipshuffle")

# select phylobeta_obs_z representing the ses phylobetadiversity per region and run the PCoA
# to select the axes with high explanatory power

# Perform PCoA - the axes are shown
pcoa_active <- cmdscale(phylobeta_active$phylobeta_obs_z, eig=TRUE, add=TRUE)

pcoa_active_axes <- pcoa_active$points

colnames(pcoa_active_axes) <- c("pcoa1", "pcoa2")

percent_explained <- 100 * pcoa_active$eig / sum(pcoa_active$eig)

pcoa_active_axes <- pcoa_active_axes %>%
  as.data.frame() %>% 
  rownames_to_column(., var = "grids") 

pcoa_passive <- cmdscale(phylobeta_passive$phylobeta_obs_z, eig=TRUE, add=TRUE)

pcoa_passive_axes <- pcoa_passive$points

colnames(pcoa_passive_axes) <- c("pcoa1", "pcoa2")

percent_explained <- 100 * pcoa_passive$eig / sum(pcoa_passive$eig)

pcoa_passive_axes <- pcoa_passive_axes %>%
  as.data.frame() %>% 
  rownames_to_column(., var = "grids") 

# combine with the sparse_com_active$poly_shp and set crs = 4326
df_pcoa_active <- rWCVPdata::wgsrpd3 %>% 
  st_join(st_set_crs(st_as_sf(sparse_com_active$poly_shp), 4326)) %>% 
  st_drop_geometry() %>% 
  left_join(as.data.frame(pcoa_active_axes), by = "grids") %>% 
  group_by(LEVEL3_NAM) %>% 
  summarise(across(c(abundance, richness, pcoa1, pcoa2), ~ mean(., na.rm=T))) %>% 
  rename(pcoa1_active = pcoa1, pcoa2_active = pcoa2) %>% 
  select(-c(abundance, richness))


df_pcoa_passive <- rWCVPdata::wgsrpd3 %>% 
  st_join(st_set_crs(st_as_sf(sparse_com_passive$poly_shp), 4326)) %>% 
  st_drop_geometry() %>% 
  left_join(as.data.frame(pcoa_passive_axes), by = "grids") %>% 
  group_by(LEVEL3_NAM) %>% 
  summarise(across(c(abundance, richness, pcoa1, pcoa2), ~ mean(., na.rm=T))) %>% 
  rename(pcoa1_passive = pcoa1, pcoa2_passive = pcoa2) %>% 
  select(-c(abundance, richness))


######
# GLMM

data_active <- data_glm %>% 
  dplyr::select(geo_entity, kingdom_WWF, mean_temp, precipatation, temp_seasonality, 
                temp_anomaly, prec_anomaly, ivs, active, rich) %>% 
  rename(precipitation = precipatation) %>% 
  drop_na() %>% 
  left_join(df_pcoa_active, by = c("geo_entity" = "LEVEL3_NAM"))


data_passive <- data_glm %>% 
  dplyr::select(geo_entity, kingdom_WWF, mean_temp, precipatation, temp_seasonality, 
                temp_anomaly, prec_anomaly, ivs, passive, rich) %>%
  rename(precipitation = precipatation) %>% 
  drop_na() %>% 
  left_join(df_pcoa_passive, by = c("geo_entity" = "LEVEL3_NAM"))


glmm_active <- data_active %>% 
  mutate(precipitation_scaled = scale(log10(precipitation+1))[,1],
         mean_temp_scaled = scale(mean_temp)[,1],
         temp_seasonality_scaled = scale(log10(temp_seasonality+1))[,1],
         ivs_scaled =scale(log10(ivs))[,1],
         temp_anomaly_scaled = scale(log10(temp_anomaly+1))[,1],
         prec_anomaly_scaled =scale(log10(prec_anomaly+1))[,1],
         pcoa1_scaled = scale(pcoa1_active)[,1],
         pcoa2_scaled = scale(pcoa2_active)[,1]
  ) %>% 
  drop_na()


m1_active <- glmer(cbind(active,rich-active) ~ 
                     precipitation_scaled + mean_temp_scaled + 
                     temp_seasonality_scaled + ivs_scaled +
                     temp_anomaly_scaled + prec_anomaly_scaled +
                     pcoa1_scaled + pcoa2_scaled +
                     (1|kingdom_WWF), data=glmm_active, family = "binomial",
                   na.action = na.exclude, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),nAGQ = 0)


summ(m1_active)
r.squaredGLMM(m1_active)


glmm_passive <- data_passive %>% 
  mutate(precipitation_scaled = scale(log10(precipitation+1))[,1],
         mean_temp_scaled = scale(mean_temp)[,1],
         temp_seasonality_scaled = scale(log10(temp_seasonality+1))[,1],
         ivs_scaled =scale(log10(ivs))[,1],
         temp_anomaly_scaled = scale(log10(temp_anomaly+1))[,1],
         prec_anomaly_scaled =scale(log10(prec_anomaly+1))[,1],
         pcoa1_scaled = scale(pcoa1_passive)[,1],
         pcoa2_scaled = scale(pcoa2_passive)[,1]
  ) %>% 
  drop_na()

m1_passive <- glmer(cbind(passive,rich-passive) ~ 
                      precipitation_scaled + mean_temp_scaled + 
                      temp_seasonality_scaled + ivs_scaled +
                      temp_anomaly_scaled + prec_anomaly_scaled +
                      pcoa1_scaled + pcoa2_scaled +
                      (1|kingdom_WWF), data=glmm_passive, family = "binomial",
                    na.action = na.exclude, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),nAGQ = 0)


summ(m1_passive)
r.squaredGLMM(m1_passive)



fig5_a <- jtools::plot_summs(m1_active,
                             m1_passive,
                             scale = TRUE,
                             coefs = c("Preciptation" = "precipitation_scaled",
                                       "Temperature" = "mean_temp_scaled",
                                       "Seasonality" = "temp_seasonality_scaled",
                                       "IVS" = "ivs_scaled",
                                       "Temp anomaly" = "temp_anomaly_scaled",
                                       "Prec anomaly" = "prec_anomaly_scaled",
                                       "PCoA.axis_1" = "pcoa1_scaled",
                                       "PCoA.axis_2" = "pcoa2_scaled"),
                             colors = c("dodgerblue3", "gold2"),
                             model.names = c("Active", "Passive")) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 11),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.position = c(0.95, 0.40),
        legend.justification = c(1,1),
        legend.box.background = element_blank(),
        plot.margin = unit(c(0.25, 0, 0.25, 0), "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5)))+
  ggtitle("(a) Climbing mechanisms model")



# realtive importance of predictor variables (for fixed effects)
r_active <- r2beta(model = m1_active, method = 'sgv', data = glmm_active) # active
r_passive <- r2beta(model = m1_passive, method = 'sgv', data = glmm_passive) # passive

total_sum_active <- r_active %>% 
  filter(Effect != "Model") %>% 
  summarise(sum(Rsq)) %>% pull(.)

total_sum_passive <- r_passive %>% 
  filter(Effect != "Model") %>% 
  summarise(sum(Rsq)) %>% pull(.)

rel_var_active <- r_active %>% 
  select(Effect, Rsq) %>% 
  filter(Effect != "Model") %>% 
  mutate(var_type = case_when(Effect == "pcoa1_scaled" ~ "phylogeny",
                              Effect == "pcoa2_scaled" ~ "phylogeny",
                              Effect == "temp_seasonality_scaled" ~ "climate",
                              Effect == "mean_temp_scaled" ~ "climate",
                              Effect == "precipitation_scaled" ~ "climate",
                              Effect == "prec_anomaly_scaled" ~ "paleoclimate",
                              Effect == "temp_anomaly_scaled" ~ "paleoclimate",
                              Effect == "ivs_scaled" ~ "vegetation")) %>% 
  group_by(var_type) %>% 
  summarise(perc = (sum(Rsq))) %>% 
  mutate(perc_rescaled = perc / total_sum_active * 100) 

rel_var_passive <- r_passive %>% 
  select(Effect, Rsq) %>% 
  filter(Effect != "Model") %>% 
  mutate(var_type = case_when(Effect == "pcoa1_scaled" ~ "phylogeny",
                              Effect == "pcoa2_scaled" ~ "phylogeny",
                              Effect == "temp_seasonality_scaled" ~ "climate",
                              Effect == "mean_temp_scaled" ~ "climate",
                              Effect == "precipitation_scaled" ~ "climate",
                              Effect == "prec_anomaly_scaled" ~ "paleoclimate",
                              Effect == "temp_anomaly_scaled" ~ "paleoclimate",
                              Effect == "ivs_scaled" ~ "vegetation")) %>% 
  group_by(var_type) %>% 
  summarise(perc = (sum(Rsq))) %>% 
  mutate(perc_rescaled = perc / total_sum_passive * 100) 


var_explained <- rel_var_active %>% 
  select(-perc) %>% 
  rename(perc_active = perc_rescaled) %>% 
  left_join(rel_var_passive, by = "var_type") %>% 
  select(-perc) %>% 
  rename(perc_passive = perc_rescaled)


fig5_b <- rel_var_active %>% 
  left_join(rel_var_passive, by = "var_type") %>% 
  select(var_type, perc_rescaled.x, perc_rescaled.y) %>% 
  rename(perc_active = perc_rescaled.x, perc_passive = perc_rescaled.y) %>% 
  pivot_longer(cols = c("perc_active", "perc_passive")) %>% 
  ggplot(aes(fill=var_type, y=value, x=name)) + 
  geom_bar(position="stack", stat="identity", width = 0.6) +
  xlab("") + ylab("Variance explained (%)") +
  scale_fill_viridis_d() +
  scale_x_discrete(labels = c("Active", "Passive")) + 
  theme_minimal() +
  theme(legend.title = element_blank()) +
  ggtitle("(b) Variable importance")



fig5 <- plot_grid(fig5_a, fig5_b,
                  ncol = 1, nrow = 2, rel_widths = c(1,1,1,1))

ggsave(fig5, filename = "fig5_2023.png", bg = "white",
       path = "Desktop/Projects_Git/lianas_GEB/figures/",
       height = 10, width = 12, scale = .75)


# Moran coefficient

xy_active <- glmm_active %>% 
  select(geo_entity) %>% 
  left_join(com_active, by = c("geo_entity" = "country")) %>% 
  as_tibble() %>% 
  select(-species) %>% 
  distinct()


# test spatial autocorrelation
testSpatialAutocorrelation(simulateResiduals(m1_active), 
                           x = xy_active$longitude, 
                           y= xy_active$latitude)

# P < 0.001 but I = 0.10


xy_passive <- glmm_passive %>% 
  select(geo_entity) %>% 
  left_join(com_passive, by = c("geo_entity" = "country")) %>% 
  as_tibble() %>% 
  select(-species) %>% 
  distinct()

# test spatial autocorrelation
testSpatialAutocorrelation(simulateResiduals(m1_passive), 
                           x = xy_passive$longitude, 
                           y= xy_passive$latitude)
# P = 0.6

# END


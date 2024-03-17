

# Load libraries

library(tidyverse)
library(here)
library(janitor)
library(grid)
library(cowplot)
library(visdat)
library(patchwork)


theme1 <- theme( axis.title = element_text(size=9),
                 axis.text.x= element_text(size=9), axis.text.y = element_text(size=9),
                 legend.text= element_text(size= 9)) 


lianas <- read.csv(here::here("Desktop", "Projects_Git", "lianas_GEB", "data_prepared", "traits_orig.csv"),
                   stringsAsFactors = FALSE) %>% 
  rename(SLA = sla, LA = la, SM = seed, Nmass = leafnitro, WD = wd)


#  a. All lianas in our study
a <- vis_miss(lianas[,c("SLA","LA","SM","Nmass","WD")]) + 
  labs(title = "(a) Liana species in this study") +
  ylab("Species") +
  theme1


dflianas_1 <- read.csv(here::here("Documents", "R_projects", "GEB_resubmission", "data_wilsonetal2022", "full_met_analysis_data_8_Feb.csv"),
                       stringsAsFactors = FALSE) %>% 
  clean_names()

dflianas_2 <- read.csv(here::here("Documents", "R_projects", "GEB_resubmission", "data_wilsonetal2022", "filtered_TRY_analysis_16-03-21.csv"),
                       stringsAsFactors = FALSE) %>% 
  clean_names()


lianas1 <- dflianas_1 %>% 
  filter(growth_form == 'liana') %>% 
  select(acc_species_name, al_as_m2_cm2, n_g_kg, wd_g_cm3, aarea_umol_m2_s, sla_cm2_g) %>%
  rename(Species = acc_species_name, sla = sla_cm2_g, leaf_nitro = n_g_kg, wd = wd_g_cm3, amax = aarea_umol_m2_s, la = al_as_m2_cm2) %>%
  as_tibble()

lianas2 <- dflianas_2 %>% 
  filter(growth_form == 'liana') %>% 
  select(name, sla_mm2_mg, nmass_mg_g, ssd_g_cm3, aarea_mmol_co2_m2_s, la_cm_2) %>% 
  rename(Species = name, sla = sla_mm2_mg, leaf_nitro = nmass_mg_g, wd = ssd_g_cm3, amax = aarea_mmol_co2_m2_s, la = la_cm_2) %>%
  #mutate(amax = amax * la) %>% 
  as_tibble()

# wilson et al. (2022) complete trait information
lianas_wilson2022 <- lianas2 %>% 
  left_join(lianas1, by = c("Species" = "Species")) %>% 
  rowwise() %>% 
  mutate(
    la = mean(c(la.x, la.y), na.rm = T),
    sla = mean(c(sla.x, sla.y), na.rm = T),
    wd = mean(c(wd.x, wd.y), na.rm = T),
    leaf_nitro = mean(c(leaf_nitro.x, leaf_nitro.y), na.rm = T)
  ) %>% 
  select(Species, sla, la, leaf_nitro, wd) %>% 
  rename(species = Species, leafnitro = leaf_nitro) %>% 
  ungroup() %>% 
  rename(SLA = sla, LA = la, Nmass = leafnitro, WD = wd)

#  b. All lianas in Wilson et al. (2022)
b <- vis_miss(lianas_wilson2022[,c("SLA","LA","Nmass","WD")]) + 
  labs(title = "(b) Liana species in Wilson et al. (2022)") +
  ylab("Species") +
  theme1

fig_S1 <- a | b 

ggsave(fig_S1, filename = "figS1_trait_coverage.png", bg = "white",
       path = "Desktop/Projects_Git/lianas_GEB/figures/",
       height = 9, width = 15, scale = .75)

 
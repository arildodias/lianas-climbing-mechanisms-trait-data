


################
# Figure 4


library(rWCVP)
library(tidyverse)
library(sf)
library(conflicted)
library(cowplot)

richness <- st_read("Desktop/Projects_Git/lianas_GEB/data_prepared/liana_countries.shp") %>% 
  st_make_valid()


# relative richness of active to passive especies

habit_names <- c(
  `ric_prop_active` = "Active",
  `ric_prop_passive` = "Passive"
)


fig4_a <- richness %>% 
  st_drop_geometry() %>% 
  rename(mechanism = mechnsm, species = specs_d, LEVEL3_NAM = LEVEL3_N) %>% 
  group_by(LEVEL3_NAM, habit) %>% 
  summarise(species_richness = n_distinct(species)) %>% 
  pivot_wider(names_from = habit, values_from = species_richness) %>% 
  mutate(
    ric_ratio = active / passive,
    ric_total = coalesce(active, 0) + coalesce(passive, 0),
    ric_prop_active = (active / ric_total) * 100, # relative abundance measured as proportion representation of each group by region
    ric_prop_passive = (passive / ric_total) * 100
  ) %>% 
  select(LEVEL3_NAM, ric_prop_active, ric_prop_passive) %>% 
  pivot_longer(cols = ric_prop_active:ric_prop_passive,
               names_to = "habit") %>% 
  dplyr::filter(habit != "NA") %>% 
  right_join(rWCVPdata::wgsrpd3, by = "LEVEL3_NAM") %>% 
  select(-c(LEVEL2_COD, LEVEL1_COD)) %>% 
  dplyr::filter(habit != "NA") %>% 
  ungroup() %>% 
  st_as_sf() %>% 
  ggplot()+
  geom_sf(data = rWCVPdata::wgsrpd3) +
  geom_sf(aes(fill = value)) +
  scale_fill_viridis_c(name = "Active\nrichness")+ 
  scale_colour_viridis_c(name = "Passive\nrichness")+ 
  #remove extra whitespace
  coord_sf(expand=FALSE) +
  theme(legend.title = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold")) +
  facet_wrap(~ habit, nrow = 1, labeller = as_labeller(habit_names)) +
  labs(title = "(a) Proportional richness of climbing mechanisms")


# relative contribution of active and passive species per botanical continent

fig4_b <- richness %>% 
  st_drop_geometry() %>% 
  rename(mechanism = mechnsm, species = specs_d, LEVEL3_NAM = LEVEL3_N) %>% 
  right_join(rWCVPdata::wgsrpd3, by = "LEVEL3_NAM") %>% 
  select(-LEVEL2_COD) %>%
  group_by(LEVEL1_COD, habit) %>% 
  summarise(species_richness = n_distinct(species)) %>% 
  pivot_wider(names_from = habit, values_from = species_richness) %>% 
  mutate(
    ric_total = coalesce(active, 0) + coalesce(passive, 0),
    ric_prop_active = (active / ric_total) * 100, # relative abundance measured as proportion representation of each group by region
    ric_prop_passive = (passive / ric_total) * 100
  ) %>% 
  select(LEVEL1_COD, ric_prop_passive, ric_prop_active) %>% 
  mutate(bot_continent = case_when(
    LEVEL1_COD == 1 ~ "Europe",
    LEVEL1_COD == 2 ~ "Africa",
    LEVEL1_COD == 3 ~ "Asia Temperate",
    LEVEL1_COD == 4 ~ "Asia Tropical",
    LEVEL1_COD == 5 ~ "Australasia",
    LEVEL1_COD == 6 ~ "Pacific",
    LEVEL1_COD == 7 ~ "Northern America",
    LEVEL1_COD == 8 ~ "Sourthern America",
    LEVEL1_COD == 9 ~ "Antartic",
  )) %>% 
  pivot_longer(cols = ric_prop_passive:ric_prop_active,
               names_to = "habit") %>% 
  dplyr::filter(!(bot_continent %in% c("Pacific", "Antartic"))) %>% 
  # colors: passive gold2, active dodgerblue3
  ggplot(aes(x = factor(bot_continent), y = value, fill = habit)) +
  geom_bar(position="stack", stat="identity", width = 0.6) +
  labs(y = "Relative (%)", x = NULL) +
  ggtitle("(b) Contribution of climbing mechanisms to botanical continents") +
  scale_fill_manual(values = c("ric_prop_active" = "dodgerblue3", "ric_prop_passive" = "gold2"),
                    labels = c("ric_prop_active" = "Active", "ric_prop_passive" = "Passive")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.y =element_text(colour=c("black"), size=10),
        axis.text.x =element_text(colour=c("black"), size=10))



fig4 <- plot_grid(fig4_a, fig4_b,
                  ncol = 1, nrow = 2, rel_widths = c(1,1,1,1))

ggsave(fig4, filename = "fig4_2023.png", bg = "white",
       path = "Desktop/Projects_Git/lianas_GEB/figures/",
       height = 9, width = 15, scale = .75)

# END
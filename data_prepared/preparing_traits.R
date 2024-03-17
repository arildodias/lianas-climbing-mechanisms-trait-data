

# This is code is adapted from M. Paola Barajas Barbosa et al (2022) Assembly of functional diversity 
# in an oceanic island flora


# NOTE: We used phylogenetic trait imputation to estimate missing trait values. 
# To run analysis and figures of this study, trait imputation has to be run first

# Loading libraries

library(here)
library(janitor)
library(tidyverse)
library(reshape2)
library(missForest)
library(adephylo)
library(phytools)
library(ape)
library(pez)


# Tratis data
lianas_plot <- read.csv(here::here("Dropbox", "Arildo", "R", "Gentry-transects", "lianascomtax.csv"),
                        stringsAsFactors = FALSE) %>% 
  rename(species = newSpecies)

# Accepted taxon names
final_matches <- read.csv(here::here("Desktop", "Projects_Git", "lianas_GEB", "data_prepared", 
                                     "lianas_taxaname_matched.csv"), stringsAsFactors = FALSE)

# Phylogenetic tree
tree_lianas <- read.tree("Desktop/Projects_Git/lianas_GEB/data_prepared/phylotree_lianas.txt")


#  ---- Check and clean the traits dataset to imputation ----

# Calculate the total number of traits without NA values for distinct species by climbing 
# mechanism (habit) but without including only species with complete information for multiple traits
lianas_plot %>%
  group_by(habit, species) %>% 
  summarise(across(c(sla:wd), ~ n_distinct(species[!is.na(.)]))) %>%  
  ungroup() %>%
  group_by(habit) %>%
  summarise(across(c(sla:wd), sum))

# Number of species with information for multiple traits per climbing mechanism
lianas_plot %>% 
  group_by(habit, species) %>% 
  dplyr::select(habit, species, sla, la, leafnitro, wd) %>% 
  drop_na() %>% 
  summarise(across(c(sla:wd), ~ n_distinct(species[!is.na(.)]))) %>% 
  ungroup() %>%
  group_by(habit) %>%
  summarise(across(c(sla:wd), sum))

# For species with multiple data available for the same trait we use the average per species
lianas_df <- lianas_plot %>% 
  dplyr::select(habit, sla, la, seed, leafnitro, wd, species, genus, family) %>% 
  group_by(habit, species) %>% 
  summarise(across(c(sla:wd), ~ mean(., na.rm = TRUE)))

# Final trait dataset
lianas_traits <- final_matches %>% 
  left_join(lianas_df, by = c("scientific_name" = "species")) %>% 
  dplyr::select(habit, accepted_taxon_name, accepted_genus_name, accepted_family_name,
                sla:wd) 

# Two species are classified with both types active and passive so we correct manually and keep only one 
# Adenopodia scelerata has spines and therefore we classified as Hook/passive mechanisms. 
# Dichapetalum sessiliflorum is a scrambler and we keep only the species classified as passive

traits <- lianas_traits %>%
  dplyr::filter(!(accepted_taxon_name == "Adenopodia scelerata" & habit == "active"),
                !(accepted_taxon_name == "Dichapetalum sessiliflorum" & habit == "active")) %>%
  # to match with the phylogenetic tree
  mutate(accepted_taxon_name = str_replace_all(accepted_taxon_name, " ", "_")) %>% 
  rename(sp = accepted_taxon_name) %>% 
  mutate(IDmaster = 1:nrow(.), .before = 1) %>% 
  mutate(across(c(sla:wd), ~ ifelse(. == 0, NA, .)))


length(tree_lianas$tip.label)
sp_tree = as.data.frame(tree_lianas$tip.label)
duplicated(sp_tree)

traits$dupes = duplicated(traits$sp)
traits = filter( traits, dupes == FALSE)

length( setdiff(sp_tree$`tree_lianas$tip.label`, traits$sp) )
length( setdiff(traits$sp, sp_tree$`tree_lianas$tip.label`) )

traits = traits %>% 
  filter(sp %in% sp_tree$`tree_lianas$tip.label`) %>% 
  select(sp, sla, la, seed, leafnitro, wd) %>%
  column_to_rownames( var = "sp")


# Calculate Moran eigenvetors to use phylogenetic correlation structure to predict traits
# Eigenvectors eliminate features that have a strong correlation between them and also help in reducing over-fitting.

# Create phylogenetic proximity table. 
prox.Ab.all <- proxTips(tree_lianas, method = "Abouheif", normalize="none")
dim(prox.Ab.all) # 702 - the number of species

prox <- prop.table(prox.Ab.all, 1) 
prox <- 0.5 * (prox + t(prox))    

ME <- me.phylo(prox = prox)  
ME <- ME[rownames(traits),]

# To choose the number of eigeinvectors to include in the analysis is needed to test what is the number
# of eigenvectors to reduce the error. The recomendation is including fewer than 15–25 variables.
# By selecting 1:30 covers this maximum of variables

trait.imp <- cbind(traits, ME[,1:30]) # Morans Eigenvectors plus species by traits matrix (with column_to_rownames). 


#  ---- Phylogenetic imputation using missForest ----
set.seed(8) 
dfk <- data.frame(matrix(NA, nrow = 30, ncol = 6)) 
colnames(dfk) <- c("k", "OOB_sla","OOB_la","OOB_seed","OOB_leafnitro","OOB_wd")  

for (n in 1:30) {
  dfimp <- trait.imp[, 1: ( 5 +n)]  
  o <- missForest(dfimp, maxiter = 25, ntree = 100 , variablewise = TRUE) 
  dfk[n, 1] <- n
  dfk[n,2] <- o$OOBerror[1]
  dfk[n,3] <- o$OOBerror[2]    
  dfk[n,4] <- o$OOBerror[3]
  dfk[n,5] <- o$OOBerror[4]
  dfk[n,6] <- o$OOBerror[5]
  
}

dfk2<-dfk %>%                                                                          
  summarize(min_sla = min(OOB_sla), k_min_sla = k[which.min(OOB_sla)], 
            min_la = min(OOB_la), k_min_la = k[which.min(OOB_la)],
            min_seed = min(OOB_seed), k_min_seed = k[which.min(OOB_seed)],
            min_leafnitro = min(OOB_leafnitro), k_min_leafnitro = k[which.min(OOB_leafnitro)],
            min_wd = min(OOB_wd), k_min_wd = k[which.min(OOB_wd)]
  )

# phylogenetically-informed imputed data have generally similar/lower error rates
# Chose the number of eigenvectors that minimizes the imputation error (i.e., OOB_trait).

set.seed(8) 

sla_ideal   <-missForest(trait.imp[, 1: (10+dfk2$k_min_sla)], maxiter = 25, ntree = 100 ,   variablewise = TRUE)  
la_ideal <-missForest(trait.imp[, 1: (10+dfk2$k_min_la)], maxiter = 25, ntree = 100 ,   variablewise = TRUE) 
seed_ideal  <-missForest(trait.imp[, 1: (10+dfk2$k_min_seed)], maxiter = 25, ntree = 100 ,   variablewise = TRUE) 
leafnitro_ideal  <-missForest(trait.imp[, 1: (10+dfk2$k_min_leafnitro)], maxiter = 25, ntree = 100 ,   variablewise = TRUE) 
wd_ideal <-missForest(trait.imp[, 1: (10+dfk2$k_min_wd)], maxiter = 25, ntree = 100 ,   variablewise = TRUE) 


best_sla     <-tibble(sla=sla_ideal$ximp$sla, species=rownames(sla_ideal$ximp))
best_la   <-tibble(la=la_ideal$ximp$la, species=rownames(la_ideal$ximp))
best_seed    <-tibble(seed=seed_ideal$ximp$seed, species=rownames(seed_ideal$ximp))
best_leafnitro    <-tibble(leafnitro=leafnitro_ideal$ximp$leafnitro, species=rownames(leafnitro_ideal$ximp))
best_wd  <-tibble(wd=wd_ideal$ximp$wd, species=rownames(wd_ideal$ximp))


impute_out<-left_join(best_sla,best_la,by="species")%>%   
  left_join(.,best_seed, by="species")%>%
  left_join(.,best_leafnitro,by="species")%>%
  left_join(.,best_wd, by="species")%>%
  select(., species,sla,la,seed,leafnitro,wd)

# orginal trait data
traits_orig <- lianas_traits %>%
  dplyr::filter(!(accepted_taxon_name == "Adenopodia scelerata" & habit == "active"),
                !(accepted_taxon_name == "Dichapetalum sessiliflorum" & habit == "active")) %>%
  mutate(IDmaster = 1:nrow(.), .before = 1) %>% 
  mutate(across(c(sla:wd), ~ ifelse(. == 0, NA, .)))


Trait_imputed_phylo <- impute_out 
Trait_imputed_phylo$accepted_taxon_name <- gsub("_" , " "    , Trait_imputed_phylo$species)
Trait_imputed_phylo$species = NULL
z <- traits_orig[, c("IDmaster" , "accepted_taxon_name")]

Trait_imputed_phylo <- left_join(Trait_imputed_phylo, z,  by = "accepted_taxon_name") 


Trait_imputed_phylo$IDmaster <- as.character(Trait_imputed_phylo$IDmaster)

#  ---- Naive trait imputation ----
# Here we do not used phylogenetic information. 

traits_1 <- traits_orig[, c("IDmaster","sla","la","seed", "leafnitro", "wd")] 
traits_1 <- column_to_rownames(traits_1, var ="IDmaster") 

naive <- missForest(traits_1, maxiter = 25, ntree = 100 , variablewise = TRUE) 

Trait_imputed_naive <- naive$ximp 

Trait_imputed_naive["IDmaster"] <- rownames(Trait_imputed_naive)

write.csv(Trait_imputed_naive, "Desktop/Projects_Git/lianas_GEB/data_prepared/Trait_imputed_naive.csv",row.names=F) 



Trait_imputed_phylo_1 = left_join(Trait_imputed_naive, Trait_imputed_phylo, by = "IDmaster" )

summary(Trait_imputed_phylo_1)

Trait_imputed_phylo_2 <- Trait_imputed_phylo_1 %>%   
  mutate(sla.y = coalesce(sla.y, sla.x)) %>%
  mutate(la.y = coalesce(la.y, la.x)) %>%
  mutate(seed.y = coalesce(seed.y, seed.x)) %>%
  mutate(leafnitro.y = coalesce(leafnitro.y, leafnitro.x)) %>%
  mutate(wd.y = coalesce(wd.y, wd.x)) 

Trait_imp_phylo <- Trait_imputed_phylo_2[, 1:11]
Trait_imp_phylo <- data.frame(sla = Trait_imp_phylo$sla.y, 
                              la = Trait_imp_phylo$la.y, 
                              seed = Trait_imp_phylo$seed.y, 
                              leafnitro = Trait_imp_phylo$leafnitro.y, 
                              wd = Trait_imp_phylo$wd.y, 
                              IDmaster = as.integer(Trait_imp_phylo$IDmaster), 
                              IDmaster = Trait_imp_phylo$IDmaster)


Trait_imp_phylo$IDmaster.1 = NULL

write.csv(Trait_imp_phylo, "Desktop/Projects_Git/lianas_GEB/data_prepared/Trait_imputed_phylo8.csv", row.names=F) 


#  ---- Imputation errors ----
# Naive errors
OOB_multi<-data.frame(rbind(naive$OOBerror))[,1:5]
colnames(OOB_multi) <- c("OOB_sla","OOB_la","OOB_seed","OOB_leafnitro","OOB_wd")
OOB_multi$trait<-"naive"   

errors_naive <-OOB_multi
errors_phylo_info <- dfk2
errors_phylo_info <- errors_phylo_info[c ("min_sla", "min_la", "min_seed", "min_leafnitro", "min_wd")]

# Phylogenetically informed errors
errors_phylo_info <- data.frame(OOB_sla  = errors_phylo_info$min_sla, 
                                OOB_la = errors_phylo_info$min_la, 
                                OOB_seed = errors_phylo_info$min_seed,
                                OOB_leafnitro = errors_phylo_info$min_leafnitro, 
                                OOB_wd = errors_phylo_info$min_wd)

errors_phylo_info$trait <- "phylo"
OOB_errors = rbind(errors_phylo_info, errors_naive)
OOB_errors

# END





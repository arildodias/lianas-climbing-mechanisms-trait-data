


# This is code is adapted from M. Paola Barajas Barbosa et al (2022) Assembly of functional diversity 
# in an oceanic island flora

require(tidyverse)
require(phytools)
require(picante)
require(pez)
require(V.PhyloMaker2)
require(TNRS)
require(caper)
library(here)

#----------------------#
# Load data #
#----------------------#

traits_orig <- read.csv("Desktop/Projects_Git/lianas_GEB/data_prepared/traits_orig.csv")

#----------------------------#
# Build phylogenies          #
#----------------------------#

trts_spp <- traits_orig %>% 
  dplyr::select(., accepted_taxon_name, accepted_genus_name, accepted_family_name) %>% 
  distinct(.) %>% 
  data.frame(.)

ten_phyo<-phylo.maker(trts_spp,tree = GBOTB.extended.TPL, scenarios="S2",r=100)

#---------------------------#
# Phylogenetic signal       #
# of each traits            #
#---------------------------#

ten_LA <- traits_orig %>% 
  dplyr::select(accepted_taxon_name, la) %>% 
  mutate(accepted_taxon_name=str_replace_all(accepted_taxon_name, " ","_")) %>% 
  drop_na(.) %>% 
  tibble::column_to_rownames(., var="accepted_taxon_name")

ten_SLA <- traits_orig %>% 
  dplyr::select(accepted_taxon_name, sla) %>% 
  mutate(accepted_taxon_name=str_replace_all(accepted_taxon_name, " ","_")) %>% 
  drop_na(.) %>% 
  tibble::column_to_rownames(., var="accepted_taxon_name")

ten_Nmass <- traits_orig %>% 
  dplyr::select(accepted_taxon_name, leafnitro) %>% 
  mutate(accepted_taxon_name=str_replace_all(accepted_taxon_name, " ","_")) %>% 
  drop_na(.) %>% 
  tibble::column_to_rownames(., var="accepted_taxon_name")

ten_SM <- traits_orig %>% 
  dplyr::select(accepted_taxon_name, seed) %>% 
  mutate(accepted_taxon_name=str_replace_all(accepted_taxon_name, " ","_")) %>% 
  drop_na(.) %>% 
  tibble::column_to_rownames(., var="accepted_taxon_name")

ten_SSD <- traits_orig %>% 
  dplyr::select(accepted_taxon_name, wd) %>% 
  mutate(accepted_taxon_name=str_replace_all(accepted_taxon_name, " ","_")) %>% 
  drop_na(.) %>% 
  tibble::column_to_rownames(., var="accepted_taxon_name")


PhySig_out <- list(); 

for (i in 1:100) {
  # LA
  phy_up_la <- drop.tip(ten_phyo$scenario.2[[i]], 
                        setdiff(ten_phyo$scenario.2[[i]]$tip.label, rownames(ten_LA)))
  
  la_match<-match.phylo.data(phy_up_la, ten_LA)
  la_phy<-phylosig(la_match$phy, la_match$data[,1], method="lambda", test=FALSE, nsim=1000)
  la_phy<-cbind.data.frame(Trait="LA",Lambda=la_phy$lambda)
  
  # SLA
  phy_up_sla <- drop.tip(ten_phyo$scenario.2[[i]], 
                         setdiff(ten_phyo$scenario.2[[i]]$tip.label, rownames(ten_SLA)))
  
  sla_match<-match.phylo.data(phy_up_sla, ten_SLA)
  sla_phy<-phylosig(sla_match$phy, sla_match$data[,1], method="lambda", test=FALSE, nsim=1000)
  sla_phy<-cbind.data.frame(Trait="SLA",Lambda=sla_phy$lambda)
  
  # Nmass
  phy_up_nmass <- drop.tip(ten_phyo$scenario.2[[i]], 
                           setdiff(ten_phyo$scenario.2[[i]]$tip.label, rownames(ten_Nmass)))
  
  nmass_match<-match.phylo.data(phy_up_nmass, ten_Nmass)
  nmass_phy<-phylosig(nmass_match$phy, nmass_match$data[,1], method="lambda", test=FALSE, nsim=1000)
  nmass_phy<-cbind.data.frame(Trait="Nmass",Lambda=nmass_phy$lambda)
  
  # SM
  phy_up_sm <- drop.tip(ten_phyo$scenario.2[[i]], 
                        setdiff(ten_phyo$scenario.2[[i]]$tip.label, rownames(ten_SM)))
  
  sm_match<-match.phylo.data(phy_up_sm, ten_SM)
  sm_phy<-phylosig(sm_match$phy, sm_match$data[,1], method="lambda", test=FALSE, nsim=1000)
  sm_phy<-cbind.data.frame(Trait="SM",Lambda=sm_phy$lambda)
  
  # SSD
  phy_up_ssd <- drop.tip(ten_phyo$scenario.2[[i]], 
                         setdiff(ten_phyo$scenario.2[[i]]$tip.label, rownames(ten_SSD)))
  
  ssd_match<-match.phylo.data(phy_up_ssd, ten_SSD)
  ssd_phy<-phylosig(ssd_match$phy, ssd_match$data[,1], method="lambda", test=FALSE, nsim=1000)
  ssd_phy<-cbind.data.frame(Trait="SSD",Lambda=ssd_phy$lambda)
  
  a<-bind_rows(la_phy, sla_phy) 
  a<-bind_rows(a, nmass_phy) 
  a<-bind_rows(a,sm_phy )  
  out<-bind_rows(a,ssd_phy) 
  out$Iteration<-i
  cat("progress", i, sep=' ','\n')
  PhySig_out[[i]]<-rbind.data.frame(out)
}

PhySig_out<-bind_rows(PhySig_out)

PhySig_out %>% 
  group_by(Trait) %>% 
  summarise(avg_lambda = mean(Lambda))

# END
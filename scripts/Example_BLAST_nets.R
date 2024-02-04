library(tidyverse)
library(igraph)
setwd('~/Documents/nothingButNet')
p2d <- file.path('data/raw')
files <- list.files(p2d)

source('scripts/functions.R')

blast <- read.csv(file.path(p2d, files[grep('tally_blast', files)])) %>% 
  group_by(Sample) %>% 
  slice_max(Prcnt_seqs, n = 10) %>% 
  select(Sample, Taxon  = taxid, Prcnt_seqs) %>% 
  mutate('Method' = 'BLAST') %>% 
  drop_na() %>% 
  ungroup() 

blast_model <- read.csv(file.path(p2d, files[grep('Post.*BLAST', files)])) %>% 
  select(Sample, Taxon = TAXON_NEW, Prcnt_seqs) %>% 
  group_by(Sample) %>% 
  slice_max(Prcnt_seqs, n = 10) %>% 
  mutate('Method' = 'Modelled') %>% 
  drop_na() %>% 
  distinct()%>% 
  ungroup() 
  
blast_expert <- read.csv(file.path(p2d, files[grep('Integrated_Corbiculae', files)])) %>% 
  select(Sample = sample.id, Taxon, Prcnt_seqs = Percent_Sample) %>% 
  mutate('Method' = 'Expert') %>% 
  drop_na() %>% 
  distinct()

arranged_plants <- read.csv(file.path(p2d, files[grep('species', files)])) %>% 
  pull(abbreviation)
plants_legend <- read.csv(file.path(p2d, files[grep('species', files)])) 

arranged_bees <- unlist(strsplit(readLines(
  file.path(p2d, files[grep('bees', files)])), split = ',')) 

seqs <- rbind(blast, blast_model, blast_expert) %>% 
  mutate(Sample = as.numeric(Sample)) 

corbiculae <- read.csv(file.path(p2d, files[grep('Corbiculae_Samples', files)])) %>% 
  select(-site, -plant.species) 

seqs <- left_join(seqs, corbiculae, by = c('Sample' = 'sample.id'))

p2d <- ('~/Documents/assoRted/floral_observations/data')
files <- list.files(p2d) 

weeks <- read.csv(
  paste0(p2d, '/', files[grep('queen_observations', files)]) )  %>% 
  select(date, week) %>% 
  distinct(date, .keep_all = T) %>% 
  mutate(period = case_when(
    week %in% 3:5 ~ 'Early',
    week %in% 6:8 ~ 'Mid',
    week %in% 9:11 ~ 'Late'
  ))

seqs <- left_join(seqs, weeks, by = 'date') %>% 
  filter(caste == 'q') %>% 
  drop_na() %>% 
  mutate(bombus.species = gsub('Bombus.|.dark|.orange',  "", bombus.species)) %>% 
  select(Taxon, period, Prcnt_seqs, bombus.species, date, Method, Sample)

rm(blast, blast_model, blast_expert, files, p2d, corbiculae, weeks)


bee_obs_wk <- seqs %>% 
  ungroup() %>% 
  group_by(bombus.species, period, Method) %>% 
  mutate( Samples_time = length(unique(Sample)) ) %>% 
  group_by(bombus.species, period, Method, Taxon) %>% 
  mutate( Interactions_total  = sum(Prcnt_seqs)/Samples_time )  %>% 
  distinct(Taxon, bombus.species, period, .keep_all = T) %>% 
  select(bombus.species, Taxon, period, Interactions_total) %>% 
  ungroup(Taxon) 

resin <- bee_obs_wk %>%  # if sp. we need to leave genus spelled out, if epithet present, than abbreviate species. 
  ungroup() %>% 
  mutate(periods = str_count(Taxon, pattern =  fixed('.'))) %>% 
  mutate(Taxon = str_replace(Taxon, " ", "."),
         Taxon = if_else(periods == 0,  paste0(substr(Taxon, 1,1), '.',  
                                               gsub("^.*?\\.", "", Taxon )),
                         Taxon)) %>% 
  select(-periods) %>% 
  split(f = .$Method) %>% 
  map(~ .x  %>% 
        pivot_wider(names_from = Taxon, values_from = Interactions_total)
  ) %>% 
  bind_rows() %>% 
  split(f = list(.$Method, .$period) ) %>% 
  lapply(., select, -period, -Method) %>% 
  map(~ .x  %>% 
        column_to_rownames(., var = 'bombus.species') %>% 
        mutate(across(.cols = everything(), ~replace_na(.x, 0)))) %>% 
  map(., function(x){t(scales::rescale(t(x))) * 100}) %>% # scale big datasets. 
  map(., arrange_nets, col_arrange = arranged_plants, # now use the ARRANGE NET fun
      row_arrange = arranged_bees) # to remove empty plants

list2env(resin,env = environment())


exp <- bee_obs_wk %>% 
  filter(Method == 'Expert') %>% 
  ungroup() %>% 
  group_by(bombus.species, Taxon) %>% 
  mutate(
         Taxon = str_replace(Taxon, " ", "."),
         Taxon = paste0(substr(Taxon, 1,1), '.',
                                       gsub("^.*\\.", "", Taxon ))) %>%
  mutate(Taxon = case_when(
    Taxon == 'M.cilita' ~ 'M.ciliata',
    Taxon == 'L.lanszwertii' ~ 'L.leucanthus',
    Taxon == 'E.Ericaceae' ~ 'Ericaceae.sp.',
    Taxon == 'S.Salix' ~ 'Salix.sp.',
    Taxon == 'E.Epilobium' ~ 'Epilobium.sp.',
    Taxon == 'F.speciosus' ~ 'F.speciosa',
    TRUE ~ Taxon
  )) %>% 
  mutate(
    Interactions_total = sum(Interactions_total),
  ) %>% 
  ungroup(Taxon) %>% 
  mutate(
    Interactions_total = ( Interactions_total / sum(Interactions_total)) * 100) %>% 
  select(-period, -Method) %>% 
  distinct() %>% 
  pivot_wider(names_from = Taxon, values_from = Interactions_total) %>% 
  column_to_rownames(., var = 'bombus.species') %>% 
  mutate(across(.cols = everything(), ~replace_na(.x, 0)) )  


out <- arrange_nets(exp, col_arrange = arranged_plants, row_arrange = arranged_bees)

graphDrawer(out, lbl_fnt = 14,
            edge_clr = 'lightseagreen',
            node_clrs  = c("#CEAB07", "deeppink2"),
            legend_items = c("Bombus", "Plant"),
            directory = 'NetworkGraphs/BLAST_nets',
            ntwrks_page = 6,
            col = 2
)




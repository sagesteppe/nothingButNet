library(tidyverse)

setwd('~/Documents/nothingButNet')
source('scripts/functions.R')

p2d <- 'data'
files <- list.files(file.path(p2d, 'rmbl'), pattern = 'csv')

flr_ranks <- read.csv(
  file.path(p2d, 'rmbl', files[grep('flower_ranks', files)] )) %>% 
  select(site, doy, week, tot.obs.length, habitat, plant.species:abun.rank)

observation_times <- read.csv(
  file.path(p2d, 'rmbl', files[grep('queen_observations', files)]) ) %>% 
  select(site, week, tot.obs.length) %>% 
  distinct(site, week, .keep_all = T)

bee_obs <- read.csv(
  file.path(p2d, 'rmbl', files[grep('queen_observations', files)]) ) %>% 
  filter(caste == 'q') %>% 
  mutate(species = gsub('\\..*$', '', species)) %>% 
  select(site, doy, week, tot.obs.length, species, plant.species,
         resource.coll:pollen.collected, fl.switch.1:fl.switch.2.resource.coll) %>% 
  mutate(across(resource.coll:fl.switch.2.resource.coll, ~ na_if(.x, ''))) %>% 
  pivot_longer(c(plant.species, fl.switch.1, fl.switch.2), 
               names_to = 'observation', values_to = 'plant.species') %>% 
  drop_na(plant.species) %>% 
  group_by(site, week, species, plant.species) %>% 
  count(name = 'interactions') %>% 
  left_join(., observation_times, by = c("site", "week")) %>% 
  mutate(period = case_when(
    week %in% 3:5 ~ 'Early',
    week %in% 6:8 ~ 'Mid',
    week %in% 9:11 ~ 'Late'
  ))


f <- list.files(file.path(p2d, 'raw'))
arranged_plants <- read.csv(file.path(p2d, 'raw', f[grep('species', f)])) %>%
  pull(abbreviation)

arranged_bees <- unlist(strsplit(readLines(
  file.path(p2d, 'raw', f[grep('bees', f)])), split = ','))

rm(p2d, files, f, observation_times)

bee_obs_wk <- bee_obs %>% 
  ungroup() %>% 
  group_by(species, plant.species, period) %>% 
  mutate(Interactions_total  = sum(interactions)) %>% 
  distinct(species, plant.species, period, .keep_all = T) %>% 
  select(species, plant.species, period, Interactions_total, tot.obs.length) %>% 
  ungroup(plant.species) %>% 
  mutate(n = sum(Interactions_total),
         Prop = Interactions_total/n)

resin <- bee_obs_wk %>% 
  select(-tot.obs.length, -n, -Prop) %>%
  mutate(plant.species = paste0(substr(plant.species, 1,1), '.',
                                gsub("^.*\\.", "", plant.species ))) %>% 
  pivot_wider(names_from = plant.species, values_from = Interactions_total) %>% 
  mutate(across(everything(), ~replace_na(.x, 0))) %>% 
  arrange(period, species) %>% 
  split(., f = .$period) %>% 
  lapply(., column_to_rownames, 'species') %>% 
  lapply(., select, -period) %>% 
  lapply(., as.matrix)

# networks

tet <- lapply(resin, arrange_nets, col_arrange = arranged_plants, 
              row_arrange = arranged_bees)
rm(bee_obs_wk, resin, bee_obs)

observations_early <- tet[['Early']]
observations_mid <- tet[['Mid']]
observations_late <- tet[['Late']]

graphDrawer(observations_late, lbl_fnt = 14,
            edge_clr = 'lightseagreen',
            node_clrs  = c("#CEAB07", "deeppink2"),
            legend_items = c("Bombus", "Plant"),
            ntwrks_page = 12,
            col = 3
)

rm(observations_early, observations_mid, observations_late, tet)


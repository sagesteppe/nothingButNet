library(here)
library(tidyverse)
library(igraph)
# set_here('~/Documents/floral_observations')  # move this to be the root of the project folder on your box
setwd('~/Documents/floral_observations')

p2d <- 'data'
files <- list.files(p2d, pattern = 'csv')

flr_ranks <- read.csv(
  paste0(p2d, '/', files[grep('flower_ranks', files)] )) %>% 
  select(site, doy, week, tot.obs.length, habitat, plant.species:abun.rank)

observation_times <- read.csv(
  paste0(p2d, '/', files[grep('queen_observations', files)]) ) %>% 
  select(site, week, tot.obs.length) %>% 
  distinct(site, week, .keep_all = T)

bee_obs <- read.csv(
  paste0(p2d, '/', files[grep('queen_observations', files)]) ) #%>% 
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
  left_join(., observation_times) %>% 
  mutate(period = case_when(
    week %in% 3:5 ~ 'Early',
    week %in% 6:8 ~ 'Mid',
    week %in% 9:11 ~ 'Late'
  ))

arranged_plants <- c('C.gunnisonii', 'E.grandiflorum',
                     
                     'A.columbianum', 'A.coerulea', "T.fendleri",
                     'D.barbeyi', 'D.nuttallianum',
                     
                     "L.leucanthus", 'L.bakeri', "V.americana",
                     
                     "E.triflora","P.pulcherrima",
                     
                     "A.lewisii",
                     "S.sp",
                     "V.praemorsa", 
                     "C.lanceolata",
                     
                     "D.pulchellum",
                     
                     "M.ciliata" ,"M.fusiformis",
                     "H.capitatum", "H.fendleri",
                     "A.urticifolia", "P.procera", "P.bracteosa",
                     
                     "C.parryi", "D.hoopesii", "H.quinquenervis", "L.bigelovii",
                     "E.speciosus", "S.integerrimus", "T.officinale", "W.arizonica",
                     
                     "D.involucrata",
                     "V.occidentalis",
                      
                     "F.speciosa",
                     'O.occidentalis'
)

arranged_bees <- c('bifarius', 'mixtus', 'rufocinctus', 'sylvicola',
                   'flavifrons',
                   'appositus', 'californicus', 'nevadensis', 
                   'unknown') 

rm(p2d, files, observation_times)

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
  mutate(plant.species = paste0(substr(plant.species, 1,1), '. ',
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

early <- tet[['Early']]
mid <- tet[['Mid']]
late <- tet[['Late']]

graphDrawer(late, lbl_fnt = 14,
            plot_name = 'observation-early',
            edge_clr = 'lightseagreen',
            node_clrs  = c("#CEAB07", "deeppink2"),
            legend_items = c("Bombus", "Plant"),
            fname = 'molecular-late', 
            ntwrks_page = 12,
            col = 3
)


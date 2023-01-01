library(tidyverse)

setwd('~/Documents/nothingButNet')
source('scripts/functions.R')

p2d <- 'data'
files <- list.files(p2d, recursive = T)

arranged_bees <- unlist(strsplit(readLines(
  file.path(p2d,  files[grep('bees', files)])), split = ',')) %>% 
  str_trim(.)

lkp <- data.frame(doy = 146:208, week = rep(3:11, each = 7)) %>% 
  mutate(date = as.Date(doy, origin = "2015-01-01"))

morpho <- read.csv(file.path(p2d, files[grep('morpho', files)])) %>%
  mutate(date = str_replace(date, '2019', '2015'),  
         date = as.Date(date),
         species = str_remove(species, 'Bombus.|.dark|.orange'), 
         morphotype = str_remove(morphotype, 'OTHER_|_[(].*$')) %>%
  left_join(., lkp, by = 'date') %>% 
  mutate(period = case_when(
    week %in% 3:5 ~ 'Early',
    week %in% 6:8 ~ 'Mid',
    week %in% 9:11 ~ 'Late'
  )) %>% 
  group_by(species, period, morphotype) %>% 
  summarize(Grains = sum(sum)) %>% 
  filter(! morphotype  %in% c('CONIFER', 'INDISTINGUISHABLE_CLUMP', 'NER', 
                              "CAN'T_DISTINGUISH", 'DEHYDRATED'),
         species != 'kirbiellus',
         Grains > 1) %>% 
  mutate(Percent = (Grains/sum(Grains)) * 100)  %>% 
  select(-Grains)

rm(lkp)

arranged_plants <- unique(morpho$morphotype)
arranged_plants <- sort(arranged_plants)


test <- morpho %>% 
  pivot_wider(names_from = morphotype, values_from = Percent) %>% 
  mutate(across(everything(), ~ if_else(.x < 1, 1, .x)),
    across(everything(), ~replace_na(.x, 0))) %>% 
  arrange(period, species) %>% 
  split(., f = .$period) %>% 
  lapply(., column_to_rownames, 'species') %>% 
  lapply(., select, -period) %>% 
  lapply(., as.matrix)

tet <- lapply(test, arrange_nets, col_arrange = arranged_plants, 
              row_arrange = arranged_bees)

palyn_early <- tet[['Early']]
palyn_mid <- tet[['Mid']]
palyn_late <- tet[['Late']]

graphDrawer(palyn_mid, lbl_fnt = 14,
            edge_clr = 'lightseagreen',
            node_clrs  = c("#CEAB07", "deeppink2"),
            legend_items = c("Bombus", "Plant"),
            ntwrks_page = 12,
            col = 3
)


blanker(palyn_mid)

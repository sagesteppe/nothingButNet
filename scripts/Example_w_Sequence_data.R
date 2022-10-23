library(here)
library(tidyverse)
library(igraph)
#set_here('~/Documents/floral_observations')  # move this to be the root of the project folder on your box
# i_am('~/Documents/floral_observations')
setwd('~/Documents/nothingButNet/data/raw')
p2d <- file.path('~/Documents/nothingButNet/data/raw')
files <- list.files(p2d, pattern = 'csv')


# BLAST IS MISSING SAMPLES!!!!!!!!!!!! LIKE 3- 12 ????!??!?!?!??!?!?!?!?

blast <- read.csv(file.path(p2d, files[grep('blast', files)])) %>% 
  select(sample.id, Taxon  = name, n = Seqs_pr_taxon) %>% 
  mutate('Method' = 'BLAST') %>% 
  drop_na()
bracken <- read.csv(file.path(p2d, files[grep('bracken', files)])) %>% 
  filter(taxonomy_lvl %in% c('S', 'G')) %>% 
  select(sample.id, Taxon, n = new_est_reads) %>% 
  mutate('Method' = 'Bracken')
kraken <- read.csv(file.path(p2d, files[grep('kraken', files)]))  %>% 
  select(sample.id = Sample, Taxon, n) %>% 
  mutate('Method' = 'Kraken')

arranged_plants <- read.csv(file.path(p2d, files[grep('species', files)])) %>% 
  pull(abbreviation)
plants_legend <- read.csv(file.path(p2d, files[grep('species', files)])) 

seqs <- rbind(blast, bracken, kraken) %>% 
  mutate(sample.id = as.numeric(sample.id)) 

corbiculae <- read.csv(file.path(p2d, files[grep('Corbiculae', files)])) %>% 
  select(-site, -plant.species) 

seqs <- left_join(seqs, corbiculae, by = 'sample.id')

p2d <- ('~/Documents/floral_observations/data')
files <- list.files(p2d, pattern = 'csv') 

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
  select(Taxon, period, n, bombus.species, date, Method)

rm(blast, bracken, kraken, files, p2d, corbiculae, weeks)

arranged_bees <- c('bifarius', 'mixtus', 'rufocinctus', 'sylvicola',
                   'flavifrons',
                   'appositus', 'californicus', 'nevadensis', 
                   'unknown') 


bee_obs_wk <- seqs %>% 
  group_by(bombus.species, Taxon, period, Method) %>% 
  mutate(Interactions_total  = sum(n)) %>% 
  distinct(Taxon, bombus.species, period, .keep_all = T) %>% 
  select(bombus.species, Taxon, period, Interactions_total) %>% 
  ungroup(Taxon) 

resin <- bee_obs_wk %>%  # if sp. we need to leave genus spelt out, if epithet present, than abbreviate species. 
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


graphDrawer(Kraken.Mid, lbl_fnt = 14,
            plot_name = 'Mid.BLAST',
            edge_clr = 'lightseagreen',
            node_clrs  = c("#CEAB07", "deeppink2"),
            legend_items = c("Bombus", "Plant"),
            ntwrks_page = 9,
            col = 3
)

 # PROBLEM CURRENTLY READING IN ACROSS INTENDED PATTERNS - 
# TOP ROW IS BLAST, AND NEXT ROW IS BRACKEN 

netPage(col_var = c('Kraken', 'Bracken', 'BLAST'), mainT = 'Comparision of Foraging
             Patterns from Three Sequence Alignment Algorithms',
              row_var = c('Early', 'Mid', 'Late'), sep = '.')


graph_dims <- function(ntwrks_page){
  
  col <- length(col_var)
  W = round(1984 / col, 0)
  H = round(2864 / (ntwrks_page/col), 0 )
  if(W > H){ W <- H} else {H <- W}
  
  dims <- list('W' = W, 
               'H' = H)
  return(dims)
}


#' Creates a simple abbreviation legend for graphs
#' 
#' 
#' @param values a vector of names 
#' @param colN number of columns ot spread legend across
#' @fname a file name for the legend, defaults to legend. 
#' @ntwrks_page the number of networks you plan to place on the page
#' 
legenDrawer <- function(values, colN, directory, fname, ntwrks_page){
  
  if(missing(directory)) { directory <- 'NetworkGraphs' }
  if(missing(fname)) {fname <- 'legend.png'}
  
  col_var <- rep('A',length(colN))
  dims <- graph_dims(ntwrks_page, col_var)

  lname <- file.path(directory, fname)
  
  grp <- ceiling(length(values)/colN)
  l <- (grp * colN) - length(values)
  padding <- rep("", l)
  
  values <- c(values, rep("", l))
  v <- matrix(data = values , nrow = grp, ncol = colN)
  tt2 <- gridExtra::ttheme_minimal(core=list(fg_params=list(hjust=0, x=0.1)),
                                   rowhead=list(fg_params=list(hjust=0, x=0)))
  
  p <- gridExtra::tableGrob(v, theme = tt2)
  
  png(lname, height = dims$H, width = 1984)
  gridExtra::grid.arrange(p)
  dev.off()
  
}

legenDrawer(plants_legend$full, colN = 8)


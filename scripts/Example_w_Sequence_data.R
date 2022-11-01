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

# list2env(resin,env = environment())

graphDrawer(Kraken.Late, lbl_fnt = 14,
            plot_name = 'Kraken.Mid',
            edge_clr = 'lightseagreen',
            node_clrs  = c("#CEAB07", "deeppink2"),
            legend_items = c("Bombus", "Plant"),
            ntwrks_page = 9,
            col = 3
)

rm(seqs, bee_obs_wk, plants_legend)

netPage(col_var = c('Kraken', 'Bracken', 'BLAST'), mainT = 'Comparision of Foraging
             Patterns from Three Sequence Alignment Algorithms',
              row_var = c('Early', 'Mid', 'Late'), sep = '.')

plsl <- plants_legend %>% 
  filter(str_detect(full, '.sp.$', negate = T)) %>% 
  pull(full) 

tableLegend(plsl, LegcolN = 8, ntwrks_page = 9, colN = 3)

sizeLegend(x = resin, y.space = c(1, 1.35, 1.25, 2.25, 1.95), fill_col = 'black',
           title = 'No. of\nInteractions', colN = 3, ntwrks_page = 9)

categoryLegend(node_clrs  = c("#CEAB07", "deeppink2"), LcolN = 1,
                       legend_items = c("Bombus", "Plant"), ntwrks_page = 9, colN = 3)


tableLegend2 <- function(x, table_items, directory, fname, colors, legend_items, node_clrs,
                        ntwrks_page, colN, LcolN, fill_col, y.space, values, LegcolN, title){
  
  if(missing(directory)) { directory <- 'NetworkGraphs' }
  if(missing(fname)) {fname <- 'TableLegend.png'}
  if(missing(fill_col)){fill_col <- 'white'}
  if(missing(y.space)){y.space <- seq(from = 1, to = 2, length.out = 5)}
  
  dims <- graph_dims(ntwrks_page, col = colN)
  
  # create the categorical legend.

  d <- data.frame(Group = c('Bombus', 'Plant'),
                  Dummy = c(1,1), Dummy2 = c(1,1))
  
  v <- c("#CEAB07", "deeppink2")
  names(v) <- c('Bombus', 'Plant')
  
  p <- ggplot2::ggplot(d, 
                       ggplot2::aes(x= Dummy, y=Dummy2, color = Group)) + 
    ggplot2::geom_point() +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size=10))) +
    ggplot2::scale_fill_manual(v) +
    ggplot2::labs(color = 'Major Interacting\nGroups ') +
    ggplot2::theme(legend.key = ggplot2::element_rect(fill = "white"),
                   legend.title.align = 0.5,
                   legend.key.size = unit(1.5, 'cm'), 
                   legend.key.height = unit(3, 'cm'), 
                   legend.key.width = unit(2, 'cm'), 
                   legend.title = ggplot2::element_text(size=16),
                   legend.text = ggplot2::element_text(size=14)) 
  
  catLeg <- cowplot::get_legend(p)
  
  #  create the legend table. 
  lname <- file.path(directory, fname)
  grp <- ceiling(length(values)/LegcolN)
  l <- (grp * LegcolN) - length(values)
  
  values <- c(values, rep("", l))
  v <- matrix(data = values , nrow = grp, ncol = LegcolN)
  tt2 <- gridExtra::ttheme_minimal(core=list(fg_params=list(hjust=0, x=0.1)),
                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                   base_size = 12)
  
  tableLeg <- gridExtra::tableGrob(v, theme = tt2)
  
  # Create the bubble size legend
  net <- lapply(x, igraph::graph_from_incidence_matrix, weight = T)
  deg <- lapply(net, igraph::centr_degree,  mode = "all")
  
  vals <- vector(mode = 'list', length = length(deg))
  for (i in 1:length(deg)){
    vals[[i]] <- deg[[i]][['res']]
  }
  
  interaction_no <- Reduce(c,vals)
  breaks <- (max(interaction_no) - min(interaction_no)) / 4
  Intervals <- c(min(interaction_no), round(breaks * 1),
                 round(breaks * 2), round(breaks * 3), max(interaction_no)
  )
  vals_size <- 2.5 * sqrt(Intervals)
  size_legend <- data.frame('Interactions' = Intervals,  'Area' = vals_size)
  
  # prepare the figure 
  p <- ggplot2::ggplot(size_legend, ggplot2::aes(Interactions, Area, size = Area)) + 
    ggplot2::geom_point() +
    ggplot2::theme(legend.key = 
                      ggplot2::element_rect(fill = "white"),
                   legend.title.align = 0.5,
                   legend.key.size = unit(1, 'cm'), 
                   legend.key.height = unit(2, 'cm'), 
                   legend.key.width = unit(2, 'cm'), 
                   legend.title = ggplot2::element_text(size=16),
                   legend.text = ggplot2::element_text(size=14)) +
    ggplot2::labs(size = 'No. of Interacting\nspecies')
  
  sizeLeg <- cowplot::get_legend(p)
  
  # now combine the three legends. 
  grobBY <- list(tableLeg, sizeLeg, catLeg)
  
  lays <- rbind(
    c(1,1,1,1,1, 2),
    c(1,1,1,1,1, 2),
    c(1,1,1,1,1, 2),
    c(1,1,1,1,1, 3),
    c(1,1,1,1,1, 3)
  )
  
   bucket <- gridExtra::grid.arrange(grobs = grobBY, layout_matrix = lays)
    png(lname,
        width = 1984, height = dims$H, units = "px", pointsize = 12)
    invisible(grid::grid.draw(bucket))
    invisible(dev.off())
  
    message(paste0("'", fname, 
                   "' has been rendered as a legend and saved to:\n ",
                   file.path(directory, fname)))
}


tableLegend2(x = resin, node_clrs = c("#CEAB07", "deeppink2"), ntwrks_page = 9,
                 colN = 3, LcolN = 1, legend_items = c("Bombus", "Plant"), 
                 table_items = arranged_plants, values = arranged_plants, fill_col = 'black',
                 LegcolN = 8, title = 'bill walton', y.space = c(1, 1.35, 1.25, 2.25, 1.95))

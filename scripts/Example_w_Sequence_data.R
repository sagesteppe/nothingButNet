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



netPage(col_var = c('Kraken', 'Bracken', 'BLAST'), mainT = 'Comparision of Foraging
             Patterns from Three Sequence Alignment Algorithms',
              row_var = c('Early', 'Mid', 'Late'), sep = '.')

plsl <- plants_legend %>% 
  filter(str_detect(full, '.sp.$', negate = T)) %>% 
  pull(full) 

legenDrawer(plsl, colN = 8, ntwrks_page = 9, LegcolN = 5)


# working on making a legend of node size 












node_sizes <- function(x){
  
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
  
  size_legend <- data.frame(
    'Interactions' = Intervals, 
             'Area' = vals_size)
  
  a <- legend('center', 
              legend = size_legend$Interactions, 
              pt.cex = size_legend$Area/100, bty = 'n', 
              y.intersp = c(1,1.35,1.25,2.25,1.95),
              col='white', pch=21, pt.bg='white', 
              title = 'No. of\nInteractions')
  
  png(filename,
      width = dims$W, height = dims$H, units = "px", pointsize = 12)
  
  plot(er_graph, vertex.label=NA, vertex.color = NA, edge.color = NA, 
       vertex.frame.color = NA)
  legend('center', legend = size_legend$Interactions, 
         pt.cex = size_legend$Interactions, col='black', bty = 'n', 
         y.intersp= c(1,1.35,1.25,2.25,1.95), title = 'No. of\nInteractions'
  )
  symbols(a$text$x + a$rect$left, a$text$y, circles = size_legend$Area/100, 
          inches = FALSE, add = TRUE, bg = NA)
  
  invisible(dev.off())
  
  message(paste0("'", plot_name, 
                 "' has been rendered as a legend and saved to:\n ", filename))
  
}


size_legend <- node_sizes(resin)


a <- legend('center', 
            legend = size_legend$Interactions, 
            pt.cex = size_legend$Area/100, bty = 'n', 
            y.intersp = c(1,1.35,1.25,2.25,1.95),
            col='white', pch=21, pt.bg='white', 
            title = 'No. of\nInteractions')

png(filename,
    width = dims$W, height = dims$H, units = "px", pointsize = 12)

plot(er_graph, vertex.label=NA, vertex.color = NA, edge.color = NA, 
     vertex.frame.color = NA)
legend('center', legend = size_legend$Interactions, 
       pt.cex = size_legend$Interactions, col='black', bty = 'n', 
       y.intersp= c(1,1.35,1.25,2.25,1.95), title = 'No. of\nInteractions'
       )
symbols(a$text$x + a$rect$left, a$text$y, circles = size_legend$Area/100, 
        inches = FALSE, add = TRUE, bg = NA)

invisible(dev.off())

message(paste0("'", plot_name, 
               "' has been rendered as a legend and saved to:\n ", filename))

#' Draw a node color legend grob 
#' 
#' This quickly draws a grob to be assembled onto the 'legend' section of a page
#' it is quite minimal and does not currently offer much flexibility.
#' 
#' @param legend_items character vector with names of node items
#' @param node_clrs character vector with node item colors
#' @param colN number of columns to split legend across
#' @example category_legend_drawer(node_clrs  = c("#CEAB07", "deeppink2"), 
#' legend_items = c("Bombus", "Plant"))
#' @seealso legenDrawer
#' 
category_legend_drawer <- function(legend_items, node_clrs, colN){
  
  if(missing(colN)) {colN <- 1}
  
  er_graph <- igraph::erdos.renyi.game(100, 5/100) 
  plot(er_graph, vertex.label=NA, vertex.color = NA, edge.color = NA, 
       vertex.frame.color = NA)
  legend('center', legend_items, 
         pch=21, col="#777777", 
         pt.bg=node_clrs, bty = "n",
         pt.cex=2, cex=.8,  ncol = colN)
  
}

category_legend_drawer(node_clrs  = c("#CEAB07", "deeppink2"), 
                       legend_items = c("Bombus", "Plant"))

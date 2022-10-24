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


graphDrawer(Kraken.Late, lbl_fnt = 14,
            plot_name = 'Kraken.Mid',
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

tableLegend(plsl, LegcolN = 8, ntwrks_page = 9, colN = 3)

sizeLegend(x = resin, y.space = c(1, 1.35, 1.25, 2.25, 1.95), fill_col = 'black',
           title = 'No. of\nInteractions', colN = 3, ntwrks_page = 9)

categoryLegend(node_clrs  = c("#CEAB07", "deeppink2"), LcolN = 1,
                       legend_items = c("Bombus", "Plant"), ntwrks_page = 9, colN = 3)




##### GROBS SEEMS ALL FUCKED. I THINK THE APPROACH IS TO CREATE THE THREE LEGEND GROBS SIMULTANEOSULY IN R AND ARRANGE FROM THERE. ALL FUCKED WITH SAVING EXTERNALLY

#' Arrange the legends into a single table grob
#' 
#' Combine the three legends into a single table to be placed at the bottom of 
#' a page of networks. If you have been using default filenames and path this function
#' can be run without any arguments.
#' @param directoryIN
#' @param filenamesIN
#' @param directoryOUT
#' @param filenameOUT
#' @example 
#' @seealso 'sizeLegend', 'abbreviaitonTable', 'categoryLegend 

grobLegend <- function(directoryIN, filenamesIN, directoryOUT, filenameOUT, ntwrks_page,
                       colN){
  
  if(missing(directoryIN)) {directoryIN <- 'NetworkGraphs' }
  if(missing(directoryOUT)) {directoryOUT <- 'NetworkGraphs' }
  if(missing(filenamesIN))
    {filenamesIN <- c('TableLegend', 'CategoricalLegend', 'SizeLegend')}
  if(missing(filenameOUT)) {filenameOUT <- 'LegendGrob'}
  
  dims <- graph_dims(ntwrks_page, col = colN)
  
  filenamesIN <- file.path(directoryIN, paste0(filenamesIN, '.png'))
  images = lapply(filenamesIN, png::readPNG)
  grob_images = lapply(images, grid::rasterGrob)
  
  layout <- matrix(data = c(1,2,1,3), ncol = 2, byrow = T)
  gzup <- gridExtra::arrangeGrob(grobs = grob_images, layout_matrix = layout)
  
  png(file.path(directoryOUT, paste0(filenameOUT, '.png')),
      width = 1984, height = dims$H, units = "px", pointsize = 12)
  print(gzup)
  invisible(dev.off())
  
  message(paste0("'", filenameOUT, 
                 "' has been rendered as the legend panel and saved to:\n ", 
                 file.path(directoryOUT, filenameOUT)))
}

#grobLegend(ntwrks_page = 9, colN = 3)


list2env(grob_images, env = environment())


dims <- graph_dims(ntwrks_page = 9, col = 3)

filenamesIN <- file.path(directoryIN, paste0(filenamesIN, '.png'))
images = lapply(filenamesIN, png::readPNG)
grob_images = lapply(images, grid::rasterGrob)

layout <- matrix(data = c(1,2,1,3), ncol = 2, byrow = T)
gzup <- gridExtra::arrangeGrob(grobs = grob_images, layout_matrix = grob2)


gridExtra::grid.arrange(grob_images[[1]], 
                        gridExtra::arrangeGrob(grob_images[[2]],
                                               grob_images[[3]]), ncol = 2,
                        heights=c(4, 1))
                  #      layout_matrix = rbind(c(1,1,2), c(1,1,3), c(1,1,3), c(1,1,3)))


message(paste0("'", filenameOUT, 
               "' has been rendered as the legend panel and saved to:\n ", 
               file.path(directoryOUT, filenameOUT)))


r <- sapply(images, dim)
r

r[2,1] / r[2,2]

?arrangeGrob

grid.arrange(grob_images[[1]], 
             arrangeGrob(grob_images[[2]], 
                         grob_images[[3]], 
                         ncol = 1, heights = c(1,3)),
             ncol=2, heights = c(1, 1), 
             widths = c(5, 1))


grid.arrange(grob_images[[1]], grob_images[[3]],
             ncol=2, heights = c(1, 1), 
             widths = c(5, 1))

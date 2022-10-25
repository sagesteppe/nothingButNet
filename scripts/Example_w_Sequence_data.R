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

grobLegend <- function(directoryIN, filenamesIN, directoryOUT, filenameOUT, ntwrks_page,colN){
  
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


grid.arrange(grob_images[[1]], 
             arrangeGrob(grob_images[[2]], 
                         grob_images[[3]], 
                         ncol = 1, heights = c(1,3)),
             ncol=2, heights = c(1, 1), 
             widths = c(5, 1))


grid.arrange(grob_images[[1]], grob_images[[3]],
             ncol=2, heights = c(1, 1), 
             widths = c(5, 1))











































tableLegend <- function(x, table_items, directory, fname, colors, legend_items, node_clrs,
                        ntwrks_page, colN, LcolN, fill_col, y.space, values, LegcolN, title){
  
  if(missing(directory)) { directory <- 'NetworkGraphs' }
  if(missing(fname)) {fname <- 'TableLegend.png'}
  if(missing(fill_col)){fill_col <- 'white'}
  if(missing(y.space)){y.space <- seq(from = 1, to = 2, length.out = 5)}
  
  dims <- graph_dims(ntwrks_page, col = colN)
  
  # create the categorical legend.
  er_graph <- igraph::erdos.renyi.game(100, 5/100) 
  catLeg <- plot(er_graph, vertex.label=NA, vertex.color = NA, edge.color = NA, 
                 vertex.frame.color = NA) 
  legend('center', legend_items, 
           pch = 21, col="#777777", 
           pt.bg = node_clrs, bty = "n", 
           pt.cex = 4, cex = 2,  ncol = LcolN)
  
  dev.off()
  
  #  create the legend table. 
  lname <- file.path(directory, fname)
  grp <- ceiling(length(values)/LegcolN)
  l <- (grp * LegcolN) - length(values)
  
  values <- c(values, rep("", l))
  v <- matrix(data = values , nrow = grp, ncol = LegcolN)
  tt2 <- gridExtra::ttheme_minimal(core=list(fg_params=list(hjust=0, x=0.1)),
                                   rowhead=list(fg_params=list(hjust=0, x=0)))
  
  table_grob <- gridExtra::tableGrob(v, theme = tt2)
  
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
  Size <- plot(er_graph, vertex.label=NA, vertex.color = NA, edge.color = NA, 
               vertex.frame.color = NA) # empty dummy plot
  legend('center', legend = size_legend$Interactions, 
              pt.cex = size_legend$Area/100, bty = 'n', 
              y.intersp = y.space, cex = 2,
              col='white', pch=21, pt.bg='white', title = title)
  a <- legend('center', legend = size_legend$Interactions, 
                        pt.cex = size_legend$Interactions, col='black', bty = 'n', 
                        y.intersp = y.space, title = title, cex = 2,
  )
  Size_grob <- Size + symbols(a$text$x + a$rect$left, a$text$y, 
                                   circles = size_legend$Area/100, 
          inches = FALSE, add = TRUE, bg = fill_col)
  
  # assemble the super legend. 
  
  #legend_grob <- gridExtra::grid.arrange(table_grob, catLeg, Size_grob,
  #             ncol = 2, heights = c(1, 1), LcolN = 1, 
  #             widths = c(5, 1))
  
  #png(lname,
  #    width = 1984, height = dims$H, units = "px", pointsize = 12)
  
  #invisible(dev.off())
  
  message(paste0("'", fname, 
                 "' has been rendered as a legend and saved to:\n ",
                 file.path(directory, fname)))
  
  return(catLeg)
}


p <- tableLegend(x = resin, node_clrs = c("#CEAB07", "deeppink2"), ntwrks_page = 9,
            colN = 3, LcolN =1,
             legend_items = c("Bombus", "Plant"), table_items = arranged_plants,
            values = arranged_plants, LegcolN = 6, title = 'bill walton')

plot(p)

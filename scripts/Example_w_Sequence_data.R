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

tableLegend2(x = resin, node_clrs = c("#CEAB07", "deeppink2"), ntwrks_page = 9,
                 colN = 3, LcolN = 1, legend_items = c("Bombus", "Plant"), 
                 table_items = arranged_plants, fill_col = 'black',
                 LegcolN = 8, title = 'bill walton', y.space = c(1, 1.35, 1.25, 2.25, 1.95))


netPage2 <- function(directory, col_var, row_var, fname, sep_char, mainT, Tlegend_fname) {
  
  if(missing(directory)) { directory <- 'NetworkGraphs' }
  if(missing(fname)) {fname <- 'mosaiced_Nets.pdf'}
  if(missing(sep_char))  {sep_char <- ''}
  if(missing(mainT))  {mainT <- 'Bill walton thinks u forgot a title'}
  if(missing(Tlegend_fname)) {Tlegend_fname <- 'TableLegend.png'}
  
  rowOrder <- matrix( # this sets up a matrix to have rows groups in order
    data = rep(row_var, each = length(col_var)),
    ncol = length(col_var), byrow = T)
  colOrder <- matrix( # this sets up a matrix to have column groups in order
    data = rep(col_var, each = length(col_var)),
    ncol = length(col_var), byrow = F)
  image_orders <- matrix(paste0(colOrder, sep_char, rowOrder, '.png'), 
                         ncol = length(col_var))
  image_orders <- as.vector(t(image_orders)) 
  image_orders <- file.path(directory, image_orders)  
  images = lapply(image_orders, png::readPNG)
  grob_images = lapply(images, grid::rasterGrob)
  
  nPlots <- (length(col_var) * length(row_var))
  colN <- length(col_var)
  rowN <- nPlots / length(col_var)
  
  # define the layout
  layout <- matrix(nrow = rowN, ncol = colN, byrow = T,
                   data = c(rep(1:(rowN -1), each = colN),
                            rowN:((rowN-1) + colN) )
  )
  layout <- rbind(layout, rep((max(layout) + 1), nrow(layout)))
  
  # Recover the top grobs!!
  top_grobs <- split(grob_images[1:(length(grob_images) - rowN)],
                     ceiling(seq_along(grob_images[1:(length(grob_images) -
                                                        rowN)]) / rowN))
  names(top_grobs) <- paste0('t', seq(1:length(top_grobs)))
  
  # title up the real top grob
  t1 <- gridExtra::arrangeGrob(grobs = top_grobs$t1, ncol = colN, 
                               padding = unit(0.0, "line"),
                               left = row_var[1], top = mainT)
  
  # this will get all of the subsequent top grobs !!!!!!!!
  top_grobs <- mapply(gridExtra::arrangeGrob, 
                      grobs = top_grobs[2:length(top_grobs)], ncol = colN, 
                      left = row_var[2:(length(row_var)-1)], 
                      padding = unit(0.0, "line"))
  
  list2env(top_grobs,env = environment())
  
  # this recovers the bottom grobs
  L <- length(grob_images)
  range <- ((L - rowN) + 1 ): L
  bottom_grobs <- split(grob_images[((L - rowN) + 1 ): L],
                        seq_along(range) / 1)
  names(bottom_grobs) <- paste0('b', seq(1:length(bottom_grobs)))
  
  b1 <- gridExtra::arrangeGrob(grobs = bottom_grobs$b1, bottom = col_var[1],
                               left = row_var[length(row_var)],
                               padding = unit(0.0, "line"))
  bottoms <- mapply(gridExtra::arrangeGrob, 
                    grobs = bottom_grobs[2:length(bottom_grobs)],
                    bottom = col_var[2:(length(col_var))],
                    padding = unit(0.0, "line"))
  
  list2env(bottoms, env = environment())
  
  
  # ensure the grobs are in the correct order for the mosaic. 
  g2p <- mget(ls(pattern = '[t|b][1-9]{1}'))
  groborder <- c(paste0('t', seq(1:9)), paste0('b', seq(1:9)))
  ord2grab <- match(groborder, names(g2p)) |> na.omit()
  g2p <- g2p[ord2grab]
  
  # load and place the legend onto the grobs2plot
  legend <- grid::rasterGrob(png::readPNG(file.path(directory, Tlegend_fname)))
  legend <- gridExtra::arrangeGrob(legend, nrow = 1)
  g2p <- c(g2p, leg = list(legend)) 
  
  
  # place on the page and print.
  ml <- gridExtra::marrangeGrob(grobs = g2p, 
                                layout_matrix = layout, top = '')
  pdf(file = file.path(directory, fname), paper = 'a4')
  print(ml)
  invisible(dev.off())
  
  message(paste0("'", fname, 
                 "' has been rendered as a pdf and saved to:\n ",
                 file.path(directory, fname)))
 return(g2p) 
}

legend <- netPage2(col_var = c('Kraken', 'Bracken', 'BLAST'), mainT = 'Comparision of Foraging Patterns from Three Sequence Alignment Algorithms',
        row_var = c('Early', 'Mid', 'Late'), sep = '.')


str(a)
str(b)


t <- grid::textGrob(".")
g <- list(t, legend, t)

r1 <- c(1, rep(2, times = 8), 3)
t <- matrix( data = rep(r1, times = 10), nrow = 10, byrow = T)


a <- gridExtra::arrangeGrob(legend, nrow = 1)

grid::grid.draw(legend)

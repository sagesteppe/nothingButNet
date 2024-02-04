setwd('~/Documents/nothingButNet/scripts')

library(tidyverse)
library(igraph)
dat <- read.csv('../data/JeremieReclassTable.csv') %>% 
  mutate(value = if_else(value == 0, 0.1, value))
source('functions.R')

# we will create Three different networks. Where each one displays a different type of information


p2d <- file.path('../data/raw')
files <- list.files(p2d)

arranged_plants <- read.csv(file.path(p2d, files[grep('species', files)])) %>% 
  mutate(Genus = str_extract(full, '^.*[.]+?'), 
         Genus = str_remove(Genus, '[.].*$')) %>% 
  drop_na(Genus) %>% 
  filter(Genus != 'P') %>% 
  distinct(Genus) %>% 
  pull(Genus)

arranged_bees <- unlist(strsplit(readLines(
  file.path(p2d, files[grep('bees', files)])), split = ',')) 

rm(p2d, files)

out <- dat %>%
  split(f = list(.$method) ) %>% 
  map(~ .x  %>% 
        pivot_wider(names_from = plant, values_from = value)
  ) %>% 
  map(~ .x  %>% 
        column_to_rownames(., var = 'bee') %>% 
        mutate(across(.cols = everything(), ~replace_na(.x, 0))) %>% 
        select(-method)) %>% 
  map(., arrange_nets, col_arrange = arranged_plants, # now use the ARRANGE NET fun
      row_arrange = arranged_bees) # to remove empty plants
  
Observations <- out[[1]]
Morphology <- out[[2]]
Molecular <- out[[3]]


simple_grapher <- function(x){
  
  obs <- rbind(x, 
               setNames(data.frame(matrix(ncol = ncol(x), nrow = 3, data = 0)), 
                        colnames(x)))
  obs <- cbind(obs, data.frame(matrix(data = 0, nrow = nrow(obs), ncol = 3)))
  
  VNames <- vertex_names(obs)
  net = igraph::graph_from_incidence_matrix(obs, weight=T)
  deg = igraph::centr_degree(net, mode="all")
  
  net <- set_colors(x = obs, net, node_clrs  = c("#CEAB07", "deeppink2"), VNames = VNames)
  template <- layout_in_circle(net)
  
  V(net)$label.color <- 'black'
  V(net)$size = 3.2*sqrt(deg$res)
  E(net)$width = E(net)$weight/5 
  E(net)$width <- ifelse( E(net)$width < 0.75, 0.75, E(net)$width)# make so lines can be seen
  
  VLPs <- vert_lab_position(obs, col = 1)
  VLPs$dist_vals <- rep(3, length(VLPs$degree_shift))
  
  return(setNames(
    list(net, template, VNames, VLPs), c('net', 'template', 'VNames', 'VLPs')
    )
  )
}


ob <- simple_grapher(Observations)
VLPs <- ob[['VLPs']] ; net <- ob[['net']]; template <- ob[['template']]; VNames <- ob[['VNames']]
VLPs$degree_shift[6] <- -0.15 # californicus
VLPs$degree_shift[7] <- -0.3 # nevadensis
VLPs$degree_shift[11] <- 3.4 # aquilegia
VLPs$degree_shift[22] <- 3.05 # agastache
VLPs$degree_shift[23] <- 2.85 # Cynoglossum
VLPs$degree_shift[24] <- 2.7 # Mertensia
VLPs$degree_shift[25] <- 2.3 # hydrophyllum
VLPs$degree_shift[26] <- 1.4 # erigeron
VLPs$degree_shift[27] <- 0.3 # scabrethia
VLPs$degree_shift[28] <- 0.1 # Symphyotrichum


png('../NetworkGraphs/Manuscript/Observations_manual.png',
    width = 720, height = 720, units = "px", pointsize = 12)
par(mar=c(3.5,6.5,3.5,6.5))
plot(net, layout=template, 
     main = 'Observations',
     edge.color = 'lightseagreen',
     vertex.label.cex = 1.25,
     vertex.label = VNames, 
     vertex.label.dist= VLPs$dist_vals, 
     vertex.label.degree = VLPs$degree_shift
) 
dev.off()


# repeat for molecular

ob <- simple_grapher(Molecular)
VLPs <- ob[['VLPs']] ; net <- ob[['net']]
template <- ob[['template']]; VNames <- ob[['VNames']]

VLPs$degree_shift[6] <- -0.15 # californicus
VLPs$degree_shift[7] <- -0.3 # nevadensis
VLPs$degree_shift[22] <- 3.05 #  hydrophyllum
VLPs$degree_shift[23] <- 2.9 #  pedicularis
VLPs$degree_shift[24] <- 2.4 #  Erigeron
VLPs$degree_shift[25] <-  0.45 # Helianthella
VLPs$degree_shift[26] <-  0.3 # senecio # 0.3
VLPs$degree_shift[27] <- 0.15 #  Taraxacum

png('../NetworkGraphs/Manuscript/Molecular_manual.png',
    width = 720, height = 720, units = "px", pointsize = 12)
par(mar=c(3.5,6.5,3.5,6.5))
plot(net, layout = template, 
     main = 'Molecular',
     edge.color = 'lightseagreen',
     vertex.label.cex = 1.25,
     vertex.label = VNames, 
     vertex.label.dist= VLPs$dist_vals, 
     vertex.label.degree = VLPs$degree_shift
) 
dev.off()


# repeat for morphological

ob <- simple_grapher(Morphology)
VLPs <- ob[['VLPs']] ; net <- ob[['net']]
template <- ob[['template']]; VNames <- ob[['VNames']]

VLPs$degree_shift[6] <- -0.15 # californicus
VLPs$degree_shift[7] <- -0.3 # nevadensis
VLPs$degree_shift[20] <- 3.05 #  hydrophyllum
VLPs$degree_shift[21] <- 2.9 # Pedicularis
VLPs$degree_shift[22] <- 0.6 # Erigeron
VLPs$degree_shift[23] <- 0.25 # helianthella

png('../NetworkGraphs/Manuscript/Morphology_manual.png',
    width = 720, height = 720, units = "px", pointsize = 12)
par(mar=c(3.5,6.5,3.5,6.5))
plot(net, layout = template, 
     main = 'Morphology',
     edge.color = 'lightseagreen',
     vertex.label.cex = 1.25,
     vertex.label = VNames, 
     vertex.label.dist= VLPs$dist_vals, 
     vertex.label.degree = VLPs$degree_shift
) 
dev.off()




graphDrawer(Observations, lbl_fnt = 14,
            fname = 'Method.Observations',
            edge_clr = 'lightseagreen',
            node_clrs  = c("#CEAB07", "deeppink2"),
            legend_items = c("Bombus", "Plant"),
            directory = '../NetworkGraphs/Manuscript',
            ntwrks_page = 3,
            col = 1
)

graphDrawer(Morphology, lbl_fnt = 14,
            fname = 'Method.Morphology',
            edge_clr = 'lightseagreen',
            node_clrs  = c("#CEAB07", "deeppink2"),
            legend_items = c("Bombus", "Plant"),
            directory = '../NetworkGraphs/Manuscript',
            ntwrks_page = 3,
            col = 1
)

graphDrawer(Molecular, lbl_fnt = 14,
            fname = 'Method.Molecular',
            edge_clr = 'lightseagreen',
            node_clrs  = c("#CEAB07", "deeppink2"),
            legend_items = c("Bombus", "Plant"),
            directory = '../NetworkGraphs/Manuscript',
            ntwrks_page = 3,
            col = 1
)


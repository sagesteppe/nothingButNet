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
  paste0(p2d, '/', files[grep('queen_observations', files)]) ) %>% 
  filter(caste == 'q') %>% 
  mutate(species = gsub('\\..*$', '', species)) %>% 
  select(site, doy, week, tot.obs.length, species, plant.species,
         resource.coll:pollen.collected, fl.switch.1:fl.switch.2.resource.coll) %>% 
  mutate(across(resource.coll:fl.switch.2.resource.coll, ~ na_if(.x, ''))) %>% 
  pivot_longer(c(plant.species, fl.switch.1, fl.switch.2), names_to = 'observation', values_to = 'plant.species') %>% 
  drop_na(plant.species) %>% 
  group_by(site, week, species, plant.species) %>% 
  count(name = 'interactions') %>% 
  left_join(., observation_times) %>% 
  mutate(period = case_when(
    week %in% 3:5 ~ 'Early',
    week %in% 6:8 ~ 'Mid',
    week %in% 9:11 ~ 'Late'
  ))

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

early <- resin[['Early']]
mid <- resin[['Mid']]
late <- resin[['Late']]

early <- early[rowSums(early) > 0, colSums(early) > 0]
mid <- mid[rowSums(mid) > 0, colSums(mid) > 0]
late <- late[rowSums(late) > 0, colSums(late) > 0]

rm(bee_obs_wk, resin, bee_obs)

# networks

net = graph_from_incidence_matrix(test, weight=T)
deg = centr_degree(net, mode="all")


colrs <- c("#CEAB07", "deeppink2") # magenta2, goldenrod3
V(net)$color = rbind(c(
  rep(colrs[1],nrow(test)),
  rep(colrs[2],ncol(test))
  )
)

V(net)$size = 5*sqrt(deg$res) # make the more abundant taxa larger
E(net)$width = E(net)$weight/3 # reduce some edge sizes.
# V(i_net)$label <- NA 

template <- layout_in_circle(net)

template_test <- nrow(template)

plot(net, layout=template, 
#            edge.color = "lightseagreen",
 #   vertex.label.dist= lab_dist_vals, 
#   vertex.label.degree = degree_shift, 
    label.font = 3)
legend(x=-2.25, y=1.25, c("Bombus", "Plant"), pch=21, col="#777777", 
       pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)


graphDrawer <- function(x, edge_clr, node_clrs = c(), bg_clr, label_fnt_size, legend_items){
  
  i_net = graph_from_incidence_matrix(x, weight=T)
  deg = centr_degree(i_net, mode="all")
  
  V(i_net)$color = rbind(c(
    rep(node_clrs[1],nrow(early)),
    rep(node_clrs[2],ncol(early))
  )
  )
}


template <- layout_in_circle(net)
template_test <- nrow(template)

plot(net, layout=template, 
     edge.color = "lightseagreen",
     vertex.label.dist= lab_dist_vals, 
     vertex.label.degree = degree_shift, 
     label.font = 3)
legend(x=-2.25, y=1.25, c("Bombus", "Plant"), pch=21, col="#777777", 
       pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)



blanker <- function(network){
  
  netL <- nrow(network) + ncol(network) 
  QnetL <- round(netL / 4 , 0 ) 

  
  if(
    nrow(network) < QnetL
    )
    
    p <- data.frame(
      network[,1:floor(QnetL)], 
      matrix(data = 0, nrow(network), ncol = 2),
      network[,(floor(QnetL) + 1):floor(ncol(network)/2)],
      matrix(data = 0, nrow(network), ncol = 2),
      network[,(floor(ncol(network)/2) + 1):ncol(network)]
      ) 
  
  else 
    
    A = network[1:floor(QnetL),] 
    B = matrix(data = 0, ncol = ncol(network), nrow = 2)
    C = network[(floor(QnetL) + 1):nrow(network),] 
    p <- rbind(A,B,C)
    
    p <- data.frame(
      p[,1:(floor(ncol(p)/2) + 2)],
      matrix(data = 0, nrow(p), ncol = 2),
      p[,(floor(ncol(p)/2) + 5):ncol(p)]
    ) 
    
  names(p) <- gsub('\\.\\.', '\\.', names(p))
  return(p)
}
test <- blanker(mid)

vertex_names <- function(x){
  
  #' After inserting the blank 'observations' into the graph to act as 
  #' vertices we will remove their row and column names to hide the 
  #' dummy names R gives them by default.
  
  vertice_names <- c(rownames(x), colnames(x))
  vertice_names <- sub('^X$', '', vertice_names)
  vertice_names <- sub('^X1$', '', vertice_names)
  vertice_names <- sub('^X2$', '', vertice_names)
  vertice_names <- sub('^X.1$', '', vertice_names)
  vertice_names <- sub('^X.2$', '', vertice_names)
}
VNames <- vertex_names(test)

vert_lab_position <- function(x){
  
  #' Takes a single input, a network, this will serve to specify the positions
  #' for the labels of objects in it. It does not modify the label values, just #
  #' positions and orientation. NOTE: some extra terms are left in the function
  #' for modifications of the upper and lower positions
  
  nodes <- nrow(x) + ncol(x)
  degrees <- 360/nodes
  positions <- degrees * 1:nodes
  
  
  degree_shift <- c(
    rep(0, length(which((positions < 57.5) == TRUE))),
    
    rep(0, length(which((positions > 57.5 & positions < 78.75) == TRUE))),
    rep(-pi/2, length(which((positions > 78.75 & positions < 101.25) == TRUE))),
    rep(-pi/2, length(which((positions > 101.25 & positions < 112.5) == TRUE))),
    
    rep(pi, length(which((positions > 112.5 & positions < 247.5) == TRUE))),
    
    rep(pi/2, length(which((positions > 247.5 & positions < 278.75) == TRUE))),
    rep(pi/2, length(which((positions > 278.75 & positions < 301.25) == TRUE))),
    rep(0, length(which((positions > 301.25 & positions < 312.5) == TRUE))), 
    
    rep(0, length(which((positions > 312.5) == TRUE)))
    # '0' right, '-pi/2' up, 'pi' left, 'pi/2' down
  )
  
  
  dist_vals <- c(
    rep(8, length(which((positions < 57.5) == TRUE))),
    
    rep(6, length(which((positions > 57.5 & positions < 78.75) == TRUE))),
    rep(c(2,3), each = 1, length.out = 
          length(which((positions > 78.75 & positions < 101.25) == TRUE))),
    rep(2, length(which((positions > 101.25 & positions < 112.5) == TRUE))),
    
    rep(8, length(which((positions > 112.5 & positions < 247.5) == TRUE))),
    
    rep(2, length(which((positions > 247.5 & positions < 278.75) == TRUE))),
    rep(c(3,2), each = 1, length.out = 
          length(which((positions > 278.75 & positions < 301.25) == TRUE))),
    rep(6, length(which((positions > 301.25 & positions < 312.5) == TRUE))), 
    
    rep(8, length(which((positions > 312.5) == TRUE)))
    # '0' right, '-pi/2' up, 'pi' left, 'pi/2' down
  )
  
  labels <- list('degree_shift' = degree_shift, 
                 'dist_vals' = dist_vals)
  
  return(labels)
}
VLPs <- vert_lab_position(test)

net = graph_from_incidence_matrix(test, weight=T)
deg = centr_degree(net, mode="all")

set_colors <- function(x, node_clr1, node_clr2, bg_clr){
  
  #' This serves to set the colors of up to three groups in the plot
  #' the first two groups being those major groups which interact in the 
  #' network. the final group is the background color of the plot. This will
  #' default to a white if left blank. 
  
  if(missing(bg_clr)){
    bg_clr <- '#FFFFFF'
  } # doesnt really change the graph... but leave anyways.
  
  my_cols = rbind(c(
    rep(node_clr1,nrow(test)),
    rep(node_clr2,ncol(test)))
  )
  
  my_cols[which(VNames == "")] <- bg_clr 
  V(x)$color <- my_cols
  return(x)
  
}

net <- set_colors(net, node_clr1 = "#CEAB07", node_clr2 = "deeppink2")


V(net)$size = 5*sqrt(deg$res) 
E(net)$width = E(net)$weight/3 
template <- layout_in_circle(net)

plot(net, layout=template, 
    edge.color = "lightseagreen",
    vertex.label = VNames,
    vertex.label.dist= VLPs$dist_vals, 
    vertex.label.degree = VLPs$degree_shift, 
    label.font = 3)


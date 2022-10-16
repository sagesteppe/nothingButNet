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
  
  "M.ciliata" ,"M.fusiformis",
  
  "A.lewisii",
  "C.lanceolata",
  "D.involucrata",
  "S.sp",
  "V.praemorsa", # ???/
  
  "H.capitatum", "H.fendleri",
  "A.urticifolia", "P.procera", "P.bracteosa",
  
  "C.parryi", "D.hoopesii", "H.quinquenervis", "L.bigelovii",
  "E.speciosus", "S.integerrimus", "T.officinale", "W.arizonica",
  
  "F.speciosa",
  "V.occidentalis",
  "D.pulchellum",
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

arrange_nets <- function(x, col_arrange, row_arrange){
  
  #' This function serves to help standardize the order of species in these networks.
  #' We use two different schema to arrange our taxa in one net! The insects we arrange
  #' based on functional parameters, in this case proboscis length. The plants we 
  #' arrange using a phylogenetic and alphabetical approach. These arguments must be 
  #' given as character vectors with names in the EXACT SAME format as the input data.
  #' 
  #' INPUTS: 
  #' x = a network in split format
  #' col_arrange = a character vector with it's elements in the order which you want the columns arranged (all column names must be in the vector, and be in the format 'GENUS.EPITHET' (note no spaces))
  #' row_arrange = a character vector with it's elements in the order which you want the rows arranged (all row names must be in the vector)
  #' 
  if(missing(row_arrange)) { row_arrange <- rownames(x) }
  
  x <- x[rowSums(x) > 0, colSums(x) > 0]
  colnames(x) <- sub(" ", "", colnames(x))
  col_arrange <- sub(" ", "", col_arrange)
  row_arrange <- sub(" ", "", row_arrange)
  
  col2grab <- match(col_arrange, colnames(x)) |> na.omit()
  x <- x[, col2grab]
  
  row2grab <- match(row_arrange, rownames(x)) |> na.omit()
  x <- x[row2grab, ]
  
  return(x)
  
}
blanker <- function(x){
  
  netL <- nrow(x) + ncol(x) 
  QnetL <- round(netL / 4 , 0 ) 
  
  if(
    nrow(x) == ncol(x)) {
    
    A = x[1:floor(QnetL),] 
    B = matrix(data = 0, ncol = ncol(x), nrow = 3)
    C = x[(floor(QnetL) + 1):nrow(x),] 
    p <- rbind(A,B,C)
    
    p <- data.frame(
      p[,1:(floor(QnetL))],
      matrix(data = 0, nrow(p), ncol = 3),
      p[,(floor(QnetL)+1):ncol(x)]
    )
    
  } else if(
    nrow(x) < QnetL
  ) {
    
    p <- data.frame(
      x[,1:floor(QnetL)], 
      matrix(data = 0, nrow(x), ncol = 3),
      x[,(floor(QnetL) + 1):floor(ncol(x)/2)],
      matrix(data = 0, nrow(x), ncol = 3),
      x[,(floor(ncol(x)/2) + 1):ncol(x)]
    ) 
    
  } else {
    
    A = x[1:floor(QnetL),] 
    B = matrix(data = 0, ncol = ncol(x), nrow = 3)
    C = x[(floor(QnetL) + 1):nrow(x),] 
    p <- rbind(A,B,C)
    
    p <- data.frame(
      p[,1:(floor(ncol(p)/2) + 2)],
      matrix(data = 0, nrow(p), ncol = 3),
      p[,(floor(ncol(p)/2) + 5):ncol(p)]
    )
    
  }
  
  names(p) <- gsub('\\.\\.', '\\.', names(p))
  return(p)
}
vertex_names <- function(x){
  
  #' After inserting the blank 'observations' into the graph to act as 
  #' vertices we will remove their row and column names to hide the 
  #' dummy names R gives them by default.
  
  vertice_names <- c(rownames(x), colnames(x))
  vertice_names <- sub('^X$', '', vertice_names)
  vertice_names <- sub('^X1$', '', vertice_names)
  vertice_names <- sub('^X3$', '', vertice_names)
  vertice_names <- sub('^X2$', '', vertice_names)
  vertice_names <- sub('^X.1$', '', vertice_names)
  vertice_names <- sub('^X.2$', '', vertice_names)
  vertice_names <- sub('^X.3$', '', vertice_names)
  
  vertice_names <- gsub("\\.(?=[A-Za-z])", ". ", vertice_names, perl = TRUE)
}
vert_lab_position <- function(x){
  
  #' Takes a single input, a network, this will serve to specify the positions
  #' for the labels of objects in it. It does not modify the label values, just
  #' positions and orientation. NOTE: some extra terms are left in the function
  #' for modifications of the upper and lower positions
  
  nodes <- nrow(x) + ncol(x)
  degrees <- 360/nodes
  positions <- degrees * 1:nodes
  
  degree_shift <- c(
    rep(0, length(which((positions < 65) == TRUE))),
    rep(-pi/2, length(which((positions > 65 & positions < 125) == TRUE))),
    
    rep(pi, length(which((positions > 125 & positions < 245) == TRUE))),
    
    rep(pi/2, length(which((positions > 245 & positions < 295) == TRUE))),
    rep(0, length(which((positions > 295) == TRUE)))
    # '0' right, '-pi/2' up, 'pi' left, 'pi/2' down
  )
  
  
  dist_vals <- c(
    rep(4, length(which((positions < 65) == TRUE))),
    rep(c(2,3), each = 1, length.out = 
          length(which((positions > 65 & positions < 125) == TRUE))),
    
    rep(4, length(which((positions > 125 & positions < 245) == TRUE))),
    
    rep(c(2,3), each = 1, length.out = 
          length(which((positions > 245 & positions < 295) == TRUE))),
    rep(4, length(which((positions > 295) == TRUE)))
    
    # '0' right, '-pi/2' up, 'pi' left, 'pi/2' down
  )
  
  labels <- list('degree_shift' = degree_shift, 
                 'dist_vals' = dist_vals)
  
  return(labels)
}
set_colors <- function(x, net, node_clrs,  bg_clr, VNames){
  
  #' This serves to set the colors of up to three groups in the plot
  #' the first two groups being those major groups which interact in the 
  #' network. the final group is the background color of the plot. This will
  #' default to a white if left blank. 
  
  if(missing(bg_clr)){
    bg_clr <- '#FFFFFF'
  } # doesn't really change the graph... but leave anyways.
  
  my_cols = rbind(c(
    rep(node_clrs[1],nrow(x)),
    rep(node_clrs[2],ncol(x)))
  )
  
  my_cols[which(VNames == "")] <- bg_clr 
  V(net)$color <- my_cols
  return(net)
  
}
graph_dims <- function(ntwrks_page, col){
  
  #' this function hopefully returns square graphs of similar size which may be
  #' somewhat readily assembled onto an A4 page with what i hope are standard
  #' page margins (2.97cm top & bottom, 2.1 cm sides). It is easily modifiable for
  #' alternative dimensions. It assumes you want to fill an entire page with graphs.
  
  W = round(1984 / col, 0)
  H = round(2864 / (ntwrks_page/col), 0 )
  if(W > H){ W <- H} else {H <- W}
  
  dims <- list('W' = W, 
               'H' = H)
  return(dims)
}

graphDrawer <- function(data, plot_name, edge_clr, node_clrs,  bg_clr, 
                        lbl_fnt, legend_items, directory, fname,
                        ntwrks_page, col, H, W){
  
  #' this function serves to draw multiples of similar graphs, it's most logical 
  #' applications are in displaying information across temporal slices, or 
  #' spatial replicates or treatments thereof.
  #' 
  #' INPUTS
  #' data = a list of dataframes, wherein one group of interacting organisms (e.g.)
  #' insects are along one axis, like rows, and the other group of organisms is 
  #' along the other axis. This will be coerced to an igraph object. Named lists may 
  #' serve as titles of the rendered plots. 
  #' edge_clr, a color for mapping the edges (interactions) between the nodes (species)
  #' node_clrs, a set of two colors for your major groups
  #' bg_clr, a background color for the network, defaults to white
  #' label_fnt_size, label font size defaults to 9
  #' legend item, the name of the two groups on the axis, e.g. c('Insects', 'Plants')
  #' directory = name of directory to save graphs to, defaults to 'NetworkGraphs'
  #' fname = filename for the output graph, defaults to name of input
  #' ntwrks_page = (numeric) how many networks per page? Defaults to one
  #' col = (numeric) how many columns of networks? Defaults to one
  #' H = (numeric) heigh of graph in pixels. If left blank is calculated to maximixe 
  #' graph sizes on an A4 page with standard margins based on 'ntwrks_page' and 'col'
  #' W = (numeric) length of graph in pixels. If left blank is calculated to maximixe 
  #' graph sizes on an A4 page with standard margins based on 'ntwrks_page' and 'col'
  
  if(missing(plot_name)) { plot_name <- substitute(data) }
  if(missing(lbl_fnt)) { lbl_fnt <- 14 }
  if(missing(directory)) { directory <- 'NetworkGraphs' }
  if(missing(fname)) { fname <- paste0(substitute(data), '.png') 
  } else {fname <- paste0(fname, '.png')}
  
  if(missing(ntwrks_page)) { ntwrks_page <- 1 } | if(missing(col)) { col <- 1 }
  if(missing(H)) { H <- NA}  |  if(missing(W)) { W <- NA}
  dims <- list('H' = H, 'W' = W)
  if(is.na(dims$H) + is.na(dims$W) == 0){
    dims 
  } else if(is.na(dims$H) & !is.na(dims$W)){
    dims$H <- dims$W
  } else if(is.na(dims$W) & !is.na(dims$H)){
    dims$W <- dims$H
  } else {
    dims <- graph_dims(ntwrks_page, col)
  }
  
  filename <- file.path(directory, fname)
  
  ifelse(!dir.exists(file.path(directory)),
         dir.create(file.path(directory)), FALSE)
  
  blanked_data <- blanker(data)
  VNames <- vertex_names(blanked_data)
  VLPs <- vert_lab_position(blanked_data)
  
  net = graph_from_incidence_matrix(blanked_data, weight=T)
  deg = centr_degree(net, mode="all")
  
  net <- set_colors(x = blanked_data, net, node_clrs, VNames = VNames)
  
  V(net)$label.color <- 'black'
  V(net)$size = 5*sqrt(deg$res) 
  E(net)$width = E(net)$weight/4
  template <- layout_in_circle(net)
  
  png(filename,
      width = dims$W, height = dims$H, units = "px", pointsize = 12)
  
  par(mar=c(5,6,5,6))
  plot(net, layout=template, 
       edge.color = edge_clr,
       vertex.label = VNames,
       vertex.label.dist= VLPs$dist_vals, vertex.label.degree = VLPs$degree_shift, 
       label.font = lbl_fnt,
       main = plot_name
       ) 
    legend(x= -0.25, y=-1.35, legend_items, 
           pch=21, col="#777777", 
           pt.bg=node_clrs, 
           pt.cex=2, cex=.8, bty="n", ncol=1)
    
 invisible(dev.off())
 
 message(paste0("'", plot_name, "' has been rendered as a graph and saved to:\n ", filename))

}

# networks

tet <- lapply(resin, arrange_nets, col_arrange = arranged_plants, 
              row_arrange = arranged_bees)
rm(bee_obs_wk, resin, bee_obs)

early <- tet[['Early']]
mid <- tet[['Mid']]
late <- tet[['Late']]

graphDrawer(mid, lbl_fnt = 14,
            plot_name = 'mid',
            edge_clr = 'lightseagreen',
            node_clrs  = c("#CEAB07", "deeppink2"),
            legend_items = c("Bombus", "Plant"),
            fname = 'floral-mid', 
            ntwrks_page = 12,
            col = 3
)

# the below looks cool...
# map(out.list, ~ network.fun(nodes = .x$nodes, edges = .x$edges))


map2(out.list, names(out.list), ~ network.fun(plotName = .y, nodes = .x$nodes, edges = .x$edges))


names(late)

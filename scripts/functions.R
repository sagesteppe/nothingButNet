#' Rearrange elements of a network thematically
#' 
#' 'arrange_nets()' helps standardize the order of species in these networks.
#' We use two different schema to arrange our taxa in one net! The insects we arrange
#' based on functional parameters, in this case proboscis length. The plants we 
#' arrange using a phylogenetic and alphabetical approach. These arguments must be 
#' given as character vectors with names in the EXACT SAME format as the input data.
#' 
#' @param x = a network in split format
#' @param col_arrange = a character vector with it's elements in the order which you want the columns arranged (all column names must be in the vector, and be in the format 'GENUS.EPITHET' (note no spaces))
#' @param row_arrange = a character vector with it's elements in the order which you want the rows arranged (all row names must be in the vector)
#' @export
arrange_nets <- function(x, col_arrange, row_arrange){
  
  if(missing(row_arrange)) { row_arrange <- rownames(x) }
  
  x <- x[rowSums(x) > 0, colSums(x) > 0, drop = F]
  colnames(x) <- sub(" ", "", colnames(x))
  col_arrange <- sub(" ", "", col_arrange)
  row_arrange <- sub(" ", "", row_arrange)
  
  col2grab <- match(col_arrange, colnames(x)) |> na.omit()
  x <- x[, col2grab, drop = F]
  
  row2grab <- match(row_arrange, rownames(x)) |> na.omit()
  x <- x[row2grab, , drop = F]
  
  return(x)
  
}


#' cleanup column and row names from blanker
#' Create information rich circular networks
#' 
#' 'blanker()' inserts several 'empty' observations into the rows and/or columns of
#' a matrix of interaction data so that nodes can be labelled at the top and bottom
#' of an igraph circular graph.
#' 
#' @param x = a network in matrix format
#' @returns = a network in matrix format with blank rows/columns inserted
#' @seealso 'vertex_names', 'graphDrawer'
#' @export
blanker <- function(x){
  
  netL <- nrow(x) + ncol(x) 
  QnetL <- round(netL / 4 , 0 ) 
  colN <- ceiling(netL/10)
  if(colN < 3){colN <- 3}
  thirdStart <- round( round(netL * 3/4) - (colN*3/4) ) 
  
  if( # net work with equal number of members in each bipartite group
    nrow(x) == ncol(x)) {
    
    A = x[1:round(QnetL),] 
    B = matrix(data = 0, ncol = ncol(x), nrow = 3)
    C = x[(round(QnetL) + 1):nrow(x),] 
    p <- rbind(A,B,C)
    
    p <- data.frame(
      p[,1:(round(QnetL))],
      matrix(data = 0, nrow(p), ncol = 3),
      p[,(round(QnetL)+1):ncol(x)]
      
    )
    
  } else if( # one bipartite group has many more members than the other
    nrow(x) < QnetL
  ) {
    
    sides <- netL/2
    if(netL %% 2 == 1){ # variables for odd numbers of members in a network
      
      Lside <- ceiling(sides) - ceiling(sides) %% 2 ; Rside <- netL - Lside
      RSUpper <- ceiling(Rside/2) - ceiling(Rside/2) %% 2
      RSLower <- Rside - RSUpper 
      
      
    }  else { # variables for even numbers of members in a network
      
      Lside <- sides; RSUpper <- ceiling(sides/2);  RSLower <- floor(sides/2)
      
    }
    
    p <- data.frame( # make all plots for networks with many more of one bipartite member than the other
      x[,1:(RSUpper - nrow(x))], 
      matrix(data = 0, nrow(x), ncol = colN),
      x[,(RSUpper - nrow(x) + 1):(ncol(x) - (RSLower))],
      matrix(data = 0, nrow(x), ncol = colN),
      x[, (ncol(x) - RSLower + 1):ncol(x)] )

  } else { # one group has only a few more members than the other
    
    pA <- rbind(
      x[1:floor((netL / 4)),],
      matrix(data = 0, nrow = colN, ncol = ncol(x)),
      x[ceiling((netL / 4)):nrow(x), , drop=FALSE] )
    
    r <- (netL / 2)  + (netL / 4) -  nrow(x)
    
    p <- data.frame(
      pA[,1:floor(r)], 
      matrix(data = 0, nrow(pA), ncol = colN),
      pA[,(floor(r) + 1):ncol(x)] )
    
  }
  
  colnames(p) <- gsub('\\.\\.', '\\.', colnames(p))
  rownames(p) <- gsub('X[.]', 'X', rownames(p))
  return(p)
}


#' Clean up element names from 'blanker'
#' 
#' After inserting the blank 'observations' into the graph to act as 
#' vertices we will remove their row and column names to hide the 
#' dummy names R gives them by default.
#' @param x = a matrix output from 'blanker'
#' @seealso 'blanker', 'graphDrawer'
#' @export
vertex_names <- function(x){
  
  #' After inserting the blank 'observations' into the graph to act as 
  #' vertices we will remove their row and column names to hide the 
  #' dummy names R gives them by default.
  
  vertice_names <- c(rownames(x), colnames(x))
  vertice_names <- sub('^X$', '', vertice_names)
  vertice_names <- sub('^X[0-9]$', '', vertice_names)
  vertice_names <- sub('^X.[0-9]$', '', vertice_names)
  vertice_names <- sub('^X[0-9].[0-9]$', '', vertice_names)
  
  vertice_names <- gsub("\\.(?=[A-Za-z])", ". ", vertice_names, perl = TRUE)
}


#' Give node labels rendering space
#' 
#' 'vert_lab_position' specifies the positions for the labels of objects in the 
#' graphs. It does not modify the label values, just positions and orientation. 
#' 
#' @param x = a matrix output from vertex_names
#' @col = the number of columns of graphs
#' @seealso 'vertex_names', 'graphDrawer'
#' @export
vert_lab_position <- function(x, col){
  
  m_dist <- if(col == 3){m_dist = 4.5} else if(col < 3){m_dist = 3} 
  
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
    rep(m_dist, length(which((positions < 65) == TRUE))),
    seq(from = 2, to = 7, length.out = 
          length(which((positions > 65 & positions < 90) == TRUE))),
    
    seq(from = 7, to = 2, length.out = 
          length(which((positions > 90 & positions < 125) == TRUE))),
    
    rep(m_dist, length(which((positions > 125 & positions < 245) == TRUE))),
    
    seq(from = 2, to = 7, length.out = 
          length(which((positions > 245 & positions < 270) == TRUE))),
    seq(from = 7, to = 2, length.out = 
          length(which((positions > 270 & positions < 295) == TRUE))),
    
    rep(m_dist, length(which((positions > 295) == TRUE)))
    
  )
  
  labels <- list('degree_shift' = degree_shift, 
                 'dist_vals' = dist_vals)
  
  return(labels)
}

#' Custom colouring for the net
#' 
#' This serves to set the colors of up to three groups in the plot
#' the first two groups being those major groups which interact in the 
#' network. the final group is the background color of the plot. This will
#' default to a white if left blank.
#' @param x = the matrix output from 'blanker'
#' @param net = network from 'igraph::graph_from_incidence_matrix()'
#' @param node_clrs = a vector of length two of hex or named colors for the major groups
#' @param bg_clr = background color for the plot
#' @param VNames = a vector of names from 'vertex_names' function
#' @seealso = 'vertex_names', 'blanker', 'graphDrawer'
#' @export
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

#' Size graphs for assembling on an A4 page
#' 
#' this function hopefully returns square graphs of similar size which may be
#' somewhat readily assembled onto an A4 page with what i hope are standard
#' page margins (2.97cm top & bottom, 2.1 cm sides). It is easily modifiable for
#' alternative dimensions. It assumes you want to fill an entire page with graphs.
#' @param ntwrks_page = the number of networks to place on the page
#' @col  the number of columns to spread the networks across
#' @seealso  'graphDrawer', 'legenDrawer'
#' @export
graph_dims <- function(ntwrks_page, col){
  
  W = round(1984 / col, 0)
  H = round(2864 / (ntwrks_page/col), 0 )
  if(W > H){ W <- H} else {H <- W}
  
  dims <- list('W' = W, 
               'H' = H)
  return(dims)
  
}


#' Render multiple graphs of similar thematic content
#' 
#' 'graph_dims' draw multiples of similar graphs, some applications are in 
#' displaying information across temporal slices, or spatial replicates or 
#' treatments thereof. It calls most of the functions in the package.
#' 
#' @param data = a list of matrcies, wherein one group of interacting organisms 
#' (e.g.) insects are along one axis, like rows, and the other group of organisms
#'  is along the other axis. This will be coerced to an igraph object. Named 
#'  lists may serve as titles of the rendered plots. 
#' @param edge_clr, a color for mapping the edges (i.e. interactions) between 
#' the nodes (i.e. species)
#' @param node_clrs, a set of two colors for your major groups
#' @param bg_clr, a background color for the network, defaults to white
#' @param lbl_fnt, label font size defaults to 14
#' @param legend item, the name of the two groups on the axis, 
#' e.g. c('Insects', 'Plants')
#' @param directory = directory to save graphs to, defaults to 'NetworkGraphs'
#' @param fname = filename for the output graph, defaults to name of input ('data'),
#' but see 'netPage' for naming conventions that allow quick assembly of plots
#' onto a single page. 
#' @param ntwrks_page = (numeric) how many networks per page? Defaults to 1
#' @param col = (numeric) how many columns of networks? Defaults to 1
#' @param H = (numeric) height of graph in pixels. If left blank is calculated 
#' to maximixe graph sizes on an A4 page with standard margins based on 
#' 'ntwrks_page' and 'col'
#'  and 'col'
#' @param W = see 'H'
#' @return a graph to a location on disk as a png 
#' @seealso 'blanker', 'vertex_names, vert_lab_position','set_colors',
#' 'graph_dims'
#' @export
graphDrawer <- function(data, edge_clr, node_clrs,  bg_clr, 
                        lbl_fnt, legend_items, directory, fname,
                        ntwrks_page, col, H, W){
  
  if(missing(lbl_fnt)) { lbl_fnt <- 14 }
  if(missing(directory)) { directory <- 'NetworkGraphs/Intermediates' }
  if(missing(fname)) { 
    fname <- paste0(substitute(data), '.png')} else{
      fname <- paste0(fname, '.png')}
  if(missing(ntwrks_page)) { ntwrks_page <- 1 } | if(missing(col)) { col <- 1 }
  if(missing(H)) { H <- NA}  |  if(missing(W)) { W <- NA}
  if(missing(edge_clr)) {edge_clr <- 'black'}
  
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
  
  filename <- file.path(directory, paste0(fname))
  ifelse(!dir.exists(file.path(directory)),
         dir.create(file.path(directory)), FALSE)
  
  blanked_data <- blanker(data)
  VNames <- vertex_names(blanked_data)
  VLPs <- vert_lab_position(blanked_data, col)
  
  net = igraph::graph_from_incidence_matrix(blanked_data, weight=T)
  deg = igraph::centr_degree(net, mode="all")
  
  net <- set_colors(x = blanked_data, net, node_clrs, VNames = VNames)
  
  V(net)$label.color <- 'black'
  V(net)$size = 2.5*sqrt(deg$res)
  E(net)$width = E(net)$weight/2
  
  template <- layout_in_circle(net)
  
  png(filename,
      width = dims$W, height = dims$H, units = "px", pointsize = 12)
  
  par(mar=c(3.5,6.5,3.5,6.5))
  plot(net, layout=template, 
       edge.color = edge_clr,
       vertex.label.cex = 1.25,
       vertex.label = VNames,
       vertex.label.dist= VLPs$dist_vals, 
       vertex.label.degree = VLPs$degree_shift, 
       label.font = lbl_fnt
       ) 
    
 invisible(dev.off())
 
 message(paste0("'", deparse(substitute(data)), 
                "' has been rendered as a graph and saved to:\n ", filename))
 
}


#' Draw three legends for the plots and save to png
#' 
#' This quickly draws three different legends. One is a table to show the 
#' abbreviations which are used in the nets. One legend shows the two major
#' categorical colors which are used in the nets. The final legend is a size legend, 
#' which is not to scale, but numerically shows the range in the number of 
#' interacting species. 
#' 
#' @param x a list of matrices to render graphs of (output of 'arrange_nets')
#' @param table_title a character, table name
#' @param table_items vector of names for legend portion of the table
#' @param directory location to size legend before assembling to final product
#' @param fname a file name for the legend, defaults to legend. 
#' @param legend_items character vector with names of node items
#' @param node_clrs character vector with node item colors
#' @param fill_col a character vector of fill colors for the size bubbles
#' @param y.space a numeric vector of 'distances' for the legend elements
#' @param filename defaults to 'SizeLegend'
#' @param ntwrks_page the number of networks you plan on placing on each page
#' @param LcolN number of columns to split legend across
#' @param LegcolN number of columns in the table
#' @param colN number of columns for the nets on the page
#' @seealso 
#' @example category_legend_drawer(node_clrs  
#' = c("#CEAB07", "deeppink2"), LcolN = 1,
#'  legend_items = c("Bombus", "Plant"), ntwrks_page = 9, colN =3)
#' @export
tableLegend <- function(x, table_title, table_items, directory, fname, legend_items, node_clrs,
                        fill_col, y.space, ntwrks_page, colN, LcolN, LegcolN){
  
  if(missing(directory)) { directory <- 'NetworkGraphs/Intermediates' }
  if(missing(fname)) {fname <- 'TableLegend.png'}
  if(missing(table_title)){table_title <- ""}
  if(missing(fill_col)){fill_col <- 'white'}
  if(missing(y.space)){y.space <- seq(from = 1, to = 2, length.out = 5)}
  
  dims <- graph_dims(ntwrks_page, col = colN)
  
  # create the categorical legend.
  
  d <- data.frame(Group = c('Bombus', 'Plant'),
                  Dummy = c(1,1), Dummy2 = c(1,1))
  
  v <- c("#CEAB07", "deeppink2")
  names(v) <- c('Bombus', 'Plant')
  
  p <- ggplot2::ggplot(d, ggplot2::aes(x= Dummy, y=Dummy2, color = Group)) + 
    ggplot2::geom_point() +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size=12))) + 
    ggplot2::scale_color_manual(values = v) + 
    ggplot2::labs(color = 'Major Interacting\nGroups ') +
    ggplot2::theme(legend.key = ggplot2::element_rect(fill = "white"),
                   legend.title.align = 0.5,
                   legend.key.size = unit(2.5, 'cm'), 
                   legend.key.height = unit(3, 'cm'), 
                   legend.key.width = unit(3, 'cm'), 
                   legend.title = ggplot2::element_text(size = 28),
                   legend.text = ggplot2::element_text(size = 24)) 
  
  catLeg <- cowplot::get_legend(p)
  
  #  create the legend table. 
  lname <- file.path(directory, fname)
  grp <- ceiling(length(table_items)/LegcolN)
  l <- (grp * LegcolN) - length(table_items)
  
  table_items <- c(table_items, rep("", l))
  v <- matrix(data = table_items , nrow = grp, ncol = LegcolN)
  tt2 <- gridExtra::ttheme_minimal(core=list(fg_params=list(hjust=0, x=0.1)),
                                   rowhead=list(fg_params=list(hjust=0, x=0)),
                                   base_size = 20)
  
  tableLeg <- gridExtra::tableGrob(v, theme = tt2)
  tableLeg <- gridExtra::grid.arrange(top = grid::textGrob(
    table_title, gp=grid::gpar(fontsize=40)),
    tableLeg)
  
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
                   legend.title = ggplot2::element_text(size = 28),
                   legend.text = ggplot2::element_text(size = 24)) +
    ggplot2::labs(size = 'No. of Interacting\nspecies')
  
  sizeLeg <- cowplot::get_legend(p)
  
  # now combine the three legends. 
  grobBY <- list(tableLeg, sizeLeg, catLeg)
  
  lays <- rbind(
    c(NA,NA,NA,NA,NA,NA),
    c(1,1,1,1,1, 2),
    c(1,1,1,1,1, 2),
    c(1,1,1,1,1, 2),
    c(1,1,1,1,1, 2),
    c(1,1,1,1,1, 2),
    c(1,1,1,1,1, 3),
    c(1,1,1,1,1, 3),
    c(1,1,1,1,1, 3),
    c(NA,NA,NA,NA,NA,NA)
  )
  
  bucket <- invisible(gridExtra::grid.arrange(grobs = grobBY, layout_matrix = lays))
  png(lname,
      width = 1984, height = dims$H, units = "px", pointsize = 12)
  invisible(grid::grid.draw(bucket))
  invisible(dev.off())
  
  message(paste0("'", fname, 
                 "' has been rendered as a legend and saved to:\n ",
                 file.path(directory, fname)))
}


#' mosaic the grids together using gridExtra
#' 
#' netPage simply uses list.files with cowplot to pull back in your plots for 
#' quick assembly into a grid. It's use is optional. 
#' @param directory the folder holding your plots output from 'graphDrawer'
#' @param col_var a character vector of grouping variable to sort columns from
#'  left to right
#' @param row_var a character vector of a grouping variable to sort rows from 
#' left to right
#' @example 
#' netPage(col_var = c('observation', 'microscopy', 'molecular'), 
#"        row_var = c('early', 'mid', 'late'))
#' @return a cowplot arranged grid.
#' @export
#' 
#' 
#' 
nets2Page <- function(directory, col_var, row_var, fname, sep_char, mainT, Tlegend_fname){
  
  if(missing(directory)) {directory <- 'NetworkGraphs/Intermediates'}
  dirIN <- file.path(directory); dirOUT <- file.path(directory, 'Output')
  if(missing(fname)) {fname <- 'mosaiced_Nets.pdf'}
  if(missing(sep_char))  {sep_char <- ''}
  if(missing(mainT))  {mainT <- ''}
  if(missing(Tlegend_fname)) {Tlegend_fname <- 'TableLegend.png'}
  
  rowOrder <- matrix( # this sets up a matrix to have rows groups in order
    data = rep(row_var, each = length(col_var)),
    ncol = length(col_var), byrow = T)
  colOrder <- matrix( # this sets up a matrix to have column groups in order
    data = rep(col_var, each = length(row_var)),
    ncol = length(col_var), byrow = F)
  image_orders <- matrix(paste0(colOrder, sep_char, rowOrder, '.png'), 
                         ncol = length(col_var))
  image_orders <- as.vector(t(image_orders)) 
  image_orders <- file.path(dirIN, image_orders)  
  images = lapply(image_orders, png::readPNG)
  grob_images = lapply(images, grid::rasterGrob)
  
  colN <- length(col_var)
  nPlots <- (colN * length(row_var))
  rowN <- nPlots / length(col_var)
  
  # define the layout
  mlayout <- matrix(nrow = rowN, ncol = colN, byrow = T,
                   data = c(rep(1:(rowN -1), each = colN),
                            rowN:((rowN-1) + colN) ) )
  mlayout <- rbind(mlayout, rep((max(mlayout) + 1), ncol(mlayout)))  

  # Recover the top grobs!!
  top_grobs <- split(grob_images[1:(length(grob_images) - rowN)],
                     ceiling(seq_along(grob_images[1:(length(grob_images) -
                                                        rowN)]) / rowN))
  names(top_grobs) <- paste0('t', seq(1:length(top_grobs)))
  
  # title up the real top grob
  t1 <- gridExtra::arrangeGrob(grobs = top_grobs$t1, ncol = colN, 
                               padding = unit(0.0, "line"),
                               left = row_var[1], top = grid::textGrob(mainT, vjust = -1.1))

  # this will get all of the subsequent top grobs !!!!!!!!
  
  if(length(top_grobs) >= 3){
  top_grobs <- mapply(gridExtra::arrangeGrob, 
                      grobs = top_grobs[(colN+1):length(top_grobs)], ncol = colN, 
                      left = row_var[(colN+1):(length(row_var)-1)],  
                      padding = unit(0.0, "line"))
  } else { 
    top_grobs <- gridExtra::arrangeGrob(
                        grobs = top_grobs[2], ncol = colN, 
                        left = row_var[2],  
                        padding = unit(0.0, "line") )
    }
    
  list2env(top_grobs, env = environment())
  
  # this recovers the bottom grobs
  L <- length(grob_images)
  range <- ((L - rowN) + 1 ):L
  bottom_grobs <- split(grob_images[((L - rowN) + 1 ): L],
                        seq_along(range) / 1)
  names(bottom_grobs) <- paste0('b', seq(1:length(bottom_grobs)))
  
  if(colN > 2){
  b1 <- gridExtra::arrangeGrob(grobs = bottom_grobs$b1, bottom = col_var[1],
                               left = row_var[length(row_var)],
                               padding = unit(0.0, "line")
                         )
  
  bottoms <- mapply(gridExtra::arrangeGrob, 
                    grobs = bottom_grobs[2:length(bottom_grobs)],
                    bottom = col_var[1:(length(col_var))],
                    padding = unit(0.0, "line")
  )
  } else {
    b1 <- gridExtra::arrangeGrob(grobs = bottom_grobs$b1, bottom = col_var[1],
                                 left = row_var[length(row_var)],
                                 padding = unit(0.0, "line") )
    
    bottoms <- mapply(gridExtra::arrangeGrob, 
                    grobs = bottom_grobs[2],
                    bottom = col_var[2],
                    padding = unit(0.0, "line") )
  }
  
  list2env(bottoms, env = environment())

  # ensure the grobs are in the correct order for the mosaic. 
  g2p <- mget(ls(pattern = '[t|b][1-9]{1}'))
  groborder <- c(paste0('t', seq(1:9)), paste0('b', seq(1:9)))
  ord2grab <- match(groborder, names(g2p)) |> na.omit()
  g2p <- g2p[ord2grab]
  
  # load and place the legend onto the grobs2plot
  legend <- grid::rasterGrob(png::readPNG(file.path(dirIN, Tlegend_fname)))
  legend <- gridExtra::arrangeGrob(legend, nrow = 1)
  g2p <- c(g2p, leg = list(legend)) 
  
  # place on the page and print.
  
  ifelse(!dir.exists(file.path(dirOUT)),
         dir.create(file.path(dirOUT)), FALSE)
 
  return(groborder)
  ml <- gridExtra::marrangeGrob(grobs = g2p, 
                                layout_matrix = mlayout, top = "")
  pdf(file = file.path(dirOUT, fname), paper = 'a4')
  print(ml)
  invisible(dev.off())

  message(paste0("'", fname, 
                 "' has been rendered as a pdf and saved to:\n ",
                 file.path(dirOUT, fname)))
  
}


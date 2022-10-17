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
  vertice_names <- sub('^X1$', '', vertice_names)
  vertice_names <- sub('^X3$', '', vertice_names)
  vertice_names <- sub('^X2$', '', vertice_names)
  vertice_names <- sub('^X.1$', '', vertice_names)
  vertice_names <- sub('^X.2$', '', vertice_names)
  vertice_names <- sub('^X.3$', '', vertice_names)
  
  vertice_names <- gsub("\\.(?=[A-Za-z])", ". ", vertice_names, perl = TRUE)
}


#' Give node labels rendering space
#' 
#' 'vert_lab_position' specifies the positions for the labels of objects in the 
#' graphs. It does not modify the label values, just positions and orientation. 
#' 
#' @param x = a matrix output from vertex_names
#' @seealso 'vertex_names', 'graphDrawer'
#' @export
vert_lab_position <- function(x){
  
  
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
    rep(4.5, length(which((positions < 65) == TRUE))),
    rep(c(2.5,3.5), each = 1, length.out = 
          length(which((positions > 65 & positions < 125) == TRUE))),
    
    rep(4.5, length(which((positions > 125 & positions < 245) == TRUE))),
    
    rep(c(2.5,3.5), each = 1, length.out = 
          length(which((positions > 245 & positions < 295) == TRUE))),
    rep(4.5, length(which((positions > 295) == TRUE)))
    
    # '0' right, '-pi/2' up, 'pi' left, 'pi/2' down
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
#' @col = the number of columns to spread the networks across
#' @seealso = 'graphDrawer'
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
graphDrawer <- function(data, plot_name, edge_clr, node_clrs,  bg_clr, 
                        lbl_fnt, legend_items, directory, fname,
                        ntwrks_page, col, H, W){
  
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
  
  net = igraph::graph_from_incidence_matrix(blanked_data, weight=T)
  deg = igraph::centr_degree(net, mode="all")
  
  net <- set_colors(x = blanked_data, net, node_clrs, VNames = VNames)
  
  V(net)$label.color <- 'black'
  V(net)$size = 5*sqrt(deg$res) 
  E(net)$width = E(net)$weight/4
  template <- layout_in_circle(net)
  
  png(filename,
      width = dims$W, height = dims$H, units = "px",pointsize = 12)
  
  par(mar=c(3,6.5,3,6.5))
  plot(net, layout=template, 
       edge.color = edge_clr,
       vertex.label.cex = 1.5,
       vertex.label = VNames,
       vertex.label.dist= VLPs$dist_vals, 
       vertex.label.degree = VLPs$degree_shift, 
       label.font = lbl_fnt#,
       #main = plot_name
       ) 
#    legend(x= -0.25, y=-1.35, legend_items, 
#           pch=21, col="#777777", 
#           pt.bg=node_clrs, 
#           pt.cex=2, cex=.8, bty="n", ncol=1)
    
 invisible(dev.off())
 
 message(paste0("'", plot_name, 
                "' has been rendered as a graph and saved to:\n ", filename))

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
#' @return a cowplot arranged grid.

netPage <- function(directory, col_var, row_var, col, fname){
  
  if(missing(directory)) { directory <- 'NetworkGraphs' }
  if(missing(fname)) {fname <- 'mosaiced_Nets.pdf'}
  
  rowOrder <- matrix( # this sets up a matrix to have rows groups in order
    data = rep(row_var, each = length(col_var)),
    ncol = length(col_var))
  colOrder <- matrix( # this sets up a matrix to have column groups in order
    data = rep(col_var, each = length(col_var)),
    ncol = length(col_var), byrow = F)
  image_orders <- matrix(paste0(colOrder, '-', rowOrder, '.png'), 
                   ncol = length(col_var))
  image_orders <- file.path(directory, image_orders)
  
  images = lapply(image_orders, png::readPNG)
  grob_images = lapply(images, grid::rasterGrob)
  
  
  # we might need to arrange each row seperately in order to add the 
  # title to the left side. 
  
  pdf(file = file.path(directory, fname), paper = 'a4')
  plot <- gridExtra::grid.arrange(grobs= grob_images, ncol = length(col_var), 
                          top="Bill Walton", left="mid", bottom = 'molecular')

#  grid::grid.text(c("right", "centre"), 3, 2)
  dev.off()
  
  return(plot)
}


r <- netPage(col_var = c('observation', 'microscopy', 'molecular'), 
        row_var = c('early', 'mid', 'late'))



# probably need to arrange each row, less the last, into their own 'grob', and
# then stack those on top of each other in order to label the left side of each.
# then add three different plots to the bottom  so that each may be labelled at
# the bottom

col_var <- c('aC', 'bC', 'cC')
row_var <- c('aR', 'bR', 'cR', 'dR')
nPlots <- (length(col_var) * length(row_var))
colN <- length(col_var)
rowN <- nPlots / length(col_var)


# need to populate this empty matrix with values 


layout <- matrix(nrow = rowN, ncol = colN, byrow = T,
                 data = c(rep(1:(rowN -1), each = colN),
                          rowN:((rowN-1) + colN) )
)


x <- 1:(nPlots - rowN) +1
res <- x[!x%%colN]




g1 <- arrangeGrob(grobs = grob_images[1:colN], ncol = colN, left = row_var)
g2 <- arrangeGrob(grobs = grob_images, ncol = nCol, left = row_var)
g3 <- arrangeGrob(grobs = grob_images, ncol = nCol, left = row_var)

grid.arrange(grobs = gs, layout_matrix = layout)





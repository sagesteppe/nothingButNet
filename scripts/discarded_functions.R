# yes you read it right... that is how much a production this has been, about a half dozen functions fell before I even initialized this repo... So...



vert_lab_position <- function(x){
  
  #' Takes a single input, a network, this will serve to specify the positions
  #' for the labels of objects in it. It does not modify the label values, just
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
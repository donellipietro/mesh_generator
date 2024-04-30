# Generates a grid (either hexagonal or square) of step h on the rectangular 
# box identified by bbox

generate_grid <- function(bbox, h, seed_point, type = "hexagonal") {
  
  ## get boundary box
  xmin <- bbox$xmin - 2*h
  xmax <- bbox$xmax + 2*h
  ymin <- bbox$ymin - 2*h
  ymax <- bbox$ymax + 2*h
  
  ## set defaults 
  if(is.null(seed_point)){
    seed_point <- st_as_sf(data.frame(x = xmin, y = ymin), coords = c("x", "y"))
  }
  if(type == "square"){
    type = "regular"
  }
  
  # create a polygon representing the rectangle
  bbox_sf <- st_polygon(list(as.matrix(cbind(c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin)))))
  
  # create a grid within the rectangle
  grid <- st_make_grid(bbox_sf, what = "centers", square = F, cellsize = h)
  
  # Find the index of the point with the minimum distance from the seed
  distances <- dist_point_from_points(seed_point, grid)
  closest_point_index <- which.min(distances)
  closest_point <- grid[closest_point_index, ]
  translation <- data.frame(st_coordinates(seed_point) - st_coordinates(closest_point))
  
  # Translation to match the seed
  grid <- data.frame(st_coordinates(grid))
  grid$X <- grid$X + translation$X
  grid$Y <- grid$Y + translation$Y
  grid <- st_as_sf(grid, coords = c("X", "Y"))

  return(grid)
}

# Generates a lattice (either hexagonal or square) on the smallest rectangular
# box containing the locations points cloud.
# The polygons that do not contain any point are then discarded and the domain 
# is generated as the union of the remaining polygons.

generate_lattice <- function(locations, h, bbox, seed_point, type = "hexagonal") {  if(is.null(bbox)){
    bbox <- attributes(locations$geometry)$bbox
  }
  
  grid <- generate_grid(bbox, h, seed_point, type)
  
  check <- c()
  polygons <- list()
  for(i in 1:nrow(st_coordinates(grid))){
    point <- grid[i,]
    if(type == "hexagonal"){
       polygon <- hex(point, h/sqrt(3) + 1e-9)
    } else if(type == "square") { 
      polygon <- square(point, h/sqrt(2) + 1e-9) 
    }
    check[i] <- any(apply(st_intersects(locations, polygon), 1, any))
    if (check[i]) {
      polygons <- c(polygons, list(polygon))
    }
  }
  
  
  # Domain
  lattice <- st_multipolygon(polygons)
  domain <- st_union(lattice, by_feature = TRUE)
  
  return(list(grid = grid,
              lattice = lattice,
              domain = domain))
}

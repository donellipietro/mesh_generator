# polygons ----

# Generates a hexagon given the central point and the radius

hex <- function(p, h) {
  alpha_vect <- seq(pi/6, 2*pi+pi/6, by = pi/3)
  
  ## p coordinates
  p_coords <- st_coordinates(p)
  
  # calculate coordinates of the hexagon vertices
  hex_coords <- data.frame(x = p_coords[, "X"] + h * cos(alpha_vect),
                           y = p_coords[, "Y"] + h * sin(alpha_vect))
  
  # create an sf object for the hexagon
  hex_sf <- st_polygon(list(as.matrix(hex_coords)))
  
  return(hex_sf)
}


# metrics ----

# returns the distance between a point p and all the points contained in points
dist_point_from_points <- function(p, points){
  
  ## p coordinates
  p_coords <- st_coordinates(p)
  points_coords <- st_coordinates(points)
  
  ## distance
  dist <- sqrt((points_coords[, "X"] - p_coords[, "X"])^2 +
               (points_coords[, "Y"] - p_coords[, "Y"])^2)
  
  return(dist)
}

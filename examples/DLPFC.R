# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Example: DLPFC mesh generation %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
VERBOSE <- FALSE

name_example <- "DLPFC"


# libraries ----
suppressMessages(library(sf))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(raster))
suppressMessages(library(rmapshaper))
suppressMessages(library(fdaPDE))


# sources ----
source("src/utils/cat.R")
source("src/utils/directories.R")
source("src/utils/plots.R")
source("src/meshing.R")
source("src/plots.R")

# paths ----

path_images <- paste("images/", name_example, "/", sep = "")
mkdir(path_images)

# global variables ----
figure_width <- 7.5
figure_height <- 8


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cat.script_title("Example: DLPFC mesh generation")


# data ----

## load data
load("data/DLPFC/DLPFC_sample9.RData")
rm(counts, true_labels)

## create sp object
locations <- SpatialPoints(locations)

# range(locations@coords[,1])
xmin <- 125
xmax <- 510 
# range(locations@coords[,2])
ymin <- -515
ymax <- -110 


# Plot
plot <- ggplot() +
  xlim(xmin, xmax) + ylim(ymin, ymax) + mesh_generation_plot_settings() + 
  geom_sf(data = st_as_sf(locations), color = "black", size = 0.1)
ggsave(paste(path_images, "0_initial_locations.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")


# domain identification ----

## grid ----

# Lattice type
type = "hexagonal"

# Why an hexagonal grid?
# https://strimas.com/post/hexagonal-grids/

# Regular hexagons are the closest shape to a circle that can be used for the
# regular tessellation of a plane and they have additional symmetries compared 
# to squares. These properties give rise to the following benefits.

# - Reduced edge effects: a hexagonal grid gives the lowest perimeter to area 
#   ratio of any regular tessellation of the plane. In practice, this means that 
#   edge effects are minimized when working with hexagonal grids. This is 
#   essentially the same reason beehives are built from hexagonal honeycomb: it 
#   is the arrangement that minimizes the amount of material used to create a 
#   lattice of cells with a given volume.
# - All neighbours are identical: square grids have two classes of neighbours, 
#   those in the cardinal directions that share an edge and those in diagonal 
#   directions that share a vertex. In contrast, a hexagonal grid cell has six 
#   identical neighbouring cells, each sharing one of the six equal length 
#   sides. Furthermore, the distance between centroids is the same for all 
#   neighbours.
# - Better fit to curved surfaces: when dealing with large areas, where the 
#   curvature of the earth becomes important, hexagons are better able to fit 
#   this curvature than squares. This is why soccer balls are constructed of 
#   hexagonal panels.
# - They look badass: it can’t be denied that hexagonal grids look way more 
#   impressive than square grids!

# Step of the grid
h <- 6.4
# It is decided by the user by looking at the initial distribution of location. 
# Lower the value of h more precise will be the domain reconstruction.
# However, it can not be too low otherwise there could be unwelcome holes in 
# the domain

# Seed Point
seed_point <- SpatialPoints(data.frame(x = 235, y = -126))
# It is the seed for the generation of the grid, the final grid is guaranteed to
# contain this point. It is used to have always the same grid.

# Plot
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(locations), color = "black", size = 0.1) +
  geom_sf(data = st_as_sf(seed_point), color = "red", size = 0.15) +
  geom_sf(data = st_as_sf(hex(seed_point, h/sqrt(3))), fill = "blue", alpha = 0.5, color = "black")
ggsave(paste(path_images, "1_initial_locations_check.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

# Grid generation
grid <- generate_grid(locations@bbox, h, seed_point, type = type)

# Plot
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(grid), color = "blue", size = 0.05) +
  geom_sf(data = st_as_sf(seed_point), color = "red", size = 1)
ggsave(paste(path_images, "2_grid.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")


## Lattice ----
## ||||||||||||

check <- c()
polygons <- list()
polygons_all <- list()
for(i in 1:nrow(grid@coords)){
  point <- grid[i,]
  polygon <- hex(point, h/sqrt(3) - 1e-9)
  check[i] <- any(!is.na(over(locations, polygon)))
  if(check[i]){
    polygons <- c(polygons, polygon@polygons[[1]]@Polygons[[1]])
  }
  polygons_all <- c(polygons_all, polygon@polygons[[1]]@Polygons[[1]])
}

lattice <- SpatialPolygons(list(Polygons(polygons, ID = "hex_lattice")))
lattice_all <- SpatialPolygons(list(Polygons(polygons_all, ID = "hex_lattice")))

# Plot all
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice_all), fill = "blue", alpha = 0.5, color = "black")
ggsave(paste(path_images, "3_lattice.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

# Plot all and locations
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice_all), fill = "blue", alpha = 0.5, color = "black") +
  geom_sf(data = st_as_sf(locations), color = "black", size = 0.1)
ggsave(paste(path_images, "4_lattice_and_locations.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

# Plot selected
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice), fill = "blue", alpha = 0.5, color = "black")
ggsave(paste(path_images, "5_lattice_only_selected.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")


## Simplification ----
## |||||||||||||||||||

# User defiend parameter
simplification <- 0.25
# This parameter represents the percentage of boundary vertices to be kept.
# The user should set a value such that the boundary is simplified enough
# but without exeeding otherwise there will be a lot of discarded points

# About holes
remove_holes <- FALSE
minimum_area_hole <- NULL
simplification_hole <- 0.5

# Lattice
lattice <- generate_lattice(locations, h, locations@bbox, seed_point, type = type)

# Plot
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice$domain), fill = "grey", alpha = 0.5, color = "black")
ggsave(paste(path_images, "6_domain.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

# Simplification
lattice_simplified <- simplify_domain(lattice, simplification,
                                      remove_holes, minimum_area_hole, simplification_hole)

# During this step the islands, namely the parts of domain that are not joined
# to the one with largest, area are discarded

# Plot
plot <- plot.discretized_domain(lattice_simplified, size = 0.1) +
  mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) + guides(color = "none")
ggsave(paste(path_images, "7_domain_simplified.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

# Discarded points
indexes.discarded_locations <- is.na(over(locations, lattice_simplified$domain))
locations.final <- locations[!indexes.discarded_locations,]

# Plot
plot <- plot.final_locations(locations, SpatialPoints(locations.final), lattice_simplified, size = 0.05) +
  mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste(path_images, "8_final_locations.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")


# Meshing ----

# User defiend parameter
maximum_area <- 50
# It is the threshold for the larger possible value for an element of the mesh.
# Is the original mesh contain elements larger than it it is refined until all
# the elements meet this constraint.
# The user should set a value such that the final number of nodes of the mesh
# has the same order of magnitude of the number of locations.

# Mesh generation
mesh <- generate_mesh(lattice_simplified, maximum_area)

# Plot
plot <- plot.fdaPDE_mesh(mesh) +
  mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste(path_images, "9_mesh.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

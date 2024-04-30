# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Example: DLPFC mesh generation %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

name_example <- "DLPFC"


# libraries ----
suppressMessages(library(sf))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(raster))


# sources ----
source("src/utils/cat.R")
source("src/utils/directories.R")
source("src/utils/plots.R")
source("src/geometry.R")
source("src/grids.R")

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

## create sf object
locations <- st_as_sf(locations, coords = c("x", "y"))

## get boundary box
bbox <- attributes(locations$geometry)$bbox
xmin <- floor(bbox$xmin) - 3 
xmax <- ceiling(bbox$xmax) + 3 
ymin <- floor(bbox$ymin) - 3 
ymax <- ceiling(bbox$ymax) + 3 

# Plot
plot <- ggplot() +
  xlim(xmin, xmax) + ylim(ymin, ymax) + mesh_generation_plot_settings() + 
  geom_sf(data = locations, color = "black", size = 0.1)
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
seed_point <- st_as_sf(data.frame(x = 235, y = -126), coords = c("x", "y"))
# It is the seed for the generation of the grid, the final grid is guaranteed to
# contain this point. It is used to have always the same grid.

# Plot
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = locations, color = "black", size = 0.1) +
  geom_sf(data = seed_point, color = "red", size = 0.1) +
  geom_sf(data = hex(seed_point, h/sqrt(3)), fill = "blue", alpha = 0.5, color = "black")
ggsave(paste(path_images, "1_initial_locations_check.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

# Grid generation
grid <- generate_grid(bbox, h, seed_point, type = type)

# Plot
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = grid, color = "blue", size = 0.05) +
  geom_sf(data = seed_point, color = "red", size = 0.05)
ggsave(paste(path_images, "2_grid.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")


## Lattice ----
## ||||||||||||

check <- c()
polygons <- list()
polygons_all <- list()
for(i in 1:nrow(st_coordinates(grid))){
  point <- grid[i,]
  polygon <- hex(point, h/sqrt(3) - 1e-9)
  check[i] <- any(apply(st_intersects(locations, polygon), 1, any))
  if (check[i]) {
    polygons <- c(polygons, list(polygon))
  }
  polygons_all <- c(polygons_all, list(polygon))
}

lattice <- st_multipolygon(polygons)
lattice_all <- st_multipolygon(polygons_all)

# Plot all
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = lattice_all, fill = "blue", alpha = 0.5, color = "black")
ggsave(paste(path_images, "3_lattice.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

# Plot all and locations
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = lattice_all, fill = "blue", alpha = 0.5, color = "black") +
  geom_sf(data = locations, color = "black", size = 0.1)
ggsave(paste(path_images, "4_lattice_and_locations.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

# Plot selected
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = lattice, fill = "blue", alpha = 0.5, color = "black")
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
lattice <- generate_lattice(locations, h, bbox, seed_point, type = type)

# Plot
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = lattice$domain, fill = "grey", alpha = 0.5, color = "black")
ggsave(paste(path_images, "6_domain.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")
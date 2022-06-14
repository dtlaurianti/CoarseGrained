################################################################################
#
# make_smoothdata.R
# script to 
#
################################################################################

#setwd("~/CoarseGrained") # to be updated later
install.packages("tidyverse", repos = "http://cran.us.r-project.org")
install.packages("geometry", repos = "http://cran.us.r-project.org")

library(tidyverse) # install.packages("tidyverse")
library(geometry)  # install.packages("geometry")

args <- commandArgs(trailingOnly = TRUE)

#' @param samples, a data frame with columns x, y, and z. 
#' @param resolution, the resolution of the rectangular grid on which to compute the surface
#' @param bandwidth, the bandwidth of the smoother. Smaller bandwidth -> more wiggly smoothers
#' @return grid, a truncated grid of points on which the loess smoother is evaluated. 
smoothed_surface <- function(samples, resolution, bandwidth){
  # compute convex hull of x and y points
  ch <- convhulln(cbind(samples$x, samples$y))	
  # create the grid
  grid <- expand.grid(x = seq(min(samples$x), max(samples$x), by = resolution), 
                      y = seq(min(samples$y), max(samples$y), by = resolution), 
                      stringsAsFactors = FALSE) %>% 
    tibble()
  # for each point, check whether it's in the convex hull
  grid["in_grid"] <- inhulln(ch, cbind(grid$x, grid$y)) 
  grid <- grid %>% filter(in_grid) 
  # fit the smoother
  fit <- loess(z ~ x + y, data = samples, span = 0.05)
  # add smoother predictions to the grid
  grid <- broom::augment(fit, newdata = data.frame(grid)) %>% 
    mutate(z_ = .fitted) 
  # return the grid
  grid %>% 
    select(x, y, z_) %>%
    mutate(z = z_)
}

filename <- paste0("data/visualization_data/", args[1])

samples <- read_csv( paste0(filename, ".csv") )
the_grid <- smoothed_surface(samples, 0.1, 0.1)

write_csv(the_grid, paste0(filename, "_smooth.csv"))


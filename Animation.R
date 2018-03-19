########################################################
#
# This script makes an animation of a Vicsek model simulation
# by reading a csv file with information about every
# particle in every time instant.
#
# Dependecies needed:
# FFmpeg
# devtools (R package)
# gganimate (Installed via devtools)
#
# Author: Ivan A. Moreno Soto
# Last updated: 17/March/2018
#
########################################################

library(ggplot2)
library(gganimate)

########################################################

# We read the simulation data.
simulation <- read.csv('sim.csv', header=TRUE, sep=',')

vicsek <- ggplot(simulation, aes(x = px, y = py, frame = t)) + geom_point()
gganimate(vicsek, "vicsek_sim.mp4", ani.width=1600, ani.height=900, interval=0.05)

########################################################

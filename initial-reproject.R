#################################################################################
#################################################################################

# Channel-centric coordinate system reprojection
# by DR 
# Jan 22 2018

#################################################################################

#################################################################################
########## Summary of steps #####################################################

# 1. load data of interest. Input data, at a minimum, includes xy coords for  
# some points of interest

# 2. Subset thalweg points from dataset (or load points used to generate
# channel centerline)

# 3. If desired, subset other points (e.g. bank breaklines) to be reprojected
# and interpolated seperately 

# 4. Interpolate or sample thalweg to desired point density

# 5. Create line object from centerline points (psp format)

# 6. Extract/calculate lengths of line segments

# 7. run primary function which locates line segment closest to each point to
# be transformed, calculates perpendicular distance from that segment to 
# the point, and calculates how far along the segment that intersection 
# is. This is the basis for the new coordinates

# 8. merge all necessary data, calculate downstream coordinate by adding 
# distance from line segment vertices that the point of interest falls,
# to the cumulative downstream distance along the centerline. 

# clear workspace
rm(list = ls(all = TRUE))



# -- load data and define variables---------------------------------------------------------------------------------



# load packages
require(spatstat)
library(sp)

# load functions - three total here
source("centerline_interp.R")
source("cartesian_to_channelcentric.R")
source("channelcentric_to_cartesian.R")

# set working directory
setwd("C:/Users/geography/OneDrive/Work/Projects/Carnation/Data/SA bed survey data/MENV work/Processed data/SA2")

# load data
original_points = read.csv("SA2 1980 checked.csv")

# load centerline points
centerline = read.csv("sa2centerline.csv")

# these are distances, in m, which the centerline will be offset
x.offset = 5
y.offset = 15

# this is the densification factor for the line interpolation
interp_factor = 20

# plot centerline and data to be transformed to inspect

plot(centerline$POINT_X, centerline$POINT_Y,
     type = "l", col = "red",
     xlim = c(min(original_points$Easting)-20, 
              max(original_points$Easting)+20),
     ylim = c(min(original_points$Northing)-20, 
              max(original_points$Northing)+20))
points(original_points$Easting, original_points$Northing,
       col = "blue")
legend(bty = "n", "topright", pch = c(NA,21),
       lty = c(1,NA), col = c("red", "blue"),
       legend = c("centerline points", "untransformed points"))



## -- interpolate centerline to desired density ------------------------------------------------------------------------



# 1. load function which takes centerline and fits a spline, and offsets
# 2. inputs are the untransformed centerline, the offsets, and the point 
#   densification factor
# 3. function returns a list containing a psp line object,
#   vector of line segment lengths, cumulative lengths, and a dataframe
#   with x y coordinate vertices

# call function
lines = centerline_interp(centerline, 
                          x.offset, 
                          y.offset, 
                          interp_factor)

# extract data from the function output list - these are all needed 
#   for the next function

thal.line = lines$thal.line
line.lengths = lines$line.lengths
line.start.cumu = lines$line.start.cumu
simple_centerline = lines$simple_centerline



## -- project to channel-centric---------------------------------------------------------------------------------------



# this function converts cartesian xy data into data projected against a channel 
# centerline. Input requirements are:

# 1. the data in cartesian form, with colums of "easting", "northing" and "elevation
# 2. psp line file contianing information of centerline segments
# 3. vector showing cumulative distance along line segments present in psp line file

# output is a dataframe with new coordinates

# call function
coords_new = cartesian_to_channelcentric(original_points, 
                                         thal.line, 
                                         line.start.cumu)


# plot reprojected data (i.e. function output)
plot(coords_new$s.coord, coords_new$n.coord, 
     xlim = c(-10,50), ylim = c(110,10),
     xlab = "across-channel coordinate (m)", 
     ylab = "downstream coordinate (m)",
     main = "transformed")
abline(v = 0, lty = 2)



# -- convert back to cartesian from channel-centric ------------------------------------------------------------------


# run the function, inputs are the following:

# 1. coordinates in channel centric form
# 2. centerline points
# 3. cumulative distances along centerline
# 4. psp line object of centerline

# output: xyz data in cartesian form
data_out = channelcentric_to_cartesian(coords_new, 
                                       thal.line, simple_centerline, 
                                       line.start.cumu,
                                       line.lengths)


# plot new data and old data for comparison                                    

plot(data_out$new_x, data_out$new_y, 
     pch = 21, bg = "red", cex = 0.5)
points(original_points$Easting, original_points$Northing, 
       col = "green", cex = 1)
lines(simple_centerline$x, simple_centerline$y)







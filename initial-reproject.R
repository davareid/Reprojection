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

# load data and define variables---------------------------------------------------------------------------------


# load packages
require(spatstat)
library(sp)

# load functions
source("centerline_interp.R")

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
windows()
plot(centerline$POINT_X, centerline$POINT_Y, type = "l", col = "red",
     xlim = c(min(original_points$Easting)-20, max(original_points$Easting)+20),
     ylim = c(min(original_points$Northing)-20, max(original_points$Northing)+20))
points(original_points$Easting, original_points$Northing, col = "blue")
legend(bty = "n", "topright", pch = c(NA,21), lty = c(1,NA), col = c("red", "blue"),
       legend = c("centerline points", "untransformed points"))

## interpolate centerline to desired density -----------------------------------------------
# this is also where the function could begin


# load function which takes centerline and fits a spline, and offsets
# inputs are the untransformed centerline, the offsets, and the point densification factor
# function returns a list containing a psp line object,
# vector of line segment lengths, cumulative lengths, and a dataframe with 
# x y coordinate vertices
lines = centerline_interp(centerline, x.offset, y.offset, interp_factor)

# extract data from the function output list - these are all needed for the next function
thal.line = lines$thal.line
line.lengths = lines$line.lengths
line.start.cumu = lines$line.start.cumu
simple_centerline = lines$simple_centerline

## reproject--------------------------------------------------------------------------------

# load functions
source("cartesian_to_channelcentric.R")

coords_new = cartesian_to_channelcentric(original_points, thal.line, line.start.cumu)


# start function here --------------------------------------------------------------------------------------------------

#cartesian_to_channelcentric = function(original_points, thal.line, line.start.cumu){

## use reprojection function to locate nearest line segment and calculate new coords

# For each bank point, this loop will:
# 1. calculate the angles between the point and vertices of each line segmet
# 2. check to see if both of those angles are < 90 degrees
# 3. if both are <90 degrees then this is probably the closest segment
# 4. identify ID of line segment closest to point
# 5. output a matrix with bank point and line segment associated
# 6. calculate perpendicular distance from point to line segment (this is the n coord)
# 6. in cases where multiple segments result in two acute angles, 
# select the segment with the shorter n coordinate
# 7. Calculate the point along closest line segment where n coord intersects
# and restate this as distance from upstream line vertices
# 8. sort data by which line segment points are closest to, and add the distance
# to nearest vertex (see previous point) to cumulative distance for upstream 
# end of nearest line segment. You now have the s coordinate


# generate some empty matrices which are populated in function
closest.line = matrix(nrow = nrow(original_points), ncol = thal.line$n)
B.side.final = matrix(nrow = nrow(original_points), ncol = thal.line$n)
ba.ang.final = matrix(nrow = nrow(original_points), ncol = thal.line$n)
n.coord.dist = matrix(nrow = nrow(original_points), ncol = thal.line$n)
counter = seq(1,thal.line$n,1)
nearest_line = matrix(nrow = nrow(original_points), ncol = 2)
sides = matrix(nrow = nrow(original_points), ncol = 2)

# loop for performing calculations
for(i in 1:nrow(original_points)){
  for(j in 1:thal.line$n){
    
    # calculate triangle side lengths - B side is between point and upstream
    # vertex of line segment
    B.side.l = sqrt((original_points$Easting[i] - thal.line$ends$x0[j])^2 + 
                      (original_points$Northing[i] - thal.line$ends$y0[j])^2)
    
    # calculate line length between point and downstream line segment vertex
    C.side.l = sqrt((original_points$Easting[i] - thal.line$ends$x1[j])^2 + 
                      (original_points$Northing[i] - thal.line$ends$y1[j])^2)
    
    # calculate angles within triangle - this is upstream angle between B and centerline segment
    ba.ang = acos((B.side.l^2 + line.lengths[j]^2 - 
                     C.side.l^2)/(2*B.side.l*line.lengths[j]))
    
    # downstream angle (between C and centerline segment)
    ac.ang = acos((line.lengths[j]^2 + C.side.l^2 - 
                     B.side.l^2 )/(2*C.side.l*line.lengths[j]))
    
    # check to see if both angles are acute (units are radians)
    # if angles check out, extract necessary side lengths and angles for calculating downstream (n)
    # coordinate
    if(ba.ang < 1.575 & ac.ang < 1.575){
      
      # add a marker so that the correct segment can be idenfitied
      closest.line[i,j] = 1
      n.coord.dist[i,j] = B.side.l * sin(ba.ang) # calculate n coord
      B.side.final[i,j] = B.side.l # save this for later
      ba.ang.final[i,j] = ba.ang   # also save this for later calc
    }
    
    # replace NA's with a big number - makes next step easier
    else{n.coord.dist[i,j] = 1000}
  }
  
  # extract index of smallest s coord distance in matrix, in case of multiple segments meetint 
  # angle criteria
  nearest_line[i,1] =  which.min(n.coord.dist[i,])
  nearest_line[i,2] = min(n.coord.dist[i,])
  sides[i,1] = B.side.final[i,nearest_line[i,1]]
  ba.ang.final[i,1] = ba.ang.final[i,nearest_line[i,1]]
  
  # calculate s coord from B side length and ba.ang.final
  sides[i,2] = sides[i,1]*cos(ba.ang.final[i,1])
}

# re-organize and calculate final coordinates ---------------------------------------------------

# combine relevant data and sort by line number of nearest centerline segment
coords_new = as.data.frame(cbind(nearest_line, sides, seq(1,nrow(original_points),1)))

# sort
coords_new = coords_new[order(coords_new[,1]),]

# merge the cumulative line breaks with reprojected data
reproj_1 = cbind(counter, line.start.cumu)
reproj_2 = merge(reproj_1, nearest_line, by.x = 1, by.y = 1)

# merge the cumulative segment distance with the 
# along-segment distance
coords_new = cbind(coords_new, reproj_2$line.start.cumu)

# rename dataframe with identifiable column names
names(coords_new) = c("nearest.seg", "s.coord", "B.seg.length", 
                      "dist.along.seg", "point id", "line.seg.dist")

# calculate final downstream coordinate
n.coord = coords_new$line.seg.dist + coords_new$dist.along.seg

# merge with rest of the data
coords_new = cbind(coords_new, n.coord)

# test plot of reprojected data

# end function here 
#}



# ---------------------------------------------------------------------------------------------------------------------

windows()
plot(coords_new$s.coord, coords_new$n.coord, xlim = c(-10,50), ylim = c(110,10),
     xlab = "across-channel coordinate (m)", ylab = "downstream coordinate (m)",
     main = "transformed")
abline(v = 0, lty = 2)
#text(coords_new$s.coord, coords_new$n.coord, cex = 0.6, pos = 4, col = "black")


# Steps for re-projection of interpolated data -----------------------------------------------------

# these are conceptual steps which essentially work backwards from those given above. 
# they probably need revision but I think the concept should work out. 

# Going point by point... 
# 1. Given the downstream coordinate, locate the thalweg segment each new point is adjacent to
# 2. isolate the line segment vertex coordinates
# 3. Calculate how far from the upstream or downstream vertex a perpendicular line drawn from the point to the
# line segment is located
# 4. We can then calculate the distance from the unknown point to one of the vertices, and then to 
# the other vertex. Most importantly, we need the distance from the unknown point to the upstream 
# centerline vertex
# 5. Once you have all sides of the right angle triangle, calculate the slope of the centerline
#   segment from the segment's vertices
# 6. Calculate the angle between the centerline segment and the hypotenuse of the triangle
# 7. Using this angle and the slope of the centerline segment, find the slope of the hypotenuse
# 8. Using this slope and the cartesian coordinates of the upstream centerline segment vertex,
#   calculate the new coordinates. 


# --- define points to be transformed back to original coordinates
# this would be the post-interpolation points. Ultimately, all that is needed is the S and N coord 
# values

dat_interp = coords_new[complete.cases(coords_new),]

# somethign needed later 
counter = seq(1,thal.line$n,1)

# --- calculate some parameters related to the centerline segments

# calculate the slope of each centerline segment
centerline_slope = matrix(nrow = nrow(simple_centerline)-1, ncol = 1)

for(i in 1:nrow(simple_centerline)-1){
  centerline_slope[i] = (simple_centerline$y[i] - simple_centerline$y[i+1])/
    (simple_centerline$x[i] - simple_centerline$x[i+1])
}

# merge line ID, counter, and line length  and coords into single file
centerline_data = as.data.frame(cbind(counter, centerline_slope, line.lengths, line.start.cumu, thal.line$ends))
cumu_upstream = c(centerline_data$line.start.cumu[2:nrow(centerline_data)],1000)

# For each point to be interpolated back, figure out which centerline segment is closest
# generate matrix to be populated in loop

nearest_line = matrix(nrow = nrow(dat_interp), ncol = nrow(centerline_data))

# this loop finds the closest centerline segment to any point, given the n coord
# for the point and the vertices of the centerline segments. This is needed to 
# transform the data back

for(i in 1:nrow(nearest_line)){
  for(j in 1:nrow(centerline_data)){
  
  nearest_line[i,j] = ifelse(dat_interp$n.coord[i] > centerline_data$line.start.cumu[j]
                             & dat_interp$n.coord[i]<=cumu_upstream[j],
                             centerline_data$counter[j], NA)
  }
}

# join nearest_line to output data file
near_line = as.data.frame(rowSums(nearest_line, na.rm = TRUE))
names(near_line) = "near_segment"
dat_interp = cbind(dat_interp, near_line)

# merge data from line segments to the interpolated data file
# we should now have everything we need to retransform the data
#  back to original coordinates. This might need to be changed 
dat_interp = merge(dat_interp, centerline_data, by.x = 8, by.y = 1)

# --- Calculate items used for retransformation ---

# step 1: calculate the distance of the n coord for each point from the upstream limit of te 
# centerline segment, and the final side of the "triangle" made up of s coord and 
# the line just calculated

segment_sublength = dat_interp$n.coord - dat_interp$line.start.cumu
hypotenuse = sqrt(segment_sublength^2 + dat_interp$s.coord^2)

# step 2: calculate the angle made between the centerline segment and the 
# hypotenuse line

theta = acos((segment_sublength^2 + hypotenuse^2 - dat_interp$s.coord^2)/
               (2*segment_sublength*hypotenuse))

# Step 3: calculate the slope of the hypotenuse line
m_hypotenuse = tan((atan(dat_interp$centerline_slope))-(pi-theta))

# step 4: calculate coordinates of transformed points

new_x = dat_interp$x0 + hypotenuse*(1/(sqrt(1+m_hypotenuse^2)))
new_y = dat_interp$y0 + hypotenuse*(m_hypotenuse/(sqrt(1+m_hypotenuse^2)))

# plot new data and old data for comparison                                    
windows()
plot(centerline_data$x0, centerline_data$y0, xlim = c(min(centerline_data$x0)-20, max(centerline_data$x0)+20),
     ylim = c(min(centerline_data$y0)-20, max(centerline_data$y0)+20))
points(new_x, new_y, pch = 21, bg = "red", cex = 0.5)
points(original_points$Easting, original_points$Northing, col = "green", cex = 1)
lines(centerline_data$x0, centerline_data$y0)

# General notes ------------------------------------------------------------------------------------

# if we want to deal with all study area data at once, it needs to be split based on which side of
# the centerline it is on. This could be done by automatically generating a mask and subsetting
# data which falls either inside or outside, and then running the reprojection separately for both
# halves. 

# In some cases where points fall really close to the intersection of two line segments, it is hard
# to find instances where two acute angles are present. I am not sure why this is, exactly, and
# there is some discussion about this in the merwade papers. It is just a small number of points
# but ideally we would find a solution to this. I've dealt with it for now by raising the threshold
# for finding closest segment to slightly more than 90 degrees in case rounding or something is the 
# problem. It doesn't fix everything but seems to help. It also slightly moves the points off their
# "nearest" segment by assigning a negative s coordinat. 






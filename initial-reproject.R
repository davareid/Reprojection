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

# load data ---------------------------------------------------------------------------------


# load packages
require(spatstat)
library(sp)

# set working directory
setwd("C:/Users/geography/OneDrive/Work/Projects/Carnation/Data/SA bed survey data/MENV work/Processed data/SA2")

# load data
dat = read.csv("SA2 1988 checked.csv")

# subset thalweg
thalweg = subset(dat, Thalweg == "TRUE")

# either subset data or use all data. 
original_points = subset(dat, Banks == "TR")


## interpolate centerline to desired density -----------------------------------------------


# this resamples to a lower point spacing - but may be unnecessary if a different
# approach is taken to generating the centerline
simple_centerline = approx(thalweg$Easting,
                           thalweg$Northing, method="linear", 
                           n=nrow(original_points)/4, rule = 1, f = 0, ties = mean)

# plot centerline and data to be transformed
windows()
plot(simple_centerline$x, simple_centerline$y, type = "l", col = "red")
points(original_points$Easting, original_points$Northing, col = "blue")
legend(bty = "n", "topright", pch = c(NA,21), lty = c(1,NA), col = c("red", "blue"),
       legend = c("centerline points", "untransformed points"))



# Prep data for reprojection function ------------------------------------------------------


# reformat line object
# extract xy from thalweg data
simple_centerline = as.data.frame(cbind(simple_centerline[[1]],
                                        simple_centerline[[2]]))

# rename dataframe
names(simple_centerline) = c("x", "y")

# define observation window for format conversion (needed as input for psp object)
thal.window = owin(c(min(thalweg$Easting)-10, 
                     max(thalweg$Easting)+10), c(min(thalweg$Northing)-10, 
                                                 max(thalweg$Northing)+10))

# create a psp object (similar to SpatialLines), inputs are coords for 
# line ends and the window object
thal.line = psp(simple_centerline$x[1:nrow(simple_centerline)-1], 
                simple_centerline$y[1:nrow(simple_centerline)-1],
                simple_centerline$x[2:nrow(simple_centerline)], 
                simple_centerline$y[2:nrow(simple_centerline)], thal.window)

# calculate lengths of individual line segments from psp centerline object
line.lengths = lengths.psp(thal.line)

# calculate cumulative segment start positions - each n coord calculated will be added
# to the appropriate break location 
# empty vector
line.start.cumu = line.lengths*0

for(i in 2:length(line.lengths)){
  line.start.cumu[i] = line.lengths[i-1]+line.start.cumu[i-1]
}


## reproject--------------------------------------------------------------------------------


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
    if(ba.ang < 1.6 & ac.ang < 1.6){
      
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

windows()
plot(coords_new$s.coord, coords_new$n.coord, xlim = c(-10,20), ylim = c(100,0),
     xlab = "across-channel coordinate (m)", ylab = "downstream coordinate (m)",
     main = "transformed")
abline(v = 0, lty = 2)
text(coords_new$s.coord, coords_new$n.coord, cex = 0.6, pos = 4, col = "black")


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
# this would be the post-interpolation points 

dat_interp = coords_new[complete.cases(coords_new),]

# --- calculate some parameters related to the centerline segments

# calculate the slope of each centerline segment
centerline_slope = matrix(nrow = nrow(simple_centerline)-1, ncol = 1)

for(i in 1:nrow(simple_centerline)-1){
  centerline_slope[i] = (simple_centerline$y[i] - simple_centerline$y[i+1])/
    (simple_centerline$x[i] - simple_centerline$x[i+1])
}

# merge line ID, counter, and line length  and coords into single file
centerline_data = as.data.frame(cbind(counter, centerline_slope, line.lengths, line.start.cumu, thal.line$ends))

# For each point to be interpolated back, figure out which centerline segment is closest
# generate matrix to be populated in loop

nearest_line = matrix(nrow = nrow(dat_interp), ncol = nrow(centerline_data))

for(i in 1:nrow(dat_interp)){
  for (j in 1:nrow(centerline_data)-1){
    if(dat_interp$n.coord[i] > centerline_data$line.start.cumu[j] & dat_interp$n.coord[i]<=centerline_data$line.start.cumu[j+1]){
      nearest_line[i,j] = centerline_data$counter[j]
    } else{nearest_line[i,j] = NA}
  }
}


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






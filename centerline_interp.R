centerline_interp = function(centerline, x.offset, y.offset, interp_factor){
  
  # this resamples to a lower point spacing - but may be unnecessary if a different
  # approach is taken to generating the centerline
  
  simple_centerline = centerline[,3:4]
  names(simple_centerline) = c("x", "y")
  
  # offset centerline so that entire channel is to one side of it.
  # offset values will need to change each time
  simple_centerline$x = simple_centerline$x - x.offset
  simple_centerline$y = simple_centerline$y - y.offset
  
  # Prep data for reprojection function ------------------------------------------------------
  
  
  # generate a spline interpolation to create short line segments
  x = simple_centerline$x
  spline = spline(simple_centerline$x, y = simple_centerline$y, n = interp_factor*length(x), method = "fmm",
                  xmin = min(x), xmax = max(x), ties = mean)
  
  simple_centerline = cbind(spline$x, spline$y)
  
  # reformat line object
  # extract xy from thalweg data
  simple_centerline = as.data.frame(cbind(simple_centerline[[1]],
                                          simple_centerline[[2]]))
  
  simple_centerline = as.data.frame(cbind(spline$x, spline$y))
  
  # rename dataframe
  names(simple_centerline) = c("x", "y")
  
  # define observation window for format conversion (needed as input for psp object)
  thal.window = owin(c(min(simple_centerline$x)-100, 
                       max(simple_centerline$x)+100), c(min(simple_centerline$y)-100, 
                                                        max(simple_centerline$y)+100))
  
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
  
  # return the following variables
  return(thal.line)
  
} # end function



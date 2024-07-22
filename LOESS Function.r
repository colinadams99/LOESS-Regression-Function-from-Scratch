# ----- LOESS function -----
myloess <- function(x, y, span = 0.5, degree = 1, show.plot = TRUE){
  
  
  # Sorting data frame by x
  data = data.frame(x,y)
  colnames(data) = c("x","y")
  data = data[order(x),]
  
  span = span # Span
  degree = degree # Degree
  N_total = length(x) # Total number of points
  n_points = round(N_total * span, digits=0) # Number of points per window (rounded down to nearest integer)
  Win_total = N_total # Number of windows
  
  # How to do LOESSplot
  
  
  fitted_value = data.frame(matrix(ncol=1, nrow=0)) # This will hold fitted value for each point
  colnames(fitted_value) = c('fitted')
  
  # Window loop with every point in data
  for (index in 1:N_total)
  {
    focal_x = data[index,1] # X-coordinate of focal point, changes every time
    Distance = NULL # Resets each time focal point changes
    
    # Calculates distance between focal point and all points
    for (index in 1:N_total)
    {
      Distance[index] = euc.dist(focal_x, data[index,1])
    }
    data_distance = cbind(data, Distance)
    data_distance = data_distance[order(Distance),] # Sorts points by distance for subsetting
    
    subset = data_distance[1:n_points,] 
    
    # Calculating Scaled Distance and Weight per subset
    max_distance = max(subset$Distance)
    subset['Scaled Distance'] = (subset$Distance / max_distance)
    subset['Weight'] = (1 - ((abs(subset$`Scaled Distance`))^3))^3
    
    if (degree == 1)
    {
      wls = lm(y~x, weights = subset$Weight, data = subset)
      local_regression_value = wls$coefficients[1] + (focal_x * wls$coefficients[2])
    }
    if (degree == 2)
    {
      wls = lm(y~(x + I(x^2)), weights = subset$Weight, data = subset) # Weighted Least-Squares for each subset
      local_regression_value = wls$coefficients[1] + (focal_x * wls$coefficients[2]) + (((focal_x)^2)  * wls$coefficients[3])
    }
    
    # Bind each local fitted value into vector of fitted values
    fitted_value = rbind(fitted_value,local_regression_value)
    colnames(fitted_value) = c('fitted')
  }
  
  linedata = data.frame(cbind(data$x, data$y, fitted_value$fitted))
  colnames(linedata) = c('x', 'y', 'fitted')
  
  residuals = linedata$y - linedata$fitted
  SSE = (t(residuals) %*% residuals)[1]
  
  loessplot = ggplot(linedata, aes(x = x, y = y)) + theme_bw() + geom_point(color = "grey") + geom_line(aes(x, fitted))
  
  return(list(
    "span" = span,
    "degree" = degree,
    "N_total" = N_total,
    "Win_total" = Win_total,
    "n_points" = n_points,
    "SSE" = SSE,
    "loessplot" = loessplot))
}

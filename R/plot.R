#----------
# PLOT                   
#----------

plot_Output <- function(res) {
  y1 = min(min(res$fitness_average),min(res$fitness_best))
  y2 = max(max(res$fitness_average),max(res$fitness_best))
  
  if (y1 >= 0) {
    ymin = y1*0.8
    ymax = y2*1.2
  }
  if (y1 <= 0 & y2 <= 0) {
    ymin = y1*1.2
    ymax = y2*0.8
  }   
  
  if (y1 <= 0 & y2 >= 0) {
    ymin = y1*1.2
    ymax = y2*1.2
  }
  
  plot(fitness_average ~ generation,
       type = "l", data = res, col = "red",
       ylim = (c(ymin, ymax)), ylab = "Fitness", xlab="Iteration"
  )
  title("QGA - Optimization")
  points(fitness_best ~ generation, type = "l", data = res, col = "blue")
  legend("bottomright", c("Average fitness","Best fitness"), fill = c("red", "blue"), cex=0.8)
}  
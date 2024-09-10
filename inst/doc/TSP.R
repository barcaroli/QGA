#--------------------------------------
# Traveller Salesman problem  with QGA
#--------------------------------------
library(QGA)


## -----------------------------------------------------------------------------------------------
TravellingSalesman <- function(solution,distance) {
  l = 0.0  
  for (i in 2:length(solution)) {
    l = l+distance[solution[i-1], solution[i]]
  }
  # Fitness function
  l = l + distance[solution[1],solution[length(solution)]]
  # Penalization
  penal <- ((nrow(distance)) - length(table(solution)))*sum(distance)/10
  cost <- -(penal+l)
  return(cost)
}


## -----------------------------------------------------------------------------------------------
cities <- read.csv("TSP_cities.csv")
ncities <- 9
cities <- cities[c(1:ncities),]
cities


## -----------------------------------------------------------------------------------------------
distance <- as.matrix(dist(cities[,c(2:3)]))
distance


## ----eval=FALSE---------------------------------------------------------------------------------
popsize = 20
Genome = nrow(cities)
nvalues_sol = nrow(cities)
set.seed(4321)
TSPsolution <- QGA(popsize,
                generation_max = 1000,
                nvalues_sol,
                Genome,
                thetainit = 3.1415926535 * 0.01,
                thetaend = 3.1415926535 * 0.01,
                # pop_mutation_rate_init = 1/(popsize + 1),
                # pop_mutation_rate_end = 1/(popsize + 1),
                # mutation_rate_init = 1/(Genome + 1),
                # mutation_rate_end = 1/(Genome + 1),
                mutation_flag = FALSE,
                plotting = TRUE,
                verbose = FALSE,
                progress = FALSE,
                eval_fitness = TravellingSalesman,
                eval_func_inputs = distance)


## -----------------------------------------------------------------------------------------------
solution <- TSPsolution[[1]]
cities$city[solution]
cities_tsp <- cities[solution,]


## ----message=F, warning=FALSE-------------------------------------------------------------------
# Plot
if (!require(maps)) install.packages("maps")
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(ggplot2)) install.packages("ggplot2")
library(maps)
library(tidyverse)
library(ggplot2)
italy_map <- map_data("italy")
lines <- data.frame(
  x = cities_tsp$x[-nrow(cities_tsp)],
  y = cities_tsp$y[-nrow(cities_tsp)],
  xend = cities_tsp$x[-1],
  yend = cities_tsp$y[-1]
)
gg <- ggplot() +
  geom_polygon(data = italy_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black") +
  geom_segment(data = lines, aes(x = x, y = y, xend = xend, yend = yend),
               color = "red", size = 1) +
  geom_point(data = cities_tsp, aes(x = x, y = y), color = "blue", size = 3) +
  geom_text(data = cities_tsp, aes(x = x, y = y, label = city), 
            hjust = 0, vjust = 0, nudge_x = 0.1, nudge_y = 0.1) +
  coord_fixed(1.3) +
  theme_minimal() +
  ggtitle("Best itinerary") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
  xlim(6, 19) + ylim(36, 47)
print(gg)


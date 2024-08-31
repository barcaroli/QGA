#---------------------------------------------
# Application of the Quantum Genetic Algorithm
# to the Traveler Salesman Problem
#---------------------------------------------

#----------------------
library(QGA)
#----------------------

#-------------------------------------------------
# Fitness evaluation for Travelling Salesman Problem
#-------------------------------------------------
TravellingSalesman <- function(solution,distance) {
  l = 0.0  
  for (i in 2:length(solution)) {
    l = l+distance[solution[i-1], solution[i]]
  }
  l = l + distance[solution[1],solution[length(solution)]]
  # Penalization for solution with repeated cities
  penal <- (nrow(distance) - length(table(solution))) / nrow(distance) * sum(distance) * 0.25
  cost <- -(penal+l)
  return(cost)
}

#-------------------------------------------------
# Prepare data for fitness evaluation
cities <- read.csv("TravelerSalesman_cities.csv")
ncities <- 8
cities <- cities[c(1:ncities),]
distance <- as.matrix(dist(cities[,c(2:3)]))
#----------------------

#----------------------
# Perform optimization
popsize = 20
Genome = nrow(cities)
nvalues_sol = nrow(cities)
set.seed(4321)
solutionQGA <- QGA(
                popsize,
                generation_max = 2000,
                nvalues_sol,
                Genome,
                thetainit = 3.1415926535 * 0.025,
                thetaend = 3.1415926535 * 0.0025,
                pop_mutation_rate_init = 5/(popsize + 1),
                pop_mutation_rate_end = 0/(popsize + 1),
                mutation_rate_init = 1/(Genome + 1),
                mutation_rate_end = 0/(Genome + 1),
                mutation_flag = TRUE,
                plotting = TRUE,
                verbose = FALSE,
                eval_fitness = TravellingSalesman,
                eval_func_inputs = distance)

#----------------------
# Analyze results
solution <- solutionQGA[[1]]
cities$city[solution]
cities_tsp <- cities[solution,]
# Plot
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

#-----------------------------------------
# Compare with classical genetic algorithm

library(genalg)
evaluate <- function(solution) {
  solution <- round(solution)
  l = 0.0  
  for (i in 2:length(solution)) {
    l = l+distance[solution[i-1], solution[i]]
  }
  l = l + distance[solution[1],solution[length(solution)]]
  penal <- ((nrow(distance)) - length(table(solution)))*sum(distance)/10
  cost <- penal+l
  return(cost)
}
solutionGA <- rbga(stringMin=c(rep(1,nrow(cities))), 
                     stringMax=c(rep(nrow(cities),nrow(cities))),
                     popSize=20, 
                     iters=1000, 
                     elitism=NA, 
                     evalFunc=evaluate)
plot(solutionGA)
filter = solutionGA$evaluations == min(solutionGA$evaluations)
bestObjectCount = sum(rep(1, solutionGA$popSize)[filter])
if (bestObjectCount > 1) {
  bestSolution = solutionGA$population[filter, ][1, 
  ]
} else {
  bestSolution = solutionGA$population[filter, ]
}
bestSolution <- round(bestSolution)
cities$city[bestSolution]
l = 0.0  
for (i in 2:length(bestSolution)) {
  l = l+distance[bestSolution[i-1], bestSolution[i]]
}
l = l + distance[bestSolution[1],bestSolution[length(bestSolution)]]
l
cities_tsp <- cities[bestSolution,]
# Plot
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
  ggtitle("Best itinerary (genalg)") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
  xlim(6, 19) + ylim(36, 47)
print(gg)



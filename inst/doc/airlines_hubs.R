## ----setup, include = FALSE---------------------------------------------------------------------
options(width = 999)
knitr::opts_chunk$set(fig.width=8, fig.height=6,
  collapse = TRUE,
  comment = "#>"
)
# if (!require(bookdown)) install.packages("bookdown", dependencies=TRUE)
# library(bookdown)
library(QGA)


## -----------------------------------------------------------------------------------------------
cities <- read.csv("Airlines_cities.csv")
print(cities)


## -----------------------------------------------------------------------------------------------
cost <- read.csv("Airlines_cost.csv",header=F)
str(cost)


## -----------------------------------------------------------------------------------------------
flow <- read.csv("Airlines_flow.csv",header=F)
str(flow)


## -----------------------------------------------------------------------------------------------
airline_hub <- function(solution,input) {
  cost <- input$cost
  flow <- input$flow
  penalization <- sum(cost)*0.0005
  obj <- -sum(apply(cost[,solution] * flow[,solution], 1, min)) / sum(cost)
  if (length(table(solution)) < length(solution)) {
    obj <- obj - penalization
  }
  # cat("\nSolution:",solution,"  obj: ",obj)
  return(obj)
}


## -----------------------------------------------------------------------------------------------
input <- list(cost=cost,flow=flow)
popsize = 20
Genome = 3 # Number of desired hubs
nvalues_sol = 25
set.seed(4321)
solutionQGA <- QGA(
                popsize,
                generation_max = 200,
                nvalues_sol,
                Genome,
                thetainit = 3.1415926535 * 0.025,
                thetaend = 3.1415926535 * 0.005,
                pop_mutation_rate_init = 3/(popsize + 1),
                pop_mutation_rate_end = 1/(popsize + 1),
                mutation_rate_init = 3/(Genome + 1),
                mutation_rate_end = 1/(Genome + 1),
                mutation_flag = TRUE,
                plotting = FALSE,
                progress = FALSE,
                verbose = FALSE,
                eval_fitness = airline_hub,
                eval_func_inputs = input)


## -----------------------------------------------------------------------------------------------
solution <- solutionQGA[[1]]
cities$City[solution]
sol <- cost[,solution] 
find_min_column <- function(row) {
  which.min(row)
}
min_column_vector <- apply(sol, 1, find_min_column)
print(min_column_vector)
cities$Hub <- cities$City[solution][min_column_vector]
print(cities[,c(2,5)])


## ----message=FALSE, warning=FALSE---------------------------------------------------------------
if (!require(tidyverse)) install.packages("tidyverse", dependencies=TRUE)
if (!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if (!require(maps)) install.packages("maps", dependencies=TRUE)
library(tidyverse)
library(ggplot2)
library(maps)
connections <- cities %>%
  left_join(cities %>% select(City, Lat, Long), by = c("Hub" = "City")) %>%
  rename(Hub_Lat = Lat.y, Hub_Long = Long.y) %>%
  select(City, Lat = Lat.x, Long = Long.x, Hub, Hub_Lat, Hub_Long)
hubs <- unique(cities$Hub)
hub_connections <- expand.grid(Hub1 = hubs, Hub2 = hubs, stringsAsFactors = FALSE) %>%
  filter(Hub1 < Hub2) %>%  # Evita duplicati e auto-connessioni
  left_join(cities %>% select(City, Lat, Long), by = c("Hub1" = "City")) %>%
  rename(Lat1 = Lat, Long1 = Long) %>%
  left_join(cities %>% select(City, Lat, Long), by = c("Hub2" = "City")) %>%
  rename(Lat2 = Lat, Long2 = Long)
usa_map <- map_data("state")
gg <- ggplot() +
  geom_polygon(data = usa_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
  coord_fixed(1.3) +
  xlim(-125, -66) + ylim(25, 50)
gg <- gg + geom_segment(data = connections, 
                        aes(x = Long, y = Lat, xend = Hub_Long, yend = Hub_Lat),
                        color = "blue", alpha = 0.5, size = 1)
gg <- gg + geom_segment(data = hub_connections,
                        aes(x = Long1, y = Lat1, xend = Long2, yend = Lat2),
                        color = "red", alpha = 0.7, size = 1)
gg <- gg + geom_point(data = cities, aes(x = Long, y = Lat), color = "blue", size = 2)
gg <- gg + geom_point(data = cities[solution,], aes(x = Long, y = Lat), color = "red", size = 3)
gg <- gg + geom_text(data = cities, aes(x = Long, y = Lat, label = City), hjust = 1, vjust = 1, size = 3)
gg <- gg + ggtitle("Airlines hubs (cost minimization)") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
print(gg)


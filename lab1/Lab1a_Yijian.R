#Load packages
library(magrittr)
library(igraph)
library(vosonSML)

#Topic : San Jose passes first U.S. law requiring gun owners to get liability insurance and pay annual fee

#Part 1
# Reddit data collection
myRedditUrls <- c("https://www.reddit.com/r/moderatepolitics/comments/sd8rdt/san_jose_passes_first_us_law_requiring_gun_owners/", "https://www.reddit.com/r/politics/comments/sdgrg6/san_jose_passes_firstofitskind_insurance/","https://www.reddit.com/r/SanJose/comments/sd61jo/san_jose_becomes_1st_city_in_nation_to_approve/","https://www.reddit.com/r/politics/comments/sddk8e/san_jose_rules_gun_owners_must_have_liability/")
# authentication
redditData <- Authenticate("reddit") %>%
  Collect(threadUrls = myRedditUrls, waitTime = 5)
View(redditData)

# create an actor network with comment text as edge attribute
actorGraph <- redditData %>% Create("actor") %>% AddText(redditData) %>% Graph

## clean up the graph data removing self-loop 
edge_cleanup <- function(graph = actorGraph){
  library(igraph)
  df <- get.data.frame(actorGraph)
  names_list <- data.frame('name' = as.character(V(actorGraph)$name),
                           'label' = as.character(V(actorGraph)$label))
  df$from <- sapply(df$from, function(x) names_list$label[match(x,names_list$name)] %>% as.character())
  df$to <- sapply(df$to, function(x) names_list$label[match(x,names_list$name)] %>% as.character())
  nodes <- data.frame(sort(unique(c(df$from,df$to))))
  links <- df[,c('from','to')]
  net <- graph.data.frame(links, nodes, directed=T)
  E(net)$weight <- 1
  net <- igraph::simplify(net,edge.attr.comb="sum")
  return(net)
}

actorGraph <- edge_cleanup() # Runs the function to remove self-loops

# check if the network is directed or undirected
is.directed(actorGraph)

vcount(actorGraph) ## the number of nodes/actors/users

ecount(actorGraph) ## the number of edges

graph.density(actorGraph)# calculate the density of the network

save.image('Lab1a_Yijian.RData')

load('Lab1a_Yijian.RData')

#Part 2
# Calculate the number of components in the graph
comp <- components(actorGraph)
comp

actorGraph %>% 
  plot(.,
       margin = c(-0.2,-0.2,-0.2,-0.2),   ## values for the size of the bottom, left, top, and right plot margins

       vertex.size = 3,                   ## node size
       vertex.color = 'brown1',           ## node color
       
       vertex.label = NA,                 ## uncomment  this line to remove node labels
       vertex.label.cex = .2,             ## node label size
       vertex.label.color = 'darkcyan',   ## node label color

       edge.arrow.size = .3,              ## arrow size
       edge.color = 'darkgray',           ## arrow color

      layout = layout_with_dh(.)          ## Davidson and Harel algorithm
  ) 

# Take out a giant component from the graph
giantGraph <- actorGraph %>% 
  induced.subgraph(., which(comp$membership == which.max(comp$csize)))
vcount(giantGraph) ## the number of nodes/actors/users
ecount(giantGraph) ## the number of edges

# Plot a graph of giant component
giantGraph %>% 
  plot(.,
       layout = layout_with_drl(.),
       edge.arrow.size = .3,
       vertex.label = NA,
       vertex.size = 3,
       vertex.color = 'green',
       vertex.label.cex = .5,
       vertex.label.color = 'black')



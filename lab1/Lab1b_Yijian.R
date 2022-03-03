
library(magrittr)
library(igraph)
library(vosonSML)


load('Lab1a_Yijian.RData')

comp <- components(actorGraph)
giantGraph <- actorGraph %>% 
  induced.subgraph(., which(comp$membership == which.max(comp$csize)))
vcount(giantGraph) 
ecount(giantGraph) 
plot(giantGraph, vertex.size = 7, vertex.label = NA) 

sna_g <- igraph::get.adjacency(giantGraph, sparse=FALSE) %>% network::as.network.matrix()

detach('package:igraph')
library(statnet)

odegScores <- degree(sna_g, cmode = 'outdegree')

centralities <- data.frame('node_name' = as.character(network.vertex.names(sna_g)),
                           'in_degree' = degree(sna_g, cmode = 'indegree'))

centralities$out_degree <- degree(sna_g, cmode = 'outdegree')


centralities$betweenness <- betweenness(sna_g)

centralities$incloseness <- igraph::closeness(giantGraph, mode = 'in')
centralities$outcloseness <- igraph::closeness(giantGraph, mode = 'out')

centralities$eigen <- igraph::eigen_centrality(giantGraph)$vector

centralities$netconstraint <- igraph::constraint(giantGraph)
help(constraint)

centralities$authority <- igraph::authority_score(giantGraph, scale = TRUE)$`vector`
centralities$hub <- igraph::hub_score(giantGraph, scale = TRUE)$`vector`

View(centralities)
detach('package:statnet', unload = TRUE)
library(igraph)

kcore <- giantGraph %>% graph.coreness(.) 

giantGraph %>% 
  plot(.,
       layout = layout_with_gem(.),
       edge.arrow.size = .3,
       vertex.size = 4,
       vertex.label = NA,
       vertex.color = adjustcolor(graph.coreness(.), alpha.f = .3),
       vertex.label.cex = .5,
       vertex.label.color = 'black',
       mark.groups = by(seq_along(graph.coreness(.)), graph.coreness(.), invisible),
       mark.shape = 1/4,
       mark.col = rainbow(length(unique(graph.coreness(.))),alpha = .1),
       mark.border = NA
  )

cluster <- giantGraph %>% cluster_walktrap() 
cluster

modularity(cluster)

membership(cluster)

length(cluster)
sizes(cluster) 
cluster %>% plot(.,giantGraph,
                 layout = layout_with_fr(giantGraph),
                 edge.arrow.size = .3,
                 vertex.size = 4,
                 vertex.label = NA,
                 vertex.color = adjustcolor(membership(.), alpha.f = .3),
                 vertex.label.cex = .5,
                 vertex.label.color = 'black',
                 mark.groups = by(seq_along(membership(.)), membership(.), invisible),
                 mark.shape = 1/4,
                 mark.col = rainbow(length(.),alpha = .1),
                 mark.border = NA
)
giantGraph %>% degree.distribution(.,mode="in") %>% 
  plot(., col = 'black', pch = 19, cex = 1.5,
       main = 'In-degree Distribution',
       ylab = 'Density',
       xlab = 'In-degree')

giantGraph %>% 
  degree.distribution(.,cumulative = TRUE,mode ='in') %>% 
  plot(1:(max(degree(giantGraph,mode='in'))+1),., #
       log='xy', type = 'l',
       main = 'Log-Log Plot of In-degree',
       ylab = 'CCDF',
       xlab = 'In-degree')
in_power <- giantGraph %>% 
  degree.distribution(., mode='in') %>%
  power.law.fit(.)
in_power

giantGraph %>% degree.distribution(.,mode="out") %>% 
  plot(., col = 'black', pch = 19, cex = 1.5,
       main = 'Out-degree Distribution',
       ylab = 'Density',
       xlab = 'Out-degree')
giantGraph %>% 
  degree.distribution(.,cumulative = TRUE,mode ='out') %>% 
  plot(1:(max(degree(giantGraph,mode='out'))+1), 
       ., log='xy', type = 'l',
       main = 'Log-Log Plot of Out-degree',
       ylab = 'CCDF',
       xlab = 'Out-degree')

out_power <- giantGraph %>% 
  degree.distribution(., mode='out') %>%
  power.law.fit(.)
out_power

ntrials <- 1000 
cl.rg <- numeric(ntrials) 
apl.rg <- numeric(ntrials) 
for (i in (1:ntrials)) {
  g.rg <- rewire(giantGraph, keeping_degseq(niter = 100))
  cl.rg[i] <- transitivity(g.rg, type = 'average')
  apl.rg[i] <- average.path.length(g.rg)
}

hist(cl.rg,
     main = 'Histogram of Clustering Coefficient',
     xlab = 'Clustering Coefficient')
par(xpd = FALSE)

abline(v = giantGraph %>% transitivity(., type = 'average'), col = 'red', lty = 2)

t.test(cl.rg, mu=giantGraph %>% transitivity(., type = 'average'),
       alternative = 'greater') 

hist(apl.rg,
     main = 'Histogram of Average Path Length',
     xlab = 'Average Path Length')

abline(v = giantGraph %>% average.path.length(), col = 'red', lty = 2)

t.test(apl.rg, mu=giantGraph %>% average.path.length(.),
       alternative = 'greater')

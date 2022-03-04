library(statnet)

adviceEdgelist <- read.csv("adviceEdgelist.csv")
head(adviceEdgelist)

advice <- as.network.matrix(adviceEdgelist, matrix.type = "edgelist") 
advice 


set.vertex.attribute(advice, "department",read.csv("departmentNode.csv",stringsAsFactors=FALSE)$department) 
set.vertex.attribute(advice, "leader",read.csv("leaderNode.csv")$leader) 
set.vertex.attribute(advice, "tenure",read.csv("tenureNode.csv")$tenure) 
set.vertex.attribute(advice, "office",read.csv("officeNode.csv")$office) 
set.vertex.attribute(advice, "female",read.csv("femaleNode.csv")$female)
advice 


get.vertex.attribute(advice,"department")
get.vertex.attribute(advice,"leader")
get.vertex.attribute(advice,"tenure")
get.vertex.attribute(advice,"office")
get.vertex.attribute(advice,"female")


messageEdgelist <- read.csv("messageEdgelist.csv")
head(messageEdgelist) 
messages <- matrix(nrow = 66, ncol = 66) 
for (i in 1:nrow(messageEdgelist)) { 
  messages[messageEdgelist$SenderId[i], messageEdgelist$ReceiverId[i] ] <- as.numeric(messageEdgelist$MessagesSent[i])
}
for (i in 1:66) { 
  messages[i,i] <- as.numeric(0)
}
hundreds_messages <- messages / 100 
                                   

summary(advice, print.adj = FALSE)       
library('igraph') 


igraph_options(vertex.size = 9, vertex.color = 'cadetblue1', 
               edge.color='red', edge.arrow.size=.4, 
               vertex.label = NA)                     

advice_igraph <- graph.adjacency(as.matrix.network(advice)) 
advice_igraph <- set_vertex_attr(advice_igraph,"female",value = read.csv("femaleNode.csv")$female)
net_layout <- layout_with_fr(advice_igraph) 
                                         
plot(advice_igraph, layout=net_layout, edge.color='black', vertex.label = V(advice_igraph))


V(advice_igraph)$color = ifelse (V(advice_igraph)$female == 1, " cadetblue1 ", " grey ")
plot(advice_igraph, layout=net_layout, edge.color='black', vertex.label = V(advice_igraph))


messages_igraph <- graph.adjacency(messages, weighted = TRUE)

plot(messages_igraph, layout=net_layout, edge.color = adjustcolor('blue',alpha=.2), vertex.label = V(advice_igraph), edge.width=log(E(messages_igraph)$weight), edge.arrow.width =.2)


detach(package:igraph) 
library(statnet)
options(ergm.loglik.warn_dyads=FALSE) 


help("ergm-terms",package = "ergm") 

summary(advice ~ edges)                     
summary(advice ~ mutual)                   
summary(advice ~ odegree(0:5))              
                                          
summary(advice ~ idegree(0:65))             
summary(advice ~ gwodegree(log(2),fixed=T)) 
summary(advice ~ gwidegree(log(2),fixed=T)) 
summary(advice ~ desp(1:5))                 
summary(advice ~ dgwesp(log(2),fixed = T))  


summary(advice ~ nodeicov("office"))             
summary(advice ~ nodeocov("office"))            
summary(advice ~ nodematch("female"))            
summary(advice ~ nodematch("department"))        
summary(advice ~ nodemix("leader",levels2=NULL)) 
summary(advice ~ diff("tenure"))                 
summary(advice ~ edgecov(hundreds_messages))    
                                               

model1 <- ergm(advice ~ edges   
               + mutual                      
               + edgecov(hundreds_messages)   
               + nodemix("leader",base = 3)
               , constraints =~ bd(maxout=5)
) 
summary(model1) 

exp(-2.73064)


model2 <- ergm(advice ~  
               
               mutual
               + gwidegree(log(2), fixed = T)                 
               + gwodegree(2, fixed = T, cutoff = 5)              
               + dgwesp(log(2), type = "OTP", fixed = T, cutoff =5)   
               + nodematch("female")                                   
               + nodemix("leader", base = 3)                            
               + nodematch("department") 
               + nodeicov("office")                                    
               + nodeocov("office")                                   
               + diff("tenure")                                        
               + edgecov(hundreds_messages)   
               , constraints =~ bd(maxout=5)                           
               , control = control.ergm(MCMC.effectiveSize = 50)
) 
summary(model2) 
exp(0.41)
library(texreg)
screenreg(list("model1"=model1,"model2"=model2))

save.image("Lab2_files.RData")


pdf('model1diagnostics.pdf')         
mcmc.diagnostics(model1) 
dev.off()                

pdf('model2diagnostics.pdf')
mcmc.diagnostics(model2) 
dev.off()               

sim1 <- simulate(model1, burnin=100000, interval=100000, nsim=100, verbose=T)  
sim1_net1 <- igraph::graph.adjacency(as.matrix.network(sim1[[1]]))
igraph::plot.igraph(sim1_net1,layout=net_layout,edge.color="brown",  
                    vertex.color = 'grey',edge.arrow.size=.4)                
sim1_net10 <- igraph::graph.adjacency(as.matrix.network(sim1[[10]]))
igraph::plot.igraph(sim1_net10,layout=net_layout,edge.color="red",  
                    vertex.color = 'grey',edge.arrow.size=.4)

sim2 <- simulate(model2, burnin=100000, interval=100000, nsim=100, verbose=T)  
sim2_net1 <- igraph::graph.adjacency(as.matrix.network(sim2[[1]]))
igraph::plot.igraph(sim2_net1,layout=net_layout,edge.color="grey",  
                    vertex.color = 'grey',edge.arrow.size=.4)         
sim2_net10 <- igraph::graph.adjacency(as.matrix.network(sim2[[10]]))
igraph::plot.igraph(sim2_net10,layout=net_layout,edge.color="purple",  
                    vertex.color = 'grey',edge.arrow.size=.4)


model1.tridist <- sapply(1:100, function(x) summary(sim1[[x]] ~triangle)) 
hist(model1.tridist,xlim=c(0,1000),breaks=10)                             
advice.tri <- summary(advice ~ triangle)                                   
advice.tri
arrows(advice.tri,20, advice.tri, 0.5, col="red", lwd=3)                      
c(obs=advice.tri,mean=mean(model1.tridist),sd=sd(model1.tridist),
  tstat=abs(mean(model1.tridist)-advice.tri)/sd(model1.tridist))


model2.tridist <- sapply(1:100, function(x) summary(sim2[[x]] ~triangle)) 
hist(model2.tridist,xlim=c(0,1000),breaks=10)                             
arrows(advice.tri,20, advice.tri, 0.5, col="red", lwd=3)                   
c(obs=advice.tri,mean=mean(model2.tridist),sd=sd(model2.tridist),
  tstat=abs(mean(model2.tridist)-advice.tri)/sd(model2.tridist))

gof1 <- gof(model1, verbose=T, burnin=1e+5, interval=1e+5, control = control.gof.ergm(nsim = 200))

dev.off()           
plot(gof1)         
                   
gof1                

gof2 <- gof(model2, verbose=T, burnin=1e+5, interval=1e+5, control = control.gof.ergm(nsim = 200))
dev.off()
plot(gof2)
gof2

library("igraph")
library("network")
library("sna")
library("ndtv")

nodes<-read.csv("./figures/hypotheses_as_graphs/hawaii_vertices.csv")

stepedges<-read.csv("./figures/hypotheses_as_graphs/step_edges.csv",header=F)
colrs<-c("black","gray50")
stepedgesNS<-read.csv("./figures/hypotheses_as_graphs/step_edges_NtoS.csv",header=F)



#island model make all possible connections

edges<-t(combn(x=nodes$vertex,m=2))
# and the reverse
edges2<-rbind(edges,edges[,c(2,1)])

#make the n-island model
island<-graph_from_data_frame(edges2,nodes,directed=T)


hi_islands<-as.matrix(nodes[,c("long","lat")])




#make a stepping-stone model
stepstone <- graph_from_data_frame(d=stepedges,vertices=nodes,directed=T)

divergence <- graph_from_data_frame(d=stepedgesNS,vertices = nodes, directed=T)

divergence2 <- graph_from_data_frame(d=stepedges,vertices = nodes, directed=T)

#panmixia
plot(stepstone,layout=hi_islands, edge.arrow.size=0.3, 
    edge.curved=0.4, vertex.color=V(stepstone)$sampled, 
    edge.color=NA,vertex.size=3, vertex.label = NA,
#    vertex.label=V(stepstone)$label, vertex.label.cex=0.7, vertex.label.color = "white"
vertex.frame.color =NA)

#stepstone
plot(stepstone,layout=hi_islands, edge.arrow.size=0.3, 
     edge.curved=0.4, vertex.color=V(stepstone)$sampled, 
     edge.color="grey50",vertex.size=3, vertex.label = NA,
     #    vertex.label=V(stepstone)$label, vertex.label.cex=0.7, vertex.label.color = "white"
     vertex.frame.color =NA)

#divergence2
plot(divergence,layout=hi_islands, edge.arrow.size=0.6, edge.width = 3, edge.lty = 1,
     edge.curved=1, vertex.color=V(stepstone)$sampled, 
     edge.color="grey50",vertex.size=3, vertex.label = NA,
     #    vertex.label=V(stepstone)$label, vertex.label.cex=0.7, vertex.label.color = "white"
     vertex.border = V(stepstone)$sector)
#divergence2
plot(divergence2,layout=hi_islands, edge.arrow.size=0.6, edge.width = 3, edge.lty = 1,
     edge.curved=1, vertex.color=V(stepstone)$sampled, 
     edge.color="grey50",vertex.size=3, vertex.label = NA,
     #    vertex.label=V(stepstone)$label, vertex.label.cex=0.7, vertex.label.color = "white"
     vertex.border = V(stepstone)$sector)

#island
plot(island,layout=hi_islands, edge.arrow.size=0.3, edge.curved=0.4, 
     vertex.color=V(island)$sampled, edge.color = "grey90",vertex.frame.color = NA,
     vertex.size=3,vertex.label=NA,vertex.shape="csquare")


plot(stepstone,layout=hi_islands, edge.arrow.size=0.3, edge.curved=0.3, vertex.label=V(stepstone)$label, vertex.color=V(stepstone)$sector, edge.color="gray",vertex.label.cex=0.7,vertex.size=8, vertex.shape="csquare")


# make simulated graphs
stepedges<-stepedges[1:18,]
nodes<-nodes[1:10,]

edges<-t(combn(x=nodes$vertex,m=2))
# and the reverse
edges2<-rbind(edges,edges[,c(2,1)])

island<-graph.data.frame(edges2,nodes,directed=T)


oneD_lattice<-cbind(x=rep(1,10), y=seq(1,10,1))

island<-graph.data.frame(edges2,nodes,directed=T)

plot(island,layout=oneD_lattice, edge.arrow.size=0.05, edge.curved=0.3, vertex.color="grey", vertex.label=NA,vertex.size=8,vertex.label.cex=0.7,vertex.shape="csquare")


stepstone<-graph.data.frame(stepedges,nodes,directed=T)

plot(stepstone,layout=oneD_lattice, edge.arrow.size=0.05, edge.curved=0.3, vertex.label=NA, vertex.color="grey", edge.color="black",vertex.label.cex=0.7,vertex.size=8)

plot(stepstone,layout=oneD_lattice, edge.arrow.size=0.05, edge.curved=0.3, vertex.label=NA, vertex.color="grey", edge.color="grey",vertex.label.cex=0.7,vertex.size=8, vertex.shape="csquare")

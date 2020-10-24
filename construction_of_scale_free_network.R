#### Required packages
list.of.packages <- c("data.table","igraph")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

#### Scale-free Network ####
set.seed(2020)
N<-10000
m<-2
g<-barabasi.game(N, power = 1, m = m, out.dist = NULL, out.seq = NULL, 
              out.pref = FALSE, zero.appeal = 0, directed = FALSE) 

V(g)$name<-seq(1:N)
edges<-get.edgelist(g, names=TRUE)
edges.all.temp.undir<-data.frame('node1'=edges[,1],'node2'=edges[,2], stringsAsFactors = FALSE)
edges.all.temp.undir <- data.table(edges.all.temp.undir)

edges.all.temp<-rbind(edges,edges[,c(2,1)])
edges.all.temp<-data.frame('node1'=edges.all.temp[,1],'node2'=edges.all.temp[,2], stringsAsFactors = FALSE)
edges.all.temp <- data.table(edges.all.temp)

edges.all.final<-vector(mode='list',length=100)
for(i in 1:100){
  edges.all.final[[i]]<-edges.all.temp
  edges.all.final[[i]]$time<-i
}

edges.all.undir.final<-vector(mode='list',length=100)
for(i in 1:100){
  edges.all.undir.final[[i]]<-edges.all.temp.undir
  edges.all.undir.final[[i]]$time<-i
}
edges.all.undir.final[101]<-N

file.name<-paste("contact_network.RData")
save(edges.all.undir.final, file=file.name)

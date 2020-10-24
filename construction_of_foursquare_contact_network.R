#### Required packages
list.of.packages <- c("data.table","matrixStats","dplyr","reshape2","tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

#### Foursquare Contact Network
tky<-read.delim("foursquare_raw_data.txt", header = TRUE, sep = "\t", dec = ".")
tky1 = transform(tky, date= colsplit(UTC_time, pattern = " ", names = c('1', 'mon', "day", "time", "2", "year")))
str(tky)
tky2<-tky1[,c(-7,-8,-9,-13)]
tky2$date.mon<-factor(tky2$date.mon, c("Jan", "Feb", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), c(1:2,4:12))

temp1<-paste(tky2$date.mon, "/", tky2$date.day, "/", tky2$date.year, " ", tky2$date.time, sep="")
temp2<-as.POSIXct(temp1, format="%m/%d/%Y %H:%M:%S")
tky3<-tky2
tky3$date<-temp2
tky3<-tky3[,c(1,2,5,6,11)]
tky3$date.num<-as.numeric(tky3$date)

tky4 <- separate(data = tky3, col = date, into = c("date", "time"), sep = " ")
colnames(tky4)<-c('User_ID','Venue_ID','Latitude','Longitude','date','time','date_num')

venue<-tky4[,c("Venue_ID","Latitude","Longitude")]
venue<-unique(venue)

venue_date<-tky4[,c("Venue_ID","date")]
venue_date<-unique(venue_date)
venue_date<-venue_date[as.integer(as.Date(venue_date$date))<min(as.integer(as.Date(venue_date$date)))+110,]

system.time({
tky5<-vector(mode='list',length=nrow(venue_date))
for(i in 1:nrow(venue_date)){
  sel.1<-venue_date$Venue_ID[i]
  sel.2<-venue_date$date[i]
  tky5[[i]]<-tky4[tky4$Venue_ID == sel.1 & tky4$date == sel.2,]
}
})


edge<-c()
for(i in 1:nrow(venue_date)){
  temp.1<-tky5[[i]]
  temp.2<-unique(temp.1$User_ID)
  temp.3<-if(length(temp.2)>=2) t(combn(temp.2,2)) else next
  temp.edge<-as.data.table(temp.3)
  temp.edge$date<-temp.1$date[1]
  edge<-rbindlist(list(edge,temp.edge))
}
edge.temp<-edge
edge<-unique(edge.temp)
names(edge)<-c('node1','node2','time')
edge[, time:=as.integer(as.Date(time))]
edge[, time:=as.integer(time-min(time)+1)]
#save(edge, file="tokyo_edgelist_time.RData")


#setwd("~/Dropbox/1. PENN - Doctorate/Research/COVID-19/Covid19-2 Non Shared/Impact of contact tracing/Updates/09-2020/09-02-20/Data/cluster 0m")
#edge<-get(load("tokyo_edgelist_time.RData"))
#tky4<-get(load("Processed_dataset_TSMC2014_TKY.RData"))

edges.all.temp.undir<-edge

temp<-as.integer(as.Date(sort(unique(tky4$date))))-min(as.integer(as.Date(sort(unique(tky4$date)))))+1
edges.all.temp.undir<-edges.all.temp.undir[time %in% temp[1:100]]
nodes<-sort(unique(c(edges.all.temp.undir$node1,edges.all.temp.undir$node2)))
nodes.temp<-data.table(node.old = nodes, node.new = 1:length(nodes))
time.temp<-data.table(time.old = temp[1:100], time.new = 1:100)

cust.1<-as.integer(nodes)
cust.2<-!is.na(cust.1)
customer<-nodes.temp[['node.new']][cust.2]

edges.all.temp.undir<-merge(edges.all.temp.undir, nodes.temp, by.x = "node2", by.y = "node.old", all.x=TRUE)
names(edges.all.temp.undir)[4]<-c('node2.new')
edges.all.temp.undir<-merge(edges.all.temp.undir, nodes.temp, by.x = "node1", by.y = "node.old", all.x=TRUE)
names(edges.all.temp.undir)[5]<-c('node1.new')
edges.all.temp.undir<-merge(edges.all.temp.undir, time.temp, by.x = "time", by.y = "time.old", all.x=TRUE)
edges.all.temp.undir<-edges.all.temp.undir[,c(5,4,6)]
names(edges.all.temp.undir)<-c('node1','node2','time')


edges.all.temp.undir.rev<-copy(edges.all.temp.undir[,c(2,1,3)])
names(edges.all.temp.undir.rev)<-c('node1','node2','time')
edges.all.temp<-rbindlist(list(edges.all.temp.undir,edges.all.temp.undir.rev))

edges.all.final<-vector(mode='list',length=max(edges.all.temp$time))
for(i in 1:max(edges.all.temp$time)){
  edges.all.final[[i]]<-edges.all.temp[edges.all.temp[[3]]==i,]
}

edges.all.undir.final<-vector(mode='list',length=max(edges.all.temp.undir$time))
for(i in 1:max(edges.all.temp.undir$time)){
  edges.all.undir.final[[i]]<-edges.all.temp.undir[edges.all.temp.undir[[3]]==i,]
}
edges.all.undir.final[max(edges.all.temp.undir$time)+1]<-length(unique(c(do.call("rbind", edges.all.undir.final)$node1,do.call("rbind", edges.all.undir.final)$node2)))


file.name<-paste("contact_network.RData")
save(edges.all.undir.final, file=file.name)


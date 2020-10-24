#### Required packages
list.of.packages <- c("data.table","matrixStats")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

#### University Student Contact Network
school.temp <- read.csv("university_student_raw_data.csv", header=FALSE, comment.char="#")
school<-school.temp
school<-school[!school$V3 %in% c(-1,-2),]
school<-unique(school)
school$day<-((school$V1) %/% (3600*24))+1

colnames(school)<-c('time','node1','node2','dB','day')
school<-as.data.table(school[order(school$day,school$node1,school$node2,school$time),])
school[ , diff := time - shift(time), by = list(day,node1,node2)] 
school[ diff!=300, diff := NA]

school[ , group := cumsum(is.na(diff))]
school[!is.na(diff) , duration := cumsum(diff), by = group]

school.final<-unique(school[duration==900, c('node1','node2','day')])

node.name.old<-sort(unique(c(school.final[[1]],school.final[[2]])))
node.name.new<-1:length(node.name.old)
node.name<-data.table('node1'=node.name.old, 'node_name_new'=node.name.new, stringsAsFactors = F)
school.final<-merge(x = school.final, y = node.name, by = "node1", all.x = TRUE)
node.name<-data.table('node2'=node.name.old, 'node_name_new'=node.name.new, stringsAsFactors = F)
school.final<-merge(x = school.final, y = node.name, by = "node2", all.x = TRUE)
school.final<-school.final[,c(4,5,3)]
names(school.final)<-c('node1','node2','time')

edges.all.temp.undir<-copy(school.final[order(time,node1,node2)])

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



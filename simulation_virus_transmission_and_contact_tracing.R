
## Required packages
list.of.packages <- c("data.table","matrixStats", "reshape2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)


## Import contact network
file.name<-paste("contact_network.RData")
edges.all.undir.final<-get(load(file.name))


## Setting parameters values
cooperativity<-0.2  # cooperativity
N<-edges.all.undir.final[[length(edges.all.undir.final)]] # Total number of nodes in erdos-renyi and scale-free networks
T.max = length(edges.all.undir.final)-1 # total number of days
edges.all.undir.final[T.max+1]<-NULL
retest.delay<-3  # those who test negative can be tested again after 3 days from the test date 

prob.inf.from.symp<-0.4 # Probability with which a symptomatic individual infects a susceptible in an interaction
prob.inf.from.presymp<-prob.inf.from.symp # Probability with which a presymptomatic individual infects a susceptible in an interaction
infness.of.asymp<-0.75 # infectiousness of asymptomatic individuals relative to symptomatic
prob.inf.from.asymp<-prob.inf.from.symp*infness.of.asymp # Probability with which a asymptomatic individual infects a susceptible in an interaction
prop.asymp<-0.4 # proportion of asymptomatic individuals among infected people.
  
I0<-1 # the number of initially infected individuals
symptom<-c('Is','RT') # compartments with symptoms
nosymptom<-c('S', 'Ia-L', 'Ip-L', 'Ip','Ia','R') # compartments without symptoms
strTmp<-as.character(1:T.max)
emtpy.table.cha<-setNames(data.table(matrix('NA' ,nrow = N, ncol = T.max)), strTmp)
emtpy.table.int<-setNames(data.table(matrix(Inf ,nrow = N, ncol = T.max)), strTmp)


## This user-defined function epidemic implements simulations of disease outbreaks and k-hop contact tracing.
epidemic<-function(iteration,hhh){ 
  hops<-hhh
  
  wdata <- copy(emtpy.table.cha)
  wdata.flu <- copy(emtpy.table.cha)
  wdata.t <- copy(emtpy.table.cha)
  wdata.tp <- rep(NA,N)
  wdata.avl <- rep(Inf,N)

  ### initial condition
  I0.node<-base::sample.int(N, I0)
  presymptomatic <- symptomatic <- asymptomatic <- rep(0,N);
  symptomatic[I0.node[1]]<-1
  presymptomatic[I0.node[2]]<-1
  asymptomatic[I0.node[3]]<-1
  susceptible = 1-(presymptomatic+symptomatic+asymptomatic)
  ready.to.test.Is = rep(0,N)
  recovered = rep(0,N)
  dead = rep(0,N)
  
  set(wdata, i=which(susceptible==1), j=1L, value='S')
  if(sum(presymptomatic==1)>0){set(wdata, i=which(presymptomatic==1), j=1L, value='Ip')}
  if(sum(symptomatic==1)>0){set(wdata, i=which(symptomatic==1), j=1L, value='Is')}
  if(sum(asymptomatic==1)>0){set(wdata, i=which(asymptomatic==1), j=1L, value='Ia')}
  
  edges.all.final.ct<-vector(mode='list',length=T.max)
  edges.all.final.transmission<-vector(mode='list',length=T.max)
    

  ### COVID-19 compartmental model with parameters listed in Table 1
  {
      inf.before.symp.onset<-1
      f.Ip_pre.to.Ip<-function(Ip_pre.list,t){
        sel<-rbinom(length(Ip_pre.list),1,1./inf.before.symp.onset)
        status<-ifelse(sel==1,'Ip','Ip-L')
        return(status)
      }
      
      inf.before.rec<-1
      f.Ia_pre.to.Ia<-function(Ia_pre.list,t){
        sel<-rbinom(length(Ia_pre.list),1,1./inf.before.rec)
        status<-ifelse(sel==1,'Ia','Ia-L')
        return(status)
      }
      
      incub.per<-5-inf.before.symp.onset
      f.Ip.to.Is<-function(Ip.list,t){
        sel<-rbinom(length(Ip.list),1,1./incub.per)
        status<-ifelse(sel==1,'Is','Ip')
        return(status)
      }
      
      wait.see.per<-4
      f.Is.to.RTs<-function(Is.list,t){
        sel<-rbinom(length(Is.list),1,1./wait.see.per)
        status<-ifelse(sel==1,'RT','Is')
        return(status)
      }
      
      rec.die.per<-14-wait.see.per
      death.rate<-0.0065
      f.RTs.to.D<-function(RT.list,t){
        status<-base::sample(c('R', 'D', 'RT'), size = length(RT.list), replace = TRUE, prob = c((1./rec.die.per)*(1.-death.rate), (1./rec.die.per)*death.rate, (1.-1./rec.die.per)))
        return(status)
      }
      
      rec.Ia.per<-8-inf.before.rec
      f.Ia.to.R<-function(Ia.list,t){
        sel<-rbinom(length(Ia.list),1,1./rec.Ia.per)
        status<-ifelse(sel==1,'R','Ia')
        return(status)
      }
  }

  test<-c()
  TP.list<-c();
  noTP.list<-c();
  Tested<-c(); temp<-c();
  ct.time<-1
  for(t in 2:T.max){
    ## Contact network on day t-1
    temp.rm<-c(); temp.edge<-c(); temp.edge.rev<-c();
    temp.rm<-unique(c(which(wdata.tp %in% c('TP','D')),test))
    temp.edge<-edges.all.undir.final[[t-1]][!(node1 %in% temp.rm | node2 %in% temp.rm)][,c(1,2)]
    temp.edge.rev<-copy(temp.edge[,c(2,1)])
    names(temp.edge.rev)<-c('node1','node2')
    edges.all.final.transmission[[t-1]]<-rbindlist(list(temp.edge,temp.edge.rev))  
    
    ## Epidemic process; status update from day t-1 to day t
    { 
      S.list<-c()
      I.list<-c()
      temp<-which(wdata[[t-1]] %in% c('S'))
      if(length(temp)>0){
        set(wdata, i=temp, j=as.integer(t), value='S')
        S.list<-data.table(node1=as.integer(temp), status='S')
      }
      
      temp<-which(wdata[[t-1]] %in% c('Ip','Is','Ia'))
      I.list<-wdata[[t-1]][temp]
      I.list<-data.table(node2=as.integer(temp), status=I.list)
      
      inf.list<-c()
      temp.neighbor.of.I<-c()
      neighbor.of.I<-c()
      edges.all.ct<-c()
      
      temp.neighbor.of.I<-edges.all.final.transmission[[t-1]][node1 %in% S.list$node1 & node2 %in% I.list$node2]
      neighbor.of.I<-temp.neighbor.of.I
      
      if(nrow(neighbor.of.I)>=1){
        neighbor.of.I<-I.list[neighbor.of.I, on="node2"]  
        neighbor.of.I<-neighbor.of.I[,.(count = .N), by=.(node1,status)]
        neighbor.of.I$status<-factor(neighbor.of.I$status, levels=c("Ip", "Is", "Ia"))
        neighbor.of.I<-setDT(dcast(neighbor.of.I, node1 ~ status, fun = sum, value.var = "count", drop=F))
        neighbor.of.I[,prob.inf:=1-((1-prob.inf.from.presymp)^Ip)*((1-prob.inf.from.symp)^Is)*((1-prob.inf.from.asymp)^Ia)]
        neighbor.of.I[,new.inf:=rbinom(length(neighbor.of.I$prob.inf), size = 1, prob=neighbor.of.I$prob.inf)]
        
        temp<-rbinom(sum(neighbor.of.I$new.inf==1), size = 1, prob=prop.asymp)
        neighbor.of.I[new.inf==1,new.inf.Ia:=temp]
        neighbor.of.I$new.inf.Ia[is.na(neighbor.of.I$new.inf.Ia)]<-0
        neighbor.of.I[new.inf==1 & new.inf.Ia==1,stt:='Ia-L']
        neighbor.of.I[new.inf==1 & new.inf.Ia!=1,stt:='Ip-L']
        neighbor.of.I<-neighbor.of.I[!is.na(stt)]
        
        if(nrow(neighbor.of.I)>1){
          set(wdata, i=as.integer(neighbor.of.I$node1), j=as.integer(t), value=(neighbor.of.I$stt))
          temp.inf.list<-c()
          temp.inf.list<-neighbor.of.I[temp.neighbor.of.I, on="node1"]  
          temp.inf.list<-temp.inf.list[I.list, on="node2"] 
          temp.inf.list[, time:=t]
          inf.list<-rbindlist(list(inf.list,temp.inf.list))
        }
      }
      
      Ip_pre.list<-which(wdata[[t-1]]=='Ip-L')
      if(length(Ip_pre.list)>=1){
        temp<-f.Ip_pre.to.Ip(Ip_pre.list,t-1)
        set(wdata, i=as.integer(Ip_pre.list), j=as.integer(t), value=(temp))
      }
      
      Ia_pre.list<-which(wdata[[t-1]]=='Ia-L')
      if(length(Ia_pre.list)>=1){
        temp<-f.Ia_pre.to.Ia(Ia_pre.list,t-1)
        set(wdata, i=as.integer(Ia_pre.list), j=as.integer(t), value=(temp))
      }
      
      Ip.list<-which(wdata[[t-1]]=='Ip')
      if(length(Ip.list)>=1){
        temp<-f.Ip.to.Is(Ip.list,t-1)
        set(wdata, i=as.integer(Ip.list), j=as.integer(t), value=(temp))
      }
      
      Is.list<-which(wdata[[t-1]]=='Is')
      if(length(Is.list)>=1){
        temp<-f.Is.to.RTs(Is.list,t-1)
        set(wdata, i=as.integer(Is.list), j=as.integer(t), value=(temp))
      }
      
      RT.list<-which(wdata[[t-1]]=='RT')
      if(length(RT.list)>=1){
        temp<-f.RTs.to.D(RT.list,t-1)
        set(wdata, i=as.integer(RT.list), j=as.integer(t), value=(temp))
      }
      
      Ia.list<-which(wdata[[t-1]]=='Ia')
      if(length(Ia.list)>=1){
        temp<-f.Ia.to.R(Ia.list,t-1)
        set(wdata, i=as.integer(Ia.list), j=as.integer(t), value=(temp))
      }
      
      R.list<-which(wdata[[t-1]]=='R')
      if(length(R.list)>=1){
        set(wdata, i=as.integer(R.list), j=as.integer(t), value='R')
      }
      
      D.list<-which(wdata[[t-1]]=='D')
      if(length(D.list)>=1){
        set(wdata, i=as.integer(D.list), j=as.integer(t), value='D')
      }
    }

    
    
    ### k-hop contact tracing and testing
    ## tracable network on day t-1; cooperativity is considered.
    temp.rm<-c(); temp.edge.all<-c(); temp.edge<-c(); temp.edge.rev<-c();
    temp.rm<-unique(c(which(wdata.tp %in% c('TP','D')),test))
    if(length(temp.rm)>0){
      temp.edge.all<-edges.all.undir.final[[t-1]][!(node1 %in% temp.rm | node2 %in% temp.rm)][,c(1,2)]
      temp.num<-round(nrow(temp.edge.all)*cooperativity)
    }else{
      temp.edge.all<-edges.all.undir.final[[t-1]][,c(1,2)]
      temp.num<-round(nrow(temp.edge.all)*cooperativity)
    }
    temp.edge<-copy(temp.edge.all[sample(.N,temp.num)])
    temp.edge.rev<-copy(temp.edge[,c(2,1)])
    names(temp.edge.rev)<-c('node1','node2')
    edges.all.final.ct[[t-1]]<-rbindlist(list(temp.edge,temp.edge.rev))
    
    
    test<-c(); 
    edges.all.ct<-c()
    g.tilde.hop<-c()
    g.tilde.hop.temp<-c()
    RT.list<-unique(which(wdata[[t]]=='RT'))
    #RT.list<-unique(c(RT.list))
    
    wdata.avl<-wdata.avl+1
    Test.available<-which(wdata.avl>=retest.delay) # available nodes that can be tested on day t
    
    tt<-ifelse(t-14>1,t-14,1):(t-1)
    edges.all.ct<-unique(rbindlist(edges.all.final.ct[tt])[,1:2])

    if(length(Test.available)==0){
      TP.list<-c();
      noTP.list<-c();
      test<-c();
      next
    }
    else if((length(RT.list)==0 & ct.time==1) & length(Test.available)>0){
      TP.list<-c();
      noTP.list<-c();
      test<-c();
      next
    }
    else if((length(RT.list)>=1 & ct.time==1) & length(Test.available)>0){
      Tested<-test<-RT.list
      TP.list<-RT.list
      set(wdata.t, i=RT.list, j=as.integer(t), value='Test.RT')
      ct.time<-ct.time+1
    }
    else{
      if(hops>0){
        if(length(TP.list)>=1){
          g.tilde.hop.temp<-c(); g.tilde.hop<-c();
          for(i in 1:hops){
            temp<-if(i==1) TP.list else g.tilde.hop.temp
            g.tilde.hop.temp<-unique(edges.all.ct[edges.all.ct$node1 %in% temp]$node2)
            g.tilde.hop.temp<-if(length(g.tilde.hop.temp)==1){g.tilde.hop.temp
            }else{base::sample(g.tilde.hop.temp)}
            g.tilde.hop<-unique(c(g.tilde.hop,g.tilde.hop.temp))
          }
          
          test1.temp<-g.tilde.hop[g.tilde.hop %in% Test.available]
          test1.temp.s<-test1.temp[test1.temp %in% unique(c(which(wdata[[t]] %in% symptom)))]
          test1.temp.ns<-test1.temp[test1.temp %in% which(wdata[[t]] %in% nosymptom)]
          
          test1.s<-c(); 
          if(length(test1.temp.s)>0){
            test1.s<-test1.temp.s
            test<-c(test,test1.s)
            set(wdata.t, i=test1.s, j=as.integer(t), value='Test.CT')
            Test.available<-Test.available[!Test.available %in% test]
          }
          
          test1.ns<-c();
          test1.temp.ns<-test1.temp.ns[test1.temp.ns %in% Test.available]
          if(length(test1.temp.ns)>0){
            test1.ns<-test1.temp.ns
            test<-c(test,test1.ns)
            set(wdata.t, i=test1.ns, j=as.integer(t), value='Test.CT')
            Test.available<-Test.available[!Test.available %in% test]
          }
          
          test2.s<-c(); test2.temp.s<-c();
          test2.temp.s<-RT.list
          test2.temp.s<-test2.temp.s[test2.temp.s %in% Test.available]
          if(length(test2.temp.s)>0){
            test2.s<-test2.temp.s
            test<-c(test,test2.s)
            set(wdata.t, i=test2.s, j=as.integer(t), value='Test.RT')
            Test.available<-Test.available[!Test.available %in% test]
          }
        }
        else{
          test<-c(); test2.s<-c(); test2.temp.s<-c();
          test2.temp.s<-RT.list
          test2.temp.s<-test2.temp.s[test2.temp.s %in% Test.available]
          if(length(test2.temp.s)>0){
            test2.s<-test2.temp.s
            test<-c(test,test2.s)
            set(wdata.t, i=test2.s, j=as.integer(t), value='Test.RT')
            Test.available<-Test.available[!Test.available %in% test]
          }
        }
      }
      else{
        test<-c(); test2.s<-c(); test2.temp.s<-c();
        test2.temp.s<-RT.list
        test2.temp.s<-test2.temp.s[test2.temp.s %in% Test.available]
        if(length(test2.temp.s)>0){
          test2.s<-test2.temp.s
          test<-c(test,test2.s)
          set(wdata.t, i=test2.s, j=as.integer(t), value='Test.RT')
          Test.available<-Test.available[!Test.available %in% test]
        }
      }
    }

    
    temp<-which(wdata[[t]] %in% c('Ia-L','Ip-L','Ip','Is','Ia','RT'))
    TP.list<-test[test %in% temp]      # nodes who newly tested positive at time t
    noTP.list<-test[!test %in% temp]   # nodes who newly tested negative at time t
    Tested<-unique(c(Tested, test))     # nodes who have tested by time t
    
    if(length(TP.list)>0){
       wdata.avl[TP.list] <- -Inf
       wdata.tp[TP.list] <- 'TP'
    }
    
    if(length(noTP.list)>0){
      wdata.avl[noTP.list] <- 0 
    }
    
    D.list<-which(wdata[[t]]=='D')
    if(length(D.list)>0){
      wdata.avl[D.list] <- -Inf
      wdata.tp[D.list] <- 'D'
    }
  }
  list(process=wdata, testednodes=wdata.t)
}


## simulation for k-hop contact tracing for k = 0,1,2,3,4,5 
num_iteration<-1000 # 1000 simulation runs
for(k in 0:5){ 
  process<-vector(mode='list',length=num_iteration)
  tested.nodes<-vector(mode='list',length=num_iteration)
  for(pp in 1:num_iteration){
    set.seed(pp)
    temp<-epidemic(pp,as.integer(k))
    process[[pp]]<- copy(temp$process)
    tested.nodes[[pp]]<- copy(temp$testednodes)
    print(paste(k, '-hop contact tracing, simulation run ',pp, sep = ''))
  }
  file.name.1<-paste("progression_of_disease_over_time_",k,"-hop_contact_tracing.RData",sep="")
  file.name.2<-paste("tests_over_time_",k,"-hop_contact_tracing.RData",sep="")
  save(process, file=file.name.1)
  save(tested.nodes, file=file.name.2)
  print(paste('### The end of simulations for ', k, '-hop contact tracing', sep = ''))
}

######## Analysis #####
time.seq<-1:T.max
max.hops<-5
cumul.result.I.time<-vector(mode='list',length=max.hops+1)
new.result.I.time<-vector(mode='list',length=max.hops+1)
for(hhh in 0:max.hops){
  cumul.result.I.temp<-matrix(, nrow = T.max, ncol = num_iteration)
  new.result.I.temp<-matrix(, nrow = T.max, ncol = num_iteration)
  file.name.1<-paste("progression_of_disease_over_time_",hhh,"-hop_contact_tracing.RData",sep="")
  process.temp<-get(load(file.name.1))
  for(pp in 1:num_iteration){
    temp<-process.temp[[pp]]
    sel.1<-temp == 'Ia-L' | temp == 'Ia' | temp == 'Ip-L' | temp == 'Ip' | temp == 'Is'
    sel.1<-as.matrix(sel.1)
    temp.1<-rowCumsums(sel.1)

    temp.1.1<-temp.1>=1
    temp.1.2<-vector(mode = "integer", length = nrow(temp.1))
    for(i in 1:nrow(temp.1)){
       temp.1.2[i]<-ifelse(sum(temp.1[i,]==1)==0,Inf,min(which(temp.1[i,]==1)))
    }
    cumul.result.I.temp[,pp]<-colSums(temp.1.1)
    new.result.I.temp[,pp]<-table(factor(temp.1.2, levels = time.seq))
  }
  cumul.result.I.time[[hhh+1]]<-cumul.result.I.temp
  new.result.I.time[[hhh+1]]<-new.result.I.temp
}
cumululative.infections.over.time<-matrix(, nrow = 1+max.hops, ncol = T.max)
daily.new.infections<-matrix(, nrow = 1+max.hops, ncol = T.max)
rownames(cumululative.infections.over.time) <- c('0-hop', '1-hop', '2-hop', '3-hop', '4-hop', '5-hop')
colnames(cumululative.infections.over.time) <- 1:T.max
rownames(daily.new.infections) <- c('0-hop', '1-hop', '2-hop', '3-hop', '4-hop', '5-hop')
colnames(daily.new.infections) <- 1:T.max
for(k in 0:max.hops){
  cumululative.infections.over.time[1+k,]<-rowMeans(cumul.result.I.time[[1+k]])*(1000/N)
  daily.new.infections[1+k,]<-rowMeans(new.result.I.time[[1+k]])*(1000/N)
}
write.table(cumululative.infections.over.time, file = "cumulative_infections_per_1000.csv",row.names=TRUE,col.names=NA, sep=",") # export a matrix to a file.
write.table(daily.new.infections, file = "daily_new_infections_per_1000.csv",row.names=TRUE,col.names=NA, sep=",") # export a matrix to a file.


test.1<-vector(mode='list',length=max.hops+1)
test.2<-vector(mode='list',length=max.hops+1)
test.all<-vector(mode='list',length=max.hops+1)
for(hhh in 0:max.hops){
  test.1.temp<-matrix(, nrow = T.max, ncol = num_iteration)
  test.2.temp<-matrix(, nrow = T.max, ncol = num_iteration)
  file.name.1<-paste("tests_over_time_",hhh,"-hop_contact_tracing.RData",sep="")
  tested.nodes<-get(load(file.name.1))
  for(pp in 1:num_iteration){
    temp<-(tested.nodes[[pp]])
    test.1.temp[,pp]<-colSums(temp == 'Test.RT')
    test.2.temp[,pp]<-colSums(temp == 'Test.CT')
  }
  test.1[[hhh+1]]<-(test.1.temp)
  test.2[[hhh+1]]<-(test.2.temp)
  test.all[[hhh+1]]<-(test.1.temp+test.2.temp)
}

test.over.time<-matrix(, nrow = 1+max.hops, ncol = T.max)
rownames(test.over.time) <- c('0-hop', '1-hop', '2-hop', '3-hop', '4-hop', '5-hop')
colnames(test.over.time) <- 1:T.max
for(k in 0:max.hops){
  test.over.time[1+k,]<-rowMeans(test.all[[1+k]])*(1000/N)
}
write.table(test.over.time, file = "test_per_1000.csv",row.names=TRUE,col.names=NA, sep=",") # export a matrix to a file.

for(k in 0:5){ 
  file.name.1<-paste("progression_of_disease_over_time_",k,"-hop_contact_tracing.RData",sep="")
  file.name.2<-paste("tests_over_time_",k,"-hop_contact_tracing.RData",sep="")
  file.remove(file.name.1)
  file.remove(file.name.2)
}

print(paste('The end of the simulations', sep = ''))



###############INPUTS############################
#set working directory (where all imput data is present, including shapefiles)
setwd("C:\\Users\\Sidonie\\Dropbox\\Carvalho_et_al_Neotropical_biome_shifts\\Data_and_analyses\\Analyses\\final_Biogeobears_analyses\\mammals")
##############################
binsize = 5 # binsize
area1 = "a+|b+"
area2 = "c+|d+"
################################################
#check number tre files(should be same as nr RDATA files)

div.plots <-function(binsize = 5, area1 = "a", area2 = "b", wd = getwd()){
  number=length(grep("clade_\\d+_area_\\d\\_updated\\.tre",dir(),perl =TRUE))
  fivemilgaps=data.frame()
  
  #load in inputs
  agegap = binsize # binsize
  
  #Load dependancies
  require(BioGeoBEARS)
  require(geiger)
  
  for (o in 1:number){
    #i = 33
    a = "clade_"
    c = "_area_1.geog"
    d = "_area_1_updated.tre"
    e = "area_1_DEC"
    f = ".Rdata"
    g = ".pdf"
    
    #generate in/output names
    inputfile = paste(a,o,e,f,sep="")
    treeobjectname=paste(a,o,d,sep="")
    geogobjectname=paste(a,o,c,sep="")
    #rdataobject=paste(a,o,e,f,sep="")
    #pdfname=paste(a,o,e,g,sep="")
    
    ##########################################
    
    #loop over the number of Rdata files
    load(inputfile) #load biogeobears run
    tr=read.tree(res$input$trfn) #read tree
    branchtimes = branching.times(tr) #Get the nodeages
    maxfivemilcatogories = ceiling(max(branchtimes)/agegap) #Get total amount of bin's in the tree
    
    #Get tipranges & total available states in the dataset
    geogfn=res$input$geogfn
    tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(res$input$geogfn))
    max_range_size=res$inputs$max_range_size
    areas = getareas_from_tipranges_object(tipranges)
    statenames = areas_list_to_states_list_new(areas, maxareas=max_range_size, include_null_range=TRUE, split_ABC=FALSE)
    
    #create table containing the areas for all nodes
    areapernode=matrix(ncol=7, nrow=length(res$ML_marginal_prob_each_state_at_branch_top_AT_node[,1]), NA)
    
    #Retrieve culumn with the maximum nodeage
    for (i in 1:length(res$ML_marginal_prob_each_state_at_branch_top_AT_node[,1])){
      areapernode[i,1] = which.max(res$ML_marginal_prob_each_state_at_branch_top_AT_node[i,])
    }
    #add the age of each node
    areapernode[,2] = 0
    for (i in 1:length(names(branchtimes))){
      num = noquote(names(branchtimes)[i])
      num <- as.numeric(as.character(num))
      areapernode[num,2] = branchtimes[names(branchtimes)[i]]
    }
    
    #get area names translation
    areas = getareas_from_tipranges_object(tipranges)
    statenames = areas_list_to_states_list_new(areas, maxareas=max_range_size, include_null_range=TRUE, split_ABC=FALSE)
    
    #change the states with the statenames
    for (i in 1:length(areapernode[,1])){
      areapernode[i,3] = statenames[[as.numeric(as.character(noquote(areapernode[i,1])))]]
    }
    
    #Check wether nodes are present in the old or new world
    for (i in 1:length(areapernode[,1])){
      areapernode[i,4] <- grepl(area1,areapernode[i,3], perl=TRUE)
      areapernode[i,5] <- grepl(area2,areapernode[i,3], perl=TRUE)
    }
    
    #note if nodes are present in area1, area2 or both
    areapernode[,6] = 0
    for (i in 1:length(areapernode[,1])){
      if (areapernode[i,4] == TRUE & areapernode[i,5] == TRUE){
        areapernode[i,6] = 3
      }
      if(areapernode[i,4] == TRUE & areapernode[i,5] == FALSE){
        areapernode[i,6] = 1
      }
      if(areapernode[i,4] == FALSE & areapernode[i,5] == TRUE){
        areapernode[i,6] = 2
      }
    }
    
    
    #create a table containing the state with the highest likelihood for each branchcorner
    areaperbranch=matrix(ncol=7, nrow=length(res$ML_marginal_prob_each_state_at_branch_bottom_below_node[,1]), NA)
    for (i in 1:length(res$ML_marginal_prob_each_state_at_branch_bottom_below_node[,1])){
      if (i != length(tr$tip.label)+1){
        areaperbranch[i,1] = which.max(res$ML_marginal_prob_each_state_at_branch_bottom_below_node[i,])
      }
    }
    
    #add the age of each branchcorner
    areaperbranch[,2] = 0
    for (i in 1:length(names(branchtimes))){
      numm = noquote(names(branchtimes)[i])
      numm <- as.numeric(as.character(numm))
      areaperbranch[numm,2] = branchtimes[names(branchtimes)[i]]
    }
    
    #change the states with the statenames
    for (i in 1:length(areaperbranch[,1])){
      if (i != length(tr$tip.label)+1){
        areaperbranch[i,3] = statenames[[as.numeric(as.character(noquote(areaperbranch[i,1])))]]
      }
    }
    
    #Check wether branch_corner are present in the area1 or area2
    for (i in 1:length(areaperbranch[,1])){
      areaperbranch[i,4] <- grepl(area1,areaperbranch[i,3], perl=TRUE)
      areaperbranch[i,5] <- grepl(area2,areaperbranch[i,3], perl=TRUE)
    }
    
    #note if branchs corners are present in the old world, new world or both
    areaperbranch[,6] = 0
    for (i in 1:length(areaperbranch[,1])){
      if (areaperbranch[i,4] == TRUE & areaperbranch[i,5] == TRUE){
        areaperbranch[i,6] = 3
      }
      if(areaperbranch[i,4] == TRUE & areaperbranch[i,5] == FALSE){
        areaperbranch[i,6] = 1
      }
      if(areaperbranch[i,4] == FALSE & areaperbranch[i,5] == TRUE){
        areaperbranch[i,6] = 2
      }
    }
    
    #note nodes that have crossed the atlantic
    areapernode[,7] = 0
    for (i in 1:length(tr$edge[,1])){
      if (((areapernode[tr$edge[i,1],6] != areaperbranch[tr$edge[i,2],6])&(areapernode[tr$edge[i,1],6] != 0)&(areaperbranch[tr$edge[i,2],6] != 0)&
             (areapernode[tr$edge[i,1],6] != 3))){
        areapernode[tr$edge[i,1],7] = 1
      }
    }
    
    #note branches on which transatlantic change occured
    areaperbranch[,7] = FALSE
    for (i in 1:length(tr$edge[,1])){
      if ((areaperbranch[tr$edge[i,2],6] != areapernode[tr$edge[i,2],6])&(areaperbranch[tr$edge[i,2],6] != 0)&(areapernode[tr$edge[i,2],6] != 0)&
            (areaperbranch[tr$edge[i,2],6] != 3)){
        if (areaperbranch[tr$edge[i,1],7] == FALSE){
          areaperbranch[tr$edge[i,1],7] = TRUE
          
        }
      }
    }
    
    #Create matrix with nodes that changed
    nodetimes = matrix(ncol=2, nrow=(length(grep(1,areapernode[,7], perl=TRUE))+(length(grep(2,areapernode[,7], perl=TRUE)))), 0)
    x = 1
    for (i in 1:length(areapernode[,1])){
      if (areapernode[i,7] != 0){
        nodetimes[x,1] = as.numeric(as.character(noquote(areapernode[i,2])))
        nodetimes[x,2] = as.numeric(as.character(noquote(areapernode[i,7])))
        x=x+1
      }
    }
    
    
    #create 0-maxage matrix for each X mil period
    fivemilgaps=matrix(ncol=9, nrow=ceiling(max(branchtimes)/agegap), 0)
    for (i in 1:ceiling(max(branchtimes)/agegap)){
      fivemilgaps[i,1]=(agegap*i)
    }
    
    if(max(areapernode[,7])!=0){
      #Add a column with the number of nodes with rangeshifts every X mil years
      for (i in 1:length(nodetimes[,1])){
        for(y in 1:length(fivemilgaps[,1])){
          if ((nodetimes[i,1] <= fivemilgaps[y,1])&(nodetimes[i,1] >= (fivemilgaps[y,1]-agegap))){
            fivemilgaps[y,3] = fivemilgaps[y,3]+1
          }
        }
      }
      
      #Add a column with the number of rangeshifts at nodes (can be up to two) every X mil years
      for (i in 1:length(nodetimes[,1])){
        for(y in 1:length(fivemilgaps[,1])){
          if ((nodetimes[i,1] <= fivemilgaps[y,1])&(nodetimes[i,1] >= (fivemilgaps[y,1]-agegap))){
            fivemilgaps[y,7] = fivemilgaps[y,7]+nodetimes[i,2]
          }
        }
      }
    }
    
    #Add a column with the number of nodes every X mil years
    for (i in 1:length(branchtimes)){
      for(y in 1:length(fivemilgaps[,1])){
        if ((branchtimes[[i]] <= fivemilgaps[y,1])&(branchtimes[[i]] >= (fivemilgaps[y,1]-agegap))){
          fivemilgaps[y,2] = fivemilgaps[y,2]+1
        }
      }
    }
    
    
    
    #Create table with branch ages
    branchevents=matrix(ncol=3, nrow=length(tr$edge[,1]), 0)
    for (i in 1:length(tr$edge[,1])){
      if(tr$edge[i,2] > length(tr$tip.label)){
        branchevents[i,1]=branchtimes[[as.character(tr$edge[i,2])]]
      }
      else{
        branchevents[i,1]=0
      }
      branchevents[i,2]=branchtimes[[as.character(tr$edge[i,1])]]
    }
    for (i in 1:length(areaperbranch[,7])){
      if(areaperbranch[i,7] == TRUE){
        branchevents[match(i,tr$edge[,2]),3] = 1
      }
    }
    #add the total amount of "time" per slice to fivemilgaps table
    for (i in 1:length(branchevents[,2])){
      for(y in 1:length(fivemilgaps[,2])){
        if ((branchevents[i,1] > (fivemilgaps[y,1]-agegap))&(branchevents[i,1] < (fivemilgaps[y,1]))&(branchevents[i,2] > fivemilgaps[y,1])){
          fivemilgaps[y,5] = fivemilgaps[y,5]+(fivemilgaps[y,1]-branchevents[i,1])
        }
        if ((branchevents[i,1] <= (fivemilgaps[y,1]-agegap))&(branchevents[i,2] < fivemilgaps[y,1])&(branchevents[i,2] > (fivemilgaps[y,1]-agegap))){
          fivemilgaps[y,5] = fivemilgaps[y,5]+(branchevents[i,2]-(fivemilgaps[y,1]-agegap))
        }
        if ((branchevents[i,1] > (fivemilgaps[y,1]-agegap))&(branchevents[i,2] < fivemilgaps[y,1])){
          fivemilgaps[y,5] = fivemilgaps[y,5]+(branchevents[i,2]-branchevents[i,1])
        }
        if ((branchevents[i,1] <= (fivemilgaps[y,1]-agegap))&(branchevents[i,2] > fivemilgaps[y,1])){
          fivemilgaps[y,5] = fivemilgaps[y,5]+(fivemilgaps[y,1]-(fivemilgaps[y,1]-agegap))
        }
      }
    }
    
    #add the absolute amount of change at branches per slice to the fivemilgaps table
    for (i in 1:length(branchevents[,2])){
      if (branchevents[i,3] == 1){
        for(y in 1:length(fivemilgaps[,2])){
          if ((branchevents[i,1] > (fivemilgaps[y,1]-agegap))&(branchevents[i,1] < (fivemilgaps[y,1]))&(branchevents[i,2] > fivemilgaps[y,1])){
            fivemilgaps[y,6] = fivemilgaps[y,6]+((fivemilgaps[y,1]-branchevents[i,1])/(branchevents[i,2]-branchevents[i,1]))
          }
          if ((branchevents[i,1] <= (fivemilgaps[y,1]-agegap))&(branchevents[i,2] < fivemilgaps[y,1])&(branchevents[i,2] > (fivemilgaps[y,1]-agegap))){
            fivemilgaps[y,6] = fivemilgaps[y,6]+(((branchevents[i,2]-(fivemilgaps[y,1]-agegap)))/(branchevents[i,2]-branchevents[i,1]))
          }
          if ((branchevents[i,1] > (fivemilgaps[y,1]-agegap))&(branchevents[i,2] < fivemilgaps[y,1])){
            fivemilgaps[y,6] = fivemilgaps[y,6]+((branchevents[i,2]-branchevents[i,1])/(branchevents[i,2]-branchevents[i,1]))
          }
          if ((branchevents[i,1] <= (fivemilgaps[y,1]-agegap))&(branchevents[i,2] > fivemilgaps[y,1])){
            fivemilgaps[y,6] = fivemilgaps[y,6]+((fivemilgaps[y,1]-(fivemilgaps[y,1]-agegap))/(branchevents[i,2]-branchevents[i,1]))
          }
        }  
      }
    }
    
    #Add column with total number of shifts to fivemilgaps (So branch+node)
    for (i in 1:length(fivemilgaps[,6])){
      fivemilgaps[,8] = fivemilgaps[,6]+fivemilgaps[,7]
    }
    
    #Add column with number of rangeshifts relative to the total number of "time" on the lineages in the timeslice
    for (i in 1:length(fivemilgaps[,6])){
      fivemilgaps[,9] = fivemilgaps[,8]/fivemilgaps[,5]
    }
    
    fivemilgaps = fivemilgaps[,c(1,5,2,3,6,8,9)]
    fivemilgaps[,1] = fivemilgaps[,1]-2.5
    colnames(fivemilgaps) <- c("Age","tot.time","tot.node","Abs.node","Abs.branch","Abs.tot","Rel.tot")
  
    if(exists("fivemilgaps1")){
      fivemilgaps1 <- rbind(fivemilgaps1,fivemilgaps)
      fivemilgaps1 = aggregate(x=fivemilgaps1[,2:7],by=list(fivemilgaps1[,"Age"]),FUN=sum)
      colnames(fivemilgaps1)[1] <- colnames(fivemilgaps)[1]
    }else{
      fivemilgaps1 = fivemilgaps
    }
    fivemilgaps = fivemilgaps1
  }
  
  pdf(file="Dispersal_plots.pdf",paper = "a4")
  par(mfrow=c(2,2))
  plot(fivemilgaps[,"Age"],fivemilgaps[,"Rel.tot"],type = 'l',xlab = "Ma",ylab = "Relative number of dispersals",main="Relative dispersal")
  plot(fivemilgaps[,"Age"],fivemilgaps[,"Abs.tot"],type = 'l',xlab = "Ma",ylab = "Number of dispersal events",main="Absolute dispersal")
  plot(fivemilgaps[,"Age"],fivemilgaps[,"tot.time"],type = 'l',xlab = "Ma",ylab = "Branch length in bin (My)",main="Total branchlength per bin")
  plot(fivemilgaps[,"Age"],fivemilgaps[,"tot.node"],type = 'l',xlab = "Ma",ylab = "Nr of Nodes",main="Number of available nodes per bin")
  dev.off()
}

div.plots(wd = wd, binsize = binsize, area1 = area1, area2 = area2)

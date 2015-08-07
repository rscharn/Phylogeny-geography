#set working directory (where all imput data is present, including shapefiles)
setwd("C:\\Users\\Ruud\\Desktop\\R working Directory\\scripts fernanda\\LASTEST VERSIONS\\Ferns")
binsize = 5 # binsize
area1 = "a+"
area2 = "c+"
area = c("a","b","c","d","e","f","g")
areanames = c("Amazon","Andes","Atlantic_forest","Caribe","Central_America","Dry_diagonal","Dry_north")
################################################
#check number tre files(should be same as nr RDATA files)

div.circle.plots <-function(binsize = 5, area1 = "a", area2 = "b", wd = getwd()){
  number=length(grep("clade_\\d+_area_\\d\\_updated\\.tre",dir(),perl =TRUE))
  fivemilgaps=data.frame()
  
  #load in inputs
  agegap = binsize # binsize
  
  #Load dependancies
  require(BioGeoBEARS)
  require(geiger)
  numberevents = 0

  for (o in 1:number){
    
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
    
    #Check wether nodes are present in the area1 or area2
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
    
    #note nodes that have jumped/expanded from area1 to area2
    areapernode[,7] = 0
    for (i in 1:length(tr$edge[,1])){
      if ((areapernode[tr$edge[i,1],6] != areaperbranch[tr$edge[i,2],6])&(areapernode[tr$edge[i,1],6] != 0)&(areaperbranch[tr$edge[i,2],6] != 0)&
             (areapernode[tr$edge[i,1],6] != "3")&
            (areapernode[tr$edge[i,1],6] == "1")){
        areapernode[tr$edge[i,1],7] = "1"
      }
    }
    
    #note branches on which transatlantic change occured
    areaperbranch[,7] = "0"
    for (i in 1:length(tr$edge[,1])){
      if ((areaperbranch[tr$edge[i,2],6] != areapernode[tr$edge[i,2],6])&(areaperbranch[tr$edge[i,2],6] != 0)&(areapernode[tr$edge[i,2],6] != 0)&
            (areaperbranch[tr$edge[i,2],6] != "3")&
            (areaperbranch[tr$edge[i,2],6] == "1")){
        if (areaperbranch[tr$edge[i,2],7] == "0"){
          areaperbranch[tr$edge[i,2],7] = "1"
          
        }
      }
    }
    numberevents=sum(as.numeric(areapernode[,7]))+sum(as.numeric(areaperbranch[,7]))+numberevents
  }
  return(numberevents)
}

require(caTools)
combinations<-combs(area,2) #create vector with all area combinations
combinations<-rbind(combinations,combinations[,c(2,1)])
events_all = 0
combinations_names<-combs(areanames,2) #create vector with all area combinations
combinations_names<-rbind(combinations_names,combinations_names[,c(2,1)])
circletable = matrix(ncol = 2,nrow =(length(combinations[,1])),combinations_names)
circletable = cbind(circletable,0)
colnames(circletable) <- c("from","to","thickness")

for (r in 1:length(combinations[,1])){ #loop function over all area combinations
  area1 = paste(combinations[r,1],"+",sep="")
  area2 = paste(combinations[r,2],"+",sep="")  
  events<- div.circle.plots(binsize = binsize, area1 = area1, area2 = area2, wd = wd)
  events_all <- events + events_all
  circletable[r,3] = events
}

write.table(circletable, file = "number_of_transitions.txt",sep = "\t",row.names = FALSE,quote = FALSE)

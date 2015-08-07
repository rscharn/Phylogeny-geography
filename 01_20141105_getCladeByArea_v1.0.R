library("speciesgeocodeR")
library("geiger")
library("picante")

setwd("C:\\Users\\Sidonie\\Dropbox\\Carvalho_et_al_Neotropical_biome_shifts\\Data_and_analyses\\Analyses\\final_Biogeobears_analyses\\ferns")


#######################################################
###define function
#######################################################
getCladeByArea<-function(input_gbif, input_tree, input_polygons, percentage_filter = 10,
                         total_occurrence_number_filter = 0, num_poly = 1, area = 1,
                         min_clade_size = 10, monophyly_threshold = .15, summary.out = T, graph.out = F,
                         plug.in.manipulated.occurrence = F, max_clade_size = length(tr$tip.label),
                         plug.in.huge.data.bisse = F, huge.data.bisse.name = F){
  if(plug.in.huge.data.bisse == F){
    pol <- speciesgeocodeR:::.ConvertPoly(input_polygons)
    pol[,2] <- as.numeric(as.character(pol[,2]))
    pol[,3] <- as.numeric(as.character(pol[,3]))
    tt <-read.table(input_gbif)
    names(tt) <- c("indentifier","XCOOR", "YCOOR")
    out <- SpGeoCod(tt, pol)
    if(graph.out ==T){
      speciesgeocodeR:::.OutPlotSpPoly(out)
      speciesgeocodeR:::.OutBarChartPoly(out)
      speciesgeocodeR:::.OutBarChartSpec(out)
      speciesgeocodeR:::.OutMapAll(out)
      speciesgeocodeR:::.OutMapPerSpecies(out)
      speciesgeocodeR:::.OutMapPerPoly(out, scale = scale)
    }
    
    #filtering
    dat <- out$spec_table
    dat2 <- subset(dat,rowSums(data.frame(dat[,-1])) > total_occurrence_number_filter)#filter species with more then 9 records
    rownames(dat2) <- dat2[,1]# setting species names = row names
    others <- aggregate(out$not_classified_samples$XCOOR, by = list(out$not_classified_samples$identifier), length)
    dat2 <- merge(dat2, others, by.x = "identifier", by.y = "Group.1", all.x = T)
    dat2[is.na(dat2)] <-0
    identifier <- as.character(dat2$identifier)
    dat2 <- data.frame(dat2[,-1])# removing species names which is a text column. the function deal only with numerical fields
    names(dat2) <- c(names(out$polygons), "others")
    
    perc <- function(x){
      x / rowSums(dat2) * 100
    }
    
    filtperc <- data.frame(lapply(dat2,perc))
    
    for(i in 1:length(names(filtperc))){
      filtperc[which(filtperc[,i] < percentage_filter),i] <- 0
      filtperc[which(filtperc[,i] >= percentage_filter),i] <- 1
    }
    
    filtperc <- data.frame(identifier, filtperc)
    
    # remove everything that has only 0s
    out2 <- filtperc[,1:(num_poly+1)]
    out2<- out2[rowSums(data.frame(out2[,-1])) > 0,]
    out_geo <- data.frame(out2[,-1])
    rownames(out_geo) <- out2$identifier
    names(out_geo) <- names(out2)[-1]
  } 
  
  #go on here to get out clade
  if(plug.in.huge.data.bisse == T){ #for files from huge databases that had to be chopped before
    geo.in <- read.table(huge.data.bisse.name, sep = "\t", header = T, row.names = 1)
    geo <- data.frame(geo.in[, 1], row.names = rownames(geo.in))
    names(geo) <- names(geo.in)[1]
    
    tree <- read.nexus(input_tree)
    data <- suppressWarnings(treedata(tree,geo))
    
    nodata <- tree$tip.label[!tree$tip.label %in% row.names(geo)]
    nodatatab <- data.frame(rep(NA, length(nodata)), row.names = nodata)
    
    notclasstab <- data.frame(subset(data$data, data$data[, 1] == 0))
    notclass <- rownames(notclasstab)
    
    names(notclasstab) <- names(nodatatab) <- names(data.frame(data$dat))
    all.classified <- rbind(data.frame(data$dat), nodatatab)
    
    geo.in <- subset(geo.in[!rownames(geo.in) == "",],)
    
    others2 <- data.frame(subset(geo.in,geo.in[,1] == 0)[,2], 
                          row.names = rownames(subset(geo.in,geo.in[,1] == 0)))
    names(others2) <- names(geo.in)[1]
  }else{
    tree <- read.tree(input_tree)
    data <- suppressWarnings(treedata(tree,out_geo))
    
    others2 <- data.frame(others[,2], row.names = others[,1])
    names(others2) <- names(out_geo)
    nodata <- tree$tip.label[!tree$tip.label %in% row.names(rbind(out_geo,others2))]
    nodatatab <- data.frame(rep(NA, length(nodata)), row.names = nodata) 
    
    notclass <- tree$tip.label[!tree$tip.label %in% data$phy$tip.label]
    notclass <- notclass[!notclass %in% nodata]
    notclasstab <- data.frame(rep(0, length(notclass)), row.names = notclass)
    
    names(notclasstab) <- names(nodatatab) <- names(data.frame(data$dat))
    all.classified <- rbind(data.frame(data$dat), nodatatab, notclasstab)
    geo <- rbind(data.frame(data$dat), notclasstab, nodatatab)#change
    names(geo) <- names(out_geo) #change
  }
  
  if(summary.out == T){write.table(all.classified,
                                   paste(strsplit(input_gbif, "\\.")[[1]], "_biSSe_table.txt", sep = "")[1], 
                                   quote = F, col.names = F, sep = "\t")}
  
  if(plug.in.manipulated.occurrence == T){
    geo <- read.table(paste(strsplit(input_gbif, "\\.")[[1]], "_biSSe_table_modified.txt", sep = "")[1],
                      sep = "\t", header = F, row.names = 1)
    names(geo) <- names(out_geo) #change
  }
  
  data2 <- treedata(tree, geo)
  tr <- data2$phy
  colo <-  geo[match(tree$tip.label, rownames(geo)),]
  colo <- gsub(1,"red",colo)
  colo <- gsub(0,"blue",colo)
  colo[is.na(colo)] <- "grey"
  
  if(summary.out == T){
    pdf(paste(strsplit(input_gbif, "\\.")[[1]], "_pruning_tree.pdf", sep = "")[1], height = 11, width = 8, paper = "special")
    plot(tree, tip.color = colo, cex = .5)#, type = "fan")
    legend("topleft",bg = "white", fill = c("blue", "red", "grey"),
           legend = c("Outside", "Inside", "No occurrence data"), cex = .7)
    dev.off()
  }
  ###################################################
  #getting monopyletic areas Ruud
  getDescendants <- function(tree,node,curr = NULL){
    if(is.null(curr)) curr <- vector()
    daughters <- tree$edge[which(tree$edge[, 1] == node), 2]
    curr <- c(curr,daughters)
    w <- which(daughters>=length(tree$tip))
    if(length(w) > 0) for(i in 1:length(w)) 
      curr <- getDescendants(tree, daughters[w[i]],curr)
    return(curr)
  } 
  
  #create matrix for each node on state 1 ("not checked")
  nodenr = matrix(nrow = length((length(tr$tip.label)+1):max(tr$edge)), ncol = 2)
  nodenr[,1] = 1 # sets all nodes to 1 ie all can be checked
  nodenr[,2] = 0
  rownames(nodenr) <- (length(tr$tip.label)+1):max(tr$edge)
  colnames(nodenr) <- c("checked_monophyletic","clade_size_mono")
  
  #check for clades monophylethic for area "X" 
  for(i in (length(tr$tip.label)+1):max(tr$edge)){
    if (nodenr[as.character(i),1] == 1){ # skip nodes identified with 0 (not monophyletic for a given area OR descendent of monophyletic clade)
      
      # __ # monophyly_threshold is NUMBER of species not in area of interest
      # __  V=geo[tips(tr,i),area]
      # __  if (length(V[V == 0]) > monophyly_threshold ){
      # __    nodenr[as.character(i),1]=0
      # __  } else{
      
      # monophyly_threshold is percentage of species not in area of interest OUT OF CLADE SIZE
      V <- geo[tips(tr,i),area]
      if (length(V[V == 0]) / length(V) >= monophyly_threshold || length(tips(tr,i)) > max_clade_size){
        nodenr[as.character(i),1] = 0
      } else{
        
        nodenr[as.character(i),2]=length(tips(tr,i)) # number of tips in the clade
        nodes <- getDescendants(tr,i) 
        nodenr[as.character(nodes[nodes>length(tr$tip.label)]),1]=0 # set all descendent noeds to zero (skip check)
        nodenr[as.character(i),1]=2 # indentifier of monophyletic group
      }
    }
  }
  
  # removes from table all nodes which are not of interest
  # clades TABLE:
  # column 1: 0-> discarded (there shouldn't be any); 1-> untested (there shouldn't be any at the end of the loop); 2-> positive hit
  # CALL RUUD for bug reports or if you want to have a good time  (ruud.scharn@gmail.com)
  # column 2: size of the clade
  clades <- subset(nodenr,nodenr[,2] >= min_clade_size,)
  
  #Write output table
  lengtab <- as.vector(0)
  subtreelist <- as.list(0)
  pdf(paste(strsplit(input_gbif, "\\.")[[1]], "_subtrees.pdf", sep = "")[1], height = 11, width = 8, paper = "special")
  for (i in 1:length(clades[,2])){
    node = as.numeric(row.names(clades)[i])
    subtree = extract.clade(tr, node, root.edge = 0, interactive = F)
    lengtab[i] <- max(node.age(subtree)$ages)
    write.tree(subtree, file = paste("clade_", i, "_area_", area,".tre", sep=""))
    
    colo <-  geo[match(subtree$tip.label, rownames(geo)),]
    colo <- gsub(1,"red",colo)
    colo <- gsub(0,"blue",colo)
    colo[is.na(colo)] <- "grey"
    plot(subtree, tip.color = colo, cex = .5, main = paste("Clade", i, sep = "_"))#, type = "fan")
    legend("topleft",bg = "white", fill = c("blue", "red", "grey"),
           legend = c("Outside", "Inside", "No occurrence data"), cex = .7)
    add.scale.bar()
  }
  dev.off()
  
  #Summary+output
  if(plug.in.huge.data.bisse == F){
    logg <- c(dim(dat)[1],
              dim(subset(dat, dat[,2] > 0))[1],
              dim(out_geo)[1],
              Ntip(tree),
              Ntip(tree)-length(nodata),
              length(nodata),
              length(clades[,2]),
              round(mean(clades[,2]),digits = 2)
              )
  }else{
    logg <- c("NA",
              "NA",
              "NA",
              Ntip(tree),
              Ntip(tree)-length(nodata),
              length(nodata),
              length(clades[,2]),
              round(mean(clades[,2]),digits = 2))
  }
  
  desc <- c("GBIF",
            "In polygon before\ngeo filter",
            "In polygon after\ngeo filter",
            "Phylogeny",
            "After matching\nwith the phylogeny",
            "The phylogeny with\nno occurrence data",
            "Number of monophyletic\nclades received\nafter filtering",
            "Mean tip number\nper tree")
  
  descwrite <- c("date",
                 "input_gbif",
                 "input_polygons",
                 "percentage_filter",
                 "total_occurrence_number_filter",
                 "min_clade_size",
                 "max_clade_size",
                 "monophyly_threshold", 
                 "GBIF",
                 "In_polygon_before_geo_filter",
                 "In_polygon_after_geo_filter",
                 "Phylogeny",
                 "After_matching_with_the_phylogeny",
                 "The_phylogeny_with_no_occurrence_data",
                 "Number_of_monophyletic_clades_received_after_filtering",
                 "Mean_tip_number_per_tree")
  
  datum <- as.character(Sys.Date())
  logwrite <- c(datum,
                input_gbif,
                input_polygons,
                percentage_filter,
                total_occurrence_number_filter,
                min_clade_size,
                max_clade_size,
                monophyly_threshold,
                logg)
  
  checklist <- list.files()
  if(!"GetAreasClades_log.txt" %in% checklist){
    write.table(as.vector(data.frame(t(descwrite))),
                "GetAreasClades_log.txt",
                col.names = F, row.names = F, quote = F, sep = "\t")
  }
  
  write.table(as.vector(as.data.frame(t(logwrite))), col.names = F, row.names = F,
              "GetAreasClades_log.txt", 
              append = T, quote = F, sep = "\t")
  if(summary.out == T){
    if(plug.in.huge.data.bisse == T){
      logg <- c(0,
                0,
                0,
                Ntip(tree),
                dim(data2$data)[1],
                length(nodata),
                length(clades[,2]),
                mean(clades[,2]))
    }
    pdf(paste(strsplit(input_gbif, "\\.")[[1]], "_summary_graphs.pdf", sep = "")[1], paper = "a4")
    par(mar = c(8,3,3,3))
    barplot(logg, names.arg = paste(desc, " (", logg, ")", sep = ""), las = 2, cex.names = .7, col = "brown",
            main = paste(strsplit(input_gbif, "\\.")[[1]], " species numbers", sep = "")[1])
    box("plot")
    hist(lengtab, main = "Node age", xlab = "Node age")
    hist(clades[,2], main = "Number of tips per clade", xlab = "Number of tips per clade")
    dev.off()
  }
  
  clades_write= clades
  clades_write <- data.frame(clades_write, 
                             node_nr = as.numeric(rownames(clades_write)))[, c(3, 1, 2)] #add node number as column
  rownames(clades_write) <- 1:length(clades_write[,1]) # add number of subtree as rowname
  
  if(summary.out == T){write.table(clades_write, paste("clade_area_", area,".txt", sep=""), sep = "\t")}
  
  if(plug.in.huge.data.bisse == F){
    cat(paste("GBIF: ", dim(dat)[1], "\n"))
    cat(paste("In polygon before geographic filter: ",dim(subset(dat, dat[,2] > 0))[1], "\n"))
    cat(paste("In polygon after geographic filter: ", dim(out_geo)[1], "\n"))
  }
  cat(paste("Tips in Phylogeny: ", Ntip(tree), "\n"))
  cat(paste("Species in the phylogeny with no occurrence data:  ", length(nodata), "\n", sep = ""))
  cat(paste("After matching with the phylogeny ", dim(geo)[1], "\n", sep = ""))
  cat(paste("Number of monophyletic clades received after filtering:  ", length(clades[,2]), "\n", sep = ""))
  cat(paste("Mean tip number per tree:  ", mean(clades[,2]), "\n", sep = ""))
}

#######################################################
#run function
#######################################################
# ######################      

getCladeByArea(input_gbif = "ferns_occurence.txt", 
               input_tree = "ferns.tre",
               input_polygons = "neotropics_modified_Ruud.txt",
               percentage_filter = 10,
               total_occurrence_number_filter = 0,
               num_poly = 1,
               area = 1,
               min_clade_size = 10,
               monophyly_threshold = .15,
               summary.out = T,
               graph.out = F,
               plug.in.manipulated.occurrence = F,
               max_clade_size = 100,
               plug.in.huge.data.bisse = F,
               huge.data.bisse.name = "input\\birds_final_coding_Amazon-Guiana_shield-200m_red.txt")

##Arguments you probably want to change
# input_gbif <- "turtles.txt"
# input_tree <- "Turtles.tre"
# input_polygons <- "South_America_polygons.txt"
# percentage_filter <- 10 # the min. percentage occurring in a determined area
# total_occurrence_number_filter <- 0 # the min number of records per species should be in the complete dataset. In this case, at least 10 records per species
# monophyly_threshold <- 0.1 # if 0.2 it means that 20% of taxa are allowed to be outside the area 
# summary.out = T produce summary graphs?
# graph.out = F plot speciesgeocoder maps?
# min_clade_size = 6# minimum size of a sub-clade
# max_clade_size set to maxim size of the output sub-trees (e.g 50)

##Arguments you probably do not want to change
# num_poly <- 1 # number of polygons might not be working DO NOT CHANGE 
# area = 1 #Set number of area for which you want to select monophyletic clades #DO NOT CHANGE
# plug.in.manipulated.occurrence = F use mdified occurrences file (if you correted occurrences in the BiSSE file), 
# to do this change the bisse file from the output folder, copy it to the working directory and add _modified to the name
#plug.in.huge.data.bisse = F # to use input from huge inputdata that has been prepared before
#huge.data.bisse.name = "sa_plants_final_filtered.txt" #file name of the huge input data

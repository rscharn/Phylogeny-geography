#require(rgdal)
#require(rgeos)
require(Rcpp)
require(speciesgeocodeR)
require(geiger)
require(maptools)

#set working directory (where all imput data is present, including shapefiles)
setwd("C:\\Users\\Ruud\\Desktop\\R working Directory\\scripts fernanda\\LASTEST VERSIONS\\Ferns")

################################################################################################
#define the main function of the script
################################################################################################
generate.biogeo <- function(wwfshape, input_table, minocc, minfraq) {
  
  wwf <- readShapeSpatial(wwfshape)
  
  input <- read.table(input_table, header = T, sep = "\t", stringsAsFactors = TRUE, comment.char = "")
  input <- unique(input)
  names(input) <- c("identifier", "XCOOR", "YCOOR")# retain only unique entries
  
  # run speciesgeocoder
  out <- SpGeoCod(input, wwf, areanames = "ECO_NAME")
  dat <- out$spec_table
  rownames(dat) <- dat[, 1]
  dat <- data.frame(dat[, -1])  # removing species names which is a text column. the function deal only with numerical fields
  
  # minimum occurence filter (at 1 removes only rows without occurence)
  dat <- subset(dat, rowSums(dat) >= minocc)
  
  # fraqtion of total presence filter
  dat_rowsums = rowSums(dat)
  for (i in 1:length(dat)) {
    dat[, i] = dat[, i]/dat_rowsums
    dat[, i][dat[, i] < minfraq] <- 0
    dat[, i] = dat[, i] * dat_rowsums
  }
  
  # Create updated trees/geography files for biogeobears
  number = length(grep("clade_\\d+_area_\\d\\.tre", dir(), perl = TRUE))
  for (i in 1:number) {
    # Load in tree
    tr <- read.tree(paste("clade_", i, "_area_", "1", ".tre", sep = ""))
    
    dat2 <- dat
    data <- suppressWarnings(treedata(tr, dat2))
    tr <- data$phy
    dat2 <- data$data
    dat2[dat2 >= 1] <- 1
    
    # create/print biogeobers command line
    length_dat <- length(dat2[, 1])
    letters_data <- letters[1:length(dat2[1, ])]
    letters_data2 <- paste(letters_data[1:length(letters_data) - 1], " ", sep = "", collapse = "")
    letters_data <- paste("(", letters_data2, letters_data[length(letters_data)], ")", sep = "")
    command_string <- paste(length_dat, length(dat2[1, ]), letters_data, sep = "\t")
    write(command_string, file = paste("clade_", i, "_area_", "1", ".geog", sep = ""))
    
    # create/print BioGeoBEARS table
    biogeo <- matrix(nrow = length_dat, ncol = 2)
    biogeo[, 1] <- rownames(dat2)
    for (o in 1:length_dat) {
      biogeo[o, 2] <- paste(dat2[o, ], sep = "", collapse = "")
    }
    write.table(biogeo, file = paste("clade_", i, "_area_", "1", ".geog", sep = ""), append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE, 
                row.names = FALSE)
    write.tree(tr, file = paste("clade_", i, "_area_", "1", "_updated", ".tre", sep = ""))
  }
}

################################################################################################
#Load input and run function
################################################################################################

input_table = "ferns_occurence.txt"  #give input table e.g lizards_occurrence.txt
wwfshape = "neotropics_margins_expanded_caribe"  #name of the wwf ecoregion edited shapefile to use

minocc = 1  #minimum number of occurence filter 1 = at least one occurence
minfraq = 0.1  # fraction of total filter, 0.1 is an area needs at least 10% of its total occurences in an area in order to be coded

generate.biogeo(input_table = input_table, minocc = minocc, minfraq = minfraq, wwfshape = wwfshape) 

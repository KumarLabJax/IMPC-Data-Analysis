#This code generates histograms, density plots, and QQ plots for IMPC data 

#This sets the directory and procedure_code variables for usage in the
#following functions, you must put the / at the end of the directory. The
#procedure code refers to the three digit code used to define whatever
#protocol you are hoping to analyze, this link http://www.mousephenotype.org/impress/procedures/7
#will take you to a list of all IMPC protocols
directory <- "C:/Users/Branden/Desktop/IMPC/"
procedure_code <- "OFD"

library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)
library(lubridate)
library(gdata)
library(readr)
library(tidyr)
library(Hmisc)
library(plyr)
library(zoo)
library(signal)

getProcedureTable <- function(procedure_code, directory){
  dataDir <- paste(directory,procedure_code,"\\ScrapedProcedureTable",sep="")
  setwd(dataDir)
  # will get ONLY THE FIRST procedure table in the directory <- Needs to be improved to take into account multiple procedure tables
  files <- list.files(path=dataDir, pattern = "*.csv",full.names=TRUE, recursive = FALSE)
  procedure_table <- read_delim(files[1], delim = ',')
  colnames(procedure_table)[1] <- "all_titles"
  procedure_table$all_titles <- make.names(procedure_table$all_titles)
  return(procedure_table)
}

#Loads the IMPC data dictionary 
getIMPCDictionary <- function(directory){
  setwd(directory)
  IMPC_dictionary <- read_delim("IMPC_DataDictionary.csv", delim = ',')
  IMPC_dictionary$Parameter <- make.names(IMPC_dictionary$Parameter)
  return(IMPC_dictionary)
}

#Gets a list of all of the finalized data files, organized by phenotyping center
getWideCSVFiles <- function(procedure_code, directory){
  dataDir <- paste(directory,procedure_code,"//WideData",sep="")
  setwd(dataDir)
  files <- list.files(path=dataDir, pattern = "*.csv",full.names=TRUE, recursive = FALSE)
  return(files)
}

#Loads in the finalized data from a specific file
getCenterData <- function(filename, procedure_code = ""){
  file_extension <- strsplit(filename, "/")[[1]][2]
  center_data <- read_delim(filename, delim = ',')
  # center_data <- read.csv(filename, stringsAsFactors = F)
  if('Start_hms_new' %in% colnames(center_data)){
    center_data <- read.csv(filename,colClasses = c("Start_hms_new" = "character"))
  }
  colnames(center_data) <- make.names(colnames(center_data))
  if(procedure_code == "OFD" && is.element("Center.distance.travelled", colnames(center_data)) && is.element("Periphery.distance.travelled", colnames(center_data))){
    center_data$Total.distance.travelled <- center_data$Center.distance.travelled + center_data$Periphery.distance.travelled
  }
  return(center_data)
}

#Gets the name of the phenotyping center for the loaded data
getCenterName <- function(filename){
  file_extension <- strsplit(filename, "//")[[1]][2]
  file_extension <- strsplit(file_extension, "/")[[1]][[2]]
  center_name <- strsplit(file_extension, "_")[[1]][1]
}

#Gets a list of plottable phenotypes from the loaded data, data dictionary, and procedure_table
getPhenotypes <- function(procedure_table, IMPC_dictionary, center_data, procedure_code){
  
  
  # IMPC_dictionary <- IMPC_dictionary[IMPC_dictionary$Protocol == procedure_code,]
  IMPC_dictionary <- IMPC_dictionary[IMPC_dictionary$Phenotype == 'yes',]
  IMPC_dictionary <- IMPC_dictionary[IMPC_dictionary$Use == 'yes',]
  # gets list of desired subset of phenotypes from the DataDictionary
  index = 1
  
  desired_phenotypes <- unique(IMPC_dictionary$Parameter)
  
  # gets list of phenotypes with data
  index = 1
  
  phenotype_list <- list()
  
  if(procedure_code == "OFD"){
    phenotype_list <- list("Total.distance.travelled")
  }
  for(i in desired_phenotypes){
    if(is.element(i,colnames(center_data))){
      unique_vals <- length(unique(center_data[[i]]))
      
      if(any(is.na(center_data[[i]]))){
        unique_vals <- unique_vals - 1
      }
      
      if(unique_vals>1 && !is.null(center_data[[i]])){
        phenotype_list[[i]] <- i
      }
    }
  }
  return(phenotype_list)
}

#Gets a list of plottable metadata from the loaded data, data dictionary, and procedure_table
getMetadata <- function(procedure_table, IMPC_dictionary, center_data, procedure_code, pipeline = "IMPC"){
  full_metadata_list <- list()
  IMPC_dictionary <- IMPC_dictionary[IMPC_dictionary$Protocol==procedure_code,]
  # extract phenotype and metadata lists from procedure_table
  index = 1
  on_metadata <- FALSE
  for(i in 1:(nrow(procedure_table)-1)){
    if(procedure_table$all_titles[i] == 'NA.'){
      index = 1
      on_metadata<-TRUE
    }
    if(on_metadata){
      full_metadata_list[[index]] <- procedure_table$all_titles[i+1]
    }      
    index <- index + 1 
  }
  full_metadata_list <- lapply(full_metadata_list, function(x) strsplit(x, paste(".",pipeline,sep=""))[[1]][1])
  
  # gets list of desired subset of metadata from the DataDictionary
  index = 2
  desired_metadata <- list('production_phenotype')
  for(i in full_metadata_list){
    if(is.element(i,IMPC_dictionary$Parameter)){
      row_num <- which(IMPC_dictionary$Parameter==i)
      # check if "use" is marked in the data dictionary
      if(!is.na(IMPC_dictionary$Use[row_num]) && IMPC_dictionary$Use[row_num]=="yes" && IMPC_dictionary$Field[row_num]=="factor"){
        desired_metadata[index] <- i
        index <- index + 1
      }
    }
  }
  
  # gets list of plottable metadata
  index = 1
  metadata_list <- list()
  for(i in desired_metadata){
    if(is.element(i,colnames(center_data))){
      unique_vals <- length(unique(center_data[[i]]))
      if(any(is.na(center_data[[i]]))){
        unique_vals <- unique_vals - 1
      }
      
      if(unique_vals>1 && unique_vals<100 && i!="Start.Time"){
        # print(paste("metadata = ",i," ; unique_vals = ", unique(center_data[[i]])[1]," ", unique(center_data[[i]])[2], sep = ""))
        metadata_list[index] <- i
        index <- index + 1
      }
    }
  }
  return(metadata_list)
}

#Defines the multiplot function, allowing for multiple ggplot plots on one page/window
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#Generates histograms, density plots, and QQ plots for a singular phenotype and center 
produceHistDensity <- function(center_data, phenotype, title, center_name){
  
  y <- center_data[[phenotype]] 
  
  Sex <- center_data$sex
  
  hist0 <- ggplot() +
    theme_bw() +
    scale_fill_brewer(guide = guide_legend(),palette = 'Set1') +
    geom_histogram(aes(x = y), bins = 50)
  
  hist1 <- ggplot() +
    geom_density(aes(x = y,y = ..density..)) +
    theme_bw() + labs(title = paste(center_name,"_ALL_",phenotype)) 
    theme(plot.title=element_text(size=8))
  
  hist2 <- ggplot() +
    geom_density(aes(x = y,y = ..density..,fill = Sex),alpha = 0.2784) +
    theme_bw() +
    scale_fill_brewer(guide = guide_legend(),palette = 'Set1')  + 
    theme(legend.position="none")
  
  hist4 <- qplot(sample = y, stat = "qq")
  
  multiplot(hist1, hist0, hist2, hist4, cols=3)
  
}
dev.off()

#Generates histograms, density plots, and QQ plots for every phenotype across all phenotyping centers
produceAllHistDensity <- function(procedure_code, directory, isControl = TRUE, pipeline = "IMPC") {
  procedure_table <- getProcedureTable(procedure_code, directory)
  IMPC_dictionary <- getIMPCDictionary(directory)
  files <- getWideCSVFiles(procedure_code, directory)
  output_directory <- paste(directory, procedure_code, sep="")
  setwd(output_directory)
  date <- format(Sys.time(),"%Y%m%d_%HH%MM_")
  if(isControl){
    control_or_not <- "Controls"
  }else{
    control_or_not <- "All"
  }
  pdf(paste(date, control_or_not, "_DENSITYHIST_",procedure_code,".pdf", sep = ""), 14,4)
  for(filename in files){
    center_data <- getCenterData(filename, procedure_code)
    center_name <- getCenterName(filename)
    phenotype_list <-getPhenotypes(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code)
    if(isControl){
      center_data <- center_data[center_data$biological_sample_group == 'control',]
    }
    if(nrow(center_data)>1){
      for(phenotype in phenotype_list){
        if(is.element(phenotype, colnames(center_data))){
              title <- paste(center_name,phenotype,"\n","(",control_or_not,";",")",sep=" ")
              titleprint <- paste("Creating",center_name,phenotype,"DensityHist Plots")
              density_hist_plot <- produceHistDensity(center_data = center_data, phenotype = phenotype, title = title, center_name = center_name)
              print(titleprint)
              print(density_hist_plot)
            }
            
          }
        
    }
  

print(paste("COMPLETED DENSITYHIST: ",control_or_not," ",center_name,sep=""))

}
dev.off()
}

#This generates two sets of densityhist plots organized by phenotyping center, 
#one for the control set of mice and one for all of the mice
produceAllHistDensity(procedure_code, directory, isControl = T, pipeline = "IMPC")
produceAllHistDensity(procedure_code, directory, isControl = F, pipeline = "IMPC")

#This code generates boxplots for IMPC data 

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

#Loads the scraped procedure table
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

#Produces a boxplot for a specific phenotype from the loaded data 
produceBoxPlot <-function(center_data, phenotype, factors, title = ""){
  x <- factors
  y <- phenotype
  center_data <- center_data[complete.cases(center_data[[x]]),]
  center_data <- center_data[complete.cases(center_data[[y]]),]
  center_mean <- mean(center_data[[y]], na.rm = TRUE)
  center_std <- sd(center_data[[y]], na.rm = TRUE)
  if(length(unique(center_data[[x]]))>1){
    print(paste(phenotype,factors,sep=" "))
    plot <- ggplot() + 
      geom_point(aes(x = as.factor(center_data[[x]]),y = center_data[[y]] ,colour = sex),data=center_data,shape = 16,size = 1.0, position = position_jitter(w = 0.2, h = 0)) +  
      scale_colour_brewer(guide = guide_legend(),palette = 'Set1') +
      geom_crossbar(aes(y = center_data[[y]],x = reorder(as.factor(center_data[[x]]),center_data[[y]], FUN = mean, na.rm=TRUE), na.rm=TRUE),data=subset(center_data, !is.na(center_data[[y]])),na.rm = TRUE, colour = "Blue", width = 0.5, fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') + 
      theme_bw() + labs(x=x, y=y) + ggtitle(title) +
      geom_hline(yintercept = center_mean, color = "grey", size = 1) + geom_hline(yintercept = center_mean + center_std, color = "grey", size = 1) + geom_hline(yintercept = center_mean - center_std, color = "grey", size = 1) +
      theme(axis.text=element_text(size=16),axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold")) 
    return(plot)
  }
}

#Produces boxplots for every phenotype and phenotyping center
produceAllBoxPlot <- function(procedure_code, directory, isControl = TRUE, pipeline = "IMPC"){
  procedure_table <- getProcedureTable(procedure_code, directory)
  IMPC_dictionary <- getIMPCDictionary(directory)
  IMPC_dictionary <- IMPC_dictionary[IMPC_dictionary$Protocol == procedure_code,]
  files <- getWideCSVFiles(procedure_code, directory)
  output_directory <- paste(directory, procedure_code, sep="")
  setwd(output_directory)
  date <- format(Sys.time(),"%Y%m%d_%HH%MM_")
  if(isControl){
    control_or_not <- "Controls"
  }else{
    control_or_not <- "All"
  }
  pdf(paste(date,control_or_not,"_BOXPLOT_",procedure_code,".pdf", sep = ""), 6,6)
  for(filename in files){
    center_name <- getCenterName(filename)
    print(paste("STARTED ", center_name))
    center_data <- getCenterData(filename,procedure_code)
    phenotype_list <- getPhenotypes(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code)
    metadata_list <- getMetadata(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code, pipeline = pipeline) 
    if(isControl){
      center_data <- center_data[center_data$biological_sample_group == 'control',]
    }
    for(metadata in metadata_list){
      if(is.element(metadata, colnames(center_data))){
        # print(paste("METADATA2: ", metadata,sep=""))
        
        for(phenotype in phenotype_list){
          if(is.element(phenotype, colnames(center_data))){
            # print(paste("METADATA2x: ", metadata,sep=""))
            
            if(length(unique(center_data[[metadata]]))>1){
              # print(paste("METADATA3: ", metadata,sep=""))
              
              row_num <- which(IMPC_dictionary$Parameter==metadata)
              
              if(!is.na(IMPC_dictionary$MetadataCategory[row_num])){
                metadata_split <- IMPC_dictionary$MetadataCategory[row_num]
                title <- paste(center_name,metadata,";",phenotype,"\n","(",control_or_not,";",metadata_split,")",sep=" ")
              }else{
                title <- paste(center_name,metadata,";",phenotype,control_or_not,sep=" ")
              }
              # print(paste("HERE IS METADATA: ", metadata))
              boxplot <- produceBoxPlot(center_data = center_data, phenotype = phenotype, factors = metadata, title = title)
              if(!is.null(boxplot)){
                print(boxplot)
              }
            }
          }
        }
      }
    }
    print(paste("COMPLETED BOXPLOT: ",control_or_not, " ",center_name,sep=""))
  }
  dev.off()
}

#This generates two sets boxplots organized by center for every phenotype and\
#metadata combination from that center, one set for the control mice and
#one set for all of the mice
produceAllBoxPlot(procedure_code, directory, isControl = T)
produceAllBoxPlot(procedure_code, directory, isControl = F)

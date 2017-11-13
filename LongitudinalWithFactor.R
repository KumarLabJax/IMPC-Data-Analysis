#This code generates the longitudinal plots for IMPC data using metadata as an added factor

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


getTimeSeries <- function(IMPC_dictionary, center_data, procedure_code){
  timeseries <- list("date_of_experiment")
  IMPC_dictionary <- IMPC_dictionary[IMPC_dictionary$Protocol==procedure_code,]
  IMPC_dictionary <- IMPC_dictionary[IMPC_dictionary$Field=="time",]
  timeseries <- c(timeseries, IMPC_dictionary$Parameter)
  return(timeseries)
}

#Produces a timeseries plot using a singular phenotype with a single category of metadata as an added factor
produceMultiTimeSeries <- function(center_data, phenotype, timeplot, factors, title = "", moving_average = 1, isLine = TRUE){
  x <- timeplot
  y <- phenotype
  z <- factors
  center_data <- center_data[complete.cases(center_data[[x]]),]
  center_data <- center_data[complete.cases(center_data[[y]]),]
  center_data <- center_data[complete.cases(center_data[[z]]),]
  center_mean <- mean(center_data[[y]], na.rm = TRUE)
  center_std <- sd(center_data[[y]], na.rm = TRUE)
  center_data[[x]] <- parse_date_time(center_data[[x]],orders = c("Ymd", "Ymd H:M:S", "Y/m/d", "m/d/Y", "mdy", "mdy H:M","mdY"))
  
  if(length(unique(center_data[[z]]))>1){
    geom <- list()
    
    temp <- data.frame()
    for(val in unique(center_data[[z]])){
      if(sum(center_data[[z]]==val)>1){
        val_data <- center_data[center_data[[z]] == val,]
        temp <- rbind(temp, data.frame(x = val_data[[x]],y = val_data[[y]],color=toString(val)))
      }
    }
    
    
    if(isLine){
      g <- geom_line(data=temp,aes(x=x,y=y,col=color,group=color),na.rm=T)
    }else{
      g <- geom_point(data=temp,aes(x=x,y=y,col=color),na.rm=T, alpha=0.5)
    }
    geom <- c(geom, g)
    
    # PREVIOUS ATTEMPT AT ROLLING MEAN
    # center_zoo <- zoo(val_data[[y]], val_data[[x]])
    # m_av <- rollmean(center_zoo, moving_average, fill=list(NA,NULL,NA))
    # temp <- data.frame(x = val_data[[x]],y = coredata(m_av), color = toString(val))
    # if(isLine){
    #   g <- geom_line(data=temp,aes(x=x,y=y,col=color),na.rm=T)
    # }else{
    #   g <- geom_point(data=temp,aes(x=x,y=y,col=color),na.rm=T, alpha=0.5)
    # }  , group = combined_df$phenotyping_center
    
    
    plot <- ggplot() + unlist(geom) + labs(x=x,y=y,colour=z) + ggtitle(title) +
      geom_hline(yintercept = center_mean, color = "grey", size = 1) + geom_hline(yintercept = center_mean + center_std, color = "grey", size = 1) + geom_hline(yintercept = center_mean - center_std, color = "grey", size = 1) +
      theme(axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"), axis.text=element_text(size=16))
    return(plot)
  }else{
    return("")
  }
}

#Produces every phenotype metadata combination multi time series plot for every phenotyping center
produceAllMultiTimeSeries <- function(procedure_code, directory, isControl = TRUE, moving_average = 1, isLine = TRUE, pipeline = "IMPC"){
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
  if(isLine){
    pdf(paste(date, control_or_not,"_MULTITIMESERIES_",procedure_code,".pdf", sep = ""), 14,4)
    line_or_not <- "ASLINE"
  }else{
    pdf(paste(date, control_or_not, "_ASPOINTS", "_MULTITIMESERIES_",procedure_code,".pdf", sep = ""), 14,4)
    line_or_not <- "ASPOINTS"
  }
  for(filename in files){
    center_data <- getCenterData(filename, procedure_code)
    center_name <- getCenterName(filename)
    phenotype_list <-getPhenotypes(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code)
    timeplot_list <- getTimeSeries(IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code)
    factor_list <- getMetadata(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code, pipeline = pipeline) 
    
    
    if(isControl){
      center_data <- center_data[center_data$biological_sample_group == 'control',]
    }
    if(nrow(center_data)>1){
      for(phenotype in phenotype_list){
        if(is.element(phenotype, colnames(center_data))){
          for(timeplot in timeplot_list){
            if(is.element(timeplot, colnames(center_data))){
              for(factors in factor_list){
                if(is.element(factors, colnames(center_data))){
                  title <- paste(center_name,phenotype,";", timeplot,";",factors,"\n","(",control_or_not,";",")",sep=" ")
                  if(isLine == TRUE){
                    multitime_series_plot <- produceMultiTimeSeries(center_data = center_data, phenotype = phenotype, timeplot = timeplot, factors = factors, title = title, moving_average = moving_average, isLine = isLine)
                    print(multitime_series_plot)
                  }else{
                    multitime_series_plot <- produceMultiTimeSeries(center_data = center_data, phenotype = phenotype, timeplot = timeplot, factors = factors, title = title, moving_average = moving_average, isLine = isLine)
                    print(multitime_series_plot)
                  }
                }
                
              }
            }
          }
        }
      }
    }
    print(paste("COMPLETED MULTITIMESERIES: ",control_or_not," ",line_or_not," ",center_name,sep=""))
    
  }
  dev.off()
}

#This generates four sets of longitudinal plots organized by phenotyping center 
#using metadata as an added factor, creating every combination of control or 
#not control and connecting the plots with a line or with points
produceAllMultiTimeSeries(procedure_code, directory, isControl = T, isLine = T)
produceAllMultiTimeSeries(procedure_code, directory, isControl = F, isLine = T)
produceAllMultiTimeSeries(procedure_code, directory, isControl = T, isLine = F)
produceAllMultiTimeSeries(procedure_code, directory, isControl = F, isline = F)

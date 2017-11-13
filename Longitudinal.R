#This code generates the longitudinal plots for IMPC data

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
  #This line only gets the first procedure table from the directory
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

#Loads the final data from a specific file, this is combined with getWideCSVFiles to load the data
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

getTimeSeries <- function(IMPC_dictionary, center_data, procedure_code){
  timeseries <- list("date_of_experiment")
  IMPC_dictionary <- IMPC_dictionary[IMPC_dictionary$Protocol==procedure_code,]
  IMPC_dictionary <- IMPC_dictionary[IMPC_dictionary$Field=="time",]
  timeseries <- c(timeseries, IMPC_dictionary$Parameter)
  return(timeseries)
}

produceTimeSeries <- function(center_data, phenotype, timeplot, title = "", moving_average = 1){
  x <- timeplot
  y <- phenotype
  center_data <- center_data[complete.cases(center_data[[x]]),]
  center_data <- center_data[complete.cases(center_data[[y]]),]
  center_mean <- mean(center_data[[y]], na.rm = TRUE)
  center_std <- sd(center_data[[y]], na.rm = TRUE)
  center_data[[x]] <- parse_date_time(center_data[[x]],orders = c("Ymd", "Ymd H:M:S", "Y/m/d", "m/d/Y", "mdy", "mdy H:M","mdY"))
  center_data[[x]] <- as.Date(center_data[[x]])
  plot <- ggplot() +
    geom_point(aes(x = center_data[[x]],y = center_data[[y]],colour = sex),data=center_data,shape = 16,size = 1.0) +
    geom_line(aes(x = center_data[[x]],y = center_data[[y]]),data=center_data,size = 1.0,fun.data = mean_sdl,stat = 'summary') +
    scale_colour_brewer(guide = guide_legend(),palette = 'Set1') + labs(x=x, y=y) + ggtitle(title) +
    geom_hline(yintercept = center_mean, color = "grey", size = 1) + geom_hline(yintercept = center_mean + center_std, color = "grey", size = 1) + geom_hline(yintercept = center_mean - center_std, color = "grey", size = 1) +
    theme(axis.text=element_text(size=16),axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"))
  return(plot)
}

#Produces a longitudinal timeplot for every phenotype for every phenotyping center
produceAllTimeSeries <- function(procedure_code, directory, isControl = TRUE, moving_average = 1){
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
  pdf(paste(date, control_or_not,"_TIMESERIES_",procedure_code,".pdf", sep = ""), 14,4)
  for(filename in files){
    center_data <- getCenterData(filename, procedure_code)
    center_name <- getCenterName(filename)
    phenotype_list <-getPhenotypes(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code)
    timeplot_list <- getTimeSeries(IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code)
    
    if(isControl){
      center_data <- center_data[center_data$biological_sample_group == 'control',]
    }
    if(nrow(center_data)>1){
      for(phenotype in phenotype_list){
        if(is.element(phenotype, colnames(center_data))){
          for(timeplot in timeplot_list){
            if(is.element(timeplot, colnames(center_data))){
              title <- paste(center_name,phenotype,";", timeplot,"\n","(",control_or_not,")",sep=" ")
              time_series_plot <- produceTimeSeries(center_data = center_data, phenotype = phenotype, timeplot = timeplot, title = title, moving_average = moving_average)
              # tryCatch(print(time_series_plot),error=function(x){print(paste(center_name, phenotype, timeplot,sep=" "))})
              print(time_series_plot)
            }
          }
        }
      }
    }
  }
  dev.off()
}

#This generates two sets of longitudinal plots organized by phenotyping center 
#one for the control set of mice and one for all of the mice
produceAllTimeSeries(procedure_code, directory, isControl = T)
produceAllTimeSeries(procedure_code, directory, isControl = F)

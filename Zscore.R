#This code generates Z-score allele plots for IMPC data 

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

#Produces an allele plot for a singular phenotype
produceAllelePlot <- function(center_data, phenotype, title = "", isSorted = FALSE, center_name = "", z_score_plots = FALSE){
  x <- "allele_symbol_zygosity"
  y <- phenotype
  
  
  center_data$allele_symbol <- as.character(center_data$allele_symbol)
  center_data$strain_name <- as.character(center_data$strain_name)
  center_data$allele_symbol[(center_data$allele_symbol == "" | is.na(center_data$allele_symbol)) & center_data$biological_sample_group == 'control'] <- center_data$strain_name[(center_data$allele_symbol == "" | is.na(center_data$allele_symbol)) & center_data$biological_sample_group == 'control']
  
  # center_data$zygosity <- " "
  center_data$allele_symbol_zygosity <- paste(center_data$allele_symbol,center_data$zygosity,sep=" ")
  controls_only <- center_data[center_data$biological_sample_group == 'control',]
  controlMean <- mean(controls_only[[y]], na.rm = T)
  controlSD <- sd(controls_only[[y]], na.rm = T)
  center_data$z_score <- (center_data[[y]]-controlMean) / controlSD
  if(z_score_plots){
    plot <- ggplot() +
      geom_pointrange(aes(x = reorder(center_data[[x]], center_data$z_score, FUN = mean,na.rm=TRUE),y = center_data$z_score, na.rm = TRUE),data= center_data,colour = '#0000ff',fill = '#ffffff',size = 0.8, shape = 20 ,fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') + 
      geom_point(aes(x = center_data[[x]],y = center_data$z_score, na.rm = TRUE),data=center_data,shape = 20,colour = '#999999',size = 1.0, alpha = 0.3) +
      theme_bw(base_size = 4.0) + labs(title = title) +
      geom_hline(yintercept = -1, color  = "red") + 
      geom_hline(yintercept = 0, color  = "red") +
      geom_hline(yintercept = 1, color  = "red") +
      labs(x=x,y=y) +
      coord_flip() +
      theme(axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"))
  }else{
    if(isSorted){
      plot <- ggplot() +
        geom_pointrange(aes(x = reorder(center_data[[x]], center_data[[y]], FUN = mean,na.rm=TRUE),y = center_data[[y]], na.rm = TRUE),data= center_data,colour = '#0000ff',fill = '#ffffff',size = 0.8, shape = 20 ,fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') + 
        geom_point(aes(x = center_data[[x]],y = center_data[[y]], na.rm = TRUE),data=center_data,shape = 20,colour = '#999999',size = 1.0, alpha = 0.3) +
        scale_colour_brewer(guide = guide_legend(),palette = 'Set1') + theme_bw(base_size = 4.0) + labs(title = title) +
        geom_hline(yintercept = controlMean, color  = "red") + 
        geom_hline(yintercept = controlMean + controlSD, color  = "red") +
        geom_hline(yintercept = controlMean - controlSD, color  = "red") +
        labs(x=x,y=y) +
        coord_flip() +
        theme(axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"))
      return(plot)
    }else{
      plot <- ggplot() +
        geom_point(aes(x = center_data[[x]],y = center_data[[y]], na.rm = TRUE),data= center_data,shape = 20,colour = '#999999',size = 1.0, alpha = 0.3) +
        geom_pointrange(aes(x = center_data[[x]],y = center_data[[y]], na.rm = TRUE),data= center_data,colour = '#0000ff',fill = '#ffffff',size = 0.8, shape = 20,fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary',na.rm=T) + 
        scale_colour_brewer(guide = guide_legend(),palette = 'Set1') + theme_bw(base_size = 4.0) + labs(title = title) + 
        geom_hline(yintercept = controlMean, color  = "red") + 
        geom_hline(yintercept = controlMean + controlSD, color  = "red") +
        geom_hline(yintercept = controlMean - controlSD, color  = "red") + 
        labs(x=x,y=y) +
        coord_flip() +
        theme(axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"))
      return(plot)
    }
  }
}

produceAllZScorePlots <- function(procedure_code, directory, plotAllFile = T){
  procedure_table <- getProcedureTable(procedure_code, directory)
  IMPC_dictionary <- getIMPCDictionary(directory)
  files <- getWideCSVFiles(procedure_code, directory)
  output_directory <- paste(directory, procedure_code, sep="")
  setwd(output_directory)
  date <- format(Sys.time(),"%Y%m%d_%HH%MM_")
  finished_phenotypes <- list()
  if(plotAllFile){
    pdf(paste(date,"_ALL_ZPlots_",procedure_code,".pdf", sep = ""), onefile = T, 10,50)
  }
  
  for(filename in files){
    center_data <- getCenterData(filename, procedure_code)
    # View(center_data[,c("Total.distance.travelled", "biological_sample_id")])
    phenotype_list <-getPhenotypes(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code)
    for(phenotype in phenotype_list){
      if(!is.element(phenotype, finished_phenotypes)){
        combined_df <- data.frame()
        #combined_df <- data.frame(biological_sample_id = integer(0), phenotyping_center = character(0), allele_symbol_zygosity = character(0), z_score = numeric(0))
        
        if(!plotAllFile){
          pdf(paste(date,"_",phenotype,"_ZPlots_",procedure_code,".pdf", sep = ""), 10,20)
        }
        for(current in files){
          current_data <- getCenterData(current, procedure_code)
          current_name <- getCenterName(current)
          current_list <-getPhenotypes(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = current_data, procedure_code = procedure_code)
          if(is.element(phenotype, current_list) & is.element(phenotype,colnames(current_data))){
            
            
            # produce graph for single center
            if(!plotAllFile){
              title <- paste(current_name,"allele_symbol",";",phenotype,("z-scores"),sep=" ")
              alleleplot_sorted <- produceAllelePlot(center_data = current_data, phenotype = phenotype, title = title, isSorted = TRUE, center_name = current_name,z_score_plots = TRUE)
              print(alleleplot_sorted)
            }
            
            
            
            ######CODE TO MODIFY/MERGE DATA FRAMES TO CREATE FINAL COMBINED ALLELE PLOT#################
            
            # replace empty cells in allele column with controls
            modified_current <- current_data
            modified_current$allele_symbol <- as.character(modified_current$allele_symbol)
            modified_current$strain_name <- as.character(modified_current$strain_name)
            modified_current$phenotyping_center <- as.character(modified_current$phenotyping_center)
            # replacing blanks w/ "control"
            modified_current$allele_symbol[(modified_current$allele_symbol == "" | is.na(modified_current$allele_symbol)) & modified_current$biological_sample_group == 'control'] <- modified_current$strain_name[(modified_current$allele_symbol == "" | is.na(modified_current$allele_symbol)) & modified_current$biological_sample_group == 'control']
            modified_current$allele_symbol_zygosity <- paste(modified_current$phenotyping_center,modified_current$allele_symbol,modified_current$zygosity,sep=" ")
            
            # get z-scores
            modified_current_controls <- modified_current[modified_current$biological_sample_group == 'control',]
            controlMean <- mean(modified_current_controls[[phenotype]], na.rm = T)
            controlSD <- sd(modified_current_controls[[phenotype]], na.rm = T)
            modified_current$z_score <- (modified_current[[phenotype]]-controlMean) / controlSD
            
            
            # subset modified dataframe
            modified_current <- modified_current[,c("biological_sample_id", "phenotyping_center", "allele_symbol_zygosity", "z_score")]
            
            # add dataframe to combined dataframe
            combined_df <- rbind(combined_df, modified_current)
            
          }
        }
        x <- "allele_symbol_zygosity"
        y <- "z_score"
        # View(modified_current)
        # View(combined_df)
        # plot <- ggplot() +
        #   geom_pointrange(aes(x = reorder(combined_df$allele_symbol_zygosity, combined_df$z_score, FUN = mean,na.rm=TRUE),y = combined_df$z_score, na.rm = TRUE),data= combined_df,colour = '#0000ff',fill = '#ffffff',size = 0.4, shape = 20 ,fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') + 
        #   geom_point(aes(x = combined_df$allele_symbol_zygosity,y = combined_df$z_score, na.rm = TRUE),data=combined_df,shape = 20,colour = '#999999',size = 1.0, alpha = 0.3) +
        #   theme_bw(base_size = 4.0) + labs(title = title) +
        #   geom_hline(yintercept = -1, color  = "red") + 
        #   geom_hline(yintercept = 0, color  = "red") +
        #   geom_hline(yintercept = 1, color  = "red") +
        #   labs(x=x,y=y) +
        #   coord_flip() +
        #   theme(axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"))
        if(plotAllFile){
          plot <- ggplot() +
            geom_pointrange(aes(x = reorder(combined_df$allele_symbol_zygosity, combined_df$z_score, FUN = mean,na.rm=TRUE),y = combined_df$z_score, color = combined_df$phenotyping_center, group = combined_df$phenotyping_center, na.rm = TRUE),data= combined_df,size = 0.4, shape = 20 ,fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') + 
            scale_colour_brewer(guide = guide_legend(),palette = 'Set1') + theme_bw(base_size = 4.0) + labs(title = paste(phenotype," z-scores")) +
            geom_hline(yintercept = -1, color  = "red") + 
            geom_hline(yintercept = 0, color  = "red") +
            geom_hline(yintercept = 1, color  = "red") +
            labs(x=x,y=y,color="phenotyping centers") +
            coord_flip() +
            theme(axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"))
          print(plot)
          print(paste("COMPLETED FULL ZPLOT: ", phenotype,sep=""))
          
        }
        
        if(!plotAllFile){
          dev.off()
        }
        finished_phenotypes <- c(finished_phenotypes, phenotype)
      }
      print(paste("COMPLETED INDIVIDUAL ZPLOT: ", phenotype, sep=""))
    }
  }
  if(plotAllFile){
    dev.off()
  }
}

# IMPC Data Analysis Code
The purpose of this code is to download KOMP metadata relating to a variety of procedures, phenotypes and metadata across a variety of phenotyping centers from the IMPC pipeline

[This](http://www.mousephenotype.org/impress/procedures/7) link can be used to find a list of the different procedures and the phenotypes/metadata that is recorded for them. The procedures all have a three character procedure code that is used in the code to determine what type of data you wish to download.


## Prerequisites 
Python 3,
Gawk, 
a unix emulator if you are on PC, and
R

### Packages

#### Python
Python packages: csv, pandas, numpy, sys, requests, bs4, os, subprocess, time, datetime, re, glob, lxml, and pprint

#### R
R packages:ggplot2, stringr, dplyr, reshape2, lubridate, gdata, readr, tidyr, Hmisc, plyr, zoo, and signal

## Step By Step Instructions For Running The Code

Make sure you run the different scripts as organized here from top to bottom

### First Python Script
1. Edit the directory variable in IMPC_A.py to the directory you want the data to be downloaded to
2. At the bottom of IMPC_A.py change the procedure code in the function DownloadAll to the 3 digit code for the procedure data you would like to download. The codes are in the format AAA. 
3. Run the IMPC_A.py script 

### Unix Code
1. Set your working directory to the location of the raw data that the first script downloaded. This is located within a folder that is titled your procedure code (such as OFD), which is within your originl directory folder you directed the data to be downloaded to.
2. Run the following line of unix code, where AAA would represent whatever procedure code you are using 
```
awk 'NR%2-1{gsub(/,/,";")}1' RS=\" ORS=\" AAA_data_2.csv > AAA_data.csv
```

### Second Python Script 
1. Edit the directory variable in IMPC_B.py to the directory you want the data to be downloaded to
2. At the bottom of IMPC_B.py change the procedure code in the function create_wide_data to the 3 digit code for the procedure data you would like to download. The codes are in the format AAA. 
3. Run the IMPC_B.py script 

### R All Plots Script
1. Change the directory variable to the directory you used in the two python scripts
2. Change the procedure_code variable to the procedure code you used in the two python scripts 
3. Run the R script and it should begin generating plots 

### R Specific Plot Script
There are multiple R files where each one is for a single plot type that is in the AllPlots.r script, if you only want to generate a specific type of plot, follow the same instructions as the R All Plots Script, except for the plot script you want to use.

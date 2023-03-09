rm(list=ls())
# Set wokring directory
setwd("~/Documents/GitHub/popgen_binp29/Data")

# Load required library
library("dplyr")
library("stringr")

# Read txt file, this file has ID and mean date
date_data=read.table(list.files(pattern = "\\.txt$"),header=FALSE)              # OPTIONAL TO TAKE USER INPUT, write if else for the correct format of txt file
# Make sure that is a numeric data
date_data$V2 <- as.numeric(date_data$V2) 

# Read fam file
fam_data=read.table(list.files(pattern = "\\.fam$"),header=FALSE)               # OPTIONAL TO TAKE USER INPUT

# Filter fam data against the ID data from txt
# Get unique IDs
filter_data = unique(date_data[date_data$V1 %in% fam_data$V2,])

# Create time step
#Round up to the hundred, or thousand, depends on the max value of date
round_digit=nchar(max(filter_data$V2))-2
end_time=round(max(filter_data$V2),-round_digit)
start_time=0
time_step=2500
time_vec=seq(start_time,end_time+time_step,time_step)

# Read bim file
bim_data=read.table(list.files(pattern = "\\.bim$"),header=FALSE)
snps=bim_data$V2

# Read ped file
ped_data=read.delim(list.files(pattern = "\\.ped$"),header=FALSE)

col_index_snp=match("rs11807848",snps)+6
col_index_allele=match("rs11807848",snps)

# The snps allele
temp1=ped_data[,col_index_snp]
temp2=temp1[!temp1=="0 0"]

# Retrieve the minor and major
minor_temp=as.character(bim_data$V5[col_index_allele])
major_temp=as.character(bim_data$V6[col_index_allele])

# Convert all the allele into one strings for further counting
allele_string=paste(temp2,collapse=' ')
minor_count=str_count(allele_string,minor_temp)
major_count=str_count(allele_string,major_temp)
maf=minor_count/(minor_count+major_count)




rm(list=ls())
# Set wokring directory
setwd("~/Documents/GitHub/popgen_binp29/Data")

# Load required library
if(!require("dplyr")) install.packages("dplyr")
if(!require("stringr")) install.packages("stringr")
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("egg")) install.packages("egg")

# Read txt file, this file has ID and mean date
date_data=read.table(list.files(pattern = "\\.txt$"),header=FALSE)              # OPTIONAL TO TAKE USER INPUT, write if else for the correct format of txt file
# Make sure that is a numeric data
date_data$V2 <- as.numeric(date_data$V2) 

# Read fam file
fam_data=read.table(list.files(pattern = "\\.fam$"),header=FALSE)               # OPTIONAL TO TAKE USER INPUT

# Filter fam data against the ID data from txt
# Get unique IDs
filter_data = unique(date_data[date_data$V1 %in% fam_data$V2,])

# Read bim file
bim_data=read.table(list.files(pattern = "\\.bim$"),header=FALSE)
snps=bim_data$V2

# Read ped file
ped_data=read.delim(list.files(pattern = "\\.ped$"),header=FALSE)

# Create time step
#Round up to the hundred, or thousand, depends on the max value of date
round_digit=nchar(max(filter_data$V2))-2
end_time=round(max(filter_data$V2),-round_digit)
start_time=0
time_step=2500
time_vec=seq(start_time,end_time+time_step,time_step)

# FUNCTION TO RUN CODE FROM LINE 30 to 36
create_time_vec <- function(step,data){
  round_digit=nchar(max(data$V2))-2
  end_time=round(max(data$V2),-round_digit)
  start_time=0
  time_vec=seq(start_time,end_time+step,step)
  return(time_vec)
}
# EXAMPLE RUN: create_time_vec(step=1500,data=filter_data)

time_data_list=list()
ped_data_list=list()
# Separate data by time_vec
for (i in 1:(length(time_vec)-1)) { 
  # Create time data, add to a time_data list
  time_temp=filter_data[filter_data$V2 >= time_vec[i] & filter_data$V2 < time_vec[i+1], ]
  time_data_list[[i]]=time_temp
  
  # Separate ped data, 
  ped_temp=ped_data[ped_data$V2 %in% time_temp$V1,]
  ped_data_list[[i]]=ped_temp
  #assign(gsub(" ","",paste("data_",as.character(i),"_time")),temp)
  #assign(gsub(" ","",paste("data_",as.character(i))),
  #       ped_data[ped_data$V2 %in% temp$V1,])
}

snp_name=readline("Input: ")

# Ask for user input using readline
col_index_snp=match(snp_name,snps)+6
col_index_allele=match(snp_name,snps)

# The snps allele
temp2=ped_data[,col_index_snp]

# Retrieve the minor and major
minor_temp=as.character(bim_data$V5[col_index_allele])
major_temp=as.character(bim_data$V6[col_index_allele])

# Convert all the allele into one strings for further counting
allele_string=paste(temp2,collapse=' ')
minor_count=str_count(allele_string,minor_temp)
major_count=str_count(allele_string,major_temp)
maf=minor_count/(minor_count+major_count)


minor_count_vec=NULL
major_count_vec=NULL
maf_vec=NULL
for (i in 1:length(time_data_list)){
  # Retrieve ped data
  sep_ped_temp=ped_data_list[[i]][,col_index_snp]
  
  # Look for minor and major allele in bim file
  minor_temp=as.character(bim_data$V5[col_index_allele])
  major_temp=as.character(bim_data$V6[col_index_allele])
  
  # Paste the allele data to a string
  allele_string=paste(sep_ped_temp, collapse=' ')
  
  # Count minor, major and calculate maf
  # Minor count and append the the current count to a vector
  minor_count=str_count(allele_string, minor_temp)
  minor_count_vec = c(minor_count_vec,minor_count)
  
  # Minor count and append the the current count to a vector
  major_count = str_count(allele_string, major_temp)
  major_count_vec = c(major_count_vec, major_count)
  
  # Calculate maf and append the maf to a vector
  maf = minor_count/(minor_count+major_count)
  maf_vec = c(maf_vec, maf) 
}

#Label minor and major allele
base_lab <- c("A","C","G","T")
# Label minor
result_minor <- tryCatch({
  as.numeric(minor_temp)
}, error = function(e) {
  return(NULL)
})
if (is.null(result_minor)){
  minor_lab=base_temp
} else (minor_lab=base_lab[result_minor])

#Label major
result_major <- tryCatch({
  as.numeric(major_temp)
}, error = function(e) {
  return(NULL)
})
if (is.null(result_major)){
  major_lab=base_temp
} else (major_lab=base_lab[result_major])




# Plotting
time_vec_temp <- c(2500,5000,7500,10000,12500,15000)

# Data for stacked bar plot
plot_data_bar <- data.frame(minor_count_vec,major_count_vec,time_vec_temp)
colnames(plot_data_bar) <- c(minor_lab,major_lab,"Year")
plot_data_bar_long <- reshape2::melt(plot_data_bar, id.vars= "Year", variable.name = "Allele", value.name = "Count")


plot_data_line <- data.frame(maf_vec,time_vec_temp)
colnames(plot_data_line) <- c("MAF","Year")

install.packages("egg")
library(egg)


egg::ggarrange(plot1, plot2, heights = c(0.25, 0.4))



install.packages("shiny")
library(shiny)

plot1 <- ggplot()+
  geom_point(data=plot_data_line, aes(x = Year, y = MAF))+
  geom_line(data=plot_data_line, aes(x = Year, y = MAF))+
  theme_classic()+
  scale_x_reverse(expand = expansion(add = c(1625, 1625))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.01)), limits = c(0,1.1)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + #remove the y axis label
  theme(axis.ticks.length.x = unit(0, "cm")) + #set the length of the divider of value in the x axis to 0
  theme(plot.margin = unit(c(0,0,0,0), "cm")) #expand the plot so when ggarange there is no gap in between two plots

# Part of the trickiness here is that padding works differently 
# for points/lines vs. bars, which have width around the data point. 
# So if we want our points to be aligned with our bars, we need to tweak that.
# To estimate what this should be, I note that the Year resolution is 2500, the
# bars will be 90% of that width (2,250), so each side of the bar extends 1,125 
# beyond the range that the points need. So I add 1,125 to the point padding.

plot2 <- ggplot()+
  geom_bar(data=plot_data_bar_long,stat = "identity", aes(x=Year,y=Count,fill=Allele)) + 
  scale_fill_manual(values = c("C"="green","T"="blue")) +
  scale_x_reverse(expand = expansion(add = c(500, 500))) + #expand the bar and also reverse the x axis value
  scale_y_reverse(expand = expansion(mult = c(0.05, 0))) + #expand the y axis upward and downward so the x-axis kinda disappear when merge two plot
  theme_classic() + theme(legend.position = "none") + #Remove legend
  theme(plot.margin = unit(c(0,0,0,0), "cm")) #Expand the plot so when merge it seamless

egg::ggarrange(plot1, plot2, heights = c(0.25, 0.4))

rm() 

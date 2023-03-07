setwd("~/Documents/GitHub/popgen_binp29/Data")
install.packages("readxl")
install.packages("dplyr")
library("dplyr")
library("readxl")

# Load the excel file (main data)
maindata <- read_excel("V50.xlsx")

# Select only the ID and mean date
filter1_maindata = maindata %>%
  select(c(`Master ID`,`Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]`))

# Load the .fam data (has the ID)
famdata <- read.delim(file="DataS1.fam", sep =" ", header=FALSE)
colnames(famdata) <- c("Country","ID","FatherID","MotherID","Sex","Phenotype")

# Filter famdata against the maindata, to get only the unique IDs presented in both files
filter2_maindata = unique(filter1_maindata[filter1_maindata$`Master ID` %in%  famdata$ID,])
colnames(filter2_maindata) <- c("MasterID","DateMean")

# Sort maindata
sort_maindata <- filter2_maindata[order(filter2_maindata$DateMean),]


data1 <-sort_maindata[sort_maindata$DateMean>10000,]
# Seperate data by time, each step is 500
end_time=(max(filter2_maindata$DateMean))
start_time=min(filter2_maindata$DateMean)
time_step=500
timeseq=seq(start_time,end_time+500,time_step)

# Generate different data base on different time step
for (i in 1:(length(timeseq)-1)) { 
  # paste generate a variable name
  # gsub remove the space in the variable name
  # assign assign the variable name with the data 
  # data was filter base on the timeseq
  assign(gsub(" ","",paste("data_",as.character(i))),filter2_maindata[filter2_maindata$DateMean>=timeseq[i] & filter2_maindata$DateMean < timeseq[i+1], ])
}

rm(list=ls())



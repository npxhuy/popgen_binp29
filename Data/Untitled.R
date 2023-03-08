setwd("~/Documents/GitHub/popgen_binp29/Data")
# install.packages("readxl")
# install.packages("dplyr")
library("dplyr")
library("readxl")
rm(list=ls())
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

# Seperate data by time, each step is 500
end_time=(max(filter2_maindata$DateMean))
start_time=min(filter2_maindata$DateMean)
time_step=1000
timeseq0=seq(start_time-1,end_time+1000,time_step)

# Seperate timeseq into 3 different time step
timeseq1=timeseq0[timeseq0<=8000]
timeseq2_1=timeseq0[timeseq0>8000 & timeseq0<=11000]
timeseq2_2=timeseq2_1[c(FALSE,FALSE,TRUE)] # Time step = 2000
timeseq3=max(timeseq0) # Take the rest

# Overall timeseq
timeseq=c(timeseq1,timeseq2_2,timeseq3)

# Read the ped file
pedfile=read.delim("DataS1.ped",header=FALSE)


# Generate different data base on different time step
for (i in 1:(length(timeseq)-1)) { 
  # paste generate a variable name
  # gsub remove the space in the variable name
  # assign assign the variable name with the data 
  # data was filter base on the timeseq
  temp=filter2_maindata[filter2_maindata$DateMean>=timeseq[i] & filter2_maindata$DateMean < timeseq[i+1], ]
  maxx=max(temp$DateMean)
  minn=min(temp$DateMean)
  assign(gsub(" ","",paste("data_",as.character(minn),"_",as.character(maxx))),
    pedfile[pedfile$V2 %in% temp$MasterID,])
  #write.csv(mydata,file=gsub(" ","",paste("data_",as.character(i))))
}



pedfile %>%
  count(V8)

assign(gsub(" ","",paste("data_",as.character(i))),
       filter2_maindata[filter2_maindata$DateMean>=timeseq[i] & filter2_maindata$DateMean < timeseq[i+1], ])






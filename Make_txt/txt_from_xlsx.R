# This R code will take the excel file that has the ID and the data of the time and put them into txt file
# Highly depends on the excel file column names
# Remember to set working directory to the same directory that you have the .xlsx file
install.packages(c("dplyr","readxl"))
library("dplyr")
library("readxl")
maindata <- read_excel("V50.xlsx") # Changing V50.xlsx to your own xlsx name "
maindata = maindata %>%
  select(c(`Master ID`,`Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]`)) # Change code in this line to the column name of ID and Date.
colnames(maindata) <- NULL
write.table(maindata,file="ID_and_date.txt",row.names = FALSE,quote=FALSE)
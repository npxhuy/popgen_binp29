# This R file is to take the excel file that has the ID and the data of the time and put them into txt file
# Highly depends on the excel file column names
# Changing code in line 8 as your excel file
# Change code in line 10 to the column name of ID and Date, or you can just provide a txt file that has the correct format
install.packages(c("dplyr","readxl"))
library("dplyr")
library("readxl")
maindata <- read_excel("V50.xlsx")
maindata = maindata %>%
  select(c(`Master ID`,`Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]`))
colnames(maindata) <- c("MasterID","DateMean")
write.table(maindata,file="ID_and_date.txt",row.names = FALSE,quote=FALSE)
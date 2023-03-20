SNPulse is an application that will produce an minor allele frequency plot overtime for a certain SNP.

This readme contains the explanation for files type, steps that involving in the population genetics project of the BINP29 course.

# USAGE OF APPLICATION
- Clone the application **SNPulse.R** and **readme.md** file to your computer, run **SNPulse.R** will **R Studio**.
- If you did not have sample data, they could be cloned from folder Data on github.
- Upload the required files, fill in the parameters and hit submit button.
- The generated speed highly depends on the size of data, be patient! 
- The information of files and parameters is written below. The download button can be used to download the plot in high definition.

# REQUIRED FILES
## *.bim*, *.fam* *.bed* & *.ped* file
### ***.bim* file**

Extended variant information file accompanying a .bed binary genotype table.
A text file with no header line, and one line per variant with the following six fields:
1. Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
2. Variant identifier
3. Position in morgans or centimorgans (safe to use dummy value of '0')
4. Base-pair coordinate (1-based; limited to 2<sup>31</sup>-2)
5. Allele 1 (corresponding to clear bits in .bed; usually minor)
6. Allele 2 (corresponding to set bits in .bed; usually major)

Allele codes can contain more than one character. Variants with negative bp coordinates are ignored by PLINK.

### ***.fam* file**
Sample information file accompanying a .bed binary genotype table. (--make-just-fam can be used to update just this file.) Also generated by "--recode lgen" and "--recode rlist".

A text file with no header line, and one line per sample with the following six fields:

1. Family ID ('FID')
2. Within-family ID ('IID'; cannot be '0')
3. Within-family ID of father ('0' if father isn't in dataset)
4. Within-family ID of mother ('0' if mother isn't in dataset)
5. Sex code ('1' = male, '2' = female, '0' = unknown)
6. Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)


### ***.bed* & *.ped* file**
Primary representation of genotype calls at biallelic variants. Must be accompanied by .bim and .fam files.

This .bed file is a binary file, that can not be used directly. Hence it will be used to generated another readable file call .ped file. *PLINK* will be used to generated it.

To use plink, we could either: 
1. Download PLINK, choose the correct operating system, follow the [website instruction](https://www.cog-genomics.org/plink/1.9/) to download

After extraction, copy the *plink* file to the folder that has the data. To use the *plink* file, we run to *plink* file as ./plink --yourcommand

2. Install plink through anaconda, which needed anaconda to be install first, read more about installation of conda in through this [link](https://docs.conda.io/en/latest/miniconda.html).

In this case, we will use the second option.
```sh
# Installing plink, this version is 1.90
conda install plink
# Convert the DataS1.bed file to ped file.
plink --bfile DataS1 --recode --tab --out DataS1
```
This will generate a *.ped* file, a *.map* file that has similar information as *.bim* file but without the allele information.

The *.ped* file contains no header line, and one line per sample with 2V+6 fields where V is the number of variants. The first six fields are the same as those in a *.fam* file. The seventh and eighth fields are allele calls for the first variant in the .map file ('0' = no call); the 9th and 10th are allele calls for the second variant; and so on.

This application required the *.ped* file that was just generated, *.bim* and *.fam* file to run.

## *.txt* file contains ID & Mean Date
The information about mean date was not inclued in the any of the previous file. In this case the mean date data was obtain through an Excel file. 

From that Excel file, the data of the ID of individuals and their mean date will be extracted into a single *.txt* file, with two columns, the first column contains the IDs and the second column contains their mean date.

The ID in this Excel file should also appear in the second column in the *.fam* file

Considering the following code that extract the ID & date information from your own Excel data:

```r
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
```
The R code file and example xlsx file are in the Make_txt folder on github
# OTHERS REQUIREMENTS
## Time step
A positive interger number that indicates the time step to generated the time vector. The default time step is set as 2500. The application will take the absolute of the input time step, so negative number will be converted to positive.

## SNP
A single name of a snp. The input snp should be from one of the snp from the *.bim* file in the second column.

# METHODOLOGY
The process can be summed up as: 
1. Retrieve data from input files & parameters
2. Make data for plotting
3. Produce plot
Specifically, the the process included the following minor steps:
1. Filter *.ped* file with input SNP
2. Filter *.txt* file with *.fam* file
3. Construct time vector from the input time step
4. From 3 data above, separate the *.ped* file base on time vector
5. Calculate the MAF and couting allele
6. Plot point-to-point graph of MAF
7. Plot stacked bar plot of allele count
8. Combine those two plot into one plot

# FAQ
1. Can I use other types of files instead of Plink file format for this app?\
-> I am sorry but no.
2. I only have the *.bed*, *.fam* and *.bim* file, what is the *.ped* file?\
-> Use Plink to generate *.ped* file from those 3 files, more information is in section about *.bed* & *.ped* file
3. What is the *.txt* file?\
-> Read section about *.txt* file contains ID & Mean Date
4. No plot was shown, what was wrong?\
-> Did you remember to click the submit button? The submit button should have a blue *aura* when click it.
5. I click the submit button, no plot was shown, what was wrong?\
-> If no error message was display, the application would have been still running. For the example data set I provide in **Data** directory, it took around 2-3 mins to generate on my computer. Please wait!
6. An error was displayed but it did not tell me what was the error?\
-> I am so sorry that you experience an error with my application. Please open a new issue and give me details about your error, I will find a way to fix it as soon as possible.
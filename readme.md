This readme file contains the explanation for files type, codes, steps that involving in the population genetics project of the BINP29 course.

To use plink, we could either: 
1. Download PLINK, choose the correct operating system, follow the website instruction to download
> https://www.cog-genomics.org/plink/1.9/
After extraction, copy the *plink* file to the folder that has the data. To use the *plink* file, we run to *plink* file as ./plink --yourcommand
2. Install plink through anaconda, which needed anaconda to be install first, read more about instsallation of conda in through this [link](https://docs.conda.io/en/latest/miniconda.html).

In this case, we will use the second option.
```sh
# Installing plink, this version is 1.90
conda install plink
# Convert the DataS1.bed file to ped file.
plink --bfile DataS1 --recode --tab --out DataS1
```
This will generate a .ped file, a .map file that has similar information as .bim file but without the allele information.

# Data explanation
## .bim files (PLINK extended MAP file)
Extended variant information file accompanying a .bed binary genotype table.
A text file with no header line, and one line per variant with the following six fields:
1. Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
2. Variant identifier
3. Position in morgans or centimorgans (safe to use dummy value of '0')
4. Base-pair coordinate (1-based; limited to 2<sup>31</sup>-2)
5. Allele 1 (corresponding to clear bits in .bed; usually minor)
6. Allele 2 (corresponding to set bits in .bed; usually major)

Allele codes can contain more than one character. Variants with negative bp coordinates are ignored by PLINK.

## .fam (PLINK sample information file)
Sample information file accompanying a .bed binary genotype table. (--make-just-fam can be used to update just this file.) Also generated by "--recode lgen" and "--recode rlist".

A text file with no header line, and one line per sample with the following six fields:

1. Family ID ('FID')
2. Within-family ID ('IID'; cannot be '0')
3. Within-family ID of father ('0' if father isn't in dataset)
4. Within-family ID of mother ('0' if mother isn't in dataset)
5. Sex code ('1' = male, '2' = female, '0' = unknown)
6. Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

In this study we do not need to care about the 3 and 4, and also 6 because we are doing population genetics.

## .ped
Contains no header line, and one line per sample with 2V+6 fields where V is the number of variants. The first six fields are the same as those in a .fam file. The seventh and eighth fields are allele calls for the first variant in the .map file ('0' = no call); the 9th and 10th are allele calls for the second variant; and so on.

# Data exploratory

```sh
cat DataS1.ped | cut -f 1 | sort | uniq | wc -l #34 countries

awk '{print NF}' DataS1.ped | sort -nu | tail -n 1 #294464 fields
# excluding first 6 fields, we have a total of 294458 variants

# but only 147235 fields, hmmm
cat DataS1.bim | wc -l #147229
# ahh because we have a total of 147229 variants, plus 6 fields of information we have 147235 fields


cat DataS1.ped | cut -f 147232


```


Could it be plink and freq to have the minor allele freq and plot against it????

Calculate minor allele freq can we use the plink or have to manually write it??

The plot show the minor allele freq overtime, and I assume the data set is one year only??? I guess

What is 1,2,3,4 stand for, I know 1 is A, 3 is C, what about 2 and 4??


```sh
#Cut the ID from the .fam file for filter in the excel file.
cat DataS1.fam | cut -d " " -f2 >> ID 

```
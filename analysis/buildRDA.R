############################################
### Creating the .rda file from raw data ###
############################################

setwd("G:/My Drive/Side projects/Lake herring colonization model/InlandCisco/")

### Create genind file with populations specified
library(adegenet)
library(strataG)

dat.usats <- import2genind("./extData/final_usat.gen", ncode=3)
pops <- c(rep("Big Manistique Lake",33),
          rep("Browns Lake",27),
          rep("Elk Lake",28),
          rep("Glen Lake",27),
          rep("Green Lake",26),
          rep("Harwood Lake",19),
          rep("Howard Lake",30),
          rep("Lime Lake",31),
          rep("Murray Lake",23),
          rep("North Lake Leelanau",30),
          rep("North Sand Lake",25),
          rep("Northern Lake Huron",26),
          rep("Little Traverse Bay",14),
          rep("Walloon Lake",23))
dat.usats@pop <- as.factor(pops) ## Assign populations to each individual, stored within genind object

### Read in latitude and longitude data
latLongs <- read.delim("./extData/latLongs.txt", sep="\t", header = FALSE)

### Read in linkage disequilibrium results from GenePop
ld <- read.delim("./extData/ldResults.txt", sep="\t", header = TRUE)

### Read in distance to Great Lake matrix
distToGL <- read.delim("./extData/distToGL.txt", sep="\t", header = FALSE)

### Read in individual haplotypes
indHaps <- read.delim("./extData/indHaps.csv", sep="\t", header = TRUE)

### Read in Arlequin mtDNA file
datDNA <- read.arlequin("./extData/final_mtDNA.arp")
tmp1 <- expandHaplotypes(datDNA)
datDNA <- labelHaplotypes(tmp1)
datDNA <- datDNA[[1]]

### Save as a .rda file
save(dat.usats, latLongs, ld, distToGL, indHaps, datDNA, file = "./data/InlandCisco.rda")

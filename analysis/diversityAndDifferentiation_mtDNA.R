############################################################
#####  Genetic diversity and differentiation - mtDNA   #####
############################################################

##### Load packages #####
library(strataG)
library(InlandCisco)
data(InlandCisco)

splitDNA <- strataSplit(datDNA)

#Nucleotide Diversity - mtDNA
nucdiv <- sapply(as.data.frame(sapply(splitDNA, nucleotideDiversity)), mean)

#Number of haplotypes (num.alleles) and haplotype diversity (heterozygosity)
summary(datDNA)$strata.smry


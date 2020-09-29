#######################################################################
#####   Checking common population genetic analysis assumptions   #####
#######################################################################

##### Load libraries and data #####
library(PopGenReport)
library(adegenet)
library(tidyverse)
library(qvalue)
library(glHerring)
data(glHerring)

######        Run PopGenReport to obtain allele counts, Hardy-Weinberg tests, Fst table,          ######
######  allele distributions by population, null allele frequencies, allelic richness estimates   ######
######             and differentiation stats (Nei's Gst, Hedrick's Gst, Jost's D)                 ######

popgenreport(dat.usats, mk.counts = TRUE, mk.map = FALSE, mk.locihz = FALSE,
             mk.hwe = TRUE, mk.fst = TRUE, mk.gd.smouse = FALSE,
             mk.gd.kosman = FALSE, mk.pcoa = FALSE, mk.spautocor = FALSE,
             mk.allele.dist = TRUE, mk.null.all = TRUE, mk.allel.rich = TRUE,
             mk.differ.stats = TRUE, mk.custom = FALSE, fname = "PopGenReport_Herring",
             foldername = "PopGenReport_Herring", path.pgr = getwd(), mk.Rcode = FALSE,
             mk.complete = FALSE, mk.pdf = TRUE)


##### FDR evaluation of Hardy-Weinberg proportions ######
hwe <- read.csv("./PopGenReport_Herring/PopGenReport_Herring-HWE_by_locus_location.csv")
hwe.vector <- c(hwe$Sfo8.Lower, hwe$Sfo8.Upper, hwe$Bwf.1, hwe$Bwf.2, hwe$C2.157, hwe$Cocl23, hwe$Cisco.59,
                hwe$Cisco.90, hwe$Cisco.126, hwe$Cisco.181, hwe$Cisco.200, hwe$Cisco.183, hwe$Cisco.179 , hwe$Cisco.106)

hwe.p.corr <- p.adjust(hwe.vector, method = "fdr")
hwe.p.corr[hwe.p.corr < 0.05]

##### FDR evaluation of linkage disequilibrium results imported from Genepop #####
ld.vec <- as.numeric(levels(ld$P.Value))[ld$P.Value]
ld.vec <- ld.vec %>% na.omit()

qval <- qvalue(ld.vec)$qvalues
min(qval)

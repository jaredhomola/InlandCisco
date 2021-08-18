#Script to run hABC model for lake herring
#Usage
#Rscript lake_loop.R <node> <nreps> <who> <group>
#Rscript lake_loop.R 1 10000 JDR N


#Command Line arguments for cluster sims
args <- commandArgs(TRUE)
i <- as.numeric(args[1])    #The "node" for the simulation, one per array job, label for outfile
nreps <- as.numeric(args[2])     #The number of replicates to simulate
lake <- as.numeric(args[3])	#
if(length(args) == 0){
    i <- 1
    nreps <- 10
    lake <- 1
 }

#Run a simulation for all N lakes in the set we're considering
## Check
#codedir = "~/Google Drive/MSU/Students/Homola/Herring/Lake Herring hABC model/development/hierarchical"
#outdir = "~/Desktop/SIMS"

#Jared's directories
#codedir = "G:/My Drive/Side projects/Lake herring colonization model/Lake Herring hABC model/development/"
#outdir = "C:/Users/HP/Desktop/SIMS"
#simdir = "C:/Users/HP/Desktop/test_tmp"

#Cluster directories
codedir = "/mnt/research/ABC/herring/SRC"
outdir = "/mnt/research/ABC/herring/OUT"
#outdir = "/mnt/gs18/scratch/users/homolaj1"
#simdir = system("echo $TMPDIR", intern = TRUE)
simdir = "/mnt/gs18/scratch/users/homolaj1/herringtmp"


setwd(codedir)
source("fxns.R")
laketab <- read.csv("laketable.csv", header = TRUE)

setwd(outdir)

library(strataG)

outfile <- paste0("LakeHerr_", lake, "_", i,".csv")
control <- read.csv("../SRC/control.csv", header = TRUE)

for(r in 1:nreps) {
	setwd(codedir)
	parms <- drawParms("control.csv", n.lakes = length(lakeIDnum))
	
	if(parms$Tadmix >= parms$Tgl){
	  parms$Tadmix = runif(1, control$min[control$parameter == "Tadmix"], parms$Tgl)
	} 

	out <- parms

	setwd(simdir)
#	for(lake in 1:length(lake.n.SSR)) {
		stats <- runrep(parms = parms, lake.n.SSR = laketab$nSSR[lake], lake.n.DNA = laketab$nDNA[lake], ind.num = lake, lakeID = laketab$lakeID[lake])
		out <- cbind(out, stats)
#	}

	setwd(outdir)
	if(!file.exists(outfile)) {
        write.table(out, outfile, sep = ",", quote = FALSE, row.names = FALSE)
    } else {
        write.table(out, outfile, quote = FALSE, row.names = FALSE, sep=",", append = TRUE, col.names = FALSE)
    }
    
    rm(out, parms)

}

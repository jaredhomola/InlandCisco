#Function to run the simulations (both markers) and calculate stats

runrep <- function(parms=NULL, lake.n.SSR = NULL, lake.n.DNA = NULL, ind.num = NULL, lakeID = NULL) {
	
  #SSR
	#fscPopInfo()
	#fscPopInfo(pop.size, sample.size, sample.times = 0, growth.rate = 0)
	piSSR <- fscPopInfo(pop.size = c(parms$N.ATL, 
									parms$N.MS, 
									parms$N.GL,
									parms$N.lake, 
									parms$N.STOCK), 
						sample.size = c(0,0,26,lake.n.SSR,0))
	
	#fscLocusParams
	#fscLocusParams(locus.type = c("dna", "msat", "snp"), sequence.length = NULL,
  		#num.loci = NULL, mut.rate = NULL, transition.rate = 1/3,
  		#gsm.param = 0, range.constraint = 0, recomb.rate = 0,
  		#chromosome = NULL, num.chrom = NULL, ploidy = 2)
	locSSR = fscLocusParams(locus.type = "msat", num.chrom = 14, num.loci = 1, mut.rate = parms$muSSR)

	#fscHistEv
	#fscHistEv(num.gen = 0, source.deme = 0, sink.deme = 0,
  		#prop.migrants = 1, new.sink.size = 1, new.sink.growth = 0,
  		#new.mig.mat = 0)

	hist.ev = fscHistEv(num.gen = c(parms$Tstock/parms$G,
									(parms$Tstock+parms$G)/parms$G,
									parms$Tadmix/parms$G,
									(parms$Tlake-parms$G)/parms$G,
									parms$Tlake/parms$G,
									parms$Tgl/parms$G,
									parms$Trefdiv/parms$G), 
						source.deme = c(3,4,2,3,3,2,1), 
						sink.deme =   c(4,2,0,3,2,1,0), 
						prop.migrants = c(parms$pSTOCK,1,parms$pATL,1,1,1,1), 
						new.sink.size = c(1,1,1,parms$Bneck,1,1,1), 
						new.sink.growth = rep(0,7), new.mig.mat = rep(0,7))

	#MTDNA
	piDNA <- piSSR
	piDNA[,1] <- piDNA[,1]/4
	piDNA[,2] <- c(0, 0, 26, lake.n.DNA, 0)

	locDNA <- fscLocusParams(locus.type = "dna", sequence.length = 525, mut.rate = parms$muDNA, ploidy = 1)

	#fastsimcoal
	#fastsimcoal(pop.info, locus.params, mig.rates = NULL, hist.ev = NULL,
		#  label = NULL, quiet = TRUE, delete.files = TRUE, exec = "fsc252",
		#  num.cores = NULL, label.haplotypes = TRUE)

	datSSR <- fastsimcoal(pop.info = piSSR, locus.params = locSSR, hist.ev = hist.ev, label = paste0("HerrSSR-",ind.num,"_",runif(1)), delete.files = TRUE, exec = "fsc26")

	datDNA <- fastsimcoal(pop.info = piDNA, locus.params = locDNA, hist.ev = hist.ev, label = paste0("HerrDNA-",ind.num,"_",runif(1)), delete.files = TRUE, exec = "fsc26")


	#Now calculate some summary statistics for each 
	split_datSSR <- strataSplit(datSSR)
	split_datDNA <- strataSplit(datDNA)

	#Nucleotide Divergence - mtDNA
	Da <- nucleotideDivergence(datDNA)$gene.1$between$dA

	#Nucleotide Diversity - mtDNA      
	nucdiv <- mean(nucleotideDiversity(split_datDNA[[2]]))

	#Fu's Fs - mtDNA
	FuF <- fusFs(split_datDNA[[2]])["gene.1"]

	#Tajima's D - mtDNA
	TajD <- tajimasD(expandHaplotypes(split_datDNA[[2]]))[1,1]

	#Number of Haplotypes
	#nhap <- sum(table(datDNA@data$gene.1[datDNA@data$strata == "Sample 2"]) > 0)
	nhap <- summary(split_datDNA[[2]])$strata.smry[3]

	hapdiv <- summary(split_datDNA[[2]])$strata.smry[5]

	#Allelic Richness - SSR
	AR <- mean(allelicRichness(split_datSSR[[2]]))

	#Number of Alleles - SSR

	#Private Alleles - SSR
	privall <- sum(privateAlleles(datSSR)[,2])

	#Fst - SSR
	Fst <- statFst(datSSR)$result["estimate"]

	#Jost's D - SSR
	JostD <- statJostD(datSSR)$result["estimate"]

	#Expected & Observed Heterozygosity - SSR
	He <- mean(exptdHet(split_datSSR[[2]]))
	#Ho <- mean(obsvdHet(split_datSSR[[2]]))

	#Fis - het based
	#F = 1 - (Ho/He)
	loc_Fis <- 1-(obsvdHet(split_datSSR[[2]])/exptdHet(split_datSSR[[2]]))
	loc_Fis[is.na(loc_Fis)] <- 0
	Fis <- mean(loc_Fis)

	#M-Ratio - SSR
	Mrat <- mean(mRatio(split_datSSR[[2]], rpt.size = 1), na.rm = TRUE)

	#DAPC
	giObj <- gtypes2genind(datSSR, type = "codom")
	dapc.out <- dapc(giObj, var.contrib = TRUE, scale = FALSE, n.pca = 25, n.da = 1)
	dapc.post <- dapc.out$posterior[giObj@pop == "Sample 2",2]
	dapc.var.exp <- dapc.out$var
	dapc.post.mean <- mean(dapc.post)
	dapc.post.var <- var(dapc.post)

	# bind it up 
	stats <- data.frame(mtDNA.Da = Da, 
					mtDNA.pi = nucdiv,
					mtDNA.FuF = FuF,
					mtDNA.TajD = TajD,
					mtDNA.nhap = nhap,
					mtDNA.Hd = hapdiv,
					SSR.AR = AR,
					SSR.priv = privall,
					SSR.Fst = Fst,
					SSR.Fis = Fis,
					SSR.JostD = JostD,
					SSR.He = He,
					SSR.Mrat = Mrat,
					row.names = NULL, 
					dapc.var.exp = dapc.var.exp,
					dapc.post.mean = dapc.post.mean, 
					dapc.post.var = dapc.post.var)
					#dapc.post.min = dapc.post.min, 
					#dapc.post.max = dapc.post.max)

	names(stats) <- paste0(lakeID, ".", names(stats))

	Fucols <- grep("FuF",names(stats))
	Dcols <- grep("TajD",names(stats))
	stats[which(is.na(stats) & names(stats) %in% names(stats)[c(Dcols,Fucols)])] <- 0
	
	stats

}





#Function to draw parameter values for a replicate - draws all parm values...


drawParms <- function(control = NULL, n.lakes = NULL) {

	priors <- read.csv(control, header = TRUE, row.names = 1)

	parms <- c()

	for(p in 1:length(priors[,1])) {
		
#Prior draws for "Fixed" Parameters - these don't vary across simulations, Generation time for instance

		if(priors$type[p] == "fixed") {
			parms[p] <- as.numeric(as.character(priors$min[p]))
			
			#Round if needed
			if(priors$is.int[p] == 1) {
				parms[p] <- round(parms[p])
			}
			names(parms)[p] <- rownames(priors)[p]
		}

#Prior draws for Parameters with flat priors - these are shared in simulations for all lakes
		
		else if(priors$type[p] == "uniform") {
			if(priors$is.cond[p] == 0) {
				min <- as.numeric(as.character(priors$min[p]))
				max <- as.numeric(as.character(priors$max[p]))
				parms[p] <- runif(1, min, max)
				rm(min)
				rm(max)
				names(parms)[p] <- rownames(priors)[p]
				#Round if needed
				if(priors$is.int[p] == 1) {
					parms[p] <- round(parms[p])
				}
			} else if(priors$is.cond[p] == 1) {
				if(length(which(rownames(priors)==as.character(priors$min[p]))) == 0) {
					min <- as.numeric(as.character(priors$min[p]))
					max <- parms[which(rownames(priors)==priors$max[p])]
					parms[p] <- runif(1, min, max)
					rm(min)
					rm(max)
					names(parms)[p] <- rownames(priors)[p]
					#Round if needed
					if(priors$is.int[p] == 1) {
						parms[p] <- round(parms[p])
					}
				} else if(length(which(rownames(priors)==as.character(priors$max[p]))) == 0) {
					min <- parms[which(rownames(priors)==priors$min[p])]
					max <- as.numeric(as.character(priors$min[p]))
					parms[p] <- runif(1, min, max)
					rm(min)
					rm(max)
					names(parms)[p] <- rownames(priors)[p]
					#Round if needed
					if(priors$is.int[p] == 1) {
						parms[p] <- round(parms[p])
					}
				}	
			}

#Prior draws for Hyper Parameters - save #, mean, var in addition to individual parameter values

		} else if(priors$type[p] == "hyper") {
			n.vals <- n.lakes
			
			if(priors$is.cond[p] == 0) {
				min <- as.numeric(as.character(priors$min[p]))
				max <- as.numeric(as.character(priors$max[p]))
				n.parms <- sample(c(1:n.vals),1,replace = FALSE)
				hyp.parms <- runif(n.parms, min = min, max = max)
				redraw <- n.vals - n.parms
				hyp.parms <- c(hyp.parms, sample(hyp.parms, redraw, replace = TRUE))
				hyp.parms <- hyp.parms[sample(c(1:length(hyp.parms)))]

				names(n.parms) <- paste0("n.",rownames(priors)[p])
				mean.parm <- mean(hyp.parms)
				names(mean.parm) <- paste0("mean.", rownames(priors)[p])
				var.parm <- var(hyp.parms)
				names(var.parm) <- paste0("var.", rownames(priors)[p])

				names(hyp.parms) <- paste0("hyp.",rownames(priors)[p], ".", c(1:n.vals))
				
				rm(min)
				rm(max)

				if(priors$is.int[p] == 1) {
					hyp.parms <- round(hyp.parms)
				}
				parms <- c(parms, n.parms, mean.parm, var.parm, hyp.parms)
			} else if(priors$is.cond[p] == 1) {
				if(length(which(rownames(priors)==as.character(priors$min[p]))) == 0) {
					min <- as.numeric(as.character(priors$min[p]))
					max <- parms[which(rownames(priors)==priors$max[p])]
					n.parms <- sample(c(1:n.vals),1,replace = FALSE)
					hyp.parms <- runif(n.parms, min = min, max = max)
					redraw <- n.vals - n.parms
					hyp.parms <- c(hyp.parms, sample(hyp.parms, redraw, replace = TRUE))
					hyp.parms <- hyp.parms[sample(c(1:length(hyp.parms)))]

					names(n.parms) <- paste0("n.",rownames(priors)[p])
					mean.parm <- mean(hyp.parms)
					names(mean.parm) <- paste0("mean.", rownames(priors)[p])
					var.parm <- var(hyp.parms)
					names(var.parm) <- paste0("var.", rownames(priors)[p])

					names(hyp.parms) <- paste0("hyp.",rownames(priors)[p], ".", c(1:n.vals))

					rm(min)
					rm(max)

					if(priors$is.int[p] == 1) {
						hyp.parms <- round(hyp.parms)
					}
					parms <- c(parms, n.parms, mean.parm, var.parm, hyp.parms)
				} else if(length(which(rownames(priors)==as.character(priors$max[p]))) == 0) {
					min <- parms[which(rownames(priors)==priors$min[p])]
					max <- as.numeric(as.character(priors$min[p]))
					n.parms <- sample(c(1:n.vals),1,replace = FALSE)
					hyp.parms <- runif(n.parms, min = min, max = max)
					redraw <- n.vals - n.parms
					hyp.parms <- c(hyp.parms, sample(hyp.parms, redraw, replace = TRUE))
					hyp.parms <- hyp.parms[sample(c(1:length(hyp.parms)))]
					
					names(n.parms) <- paste0("n.",rownames(priors)[p])
					mean.parm <- mean(hyp.parms)
					names(mean.parm) <- paste0("mean.", rownames(priors)[p])
					var.parm <- var(hyp.parms)
					names(var.parm) <- paste0("var.", rownames(priors)[p])

					names(hyp.parms) <- paste0("hyp.",rownames(priors)[p], ".", c(1:n.vals))

					rm(min)
					rm(max)

					if(priors$is.int[p] == 1) {
						hyp.parms <- round(hyp.parms)
					}
					parms <- c(parms, n.parms, mean.parm, var.parm, hyp.parms)
				}
			}		
#Prior draws for "Multiple" Parameters - these are lake specific, non-hyper parameters

		} else if(priors$type[p] == "multiple") {
			if(priors$is.cond[p] == 0) {
				min <- as.numeric(as.character(priors$min[p]))
				max <- as.numeric(as.character(priors$max[p]))
				n.parms <- n.lakes
				mult.parms <- runif(n.parms, min = min, max = max)
				names(mult.parms) <- paste0(rownames(priors)[p],".",c(1:n.parms))
				rm(min)
				rm(max)
				if(priors$is.int[p] == 1) {
					mult.parms <- round(mult.parms)
				}
				parms <- c(parms, mult.parms)
			} else if(priors$is.cond[p] == 1) {
				if(length(which(rownames(priors)==as.character(priors$min[p]))) == 0) {
					min <- as.numeric(as.character(priors$min[p]))
					max <- parms[which(rownames(priors)==priors$max[p])]
					n.parms <- n.lakes
					mult.parms <- runif(n.parms, min = min, max = max)
					names(mult.parms) <- paste0(rownames(priors)[p],".",c(1:n.parms))
					rm(min)
					rm(max)
					if(priors$is.int[p] == 1) {
						mult.parms <- round(mult.parms)
					}
					parms <- c(parms, mult.parms)
				} else if(length(which(rownames(priors)==as.character(priors$max[p]))) == 0) {
					min <- parms[which(rownames(priors)==priors$min[p])]
					max <- as.numeric(as.character(priors$min[p]))
					n.parms <- n.lakes
					mult.parms <- runif(n.parms, min = min, max = max)
					names(mult.parms) <- paste0(rownames(priors)[p],".",c(1:n.parms))
					rm(min)
					rm(max)
					if(priors$is.int[p] == 1) {
						mult.parms <- round(mult.parms)
					}
					parms <- c(parms, mult.parms)
				}
			}
		}
	}
		
	parms <- as.data.frame(t(parms))
	parms
}


#Making some figures... herring work
dir = "~/Desktop"

###1 - CV4MS for mnlogistic model selection
##1a - N
cvdat <- read.csv("LakeHerr_N_CV4MS_mnlog0.0005_allclean.csv")
par(mfrow = c(2,3))
for(x in unique(cvdat$true.model)) {
	cols = rep("white",length(unique(cvdat$true.model)))
	cols[x] = "light grey"
	boxplot(cvdat[cvdat$true.model == x,-1], main = paste("N - mnlog - model", x),  ylab = "Posterior Probability", names = paste0("mod",unique(cvdat$true.model)), col = cols)
}
rm(cvdat)

##1b - S
cvdat <- read.csv("LakeHerr_S_CV4MS_mnlog0.0005_allclean.csv")
par(mfrow = c(3,3))
for(x in unique(cvdat$true.model)) {
	cols = rep("white",length(unique(cvdat$true.model)))
	cols[x] = "light grey"
	boxplot(cvdat[cvdat$true.model == x,-1], main = paste("S - mnlog - model", x),  ylab = "Posterior Probability", names = paste0("mod",unique(cvdat$true.model)), col = cols)
}
rm(cvdat)


##################################
##################################
##################################
##################################
##################################


###2 - CV4PE plots for loclinear, neural net
##2a - N loclin
setwd("~/Desktop")
load("LakeHerr_N_CV4PEloclin_mod5.Rws")
par(mfrow = c(5,6))
for(x in 1:dim(CV4PE$true)[2]){
	plot(CV4PE$true[,x], CV4PE$estim$tol0.001[,x], xlim = range(CV4PE$true[,x]), ylim = range(CV4PE$true[,x]), pch = 19, xlab = "Simulated", ylab = "Estimated", main = paste("loclin - N -", names(CV4PE$true)[x]))
	abline(a = 0, b = 1, lwd = 1, col = "red")
	legend("topleft", legend = paste("Pearson's r =", round(cor(CV4PE$true[,x], CV4PE$estim$tol0.001[,x]), 3)), bty = "n")
}

##2b - N nnet
load("LakeHerr_N_CV4PEnnet_mod5.Rws")
par(mfrow = c(5,6))
for(x in 1:dim(CV4PE$true)[2]){
	plot(CV4PE$true[,x], CV4PE$estim$tol0.01[,x], xlim = range(CV4PE$true[,x]), ylim = range(CV4PE$true[,x]), pch = 19, xlab = "Simulated", ylab = "Estimated", main = paste("nnet - N -", names(CV4PE$true)[x]))
	abline(a = 0, b = 1, lwd = 1, col = "red")
	legend("topleft", legend = paste("Pearson's r =", round(cor(CV4PE$true[,x], CV4PE$estim$tol0.01[,x]), 3)), bty = "n")
}


##2c - S loclin
load("LakeHerr_S_CV4PEloclin_mod5.Rws")
par(mfrow = c(6,6))
for(x in 1:dim(CV4PE$true)[2]){
	plot(CV4PE$true[,x], CV4PE$estim$tol0.001[,x], xlim = range(CV4PE$true[,x]), ylim = range(CV4PE$true[,x]), pch = 19, xlab = "Simulated", ylab = "Estimated", main = paste("loclin - S -", names(CV4PE$true)[x]))
	abline(a = 0, b = 1, lwd = 1, col = "red")
	legend("topleft", legend = paste("Pearson's r =", round(cor(CV4PE$true[,x], CV4PE$estim$tol0.001[,x]), 3)), bty = "n", cex = 0.8)
}

##2d - S nnet
load("LakeHerr_S_CV4PEnnet_mod5.Rws")
par(mfrow = c(6,6))
for(x in 1:dim(CV4PE$true)[2]){
	plot(CV4PE$true[,x], CV4PE$estim$tol0.01[,x], xlim = range(CV4PE$true[,x]), ylim = range(CV4PE$true[,x]), pch = 19, xlab = "Simulated", ylab = "Estimated", main = paste("nnet - S -", names(CV4PE$true)[x]))
	abline(a = 0, b = 1, lwd = 1, col = "red")
	legend("topleft", legend = paste("Pearson's r =", round(cor(CV4PE$true[,x], CV4PE$estim$tol0.01[,x]), 3)), bty = "n", cex = 0.8)
}


##################################
##################################
##################################
##################################
##################################


###3 - Parm Posteriors
##3a - N 
setwd("~/Desktop")
library(abc)
load("LakeHerr_N_abcPE_mod5.Rws")
par(mfrow = c(5,6), mar = c(4,4,2,2))
for(x in 1:ncol(parms)) {
	plot(density(parms[,x]), lwd = 0.5, lty = "dotted", main = paste("N mod 5 -", names(parms)[x]), xlim = range(parms[,x]), ylim = c(0, max(density(parms[,x])$y, density(PE_loclin_0.001$adj.values[,x])$y, density(PE_loclin_0.005$adj.values[,x])$y)), xlab = names(parms)[x],  ylab = "Posterior Density")
	lines(density(PE_loclin_0.005$adj.values[,x]), lwd = 1)
	lines(density(PE_loclin_0.001$adj.values[,x]), lwd = 2)
	lines(density(PE_nnet_0.05$adj.values[,x]), lwd = 1, col = "red")
	lines(density(PE_nnet_0.01$adj.values[,x]), lwd = 2, col = "red")
}
plot.new()
legend("center", lty = c("dotted", "solid", "solid", "solid", "solid"), lwd = c(0.5, 1, 2, 1, 2), col = c(rep("black",3), rep("red",2)), legend = c("Prior", "LocLin, Tol = 0.005", "LocLin, Tol = 0.001", "NNet, Tol = 0.05", "NNet, Tol = 0.01"), bty = "n", cex = 1.5)


##3b - S 
setwd("~/Desktop")
library(abc)
load("LakeHerr_S_abcPE_mod5.Rws")
par(mfrow = c(6,6), mar = c(4,4,2,2))
for(x in 1:ncol(parms)) {
	plot(density(parms[,x]), lwd = 0.5, lty = "dotted", main = paste("S mod5 -", names(parms)[x]), xlim = range(parms[,x]), ylim = c(0, max(density(parms[,x])$y, density(PE_loclin_0.001$adj.values[,x])$y, density(PE_loclin_0.005$adj.values[,x])$y)), xlab = names(parms)[x],  ylab = "Posterior Density")
	lines(density(PE_loclin_0.005$adj.values[,x]), lwd = 1)
	lines(density(PE_loclin_0.001$adj.values[,x]), lwd = 2)
	lines(density(PE_nnet_0.05$adj.values[,x]), lwd = 1, col = "red")
	lines(density(PE_nnet_0.01$adj.values[,x]), lwd = 2, col = "red")
}
plot.new()
legend("center", lty = c("dotted", "solid", "solid", "solid", "solid"), lwd = c(0.5, 1, 2, 1, 2), col = c(rep("black",3), rep("red",2)), legend = c("Prior", "LocLin, Tol = 0.005", "LocLin, Tol = 0.001", "NNet, Tol = 0.05", "NNet, Tol = 0.01"), bty = "n", cex = 1.5)

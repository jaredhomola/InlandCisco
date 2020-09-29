############################################################
#####      Genetic diversity and differentiation       #####
############################################################

##### Load packages #####
library(hierfstat)
library(adegenet)
library(mmod)
library(tidyverse)
library(ggrepel)
library(gridExtra)
library(genepop)
library(glHerring)
data(glHerring)

##### Calculate diversity measures #####
gd <- basic.stats(dat.usats)

he <- as.data.frame.matrix(gd$Hs)
he.pop <- as.data.frame(colMeans(he))

fis <- as.data.frame(gd$Fis)
fis[is.na(fis)] <- 0
fis.pop <- colMeans(fis)

ar <- as.data.frame(allelic.richness(dat.usats))[,2:15]
ar.pop <- as.data.frame(colMeans(ar))
names(ar.pop) <- names(he.pop)

##### Assemble master dataframe for regression analyses #####
gd.df <- cbind(ar.pop, he.pop, latLongs$V3, distToGL$V1)
rownames(gd.df) <- rownames(he.pop)
names(gd.df) <- c("AR", "He", "Lat", "distToGL")

##### Calculate pairwise genetic differentiation #####
gst.nei.paired <- pairwise_Gst_Nei(dat.usats)
gst.hedrick.paired <- pairwise_Gst_Hedrick(dat.usats)

Gst_Nei(dat.usats)
Gst_Hedrick(dat.usats)

##### Exact test for differentiation ######
test_diff("./extData/final_usat.gen",
          pairs = TRUE,
          outputFile = 'herring.txt.GE')


#write.csv(as.matrix(gst.nei.paired), "./nei.paired.csv")
#write.csv(as.matrix(gst.hedrick.paired), "./hedrick.paired.csv")

##### Genetic diverity vs. latitude #####
## Models
Lat.He.mod <- lm(He ~ Lat, data = gd.df[-c(12,13),])
summary(Lat.He.mod)

GL.He.mod <- lm(He ~ distToGL, data = gd.df[-c(12,13),])
summary(GL.He.mod)

Lat.AR.mod <- lm(AR ~ Lat, data = gd.df[-c(12,13),])
summary(Lat.AR.mod)

GL.AR.mod <- lm(AR ~ distToGL, data = gd.df[-c(12,13),])
summary(GL.AR.mod)


## Plotting
p1 <- ggplot(gd.df[-c(12,13),], aes(x = Lat, y = He)) +
  geom_smooth(method='lm') +
  geom_point() +
  geom_text_repel(data=gd.df[-c(12,13),],
                  #fontface = "bold",
                  aes(label=rownames(gd.df[-c(12,13),])),
                  size=5) +
  ylab("Expected heterozygosity") +
#  xlab("Latitude") +
  ylim(0.4, 0.8) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  annotate("text", x = 45.5, y = 0.42, size = 6,
           label = "paste(Adj., italic(R) ^ 2, \" = 0.670\")", parse = TRUE) +
  annotate("text", x = 42, y = 0.8, size = 10,
           label = "bold(A.)", parse = TRUE) +
  NULL

p2 <- ggplot(gd.df[-c(12,13),], aes(x = distToGL, y = He)) +
  geom_smooth(method='lm') +
  geom_point() +
  geom_text_repel(data=gd.df[-c(12,13),],
                  #fontface = "bold",
                  aes(label=rownames(gd.df[-c(12,13),])),
                  size=5) +
#  ylab("Expected Heterozygosity") +
#  xlab("Distance to Nearest Great Lake") +
  ylim(0.4, 0.8) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  annotate("text", x = 20, y = 0.42, size = 6,
           label = "paste(Adj., italic(R) ^ 2, \" = 0.859\")", parse = TRUE) +
  annotate("text", x = 110, y = 0.8, size = 10,
           label = "bold(B.)", parse = TRUE) +
  NULL

p3 <- ggplot(gd.df[-c(12,13),], aes(x = Lat, y = AR)) +
  geom_smooth(method='lm') +
  geom_point() +
  geom_text_repel(data=gd.df[-c(12,13),],
                  #fontface = "bold",
                  aes(label=rownames(gd.df[-c(12,13),])),
                  size=5) +
  ylab("Allelic richness") +
  xlab("Latitude") +
  ylim(0, 10) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  annotate("text", x = 45.5, y = 0.8, size = 6,
           label = "paste(Adj., italic(R) ^ 2, \" = 0.435\")", parse = TRUE) +
  annotate("text", x = 42, y = 9.5, size = 10,
           label = "bold(C.)", parse = TRUE) +
  NULL

p4 <- ggplot(gd.df[-c(12,13),], aes(x = distToGL, y = AR)) +
  geom_smooth(method='lm') +
  geom_point() +
  geom_text_repel(data=gd.df[-c(12,13),],
                  #fontface = "bold",
                  aes(label=rownames(gd.df[-c(12,13),])),
                  size=5) +
#  ylab("Allelic richness") +
  xlab("Distance to Nearest Great Lake (km)") +
  ylim(0, 10) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  annotate("text", x = 20, y = 0.8, size = 6,
           label = "paste(Adj., italic(R) ^ 2, \" = 0.678\")", parse = TRUE) +
  annotate("text", x = 110, y = 9.5, size = 10,
           label = "bold(D.)", parse = TRUE) +
  NULL

grid.arrange(p1, p2, p3, p4, nrow = 2)

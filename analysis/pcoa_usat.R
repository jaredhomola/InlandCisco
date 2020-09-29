################################################################
######                  Lake Herring PCoA                 ######
################################################################

##### Load data and libraries #####
library(adegenet)
library(ggrepel)
library(tidyverse)
library(PopGenReport)
library(glHerring)
library(cowplot)
library(magick)
data(glHerring)

### Map file
map_file <- system.file("extData", "MichiganMap-01.png", package = "cowplot")

######## Smouse & Peakall Genetic Distance PCoA ##########
sp <- gd.smouse(dat.usats, verbose = TRUE)
pcoa.sp <- dudi.pco(sp, scannf=FALSE, nf=3)
summary(pcoa.sp)

pcoa.sp.df <- cbind(pcoa.sp$li, dat.usats@pop)
pcoa.sp.df[,1:2] <- pcoa.sp.df[,1:2] * -1
names(pcoa.sp.df) <- c("A1", "A2", "A3", "pop")
centroids.sp <- aggregate(cbind(A1, A2)~pop, pcoa.sp.df, mean)

pcoa.plot <- ggplot(pcoa.sp.df, aes(x = A1, y = A2, color = pop)) +
  geom_point(alpha = 0.25) +
  geom_point(data = centroids.sp, size = 5) +
  geom_text_repel(data=centroids.sp, fontface = "bold", aes(label=pop), size=6) +
  #ylim(-15, 15) +
  #xlim(-13.5, 20) +
  ylab("Axis 2 (8.0%)") +
  xlab("Axis 1 (9.4%)") +
  labs(color = "Population") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  NULL

ggdraw(pcoa.plot) +
  draw_image("./extData/MichiganMap.png", x = .32, y = .27, scale = 0.5)




ggplot(pcoa.sp.df, aes(x = A1, y = A2, group = pop)) +
  stat_ellipse() +
  #geom_point(alpha = 0.25) +
  geom_point(data = centroids.sp, size = 5) +
  geom_text_repel(data=centroids.sp, fontface = "bold", aes(label=pop), size=6) +
  #ylim(-15, 15) +
  #xlim(-13.5, 20) +
  ylab("Axis 2 (8.0%)") +
  xlab("Axis 1 (9.4%)") +
  labs(color = "Population") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  NULL

## Load data and libraries
# Jared's setwd
library(ggmap)
library(tidyverse)
library(ggsn)
library(ggrepel)

# Set Jared's Google API Key
register_google(key = "AIzaSyDNBsBTypZFWTPTKfyyE9Ul_6lPg47idzo")

# Create reference tibbles
latLong <- read.delim("G:/My Drive/Side projects/Lake herring colonization model/Data analysis/manuscript_analyses/glHerring/extData/latLongs.txt",
                      sep = "\t",
                      header = FALSE) %>%
  as_tibble() %>%
  rename(Lake = V4)


# Calculate median stocking estimates
#### Can we believe our model regarding differential stocking rates? #####

library(Hmisc)
library(corrplot)
setwd("G:/My Drive/Side projects/Lake herring colonization model/Lake Herring hABC model/")

### Load data
N.PE <- read_csv("./Plots_v1/herringFigs/PE_loclin_005.N.adjvalues.csv")
S.PE <- read_csv("./Plots_v1/herringFigs/PE_loclin_005.S.adjvalues.csv")

### Median estimated stocking rates
N.median <- N.PE %>%
  select(starts_with("pSTOCK")) %>%
  rename(`Big Manistique Lake` = pSTOCK.1,
         `Elk Lake` = pSTOCK.2,
         `Glen Lake` = pSTOCK.3,
         `North Lake Leelanau` = pSTOCK.4,
         `Walloon Lake` = pSTOCK.5
  ) %>%
  pivot_longer(cols = everything(),
               names_to = "Lake",
               values_to = "value") %>%
  group_by(Lake) %>%
  dplyr::summarize(median.sim = median(value))

all.median <- S.PE %>%
  select(starts_with("pSTOCK")) %>%
  rename(`Browns Lake` = pSTOCK.1,
         `Green Lake` = pSTOCK.2,
         `Harwood Lake` = pSTOCK.3,
         `Howard Lake` = pSTOCK.4,
         `Lime Lake` = pSTOCK.5,
         `Murray Lake` = pSTOCK.6,
         `North Sand Lake` = pSTOCK.7
  ) %>%
  pivot_longer(cols = everything(),
               names_to = "Lake",
               values_to = "value") %>%
  group_by(Lake) %>%
  dplyr::summarize(median.sim = median(value)) %>%
  bind_rows(N.median) %>%
  rename(simStock.median = median.sim)

plotting.tib <- left_join(all.median, latLong) %>%
  rename(Est_Admix.Stock = simStock.median)

source.tib <- latLong %>% filter(Lake == "Northern Lake Huron")

### Create bar chart


### Create Map
ggmap(get_googlemap(
  center = c(lon = -84.7633, lat = 44.1439),
  zoom = 7,
  scale = 2,
  maptype = 'satellite',
  color = 'color'
)) +
  geom_point(data = plotting.tib,
             aes(x = V2, y = V3, color = Est_Admix.Stock),
             size = 6) +
  #geom_point(data = source.tib,
  #           aes(x = V2, y = V3),
  #           color = "white",
  #           size = 7,
  #           shape = 17) +
  scale_color_gradient(low="blue", high="red") +
  geom_label_repel(
    aes(V2, V3, label = Lake),
    data = plotting.tib,
    size = 4,
    point.padding = 0.6,
    segment.color = 'white',
    segment.size = 1) +
  #geom_label_repel(
  #  aes(x = V2, y = V3, label = "Source"),
  #  data = source.tib,
  #  size = 4,
  #  point.padding = 0.6,
  #  segment.color = 'white',
  #  segment.size = 1) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Lake Herring Sites")

library(readxl)
library(tidyverse)

dat <- read_xlsx("../../Historic Whitefish Stocking Database.xlsx")

targets <- c("Manistique",
            "Browns",
            "Elk",
            "Glen",
            "Green",
            "Harwood",
            "Howard",
            "Lime",
            "Murray",
            "Leelanau",
            "North Sand",
            "Walloon")

waters <- dat %>% filter(str_detect(Waters_Name, paste(c("Manistique",
                                       "Browns",
                                       "Elk",
                                       "Glen",
                                       "Green",
                                       "Harwood",
                                       "Howard",
                                       "Lime",
                                       "Murray",
                                       "Leelanau",
                                       "North Sand",
                                       "Walloon"), collapse = '|')))

sites <- dat %>% filter(str_detect(Site_Name, paste(c("Manistique",
                                                         "Browns",
                                                         "Elk",
                                                         "Glen",
                                                         "Green",
                                                         "Harwood",
                                                         "Howard",
                                                         "Lime",
                                                         "Murray",
                                                         "Leelanau",
                                                         "North Sand",
                                                         "Walloon"), collapse = '|')))

cleaned <- bind_rows(waters, sites) %>% distinct(Stocking_Date, Waters_Name, .keep_all = TRUE)
write.csv(cleaned, "./LakeStockingRecords.csv")

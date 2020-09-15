## code to prepare `DATASET` dataset goes here
library(tidyverse)

PISA15MATH <- read_csv("PISA_2015_M1_M2_complete.csv") %>% 
        select(starts_with("CM") & ends_with("T"))
usethis::use_data(PISA15MATH, overwrite = TRUE, internal = T)



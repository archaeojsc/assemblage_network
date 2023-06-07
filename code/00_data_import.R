require(tidyverse)

# Import data from file ---------------------------------------------------

# Data for Site A
dat_raw <- read_csv("./data/Catalog_SiteA.csv",
                    col_select = c(LEVEL_ID, CODE))

# Data for Sites A and C
# dat_raw <- read_csv("./data/Catalog_AC.csv",
#                     col_select = c(LEVEL_ID, CODE))

# Artifacts to exclude from analysis
exclude_artifacts <-
  as.vector(read.csv("./data/code_exclude.csv", header = T)$x)

# Remove Bakelite and plastic buttons from exclusion list
exclude_artifacts <-
  exclude_artifacts[!exclude_artifacts %in% c("BKLT", "PB")]

# Import list of artifact codes
artifact_codes <-
  read.csv("./data/code_list.csv",
           header = TRUE,
           stringsAsFactors = FALSE)

# Filter data set for excluded artifact types
dat <- dat_raw %>% filter(!(CODE %in% exclude_artifacts))
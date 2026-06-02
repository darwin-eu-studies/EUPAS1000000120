# Copyright 2026 European Medicines Agency
#
# This software is developed by the DARWIN EU Coordination Centre
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# This study code contains unmodified open-source dependencies.

# # load packages -----
library(shiny)
library(shinydashboard)
library(dplyr)
library(shinyWidgets)
library(pool)

# options(shiny.maxRequestSize = 100000*1024^2)

source("functions.R")

outputPath <- here::here("results")

dbdir <- file.path(outputPath, "results_combined.duckdb")
con <- pool::dbPool(duckdb::duckdb(), dbdir = dbdir)
# con <- DBI::dbConnect(duckdb::duckdb(), dbdir = dbdir)

DBI::dbListTables(con)

cohorts <- tbl(con, "cohort_counts") %>% collect
demo <- tbl(con, "demographics")  #%>% collect
prev <- tbl(con, "prevalence") #%>% collect
hosp <- tbl(con, "hospitalization_death") %>% collect
cohort_attrition <- tbl(con, "cohort_attrition") %>% collect
drugs <- tbl(con, "ls_drugs") #%>% collect 
conditions <- tbl(con, "ls_conditions") #%>% collect 

demo %>% 
  filter(strata_name == "frailty_category") %>% 
  count(strata_level)

demo %>% 
  count(frailty_category)

prev %>% 
  count(cdm_name)

# cdm snapshot ------
cdm_snapshot_files <- list.files(outputPath, pattern = ".csv", full.names = T) %>% 
  stringr::str_subset("snapshot")

cdm_snapshot <- list()
for (i in seq_along(cdm_snapshot_files)) {
  cdm_snapshot[[i]] <- readr::read_csv(cdm_snapshot_files[[i]], show_col_types = FALSE) %>% 
    select("cdm_name", "person_count", "observation_period_count", "vocabulary_version")
}


# age_group <chr>, sex <chr>, polypharm_gte_5 <chr>, polypharm_gte_10 <chr>, 
# frailty_category <chr>, pre_2020 

cdm_snapshot <- dplyr::bind_rows(cdm_snapshot) %>% 
  rename("Database name" = "cdm_name",
         "Persons in the database" = "person_count",
         "Number of observation periods" = "observation_period_count",
         "OMOP CDM vocabulary version" = "vocabulary_version")





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

# this script extracts the cohort json files from Atlas.
# The cancer cohorts were created by Julieta and reviewed by Adam on March 4, 2024
# reloaded concept sets on May 14

library(ROhdsiWebApi)
library(dplyr)

baseUrl <- "https://atlas-dev.darwin-eu.org/WebAPI"

# remove the inst folder if it exists and recreate it
dir.create(here::here("inst"))

# this token was copied from the browser 
token <- ""
setAuthHeader(baseUrl, authHeader = token)

md0 <- getCohortDefinitionsMetaData(baseUrl)

md <- md0 %>% 
  filter(stringr::str_detect(name, "C1-009")) %>% 
  select(name, id) %>% 
  mutate(name2 = snakecase::to_snake_case(tolower(stringr::str_remove(name, "C1-009 ")))) %>% 
  mutate(name2 = stringr::str_remove(name2, "_jp$"))

cohortsToCreate <- dplyr::tibble(atlasId = md$id, cohortId = md$id, cohortName = md$name2)
readr::write_csv(cohortsToCreate, here::here("inst", "cohortsToCreate.csv"))

insertCohortDefinitionSetInPackage(
  fileName = here::here("inst/cohortsToCreate.csv"),
  baseUrl = "https://atlas-dev.darwin-eu.org/WebAPI",
  packageName = "FrailtyPolypharmacyCancerDiagnostics",
  insertCohortCreationR = F
)

for (i in seq_len(nrow(cohortsToCreate))) {
  file.rename(here::here("inst/cohorts", paste0(cohortsToCreate[i, "cohortId"], ".json")),
              here::here("inst/cohorts", paste0(cohortsToCreate[i, "cohortName"], ".json")))
}

unlink(here::here("inst", "sql"), recursive = T)
file.remove(here::here("inst", "cohortsToCreate.csv"))

df0 <- getConceptSetDefinitionsMetaData(baseUrl)

df <- df0 %>% 
  filter(stringr::str_detect(name, "C1-009")) %>% 
  filter(stringr::str_detect(name, "frailty|Frailty")) %>% 
  select(id, name) %>% 
  mutate(name = stringr::str_remove(name, "C1-009 ")) %>% 
  mutate(name = stringr::str_remove(name, "_JP$")) %>% 
  mutate(name = snakecase::to_snake_case(name)) %>% 
  mutate(name = tolower(name)) %>% 
  filter(name != "frailty")

print(df, n=100)

dir.create(here::here("inst/frailty_concept_sets"))

for (i in cli::cli_progress_along(1:nrow(df))) {
  cpt <- getConceptSetDefinition(df$id[i], baseUrl)
  
  jsonlite::write_json(cpt$expression, 
                       path = here::here("inst", "frailty_concept_sets", paste0(df$name[i], ".json")),
                       pretty = T,
                       auto_unbox = T)
}



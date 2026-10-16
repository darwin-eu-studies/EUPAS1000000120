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

# 
# # returns an example cdm to use for tests
# mock_cdm <- function(prefix = "c1009_") {
#   
#   checkmate::assertCharacter(prefix, any.missing = F, min.chars = 1)
#   
#   require(CDMConnector)
#   
#   con <- DBI::dbConnect(duckdb::duckdb(eunomia_dir("synthea-lung_cancer-10k")))
#   cdm <- cdm_from_con(con, 
#                       cdm_schema = "main", 
#                       write_schema = c(schema = "main", prefix = prefix))
#   
#   concept_set_df <- readr::read_csv(here::here("data", "concept_sets.csv"), col_types = "cci")
#   
#   # insert records so we have some test data
#   co <- cdm$condition_occurrence %>% 
#     dplyr::collect()
#   
#   max_id <- max(co$condition_occurrence_id)
#   nrows_to_add <- 100000
#   set.seed(1)
#   df <- dplyr::tibble(
#     condition_occurrence_id = (max_id+1):(max_id+nrows_to_add),
#     person_id = sample(co$person_id, nrows_to_add, replace = T),
#     condition_concept_id = sample(concept_set_df$concept_id, nrows_to_add, replace = T),
#     condition_start_date = sample(co$condition_start_date, nrows_to_add, replace = T),
#     condition_end_date = condition_start_date + 1)
#   
#   DBI::dbAppendTable(con, "condition_occurrence", df)
#   
#   cdm
# }

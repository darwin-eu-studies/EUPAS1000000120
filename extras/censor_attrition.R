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

library(dplyr)

# in iqvia we needed to censor attrition

cohort_attrition <- tibble::tribble(
           ~cdm_name, ~cohort_definition_id, ~cohort_name, ~number_records, ~number_subjects, ~reason_id,                                    ~reason, ~excluded_records, ~excluded_subjects,
              "cdm",                    1L,       "chrt1",            249L,             249L,         1L,               "Qualifying initial records",                0L,                 0L,
              "cdm",                    1L,       "chrt1",            196L,             196L,         2L,                          "no prior cancer",               53L,                53L,
              "cdm",                    1L,       "chrt1",            193L,             193L,         3L,                                "Age >= 18",                3L,                 3L,
              "cdm",                    8L,       "chrt2",           1278L,            1278L,         1L,               "Qualifying initial records",                0L,                 0L,
              "cdm",                    8L,       "chrt2",           1079L,            1079L,         2L,                          "no prior cancer",              199L,               199L,
              "cdm",                    8L,       "chrt2",           1029L,            1029L,         3L,                                "Age >= 18",               50L,                50L,
  )

supress_cohort <- function(cohort_attrition, min_cell_count = 5) {
  cohort_attrition %>% 
    group_by(cohort_name, cdm_name) %>% 
    mutate(
      max_reason_id = max(reason_id),
      censor_group = ifelse(
        dplyr::between(excluded_subjects, 1, min_cell_count-1) | 
          dplyr::between(excluded_records, 1, min_cell_count-1) | 
          dplyr::between(number_records, 1, min_cell_count-1) |
          dplyr::between(number_subjects, 1, min_cell_count-1), T, F)) %>% 
    mutate(censor_group = max(censor_group, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(
      number_records    = ifelse(censor_group == 1 & reason_id < max_reason_id, NA_real_, number_records), 
      number_subjects   = ifelse(censor_group == 1  & reason_id < max_reason_id, NA_real_, number_subjects), 
      excluded_records  = ifelse(censor_group == 1 , NA_real_, excluded_records), 
      excluded_subjects = ifelse(censor_group == 1 , NA_real_, excluded_subjects)) %>% 
    select(-max_reason_id, -censor_group)
}

# before censoring
cohort_attrition %>% 
  select(-cohort_definition_id) %>% 
  rename(n_rec = number_records, n_sub = number_subjects, 
         ex_rec = excluded_records, ex_subj = excluded_subjects) %>% 
  print(n=1e6)

# after censoring
supress_cohort(cohort_attrition, 5) %>% 
  select(-cohort_definition_id) %>% 
  rename(n_rec = number_records, n_sub = number_subjects, 
         ex_rec = excluded_records, ex_subj = excluded_subjects) %>% 
  print(n=1e6)

  
  

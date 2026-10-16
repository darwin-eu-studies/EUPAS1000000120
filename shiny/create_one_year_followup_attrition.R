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



df <- demo %>% 
  select(cdm_name, cohort = group_level, strata_level, variable_name, estimate_value) %>% 
  filter(strata_level == "overall", variable_name == "Number subjects") %>% 
  transmute(cdm_name, cohort, number_subjects_before = as.numeric(estimate_value))
  

df2 <- hosp %>% 
  select(cdm_name, cohort = group_level, strata_level, number_subjects_count) %>% 
  filter(strata_level == "overall") %>% 
  transmute(cdm_name, cohort, number_subjects_after = as.numeric(number_subjects_count))

df3 <- left_join(df, df2, by = c("cdm_name", "cohort")) %>% 
  mutate(number_subjects_dropped = number_subjects_before - number_subjects_after)


readr::write_csv(df3, "one_year_followup_attrition.csv")



# Copyright 2025 European Medicines Agency
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


purrr::walk(list.files(here::here("R"), full.names = T), source)

prefix <- "c1009_"

cdm <- mock_cdm(prefix)


# generate the cohorts

library(CDMConnector)
library(here)
library(dplyr)

cohort_set <- read_cohort_set(here::here("inst", "cohorts"))
cohorts <- readr::read_csv(here("inst", "cohortsToCreate.csv"), col_types = "iic")

con <- attr(cdm, "dbcon")
write_schema <- attr(cdm, "write_schema")

tables_to_drop <- list_tables(con, schema = write_schema) %>% 
  stringr::str_subset("cancer_cohort")

purrr::walk(tables_to_drop, ~dropTable(cdm, name = .))

# based on behavior of CDMConnector 1.3 
cohort_set <- cohort_set %>% 
  mutate(cohort_name = ifelse(cohort_name == "45", "450", cohort_name)) %>% 
  mutate(cohort_definition_id = as.integer(cohort_name)) %>% 
  left_join(cohorts, by = c("cohort_definition_id" = "cohortId")) %>% 
  select(cohort_definition_id,
         cohort_name = cohortName,
         cohort_name_snakecase = cohortName,
         cohort, 
         json) %>% 
  filter(cohort_name != "frailty")

cdm <- generate_cohort_set(cdm, 
                           cohort_set, 
                           name = "cancer_cohort",
                           compute_attrition = TRUE,
                           overwrite = TRUE)

# confirm that each person is in the cohort at most one time

cdm$cancer_cohort %>% 
  count(subject_id, name = "n_records") %>% 
  count(n_records, name = "n_persons")


# compute polypharmacy

# Drug era 

tmp1 <- cdm$cancer_cohort %>% 
  mutate(window_start = !!dateadd("cohort_start_date", -90), window_end = cohort_start_date) %>% 
  select(cohort_definition_id, subject_id, window_start, window_end) %>% 
  inner_join(cdm$drug_era, by = c("subject_id" = "person_id")) %>% 
  filter(window_start <= drug_era_end_date, drug_era_start_date <= window_end) %>% 
  compute(temporary = FALSE, overwrite = TRUE)
  
n_drugs <- union_all(
    transmute(tmp1, cohort_definition_id, subject_id, date = drug_era_start_date, flag = 1L),
    transmute(tmp1, cohort_definition_id, subject_id, date = drug_era_end_date, flag = -1L)
  ) %>% 
  group_by(cohort_definition_id, subject_id) %>% 
  arrange(date) %>% 
  # dbplyr::window_order(date) %>% 
  mutate(cumsum_flag = cumsum(flag)) %>% 
  # dplyr::show_query() # check that cumulative sum sql is correct
  group_by(cohort_definition_id, subject_id) %>% 
  summarise(n_drug_ingredients = max(cumsum_flag, na.rm = T)) %>% 
  compute(name = "polypharamacy", temporary = FALSE, overwrite = TRUE)
  # compute_query(name = "polypharamacy", temporary = FALSE, overwrite = TRUE, schema = write_schema)


# check that this includes descendants
concepts <- CodelistGenerator::codesFromConceptSet(cdm, path = here("inst", "frailtyConceptSets"))

names(concepts) <- snakecase::to_snake_case(names(concepts))

n_concept_sets <- length(concepts)

# compute frailty
frailty <- cdm$cancer_cohort %>% 
  PatientProfiles::addSex() %>% 
  PatientProfiles::addAge(ageGroup = list(c(18,44), c(45,64), c(65,74), c(75,84), c(85,120))) %>% 
  PatientProfiles::addConceptIntersectFlag(conceptSet = concepts, 
                                           window = c(-99999, 0),
                                           nameStyle = "{concept_name}") 


analytic_dataset <- frailty %>% 
  mutate(frailty_score = 
    frailty_hearing_impairment +
    frailty_mobility_transfer +
    frailty_fragility_fracture +
    frailty_anemia +
    frailty_c_kidney_dx +
    frailty_osteoporosis +
    frailty_visual_impairment +
    frailty_parkinsonism_tremor +
    frailty_cerebrovascular_dx +
    frailty_activity_limitation +
    frailty_a_fib +
    frailty_sleep_disturbance +
    frailty_social_vulnerability +
    frailty_diabetes +
    frailty_foot_problem +
    frailty_peripheral_vasc_dx +
    frailty_hypertension +
    frailty_heart_valv_dx +
    frailty_dyspnea +
    frailty_falls +
    frailty_hypotension +
    frailty_respiratory_dx +
    frailty_peptic_ulc +
    frailty_thyroid_dx +
    frailty_memeroy_cogn_problem +
    frailty_heart_failure +
    frailty_req_care +
    frailty_urinary_incontinence +
    frailty_weight_loss_anorexia_jp +
    frailty_ischemic_heart_dx +
    frailty_artrtitis +
    frailty_housebound +
    frailty_skin_ulcer
  ) %>% 
  left_join(n_drugs, by = c("cohort_definition_id", "subject_id")) %>% 
  mutate(n_drug_ingredients = coalesce(n_drug_ingredients, 0L)) %>% 
  mutate(
    polypharm_gte5 = ifelse(n_drug_ingredients >= 5L, 1L, 0L),
    polypharm_gte10 = ifelse(n_drug_ingredients >= 10L, 1L, 0L)
  ) %>% 
  compute(name = "analytic_dataset", overwrite = TRUE, temporary = FALSE)

# analytic_dataset$n_drug_ingredients

# compute counts
counts_by_cohort <- analytic_dataset %>% 
  group_by(cohort_definition_id) %>% 
  summarise(n_person = n_distinct(subject_id), n_records = n()) %>% 
  collect() %>% 
  mutate(
    n_person = ifelse(n_person < 5 & n_person > 0, -1, n_person),
    n_records = ifelse(n_records < 5 & n_records > 0, -1, n_records),
  )

counts_by_age <- analytic_dataset %>% 
  group_by(cohort_definition_id, age_group) %>% 
  summarise(n_person = n_distinct(subject_id), n_records = n()) %>% 
  collect() %>% 
  mutate(
    n_person = ifelse(n_person < 5 & n_person > 0, -1, n_person),
    n_records = ifelse(n_records < 5 & n_records > 0, -1, n_records),
  )

counts_by_sex <- analytic_dataset %>% 
  group_by(cohort_definition_id, sex) %>% 
  summarise(n_person = n_distinct(subject_id), n_records = n()) %>% 
  collect() %>% 
  mutate(
    n_person = ifelse(n_person < 5 & n_person > 0, -1, n_person),
    n_records = ifelse(n_records < 5 & n_records > 0, -1, n_records),
  )

counts_by_age_sex <- analytic_dataset %>% 
  group_by(cohort_definition_id, age_group, sex) %>% 
  summarise(n_person = n_distinct(subject_id), n_records = n()) %>% 
  collect() %>% 
  mutate(
    n_person = ifelse(n_person < 5 & n_person > 0, -1, n_person),
    n_records = ifelse(n_records < 5 & n_records > 0, -1, n_records),
  )
  
fs::dir_create(here::here("output"))
readr::write_csv(counts_by_cohort, here::here("output", "counts_by_cohort.csv"))
readr::write_csv(counts_by_age, here::here("output", "counts_by_age.csv"))
readr::write_csv(counts_by_sex, here::here("output", "counts_by_sex.csv"))
readr::write_csv(counts_by_age_sex, here::here("output", "counts_by_age_sex.csv"))
readr::write_csv(settings(cdm$cancer_cohort), here::here("output", "cohort_names.csv"))
cli::cat_line("saved cohort counts")

# get frailty score quantiles overall and by groups

probs <- seq(0, 1, by = .1)

frailty_by_cohort <- analytic_dataset %>% 
  group_by(cohort_definition_id) %>% 
  summarise_quantile(frailty_score, probs = probs) %>% 
  collect() %>% 
  tidyr::gather("quantile", "value", matches("^p")) %>% 
  mutate(measure = "frailty") 

readr::write_csv(frailty_by_cohort, here::here("output", "frailty_by_cohort.csv"))
cli::cat_line("saved frailty_by_cohort.csv")


frailty_by_age <- analytic_dataset %>% 
  group_by(cohort_definition_id, age_group) %>% 
  summarise_quantile(frailty_score, probs = probs) %>% 
  collect() %>% 
  tidyr::gather("quantile", "value", matches("^p")) %>% 
  mutate(measure = "frailty") 

readr::write_csv(frailty_by_age, here::here("output", "frailty_by_age.csv"))
cli::cat_line("saved frailty_by_age.csv")

frailty_by_sex <- analytic_dataset %>% 
  group_by(cohort_definition_id, sex) %>% 
  summarise_quantile(frailty_score, probs = probs) %>% 
  collect() %>% 
  tidyr::gather("quantile", "value", matches("^p")) %>% 
  mutate(measure = "frailty") 

readr::write_csv(frailty_by_sex, here::here("output", "frailty_by_sex.csv"))
cli::cat_line("saved frailty_by_sex.csv")

frailty_by_age_sex <- analytic_dataset %>% 
  group_by(cohort_definition_id, age_group, sex) %>% 
  summarise_quantile(frailty_score, probs = probs) %>% 
  collect() %>% 
  tidyr::gather("quantile", "value", matches("^p")) %>% 
  mutate(measure = "frailty") 

readr::write_csv(frailty_by_sex, here::here("output", "frailty_by_age_sex.csv"))
cli::cat_line("saved frailty_by_age_sex.csv")

# polypharmacy 

polypharm <- analytic_dataset %>% 
  select(cohort_definition_id, 
         subject_id, 
         cohort_start_date, 
         age_group, 
         sex, 
         polypharm_gte5,
         polypharm_gte10)


polypharm_by_age_sex <- polypharm %>% 
  group_by(cohort_definition_id, age_group, sex) %>% 
  summarise(
    denom = n(),
    polypharm_gte5 = sum(polypharm_gte5, na.rm = T),
    polypharm_gte10 = sum(polypharm_gte10, na.rm = T)
  ) %>% 
  collect() %>% 
  tidyr::gather("measure", "value", matches("^poly")) %>% 
  mutate(denom = ifelse(denom < 5 & denom > 0, -1, denom),
         value = ifelse(value < 5 & value > 0, -1, value)) %>% 
  mutate(rate = ifelse(denom >= 0 & value >= 0, value / denom, -1))

readr::write_csv(polypharm_by_age_sex, here::here("output", "polypharm_by_age_sex.csv"))
cli::cat_line("saved polypharm_by_age_sex.csv")

polypharm_by_age <- polypharm %>% 
  group_by(cohort_definition_id, age_group) %>% 
  summarise(
    denom = n(),
    polypharm_gte5 = sum(polypharm_gte5, na.rm = T),
    polypharm_gte10 = sum(polypharm_gte10, na.rm = T)
  ) %>% 
  collect() %>% 
  tidyr::gather("measure", "value", matches("^poly")) %>% 
  mutate(denom = ifelse(denom < 5 & denom > 0, -1, denom),
         value = ifelse(value < 5 & value > 0, -1, value)) %>% 
  mutate(rate = ifelse(denom >= 0 & value >= 0, value / denom, -1))

readr::write_csv(polypharm_by_age, here::here("output", "polypharm_by_age.csv"))
cli::cat_line("saved polypharm_by_age.csv")

polypharm_by_sex <- polypharm %>% 
  group_by(cohort_definition_id, sex) %>% 
  summarise(
    denom = n(),
    polypharm_gte5 = sum(polypharm_gte5, na.rm = T),
    polypharm_gte10 = sum(polypharm_gte10, na.rm = T)
  ) %>% 
  collect() %>% 
  tidyr::gather("measure", "value", matches("^poly")) %>% 
  mutate(denom = ifelse(denom < 5 & denom > 0, -1, denom),
         value = ifelse(value < 5 & value > 0, -1, value)) %>% 
  mutate(rate = ifelse(denom >= 0 & value >= 0, value / denom, -1))

readr::write_csv(polypharm_by_sex, here::here("output", "polypharm_by_sex.csv"))
cli::cat_line("saved polypharm_by_sex.csv")

polypharm_by_cohort <- polypharm %>% 
  group_by(cohort_definition_id) %>% 
  summarise(
    denom = n(),
    polypharm_gte5 = sum(polypharm_gte5, na.rm = T),
    polypharm_gte10 = sum(polypharm_gte10, na.rm = T)
  ) %>% 
  collect() %>% 
  tidyr::gather("measure", "value", matches("^poly")) %>% 
  mutate(denom = ifelse(denom < 5 & denom > 0, -1, denom),
         value = ifelse(value < 5 & value > 0, -1, value)) %>% 
  mutate(rate = ifelse(denom >= 0 & value >= 0, value / denom, -1))

readr::write_csv(polypharm_by_cohort, here::here("output", "polypharm_by_cohort.csv"))
cli::cat_line("saved polypharm_by_cohort.csv")


analytic_dataset

cdm <- insertTable(cdm, "analytic_dataset", analytic_dataset, overwrite = T)


cdm$cancer_cohort <- cdm$cancer_cohort %>% 
  PatientProfiles::addSex() %>% 
  PatientProfiles::addAge(ageGroup = list(c(18,44), c(45,64), c(65,74), c(75,84), c(85,120))) %>% 
  dplyr::compute(name = "cancer_cohort", temporary = FALSE, overwrite = TRUE) %>% 
  newCohortTable()


settings(cdm$cancer_cohort)

# large scale --------
cli::cli_alert_info("Running large scale characterisation")

ls_condition_occurrence <- cdm$cancer_cohort %>%
  PatientProfiles::summariseLargeScaleCharacteristics(
    strata =  list(c("age_group"),
                   c("sex")),
    window = list(c(-999999, -366),
                  c(-365, -31),
                  c(-30, -1),
                  c(0, 0)),
    episodeInWindow = "condition_occurrence")

ls_drug_exposure <- cdm$cancer_cohort %>%
  PatientProfiles::summariseLargeScaleCharacteristics(
    strata =  list(c("age_group"),
                   c("sex")),
    window = list(c(-999999, -366),
                  c(-365, -31),
                  c(-30, -1),
                  c(0, 0)),
    episodeInWindow = "drug_exposure")

ls_characteristics <- ls_condition_occurrence %>%
  bind_rows(ls_drug_exposure)

readr::write_csv(ls_characteristics, 
                 here("export", paste0("ls_characteristics_", cdmName(cdm), ".csv")))



cdm_disconnect(cdm)







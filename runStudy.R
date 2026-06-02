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


# generate the study cancer cohorts
# add all needed extra covariates

library(CDMConnector)
library(here)
library(dplyr)

# save results to a duckdb database
res <- DBI::dbConnect(duckdb::duckdb(), file.path(output_folder, paste0("results_", snakecase::to_snake_case(cdmName(cdm)), ".duckdb")))

DBI::dbListTables(res)

cohort_set <- read_cohort_set(here::here("inst", "cohorts")) %>% 
  filter(cohort_name != "frailty")

log4r::info(logger = logger, "Generating cohorts")
cdm <- generate_cohort_set(cdm,
                           cohort_set,
                           name = "cancer_cohort0",
                           compute_attrition = TRUE,
                           overwrite = TRUE)

log4r::info(logger = logger, "saving cohort counts")

cohort_count(cdm$cancer_cohort0) %>% 
  inner_join(settings(cdm$cancer_cohort0), by = "cohort_definition_id") %>% 
  mutate(cdm_name = cdmName(cdm)) %>% 
  select(cdm_name, cohort_definition_id, cohort_name, everything()) %>% 
  {DBI::dbWriteTable(res, name = "cohort_counts", value = ., overwrite = T)}

attrition(cdm$cancer_cohort0) %>% 
  inner_join(settings(cdm$cancer_cohort0), by = "cohort_definition_id") %>% 
  mutate(cdm_name = cdmName(cdm)) %>% 
  select(cdm_name, cohort_definition_id, cohort_name, everything()) %>% 
  as.data.frame() %>% 
  {DBI::dbWriteTable(res, name = "cohort_attrition", value = ., overwrite = T)}

if (all(cohort_count(cdm$cancer_cohort0)$number_subjects == 0)) {
  cli::cli_abort("None of the study cohorts matched any persons in your cdm database. You cannot run the study.")
}

# confirm that each person is in the cohort at most one time
check <- cohort_count(cdm$cancer_cohort0) %>%
  mutate(number_records = as.numeric(number_records), number_subjects = as.numeric(number_subjects)) %>% 
  mutate(check = (number_records == number_subjects))

if (!all(check$check)) {
  print(check)
  cli::cli_abort("Some persons are in the cohort table more than once. We only expect a single record per person.")
}


# compute polypharmacy ----

# Get all drug era records that intersect with 90 days prior to cancer index
tmp1 <- cdm$cancer_cohort0 %>%
  mutate(window_start = !!dateadd("cohort_start_date", -90), window_end = cohort_start_date) %>%
  select(cohort_definition_id, subject_id, window_start, window_end) %>%
  inner_join(cdm$drug_era, by = c("subject_id" = "person_id")) %>%
  filter(window_start <= drug_era_end_date, drug_era_start_date <= window_end) %>%
  compute(temporary = FALSE, overwrite = TRUE, name = "tmp1")

# Stack the drug era records in a long format and use a +1, -1 encoding.
# Then find the max number of drugs someone is on at any one time
n_drugs <- union_all(
    transmute(tmp1, cohort_definition_id, subject_id, date = drug_era_start_date, flag = 1L),
    transmute(tmp1, cohort_definition_id, subject_id, date = drug_era_end_date, flag = -1L)
  ) %>%
  group_by(cohort_definition_id, subject_id) %>%
  arrange(date) %>%
  mutate(cumsum_flag = cumsum(flag)) %>%
  group_by(cohort_definition_id, subject_id) %>%
  summarise(n_drug_ingredients = max(cumsum_flag, na.rm = T)) %>%
  mutate(
    polypharm_gte5 = ifelse(n_drug_ingredients >= 5L, 1L, 0L),
    polypharm_gte10 = ifelse(n_drug_ingredients >= 10L, 1L, 0L)
  ) %>%
  compute(name = "polypharamacy", temporary = FALSE, overwrite = TRUE)


# compute frailty ----
concepts <- CodelistGenerator::codesFromConceptSet(cdm, path = here("inst", "frailty_concept_sets"))
names(concepts) <- snakecase::to_snake_case(names(concepts))

# Remove concepts in the measurement value domain
concepts <- purrr::map(concepts, ~.[!(. %in% c(1621081, 45878235, 45877743, 45879223, 45878557))])


if (length(concepts) != 35) {
  stop('There should be 35 frailty concept sets. Something is wrong.')
}

log4r::info(logger, "Add frailty concepts to cohorts as flags")

# suppress warning about long sql query
suppressWarnings({
  frailty <- cdm$cancer_cohort0 %>%
    PatientProfiles::addDemographics(
      ageGroup = list(c(18,44), c(45,64), c(65,74), c(75,84), c(85,120)),
      priorObservation = FALSE,
      futureObservation = FALSE
    ) %>% 
    PatientProfiles::addConceptIntersectFlag(conceptSet = concepts,
                                             window = c(-99999, 0),
                                             nameStyle = "{concept_name}")
})

# denominator for frailty score. Add 1 for polypharmacy (>=5)
score_denominator <- length(concepts) + 1 

log4r::info(logger, "add death date")
# get first death date for each person after index
death_date <- cdm$cancer_cohort0 %>%
  inner_join(select(cdm$death, "person_id", "death_date"), by = c("subject_id" = "person_id")) %>%
  filter(death_date >= cohort_start_date) %>%
  group_by(cohort_definition_id, subject_id) %>%
  summarize(death_date = min(death_date, na.rm = TRUE)) %>%
  compute(temporary = FALSE, name = "death_date", overwrite = TRUE)

log4r::info(logger, "Get hospitalizations")
n_hospitalizations <- cdm$cancer_cohort0 %>%
  inner_join(select(cdm$visit_occurrence, "person_id", "visit_start_date", "visit_concept_id"), by = c("subject_id" = "person_id")) %>%
  # concept id 9201 = inpatient visit; concept id 262 = emergency room + inpatient visit
  filter(visit_concept_id %in% c(9201, 262),
         visit_start_date >= cohort_start_date,
         visit_start_date <= !!dateadd("cohort_start_date", 356)) %>%
  count(cohort_definition_id, subject_id, name = "n_hospitalizations_in_year_after_dx") %>%
  compute(temporary = FALSE, name = "n_hospitalizations", overwrite = TRUE)


log4r::info(logger, "Creating final analytic cohort table")

cancer_cohort <- frailty %>%
  left_join(n_drugs, by = c("cohort_definition_id", "subject_id")) %>%
  mutate(n_drug_ingredients = coalesce(n_drug_ingredients, 0L),
         polypharm_gte5 = coalesce(polypharm_gte5, 0L),
         polypharm_gte10 = coalesce(polypharm_gte10, 0L),
         index_year = !!datepart("cohort_start_date", "year")) %>%
  mutate(
    pre_2020 = ifelse(index_year < 2020L, TRUE, FALSE),
    score = (frailty_activity_limitation + frailty_anemia + frailty_arthritis + 
       frailty_atrial_fibrillation + frailty_c_kidney_dx + frailty_cerebrovascular_dx + 
       frailty_diabetes + frailty_dyspnea + frailty_falls + frailty_foot_problem + 
       frailty_fragility_fracture + frailty_hearing_impairment + frailty_heart_failure + 
       frailty_heart_valve_dx + frailty_housebound + frailty_hypertension + 
       frailty_hypotension + frailty_ischemic_heart_dx + frailty_memeroy_cogn_problem + 
       frailty_mobility_transfer + frailty_osteoporosis + frailty_parkinsonism_tremor + 
       frailty_peptic_ulcer + frailty_peripheral_vascular_dx + frailty_req_care + 
       frailty_respiratory_dx + frailty_skin_ulcer + frailty_sleep_disturbance + 
       frailty_social_vulnerability + frailty_thyroid_dx + frailty_urinary_incontinence + 
       frailty_visual_impairment + frailty_weight_loss_anorexia_jp + frailty_dizziness + 
       frailty_urinary_system_disease + polypharm_gte5) / (1.0* local(score_denominator))) %>%
  # select(-matches("frailty_")) %>%
  left_join(death_date, by = c("cohort_definition_id", "subject_id")) %>%
  left_join(n_hospitalizations, by = c("cohort_definition_id", "subject_id")) %>%
  mutate(n_hospitalizations_in_year_after_dx = coalesce(n_hospitalizations_in_year_after_dx, 0L),
         cutoff_date = !!as.Date("2021-12-31")) %>%
  mutate(frailty_category = case_when(
    0 <= score & score < .12 ~ "fit",
    .12 <= score & score < .24 ~ "mild",
    .24 <= score & score < .36 ~ "moderate",
    .36 <= score ~ "severe",
    TRUE ~ "NA"),
    days_to_cutoff = !!datediff("cohort_start_date", "cutoff_date")) %>%
  mutate(has_possible_followup_flag = ifelse(days_to_cutoff > 0, 1L, 0L),
         death_flag = ifelse(is.na(death_date), 0L, 1L),
         hospitalization_flag = ifelse(n_hospitalizations_in_year_after_dx > 0, 1L, 0L)) %>%
  PatientProfiles::addCohortName() %>% 
  select("cohort_definition_id", "cohort_name", "subject_id", "cohort_start_date",
         "cohort_end_date", 
         "sex", "age", "age_group",
         "has_possible_followup_flag",
         "n_drug_ingredients",
         polypharm_gte_5 = "polypharm_gte5",
         polypharm_gte_10 = "polypharm_gte10",
         "score", "frailty_category",
         "n_hospitalizations_in_year_after_dx",
         "death_flag", "hospitalization_flag", "pre_2020", starts_with("frailty")) %>%
  compute(name = "cancer_cohort", overwrite = TRUE, temporary = FALSE)

cdm$cancer_cohort <- newCohortTable(cancer_cohort)

log4r::info(logger, "Cohort table created")

strata <- c("age_group", 
            "sex",
            "polypharm_gte_5", 
            "polypharm_gte_10", 
            "frailty_category",
            "pre_2020")

# need to include all intersections of all strata
create_all_strata <- function(strata) {
  
  result <- as.list(strata)
  
  if (length(strata) < 2) return(result)
  
  for (m in 2:length(strata)) {
    x <- combn(strata, m = m)
    
    for (i in 1:ncol(x)) {
      result[[length(result) + 1]] <- x[, i]
    }
  }
  result
}

all_strata_combinations <- create_all_strata(strata)

# we don't need to include both polypharm flags in the same strata
all_strata_combinations <- purrr::discard(all_strata_combinations, ~("polypharm_gte_5" %in% . && "polypharm_gte_10" %in% .))


# demographics ----
log4r::info(logger, "Summarising demographics")
demographics <- cdm$cancer_cohort %>%
  PatientProfiles::summariseCharacteristics(strata = all_strata_combinations)

demographics %>% 
  PatientProfiles::suppress(minCellCount) %>% 
  as.data.frame() %>% 
  mutate(strata_nm = stringr::str_split(strata_name, " and "),
         strata_lvl = stringr::str_split(strata_level, " and "),
         strata = purrr::map2(strata_lvl, strata_nm, setNames)) %>% 
  mutate(age_group = purrr::map_chr(strata, ~.["age_group"]),
         sex = purrr::map_chr(strata, ~.["sex"]),
         polypharm_gte_5 = purrr::map_chr(strata, ~.["polypharm_gte_5"]),
         polypharm_gte_10 = purrr::map_chr(strata, ~.["polypharm_gte_10"]),
         frailty_category = purrr::map_chr(strata, ~.["frailty_category"]),
         pre_2020 = purrr::map_chr(strata, ~.["pre_2020"])) %>% 
  mutate(across(c("sex", "age_group", "polypharm_gte_5", "polypharm_gte_10", "frailty_category", "pre_2020"), ~tidyr::replace_na(., "NA"))) %>% 
  select(-strata_nm, -strata_lvl, -strata) %>% 
  {DBI::dbWriteTable(res, name = "demographics", value = ., overwrite = T)}

rm(demographics)

# Prevalence ----
log4r::info(logger, "Computing prevalence")

binaryVariables <- c("frailty_activity_limitation","frailty_anemia","frailty_arthritis",
       "frailty_atrial_fibrillation","frailty_c_kidney_dx","frailty_cerebrovascular_dx",
       "frailty_diabetes","frailty_dyspnea","frailty_falls","frailty_foot_problem",
       "frailty_fragility_fracture","frailty_hearing_impairment","frailty_heart_failure",
       "frailty_heart_valve_dx","frailty_housebound","frailty_hypertension",
       "frailty_hypotension","frailty_ischemic_heart_dx","frailty_memeroy_cogn_problem",
       "frailty_mobility_transfer","frailty_osteoporosis","frailty_parkinsonism_tremor",
       "frailty_peptic_ulcer","frailty_peripheral_vascular_dx","frailty_req_care",
       "frailty_respiratory_dx","frailty_skin_ulcer","frailty_sleep_disturbance",
       "frailty_social_vulnerability","frailty_thyroid_dx","frailty_urinary_incontinence",
       "frailty_visual_impairment","frailty_weight_loss_anorexia_jp", "frailty_dizziness", 
       "frailty_urinary_system_disease", "polypharm_gte_5", "polypharm_gte_10")

variables <- list(
  numericVariables = c("age", "n_drug_ingredients", "score"),
  binaryVariables = binaryVariables,
  categoricalVariables = c("sex", "age_group", "frailty_category")
)

# check that variable names are in the cohort table
# unname(unlist(variables))[which(!unname(unlist(variables)) %in% colnames(cdm$cancer_cohort))]
# strata[which(!strata %in% colnames(cdm$cancer_cohort))]

prevalence <- PatientProfiles::summariseResult(
  cdm$cancer_cohort,
  group = list("cohort_name"),
  strata = all_strata_combinations,
  variables = variables,
  functions  = list(numericVariables = c("mean", "sd", "median", "min", "q05", "q25", "q75", "q95", "max"),
                    binaryVariables = c("count", "percentage"),
                    categoricalVariables = c("count", "percentage"))
)

log4r::info(logger, "saving prevalence results")
prevalence %>% 
  PatientProfiles::suppress(minCellCount) %>% 
  as.data.frame() %>% 
  mutate(strata_nm = stringr::str_split(strata_name, " and "),
         strata_lvl = stringr::str_split(strata_level, " and "),
         strata = purrr::map2(strata_lvl, strata_nm, setNames)) %>% 
  mutate(age_group = purrr::map_chr(strata, ~.["age_group"]),
         sex = purrr::map_chr(strata, ~.["sex"]),
         polypharm_gte_5 = purrr::map_chr(strata, ~.["polypharm_gte_5"]),
         polypharm_gte_10 = purrr::map_chr(strata, ~.["polypharm_gte_10"]),
         frailty_category = purrr::map_chr(strata, ~.["frailty_category"]),
         pre_2020 = purrr::map_chr(strata, ~.["pre_2020"])) %>% 
  mutate(across(c("sex", "age_group", "polypharm_gte_5", "polypharm_gte_10", "frailty_category", "pre_2020"), ~tidyr::replace_na(., "NA"))) %>% 
  select(-strata_nm, -strata_lvl, -strata) %>% 
  {DBI::dbWriteTable(res, name = "prevalence", value = ., overwrite = T)}

rm(prevalence)

log4r::info(logger, "Prevalence calculation complete")

# Hospitalization and mortality rates ------
hospitalization_death <- cdm$cancer_cohort %>% 
  filter(has_possible_followup_flag == 1L) %>% 
  PatientProfiles::summariseResult(
    group = list("cohort_name"),
    strata = all_strata_combinations,
    variables = list(
      numericVariables = "n_hospitalizations_in_year_after_dx",
      binaryVariables = "death_flag"
    ),
    functions  = list(numericVariables = "sum", binaryVariables = c("count", "percentage"))
  )

log4r::info(logger, "saving hospitalization and death results")

hospitalization_death %>% 
  PatientProfiles::suppress(minCellCount) %>% 
  as.data.frame() %>% 
  {DBI::dbWriteTable(res, name = "hospitalization_death_raw", value = ., overwrite = T)}
                      
 hospitalization_death %>% 
   PatientProfiles::suppress(minCellCount) %>% 
   mutate(variable_name = stringr::str_replace_all(variable_name, " ", "_")) %>% 
   tidyr::unite(col = "var", variable_name, estimate_name) %>% 
   select(-estimate_type, -variable_level, -additional_level, -additional_name) %>% 
   tidyr::pivot_wider(names_from = "var", values_from = "estimate_value") %>% 
   mutate(hospitalization_rate = as.numeric(n_hospitalizations_in_year_after_dx_sum)/as.numeric(number_subjects_count)) %>% 
   mutate(strata_nm = stringr::str_split(strata_name, " and "),
          strata_lvl = stringr::str_split(strata_level, " and "),
          strata = purrr::map2(strata_lvl, strata_nm, setNames)) %>% 
   mutate(age_group = purrr::map_chr(strata, ~.["age_group"]),
          sex = purrr::map_chr(strata, ~.["sex"]),
          polypharm_gte_5 = purrr::map_chr(strata, ~.["polypharm_gte_5"]),
          polypharm_gte_10 = purrr::map_chr(strata, ~.["polypharm_gte_10"]),
          frailty_category = purrr::map_chr(strata, ~.["frailty_category"]),
          pre_2020 = purrr::map_chr(strata, ~.["pre_2020"])) %>% 
   mutate(across(c("sex", "age_group", "polypharm_gte_5", "polypharm_gte_10", "frailty_category", "pre_2020"), ~tidyr::replace_na(., "NA"))) %>% 
   select(-strata_nm, -strata_lvl, -strata) %>% 
   as.data.frame() %>% 
   {DBI::dbWriteTable(res, name = "hospitalization_death", value = ., overwrite = T)}

# large scale characterization --------
# this could be done in one step but the result set is quite large so I split it into three files
log4r::info(logger, "Running large scale condition characterisation")

ls_conditions <- cdm$cancer_cohort %>%
  PatientProfiles::summariseLargeScaleCharacteristics(
    strata = list("age_group", 
                  "sex",
                  "polypharm_gte_5", 
                  "polypharm_gte_10", 
                  "frailty_category",
                  "pre_2020"),
    minimumFrequency = 0.01,
    window = list(c(-999999, -366),
                  c(-365, -31),
                  c(-30, -1),
                  c(0, 0)),
    episodeInWindow = c("condition_occurrence"))

log4r::info(logger, "saving large scale condition results")

ls_conditions %>% 
  PatientProfiles::suppress(minCellCount) %>% 
  as.data.frame() %>% 
  mutate(strata_nm = stringr::str_split(strata_name, " and "),
         strata_lvl = stringr::str_split(strata_level, " and "),
         strata = purrr::map2(strata_lvl, strata_nm, setNames)) %>% 
  mutate(age_group = purrr::map_chr(strata, ~.["age_group"]),
         sex = purrr::map_chr(strata, ~.["sex"]),
         polypharm_gte_5 = purrr::map_chr(strata, ~.["polypharm_gte_5"]),
         polypharm_gte_10 = purrr::map_chr(strata, ~.["polypharm_gte_10"]),
         frailty_category = purrr::map_chr(strata, ~.["frailty_category"]),
         pre_2020 = purrr::map_chr(strata, ~.["pre_2020"])) %>% 
  mutate(across(c("sex", "age_group", "polypharm_gte_5", "polypharm_gte_10", "frailty_category", "pre_2020"), ~tidyr::replace_na(., "NA"))) %>% 
  select(-strata_nm, -strata_lvl, -strata) %>% 
  {DBI::dbWriteTable(res, name = "ls_conditions", value = ., overwrite = T)}

rm(ls_conditions)

log4r::info(logger, "Running large scale drug characterisation")

ls_drugs <- cdm$cancer_cohort %>%
  PatientProfiles::summariseLargeScaleCharacteristics(
    strata = list("age_group", 
                  "sex",
                  "polypharm_gte_5", 
                  "polypharm_gte_10", 
                  "frailty_category",
                  "pre_2020"),
    minimumFrequency = 0.01,
    window = list(c(-999999, -366),
                  c(-365, -31),
                  c(-30, -1),
                  c(0, 0)),
    episodeInWindow = c("drug_exposure"))

log4r::info(logger, "saving large scale drug characterization results")

ls_drugs %>% 
  PatientProfiles::suppress(minCellCount) %>% 
  as.data.frame() %>% 
  mutate(strata_nm = stringr::str_split(strata_name, " and "),
         strata_lvl = stringr::str_split(strata_level, " and "),
         strata = purrr::map2(strata_lvl, strata_nm, setNames)) %>% 
  mutate(age_group = purrr::map_chr(strata, ~.["age_group"]),
         sex = purrr::map_chr(strata, ~.["sex"]),
         polypharm_gte_5 = purrr::map_chr(strata, ~.["polypharm_gte_5"]),
         polypharm_gte_10 = purrr::map_chr(strata, ~.["polypharm_gte_10"]),
         frailty_category = purrr::map_chr(strata, ~.["frailty_category"]),
         pre_2020 = purrr::map_chr(strata, ~.["pre_2020"])) %>% 
  mutate(across(c("sex", "age_group", "polypharm_gte_5", "polypharm_gte_10", "frailty_category", "pre_2020"), ~tidyr::replace_na(., "NA"))) %>% 
  select(-strata_nm, -strata_lvl, -strata) %>% 
  {DBI::dbWriteTable(res, name = "ls_drugs", value = ., overwrite = T)}

rm(ls_drugs)

log4r::info(logger, "Large scale characterization complete")

log4r::info(logger, "Results database contents ------")
for (t in DBI::dbListTables(res)) {
  n <- tbl(res, t) %>% tally() %>% pull()
  log4r::info(logger, paste(t, "table has", n, "rows"))
}

DBI::dbDisconnect(res, shutdown = T)

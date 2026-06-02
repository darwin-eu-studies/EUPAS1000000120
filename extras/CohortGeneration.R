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

# generate MM cohort ----
cli::cli_alert_info("Getting multiple myeloma cohort")
mm_codes <- CodelistGenerator::codesFromCohort(path = here("cohorts", "exposures"),
                                               cdm = cdm)
mm_codes <- unlist(mm_codes)

cdm <- CDMConnector::generateConceptCohortSet(cdm = cdm,
                                              conceptSet = list("Multiple myeloma" = mm_codes,
                                                                "Multiple myeloma (365 days prior observation)" = mm_codes,
                                                                "Multiple myeloma (30 days follow-up)" = mm_codes,
                                                                "Multiple myeloma (1 year of potential follow-up)" = mm_codes),
                                              end = "observation_period_end_date",
                                              limit = "first",
                                              name = "mm_cohort",
                                              overwrite = TRUE)

cohort_id <- cohortSet(cdm$mm_cohort) %>% collect()
mm_cohort_id <- cohort_id %>%
  filter(cohort_name == "Multiple myeloma") %>%
  pull("cohort_definition_id")
mm_p_cohort_id <- cohort_id %>%
  filter(cohort_name == "Multiple myeloma (365 days prior observation)") %>%
  pull("cohort_definition_id")
mm_f_cohort_id <- cohort_id %>%
  filter(cohort_name == "Multiple myeloma (30 days follow-up)") %>%
  pull("cohort_definition_id")
mm_surv_cohort_id <- cohort_id %>%
  filter(cohort_name == "Multiple myeloma (1 year of potential follow-up)") %>%
  pull("cohort_definition_id")

# MM exclusions  ----
cli::cli_alert_info("Running inclusion criteria for multiple myeloma cohort")
cdm$mm_cohort <- cdm$mm_cohort %>%
  addDemographics(ageGroup = list(c(0,17),
                                  c(18,44),
                                  c(45,59),
                                  c(60,69),
                                  c(70,150))) %>%
  mutate(study_period = if_else(cohort_start_date < as.Date("2018-01-01"),
                                "2017 or earlier", "2018 or later"))

# ensure first record, age and sex present
cdm$mm_cohort <- cdm$mm_cohort %>%
  dplyr::slice_min(n = 1,
                   order_by = cohort_start_date,
                   by = c("cohort_definition_id", "subject_id")) %>%
  compute_query() %>%
  filter(age >= 0 & age <=150,
         !is.na(sex),
         cohort_start_date >= as.Date("2012-01-01")) %>%
  compute_query() %>%
  record_cohort_attrition("MM cohort")

attr(cdm$mm_cohort, "cohort_attrition") <- attr(cdm$mm_cohort, "cohort_attrition") %>%
  filter(reason == "MM cohort") %>%
  mutate(reason_id = 1L) %>%
  compute_query()

# exclude history of cancer
exc_codes <- CodelistGenerator::codesFromConceptSet(path = here("cohorts", "exclusions"),
                                                    cdm = cdm)
exc_codes <- unlist(exc_codes)

cdm <- CDMConnector::generateConceptCohortSet(cdm = cdm,
                                              conceptSet = list("malignant_neoplasm" = exc_codes),
                                              end = "observation_period_end_date",
                                              limit = "first",
                                              name = "mm_exc_malignant_neoplasm",
                                              overwrite = TRUE)


cdm$mm_cohort <- cdm$mm_cohort %>%
  addCohortIntersectFlag(targetCohortTable = "mm_exc_malignant_neoplasm",
                         window = c(-Inf, -1), nameStyle = "mm_exc") %>%
  filter(mm_exc == 0) %>%
  select(!"mm_exc")%>%
  record_cohort_attrition("No prior diagnosis of cancer")
attr(cdm$mm_cohort, "cohort_attrition")

cdm$mm_cohort <- cdm$mm_cohort %>%
  filter(cohort_definition_id %in%
           c(mm_cohort_id, mm_f_cohort_id, mm_surv_cohort_id) |
           prior_observation > 365) %>%
  record_cohort_attrition("Year of prior observation")

cdm$mm_cohort <- cdm$mm_cohort %>%
  filter(cohort_definition_id %in%
           c(mm_cohort_id, mm_p_cohort_id, mm_surv_cohort_id) |
           future_observation >= 30) %>%
  record_cohort_attrition("30 days of follow-up")


db_year_left <- as.Date(snapshot(cdm)$latest_observation_period_end_date)  - years(1)
cdm$mm_cohort <- cdm$mm_cohort %>%
  filter(cohort_definition_id %in%
           c(mm_cohort_id, mm_f_cohort_id, mm_p_cohort_id) |
           cohort_start_date <=  db_year_left) %>%
  record_cohort_attrition("1 year of potential follow-up")

mm_attrition <- cohort_attrition(cdm$mm_cohort) %>%
  collect()
mm_set <- cohortSet(cdm$mm_cohort) %>%
  collect()

write_csv(mm_attrition %>%
            left_join(mm_set,
                      by = "cohort_definition_id") %>%
            mutate(cdm_name = cdmName(cdm)),
          here("Results", paste0("cohort_attrition_mm_cohort_",
                                 cdmName(cdm), ".csv")))

# generate comorbidity cohorts ---------
cli::cli_alert_info("Getting comorbidity cohorts")
comorbidity_cs <- CodelistGenerator::codesFromConceptSet(path = here("cohorts","comorbidity"),
                                                        cdm = cdm)
cdm <- CDMConnector::generateConceptCohortSet(cdm = cdm,
                                              name = "mm_cohort_cond",
                                              conceptSet = comorbidity_cs,
                                              limit = "first",
                                              end = "observation_period_end_date",
                                              overwrite = TRUE)
gc()
rm(comorbidity_cs)

# generate treatment cohorts -----
cli::cli_alert_info("Getting treatment cohorts")
treatments_cs <- CodelistGenerator::codesFromConceptSet(path = here("cohorts","treatments"),
                                                        cdm = cdm)
cdm <- DrugUtilisation::generateDrugUtilisationCohortSet(cdm = cdm,
                                                 name = "mm_cohort_treatments",
                                                 conceptSetList = treatments_cs)
gc()
rm(treatments_cs)

# generate medication cohorts -----
cli::cli_alert_info("Getting medication cohorts")
medication_cs <- CodelistGenerator::codesFromConceptSet(path = here("cohorts","medications"),
                                                        cdm = cdm)
cdm <- DrugUtilisation::generateDrugUtilisationCohortSet(cdm = cdm,
                                                         name = "mm_cohort_meds",
                                                         conceptSetList = medication_cs,
                                                         gapEra = 7)
gc()
rm(medication_cs)

cli::cli_alert_success("Cohorts created")

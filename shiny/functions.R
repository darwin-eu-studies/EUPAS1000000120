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


getFilters <- function(demo) {
  
  require(dplyr)
  databases <- demo %>% 
    distinct(.data$cdm_name) %>% 
    pull()
  
  cancers <- demo %>% 
    distinct(.data$group_level) %>% 
    pull()
  
  age_groups <- demo %>%
    filter(strata_name == "age_group") %>%
    distinct(strata_level) %>%
    filter(strata_level != "None") %>%
    pull() %>%
    {setNames(., .)} %>%
    sort() %>% 
    {c("All ages" = "NA", .)}
  
  sex <- demo %>%
    filter(strata_name == "sex") %>%
    distinct(strata_level) %>%
    pull() %>%
    {setNames(., .)} %>%
    {c("Both sexes" = "NA", .)}
  
  frailty <- demo %>%
    filter(strata_name == "frailty_category") %>%
    distinct(strata_level) %>%
    pull() %>%
    {setNames(., .)} %>%
    {c("All persons" = "NA", .)}
  
  polypharm <- c("All persons" = "NA", "0-4 drugs", ">=5 drugs", ">=10 drugs")
  pre_2020 <- c("All time" = "NA", "Pre 2020" = "TRUE", "2020 and later" = "FALSE")
  
  dplyr::lst(databases,
             cancers,
             age_groups,
             sex,
             frailty,
             polypharm,
             pre_2020)
}


# all distinct combinations of filters
getAllFilterCombinations <- function(demo) {
  filters <- getFilters(demo)
  
  expand.grid(databases = filters$database,
              cancers = filters$cancers,
              age_groups = filters$age_groups, 
              sex = filters$sex,
              frailty = filters$frailty,
              polypharm = filters$polypharm,
              pre_2020 = filters$pre_2020) %>% 
    tibble() %>% 
    mutate_all(as.character)
}

# get table 1
# cancer = "breast_cancer"
# polypharm = "NA"
# age = "NA"
# sex = "NA"
# frailty = "NA"
# pre_2020 = "NA"

add_value <- function(df, ci = T) {
  stopifnot(all(c("count", "percentage") %in% colnames(df)))
  
  suppressWarnings({
    df <- df %>% 
      mutate(count = as.numeric(.data$count), percentage = as.numeric(.data$percentage), p = percentage/100) %>% 
      mutate(denom = .data$count/.data$p, 
             z = 1.96*sqrt((.data$p*(1-.data$p))/denom),
             lcl = .data$p - z,
             ucl = .data$p + z
      ) %>% 
      mutate(value = case_when(
        is.na(.data$percentage) | is.na(.data$count) ~ "<supressed>",
        .data$percentage == "0" | isFALSE(ci) ~ as.character(glue::glue("{.data$count} (pct: {round(.data$percentage, 1)}%)")),
        TRUE ~ as.character(glue::glue("{.data$count} (pct: {round(.data$percentage, 1)}%, lcl: {round(.data$lcl*100, 1)}%, ucl: {round(.data$ucl*100, 1)}%)")))
      ) %>%
      select(-z, -lcl, -ucl, -denom, -p)
  })
  
  return(df)
}

# cancer = "breast_cancer"
# polypharm = "NA"
# age = "NA"
# sex = "NA"
# frailty = "NA"
# pre_2020 = "NA"

# table 1 ----

table1_filtered <- function(
    prev, 
    cancer = "breast_cancer",
    polypharm = "NA",
    age = "NA",
    sex = "NA",
    frailty = "NA",
    pre_2020 = "NA") {
  
  # polypharmacy stratification options.
  # Note that the polypharmacy strata do not partition the study population like other strata do.
  if (polypharm == "NA") {
    polypharm_gte_5 <- "NA"
    polypharm_gte_10 <- "NA"
  } else if (polypharm == "0-4 drugs") {
    polypharm_gte_5 <- "0"
    polypharm_gte_10 <- "NA"
  } else if (polypharm == ">=5 drugs") {
    polypharm_gte_5 <- "1"
    polypharm_gte_10 <- "NA"
  } else if (polypharm == ">=10 drugs") {
    polypharm_gte_5 <- "NA"
    polypharm_gte_10 <- "1"
  } else {
    stop("polypharm filter is not valid")
  }

  df <- prev %>%
    filter(
      group_level == !!cancer,
      age_group == !!age,
      sex == !!sex,
      polypharm_gte_5 == !!polypharm_gte_5,
      polypharm_gte_10 == !!polypharm_gte_10,
      frailty_category == !!frailty,
      pre_2020 == !!pre_2020
    ) %>%
    collect() 
  
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  median_age <- df %>% 
    filter(variable_name == "age", estimate_name == "median") %>%
    mutate(value = tidyr::replace_na(estimate_value, "NA"), variable = "Median Age") %>%
    select("cdm_name", "variable", "value")
  
  iqr_age <- df %>% 
    filter(variable_name == "age", estimate_name %in% c("q25", "q75")) %>% 
    select("cdm_name", "variable_name", "estimate_name", "estimate_value") %>% 
    tidyr::pivot_wider(names_from = "estimate_name", values_from = "estimate_value")
    
    if (all(c("q75", "q25") %in% colnames(iqr_age))) {
      # casting can produce warnings about type conversion
      suppressWarnings({
        iqr_age <- mutate(iqr_age, iqr = as.character(as.numeric(q75) - as.numeric(q25))) %>% 
          select("cdm_name", "iqr")
      })
    } else {
      iqr_age <- mutate(iqr_age, iqr = "NA") %>% 
        select("cdm_name", "iqr")
    }
  
  median_age <- median_age %>% 
    left_join(iqr_age, by = "cdm_name") %>% 
    mutate(value = paste0(value, " (IQR=", iqr, ")")) %>% 
    select("cdm_name", "variable", "value")
    
  # age <- df %>%
  #   filter(variable_name == "age_group", variable_level != "None") %>% 
  #   select(cdm_name, variable_name, variable_level, estimate_name, estimate_value) %>% 
  #   # dplyr::summarise(n = dplyr::n(), .by = c(cdm_name, variable_name, variable_level, estimate_name)) |>
  #   # dplyr::filter(n > 1L) 
  #   # add_count(cdm_name, variable_name, variable_level) %>% 
  #   # filter(n >1) %>% 
  #   # arrange(cdm_name, variable_name, variable_level) %>% View
  #   tidyr::pivot_wider(names_from = "estimate_name", values_from = "estimate_value") 
  # 
  # age %>% 
  #   count()
  
  age <- df %>%
    filter(variable_name == "age_group", variable_level != "None") %>%
    {if (nrow(.) > 0) {
      select(., cdm_name, variable_name, variable_level, estimate_name, estimate_value) %>%
        tidyr::pivot_wider(names_from = "estimate_name", values_from = "estimate_value") %>%
        add_value(ci = F) %>% 
        mutate(variable = glue::glue("{variable_name}: {variable_level}"))
    } else . } %>% 
    select("cdm_name", "variable", "value")
  
  sex <- df %>%
    filter(variable_name == "sex") %>%
    {if (nrow(.) > 0) {
      select(., cdm_name, variable_name, variable_level, estimate_name, estimate_value) %>%
        tidyr::pivot_wider(names_from = "estimate_name", values_from = "estimate_value") %>%
        add_value(ci = F) %>% 
        mutate(variable = glue::glue("{variable_name}: {variable_level}"))
        # mutate(value = as.character(glue::glue("{.data$count} ({round(as.numeric(.data$percentage), 1)}%)")),
        #        variable = glue::glue("{variable_name}: {variable_level}"))
    } else . } %>% 
    select("cdm_name", "variable", "value")
  
  frailty_df <- df %>%
    filter(variable_name == "score") %>% 
    select(cdm_name, variable_name, estimate_name, estimate_value) 
  
  if (nrow(frailty_df) > 0) {

    frailty_df_iqr <- df %>%
      filter(variable_name == "score", estimate_name %in% c("q25", "q75")) %>% 
      select(cdm_name, estimate_name, estimate_value)  %>%
      tidyr::pivot_wider(names_from = "estimate_name", values_from = "estimate_value") %>% 
      mutate(iqr = as.character(round(as.numeric(q75) - as.numeric(q25), 3))) %>% 
      select(cdm_name, iqr)
    
    frailty_df <- frailty_df %>%
      filter(estimate_name == "median") %>%
      transmute(cdm_name = cdm_name,
                variable = "Median frailty score",
                value = as.character(round(as.numeric(estimate_value), 3))) %>% 
      left_join(frailty_df_iqr, by = "cdm_name") %>% 
      mutate(value = paste0(value, " (IQR=", iqr, ")")) %>% 
      select(cdm_name, variable, value)
    
  } else {
    frailty_df <- tibble(cdm_name = character(),
                         variable = character(), 
                         value = character())
  }
  
  
  frailty_category <- df %>%
    filter(variable_name == "frailty_category") %>%
    {if (nrow(.) > 0) {
      select(., cdm_name, variable_name, variable_level, estimate_name, estimate_value) %>%
        tidyr::pivot_wider(names_from = "estimate_name", values_from = "estimate_value") %>%
        add_value(ci = T) %>% 
        # mutate(value = as.character(glue::glue("{.data$count} ({round(as.numeric(.data$percentage), 1)}%)")),
        mutate(variable = glue::glue("{variable_name}: {variable_level}"))
    } else .} %>% 
    select("cdm_name", "variable", "value")
  
  polypharm <- df %>%
    filter(variable_name != "frailty_category", stringr::str_detect(variable_name, "^polypharm_")) %>%
    {if (nrow(.) > 0) {
      select(., cdm_name, variable_name, estimate_name, estimate_value) %>%
        tidyr::pivot_wider(names_from = "estimate_name", values_from = "estimate_value") %>%
        add_value(ci = T) %>% 
        mutate(variable = variable_name)
        # mutate(value = as.character(glue::glue("{.data$count} ({round(as.numeric(.data$percentage), 1)}%)")),
               # variable = variable_name) 
    } else .} %>% 
    select("cdm_name", "variable", "value") 
  
  median_meds <- df %>% 
    filter(variable_name == "n_drug_ingredients", estimate_name == "median") %>%
    mutate(value = tidyr::replace_na(estimate_value, "NA"), variable = "Median number of medications") %>%
    select("cdm_name", "variable", "value")
  
  iqr_meds <- df %>% 
    filter(variable_name == "n_drug_ingredients", estimate_name %in% c("q75", "q25")) %>%
    select(cdm_name, estimate_name, estimate_value) %>% 
    tidyr::pivot_wider(names_from = "estimate_name", values_from = "estimate_value") %>% 
    mutate(iqr = as.character(as.numeric(q75) - as.numeric(q25))) %>% 
    select(cdm_name, iqr)
  
  median_meds <- median_meds %>% 
    left_join(iqr_meds, by = "cdm_name") %>% 
    mutate(value = paste0(value, " (IQR=", iqr, ")")) %>% 
    select("cdm_name", "variable", "value")
    
  
  # control sort order of comorbs.
  disease_state <- c("arthritis", "atrial_fibrillation", "c_kidney_dx", "cerebrovascular_dx", 
                     "diabetes", "foot_problem", "fragility_fracture", "heart_failure", "heart_valve_dx", "hypertension",
                     "hypotension", "ischemic_heart_dx", "osteoporosis", "parkinsonism_tremor", "peptic_ulcer", 
                     "peripheral_vascular_dx", "respiratory_dx", "skin_ulcer", "thyroid_dx", "urinary_system_disease")
  
  symptoms <- c("dyspnea", "dizziness", "falls", "memeroy_cogn_problem", "sleep_disturbance", 
                "urinary_incontinence", "weight_loss_anorexia_jp")
  
  disability <- c("activity_limitation", "hearing_impairment", "housebound", "mobility_transfer", 
                  "req_care", "social_vulnerability", "visual_impairment")
  
  abnormal_lab_value <- c("anemia")
  
  factor_levels <- paste0("frailty_", c(disease_state, symptoms, disability, abnormal_lab_value))
  
  comorbidities <- df %>%
    filter(variable_name != "frailty_category", stringr::str_detect(variable_name, "^frailty_")) %>%
    {if (nrow(.) > 0) {
      select(., cdm_name, variable_name, estimate_name, estimate_value) %>%
        mutate(variable_name = factor(variable_name, levels = factor_levels)) %>% 
        arrange(variable_name) %>% 
        tidyr::pivot_wider(names_from = "estimate_name", values_from = "estimate_value") %>%
        add_value(ci = F) %>% 
        # mutate(value = as.character(glue::glue("{.data$count} ({round(as.numeric(.data$percentage), 1)}%)")),
        mutate(variable = stringr::str_remove(as.character(variable_name), "frailty_")) %>% 
        mutate(variable = case_when(
          stringr::str_detect(variable_name, paste(disease_state, collapse = "|")) ~ paste("Disease State:", variable_name),
          stringr::str_detect(variable_name, paste(symptoms, collapse = "|")) ~ paste("Symptoms/Signs:", variable_name),
          stringr::str_detect(variable_name, paste(disability, collapse = "|")) ~ paste("Disability:", variable_name),
          stringr::str_detect(variable_name, paste(abnormal_lab_value, collapse = "|")) ~ paste("Abnormal Lab Value:", variable_name),
          T ~ variable_name
        ))
    } else .} %>% 
    select("cdm_name", "variable", "value") %>% 
    mutate(variable = stringr::str_replace_all(variable, ": frailty_", ": "))
  
  out <- dplyr::bind_rows(median_age, age, sex, frailty_df, frailty_category, median_meds, polypharm, comorbidities)
  
  if (nrow(out) == 0) NULL else tidyr::pivot_wider(out, names_from = "cdm_name", values_from = "value") 
}

# table1_filtered(prev)

# wrapper to handle selections of multiple sex groups
table1_filtered2 <- function(
    prev, 
    cancer = "breast_cancer",
    polypharm = "NA",
    age = "NA",
    sex = c("Male", "Female"),
    frailty = "NA",
    pre_2020 = "NA") {
  
  if (length(sex) > 1) {
    l <- list()
    for (i in seq_along(sex)) {
      l[[i]] <- table1_filtered(
        prev, 
        cancer = cancer,
        polypharm = polypharm,
        age = age,
        sex = sex[i],
        frailty = frailty,
        pre_2020 = pre_2020) %>% 
        rename_at(-1, ~paste0(., "_", ifelse(sex[i] == "NA", "All", sex[i])))
    }
    out <- purrr::reduce(l, left_join, by = "variable")
    out <- out[,order(colnames(out))] %>% 
      select(variable, everything())
    
  } else {
    out <- table1_filtered(
      prev, 
      cancer = cancer,
      polypharm = polypharm,
      age = age,
      sex = sex,
      frailty = frailty,
      pre_2020 = pre_2020)
  }
  out
}

# table1_filtered2(
#     prev, 
#     cancer = "breast_cancer",
#     polypharm = "NA",
#     age = "NA",
#     sex = c("Male", "Female"),
#     frailty = "NA",
#     pre_2020 = "NA")



# wrapper to handle selections of multiple age groups
table1_filtered3 <- function(
    prev, 
    cancer = "breast_cancer",
    polypharm = "NA",
    age = c("18 to 44", "45 to 64"),
    sex = c("Male", "Female"),
    frailty = "NA",
    pre_2020 = "NA") {
  
  if (length(age) == 0 || length(sex) == 0) return(NULL)
  
  if (length(age) > 1) {
    l <- list()
    for (i in seq_along(age)) {
      l[[i]] <- table1_filtered2(
        prev, 
        cancer = cancer,
        polypharm = polypharm,
        age = age[i],
        sex = sex,
        frailty = frailty,
        pre_2020 = pre_2020) %>% 
        mutate(age_group = ifelse(age[i] == "NA", "All ages", age[i])) %>% 
        select(age_group, everything())
    }
    out <- bind_rows(l) 
    
  } else {
    out <- table1_filtered2(
      prev, 
      cancer = cancer,
      polypharm = polypharm,
      age = age,
      sex = sex,
      frailty = frailty,
      pre_2020 = pre_2020)
  }
  out
}

# check that it works
# table1_filtered3(    prev, 
#                      cancer = "breast_cancer",
#                      polypharm = "NA",
#                      age = c("18 to 44", "45 to 64"),
#                      sex = c("Male", "Female"),
#                      frailty = "NA",
#                      pre_2020 = "NA")

# table1_filtered3(
#   prev,
#   cancer = "breast_cancer",
#   polypharm = "NA",
#   age = "NA",
#   sex = c("NA"),
#   frailty = "NA",
#   pre_2020 = "NA")

testTable1 <- function(prev) {
  
  df <- getAllFilterCombinations(demo)
  stopifnot(nrow(df) > 1)
  for (i in cli::cli_progress_along(1:nrow(df))) {
    table1_filtered(
      prev, 
      cancer = df$cancers[i],
      polypharm = df$polypharm[i],
      age = df$age_groups[i],
      sex = df$sex[i],
      frailty = df$frailty[i],
      pre_2020 = df$pre_2020[i])
  }
  
}

# testTable1(prev)

cdm_name = "ipci"
cancer = "breast_cancer"
age = "NA"
sex = "Male"
pre_2020 = "NA"


# table 2 -----
table2_filtered <- function(prev,
                            cdm_name = "ipci",
                            cancer = "breast_cancer",
                            age = "NA",
                            sex = "Male",
                            pre_2020 = "NA") {
    
  df <- prev %>%
    filter(
      cdm_name == !!cdm_name,
      group_level == !!cancer,
      age_group == !!age,
      sex == !!sex,
      pre_2020 == !!pre_2020
    ) %>%
    collect() 
  
  df2 <- df %>% 
    mutate(variable_level = tidyr::replace_na(variable_level, "")) %>% 
    mutate(variable = ifelse(variable_level == "", variable_name, glue::glue("{variable_name}: {variable_level}"))) %>% 
    mutate(variable = stringr::str_remove(variable, "^frailty_")) %>% 
    filter(stringr::str_detect(variable_name, "frailty_|age_group|sex"), 
           variable_level != "None", 
           variable_name != "frailty_category")
  
  if (nrow(df2) == 0) {
    return(NULL)
  }
  
  by_frailty <- df2 %>%
    filter(frailty_category != "NA", polypharm_gte_5 == "NA", polypharm_gte_10 == "NA") %>% 
    select(cdm_name, variable, estimate_name, estimate_value, frailty_category) %>%
    mutate(frailty_category = paste0("frailty_", frailty_category)) %>% 
    tidyr::pivot_wider(names_from = "estimate_name", values_from = "estimate_value") %>%
    add_value(ci = T) %>% 
    # mutate(value = as.character(glue::glue("{.data$count} ({round(as.numeric(.data$percentage), 1)}%)"))) %>% 
    select(-count, -percentage, -cdm_name) %>% 
    tidyr::pivot_wider(names_from = "frailty_category", values_from = "value")
  
  by_polypharm5 <- df2 %>% 
    filter(frailty_category == "NA", polypharm_gte_5 == "1", polypharm_gte_10 == "NA") %>% 
    select(cdm_name, variable, estimate_name, estimate_value, polypharm_gte_5) %>%
    tidyr::pivot_wider(names_from = "estimate_name", values_from = "estimate_value") %>%
    {
      if (nrow(.) > 0) {
        add_value(., ci = T) %>% 
        # mutate(., value = as.character(glue::glue("{.data$count} ({round(as.numeric(.data$percentage), 1)}%)"))) %>%
          select(variable, polypharm_gte5 = value)
      } else {
        tibble(variable = by_frailty$variable, polypharm_gte5 = "no data")}
    }
  
  by_polypharm10 <- df2 %>% 
    filter(frailty_category == "NA", polypharm_gte_5 == "NA", polypharm_gte_10 == "1") %>% 
    select(cdm_name, variable, estimate_name, estimate_value, polypharm_gte_5) %>%
    tidyr::pivot_wider(names_from = "estimate_name", values_from = "estimate_value") %>%
    {
      if (nrow(.) > 0) {
        add_value(., ci = T) %>% 
        # mutate(., value = as.character(glue::glue("{.data$count} ({round(as.numeric(.data$percentage), 1)}%)"))) %>%
          select(variable, polypharm_gte10 = value)
      } else {tibble(variable = by_frailty$variable, polypharm_gte10 = "no data")}
    }
  
  out <- by_frailty %>% 
    left_join(by_polypharm5, by = "variable") %>% 
    left_join(by_polypharm10, by = "variable")
  
  out
}

# table2_filtered(prev, 
#                 cdm_name = "ebb",
#                 cancer = "breast_cancer",
#                 age = "NA",
#                 sex = "Male",
#                 pre_2020 = "NA")

# check that all filter combinations work (don't throw error)
testTable2 <- function(prev, allFilterCombinations) {
  
  df <- getAllFilterCombinations(demo) %>% 
    distinct(databases, cancers, age_groups, sex, pre_2020)
  
  
  stopifnot(nrow(df) > 1)
  for (i in cli::cli_progress_along(1:nrow(df))) {
    
    table2_filtered(prev,
                    cdm_name = df$databases[i],
                    cancer = df$cancers[i],
                    age = df$age_groups[i],
                    sex = df$sex[i],
                    pre_2020 = df$pre_2020[i])
  }
  
}

# table 3 ----
cancer = "breast_cancer"
age = "18 to 44"
sex = "Male"
pre_2020 = "NA"

table3_filtered <- function(hosp,
                            cancer = "breast_cancer",
                            age = "18 to 44",
                            sex = "Male",
                            pre_2020 = "NA") {
  
  df <- hosp %>% 
    filter(group_level == !!cancer,
           age_group == !!age,
           sex == !!sex,
           pre_2020 == !!pre_2020) %>% 
    collect() %>% 
    mutate(
      mortality_pct = round(as.numeric(death_flag_percentage)/100, 3),
      mortality_se = sqrt(mortality_pct*(1-mortality_pct)/as.numeric(number_subjects_count)),
      mortality_lcl = round(mortality_pct - 1.96*mortality_se, 3),
      mortality_ucl = round(mortality_pct + 1.96*mortality_se, 3),
      hospitalization_rate = round(as.numeric(hospitalization_rate), 3),
      hospitalization_se = sqrt(as.numeric(n_hospitalizations_in_year_after_dx_sum)/(as.numeric(number_subjects_count)^2)),
      hospitalization_lcl = round(hospitalization_rate + (1.96*hospitalization_se), 3), 
      hospitalization_ucl = round(hospitalization_rate + (1.96*hospitalization_se), 3)
    ) %>% 
    select(cdm_name, 
           cancer = group_level, 
           age_group, 
           sex, 
           polypharm_gte_5,
           polypharm_gte_10, 
           frailty_category, 
           pre_2020, 
           mortality_pct,
           mortality_lcl, 
           mortality_ucl,
           hospitalization_rate,
           hospitalization_lcl,
           hospitalization_ucl,
           strata_name, 
           strata_level) %>% 
    mutate(one_year_mortality_risk = ifelse(is.na(mortality_pct), "", glue::glue("{mortality_pct} ({mortality_lcl}, {mortality_ucl})")), 
           one_year_hospitalzation_rate = ifelse(is.na(hospitalization_rate), "", glue::glue("{hospitalization_rate} ({hospitalization_lcl}, {hospitalization_lcl})"))) %>% 
    select(-hospitalization_rate, -age_group, -cancer, -sex, -pre_2020, -strata_name, -strata_level)
  
  by_frailty <- df %>% 
    filter(frailty_category != "NA", polypharm_gte_5 == "NA", polypharm_gte_10 == "NA") %>% 
    select(-polypharm_gte_5, -polypharm_gte_10) %>% 
    tidyr::pivot_longer(cols = c("one_year_mortality_risk", "one_year_hospitalzation_rate"), names_to = "rate", values_to = "value") %>% 
    mutate(frailty_category = paste0(frailty_category, "_frailty")) %>% 
    select(cdm_name, frailty_category, rate, value) %>% 
    tidyr::pivot_wider(names_from = "frailty_category", values_from = "value") 
  
  by_polypharm_gte_5 <- df %>% 
    filter(frailty_category == "NA", polypharm_gte_5 == "1", polypharm_gte_10 == "NA") %>% 
    select(cdm_name, one_year_mortality_risk, one_year_hospitalzation_rate) %>% 
    tidyr::pivot_longer(cols = c("one_year_mortality_risk", "one_year_hospitalzation_rate"), names_to = "rate", values_to = "polypharm_gte_5") 
  
  by_polypharm_gte_10 <- df %>% 
    filter(frailty_category == "NA", polypharm_gte_5 == "NA", polypharm_gte_10 == "1") %>% 
    select(cdm_name, one_year_mortality_risk, one_year_hospitalzation_rate) %>% 
    tidyr::pivot_longer(cols = c("one_year_mortality_risk", "one_year_hospitalzation_rate"), names_to = "rate", values_to = "polypharm_gte_10")
  
  result <- bind_rows(
              select(by_frailty, cdm_name, rate),
              select(by_polypharm_gte_5, cdm_name, rate),
              select(by_polypharm_gte_10, cdm_name, rate)) %>% 
    distinct() %>% 
    left_join(by_frailty, by = c("cdm_name", "rate")) %>% 
    left_join(by_polypharm_gte_5, by = c("cdm_name", "rate")) %>% 
    left_join(by_polypharm_gte_10, by = c("cdm_name", "rate")) %>% 
    # suppression filtering. 
    # I'm going to suppress hospitalization results for iqvia germany, belgium, 
    # and ipci because we don't have a good denominator. I'll also suppress mortality for 
    # iqvia germany and belgium since they do not have death capture. 
    filter(!stringr::str_detect(cdm_name, "iqvia")) %>% 
    filter(!(cdm_name == "ipci" & rate == "one_year_hospitalzation_rate"))
    
    # result %>% 
    #   tidyr::pivot_longer(cols = c(-1,-2)) %>% 
    #   mutate(new_value = tidyr::replace_na(value, "")) %>% 
    #   group_by(cdm_name, rate) %>% 
    #   mutate(max_value = max(new_value)) %>% 
    #   filter(max_value > 0) %>% 
    #   select(-new_value, -max_value) %>% 
    #   tidyr::pivot_wider(names_from = "name", values_from = "value")
    
  result
}

# hosp
# cancer = "breast_cancer"
# age = "18 to 44"
# sex = "Male"
# pre_2020 = "NA"
# 
# table3_filtered(
#   hosp,
#   cancer = "breast_cancer",
#   age = "18 to 44",
#   sex = "Male",
#   pre_2020 = "NA"
# )

# testTable3 <- function(hosp, allFilterCombinations) {
#   df <- getAllFilterCombinations(demo) %>% 
    # distinct(cancers, age_groups, sex, pre_2020)
#   
#   stopifnot(nrow(df) > 1)
#   for (i in cli::cli_progress_along(1:nrow(df))) {
#     table3_filtered(hosp,
#                     cancer = df$cancers[i],
#                     age = df$age_groups[i],
#                     sex = df$sex[i],
#                     pre_2020 = df$pre_2020[i])
#   }
# }

# top 20 from large scale characterization -----

top_20_from_ls <- function(df) {
  
  combos <- tidyr::expand_grid(rank = 1:20, 
                               cdm_name = unique(df$cdm_name), 
                               time_window = unique(df$timeframe), 
                               cohort = unique(df$cohort))
  
  df <- df %>%
    collect() %>% 
    mutate(percent = as.numeric(percent)) %>% 
    filter(!is.na(percent), covariate != "No matching concept") %>% 
    select(cdm_name, cohort, covariate, timeframe, count, percent) %>% 
    group_by(cdm_name, cohort, timeframe) %>% 
    arrange(cdm_name, cohort, timeframe, desc(percent)) %>% 
    mutate(rank = row_number()) %>% 
    mutate(text = paste0(covariate, " [n=", count, ", ", percent, "%]")) %>% 
    select(rank, cdm_name, cohort, time_window = timeframe, text) %>% 
    left_join(combos, ., by = c("rank", "cdm_name", "cohort", "time_window")) %>% 
    mutate(text = tidyr::replace_na(text, "")) %>% 
    # add_count(rank, cdm_name, time_window, cohort, cdm_name) |>
    # dplyr::filter(n > 1L) %>% 
    arrange(time_window, cohort, rank, cdm_name, text) %>% 
    tidyr::pivot_wider(names_from = "cdm_name", values_from = "text") %>% 
    mutate(rank = factor(rank), 
           cohort = factor(cohort), 
           time_window = factor(time_window)) %>% 
    arrange(desc(rank))
  
  if (nrow(df) == 0) return(NULL) else return(df)
}

# plotting ----
get_plot <- function(prev, plot_number = 1, age_group = "NA") {
  
  stopifnot(length(age_group) == 1,
                   unname(age_group) %in% c("NA", "18 to 44", "45 to 64", "65 to 74", "75 to 84", "85 to 120"))
  
  age_group <- unname(age_group)
  
  # 1. Overall prevalence of frailty conditions by frailty category among individuals with selected cancer types
  # 2. Prevalence of frailty categories by cancer type
  # 3. Prevalence of polypharmacy by cancer type
  # 4. Prevalence of frailty conditions by frailty category and Solid cancer type 
  # 5. Prevalence of frailty conditions by frailty category and Blood cancer type 
  
  if (age_group == "NA") {
    print('overall age group')
    df <- prev %>%
      filter(strata_name == "frailty_category") %>% 
      select(cdm_name, variable_name, estimate_name, estimate_value, cohort = group_level, frailty_category) %>%
      collect() %>% 
      filter(stringr::str_detect(variable_name, "frailty"), variable_name != "frailty_category") %>% 
      tidyr::pivot_wider(names_from = "estimate_name", values_from = "estimate_value") %>% 
      mutate(count = as.numeric(count), percentage = as.numeric(percentage)) %>% 
      mutate(variable_name = stringr::str_remove(variable_name, "frailty_"))
    
    age_group <- "Overall"
    
  } else {
    print(glue::glue("age group: {age_group}"))
    # age_group = "65 to 74"
    df <- prev %>%
      filter(strata_name == "age_group and frailty_category", age_group == .env$age_group) %>% 
      select(cdm_name, variable_name, estimate_name, estimate_value, cohort = group_level, frailty_category) %>%
      collect() %>% 
      filter(stringr::str_detect(variable_name, "frailty"), variable_name != "frailty_category") %>% 
      tidyr::pivot_wider(names_from = "estimate_name", values_from = "estimate_value") %>% 
      mutate(count = as.numeric(count), percentage = as.numeric(percentage)) %>% 
      mutate(variable_name = stringr::str_remove(variable_name, "frailty_"))
    
  }
  
  # df %>% distinct(variable_name) %>% pull %>% sort
  
  # cancer groups
  solid <- c("breast_cancer", "colorectal_cancer", "endometrial_cancer", "lung_cancer", 
             "ovarian_cancer", "pancreatic_cancer", "prostate_cancer")
  
  blood <- c("leukemia","lymphoma", "multiple_myeloma")
  
  # frailty groups
  disease_state <- c("arthritis", "atrial_fibrillation", "c_kidney_dx", "cerebrovascular_dx", 
                     "diabetes", "foot_problem", "fragility_fracture", "heart_failure", "heart_valve_dx", "hypertension",
                     "hypotension", "ischemic_heart_dx", "osteoporosis", "parkinsonism_tremor", "peptic_ulcer", 
                     "peripheral_vascular_dx", "respiratory_dx", "skin_ulcer", "thyroid_dx", "urinary_system_disease")
  
  symptoms <- c("dyspnea", "dizziness", "falls", "memeroy_cogn_problem", "sleep_disturbance", 
                "urinary_incontinence", "weight_loss_anorexia_jp")
  
  disability <- c("activity_limitation", "hearing_impairment", "housebound", "mobility_transfer", 
                  "req_care", "social_vulnerability", "visual_impairment")
  
  abnormal_lab_value <- c("anemia")
  
  library(ggplot2)
  
  # factor_levels <- df %>%
  #   group_by(variable_name) %>% 
  #   summarise(pct = mean(percentage, na.rm = T)) %>% 
  #   mutate(pct = case_when(
  #     variable_name %in% disease_state ~ pct + 100,
  #     variable_name %in% symptoms ~ pct + 10,
  #     variable_name %in% symptoms ~ pct + 1,
  #     TRUE ~ pct
  #   )) %>% 
  #   arrange(desc(pct)) %>% 
  #   pull(variable_name)
  
  factor_levels <- c(disease_state, symptoms, disability, abnormal_lab_value)
  
  cdm_names <- unique(df$cdm_name)
  
  if (plot_number == 1) {
    
    plt1 <- df %>% 
      filter(cdm_name %in% cdm_names[1:3]) %>% 
      mutate(variable_name = factor(variable_name, levels = factor_levels)) %>% 
      filter(cohort == "overall_ca_all_cancer_types_included_in_study") %>% 
      ggplot(aes(x = percentage, y = variable_name)) +
      geom_point() +
      facet_grid(cdm_name ~ frailty_category) +
      ggtitle(glue::glue("All cancers combined, age group {age_group}")) +
      labs(y = "")
    
    plt2 <- df %>% 
      filter(cdm_name %in% cdm_names[-1:-3]) %>% 
      mutate(variable_name = factor(variable_name, levels = factor_levels)) %>% 
      filter(cohort == "overall_ca_all_cancer_types_included_in_study") %>% 
      ggplot(aes(x = percentage, y = variable_name)) +
      geom_point() +
      facet_grid(cdm_name ~ frailty_category) +
      ggtitle(glue::glue("All cancers combined, age group {age_group}")) +
      labs(y = "")
    
    plt <- gridExtra::grid.arrange(plt1, plt2, ncol = 1, nrow = 2)
    
    return(plt)
  } 
  
  if (plot_number == 2) {
    solid <- c("breast_cancer", "colorectal_cancer", "endometrial_cancer", 
               "ovarian_cancer", "pancreatic_cancer", "prostate_cancer", 
               "lung_cancer", "sclc", "nsclc")
    
    blood <- c("lymphoma", "non_hk_lymphoma", "hk_lymphoma", "leukemia", "aml", "cml", "lla", "llc", "multiple_myeloma")
    
    cancer_levels <- c(solid, blood)
    
    if (age_group == "Overall") {
      p2_df <- prev %>% 
        filter(strata_name == "overall", variable_name %in% c("frailty_category")) %>% 
        collect()
    } else {
      p2_df <- prev %>% 
        filter(strata_name == "age_group", 
               strata_level == .env$age_group, 
               variable_name == "frailty_category") %>% 
        collect()
    }
    
    plt <- p2_df %>% 
      select(cdm_name, cohort = group_level, variable_name, variable_level, estimate_name, estimate_value) %>% 
      tidyr::pivot_wider(names_from = "estimate_name", values_from = "estimate_value") %>% 
      filter(cohort %in% cancer_levels) %>% 
      mutate(percentage = as.numeric(percentage),
             `Frailty category` = factor(variable_level, levels = rev(c("fit", "mild", "moderate", "severe"))),
             cohort = factor(cohort, levels = rev(cancer_levels))) %>% 
      # group_by(cohort, cdm_name) %>% 
      # summarise(pct = sum(percentage))
      ggplot(aes(y = percentage, x = cohort, fill = `Frailty category`)) + 
      geom_bar(position="stack", stat="identity") +
      # scale_fill_brewer(direction = 1) +
      facet_wrap("cdm_name") +
      coord_flip() +
      theme_bw() +
      ggtitle(glue::glue("Frailty category, age group {age_group}"))
    
    return(plt)
  }
  
  if (plot_number == 3) {
    solid <- c("breast_cancer", "colorectal_cancer", "endometrial_cancer", 
               "ovarian_cancer", "pancreatic_cancer", "prostate_cancer", 
               "lung_cancer", "sclc", "nsclc")
    
    blood <- c("lymphoma", "non_hk_lymphoma", "hk_lymphoma", "leukemia", "aml", 
               "cml", "lla", "llc", "multiple_myeloma")
    
    cancer_levels <- c(solid, blood)
    
    if (age_group == "Overall") {
      p3_df <- prev %>% 
        filter(strata_name %in% c("polypharm_gte_5", "polypharm_gte_10"), 
               variable_name == "number subjects") %>% 
        collect()
        
    } else {
      p3_df <- prev %>% 
        filter(strata_name %in% c("age_group and polypharm_gte_5", "age_group and polypharm_gte_10"), 
               variable_name == "number subjects",
               age_group == .env$age_group) %>% 
        collect()
    }
    
    plt <- p3_df %>% 
      mutate(strata_name = stringr::str_remove(strata_name, "age_group and "),
             strata_level = stringr::str_extract(strata_level, "[01]$")) %>% 
      tidyr::unite(col = "strata", strata_name, strata_level) %>% 
      select(cdm_name, cohort = group_level, strata, estimate_value) %>% 
      mutate(estimate_value = as.numeric(estimate_value)) %>% 
      tidyr::pivot_wider(names_from = "strata", values_from = "estimate_value") %>% 
      mutate(pct_polypharm_gte_5 = 100 * polypharm_gte_5_1 / (polypharm_gte_5_1 + polypharm_gte_5_0)) %>% 
      mutate(pct_polypharm_gte_10 = 100 * polypharm_gte_10_1 / (polypharm_gte_10_1 + polypharm_gte_10_0)) %>% 
      select(cdm_name, cohort, pct_polypharm_gte_5, pct_polypharm_gte_10) %>% 
      tidyr::pivot_longer(cols = c(pct_polypharm_gte_5, pct_polypharm_gte_10), names_to = "Polypharmacy", values_to = "Percent") %>% 
      mutate(Polypharmacy = case_when(
        Polypharmacy == "pct_polypharm_gte_5" ~ ">= 5 drugs",
        Polypharmacy == "pct_polypharm_gte_10" ~ ">= 10 drugs",
        T ~ "missing"
      )) %>% 
      filter(cohort %in% cancer_levels) %>% 
      mutate(cohort = factor(cohort, levels = rev(cancer_levels))) %>% 
      ggplot(aes(x = cohort, y = Percent, fill = Polypharmacy)) +
      geom_bar(stat = "identity", position = "dodge") +
      coord_flip() +
      facet_wrap("cdm_name") +
      theme_bw() +
      ggtitle(glue::glue("Polypharmacy, age group {age_group}"))
    
    return(plt)
  }
  
  if (plot_number == 4) {
    
    plt1 <- df %>% 
      filter(cdm_name %in% cdm_names[1:3]) %>% 
      mutate(variable_name = factor(variable_name, levels = factor_levels)) %>% 
      filter(cohort %in% solid) %>% 
      ggplot(aes(x = percentage, y = variable_name, color = cohort)) +
      geom_point() +
      facet_grid(cdm_name ~ frailty_category) +
      ggtitle("Solid tumors", subtitle = glue::glue("age group {age_group}")) +
      labs(y = "")
    
    
    plt2 <- df %>% 
      filter(cdm_name %in% cdm_names[-1:-3]) %>% 
      mutate(variable_name = factor(variable_name, levels = factor_levels)) %>% 
      filter(cohort %in% solid) %>% 
      ggplot(aes(x = percentage, y = variable_name, color = cohort)) +
      geom_point() +
      facet_grid(cdm_name ~ frailty_category) +
      ggtitle("Solid tumors", subtitle = glue::glue("age group {age_group}")) +
      labs(y = "")
    
    plt <- gridExtra::grid.arrange(plt1, plt2, ncol = 1, nrow = 2)
    
    return(plt)
  }
  
  if (plot_number == 5) {
    plt1 <- df %>% 
      filter(cdm_name %in% cdm_names[1:3]) %>% 
      mutate(variable_name = factor(variable_name, levels = factor_levels)) %>% 
      filter(cohort %in% blood) %>% 
      ggplot(aes(x = percentage, y = variable_name, color = cohort)) +
      geom_point() +
      facet_grid(cdm_name ~ frailty_category) +
      ggtitle("Blood cancers", subtitle = glue::glue("age group {age_group}")) +
      labs(y = "")
    
    plt2 <- df %>% 
      filter(cdm_name %in% cdm_names[-1:-3]) %>% 
      mutate(variable_name = factor(variable_name, levels = factor_levels)) %>% 
      filter(cohort %in% blood) %>% 
      ggplot(aes(x = percentage, y = variable_name, color = cohort)) +
      geom_point() +
      facet_grid(cdm_name ~ frailty_category) +
      ggtitle("Blood cancers", subtitle = glue::glue("age group {age_group}")) +
      labs(y = "")
    
    plt <- gridExtra::grid.arrange(plt1, plt2, ncol = 1, nrow = 2)
    
    return(plt)
  } 
  
  
  solid <- c("breast_cancer", "colorectal_cancer", "endometrial_cancer", 
             "ovarian_cancer", "pancreatic_cancer", "prostate_cancer", 
             "lung_cancer")
  blood <- c("lymphoma", "leukemia", "multiple_myeloma")
  cancer_levels <- c(solid, blood)
  databases <- c("ebb", "sidiap", "ipci", "cdm_gold_202307")
  
  if (plot_number == 6) {
    
    if (age_group == "Overall") {
      p6_df <- hosp %>% 
        filter(strata_name == "frailty_category") %>% 
        collect()
    } else {
      p6_df <- hosp %>% 
        filter(strata_name == "age_group and frailty_category",
               age_group == .env$age_group) %>% 
        collect()
    }
    
    # hospitalization
    plt <- p6_df %>% 
      select(cdm_name, cohort = group_level, 
             strata_name, 
             strata_level, 
             frailty_category,
             percent_hospitalized = hospitalization_rate) %>% 
      filter(cdm_name %in% c("sidiap", "ebb")) %>%
      mutate(cdm_name = factor(cdm_name, levels = databases)) %>% 
      filter(cohort %in% cancer_levels) %>% 
      mutate(frailty_category = factor(frailty_category, levels = c("fit", "mild", "moderate", "severe")),
             cohort = factor(cohort, levels = rev(cancer_levels))) %>% 
      ggplot(aes(y = percent_hospitalized, x = cohort, fill = cdm_name)) + 
      geom_bar(position="dodge", stat="identity") +
      facet_wrap("frailty_category") +
      coord_flip() +
      theme_bw() +
      labs(y = "Hospitalization rate within one year of index", x = "") +
      ggtitle("Hospitalization rate by frailty category", subtitle = glue::glue("age group {age_group}"))
    
    # ggsave("hosp_by_frailty.png", plt, height = 10, width = 10)
    return(plt)
  }
  
  if (plot_number == 7) {
    
    if (age_group == "Overall") {
      p7_df <- hosp %>% 
        filter(strata_name %in% c("polypharm_gte_10", "polypharm_gte_5"))
    } else {
      p7_df <- hosp %>% 
        filter(strata_name %in% c("age_group and polypharm_gte_10", "age_group and polypharm_gte_5"),
               age_group == .env$age_group)
    }
    
    plt <- p7_df %>% 
      select(cdm_name, cohort = group_level, strata_name, strata_level, percent_hospitalized = hospitalization_rate) %>% 
      filter(cdm_name %in% c("sidiap", "ebb"), stringr::str_detect(strata_level, "1$")) %>%
      filter(cohort %in% cancer_levels) %>% 
      mutate(strata_name = stringr::str_extract(strata_name, "polypharm_gte_5|polypharm_gte_10")) %>% 
      mutate(strata_name = factor(strata_name, levels = c("polypharm_gte_5", "polypharm_gte_10"))) %>% 
      mutate(cdm_name = factor(cdm_name, levels = databases)) %>% 
      mutate(cohort = factor(cohort, levels = rev(cancer_levels))) %>% 
      ggplot(aes(y = percent_hospitalized, x = cohort, fill = cdm_name)) + 
      geom_bar(position="dodge", stat="identity") +
      facet_wrap("strata_name") +
      coord_flip() +
      theme_bw() +
      labs(y = "Hospitalization rate within one year of index", x = "",) +
      ggtitle("Hospitalization rate by polypharmacy category", subtitle = glue::glue("age group {age_group}"))
    # ggsave("hosp_by_polypharm.png", plt, height = 10, width = 10)
    return(plt)
  }
  
  if (plot_number == 8) {
    
    if (age_group == "Overall") {
      p8_df <- hosp %>% 
        filter(strata_name == "frailty_category")
    } else {
      p8_df <- hosp %>% 
        filter(strata_name == "age_group and frailty_category",
               age_group == .env$age_group)
    }
    
    plt <- p8_df %>%
      select(cdm_name, cohort = group_level, strata_name, strata_level, percent_died = death_flag_percentage) %>% 
      mutate(percent_died = as.numeric(percent_died)) %>% 
      filter(cdm_name %in% c("sidiap", "ebb", "ipci", "cdm_gold_202307")) %>%
      mutate(cdm_name = factor(cdm_name, levels = databases)) %>% 
      filter(cohort %in% cancer_levels) %>%
      mutate(strata_level = stringr::str_extract(strata_level, "fit|mild|moderate|severe")) %>% 
      mutate(frailty_category = factor(strata_level, levels = c("fit", "mild", "moderate", "severe"))) %>% 
      mutate(cohort = factor(cohort, levels = rev(cancer_levels))) %>%
      ggplot(aes(y = percent_died, x = cohort, fill = cdm_name)) + 
      geom_bar(position="dodge", stat="identity") +
      facet_wrap("frailty_category") +
      coord_flip() +
      theme_bw() +
      labs(y = "Mortality risk within one year of cancer diagnosis", x = "") +
      ggtitle("Mortality risk by frailty category", subtitle = glue::glue("age group {age_group}"))
    
    # ggsave("mortality_by_frailty.png", plt, height = 10, width = 10)
    return(plt)
  }
  
  if (plot_number == 9) {
    
    if (age_group == "Overall") {
      p9_df <- hosp %>% 
        filter(strata_name %in% c("polypharm_gte_10", "polypharm_gte_5"))
    } else {
      p9_df <- hosp %>% 
        filter(strata_name %in% c("age_group and polypharm_gte_10", "age_group and polypharm_gte_5"),
               age_group == .env$age_group)
    }
    
    plt <- p9_df %>%
      select(cdm_name, cohort = group_level, strata_name, strata_level, percent_died = death_flag_percentage) %>% 
      mutate(percent_died = as.numeric(percent_died)) %>% 
      filter(cdm_name %in% c("sidiap", "ebb", "ipci", "cdm_gold_202307"), stringr::str_detect(strata_level, "1$")) %>%
      mutate(cdm_name = factor(cdm_name, levels = databases)) %>% 
      filter(cohort %in% cancer_levels) %>%
      mutate(polypharm = stringr::str_extract(strata_name, "polypharm_gte_10|polypharm_gte_5")) %>% 
      mutate(polypharm = factor(polypharm, levels = c("polypharm_gte_10", "polypharm_gte_5"))) %>% 
      mutate(cohort = factor(cohort, levels = rev(cancer_levels))) %>%
      ggplot(aes(y = percent_died, x = cohort, fill = cdm_name)) + 
      geom_bar(position="dodge", stat="identity") +
      facet_wrap("polypharm") +
      coord_flip() +
      theme_bw() +
      labs(y = "Mortality risk within one year of cancer diagnosis", x = "") +
      ggtitle("Mortality risk by polypharmacy category", subtitle = glue::glue("age group {age_group}"))
    
    # ggsave("mortality_by_polypharm.png", plt, height = 10, width = 10)
    return(plt)
  }
}


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

 con <- DBI::dbConnect(duckdb::duckdb(), here::here("results", "results_combined.duckdb"))

DBI::dbListTables(con)

library(dplyr)

demo <- tbl(con, "demographics") |> collect()
prev <- tbl(con, "prevalence") 

# need counts by cdm, age_group, frailty category

demo |> 
  distinct(strata_name, strata_level) |> 
  print(n=100)

# Figure 1: frailty by age, cdm ----  

frailty_by_age <- demo |> 
  filter(strata_name == "age_group and frailty_category") |> 
  filter(variable_name == "Number subjects") |> 
  filter(group_level == "overall_ca_all_cancer_types_included_in_study") |> 
  # distinct(group_level)
  select(cdm_name, strata_level, count = estimate_value) |> 
  tidyr::separate(col = "strata_level", into = c("age_group", "frailty_category"), sep = " and ") |> 
  mutate(count = as.numeric(count)) |> 
  filter(age_group != "None") |> 
  mutate(frailty_category = stringr::str_to_title(frailty_category)) |> 
  mutate(cdm_name = case_when(
    stringr::str_detect(cdm_name, "gold") ~ "CPRD-GOLD",
    stringr::str_detect(cdm_name, "germany") ~ "IQVIA Germany",
    stringr::str_detect(cdm_name, "belgium") ~ "IQVIA Belgium",
    T ~ toupper(cdm_name)
  )) |> 
  # filter(stringr::str_detect(age_group, "18 to|45 to", negate = T)) |> 
  group_by(cdm_name, age_group) |> 
  mutate(percent = count/sum(count, na.rm = T))

frailty_overall <- demo |>
  filter(strata_name == "frailty_category") |>
  filter(variable_name == "Number subjects") |>
  filter(group_level == "overall_ca_all_cancer_types_included_in_study") |>
  select(cdm_name, frailty_category = strata_level, count = estimate_value) |>
  mutate(count = as.numeric(count)) |>
  mutate(frailty_category = stringr::str_to_title(frailty_category)) |>
  mutate(cdm_name = case_when(
    stringr::str_detect(cdm_name, "gold") ~ "CPRD-GOLD",
    stringr::str_detect(cdm_name, "germany") ~ "IQVIA Germany",
    stringr::str_detect(cdm_name, "belgium") ~ "IQVIA Belgium",
    T ~ toupper(cdm_name)
  )) |>
  group_by(cdm_name) |>
  mutate(percent = count/sum(count, na.rm = T)) |>
  mutate(age_group = "Overall") |>
  select(cdm_name, age_group, frailty_category, count, percent)


# frailty_by_age <- bind_rows(frailty_overall, frailty_by_age)
readr::write_csv(frailty_by_age, "code_for_publication/csv/frailty_by_age.csv")

library(ggplot2)
frailty_by_age_plot <- frailty_by_age |> 
  mutate(frailty_category = factor(frailty_category, levels = c("Fit", "Mild", "Moderate", "Severe"))) |>
  # mutate(age_group = factor(age_group, levels = c("Overall", "65 to 74", "75 to 84", "85 to 120"), labels = c("Overall", "65 to 74", "75 to 84", ">=85"))) |> 
  mutate(age_group = factor(age_group, levels = c("18 to 44", "45 to 64", "65 to 74", "75 to 84", "85 to 120"), labels = c("18 to 44", "45 to 64", "65 to 74", "75 to 84", ">=85"))) |> 
  ggplot(aes(x = age_group, y = percent, fill = frailty_category)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~cdm_name, scales = "free_y") +
  coord_flip() +
  labs(
    x = "Age Group",
    y = "Percent",
    title = "Frailty Category by Age Groups & Database",
    fill = "Frailty Category"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA)    
  )

ggsave("code_for_publication/plots/frailty_by_age_plot.jpg", frailty_by_age_plot, height = 8, width = 8)


# figure 1b: polypharmacy by age

polypharm_by_age <- demo |> 
  filter(strata_name %in% c("age_group and polypharm_gte_5", "age_group and polypharm_gte_10")) |> 
  filter(variable_name == "Number subjects") |> 
  filter(group_level == "overall_ca_all_cancer_types_included_in_study") |> 
  # distinct(group_level)
  select(cdm_name, strata_name, strata_level, count = estimate_value) |> 
  tidyr::separate(col = "strata_level", into = c("age_group", "polypharm_level"), sep = " and ") |> 
  tidyr::separate(col = "strata_name", into = c("age_group_name", "polypharm_name"), sep = " and ") |> 
  mutate(polypharmacy = case_when(
    polypharm_name == "polypharm_gte_5" & polypharm_level == "0" ~ "<5",
    polypharm_name == "polypharm_gte_5" & polypharm_level == "1" ~ ">=5",
    polypharm_name == "polypharm_gte_10" & polypharm_level == "1" ~ ">=10",
    T ~ "drop"
  )) |> 
  filter(polypharmacy != "drop") |> 
  select(-polypharm_name, -polypharm_level, -age_group_name) |> 
  mutate(count = as.numeric(count)) |> 
  filter(age_group != "None") |> 
  mutate(cdm_name = case_when(
    stringr::str_detect(cdm_name, "gold") ~ "CPRD-GOLD",
    stringr::str_detect(cdm_name, "germany") ~ "IQVIA Germany",
    stringr::str_detect(cdm_name, "belgium") ~ "IQVIA Belgium",
    T ~ toupper(cdm_name)
  )) |> 
  filter(stringr::str_detect(age_group, "18 to|45 to", negate = T)) |> 
  group_by(cdm_name, age_group) |> 
  mutate(percent = count/sum(count, na.rm = T)) |> 
  mutate(polypharmacy = factor(polypharmacy, levels = c("<5", ">=5", ">=10")))


library(ggplot2)
polypharmacy_by_age_plot <- polypharm_by_age |> 
  ggplot(aes(x = age_group, y = percent, fill = polypharmacy)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~cdm_name, scales = "free_y") +
  coord_flip() +
  labs(
    x = "Age Group",
    y = "Percent",
    title = "Polypharmacy by Age Groups & Database",
    fill = "Polypharmacy"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA)    
  )

ggsave("code_for_publication/plots/polypharmacy_by_age_plot.jpg", polypharmacy_by_age_plot, height = 8, width = 8)

# Figure 2a: One year hospitalization rate by frailty, age group, and cdm ----

DBI::dbListTables(con)

library(dplyr)

hosp <- tbl(con, "hospitalization_death") |> collect()

poisson_ci <- function(events, people, conf.level = 0.95) {
  if (is.na(events)) return(c(NA_real_, NA_real_))
  ci <- poisson.test(events, T = people, conf.level = conf.level)$conf.int
  return(c(lower = ci[1], upper = ci[2]))
}

# hosp |> 
#   filter(strata_name == "polypharm_gte_5") |> 
#   filter(group_level == "overall_ca_all_cancer_types_included_in_study") |> 
#   select(cdm_name, strata_level, number_subjects_count, n_hosp = n_hospitalizations_in_year_after_dx_sum, death_flag_count) 
  
hosp_death_by_frailty_age_overall <- hosp |> 
  filter(strata_name == "frailty_category") |> 
  filter(group_level == "overall_ca_all_cancer_types_included_in_study") |> 
  select(cdm_name, strata_level, number_subjects_count, n_hosp = n_hospitalizations_in_year_after_dx_sum, death_flag_count) |> 
  mutate(age_group = "Overall", frailty_category = strata_level) |> 
  mutate(across(matches("number|n_|count"), as.numeric)) |> 
  filter(!is.na(number_subjects_count)) |> 
  mutate(frailty_category = stringr::str_to_title(frailty_category)) |> 
  mutate(cdm_name = case_when(
    stringr::str_detect(cdm_name, "gold") ~ "CPRD-GOLD",
    stringr::str_detect(cdm_name, "germany") ~ "IQVIA Germany",
    stringr::str_detect(cdm_name, "belgium") ~ "IQVIA Belgium",
    T ~ toupper(cdm_name)
  )) |> 
  select(cdm_name, age_group, frailty_category, number_subjects_count, n_hosp, death_flag_count)

hosp_death_by_frailty_age <- hosp |> 
  filter(strata_name == "age_group and frailty_category") |> 
  filter(group_level == "overall_ca_all_cancer_types_included_in_study") |> 
  select(cdm_name, strata_level, number_subjects_count, n_hosp = n_hospitalizations_in_year_after_dx_sum, death_flag_count) |> 
  tidyr::separate(col = "strata_level", into = c("age_group", "frailty_category"), sep = " and ") |> 
  mutate(across(matches("number|n_|count"), as.numeric)) |> 
  filter(!is.na(number_subjects_count)) |> 
  filter(age_group != "None") |> 
  mutate(frailty_category = stringr::str_to_title(frailty_category)) |> 
  mutate(cdm_name = case_when(
    stringr::str_detect(cdm_name, "gold") ~ "CPRD-GOLD",
    stringr::str_detect(cdm_name, "germany") ~ "IQVIA Germany",
    stringr::str_detect(cdm_name, "belgium") ~ "IQVIA Belgium",
    T ~ toupper(cdm_name)
  )) |> 
  filter(stringr::str_detect(age_group, "18 to|45 to", negate = T)) |> 
  bind_rows(hosp_death_by_frailty_age_overall)

# add ci
hosp_death_by_frailty_age <- hosp_death_by_frailty_age |> 
  mutate(
    death_pct = 100 * death_flag_count / number_subjects_count,
    ci = purrr::map2(death_flag_count, number_subjects_count, function(x, n) {
      if (is.na(x)) return(c(NA_real_, NA_real_))
      res <- binom::binom.confint(x, n, conf.level = 0.95, methods = "wilson")
      c(res$lower * 100, res$upper * 100)
    }),
    death_pct_lcl = purrr::map_dbl(ci, 1),
    death_pct_ucl = purrr::map_dbl(ci, 2)
  ) %>%
  select(-ci) |> 
  mutate(
    hosp_rate = n_hosp / number_subjects_count,  
    ci = purrr::map2(n_hosp, number_subjects_count, poisson_ci),
    hosp_rate_lcl = purrr::map_dbl(ci, 1),
    hosp_rate_ucl = purrr::map_dbl(ci, 2)
  ) %>%
  select(-ci)

death_by_age_plot <- hosp_death_by_frailty_age |> 
  mutate(age_group = factor(age_group, levels = c("Overall", "65 to 74", "75 to 84", "85 to 120"))) |> 
  mutate(frailty_category = factor(frailty_category, levels = c("Fit", "Mild", "Moderate", "Severe"))) |>
  filter(stringr::str_detect(cdm_name, "IQVIA", negate = T)) |> 
  ggplot(aes(x = age_group, y = death_pct, fill = frailty_category)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin = death_pct_lcl, ymax = death_pct_ucl), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~cdm_name, scales = "free_y") +
  coord_flip() +
  labs(
    x = "Age Group",
    y = "Percent",
    title = "One-year Mortality by Frailty Category by Age Group and Database",
    fill = "Frailty Category"
  ) +
  theme_minimal()

ggsave("code_for_publication/plots/death_by_frailty_age.jpg", death_by_age_plot, height = 8, width = 8)


hosp_death_by_frailty_age$frailty_category |> unique()

hosp_by_age_plot <- hosp_death_by_frailty_age |> 
  mutate(frailty_category = factor(frailty_category, levels = c("Fit", "Mild", "Moderate", "Severe"))) |> 
  mutate(age_group = factor(age_group, levels = c("Overall", "65 to 74", "75 to 84", "85 to 120"))) |> 
  filter(stringr::str_detect(cdm_name, "EBB|SIDIAP")) |>
  ggplot(aes(x = age_group, y = hosp_rate, fill = frailty_category)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin = hosp_rate_lcl, ymax = hosp_rate_ucl), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~cdm_name, scales = "free_y") +
  coord_flip() +
  labs(
    x = "Age Group",
    y = "Rate",
    title = "One-year Hospitalisation Rate by Frailty Category by Age Group and Database",
    fill = "Frailty Category"
  ) +
  theme_minimal()

ggsave("code_for_publication/plots/hosp_by_frailty_age_plot.jpg", hosp_by_age_plot, height = 8, width = 8)


# hospitalization rate by polypharmacy, age, cdm ----
hosp |> distinct(strata_name)

hosp_death_by_polypharm_age_overall <- hosp |> 
  filter(strata_name %in% c("polypharm_gte_5", "polypharm_gte_10")) |> 
  filter(group_level == "overall_ca_all_cancer_types_included_in_study") |> 
  mutate(strata_level = paste0(strata_name, '=', strata_level)) |> 
  select(cdm_name, strata_level, number_subjects_count, n_hosp = n_hospitalizations_in_year_after_dx_sum, death_flag_count) |> 
  mutate(age_group = "Overall", polypharmacy = case_when(
    strata_level == "polypharm_gte_5=0" ~ "<5",
    strata_level == "polypharm_gte_5=1" ~ ">=5",
    strata_level == "polypharm_gte_10=1" ~ ">=10",
    T ~ "drop"
  )) |> 
  filter(polypharmacy != "drop") |> 
  mutate(across(matches("number|n_|count"), as.numeric)) |> 
  filter(!is.na(number_subjects_count)) |> 
  mutate(cdm_name = case_when(
    stringr::str_detect(cdm_name, "gold") ~ "CPRD-GOLD",
    stringr::str_detect(cdm_name, "germany") ~ "IQVIA Germany",
    stringr::str_detect(cdm_name, "belgium") ~ "IQVIA Belgium",
    T ~ toupper(cdm_name)
  )) |> 
  select(cdm_name, age_group, polypharmacy, number_subjects_count, n_hosp, death_flag_count)

hosp_death_by_polypharm_age <- hosp |> 
  filter(strata_name %in% c("age_group and polypharm_gte_5", "age_group and polypharm_gte_10")) |> 
  filter(group_level == "overall_ca_all_cancer_types_included_in_study") |> 
  select(cdm_name, strata_name, strata_level, number_subjects_count, n_hosp = n_hospitalizations_in_year_after_dx_sum, death_flag_count) |> 
  tidyr::separate(col = "strata_level", into = c("age_group", "polypharm_level"), sep = " and ") |> 
  tidyr::separate(col = "strata_name", into = c("age_group_name", "polypharm_name"), sep = " and ") |> 
  mutate(polypharmacy = case_when(
    polypharm_name == "polypharm_gte_5" & polypharm_level == "0" ~ "<5",
    polypharm_name == "polypharm_gte_5" & polypharm_level == "1" ~ ">=5",
    polypharm_name == "polypharm_gte_10" & polypharm_level == "1" ~ ">=10",
    T ~ "drop"
  )) |> 
  filter(polypharmacy != "drop") |> 
  select(-polypharm_name, -polypharm_level, -age_group_name) |> 
  mutate(across(matches("number|n_|count"), as.numeric)) |> 
  filter(!is.na(number_subjects_count)) |> 
  filter(age_group != "None") |> 
  mutate(cdm_name = case_when(
    stringr::str_detect(cdm_name, "gold") ~ "CPRD-GOLD",
    stringr::str_detect(cdm_name, "germany") ~ "IQVIA Germany",
    stringr::str_detect(cdm_name, "belgium") ~ "IQVIA Belgium",
    T ~ toupper(cdm_name)
  )) |> 
  filter(stringr::str_detect(age_group, "18 to|45 to", negate = T)) |> 
  select(cdm_name, age_group, polypharmacy, number_subjects_count, n_hosp, death_flag_count) |> 
  bind_rows(hosp_death_by_polypharm_age_overall)

# add ci
hosp_death_by_polypharm_age <- hosp_death_by_polypharm_age |> 
  mutate(
    death_pct = 100 * death_flag_count / number_subjects_count,
    ci = purrr::map2(death_flag_count, number_subjects_count, function(x, n) {
      if (is.na(x)) return(c(NA_real_, NA_real_))
      res <- binom::binom.confint(x, n, conf.level = 0.95, methods = "wilson")
      c(res$lower * 100, res$upper * 100)
    }),
    death_pct_lcl = purrr::map_dbl(ci, 1),
    death_pct_ucl = purrr::map_dbl(ci, 2)
  ) %>%
  select(-ci) |> 
  mutate(
    hosp_rate = n_hosp / number_subjects_count,  
    ci = purrr::map2(n_hosp, number_subjects_count, poisson_ci),
    hosp_rate_lcl = purrr::map_dbl(ci, 1),
    hosp_rate_ucl = purrr::map_dbl(ci, 2)
  ) %>%
  select(-ci)

death_by_polypharm_age_plot <- hosp_death_by_polypharm_age |> 
  mutate(polypharmacy = factor(polypharmacy, levels = c("<5", ">=5", ">=10")),
         age_group = factor(age_group, levels = c("Overall", "65 to 74", "75 to 84", "85 to 120"))) |> 
  filter(stringr::str_detect(cdm_name, "IQVIA", negate = T)) |> 
  ggplot(aes(x = age_group, y = death_pct, fill = polypharmacy)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin = death_pct_lcl, ymax = death_pct_ucl), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~cdm_name, scales = "free_y") +
  coord_flip() +
  labs(
    x = "Age Group",
    y = "Percent",
    title = "One-year Mortality by Polypharmacy, Age Group, and Database",
    fill = "Polypharmacy"
  ) +
  theme_minimal()

ggsave("code_for_publication/plots/death_by_polypharm_age.jpg", death_by_polypharm_age_plot, height = 8, width = 8)

hosp_by_polypharm_age_plot <- hosp_death_by_polypharm_age |> 
  mutate(polypharmacy = factor(polypharmacy, levels = c("<5", ">=5", ">=10")),
         age_group = factor(age_group, levels = c("Overall", "65 to 74", "75 to 84", "85 to 120"))) |> 
  filter(stringr::str_detect(cdm_name, "EBB|SIDIAP")) |>
  ggplot(aes(x = age_group, y = hosp_rate, fill = polypharmacy)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin = hosp_rate_lcl, ymax = hosp_rate_ucl), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~cdm_name, scales = "free_y") +
  coord_flip() +
  labs(
    x = "Age Group",
    y = "Rate",
    title = "One-year Hospitalisation Rate by Polypharmacy, Age Group, and Database",
    fill = "Polypharmacy"
  ) +
  theme_minimal()

ggsave("code_for_publication/plots/hosp_by_polypharm_age_plot.jpg", hosp_by_polypharm_age_plot, height = 8, width = 8)


# polypharmacy by frailty category ----

polypharm_by_frailty <- demo |> 
  filter(strata_name %in% c("polypharm_gte_5 and frailty_category", "polypharm_gte_10 and frailty_category")) |>
  filter(variable_name == "Number subjects") |> 
  tidyr::separate(col = "strata_level", into = c("polypharm_level", "frailty_category"), sep = " and ") |> 
  tidyr::separate(col = "strata_name", into = c("polypharm_name", "frailty_name"), sep = " and ") |> 
  mutate(polypharmacy = case_when(
    polypharm_name == "polypharm_gte_5" & polypharm_level == "0" ~ "<5",
    polypharm_name == "polypharm_gte_5" & polypharm_level == "1" ~ ">=5",
    polypharm_name == "polypharm_gte_10" & polypharm_level == "1" ~ ">=10",
    T ~ "drop"
  )) |> 
  filter(polypharmacy != "drop") |> 
  filter(group_level == "overall_ca_all_cancer_types_included_in_study") |> 
  select(cdm_name, frailty_category, polypharmacy, n = estimate_value) |> 
  mutate(
    n = as.numeric(n),
    polypharmacy = factor(polypharmacy, levels = c("<5", ">=5", ">=10")),
    frailty_category = factor(stringr::str_to_title(frailty_category), levels = c("Fit", "Mild", "Moderate", "Severe"))
  ) |> 
  mutate(cdm_name = case_when(
    stringr::str_detect(cdm_name, "gold") ~ "CPRD-GOLD",
    stringr::str_detect(cdm_name, "germany") ~ "IQVIA Germany",
    stringr::str_detect(cdm_name, "belgium") ~ "IQVIA Belgium",
    T ~ toupper(cdm_name)
  )) |> 
  arrange(cdm_name, frailty_category, polypharmacy) |> 
  group_by(cdm_name, frailty_category) |> 
  mutate(pct = n/sum(n, na.rm = T)) |> 
  ungroup() |> 
  rename(cdm = cdm_name)

library(tidyr)
polypharm_by_frailty |> 
  tidyr::complete(.data$cdm, .data$frailty_category, .data$polypharmacy, fill = list(n = 0, pct = 0)) |> 
  mutate(polypharmacy = factor(polypharmacy, levels = c("<5", ">=5", ">=10"))) |> 
  ggplot(aes(x = frailty_category, y = pct, fill = polypharmacy)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~cdm, scales = "free_y") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Age Group",
    y = "Percent",
    title = "One-year Hospitalisation Rate by Polypharmacy, Age Group, and Database",
    fill = "Polypharmacy"
  ) +
  theme_minimal()

ggsave("code_for_publication/plots/hosp_by_polypharm_age_plot.jpg", hosp_by_polypharm_age_plot, height = 8, width = 8)



# figure solid tumors --

# need a dataset with frailty as strata with each cancer

prev |> distinct(strata_name)
prev2 <- prev |> 
  filter(strata_name %in% c("frailty_category")) |> 
  filter(group_level == "overall_ca_all_cancer_types_included_in_study") |> 
  collect()

condition_plot_df <- prev2 |> 
  filter(stringr::str_detect(variable_name, "frailty")) |> 
  filter(variable_name != "frailty_category") |> 
  filter(estimate_name == "percentage") |> 
  transmute(cdm_name, frailty_category = strata_level, variable_name, percent = as.numeric(estimate_value)) |> 
  mutate(variable_name = stringr::str_replace_all(variable_name, "frailty_|_jp", "")) |> 
  mutate(variable_name = stringr::str_replace_all(variable_name, "_", " ")) |> 
  mutate(variable_name = stringr::str_to_title(variable_name)) |> 
  mutate(variable_name = case_when(
    variable_name == "Memeroy Cogn Problem" ~ "Memory-Cognition Problem",
    T~variable_name
  )) |> 
  mutate(variable_name = stringr::str_replace_all(variable_name, "Dx", "DX")) |>
  mutate(variable_name = stringr::str_replace_all(variable_name, "C Kidney", "Chronic Kidney")) |>
  # filter(stringr::str_detect(cdm_name, "ebb|iqvia")) |> 
  mutate(cdm_name = case_when(
    stringr::str_detect(cdm_name, "gold") ~ "CPRD-GOLD",
    stringr::str_detect(cdm_name, "germany") ~ "IQVIA Germany",
    stringr::str_detect(cdm_name, "belgium") ~ "IQVIA Belgium",
    T ~ toupper(cdm_name)
  )) 
  

lvl <- condition_plot_df |> 
  group_by(variable_name) |> 
  summarise(m = mean(percent, na.rm = T)) |> 
  arrange(m) |> 
  pull(variable_name)

distinct(condition_plot_df, cdm_name)

condition_plot1 <- condition_plot_df |> 
  filter(stringr::str_detect(tolower(cdm_name), "ebb|iqvia")) |>
  mutate(percent = percent/100,
         frailty_category = stringr::str_to_title(frailty_category)) |> 
  mutate(variable_name = factor(variable_name, levels = lvl)) |> 
  ggplot(aes(x = percent, y = variable_name)) +
  geom_point() + 
  scale_x_continuous(labels = scales::percent) +
  labs(
    x = "Percent",
    y = "",
    title = ""
  ) +
  facet_grid(cdm_name ~ frailty_category) +
  theme_minimal()

condition_plot2 <- condition_plot_df |> 
  filter(!stringr::str_detect(tolower(cdm_name), "ebb|iqvia")) |>
  mutate(percent = percent/100,
         frailty_category = stringr::str_to_title(frailty_category)) |> 
  mutate(variable_name = factor(variable_name, levels = lvl)) |> 
  ggplot(aes(x = percent, y = variable_name)) +
  geom_point() + 
  scale_x_continuous(labels = scales::percent) +
  labs(
    x = "Percent",
    y = "",
    title = ""
  ) +
  facet_grid(cdm_name ~ frailty_category) +
  theme_minimal()

# gridExtra::grid.arrange(condition_plot1, condition_plot2)

library(patchwork)
condition_plot <- (condition_plot1 + condition_plot2) +
 plot_annotation(
  title = 'Condition Prevalence by Frailty Category'
 )

ggsave("code_for_publication/plots/condition_prevalence_by_frailty.jpg", condition_plot, height = 20, width = 20)




# demo |> 
#   distinct(strata_name)
# 
# polypharm_gte_5_by_age <- demo |>  
#   filter(strata_name == "age_group and polypharm_gte_5") |> 
#   filter(variable_name == "Number subjects") |> 
#   filter(group_level == "overall_ca_all_cancer_types_included_in_study") |> 
#   # distinct(group_level)
#   select(cdm_name, strata_level, count = estimate_value) |> 
#   tidyr::separate(col = "strata_level", into = c("age_group", "polypharm_gte_5"), sep = " and ") |> 
#   mutate(count = as.numeric(count)) |> 
#   filter(age_group != "None") |> 
#   mutate(cdm_name = case_when(
#     stringr::str_detect(cdm_name, "gold") ~ "CPRD-GOLD",
#     stringr::str_detect(cdm_name, "germany") ~ "IQVIA Germany",
#     stringr::str_detect(cdm_name, "belgium") ~ "IQVIA Belgium",
#     T ~ toupper(cdm_name)
#   )) |> 
#   filter(stringr::str_detect(age_group, "18 to|45 to", negate = T)) |> 
#   group_by(cdm_name, age_group) |> 
#   mutate(percent = count/sum(count, na.rm = T)) |> 
#   mutate(category = "polypharmacy >= 5") |> 
#   filter(polypharm_gte_5 == "1") |> 
#   select(cdm_name, age_group, category, count, percent)
# 
# 
# polypharm_gte_10_by_age <- demo |> 
#   filter(strata_name == "age_group and polypharm_gte_10") |> 
#   filter(variable_name == "Number subjects") |> 
#   filter(group_level == "overall_ca_all_cancer_types_included_in_study") |> 
#   # distinct(group_level)
#   select(cdm_name, strata_level, count = estimate_value) |> 
#   tidyr::separate(col = "strata_level", into = c("age_group", "polypharm_gte_10"), sep = " and ") |> 
#   mutate(count = as.numeric(count)) |> 
#   filter(age_group != "None") |> 
#   mutate(cdm_name = case_when(
#     stringr::str_detect(cdm_name, "gold") ~ "CPRD-GOLD",
#     stringr::str_detect(cdm_name, "germany") ~ "IQVIA Germany",
#     stringr::str_detect(cdm_name, "belgium") ~ "IQVIA Belgium",
#     T ~ toupper(cdm_name)
#   )) |> 
#   filter(stringr::str_detect(age_group, "18 to|45 to", negate = T)) |> 
#   group_by(cdm_name, age_group) |> 
#   mutate(percent = count/sum(count, na.rm = T)) |> 
#   mutate(category = "polypharmacy >= 10") |> 
#   filter(polypharm_gte_10 == "1") |> 
#   select(cdm_name, age_group, category, count, percent)
# 
# 
# polypharm_by_age <- bind_rows(polypharm_gte_5_by_age, polypharm_gte_10_by_age)
# 
# library(ggplot2)
# polypharm_by_age |> 
#   ggplot(aes(x = age_group, y = percent, fill = category)) +
#   geom_col(position = position_dodge()) +
#   facet_wrap(~cdm_name, scales = "free_y") +
#   coord_flip() +
#   labs(
#     x = "Age Group",
#     y = "Percent", 
#     title = "Polypharmacy by Age Groups & Database",
#     fill = "Polypharmcy Category"
#   ) +
#   scale_y_continuous(labels = scales::percent) +
#   theme_minimal() +
#   theme(strip.text.y = element_text(angle = 0),
#         axis.text.y.right = element_blank(),
#         axis.ticks.y.right = element_blank())
# 
# polypharm_by_age |> 
#   ggplot(aes(x = age_group, y = percent, fill = category)) +
#   geom_col(position = position_dodge()) +
#   facet_wrap(~cdm_name, scales = "free_y") +
#   coord_flip() +
#   labs(
#     x = "Age Group",
#     y = "Percent", 
#     title = "Polypharmacy by Age Groups & Database",
#     fill = "Polypharmacy Category"
#   ) +
#   scale_y_continuous(labels = scales::percent) +
#   theme_minimal() +
#   theme(
#     axis.text.y.left = element_blank(),  # remove axis labels
#   )
# 
# 

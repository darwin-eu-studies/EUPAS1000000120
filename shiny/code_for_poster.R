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

# code for poster

con <- DBI::dbConnect(duckdb::duckdb(), here::here("results", "results_combined.duckdb"))

DBI::dbListTables(con)

library(dplyr)
demo <- tbl(con, "demographics") |> collect()
prev <- tbl(con, "prevalence") 
hosp <- tbl(con, "hospitalization_death") |> collect()


poisson_ci <- function(events, people, conf.level = 0.95) {
  if (is.na(events)) return(c(NA_real_, NA_real_))
  ci <- poisson.test(events, T = people, conf.level = conf.level)$conf.int
  return(c(lower = ci[1], upper = ci[2]))
}

hosp_death_by_frailty <- hosp |> 
  filter(strata_name == "frailty_category") |> 
  filter(group_level == "overall_ca_all_cancer_types_included_in_study") |> 
  select(cdm_name, strata_level, number_subjects_count, n_hosp = n_hospitalizations_in_year_after_dx_sum, death_flag_count) |> 
  rename(frailty_category  = strata_level) |> 
  mutate(across(matches("number|n_|count"), as.numeric)) |> 
  filter(!is.na(number_subjects_count)) |> 
  mutate(frailty_category = stringr::str_to_title(frailty_category)) |> 
  mutate(cdm_name = case_when(
    stringr::str_detect(cdm_name, "gold") ~ "CPRD-GOLD",
    stringr::str_detect(cdm_name, "germany") ~ "IQVIA Germany",
    stringr::str_detect(cdm_name, "belgium") ~ "IQVIA Belgium",
    T ~ toupper(cdm_name)
  )) |> 
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


plot_data <- hosp_death_by_frailty |> 
  select(cdm_name, frailty_category, number_subjects_count, starts_with("hosp"))

hosp_by_frailty_plot_data <- hosp_death_by_frailty |> 
  select(cdm_name, frailty_category, number_subjects_count, matches("hosp")) |> 
  mutate(frailty_category = factor(frailty_category, levels = c("Fit", "Mild", "Moderate", "Severe"))) |> 
  filter(stringr::str_detect(cdm_name, "EBB|SIDIAP"))

# here is the data for the plot
# datapasta::tribble_paste(hosp_by_frailty_plot_data)
# hosp_by_frailty_plot_data <- tibble::tribble(
#  ~cdm_name, ~frailty_category, ~number_subjects_count, ~n_hosp,       ~hosp_rate,   ~hosp_rate_lcl,   ~hosp_rate_ucl,
#   "SIDIAP",             "Fit",                  49667,   58627, 1.18040147381561, 1.17086559540574, 1.18999570789186,
#   "SIDIAP",            "Mild",                  25529,   35409, 1.38701085040542, 1.37260120521151, 1.40153407270355,
#   "SIDIAP",        "Moderate",                   9998,   15181, 1.51840368073615, 1.49434477176074, 1.54275287281818,
#   "SIDIAP",          "Severe",                   3529,    5443,  1.5423632757155, 1.50165764117835, 1.58389281604653,
#      "EBB",             "Fit",                   1171,    3903,   3.333048676345, 3.22929385868331, 3.43928851242882,
#      "EBB",            "Mild",                   1021,    3478, 3.40646425073457, 3.29418414720899, 3.52159536598681,
#      "EBB",        "Moderate",                    448,    1478, 3.29910714285714, 3.13303844483964, 3.47169303191947,
#      "EBB",          "Severe",                    129,     479, 3.71317829457364, 3.38805393983496, 4.06108391270087
#  )

hosp_by_frailty_plot <- hosp_by_frailty_plot_data |> 
  ggplot(aes(x = hosp_rate, y = frailty_category, fill = frailty_category)) +
  geom_col(width = 0.9) +
  geom_errorbarh(aes(xmin = hosp_rate_lcl, xmax = hosp_rate_ucl), height = 0.25, linewidth = 0.6) +
  facet_wrap(~ cdm_name, nrow = 1) +
  scale_x_continuous(limits = c(0, 5), breaks = 0:5) +
  labs(x = "One-Year Hospitalisation Rate", y = NULL, fill = "Frailty") +
  theme_minimal(base_size = 13) +
  scale_fill_manual(
    name = "Frailty Category",
    values = c(
      Fit      = "#FF6F6F", # red
      Mild     = "#B38A1E", # brownish-gold
      Moderate = "#00A850", # green
      Severe   = "#C77CFF"  # purple
    )
  ) +
  guides(fill = guide_legend(reverse = TRUE)) + 
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 14),
    axis.text.y = element_blank()
  )

ggsave("hosp_by_frailty_plot_for_poster.jpg", hosp_by_frailty_plot, width = 12, height = 3, dpi = 1200)
 

mort_by_frailty_plot_data <- hosp_death_by_frailty |> 
  select(cdm_name, frailty_category, number_subjects_count, matches("death")) |> 
  mutate(frailty_category = factor(frailty_category, levels = c("Fit", "Mild", "Moderate", "Severe"))) |> 
  filter(stringr::str_detect(cdm_name, "EBB|SIDIAP|CPRD-GOLD|IPCI")) |> 
  mutate(cdm_name = factor(cdm_name, levels = c("EBB", "SIDIAP","IPCI", "CPRD-GOLD")))

# Data for the mortality plot
# datapasta::tribble_paste(mort_by_frailty_plot_data)
# tibble::tribble(
#       ~cdm_name, ~frailty_category, ~number_subjects_count, ~death_flag_count,       ~death_pct,    ~death_pct_lcl,    ~death_pct_ucl,
#          "IPCI",             "Fit",                  15688,              3966, 25.2804691483937,  24.6064751524279,  25.9665661292236,
#          "IPCI",            "Mild",                   4115,              1675, 40.7047387606318,  39.2130314337522,  42.2137846345935,
#          "IPCI",        "Moderate",                    855,               431, 50.4093567251462,  47.0636732989355,  53.7513781786226,
#          "IPCI",          "Severe",                     97,                56, 57.7319587628866,  47.7914606447473,  67.0833737491813,
#     "CPRD-GOLD",             "Fit",                  32416,              8014, 24.7223593287266,  24.2557530414239,  25.1949559609139,
#     "CPRD-GOLD",            "Mild",                  11012,              4635, 42.0904467853251,   41.171253523853,  43.0151565065389,
#     "CPRD-GOLD",        "Moderate",                   2753,              1421, 51.6164184525972,  49.7487081920148,  53.4796239891934,
#     "CPRD-GOLD",          "Severe",                    527,               299, 56.7362428842505,  52.4725972199635,  60.9023942808383,
#        "SIDIAP",             "Fit",                  49667,             10637,  21.416634787686,  21.0580620892555,  21.7796286643046,
#        "SIDIAP",            "Mild",                  25529,              9455, 37.0363116455795,  36.4459373906626,  37.6305867178072,
#        "SIDIAP",        "Moderate",                   9998,              5615, 56.1612322464493,   55.186439165168,  57.1312925752555,
#        "SIDIAP",          "Severe",                   3529,              2546, 72.1450835930859,  70.6425836467076,  73.5994243295503,
#           "EBB",             "Fit",                   1171,               157, 13.4073441502989,  11.5749580308059,  15.4790292484431,
#           "EBB",            "Mild",                   1021,               193, 18.9030362389814,  16.6196568063487,  21.4195399362668,
#           "EBB",        "Moderate",                    448,               116, 25.8928571428571,  22.0536028239038,  30.1420189404735,
#           "EBB",          "Severe",                    129,                45, 34.8837209302326,   27.204358425833,  43.4373370612354
# )


mort_by_frailty_plot <- mort_by_frailty_plot_data |> 
  ggplot(aes(x = death_pct, y = frailty_category, fill = frailty_category)) +
  geom_col(width = 0.9) +
  geom_errorbarh(aes(xmin = death_pct_lcl, xmax = death_pct_ucl), height = 0.25, linewidth = 0.6) +
  facet_wrap(~ cdm_name, nrow = 1) +
  # scale_x_continuous(limits = c(0, 5), breaks = 0:5) +
  labs(x = "One-Year Mortality Risk", y = NULL, fill = "Frailty") +
  theme_minimal(base_size = 13) +
  scale_fill_manual(
    name = "Frailty Category",
    values = c(
      Fit      = "#FF6F6F", # red
      Mild     = "#B38A1E", # brownish-gold
      Moderate = "#00A850", # green
      Severe   = "#C77CFF"  # purple
    )
  ) +
  guides(fill = guide_legend(reverse = TRUE)) + 
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 14),
    axis.text.y = element_blank()
  )

ggsave("mortality_by_frailty_plot_for_poster.jpg", mort_by_frailty_plot, width = 12, height = 3, dpi = 1200)



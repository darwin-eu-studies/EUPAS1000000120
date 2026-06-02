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

library(ggplot2)
library(dplyr)

source("global.R")

solid <- c("breast_cancer", "colorectal_cancer", "endometrial_cancer", 
           "ovarian_cancer", "pancreatic_cancer", "prostate_cancer", 
           "lung_cancer", "sclc", "nsclc")

blood <- c("lymphoma", "non_hk_lymphoma", "hk_lymphoma", "leukemia", "aml", "cml", "lla", "llc", "multiple_myeloma")

cancer_levels <- c(solid, blood)

# update for publiction: 
# Please remove 
# SCLC
# NSCLC
# non-HK lymphoma
# HK-lymphoma
# AML
# CML
# LLA
# LLC

solid <- c("breast_cancer", "colorectal_cancer", "endometrial_cancer", 
           "ovarian_cancer", "pancreatic_cancer", "prostate_cancer", 
           "lung_cancer")

blood <- c("lymphoma", "leukemia", "multiple_myeloma")

cancer_levels <- c(solid, blood)


prev %>% distinct(group_level)

prev %>% distinct(cdm_name)

plt1 <- prev %>% 
  filter(strata_name == "overall", variable_name %in% c("frailty_category")) %>%
  collect() %>% 
  mutate(cdm_name = replace_values(.data$cdm_name,
    "sidiap"            ~ "SIDIAP",
    "cdm_gold_202307"   ~ "CPRD GOLD",
    "iqvia_germany_da"  ~ "IQVIA DA Germany",
    "iqvia_belgium_lpd" ~ "IQVIA PD Belgium",
    "ipci"              ~ "IPCI",
    "ebb"  ~ "EBB")) %>% 
  # filter(cdm_name == "ipci") %>%
  # distinct(variable_name)
  select(cdm_name, cohort = group_level, variable_name, variable_level, estimate_name, estimate_value) %>% 
  tidyr::pivot_wider(names_from = "estimate_name", values_from = "estimate_value") %>% 
  filter(cohort %in% cancer_levels) %>% 
  mutate(percentage = as.numeric(percentage),
         `Frailty category` = factor(variable_level, levels = rev(c("fit", "mild", "moderate", "severe"))),
         cohort = factor(cohort, levels = rev(cancer_levels), labels = stringr::str_to_title(stringr::str_replace_all(rev(cancer_levels), "_", " ")))) %>% 
  # group_by(cohort, cdm_name) %>% 
  # summarise(pct = sum(percentage))
  ggplot(aes(y = percentage, x = cohort, fill = `Frailty category`)) + 
  geom_bar(position="stack", stat="identity") +
  # scale_fill_brewer(direction = 1) +
  facet_wrap("cdm_name") +
  coord_flip() +
  theme_bw()

# ggsave("frailty_barplot.png", plot = plt1, height = 10, width = 10)
ggsave("frailty_barplot_updated.png", plot = plt1, height = 10, width = 10)


prev %>% 
  filter(strata_name == "overall", variable_name %in% c("polypharm_gte_5", "polypharm_gte_10")) %>% 
  View()

plt2 <- prev %>% 
  filter(strata_name %in% c("polypharm_gte_5", "polypharm_gte_10"), variable_name == "number subjects") %>%
  collect() %>% 
  mutate(cdm_name = replace_values(.data$cdm_name,
                                   "sidiap"            ~ "SIDIAP",
                                   "cdm_gold_202307"   ~ "CPRD GOLD",
                                   "iqvia_germany_da"  ~ "IQVIA DA Germany",
                                   "iqvia_belgium_lpd" ~ "IQVIA PD Belgium",
                                   "ipci"              ~ "IPCI",
                                   "ebb"  ~ "EBB")) %>% 
  tidyr::unite(col = "strata", strata_name, strata_level) %>% 
  select(cdm_name, cohort = group_level, strata, estimate_value) %>% 
  filter(cohort %in% cancer_levels) %>% 
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
  mutate(cohort = factor(cohort, levels = rev(cancer_levels), labels = stringr::str_to_title(stringr::str_replace_all(rev(cancer_levels), "_", " ")))) %>% 
  ggplot(aes(x = cohort, y = Percent, fill = Polypharmacy)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  facet_wrap("cdm_name") +
  theme_bw()
  
ggsave("polypharm_barplot.png", plot = plt2, height = 10, width = 10)








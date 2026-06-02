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




solid <- c("breast_cancer", "colorectal_cancer", "endometrial_cancer", 
           "ovarian_cancer", "pancreatic_cancer", "prostate_cancer", 
           "lung_cancer")

blood <- c("lymphoma", "leukemia", "multiple_myeloma")

cancer_levels <- c(solid, blood)

databases <- c("ebb", "sidiap", "ipci", "cdm_gold_202307")

death_flag_percentage

# hospitalization
plt1 <- hosp %>% 
  filter(strata_name == "frailty_category") %>% 
  select(cdm_name, cohort = group_level, strata_name, strata_level, percent_hospitalized = hospitalization_rate) %>% 
  filter(cdm_name %in% c("sidiap", "ebb")) %>%
  mutate(cdm_name = factor(cdm_name, levels = databases)) %>% 
  filter(cohort %in% cancer_levels) %>% 
  mutate(frailty_category = factor(strata_level, levels = rev(c("fit", "mild", "moderate", "severe"))),
         cohort = factor(cohort, levels = rev(cancer_levels))) %>% 
  ggplot(aes(y = percent_hospitalized, x = cohort, fill = cdm_name)) + 
  geom_bar(position="dodge", stat="identity") +
  facet_wrap("frailty_category") +
  coord_flip() +
  theme_bw() +
  labs(y = "Hospitalization rate within one year of index", x = "", title = "Hospitalization rate by frailty category")

ggsave("hosp_by_frailty.png", plt1, height = 10, width = 10)

plt2 <- hosp %>% 
  filter(strata_name %in% c("polypharm_gte_10", "polypharm_gte_5")) %>% 
  select(cdm_name, cohort = group_level, strata_name, strata_level, percent_hospitalized = hospitalization_rate) %>% 
  filter(cdm_name %in% c("sidiap", "ebb"), strata_level == "1") %>%
  filter(cohort %in% cancer_levels) %>% 
  mutate(cdm_name = factor(cdm_name, levels = databases)) %>% 
  mutate(cohort = factor(cohort, levels = rev(cancer_levels))) %>% 
  ggplot(aes(y = percent_hospitalized, x = cohort, fill = cdm_name)) + 
  geom_bar(position="dodge", stat="identity") +
  facet_wrap("strata_name") +
  coord_flip() +
  theme_bw() +
  labs(y = "Hospitalization rate within one year of index", x = "", title = "Hospitalization rate by polypharmacy category")


ggsave("hosp_by_polypharm.png", plt2, height = 10, width = 10)


# mortality
plt3 <- hosp %>%
  filter(strata_name == "frailty_category") %>% 
  select(cdm_name, cohort = group_level, strata_name, strata_level, percent_died = death_flag_percentage) %>% 
  mutate(percent_died = as.numeric(percent_died)) %>% 
  filter(cdm_name %in% c("sidiap", "ebb", "ipci", "cdm_gold_202307")) %>%
  mutate(cdm_name = factor(cdm_name, levels = databases)) %>% 
  filter(cohort %in% cancer_levels) %>%
  mutate(frailty_category = factor(strata_level, levels = c("fit", "mild", "moderate", "severe"))) %>% 
  mutate(cohort = factor(cohort, levels = rev(cancer_levels))) %>%
  ggplot(aes(y = percent_died, x = cohort, fill = cdm_name)) + 
  geom_bar(position="dodge", stat="identity") +
  facet_wrap("frailty_category") +
  coord_flip() +
  theme_bw() +
  labs(y = "Mortality risk within one year of cancer diagnosis", x = "", title = "Mortality rate by frailty category")

ggsave("mortality_by_frailty.png", plt3, height = 10, width = 10)


plt4 <- hosp %>%
  filter(strata_name %in% c("polypharm_gte_10", "polypharm_gte_5")) %>% 
  select(cdm_name, cohort = group_level, strata_name, strata_level, percent_died = death_flag_percentage) %>% 
  mutate(percent_died = as.numeric(percent_died)) %>% 
  filter(cdm_name %in% c("sidiap", "ebb", "ipci", "cdm_gold_202307"), strata_level == "1") %>%
  mutate(cdm_name = factor(cdm_name, levels = databases)) %>% 
  filter(cohort %in% cancer_levels) %>%
  mutate(frailty_category = factor(strata_level, levels = c("fit", "mild", "moderate", "severe"))) %>% 
  mutate(cohort = factor(cohort, levels = rev(cancer_levels))) %>%
  ggplot(aes(y = percent_died, x = cohort, fill = cdm_name)) + 
  geom_bar(position="dodge", stat="identity") +
  facet_wrap("strata_name") +
  coord_flip() +
  theme_bw() +
  labs(y = "Mortality risk within one year of cancer diagnosis", x = "", title = "Mortality rate by polypharmacy category")

ggsave("mortality_by_polypharm.png", plt4, height = 10, width = 10)




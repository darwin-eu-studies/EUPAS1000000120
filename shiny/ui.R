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


library(shiny)
library(shinydashboard)
library(dplyr)
library(dbplyr)

databases <- demo %>% 
  distinct(.data$cdm_name) %>% 
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

cohort_names <- unique(cohorts$cohort_name)

polypharm <- c("All persons" = "NA", "0-4 drugs", ">=5 drugs", ">=10 drugs")
pre_2020 <- c("All time" = "NA", "Pre 2020" = "TRUE", "2020 and later" = "FALSE")

time_windows <- c("-999999 to -366", "-365 to -31", "-30 to -1", "0 to 0")

available_plots <- c(
  "1. Overall prevalence of frailty conditions by frailty category among individuals with selected cancer types" = 1,
  "2. Prevalence of frailty categories by cancer type" = 2,
  "3. Prevalence of polypharmacy by cancer type" = 3,
  "4. Prevalence of frailty conditions by frailty category and Solid cancer type" = 4,
  "5. Prevalence of frailty conditions by frailty category and Blood cancer type" = 5,
  "6. Hospitalization rate by frailty category" = 6,
  "7. Hospitalization rate by polypharmacy category" = 7,
  "8. Mortality risk by frailty category" = 8,
  "9. Mortality risk by polypharmacy category" = 9
)

ui <- dashboardPage(
  dashboardHeader(title = "Study P2-C1-009"),
  ## menu ----
  dashboardSidebar(
    sidebarMenu(
      menuItem(text = "Background", tabName = "background"),
      menuItem(text = "Databases", tabName = "dbs",
        menuSubItem(text = "Database details", tabName = "cdm_snapshot")
      ),
      menuItem(text = "Cohorts", tabName = "cohorts",
        menuSubItem(text = "Cohort attrition", tabName = "cohort_attrition")
      ),
      menuItem(text = "Polypharmacy and Frailty", tabName = "tables1_2",
        menuSubItem(text = "Obj 1 - Prevalence", tabName = "table1"),
        menuSubItem(text = "Obj 2 - Population characteristics", tabName = "table2")
      ),
    
      menuItem(text = "Hospitalization and Death", tabName = "table3"),
      
      menuItem(text = "Large Scale Characterization", tabName = "ls_characterization",
        menuSubItem("Conditions", tabName = "ls_conditions"),
        menuSubItem("Drugs", tabName = "ls_drugs")     
      ),
      menuItem(text = "Plots", tabName = "plots")
    )
  ),
  
  ## body ----
  dashboardBody(
    tabItems(
      tabItem(tabName = "background", includeMarkdown("abstract.md")),
      tabItem(
        tabName = "cdm_snapshot",
        htmlOutput('tbl_cdm_snaphot'),
        tags$hr(),
        div(style="display:inline-block",
            downloadButton(
              outputId = "gt_cdm_snaphot_word",
              label = "Download table as word"
            ), 
            style="display:inline-block; float:right")
      ),
      
      tabItem(tabName = "cohort_attrition",
              DT::dataTableOutput("cohort_attrition"),
              downloadButton("cohort_attrition_download", "Download csv")),
      
      tabItem(tabName = "table1",
        pickerInput("table1_cancer","Cancer", cohorts$cohort_name, inline = T),
        pickerInput("table1_age","Age Group", age_groups, inline = T, multiple = T, selected = "NA"),
        pickerInput("table1_sex", "Sex", sex, inline = T, multiple = T, selected = c("Male", "Female")),
        pickerInput("table1_frailty", "Frailty", frailty, inline = T),
        pickerInput("table1_polypharm", "polypharm", polypharm, inline = T),
        pickerInput("table1_pre_2020", "pre_2020", pre_2020, inline = T),
        tags$p("Note: NA implies data was supressed due to low cell counts."),
        DT::dataTableOutput("table1"),
        textOutput("table1_no_data_message"),
        downloadButton("table1_download", "Download csv")
      ),
      
      tabItem(tabName = "table2",
              pickerInput("table2_cdm_name","Database", databases, inline = T),
              pickerInput("table2_cancer","Cancer", cohorts$cohort_name, inline = T),
              pickerInput("table2_age","Age Group", age_groups, inline = T),
              pickerInput("table2_sex", "Sex", sex, inline = T),
              pickerInput("table2_pre_2020", "pre_2020", pre_2020, inline = T),
              DT::dataTableOutput("table2"),
              textOutput("table2_no_data_message"),
              downloadButton("table2_download", "Download csv")
      ),
      
      tabItem(tabName = "table3",
              pickerInput("table3_cancer","Cancer", cohorts$cohort_name, inline = T),
              pickerInput("table3_age","Age Group", age_groups, inline = T),
              pickerInput("table3_sex", "Sex", sex, inline = T),
              pickerInput("table3_pre_2020", "pre_2020", c("All time" = "NA", "Pre 2020" = "true", "2020 and later" = "false"), inline = T),
              DT::dataTableOutput("table3"),
              tags$p("Databases with zero counts or no clear denominator have been removed from the results."),
              tags$p("95% confidence interval is shown in parentheses"),
              textOutput("table3_no_data_message"),
              downloadButton("table3_download", "Download csv"),
              downloadButton("table3_raw_download", "Download full results as csv")
      ),
      
      tabItem(tabName = "ls_drugs",
              pickerInput("ls_drugs_cdm_name", "Database", databases, inline = T, multiple = T, selected = databases[1]),
              # pickerInput("ls_drugs_time_window", "Timeframe", time_windows, inline = T, multiple = F, selected = time_windows[1]),
              # pickerInput("ls_drugs_cohort", "Cohort", cohort_names, inline = T, multiple = F, selected = cohort_names[1]),
              pickerInput("ls_drugs_age", "Age Group", age_groups, inline = T, multiple = F, selected = age_groups[1]),
              pickerInput("ls_drugs_sex", "Sex", sex, inline = T, multiple = F, selected = sex[1]),
              pickerInput("ls_drugs_frailty", "Frailty", frailty, inline = T, multiple = F, selected = frailty[1]),
              pickerInput("ls_drugs_polypharm", "polypharm", polypharm, inline = T, multiple = F, selected = polypharm[1]),
              pickerInput("ls_drugs_pre_2020", "pre_2020", pre_2020, inline = T, multiple = F, selected = pre_2020[1]),
              DT::dataTableOutput("ls_drugs"),
              textOutput("ls_drugs_no_data_message"),
              downloadButton("ls_drugs_download", "Download csv"),
              tags$h3("Top 20 drugs by database"),
              DT::dataTableOutput("ls_drugs_top"),
              textOutput("ls_drugs_top_no_data_message"),
              downloadButton("ls_drugs_download_top_n", "Download csv")
      ),
      
      tabItem(tabName = "ls_conditions",
              pickerInput("ls_conditions_cdm_name", "Database", databases, inline = T, multiple = T, selected = databases[1]),
              # pickerInput("ls_conditions_time_window", "Timeframe", time_windows, inline = T, multiple = F, selected = time_windows[1]),
              # pickerInput("ls_conditions_cohort", "Cohort", cohort_names, inline = T, multiple = F, selected = cohort_names[1]),
              pickerInput("ls_conditions_age", "Age Group", age_groups, inline = T, multiple = F, selected = age_groups[1]),
              pickerInput("ls_conditions_sex", "Sex", sex, inline = T, multiple = F, selected = sex[1]),
              pickerInput("ls_conditions_frailty", "Frailty", frailty, inline = T, multiple = F, selected = frailty[1]),
              pickerInput("ls_conditions_polypharm", "polypharm", polypharm, inline = T, multiple = F, selected = polypharm[1]),
              pickerInput("ls_conditions_pre_2020", "pre_2020", pre_2020, inline = T, multiple = F, selected = pre_2020[1]),
              DT::dataTableOutput("ls_conditions"),
              textOutput("ls_conditions_no_data_message"),
              downloadButton("ls_conditions_download", "Download csv"),
              tags$h3("Top 20 conditions by database"),
              DT::dataTableOutput("ls_conditions_top"),
              textOutput("ls_conditions_top_no_data_message"),
              downloadButton("ls_conditions_download_top_n", "Download csv")
              
      ),
      tabItem(tabName = "plots", 
              tags$h3("Plots for the study report"), tags$br(),
              pickerInput("plot_age_group", "Age Group for plot", age_groups, inline = T, multiple = F, selected = age_groups[1]), tags$br(),
              pickerInput("plot_selection", "Select plot", available_plots, inline = F, multiple = F, selected = 1, width = "100%"), tags$br(),
              downloadButton("download_plot", "Download plot"), tags$br(),
              plotOutput("plot_to_display")
      )
    )
  )
)

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

# server shiny ----
server <- function(input, output, session) {
  
  # session$onSessionEnded(function() pool::dbDisconnect(con, shutdown = T))
  
  # CDM Snapshot -----
  output$tbl_cdm_snaphot <- renderText({
    kableExtra::kable(cdm_snapshot) %>% 
      kableExtra::kable_styling("striped", full_width = F)
  })
  
  output$gt_cdm_snaphot_word <- downloadHandler(
    filename = function() {"cdm_snapshot.docx"},
    content = function(file) {gt::gtsave(gt::gt(cdm_snapshot), file)}
  )
  
  # cohort_attrition ----
  output$cohort_attrition <- DT::renderDataTable({
    cohort_attrition %>% 
      select(-cohort_definition_id) %>% 
      mutate_at(c("cdm_name", "cohort_name", "reason"), as.factor) %>% 
      DT::datatable(
        rownames = FALSE,
        selection = "single",
        filter = "top")
  })
  
  output$cohort_attrition_download <- downloadHandler(
    filename = function() { "cohort_attrition.csv" },
    content = function(file) { readr::write_csv(cohort_attrition, file) }
  )
  
  # Table 1 -----
  table1 <- reactive({
    table1_filtered3(prev,
           cancer = input$table1_cancer,
           polypharm = input$table1_polypharm,
           age = input$table1_age,
           sex = input$table1_sex,
           frailty = input$table1_frailty,
           pre_2020 = input$table1_pre_2020)
  })
  
  output$table1 <- DT::renderDataTable(table1(), rownames = FALSE, options = list(
    pageLength = 100, scrollX = TRUE, ordering=F))
  output$table1_no_data_message <- renderText(ifelse(is.null(table1()), "There is no data for the current filter selection", ""))
  output$table1_download <- downloadHandler(
    filename = function() { "table1.csv" },
    content = function(file) { readr::write_csv(table1(), file) }
  )
  
  # Table 2 -----
  table2 <- reactive({
    table2_filtered(prev,
                    cdm_name = input$table2_cdm_name,
                    cancer = input$table2_cancer,
                    age = input$table2_age,
                    sex = input$table2_sex,
                    pre_2020 = input$table2_pre_2020) 
  })
  
  output$table2 <- DT::renderDataTable(table2(), rownames = FALSE, options = list(pageLength = 100))
  output$table2_no_data_message <- renderText(ifelse(is.null(table2()), "There is no data for the current filter selection", ""))
  output$table2_download <- downloadHandler(
    filename = function() { "table2.csv" },
    content = function(file) { readr::write_csv(table2(), file) }
  )
  
  # Table 3 -----
  
  # observe({
  #   print(paste("table3_cdm_name:", input$table3_cdm_name))
  #   print(paste("table3_age:", input$table3_age))
  #   print(paste("table3_sex:", input$table3_sex))
  #   print(paste("table3_pre_2020:", input$table3_pre_2020))
  # })
  
  table3 <- reactive({
    table3_filtered(
      hosp,
      cancer = input$table3_cancer,
      age = input$table3_age,
      sex = input$table3_sex,
      pre_2020 = input$table3_pre_2020
    )
  })
  
  output$table3 <- DT::renderDataTable(table3(), rownames = FALSE, options = list(pageLength = 100))
  output$table3_no_data_message <- renderText(ifelse(is.null(table3()), "There is no data for the current filter selection", ""))
  output$table3_download <- downloadHandler(
    filename = function() { "table3.csv" },
    content = function(file) { readr::write_csv(table3(), file) }
  )
  output$table3_raw_download <- downloadHandler(
    filename = function() { "hospitalization_death.csv" },
    content = function(file) { readr::write_csv(collect(hosp), file) }
  )
  
  # Large scale drugs ------
  ls_drugs <- reactive({
    
    # polypharmacy stratification options.
    # Note that the polypharmacy strata do not partition the study population like other strata do.
    if (input$ls_drugs_polypharm == "NA") {
      polypharm_gte_5 <- "NA"
      polypharm_gte_10 <- "NA"
    } else if (input$ls_drugs_polypharm == "0-4 drugs") {
      polypharm_gte_5 <- "0"
      polypharm_gte_10 <- "NA"
    } else if (input$ls_drugs_polypharm == ">=5 drugs") {
      polypharm_gte_5 <- "1"
      polypharm_gte_10 <- "NA"
    } else if (input$ls_drugs_polypharm == ">=10 drugs") {
      polypharm_gte_5 <- "NA"
      polypharm_gte_10 <- "1"
    } else {
      stop("ls_drugs_polypharm filter is not valid")
    }
    
    df <- drugs %>% 
      filter(
        cdm_name %in% !!input$ls_drugs_cdm_name,
        # group_name %in% !!input$ls_drugs_cohort,
        # variable_level %in% !!input$ls_drugs_time_window,
        age_group %in% !!input$ls_drugs_age,
        sex %in% !!input$ls_drugs_sex,
        frailty_category %in% !!input$ls_drugs_frailty,
        polypharm_gte_5 %in% !!polypharm_gte_5,
        polypharm_gte_10 %in% !!polypharm_gte_10,
        pre_2020 %in% !!input$ls_drugs_pre_2020) %>% 
      collect() %>% 
      transmute(
        cdm_name = as.factor(cdm_name),
        cohort = as.factor(group_level),
        timeframe = as.factor(variable_level), 
        # strata_name = as.factor(strata_name), 
        # strata_level = as.factor(strata_level),
        covariate = as.factor(variable_name), 
        estimate_type = as.factor(estimate_type),
        estimate_value) %>% 
      group_by(cdm_name, cohort, timeframe, covariate, estimate_type) %>% 
      # ensure there is one row per group. In the belgium data, for some reason we got a few duplicate rows.
      summarise(estimate_value = max(estimate_value, na.rm = T), .groups = "drop")
    
    if (nrow(df) == 0) {
      return(NULL)
    } else {
      df <- df %>% 
        mutate(
          estimate_type = case_when(
            estimate_type == "numeric" ~ "count",
            estimate_type == "percentage" ~ "percent"
          ),
          estimate_value = round(as.numeric(estimate_value), 1)) %>% 
        tidyr::pivot_wider(names_from = "estimate_type", values_from = "estimate_value")
      return(df)
    }
  })
  
  output$ls_drugs <- DT::renderDataTable(ls_drugs(), rownames = FALSE, options = list(pageLength = 10), filter = "top")
  output$ls_drugs_no_data_message <- renderText(ifelse(is.null(ls_drugs()), "There is no data for the current filter selection", ""))
  output$ls_drugs_download <- downloadHandler(
    filename = function() { "large_scale_drug_characterization.csv" },
    content = function(file) { readr::write_csv(ls_drugs(), file) }
  )
  
  ls_drugs_top <- reactive({
   top_20_from_ls(ls_drugs())
  })
  
  output$ls_drugs_top_no_data_message <- renderText(ifelse(is.null(ls_drugs_top()), "There is no data for the current filter selection", ""))
  output$ls_drugs_top <- DT::renderDataTable(ls_drugs_top(), rownames = F, filter = 'top')
  output$ls_drugs_top_download <- downloadHandler(
    filename = function() { "top_drugs_by_database.csv" },
    content = function(file) { readr::write_csv(ls_drugs_top(), file) }
  )
  
  # observe({
  #   print(paste("ls_drugs_cdm_name", input$ls_drugs_cdm_name))
  #   print(paste("ls_drugs_age", input$ls_drugs_age))
  #   print(paste("ls_drugs_sex", input$ls_drugs_sex))
  #   print(paste("ls_drugs_frailty", input$ls_drugs_frailty))
  #   print(paste("ls_drugs_polypharm", input$ls_drugs_polypharm))
  #   print(paste("ls_drugs_pre_2020", input$ls_drugs_pre_2020))
  # })
  

  # Large scale conditions ------
  ls_conditions <- reactive({
    if (input$ls_conditions_polypharm == "NA") {
      polypharm_gte_5 <- "NA"
      polypharm_gte_10 <- "NA"
    } else if (input$ls_conditions_polypharm == "0-4 drugs") {
      polypharm_gte_5 <- "0"
      polypharm_gte_10 <- "NA"
    } else if (input$ls_conditions_polypharm == ">=5 drugs") {
      polypharm_gte_5 <- "1"
      polypharm_gte_10 <- "NA"
    } else if (input$ls_conditions_polypharm == ">=10 drugs") {
      polypharm_gte_5 <- "NA"
      polypharm_gte_10 <- "1"
    } else {
      stop("ls_drugs_polypharm filter is not valid")
    }

    df <- conditions %>%
      filter(
        cdm_name %in% !!input$ls_conditions_cdm_name,
        # group_name %in% !!input$ls_conditions_cohort,
        # variable_level %in% !!input$ls_conditions_time_window,
        age_group %in% !!input$ls_conditions_age,
        sex %in% !!input$ls_conditions_sex,
        frailty_category %in% !!input$ls_conditions_frailty,
        polypharm_gte_5 %in% !!polypharm_gte_5,
        polypharm_gte_10 %in% !!polypharm_gte_10,
        pre_2020 %in% !!input$ls_conditions_pre_2020) %>%
      collect() %>%
      transmute(
        cdm_name = as.factor(cdm_name),
        cohort = as.factor(group_level),
        timeframe = as.factor(variable_level), 
        # strata_name = as.factor(strata_name), 
        # strata_level = as.factor(strata_level),
        covariate = as.factor(variable_name), 
        estimate_type = as.factor(estimate_type),
        estimate_value) %>% 
      group_by(cdm_name, cohort, timeframe, covariate, estimate_type) %>% 
      # ensure there is one row per group. In the Belgium data, for some reason we got a few duplicate rows.
      summarise(estimate_value = max(estimate_value, na.rm = T), .groups = "drop")

    if (nrow(df) == 0) {
      return(NULL)
    } else {
      df <- df %>% 
        mutate(
          estimate_type = case_when(
            estimate_type == "numeric" ~ "count",
            estimate_type == "percentage" ~ "percent"
          ),
          estimate_value = round(as.numeric(estimate_value), 1)) %>% 
        tidyr::pivot_wider(names_from = "estimate_type", values_from = "estimate_value")
      return(df)
    }
  })
  
  output$ls_conditions <- DT::renderDataTable(ls_conditions(), rownames = FALSE, options = list(pageLength = 10), filter = "top")
  output$ls_conditions_no_data_message <- renderText(ifelse(is.null(ls_conditions()), "There is no data for the current filter selection", ""))
  output$ls_conditions_download <- downloadHandler(
    filename = function() { "large_scale_condition_characterization.csv" },
    content = function(file) { readr::write_csv(ls_conditions(), file) }
  )
  
  ls_conditions_top <- reactive({
    top_20_from_ls(ls_conditions())
  })
  
  output$ls_conditions_top_no_data_message <- renderText(ifelse(is.null(ls_conditions_top()), "There is no data for the current filter selection", ""))
  output$ls_conditions_top <- DT::renderDataTable(ls_conditions_top(), rownames = F, filter = 'top')
  output$ls_conditions_top_download <- downloadHandler(
    filename = function() { "top_conditions_by_database.csv" },
    content = function(file) { readr::write_csv(ls_conditions_top(), file) }
  )
  
  selected_plot <- reactive({
    get_plot(prev, input$plot_selection, input$plot_age_group)
  })
  
  output$plot_to_display <- renderPlot(selected_plot(), height = 400L)
  
  plot_save_data <- reactive({
    plot_number <- as.integer(input$plot_selection)
    age <- input$plot_age_group
    width <- if (plot_number %in% c(1,2,4,5)) 20 else 10
    height <- if (plot_number %in% c(1,2,4,5)) 25 else 15
    return(list(num = plot_number,
                age = age,
                width = width,
                height = height))
  })
  
  plot_filename <- reactive({
    if (plot_save_data()[['age']] == "NA") {
      a <- ""
    } else {
      a <- paste0(" age ", plot_save_data()[['age']])
    }
    glue::glue("plot{plot_save_data()[['num']]}{a}.png")
  })
  
  output$download_plot <- downloadHandler(
    filename = function() { plot_filename() },
    content = function(file) { ggplot2::ggsave(file, selected_plot(), 
                                               width = plot_save_data()[['width']], 
                                               height = plot_save_data()[['height']]) }
  )
  
  # output$plot1 <- downloadHandler(
  #   filename = function() { "plot1.png" },
  #   content = function(file) { ggplot2::ggsave(file, get_plot(prev, 1, input$plot_age_group), width = 20, height = 25) }
  # )
  # 
  # output$plot2 <- downloadHandler(
  #   filename = function() { "plot2.png" },
  #   content = function(file) { ggplot2::ggsave(file, get_plot(prev,2, input$plot_age_group), width = 10, height = 15) }
  # )
  # 
  # output$plot3 <- downloadHandler(
  #   filename = function() { "plot3.png" },
  #   content = function(file) { ggplot2::ggsave(file, get_plot(prev,3, input$plot_age_group), width = 10, height = 15) }
  # )
  # 
  # output$plot4 <- downloadHandler(
  #   filename = function() { "plot4.png" },
  #   content = function(file) { ggplot2::ggsave(file, get_plot(prev,4, input$plot_age_group), width = 20, height = 25) }
  # )
  # 
  # output$plot5 <- downloadHandler(
  #   filename = function() { "plot5.png" },
  #   content = function(file) { ggplot2::ggsave(file, get_plot(prev,5, input$plot_age_group), width = 20, height = 25) }
  # )
  # 
  # output$plot6 <- downloadHandler(
  #   filename = function() { "plot6.png" },
  #   content = function(file) { ggplot2::ggsave(file, get_plot(prev,6, input$plot_age_group), width = 10, height = 15) }
  # )
  # 
  # output$plot7 <- downloadHandler(
  #   filename = function() { "plot7.png" },
  #   content = function(file) { ggplot2::ggsave(file, get_plot(prev,7, input$plot_age_group), width = 10, height = 15) }
  # )
  # 
  # output$plot8 <- downloadHandler(
  #   filename = function() { "plot8.png" },
  #   content = function(file) { ggplot2::ggsave(file, get_plot(prev,8, input$plot_age_group), width = 10, height = 15) }
  # )
  # 
  # output$plot9 <- downloadHandler(
  #   filename = function() { "plot9.png" },
  #   content = function(file) { ggplot2::ggsave(file, get_plot(prev,9, input$plot_age_group), width = 10, height = 15) }
  # )
}


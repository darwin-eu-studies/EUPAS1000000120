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

# Dependency management -----
# install.packages("renv") # if not already installed, install renv from CRAN
# renv::restore() # this should prompt you to install the various packages required for the study
# renv::activate()

# Load packages ------
library(CDMConnector)
library(here)
library(dplyr)

minCellCount <- 5

connection_args <- list(
  dbname = Sys.getenv("DBMS_NAME"),
  host = Sys.getenv("DBMS_SERVER"),
  user = Sys.getenv("DBMS_USERNAME")
)
connection_args[["password"]] <- Sys.getenv("DBMS_PASSWORD")
con <- do.call(DBI::dbConnect, c(list(drv = RPostgres::Postgres()), connection_args))

cdm_schema <- Sys.getenv("CDM_SCHEMA")
# we recommend using a prefix here
write_schema <- c(schema = Sys.getenv("WRITE_SCHEMA"), prefix = "c1009_")
cdm_name <- "ipci" # please do not use special characters or spaces in the cdm name. Underscore _ is fine.

output_folder <- here::here(glue::glue("results_{cdm_name}"))

cdm <- cdm_from_con(con,
                    cdm_schema = cdm_schema,
                    write_schema = write_schema,
                    cdm_name = cdm_name)

fs::dir_create(output_folder)
logger <- log4r::logger(threshold = "INFO", 
                        appenders = list(log4r::console_appender(), log4r::file_appender(here(output_folder, "log.txt"))))

log4r::info(logger,"- Getting cdm snapshot")
readr::write_csv(snapshot(cdm), here(output_folder, glue::glue("cdm_snapshot_{cdmName(cdm)}.csv")))

# run the study script
source(here("runStudy.R"))

zip(zipfile = here::here(paste0("results_", cdm_name, ".zip")),
    files = list.files(output_folder, full.names = TRUE, recursive = F), extras = "-j") # the -j flag is used to not include the full paths in the zip file

# drop all the prefixed tables in the cdm. Make sure you use a prefix. Otherwise this will drop all tables in the write schema.
dropTable(cdm, dplyr::everything())
cdm_disconnect(cdm)

# Optional step if you would prefer to send csv instead of a duckdb file
con <- DBI::dbConnect(duckdb::duckdb(), file.path(output_folder, paste0("results_", snakecase::to_snake_case(cdmName(cdm)), ".duckdb")))

csv_output_folder <- here::here(output_folder, paste0("results_", cdm_name, "csv_output"))

fs::dir_create(csv_output_folder)

for (table_name in DBI::dbListTables(con)) {
  dplyr::tbl(con, table_name) |> 
    dplyr::collect() |>
    readr::write_csv(file.path(csv_output_folder, paste0(table_name, ".csv")))
}

zip(zipfile = here::here(paste0("results_", cdm_name, "_csv", ".zip")),
    files = list.files(csv_output_folder, full.names = TRUE, recursive = F), extras = "-j") 

DBI::dbDisconnect(con, shutdown = T)

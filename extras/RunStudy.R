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

# cdm snapshot ----

require(here)
fs::dir_create(here("results"))

cli::cli_text("- Getting cdm snapshot")
readr::write_csv(snapshot(cdm), here("results", paste0(
  "cdm_snapshot_", cdmName(cdm), ".csv"
)))

# Cohort generation ----
cli::cli_text("- Cohort generation")
source(here("study_scripts", "01_cohort_generation.R"))

stopifnot("cancer_cohort" %in% names(cdm))

# characterization analysis -----
if (isTRUE(run_characterisation)){
cli::cli_text("- Running characterization")
  tryCatch({
    source(here("study_scripts", "02_characterization.R"))
  }, error = function(e) {
    writeLines(as.character(e), here("results", "error_characterization.txt"))
  })
}

# prevalence ----
if (isTRUE(run_counts)){
  cli::cli_text("- Running prevalence")
  tryCatch({
    source(here("study_scripts", "03_prevalence.R"))
  }, error = function(e) {
    writeLines(as.character(e), here("results", "error_counts.txt"))
  })
}

# basic counts ----
if (isTRUE(run_counts)){
  cli::cli_text("- Running survival analysis")
  tryCatch({
    source(here("study_scripts", "04_counts.R"))
  }, error = function(e) {
    writeLines(as.character(e), here("results", "error_counts.txt"))
  })
}

# survival analysis ----
if (isTRUE(run_survival)){
cli::cli_text("- Running survival analysis")
  tryCatch({
    source(here("study_scripts", "05_survival.R"))
  }, error = function(e) {
    writeLines(as.character(e), here("results", "error_survival.txt"))
  })
}


# zip results ----
# zip all results
zip(zipfile = file.path(here("results", paste0("results_", db_name, ".zip"))),
files = list.files(here("results"), full.names = TRUE, recursive = TRUE))

print("Done!")
print("-- If all has worked, there should now be a zip folder with your results in the results folder to share")
print("-- Thank you for running the study!")

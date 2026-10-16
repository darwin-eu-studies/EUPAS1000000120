# Instructions

## Prerequisites

Install R and the system requirements needed by the packages in `renv.lock`. The study is designed to run against an OMOP CDM database.

## Configure the database connection

Before running the study, set these environment variables in your local R environment:

- `DBMS_NAME`
- `DBMS_SERVER`
- `DBMS_USERNAME`
- `DBMS_PASSWORD`
- `CDM_SCHEMA`
- `WRITE_SCHEMA`

Use a write schema where temporary and study output tables can be created. The study code uses the prefix `c1009_` for write-schema tables.

## Run the study

1. Clone this repository.
2. Open `CodeToRun.R` in RStudio or another R environment.
3. Restore the package environment with `renv::restore()` if needed.
4. Review the database connection variables and `cdm_name` value in `CodeToRun.R`.
5. Run `CodeToRun.R`.

The script writes study outputs to a local folder named for the configured CDM and creates zipped output files for sharing through the agreed DARWIN EU study process.

## Clean-up

After execution, `CodeToRun.R` drops the study tables with the configured write-schema prefix and disconnects from the database.

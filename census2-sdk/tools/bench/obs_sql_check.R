# The SQL path of census2_sdk.obs_sql from R, through the duckdb package.
#
# The SDK writes the catalog test and the filter as SQL. This script runs the same two
# statements in R: the catalog test over the catalog's Parquet file gives the fragments
# of the datasets that can match, and one query counts the selected cells of those
# fragments by donor. It also counts all cells with the file list in a variable, which
# is how a SQL script can give a table function the catalog's list.
#
#   python -c "..."  # write the SQL of obs_sql.parse(text).text and obs_sql.catalog_test
#   Rscript obs_sql_check.R s3://my-bucket/spike/catalogs/catalog-103.parquet queries.json
#
# queries.json maps a name to {"where": <filter SQL>, "catalog": <catalog test SQL>}.
# The region comes from AWS_REGION, and the credentials from DuckDB's credential chain.
args <- commandArgs(trailingOnly = TRUE)
catalog <- args[1]
sqlfile <- args[2]
region <- Sys.getenv("AWS_REGION", "us-east-1")
t0 <- Sys.time()
library(duckdb)
library(jsonlite)
cat("duckdb R", as.character(packageVersion("duckdb")), "\n")
con <- dbConnect(duckdb())
invisible(dbExecute(con, "INSTALL httpfs; LOAD httpfs; INSTALL aws; LOAD aws"))
invisible(dbExecute(con, sprintf(
  "CREATE SECRET (TYPE s3, PROVIDER credential_chain, REGION '%s')", region)))
cat("connect + extensions", format(Sys.time() - t0), "\n")
queries <- fromJSON(sqlfile)
for (name in names(queries)) {
  q <- queries[[name]]
  t1 <- Sys.time()
  kept <- dbGetQuery(con, sprintf(paste(
    "SELECT replace(store_uri, '.zarr', '.obs.parquet') AS f FROM read_parquet('%s')",
    "WHERE organism = 'homo_sapiens' AND %s"), catalog, q$catalog))
  if (nrow(kept) == 0) {
    cat(name, ": no candidates\n")
    next
  }
  files <- paste0("'", kept$f, "'", collapse = ", ")
  counts <- dbGetQuery(con, sprintf(paste(
    "SELECT donor_id, count(*) AS n FROM read_parquet([%s], union_by_name = true)",
    "WHERE %s GROUP BY 1"), files, q$where))
  cat(name, ": candidates", nrow(kept), "donors", nrow(counts), "cells", sum(counts$n),
      "in", format(Sys.time() - t1), "\n")
}
t2 <- Sys.time()
invisible(dbExecute(con, sprintf(paste(
  "SET VARIABLE files = (SELECT list(replace(store_uri, '.zarr', '.obs.parquet'))",
  "FROM read_parquet('%s') WHERE organism = 'homo_sapiens')"), catalog)))
all <- dbGetQuery(con, "SELECT count(*) AS n FROM read_parquet(getvariable('files'))")
cat("all cells, file list from SET VARIABLE:", all$n, "in", format(Sys.time() - t2), "\n")

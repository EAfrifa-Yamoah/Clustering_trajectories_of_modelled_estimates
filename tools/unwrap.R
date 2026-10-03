#!/usr/bin/env Rscript
# usage: Rscript unwrap.R <outname> [file]  -- unwrap a connector result (JSON {result: "<json>"}) to a slim CSV
args <- commandArgs(trailingOnly = TRUE); out <- args[1]
dir <- "/root/.claude/projects/-home-claude/57d1fe51-c9de-54fe-bff5-075ff6ab1088/tool-results"
done <- file.path("/home/claude/gbd", "processed.txt")
seen <- if (file.exists(done)) readLines(done) else character()
f <- if (length(args) > 1) args[2] else {
  fs <- list.files(dir, pattern = "get_multi_location_trend", full.names = TRUE)
  fs <- setdiff(fs, seen); fs[which.min(file.info(fs)$mtime)] }
txt <- paste(readLines(f, warn = FALSE), collapse = "")
# strip outer {"result": "..."} wrapper without a JSON parser: unescape the inner string
inner <- sub('^\\s*\\{"result"\\s*:\\s*"', "", txt); inner <- sub('"\\s*\\}\\s*$', "", inner)
inner <- gsub('\\\\"', '"', inner); inner <- gsub("\\\\n", "", inner); inner <- gsub('\\\\\\\\', '\\\\', inner)
# records are objects inside "data": [ ... ]; extract fields by regex
recs <- regmatches(inner, gregexpr("\\{[^{}]*\"location_id\"[^{}]*\\}", inner))[[1]]
getf <- function(k) { v <- sub(paste0('.*"', k, '":\\s*"?([^",}]*)"?.*'), "\\1", recs, perl = TRUE); v[!grepl(paste0('"', k, '"'), recs)] <- NA; v }
keep <- c("location_id", "location_name", "year", "cause_name", "measure", "rate_age_std_combined", "rate_combined", "rate_lower_combined", "rate_upper_combined")
df <- as.data.frame(lapply(keep, getf), stringsAsFactors = FALSE); names(df) <- keep
df$location_name <- gsub("\\\\u00f4", "\u00f4", df$location_name)
write.csv(df, file.path("/home/claude/gbd", paste0(out, ".csv")), row.names = FALSE)
cat(sprintf("%s: %d rows, %d locations, years %s-%s\n", out, nrow(df), length(unique(df$location_id)), min(df$year), max(df$year)))
cat(f, file = done, sep = "\n", append = TRUE)

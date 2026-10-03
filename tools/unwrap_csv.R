#!/usr/bin/env Rscript
# usage: Rscript unwrap_csv.R <outname> <file>  -- unwrap a fmt='csv' connector result to a slim CSV
args <- commandArgs(trailingOnly = TRUE); out <- args[1]; f <- args[2]
txt <- paste(readLines(f, warn = FALSE), collapse = "")
inner <- sub('^\\s*\\{"result"\\s*:\\s*"', "", txt); inner <- sub('"\\s*\\}\\s*$', "", inner)
inner <- gsub("\\\\n", "\n", inner); inner <- gsub('\\\\"', '"', inner); inner <- gsub("\\\\u00f4", "ô", inner)
d <- read.csv(text = inner, stringsAsFactors = FALSE)
keep <- c("location_id","location_name","year","cause_name","measure","rate_age_std_combined","rate_combined","rate_lower_combined","rate_upper_combined")
d <- d[, keep]; write.csv(d, file.path("/home/claude/gbd", paste0(out, ".csv")), row.names = FALSE)
cat(sprintf("%s: %d rows, %d locations, years %s-%s\n", out, nrow(d), length(unique(d$location_id)), min(d$year), max(d$year)))

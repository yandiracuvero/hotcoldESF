# Clean the covariate file TODO_MADRID119m.txt ---------------------------------------------
#
# 119,940 rows (one per building / address point, `idcasa`), 138 columns.
# Removes:
#   - exact duplicate columns (identical in every row): the second name is dropped
#       sola16_64 (= Solo1664), Tot (= SoloTot), solaM64 (= SoloM64), solaM64H (= SoloM64H),
#       solaM64M (= SoloM64M), D_VPSuc (= D_AlcRes), ProSEnvej (= v15a19)
#   - constant columns: PerSolas
# Notes written to application/results/covariates_cleaning.csv:
#   - ProSEnvej is identical to v15a19 (a count of people aged 15-19, values 0-285), so
#     it is not an over-ageing proportion: the column was probably copied by mistake.
#   - D_AlcRes and D_VPSuc are identical: one of the two labels is wrong; kept as D_AlcRes.
#   - 12 rows have no census variables (NA in all 77 census columns).
#
# Output: application/data/covariates_clean.rds (data frame)
# Run from the repo root: Rscript application/00_clean_covariates.R

h <- read.csv("application/data/TODO_MADRID119m.txt", sep = ";", dec = ",", check.names = FALSE)
key  <- vapply(h, function(x) paste(x, collapse = "|"), "")
grp  <- split(names(h), match(key, unique(key)))
grp  <- lapply(grp, function(g) if (setequal(g, c("ProSEnvej", "v15a19"))) c("v15a19", "ProSEnvej") else g)
dups <- unlist(lapply(grp[lengths(grp) > 1], `[`, -1))
const <- names(h)[vapply(h, function(x) length(unique(x)) == 1, TRUE)]
drop <- unique(c(dups, const))
notes <- rbind(
  do.call(rbind, lapply(grp[lengths(grp) > 1], function(g)
    data.frame(column = g[-1], action = "dropped",
               reason = if (g[1] == "v15a19") "identical to v15a19 (a count, 0-285), so not a proportion: probably copied by mistake"
                        else paste("identical to", g[1])))),
  data.frame(column = const, action = "dropped", reason = "constant"),
  data.frame(column = "D_AlcRes", action = "kept",
             reason = "identical to D_VPSuc: one of the two labels is wrong; meaning to verify"),
  data.frame(column = "(rows)", action = "kept",
             reason = sprintf("%d rows have NA in all census columns", sum(!complete.cases(h)))))
dir.create("application/results", showWarnings = FALSE)
write.csv(notes, "application/results/covariates_cleaning.csv", row.names = FALSE)
clean <- h[, setdiff(names(h), drop)]
saveRDS(clean, "application/data/covariates_clean.rds")
message(ncol(h), " -> ", ncol(clean), " columns; dropped: ", paste(drop, collapse = ", "))

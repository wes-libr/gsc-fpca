# GSC-FPCA refits with a separate level for each site x study year, estimated from controls: absorbs any change at a site
# in any year, recorded or not (software upgrades, hardware, protocol). Cells with fewer than 10 control scans are merged
# with the previous year at the same site. Usage: Rscript scanner_check/refit_siteyear.R <worker 0|1>
args <- commandArgs(TRUE); wk <- as.numeric(args[1])
source("app2_core.R")
M$wave <- M$visit - 1
o0 <- modifyList(defaults, list(outcome = "sri24_icv_brain_volume")); d0 <- prep_data(o0)
ever <- tapply(d0$binge12, d0$id, max); ctl0 <- d0[ever[as.character(d0$id)] == 0, ]
cells <- list()
for (s in c("A", "B", "C", "D", "E")) {
  w <- sort(unique(d0$visit[d0$site == s] - 1)); n <- sapply(w, function(k) sum(ctl0$site == s & ctl0$visit - 1 == k))
  grp <- integer(length(w)); g <- 0
  for (j in seq_along(w)) { if (j == 1 || n[j] >= 10) g <- g + 1; grp[j] <- g }
  cells[[s]] <- data.frame(site = s, wave = w, n_ctl = n, cell = paste0("sy_", s, grp))
}
cells <- do.call(rbind, cells); print(cells)
key <- paste(M$site, M$wave); ck <- paste(cells$site, cells$wave)
cell_of <- cells$cell[match(key, ck)]
lv <- unique(cells$cell); lv <- setdiff(lv, "sy_E1")   # reference: site E, baseline
for (c in lv) M[[c]] <- as.numeric(!is.na(cell_of) & cell_of == c)
assign("M", M, envir = globalenv())
dir.create("scanner_check/out", recursive = TRUE, showWarnings = FALSE)
write.csv(cells, file.path("scanner_check/out", "site_year_cells.csv"), row.names = FALSE)
icv_t <- c("sexF", "race_labelBlack", "race_labelAsian", "race_labelWhite")
specs <- list(list(label = "icv_siteyear", Q = 3, outcome = "sri24_icv_brain_volume", covs_t = icv_t, covs_c = lv),
              list(label = "sf_siteyear", Q = 3, covs_c = lv),
              list(label = "icv_siteyear_placebo", Q = 3, outcome = "sri24_icv_brain_volume", covs_t = icv_t, covs_c = lv, placebo = TRUE),
              list(label = "gm_siteyear", Q = 3, outcome = "sri24_icv_brain_gm", covs_c = lv))
for (i in seq_along(specs)) {
  if ((i - 1) %% 2 != wk) next
  sp <- specs[[i]]; t0 <- Sys.time()
  r <- tryCatch(run_app(sp), error = function(e) { cat("ERROR", sp$label, conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) { a <- r$att[r$att$method == "gsc1" & r$att$bin %in% c("all", "placebo"), ]
    cat(sp$label, "ncov", length(sp$covs_c), "secs", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "ATT", round(a$est, 4), "(", round(a$lo, 4), round(a$hi, 4), ") mm3", round(a$est * a$sc), "\n") }
}

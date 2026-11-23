# Empirical example with the first-onset treatment definition (the default in app2_core.R): fits with k = 1 and k = 2,
# used for the comparison across k and the two-FPC column of the sensitivity table. Run with: Rscript queue_first_onset.R <worker>.
source("app2_core.R"); stopifnot(defaults$discont == "itt")
dir.create("claims_first_onset", showWarnings = FALSE); dir.create("res2", showWarnings = FALSE)
sites <- c("siteB", "siteC", "siteD", "siteE")
base_t <- c("icv_baseline", "sexF", "race_labelBlack", "race_labelAsian", "race_labelWhite")
nc_t <- c("sexF", "race_labelBlack", "race_labelAsian", "race_labelWhite")
long <- list(nloop = 5500, burnin = 500); short <- list(nloop = 3500, burnin = 500)
specs <- list(); add <- function(...) specs[[length(specs) + 1]] <<- list(...)
for (s in 1:4) do.call(add, c(list(label = sprintf("k2_s%d", s), Q = 2, seed = 20261104 + 2000 * s, psloo = as.numeric(s <= 2)), long))
do.call(add, c(list(label = "primary_s1", Q = 1, seed = 20262104, psloo = 1), long))
do.call(add, c(list(label = "k3_s1", Q = 3, seed = 20264104, psloo = 1), long))
for (q in 3:5) do.call(add, c(list(label = paste0("kk", q), Q = q, seed = 20261104, psloo = 1), short))
add(label = "k2_placebo", Q = 2, placebo = TRUE); add(label = "placebo", placebo = TRUE); add(label = "kk3_placebo", Q = 3, placebo = TRUE)
sens <- list(left = list(outcome = "sf_left"), right = list(outcome = "sf_right"), disc_drop = list(discont = "drop"), disc_censor = list(discont = "censor"),
  ctl_nolater = list(controls = "no_later"), ctl_zero = list(controls = "zero"), site = list(covs_c = sites),
  rich = list(covs_t = c(base_t, "ses_yoe_b", "fh_alc_b", "pds_b"), covs_c = sites),
  rich_sub = list(covs_t = c(base_t, "ses_yoe_b", "fh_alc_b", "pds_b", "mj_life_b", "tob_life_b"), covs_c = sites),
  couse = list(covs_c = c(sites, "mj_life_t", "tob_life_t")), win25 = list(amax = 25), win25_censor = list(amax = 25, discont = "censor"))
for (n in names(sens)) { do.call(add, c(list(label = paste0("k2_", n), Q = 2), sens[[n]])); do.call(add, c(list(label = n, Q = 1), sens[[n]])) }
add(label = "k2_age_bins", Q = 2, age_bins = TRUE)
add(label = "k2_negctl_icv", Q = 2, outcome = "sri24_icv_brain_volume", covs_t = nc_t)
add(label = "negctl_icv", outcome = "sri24_icv_brain_volume", covs_t = nc_t)
add(label = "k3_negctl_icv", Q = 3, outcome = "sri24_icv_brain_volume", covs_t = nc_t)
add(label = "k2_negctl_icv_site", Q = 2, outcome = "sri24_icv_brain_volume", covs_t = nc_t, covs_c = sites)
add(label = "k2_total_gm", Q = 2, outcome = "sri24_icv_brain_gm")
for (r in readRDS("rois.rds")) add(label = paste0("k2_roi_", r), Q = 2, outcome = r, nloop = 1500, burnin = 500)
for (r in 1:50) add(label = sprintf("k2_fake_%02d", r), Q = 2, fake = TRUE, seed = 9000 + r, nloop = 1000, burnin = 500)
for (sp in specs) {
  f <- file.path("res2", paste0(sp$label, ".rds"))
  if (file.exists(f)) next
  if (!dir.create(file.path("claims_first_onset", sp$label), showWarnings = FALSE)) next
  t0 <- Sys.time()
  r <- tryCatch(run_app(sp), error = function(e) { cat("ERROR", sp$label, conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) {
    a <- r$att[r$att$method == "gsc1" & r$att$bin %in% c("all", "placebo"), ]
    cat(format(Sys.time(), "%H:%M:%S"), sp$label, "N", r$meta$Ntr, r$meta$Nco, "post", a$n, "secs", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "ATT", round(a$est, 4), "(", round(a$lo, 3), round(a$hi, 3), ")\n")
  }
}
cat("queue finished\n")

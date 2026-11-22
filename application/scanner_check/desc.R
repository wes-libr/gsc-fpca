# Scanner checks for the ICV negative control. Writes aggregate tables only.
source("app2_core.R")
suppressPackageStartupMessages({library(sandwich); library(lmtest)})
OUT <- "scanner_check/out"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
o <- modifyList(defaults, list(outcome = "sri24_icv_brain_volume", covs_t = c("sexF", "race_labelBlack", "race_labelAsian", "race_labelWhite")))
d <- prep_data(o); d <- add_covariates(d, o)
d$wave <- d$visit - 1
d$model <- ifelse(is.na(d$scanner_model), "NA", d$scanner_model)
d$unit <- paste(d$site, d$model)
# new-scanner indicators, within site
d$A_prisma <- as.numeric(d$site == "A" & grepl("Prisma", d$model))
d$B_creator <- as.numeric(d$site == "B" & d$model == "SIGNA_Creator")
d$C_uhp <- as.numeric(d$site == "C" & d$model == "SIGNA_UHP")
d$D_prisma <- as.numeric(d$site == "D" & grepl("Prisma", d$model))
d$D_xa30 <- as.numeric(d$site == "D" & d$model == "Prisma_Fit_XA30")
ind <- c("A_prisma", "B_creator", "C_uhp", "D_prisma", "D_xa30")
d$new_any <- as.numeric(rowSums(d[, ind]) > 0)
ever <- tapply(d$binge12, d$id, max); d$ever <- as.numeric(ever[as.character(d$id)])
d$post <- d$binge12 == 1
lastpre <- tapply(d$visit_age[!d$post], d$id[!d$post], max)
d$yr <- floor(d$visit_age - lastpre[as.character(d$id)] + .5)
has_pre <- d$id %in% as.numeric(names(lastpre))
cat("rows", nrow(d), "units", length(unique(d$id)), "exposed units with pre", length(unique(d$id[d$ever == 1 & has_pre])), "post rows", sum(d$post & has_pre), "\n")
sc <- sd(d$sri24_icv_brain_volume); cat("ICV SD (scale of the ATT)", sc, "\n")

# 1. Scanner model by site and study year, analysis rows
tab1 <- as.data.frame(table(site = d$site, model = d$model, wave = d$wave)); tab1 <- tab1[tab1$Freq > 0, ]
write.csv(tab1, file.path(OUT, "scanner_by_site_wave.csv"), row.names = FALSE)

# 2. ICV step at each scanner change, within person, never-treated units only
ctl <- d[d$ever == 0, ]
ctl$sexF <- as.numeric(ctl$sex == "F")
f0 <- lm(sri24_icv_brain_volume ~ factor(id) + ns(visit_age, 4) + ns(visit_age, 4):sexF, data = ctl)
ind <- ind[sapply(ind, function(v) sum(ctl[[v]]) > 0)]; cat("indicators with control scans:", ind, "\n")
f1 <- update(f0, as.formula(paste(". ~ . +", paste(ind, collapse = " + "))))
V <- vcovCL(f1, cluster = ~id)
ct <- coeftest(f1, vcov. = V)[ind, ]
nsw <- sapply(ind, function(v) { s <- tapply(ctl[[v]], ctl$id, function(x) length(unique(x)) > 1); sum(s) })
nrow_new <- sapply(ind, function(v) sum(ctl[[v]]))
lab <- c(A_prisma = "A: Trio to Prisma", B_creator = "B: MR750 to SIGNA Creator (temporary)", C_uhp = "C: MR750 to SIGNA UHP", D_prisma = "D: Trio to Prisma", D_xa30 = "D: Prisma software XA30")
step <- data.frame(change = lab[ind],
  term = ind, est_mm3 = ct[, 1], se_mm3 = ct[, 2], lo_mm3 = ct[, 1] - 1.96 * ct[, 2], hi_mm3 = ct[, 1] + 1.96 * ct[, 2], p = ct[, 4],
  est_sd = ct[, 1] / sc, lo_sd = (ct[, 1] - 1.96 * ct[, 2]) / sc, hi_sd = (ct[, 1] + 1.96 * ct[, 2]) / sc, scans_on_new = nrow_new, units_switching = nsw)
print(step, digits = 3)
write.csv(step, file.path(OUT, "icv_step_at_scanner_change.csv"), row.names = FALSE)
cat("controls rows", nrow(ctl), "units", length(unique(ctl$id)), "\n")
# same for total gray matter and superior frontal, for comparison
for (y in c("sri24_icv_brain_gm", "sf_bilat")) {
  ff <- lm(as.formula(paste(y, "~ factor(id) + ns(visit_age, 4) + ns(visit_age, 4):sexF +", paste(ind, collapse = " + "))), data = ctl)
  cc <- coeftest(ff, vcov. = vcovCL(ff, cluster = ~id))[ind, ]; s2 <- sd(d[[y]], na.rm = TRUE)
  st <- data.frame(outcome = y, term = ind, est = cc[, 1], lo = cc[, 1] - 1.96 * cc[, 2], hi = cc[, 1] + 1.96 * cc[, 2], p = cc[, 4], est_sd = cc[, 1] / s2, lo_sd = (cc[, 1] - 1.96 * cc[, 2]) / s2, hi_sd = (cc[, 1] + 1.96 * cc[, 2]) / s2)
  print(st, digits = 3); write.csv(st, file.path(OUT, paste0("step_", y, ".csv")), row.names = FALSE)
}

# 3. Residual ICV by study year and site (never-treated units), from the model without scanner terms
ctl$res <- resid(f0)
r3 <- do.call(rbind, lapply(split(ctl, list(ctl$site, ctl$wave), drop = TRUE), function(x) {
  data.frame(site = x$site[1], wave = x$wave[1], n = nrow(x), mean_res = mean(x$res), se = sd(x$res) / sqrt(nrow(x)), prop_new = mean(x$new_any), main_model = names(sort(table(x$model), decreasing = TRUE))[1])
}))
r3 <- r3[order(r3$site, r3$wave), ]
write.csv(r3, file.path(OUT, "residual_icv_by_site_wave.csv"), row.names = FALSE)

# 4. Scanner mix of the treated: own pre vs post visits, and controls at the same site and age (+/- 1 year)
ex <- d[d$ever == 1 & has_pre, ]
pc <- sapply(seq_len(nrow(ex)), function(i) { w <- ctl$site == ex$site[i] & abs(ctl$visit_age - ex$visit_age[i]) <= 1; if (sum(w) < 5) NA else mean(ctl$new_any[w]) })
ex$p_ctl <- pc
bins <- c(sprintf("pre %d", -3:0), sprintf("year %d", 1:5), "all post")
mix <- do.call(rbind, lapply(bins, function(b) {
  w <- if (b == "all post") ex$post else if (grepl("^pre", b)) !ex$post & ex$yr == as.numeric(sub("pre ", "", b)) else ex$post & ex$yr == as.numeric(sub("year ", "", b))
  data.frame(bin = b, n = sum(w), prop_new_exposed = mean(ex$new_any[w]), prop_new_controls_same_site_age = mean(ex$p_ctl[w], na.rm = TRUE), mean_wave = mean(ex$wave[w]))
}))
mix$excess <- mix$prop_new_exposed - mix$prop_new_controls_same_site_age
print(mix, digits = 3)
write.csv(mix, file.path(OUT, "scanner_mix_exposed_vs_controls.csv"), row.names = FALSE)
# study year of treated post visits vs control visits at the same site and age
ex$wave_ctl <- sapply(seq_len(nrow(ex)), function(i) { w <- ctl$site == ex$site[i] & abs(ctl$visit_age - ex$visit_age[i]) <= 1; if (sum(w) < 5) NA else mean(ctl$wave[w]) })
wv <- data.frame(group = c("exposed pre-exposure visits", "exposed post-exposure visits"), mean_wave = c(mean(ex$wave[!ex$post]), mean(ex$wave[ex$post])),
  mean_wave_controls_same_site_age = c(mean(ex$wave_ctl[!ex$post], na.rm = TRUE), mean(ex$wave_ctl[ex$post], na.rm = TRUE)))
print(wv); write.csv(wv, file.path(OUT, "study_year_exposed_vs_controls.csv"), row.names = FALSE)
# treated units by site
print(table(site = ex$site[!duplicated(ex$id)])); print(table(site = ctl$site[!duplicated(ctl$id)]))
saveRDS(d, "scanner_check/d_icv.rds")

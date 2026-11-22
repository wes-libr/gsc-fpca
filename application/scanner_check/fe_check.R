# Simple within-person check: change after onset relative to each person's own earlier scans, with a smooth age trend
# estimated from all participants. Rows and units are those of the ICV negative control. Writes aggregate results only.
source("app2_core.R")
suppressPackageStartupMessages({library(sandwich); library(lmtest)})
OUT <- "scanner_check/out"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
o <- modifyList(defaults, list(outcome = "sri24_icv_brain_volume"))
d <- prep_data(o); d$sexF <- as.numeric(d$sex == "F")
md <- ifelse(is.na(d$scanner_model), "", d$scanner_model)
d$A_prisma <- as.numeric(d$site == "A" & grepl("Prisma", md)); d$B_creator <- as.numeric(d$site == "B" & md == "SIGNA_Creator")
d$C_uhp <- as.numeric(d$site == "C" & md == "SIGNA_UHP"); d$D_prisma <- as.numeric(d$site == "D" & grepl("Prisma", md))
pre_ids <- unique(d$id[d$binge12 == 0]); d <- d[d$id %in% pre_ids, ]
d$post <- d$binge12; d$siteyear <- factor(paste(d$site, d$visit - 1))
res <- list()
for (y in c("sri24_icv_brain_volume", "sri24_icv_brain_gm", "sri24_icv_brain_wm", "sri24_icv_brain_csf", "sf_bilat")) {
  dd <- d[!is.na(d[[y]]), ]; s <- sd(dd[[y]])
  for (adj in c("none", "scanner", "site_year")) {
    fm <- paste(y, "~ factor(id) + ns(visit_age, 4) + ns(visit_age, 4):sexF + post", switch(adj, none = "", scanner = "+ A_prisma + B_creator + C_uhp + D_prisma", site_year = "+ siteyear"))
    f <- lm(as.formula(fm), data = dd); ct <- coeftest(f, vcov. = vcovCL(f, cluster = ~id))["post", ]
    res[[length(res) + 1]] <- data.frame(outcome = y, adjust = adj, rows = nrow(dd), units = length(unique(dd$id)), post_rows = sum(dd$post), est_mm3 = ct[1], lo_mm3 = ct[1] - 1.96 * ct[2], hi_mm3 = ct[1] + 1.96 * ct[2], p = ct[4], est_sd = ct[1] / s, lo_sd = (ct[1] - 1.96 * ct[2]) / s, hi_sd = (ct[1] + 1.96 * ct[2]) / s)
  }
}
res <- do.call(rbind, res); rownames(res) <- NULL
print(res, digits = 3)
write.csv(res, file.path(OUT, "within_person_post_change.csv"), row.names = FALSE)

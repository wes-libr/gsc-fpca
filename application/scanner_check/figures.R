# Figures and summary table for the ICV scanner check. Reads aggregate results only.
suppressPackageStartupMessages(library(ggplot2))
R2 <- "res2"; PUB <- "res2"; OUT <- "scanner_check/out"; FIG <- "scanner_check/fig"
for (p in c(OUT, FIG)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid <- "#e6e5e1"
pal <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", yellow = "#eda100", magenta = "#e87ba4")
th <- theme_minimal(base_size = 11) + theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2), panel.grid.minor = element_blank(),
  panel.grid.major = element_line(colour = grid, linewidth = .3), plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(colour = ink2, size = 9.5),
  plot.caption = element_text(colour = ink2, size = 8, hjust = 0), legend.position = "top", legend.justification = "left", strip.text = element_text(face = "bold", hjust = 0),
  plot.background = element_rect(fill = "white", colour = NA))
wr <- function(x, w = 150) paste(strwrap(x, w), collapse = "\n")
att <- function(f, dir = R2) { a <- read.csv(file.path(dir, paste0(f, "_att.csv"))); a <- a[a$method == "gsc1", ]; a$mm3 <- a$est * a$sc; a$lo3 <- a$lo * a$sc; a$hi3 <- a$hi * a$sc; a }

# A. ICV: observed minus synthetic control by year relative to onset
bx <- function(b) ifelse(grepl("^pre", b), as.numeric(sub("pre ", "", b)), as.numeric(sub("year ", "", b)))
lvA <- c("As in the paper", "+ site + scanner model", "+ a level for each site and study year")
A <- rbind(transform(att("k3_negctl_icv"), model = lvA[1]), transform(att("icv_scan"), model = lvA[2]), transform(att("icv_siteyear"), model = lvA[3]))
A <- A[A$bin != "all", ]; A$x <- bx(A$bin); A$period <- ifelse(A$x <= 0, "pre", "post")
A$model <- factor(A$model, levels = lvA)
pl <- rbind(transform(att("icv_placebo"), model = lvA[1]), transform(att("icv_scan_placebo"), model = lvA[2]), transform(att("icv_siteyear_placebo"), model = lvA[3]))
pl$x <- 0; pl$model <- factor(pl$model, levels = levels(A$model))
nlab <- A[A$model == "As in the paper", c("x", "n")]
pA <- ggplot(A, aes(x, mm3, colour = model)) +
  annotate("rect", xmin = -3.5, xmax = .5, ymin = -Inf, ymax = Inf, fill = "#f3f2ef") +
  geom_hline(yintercept = 0, colour = ink2, linewidth = .4) +
  geom_linerange(aes(ymin = lo3, ymax = hi3), position = position_dodge(.6), linewidth = .8) +
  geom_point(position = position_dodge(.6), size = 2.4) +
  geom_point(data = pl, aes(x = -.5 + .11 * (as.numeric(model) - 2), y = mm3), shape = 5, size = 2.4, stroke = 1, show.legend = FALSE) +
  geom_linerange(data = pl, aes(x = -.5 + .11 * (as.numeric(model) - 2), ymin = lo3, ymax = hi3), linewidth = .6, linetype = "22", show.legend = FALSE) +
  geom_text(data = nlab, aes(x, y = -26500, label = paste0("n=", n)), inherit.aes = FALSE, colour = ink2, size = 3) +
  annotate("text", x = -1.5, y = 8500, label = "Before onset (fit to these scans)", colour = ink2, size = 3.2) +
  annotate("text", x = 3, y = 8500, label = "After onset (predicted, not fit)", colour = ink2, size = 3.2) +
  scale_colour_manual(values = unname(pal[c("blue", "orange", "aqua")]), name = NULL) +
  scale_x_continuous(breaks = -3:5, labels = c("-3", "-2", "-1", "0\n(last scan\nbefore onset)", 1:5)) +
  scale_y_continuous(labels = function(v) format(v, big.mark = ","), limits = c(-27500, 9500)) +
  labs(x = "Years from the last scan before onset of heavy binge drinking", y = expression("Observed minus synthetic control ICV (mm"^3*")"),
       title = "Intracranial volume: treated participants vs. their synthetic controls",
       subtitle = wr("Mean difference with 95% posterior predictive interval. Hollow diamonds between -1 and 0: the last pre-onset scan held out of the fit (placebo).", 120),
       caption = wr("115 treated and 513 control participants, ages 12-21. GSC-FPCA with two components. Scanner model: indicators for the recorded scanner changes at sites A, B, C and D. Site and study year: a separate level for every site in every study year (36 levels, sparse late years pooled), estimated from the controls, which absorbs any change at a site in any year, recorded or not.", 160)) + th
ggsave(file.path(FIG, "1_icv_effect_by_year.png"), pA, width = 9, height = 5.6, dpi = 200)

# B. Scanner model by site and study year (analysis rows)
B <- read.csv(file.path(OUT, "scanner_by_site_wave.csv"))
B$model <- factor(B$model, levels = c("MR750", "SIGNA_Creator", "SIGNA_UHP", "TrioTim", "Prisma_Fit"),
  labels = c("GE MR750", "GE SIGNA Creator", "GE SIGNA UHP", "Siemens Trio", "Siemens Prisma"))
pB <- ggplot(B, aes(factor(wave), Freq, fill = model)) + geom_col(width = .8, colour = "white", linewidth = .3) + facet_wrap(~ paste("Site", site), nrow = 1) +
  scale_fill_manual(values = unname(pal[c("blue", "yellow", "magenta", "orange", "aqua")]), name = NULL) +
  labs(x = "Study year (0 = baseline)", y = "Scans in the analysis", title = "Scanner model by site and study year",
       subtitle = wr("Scans of the 628 participants in the ICV analysis (ages 12-21). Site D's Prisma XA30 software upgrade falls after age 21 for everyone here."),
       caption = "Scan dates are almost entirely missing in FullData_Y9 (visit_date is filled for under 10% of rows, as month and year only), so study year stands in for calendar time.") + th +
  theme(legend.position = "top")
ggsave(file.path(FIG, "2_scanner_by_site_year.png"), pB, width = 9.5, height = 4.2, dpi = 200)

# C. Within-person residual ICV by site and study year, never-treated participants
C <- read.csv(file.path(OUT, "residual_icv_by_site_wave.csv")); C <- C[C$n >= 5, ]
C$main_model <- factor(C$main_model, levels = c("MR750", "SIGNA_UHP", "TrioTim", "Prisma_Fit"), labels = c("GE MR750", "GE SIGNA UHP", "Siemens Trio", "Siemens Prisma"))
pC <- ggplot(C, aes(wave, mean_res)) + geom_hline(yintercept = 0, colour = ink2, linewidth = .4) +
  geom_line(colour = grid, linewidth = .6) +
  geom_linerange(aes(ymin = mean_res - 1.96 * se, ymax = mean_res + 1.96 * se, colour = main_model), linewidth = .7) +
  geom_point(aes(colour = main_model), size = 2.2) + facet_wrap(~ paste("Site", site), nrow = 1) +
  scale_colour_manual(values = unname(pal[c("blue", "magenta", "orange", "aqua")]), name = "Most scans that year on") +
  scale_x_continuous(breaks = 0:8) + scale_y_continuous(labels = function(v) format(v, big.mark = ",")) +
  labs(x = "Study year (0 = baseline)", y = expression("Mean residual ICV (mm"^3*")"), title = "ICV over the study within never-treated participants",
       subtitle = wr("Residuals from a model with a separate level for each person and a smooth age curve by sex. Bars: 95% CI. Years with fewer than 5 scans left out.", 160),
       caption = wr("513 never-treated participants, 2,133 scans. A software or hardware change would show as a step within a site. A drift common to all sites is not visible here, because the age curve absorbs it. Late years have few scans, from the participants enrolled youngest.", 190)) + th
ggsave(file.path(FIG, "3_residual_icv_by_site_year.png"), pC, width = 10, height = 4.4, dpi = 200)

# D. Overall estimates under each adjustment
fe <- read.csv(file.path(OUT, "within_person_post_change.csv"))
row <- function(outc, lab, f, dir = R2) { a <- att(f, dir); a <- a[a$bin == "all", ]; data.frame(outcome = outc, spec = lab, est = a$mm3, lo = a$lo3, hi = a$hi3) }
ferow <- function(outc, y, adj, lab) { x <- fe[fe$outcome == y & fe$adjust == adj, ]; data.frame(outcome = outc, spec = lab, est = x$est_mm3, lo = x$lo_mm3, hi = x$hi_mm3) }
D <- rbind(
  row("Intracranial volume", "As in the paper", "k3_negctl_icv", PUB), row("Intracranial volume", "+ site", "k3_negctl_icv_site", PUB), row("Intracranial volume", "+ site + scanner model", "icv_scan"),
  row("Intracranial volume", "+ site + scanner + study year", "icv_scan_wave"), row("Intracranial volume", "+ site + scanner + study year by site", "icv_scan_sitewave"), row("Intracranial volume", "+ a level for each site and study year", "icv_siteyear"),
  ferow("Intracranial volume", "sri24_icv_brain_volume", "none", "Within-person regression"), ferow("Intracranial volume", "sri24_icv_brain_volume", "scanner", "Within-person regression + scanner"), ferow("Intracranial volume", "sri24_icv_brain_volume", "site_year", "Within-person regression + site and study year"),
  row("Superior frontal gray matter", "As in the paper", "k3_s1", PUB), row("Superior frontal gray matter", "+ site", "k3_site", PUB), row("Superior frontal gray matter", "+ site + scanner model", "sf_scan"),
  row("Superior frontal gray matter", "+ site + scanner + study year", "sf_scan_wave"), row("Superior frontal gray matter", "+ a level for each site and study year", "sf_siteyear"),
  ferow("Superior frontal gray matter", "sf_bilat", "none", "Within-person regression"), ferow("Superior frontal gray matter", "sf_bilat", "scanner", "Within-person regression + scanner"), ferow("Superior frontal gray matter", "sf_bilat", "site_year", "Within-person regression + site and study year"))
lv <- c("As in the paper", "+ site", "+ site + scanner model", "+ site + scanner + study year", "+ site + scanner + study year by site", "+ a level for each site and study year", "Within-person regression", "Within-person regression + scanner", "Within-person regression + site and study year")
D$spec <- factor(D$spec, levels = rev(lv)); D$kind <- ifelse(grepl("Within", D$spec), "Within-person regression", "GSC-FPCA")
write.csv(D, file.path(OUT, "estimates_by_adjustment.csv"), row.names = FALSE)
pD <- ggplot(D, aes(est, spec, colour = kind)) + geom_vline(xintercept = 0, colour = ink2, linewidth = .4) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = .8) + geom_point(size = 2.6) +
  facet_wrap(~ outcome, scales = "free_x") + scale_colour_manual(values = unname(pal[c("blue", "aqua")]), name = NULL) +
  scale_x_continuous(labels = function(v) format(v, big.mark = ",")) +
  labs(x = expression("Effect after onset (mm"^3*")"), y = NULL, title = "Overall estimates with more adjustment for scanner and time",
       subtitle = wr("GSC-FPCA: mean of all 182 post-onset scans, 95% posterior predictive interval. Within-person regression: change after onset relative to each person's own earlier scans, with a smooth age curve by sex; 95% CI clustered by person.", 115),
       caption = "Superior frontal 'As in the paper' and '+ site' are the paper's runs (another seed); the others use the default seed.") + th
ggsave(file.path(FIG, "4_estimates_by_adjustment.png"), pD, width = 10.5, height = 5.4, dpi = 200)

# E. ICV split into gray matter, white matter and CSF
E <- rbind(row("ICV", "GSC-FPCA", "k3_negctl_icv"), row("Gray matter", "GSC-FPCA", "comp_gm"), row("White matter", "GSC-FPCA", "comp_wm"), row("CSF", "GSC-FPCA", "comp_csf"),
  row("ICV", "GSC-FPCA + site + scanner", "icv_scan"), row("Gray matter", "GSC-FPCA + site + scanner", "comp_gm_scan"), row("White matter", "GSC-FPCA + site + scanner", "comp_wm_scan"), row("CSF", "GSC-FPCA + site + scanner", "comp_csf_scan"),
  ferow("ICV", "sri24_icv_brain_volume", "none", "Within-person regression"), ferow("Gray matter", "sri24_icv_brain_gm", "none", "Within-person regression"),
  ferow("White matter", "sri24_icv_brain_wm", "none", "Within-person regression"), ferow("CSF", "sri24_icv_brain_csf", "none", "Within-person regression"))
names(E)[1:2] <- c("tissue", "method"); E$tissue <- factor(E$tissue, levels = rev(c("ICV", "Gray matter", "White matter", "CSF")))
E$method <- factor(E$method, levels = c("GSC-FPCA", "GSC-FPCA + site + scanner", "Within-person regression"))
write.csv(E, file.path(OUT, "icv_components.csv"), row.names = FALSE)
pE <- ggplot(E, aes(est, tissue, colour = method)) + geom_vline(xintercept = 0, colour = ink2, linewidth = .4) +
  geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(-.55), linewidth = .8) + geom_point(position = position_dodge(-.55), size = 2.4) +
  scale_colour_manual(values = unname(pal[c("blue", "orange", "aqua")]), name = NULL) + guides(colour = guide_legend(nrow = 1)) +
  scale_x_continuous(labels = function(v) format(v, big.mark = ",")) +
  labs(x = expression("Effect after onset (mm"^3*")"), y = NULL, title = "Where the ICV difference comes from",
       subtitle = wr("SRI24 ICV equals gray matter + white matter + CSF (to within 0.4%). If the ICV mask were unaffected, gray matter loss would be matched by a CSF gain.", 110),
       caption = wr("Same participants, scans and covariates as the ICV negative control (sex and race; no baseline ICV). Each tissue is fit separately, so the GSC-FPCA estimates need not add up exactly.", 150)) + th
ggsave(file.path(FIG, "5_icv_components.png"), pE, width = 9, height = 4.8, dpi = 200)
print(D, digits = 3); print(E, digits = 3)
P <- rbind(transform(att("icv_placebo"), model = "paper"), transform(att("icv_scan_placebo"), model = "scanner"), transform(att("icv_siteyear_placebo"), model = "site_year"))
write.csv(P[, c("model", "bin", "n", "units", "mm3", "lo3", "hi3", "est", "lo", "hi")], file.path(OUT, "icv_placebo.csv"), row.names = FALSE); print(P[, c("model", "bin", "n", "mm3", "lo3", "hi3", "est", "lo", "hi")], digits = 3)

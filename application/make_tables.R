# Sensitivity table, ROI table, and a CSV of numbers for Section 4, from the application results.
res <- function(l) { f <- file.path("res2", paste0(l, ".rds")); if (file.exists(f)) readRDS(f) else NULL }
neg <- function(x, d = 3) { v <- sprintf(paste0("%.", d, "f"), abs(x)); ifelse(x < 0 & as.numeric(v) != 0, paste0("$-$", v), v) }
ci <- function(a, b, d = 3) sprintf("[%s,\\;%s]", neg(a, d), neg(b, d))
pick <- function(r, meth = "gsc1", bin = "all") { a <- r$att; a[a$method == meth & a$bin == bin, ] }
cell <- function(r, meth = "gsc1", bin = "all") { if (is.null(r)) return("--"); x <- pick(r, meth, bin); sprintf("%s %s", neg(x$est), ci(x$lo, x$hi)) }
cell2 <- function(r, meth = "gsc1", bin = "all") cell(r, meth, bin)
specs <- list(
  list("Primary specification ($k=3$)", "k3_s1", "k2_s1"),
  list("Prior mean zero", "k3_s1", "k2_s1", "gsc0"),
  list("\\emph{Number of FPCs}", NULL),
  list("$k=1$", "primary_s1", NULL), list("$k=2$", "k2_s1", NULL), list("$k=4$", "kk4", NULL), list("$k=5$", "kk5", NULL),
  list("\\emph{Outcome}", NULL),
  list("Left hemisphere only", "k3_left", "k2_left"), list("Right hemisphere only", "k3_right", "k2_right"),
  list("\\emph{Treatment definition}", NULL),
  list("Persistent heavy drinking only", "k3_disc_drop", "k2_disc_drop"), list("Censored at discontinuation", "k3_disc_censor", "k2_disc_censor"),
  list("\\emph{Controls}", NULL),
  list("Excluding later treated", "k3_ctl_nolater", "k2_ctl_nolater"), list("No binge drinking", "k3_ctl_zero", "k2_ctl_zero"),
  list("\\emph{Covariates added}", NULL),
  list("Site", "k3_site", "k2_site"), list("Site, education, family history, puberty", "k3_rich", "k2_rich"), list("Above, baseline cannabis, tobacco", "k3_rich_sub", "k2_rich_sub"), list("Site, cumulative cannabis, tobacco", "k3_couse", "k2_couse"),
  list("\\emph{Age range}", NULL),
  list("Post-treatment visits to age 20.5", "k3_age_bins", NULL, "gsc1", "all age<=20.5"),
  list("Ages to 25", "k3_win25", "k2_win25"), list("Ages to 25, censored at discontinuation", "k3_win25_censor", "k2_win25_censor"))
rows <- character()
for (sp in specs) {
  if (is.null(sp[[2]])) { rows <- c(rows, sprintf("\\multicolumn{4}{l}{%s} \\\\", sp[[1]])); next }
  r2 <- res(sp[[2]]); r1 <- if (is.null(sp[[3]])) NULL else res(sp[[3]]); m <- if (length(sp) >= 4) sp[[4]] else "gsc1"; b <- if (length(sp) >= 5) sp[[5]] else "all"
  nn <- if (is.null(r2)) "--" else sprintf("%d/%d/%d", r2$meta$Ntr, r2$meta$Nco, pick(r2, m, b)$n)
  rows <- c(rows, sprintf("%s & %s & %s & %s \\\\", sp[[1]], nn, cell(r2, m, b), if (is.null(sp[[3]])) "--" else cell(r1, m, "all")))
}
af <- read.csv("out_final_k3/att_table.csv"); xf <- af[af$label == "all post-exposure", ]
r0 <- res("k3_s1"); rows[1] <- sprintf("Primary specification ($k=3$) & %d/%d/%d & %s %s & %s \\\\", r0$meta$Ntr, r0$meta$Nco, xf$n, neg(xf$est), ci(xf$pred_lo, xf$pred_hi), cell(res("k2_s1")))
writeLines(rows, "tab_sens_rows.tex")
# ROI table
rois <- readRDS("rois.rds")
tab <- do.call(rbind, lapply(rois, function(ro) { r <- if (ro == "superiorfrontal_grayvol") res("k3_s1") else res(paste0("k3_roi_", ro)); if (is.null(r)) return(NULL); x <- pick(r); data.frame(roi = ro, est = x$est, lo = x$lo, hi = x$hi, pred_sd = x$pred_sd, p_neg = x$p_neg, sc = r$meta$sc) }))
if (!is.null(tab)) {
  af <- read.csv("out_final_k3/att_table.csv"); xf <- af[af$label == "all post-exposure", ]; js <- tab$roi == "superiorfrontal_grayvol"
  tab$est[js] <- xf$est; tab$lo[js] <- xf$pred_lo; tab$hi[js] <- xf$pred_hi; tab$p_neg[js] <- xf$post_prob_neg; tab$pred_sd[js] <- xf$pred_sd
  tab <- tab[order(tab$est), ]; tab$rank <- seq_len(nrow(tab))
  write.csv(tab, "roi_table.csv", row.names = FALSE)
  # rows of the regional table are written by hier_regions.R, which adds the shrunken estimates
}
# calibration with pseudo-treatments
fk <- list.files("res2", pattern = "^k3_fake_.*rds$", full.names = TRUE)
if (length(fk)) {
  z <- do.call(rbind, lapply(fk, function(f) { x <- pick(readRDS(f)); data.frame(est = x$est, lo = x$lo, hi = x$hi, n = x$n, units = x$units) }))
  z0 <- do.call(rbind, lapply(fk, function(f) { x <- pick(readRDS(f), "gsc0"); data.frame(est = x$est) }))
  cat(sprintf("fake runs %d mean %.4f sd %.4f min %.3f max %.3f cover %d avgwidth %.3f mean_n %.1f units %.1f | gsc0 mean %.4f\n", nrow(z), mean(z$est), sd(z$est), min(z$est), max(z$est), sum(z$lo <= 0 & z$hi >= 0), mean(z$hi - z$lo), mean(z$n), mean(z$units), mean(z0$est)))
  write.csv(z, "fake_table.csv", row.names = FALSE)
}

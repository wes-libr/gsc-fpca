# Prints the numbers quoted in Section 4 and Appendix B for the primary specification (first onset, k = 3). Run after make_tables.R.
RES <- "res2"; res <- function(l) readRDS(file.path(RES, paste0(l, ".rds")))
pick <- function(r, meth = "gsc1", bin = "all") { a <- r$att; a[a$method == meth & a$bin == bin, ] }
f <- function(x, meth = "gsc1", bin = "all", u = FALSE) { r <- res(x); p <- pick(r, meth, bin); s <- if (u) r$meta$sc else 1; d <- if (u) 1 else 3
  sprintf("%s: %s [%s, %s] n=%d Ntr=%d Nco=%d", x, formatC(p$est * s, format = "f", digits = d), formatC(p$lo * s, format = "f", digits = d), formatC(p$hi * s, format = "f", digits = d), p$n, r$meta$Ntr, r$meta$Nco) }
a <- read.csv("out_final_k3/att_table.csv"); sc <- a$sc[1]; cat("scale", sc, "\n")
for (i in seq_len(nrow(a))) cat(sprintf("%-18s n=%3d units=%3d est %.3f PI [%.3f, %.3f] fit [%.3f, %.3f] | mm3 %.1f [%.1f, %.1f]\n", a$label[i], a$n[i], a$units[i], a$est[i], a$pred_lo[i], a$pred_hi[i], a$fit_lo[i], a$fit_hi[i], a$est[i] * sc, a$pred_lo[i] * sc, a$pred_hi[i] * sc))
cat("chains:", sapply(1:4, function(s) round(pick(res(sprintf("k3_s%d", s)))$est, 4)), "\n")
for (x in c("k3_placebo", "k2_placebo", "placebo")) { r <- res(x); p <- pick(r, "gsc1", "placebo"); cat(x, sprintf("%.3f [%.3f, %.3f] n=%d units=%d mm3 %.1f [%.1f, %.1f]\n", p$est, p$lo, p$hi, p$n, p$units, p$est * r$meta$sc, p$lo * r$meta$sc, p$hi * r$meta$sc)) }
cat(f("k3_s1", "gsc0"), "\n")
for (x in c("primary_s1", "k2_s1", "k3_s1", "kk4", "kk5", "k3_negctl_icv", "k2_negctl_icv", "negctl_icv", "k3_negctl_icv_site", "k3_total_gm", "k2_total_gm")) cat(f(x), "\n")
for (x in c("k3_negctl_icv", "k3_total_gm", "k3_negctl_icv_site")) cat(f(x, u = TRUE), "| year3", f(x, bin = "year 3"), "\n")
r1 <- res("primary_s1"); r3 <- res("k3_s1"); for (b in c("pre -3", "pre 0")) cat(b, "k1", round(pick(r1, "gsc1", b)$est, 3), "k3", round(pick(r3, "gsc1", b)$est, 3), "\n")
nf <- do.call(rbind, lapply(1:4, function(s) res(sprintf("k3_s%d", s))$NuF)); ag <- r3$age_grid
print(round(rbind(age = ag, mean = colMeans(nf), lo = apply(nf, 2, quantile, .025), hi = apply(nf, 2, quantile, .975)), 3))
Z <- read.csv("fake_table.csv"); cat(sprintf("fake n=%d mean %.3f sd %.3f cover %d width %.3f below %d above %d\n", nrow(Z), mean(Z$est), sd(Z$est), sum(Z$lo <= 0 & Z$hi >= 0), mean(Z$hi - Z$lo), sum(Z$hi < 0), sum(Z$lo > 0)))
R <- read.csv("roi_table.csv"); cat("roi neg", sum(R$est < 0), "of", nrow(R), "excl", sum(R$hi < 0 | R$lo > 0), "sf rank", R$rank[R$roi == "superiorfrontal_grayvol"], "median", round(median(R$est), 3), "range", round(range(R$est), 3), "\n"); print(head(R[, c("roi", "est", "lo", "hi")], 8), digits = 3)
set.seed(20261105); ch <- lapply(1:4, function(s) res(sprintf("k3_s%d", s))); nb <- function(k) a$n[a$label == sprintf("year %d", k)]
yb <- function(b, n) unlist(lapply(ch, function(r) { x <- r$draws[[paste("gsc1", b)]]; x + rnorm(length(x), 0, sqrt(r$sigma2 / n)) }))
Y1 <- yb("year 1", nb(1)); Y2 <- yb("year 2", nb(2)); Y3 <- yb("year 3", nb(3)); q <- function(x) sprintf("%.3f [%.3f, %.3f]", mean(x), quantile(x, .025), quantile(x, .975)); cat("y3-y2", q(Y3 - Y2), " y3-y1", q(Y3 - Y1), "\n")
m <- read.csv("out_final_k3/lme.csv"); print(m[, 1:3], digits = 3); print(read.csv("out_final_k3/fpc_eigen.csv"), digits = 3); print(read.csv("out_final_k3/demo.csv")); print(read.csv("out_final_k3/visits.csv")); print(read.csv("out_final_k3/coef_bands.csv"), digits = 3)
dg <- read.csv("diag_k3.csv"); cat("rhat max", max(dg$rhat), "ess min", min(dg$ess_total), "\n"); print(dg[dg$rhat > 1.005, 1:6], digits = 4)

if (file.exists("hier_summary.csv")) { cat("hierarchical model for the regions\n"); print(round(read.csv("hier_summary.csv", row.names = 1), 4)); H <- read.csv("hier_regions.csv"); cat("shrunken range", round(range(H$shrunk), 3), "negative", sum(H$shrunk < 0), "interval excludes zero", sum(H$s_hi < 0 | H$s_lo > 0), "\n") }

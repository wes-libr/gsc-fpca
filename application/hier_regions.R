# Two-stage hierarchical model for the 34 regional ATT estimates (Section 4 and the regional table of the supplement).
# Run after queue_first_onset_k3.R and make_tables.R. Reads res2/k3_roi_*.rds and res2/k3_units_sf.rds (participant level, private)
# and roi_table.csv; writes hier_regions.csv, hier_summary.csv, hier_corr.csv (aggregate) and tab_roi_rows.tex.
# Level 1: theta_hat ~ N(theta, Sigma), Sigma = D R D with D the posterior predictive SDs and R the correlation of the
#          regional estimates, estimated from participant-level post-treatment differences (clustered by participant).
# Level 2: theta_j ~ N(mu, tau^2). Flat prior on mu; half-normal(0, 0.1) prior on tau (others for sensitivity).
# Superior frontal uses the primary (bilateral, four-chain) estimate, as in the regional table.
res <- function(l) readRDS(file.path("res2", paste0(l, ".rds")))
rois <- readRDS("rois.rds"); J <- length(rois); nm <- sub("_grayvol$", "", rois)
lab <- ifelse(nm == "superiorfrontal", "k3_units_sf", paste0("k3_roi_", rois))
U <- lapply(lab, function(l) res(l)$unit_post)
stopifnot(all(sapply(U, function(u) identical(u$id, U[[1]]$id))))
Dm <- sapply(U, function(u) u$diff); id <- U[[1]]$id; n <- nrow(Dm); G <- length(unique(id))
est_u <- colMeans(Dm)
Uc <- rowsum(sweep(Dm, 2, est_u), id)                      # participant sums of centered differences
Vcl <- crossprod(Uc) / n^2 * G / (G - 1); R <- cov2cor(Vcl)
tab <- read.csv("roi_table.csv"); tab <- tab[match(rois, tab$roi), ]   # estimates and predictive SDs used in the paper's table
th <- tab$est; s <- tab$pred_sd; names(th) <- nm
cat("check: unit-level means vs table estimates, max abs diff", round(max(abs(est_u - th)), 4), "\n")
cat("cluster SE / predictive SD: median", round(median(sqrt(diag(Vcl)) / s), 2), "range", round(range(sqrt(diag(Vcl)) / s), 2), "\n")
cat("correlations among regional estimates: median", round(median(R[upper.tri(R)]), 2), "range", round(range(R[upper.tri(R)]), 2), "\n")
fit <- function(Sigma, prior = function(t) dnorm(t, 0, 0.1), S = 20000, seed = 1) {
  set.seed(seed); one <- rep(1, J); tg <- seq(0.0005, 0.4, by = 0.0005)
  lp <- sapply(tg, function(t) { V <- Sigma + diag(t^2, J); Vi <- solve(V); a <- sum(Vi); mh <- sum(Vi %*% th) / a; r <- th - mh
    log(prior(t)) - 0.5 * as.numeric(determinant(V)$modulus) - 0.5 * log(a) - 0.5 * sum(r * (Vi %*% r)) })
  w <- exp(lp - max(lp)); w <- w / sum(w); ts <- sample(tg, S, TRUE, w)
  mu <- numeric(S); TH <- matrix(NA, S, J)
  for (t in unique(ts)) { ix <- which(ts == t); V <- Sigma + diag(t^2, J); Vi <- solve(V); a <- sum(Vi); mh <- sum(Vi %*% th) / a
    m <- rnorm(length(ix), mh, sqrt(1 / a)); mu[ix] <- m
    W <- diag(t^2, J) - t^4 * Vi; W <- (W + t(W)) / 2; e <- eigen(W, symmetric = TRUE); L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), J)
    B <- t^2 * Vi
    for (k in seq_along(ix)) TH[ix[k], ] <- m[k] + as.vector(B %*% (th - m[k])) + as.vector(L %*% rnorm(J)) }
  list(mu = mu, tau = ts, TH = TH)
}
q <- function(x) c(mean = mean(x), lo = quantile(x, .025, names = FALSE), hi = quantile(x, .975, names = FALSE))
Sig <- diag(s) %*% R %*% diag(s)
f <- fit(Sig)
cat("mu", round(q(f$mu), 4), "P(mu<0)", round(mean(f$mu < 0), 4), "\n"); cat("tau", round(q(f$tau), 4), "median", round(median(f$tau), 4), "\n")
out <- data.frame(region = nm, est = th, lo = tab$lo, hi = tab$hi, shrunk = colMeans(f$TH), s_lo = apply(f$TH, 2, quantile, .025), s_hi = apply(f$TH, 2, quantile, .975), p_neg = colMeans(f$TH < 0))
out <- out[order(out$shrunk), ]; write.csv(out, "hier_regions.csv", row.names = FALSE)
print(out, digits = 3, row.names = FALSE)
cat("shrunken: negative", sum(out$shrunk < 0), "interval excludes zero", sum(out$s_hi < 0 | out$s_lo > 0), "P(neg)>0.975", sum(out$p_neg > .975), "\n")
cat("contrast with mean: P(theta_j < mu) > .975 for", paste(nm[colMeans(f$TH < f$mu) > .975], collapse = ","), "| P(theta_j > mu) > .975 for", paste(nm[colMeans(f$TH > f$mu) > .975], collapse = ","), "\n")
cat("P(all 34 negative)", round(mean(apply(f$TH < 0, 1, all)), 3), " mean number negative", round(mean(rowSums(f$TH < 0)), 1), "\n")
sens <- list(halfnormal_0.05 = function(t) dnorm(t, 0, 0.05), halfcauchy_0.05 = function(t) dcauchy(t, 0, 0.05), uniform = function(t) rep(1, length(t)))
rows <- list(primary = c(q(f$mu), q(f$tau), excl = sum(out$s_hi < 0 | out$s_lo > 0)))
for (k in names(sens)) { g <- fit(Sig, sens[[k]]); lo <- apply(g$TH, 2, quantile, .025); hi <- apply(g$TH, 2, quantile, .975); rows[[k]] <- c(q(g$mu), q(g$tau), excl = sum(hi < 0 | lo > 0)) }
g <- fit(diag(s^2)); lo <- apply(g$TH, 2, quantile, .025); hi <- apply(g$TH, 2, quantile, .975); rows[["independent_errors"]] <- c(q(g$mu), q(g$tau), excl = sum(hi < 0 | lo > 0))
g <- fit(Vcl); lo <- apply(g$TH, 2, quantile, .025); hi <- apply(g$TH, 2, quantile, .975); rows[["cluster_covariance"]] <- c(q(g$mu), q(g$tau), excl = sum(hi < 0 | lo > 0))
S <- do.call(rbind, rows); colnames(S) <- c("mu", "mu_lo", "mu_hi", "tau", "tau_lo", "tau_hi", "n_excl_zero"); print(round(S, 4)); write.csv(S, "hier_summary.csv")
write.csv(round(R, 3), "hier_corr.csv")
# rows of the regional table: unshrunk estimate and 95% PI, shrunken estimate and 95% interval, ordered by the unshrunk estimate
neg <- function(x, d = 3) { v <- sprintf(paste0("%.", d, "f"), abs(x)); ifelse(x < 0 & as.numeric(v) != 0, paste0("$-$", v), v) }
ci <- function(a, b, d = 3) sprintf("[%s,\\;%s]", neg(a, d), neg(b, d))
pretty <- c(bankssts = "Banks of superior temporal sulcus", caudalanteriorcingulate = "Caudal anterior cingulate", caudalmiddlefrontal = "Caudal middle frontal", cuneus = "Cuneus", entorhinal = "Entorhinal", fusiform = "Fusiform", inferiorparietal = "Inferior parietal", inferiortemporal = "Inferior temporal", isthmuscingulate = "Isthmus cingulate", lateraloccipital = "Lateral occipital", lateralorbitofrontal = "Lateral orbitofrontal", lingual = "Lingual", medialorbitofrontal = "Medial orbitofrontal", middletemporal = "Middle temporal", parahippocampal = "Parahippocampal", paracentral = "Paracentral", parsopercularis = "Pars opercularis", parsorbitalis = "Pars orbitalis", parstriangularis = "Pars triangularis", pericalcarine = "Pericalcarine", postcentral = "Postcentral", posteriorcingulate = "Posterior cingulate", precentral = "Precentral", precuneus = "Precuneus", rostralanteriorcingulate = "Rostral anterior cingulate", rostralmiddlefrontal = "Rostral middle frontal", superiorfrontal = "Superior frontal", superiorparietal = "Superior parietal", superiortemporal = "Superior temporal", supramarginal = "Supramarginal", frontalpole = "Frontal pole", temporalpole = "Temporal pole", transversetemporal = "Transverse temporal", insula = "Insula")
o2 <- out[order(out$est), ]
writeLines(sprintf("%s & %s & %s & %s & %s \\\\", pretty[o2$region], neg(o2$est), ci(o2$lo, o2$hi), neg(o2$shrunk), ci(o2$s_lo, o2$s_hi)), "tab_roi_rows.tex")

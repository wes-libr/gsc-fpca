# Placebo for the ICV negative control: each treated participant's last pre-onset scan is held out and predicted,
# without and with the recorded scanner changes. Usage: Rscript scanner_check/placebo.R <worker 0|1>
args <- commandArgs(TRUE); wk <- as.numeric(args[1])
source("app2_core.R")
M$wave <- M$visit - 1
md <- ifelse(is.na(M$scanner_model), "", M$scanner_model)
M$A_prisma <- as.numeric(M$site == "A" & grepl("Prisma", md)); M$B_creator <- as.numeric(M$site == "B" & md == "SIGNA_Creator")
M$C_uhp <- as.numeric(M$site == "C" & md == "SIGNA_UHP"); M$D_prisma <- as.numeric(M$site == "D" & grepl("Prisma", md))
for (s in c("A", "B", "C", "D", "E")) M[[paste0("wave_", s)]] <- M$wave * (M$site == s)
assign("M", M, envir = globalenv())
sites <- c("siteB", "siteC", "siteD", "siteE"); scan <- c("A_prisma", "B_creator", "C_uhp", "D_prisma")
icv_t <- c("sexF", "race_labelBlack", "race_labelAsian", "race_labelWhite")
specs <- list(); add <- function(...) specs[[length(specs) + 1]] <<- list(...)
add(label = "icv_placebo", Q = 3, outcome = "sri24_icv_brain_volume", covs_t = icv_t, placebo = TRUE)
add(label = "icv_scan_placebo", Q = 3, outcome = "sri24_icv_brain_volume", covs_t = icv_t, covs_c = c(sites, scan), placebo = TRUE)
dir.create("res2", showWarnings = FALSE)
for (i in seq_along(specs)) {
  if ((i - 1) %% 2 != wk) next
  sp <- specs[[i]]; f <- file.path("res2", paste0(sp$label, ".rds")); if (file.exists(f)) next
  t0 <- Sys.time()
  r <- tryCatch(run_app(sp), error = function(e) { cat("ERROR", sp$label, conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) { a <- r$att[r$att$method == "gsc1" & r$att$bin == "all", ]
    cat(sp$label, "N", r$meta$Ntr, r$meta$Nco, "post", r$meta$n_post, "secs", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "ATT", round(a$est, 4), "(", round(a$lo, 4), round(a$hi, 4), ") sc", round(r$meta$sc), "\n") }
}
cat("worker", wk, "done\n")

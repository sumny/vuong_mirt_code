###################################################################################################
# Packages and versions, original results or rerun the simulations?
###################################################################################################

### for exact results use the versions specified below
### figures will only look nice if exported as .eps or .pdf

# R version 3.5.1
require(MASS) # version 7.3-51.1
require(mirt) # version 1.29
require(nonnest2) # version 0.5-2
require(SimDesign) # version 1.13

### work with the original results or rerun the simulations?
#setwd("~/vuong_mirt_code/orig") # work with the original results
#setwd("~/vuong_mirt_code/repl") # rerun the simulations

# Last mod: March/23/2018, LS (revision updates)



###################################################################################################
# Simulation 1.1: Code and results or rerun
###################################################################################################

if(length(grep("vuong_mirt_code/orig$", getwd()))) {
  cat("Working with the original simulation result files.\n")
} else if(length(grep("vuong_mirt_code/repl$", getwd()))) {
  if("sim_graded_gpcm_base.rds" %in% dir() | "sim_graded_gpcm_base-results" %in% dir()) {
    stop("New result files exist. You should inspect them prior to a rerun.\n")
  } else {
    cat("Rerunning the simulation. This may take a very long time.\n")
  }
  ### Settings
  N <- c(500, 1000, 2000) # persons
  M <- 10 # items
  K <- 4 # categories
  Design <- expand.grid(N = N, M = M, K = K)
  rm(N, M, K)

  ### SimDesign functions
  Generate <- function(condition, fixed_objects = NULL) {
    Attach(condition)
    dat <- matrix(rbinom(N * M, K - 1, 0.5), N, M)
    stopifnot(all(sapply(1:M, function(x) {
      length(unique(dat[, x]))
    }) == K))
    colnames(dat) <- 1:M
    return(dat)
  }

  Analyse <- function(condition, dat, fixed_objects = NULL) {
    ## EM, BFGS, Ramsay, Oakes, 61 quadpts, 1e-4 TOL
    graded <- mirt(dat, 1, "graded", SE = TRUE, technical = list(NCYCLES = 5000))
    gpcm <- mirt(dat, 1, "gpcm", SE = TRUE, technical = list(NCYCLES = 5000))
    return(list(vt = vuongtest(graded, gpcm),
      AIC_graded = extract.mirt(graded, "AIC"),
      AIC_gpcm = extract.mirt(gpcm, "AIC"),
      BIC_graded = extract.mirt(graded, "BIC"),
      BIC_gpcm = extract.mirt(gpcm, "BIC"), M2_graded = M2(graded),
      M2_gpcm = M2(gpcm), graded_check = c(extract.mirt(graded, "converged"),
      extract.mirt(graded, "condnum"), extract.mirt(graded, "secondordertest")),
      gpcm_check = c(extract.mirt(gpcm, "converged"),
      extract.mirt(gpcm, "condnum"), extract.mirt(gpcm, "secondordertest"))))
  }

  Summarise <- function(condition, results, fixed_objects = NULL) {
    return(list(vt_dist = mean(sapply(results, function(x) x$vt$p_omega < 0.05)),
      vt_graded = mean(sapply(results, function(x) x$vt$p_LRT$A < 0.05)),
      vt_gpcm = mean(sapply(results, function(x) x$vt$p_LRT$B < 0.05)),
      AIC = mean(sapply(results, function(x) (x$AIC_graded - x$AIC_gpcm) < 0)),
      BIC = mean(sapply(results, function(x) (x$BIC_graded - x$BIC_gpcm) < 0)),
      M2_graded = mean(sapply(results, function(x) x$M2_graded$p < 0.05)),
      M2_gpcm = mean(sapply(results, function(x) x$M2_gpcm$p < 0.05)),
      graded_cn = mean(sapply(results, function(x) x$graded_check[1])),
      graded_so = mean(sapply(results, function(x) x$graded_check[3])),
      gpcm_cn = mean(sapply(results, function(x) x$gpcm_check[1])),
      gpcm_so = mean(sapply(results, function(x) x$gpcm_check[3]))))
  }

  ### Run
  results <- runSimulation(Design, replications = 1000,
    generate = Generate, analyse = Analyse, summarise = Summarise,
    packages = c("mirt", "nonnest2"), save = TRUE, save_results = TRUE,
    save_seeds = FALSE, parallel = TRUE, progress = TRUE, verbose = FALSE,
    filename = "sim_graded_gpcm_base",
    save_details = list(save_results_dirname = "sim_graded_gpcm_base-results"))

  ### clean up
  #rm(list = c("Analyse", "Design", "Generate", "results", "Summarise"))
} else {
  stop("Please set the working directory as described above.\n")
}



###################################################################################################
# Simulation 1.1: Table 1, Figure 1 and text results
###################################################################################################

res11 <- readRDS("sim_graded_gpcm_base.rds")
#attributes(res11)

### check convergence
min(res11$graded_cn) # 0.924
min(res11$gpcm_cn)   # 0.901
min(res11$graded_so) # 1
min(res11$gpcm_so)   # 1
res11$AIC_graded <- res11$AIC
res11$AIC_gpcm <- 1 - res11$AIC
results11 <- res11[, 1:14]
results11 <- cbind(results11, "AIC_graded" = res11$AIC_graded, "AIC_gpcm" = res11$AIC_gpcm)

### results of replications in which both models converged
res11_oc <- res11
res11_oc$cn_both <- NULL
AIC_results <- list()
for(i in 1:dim(res11_oc)[1]) {
  pth <- gsub("-i", paste0("-", i), "sim_graded_gpcm_base-results/results-row-i.rds")
  temp <- readRDS(pth)
  cn <- which(sapply(temp[[2]], function(x) x$graded_check[1] == 1 & x$gpcm_check[1] == 1))
  res11_oc$cn_both[i] <- length(cn) / length(temp[[2]])
  ## vt_dist
  res11_oc[i, 4] <- mean(sapply(temp[[2]], function(x) x$vt$p_omega)[cn] < 0.05)
  ## vt_graded
  res11_oc[i, 5] <- mean(sapply(temp[[2]], function(x) x$vt$p_LRT$A)[cn] < 0.05)
  ## vt_gpcm
  res11_oc[i, 6] <- mean(sapply(temp[[2]], function(x) x$vt$p_LRT$B)[cn] < 0.05)
  ## AIC
  res11_oc[i, 7] <- mean(sapply(temp[[2]], function(x) (x$AIC_graded - x$AIC_gpcm) < 0)[cn])
  AIC_results[[i]] <- summary(sapply(temp[[2]], function(x) abs(x$AIC_graded - x$AIC_gpcm)))
  ## BIC
  res11_oc[i, 8] <- mean(sapply(temp[[2]], function(x) (x$BIC_graded - x$BIC_gpcm) < 0)[cn])
  ## M2_graded
  res11_oc[i, 9] <- mean(sapply(temp[[2]], function(x) x$M2_graded$p)[cn] < 0.05)
  ## M2_gpcm
  res11_oc[i, 10] <- mean(sapply(temp[[2]], function(x) x$M2_gpcm$p)[cn] < 0.05)
  ## graded_cn
  res11_oc[i, 11] <- mean(sapply(temp[[2]], function(x) x$graded_check[1] == 1)[cn])
  ## graded_so
  res11_oc[i, 12] <- mean(sapply(temp[[2]], function(x) x$graded_check[3] == 1)[cn])
  ## gpcm_cn
  res11_oc[i, 13] <- mean(sapply(temp[[2]], function(x) x$graded_check[1] == 1)[cn])
  ## gpcm_so
  res11_oc[i, 14] <- mean(sapply(temp[[2]], function(x) x$graded_check[3] == 1)[cn])
}
res11_oc$AIC_graded <- res11_oc$AIC
res11_oc$AIC_gpcm <- 1 - res11_oc$AIC
results11_oc <- res11_oc[, 1:14]
results11_oc <- cbind(results11_oc, "AIC_graded" = res11_oc$AIC_graded,
  "AIC_gpcm" = res11_oc$AIC_gpcm, "cn_both" = res11_oc$cn_both)
which.min(results11_oc$cn_both) # 1 N = 500

### vt_graded / vt_gpcm only if vt_dist significant
lrt_v_ds <- t(sapply(1:dim(res11_oc)[1], function(i) {
  pth <- gsub("-i", paste0("-", i), "sim_graded_gpcm_base-results/results-row-i.rds")
  temp <- readRDS(pth)
  cn <- which(sapply(temp[[2]], function(x) x$graded_check[1] == 1 & x$gpcm_check[1] == 1))
  dist <- which(sapply(temp[[2]], function(x) x$vt$p_omega < 0.05)[cn])
  vt_graded <- mean(sapply(temp[[2]][dist], function(x) x$vt$p_LRT$A) < 0.05)
  vt_gpcm <- mean(sapply(temp[[2]][dist], function(x) x$vt$p_LRT$B) < 0.05)
  c(vt_graded, vt_gpcm)
}))
rownames(lrt_v_ds) <- results11_oc$N
colnames(lrt_v_ds) <- c("vt_graded", "vt_gpcm")
results11_oc$vt_graded_ds <- lrt_v_ds[, 1]
results11_oc$vt_gpcm_ds <- lrt_v_ds[, 2]

### absolute differences in AIC values
temp <- readRDS("sim_graded_gpcm_base-results/results-row-1.rds")
cn <- which(sapply(temp[[2]], function(x) x$graded_check[1] == 1 & x$gpcm_check[1] == 1))
abs_aic11_1 <- sapply(temp[[2]], function(x) abs(x$AIC_graded - x$AIC_gpcm))[cn]

#setEPS()
#postscript("../figure1.eps", width = 4.5, height = 2.4, pointsize = 10)
#pdf("../figure1.pdf", width = 4.5, height = 2.4, pointsize = 10)

par(mai = c(0.45, 0.45, 0.1, 0.1), mgp = c(2.7, 0.7, 0), lwd = 1)
layout(matrix(c(1, 2, 2), 1, 3))
boxplot(abs_aic11_1, ylab = expression(paste("|", Delta["AIC"], "|")), ylim = c(0, 8))
hist(abs_aic11_1, ylab = "Frequency", xlab = expression(paste("|", Delta["AIC"], "|")), main = "",
  ylim = c(0, 400), col = "darkgrey")
box(bty = "L")

#dev.off()

### clean up
#rm(list = c("abs_aic11_1", "AIC_results", "cn", "i", "lrt_v_ds", "pth", "res11", "res11_oc",
#  "results11", "results11_oc", "temp"))



###################################################################################################
# Simulation 1.2: Code and results or rerun
###################################################################################################

if(length(grep("vuong_mirt_code/orig$", getwd()))) {
  cat("Working with the original simulation result files.\n")
} else if(length(grep("vuong_mirt_code/repl$", getwd()))) {
  if("sim_graded_gpcm.rds" %in% dir() | "sim_graded_gpcm-results" %in% dir()) {
    stop("New result files exist. You should inspect them prior to a rerun.\n")
  } else {
    cat("Rerunning the simulation. This may take a very long time.\n")
  }

  ### Settings
  N <- c(500, 1000, 2000) # persons
  M <- 10 # items
  W <- 0:10 # mismatching items; reference = graded
  K <- 4 # categories
  Design <- expand.grid(N = N, M = M, W = W, K = K)
  rm(N, M, W, K)

  ### SimDesign functions
  Generate <- function(condition, fixed_objects = NULL) {
    Attach(condition)
    theta <- as.matrix(rnorm(N, 0, 1))
    a <- exp(rnorm(M, 0, 0.25))
    D <- matrix(c(1, 0, -1), M, K - 1, byrow = TRUE)
    e <- rnorm(M, 0, 1)
    ## W = 0 --> graded, W = 10 --> gpcm, otherwise mixed itemtypes
    dat <- if(W > 0 && W < 10) {
      W_ind <- tail(1:M, W) # mismatching items
      d_graded <- (D + e)[-W_ind, ]
      d_gpcm <- cbind(0, D + e)[W_ind, ]
      ## force intercepts to always be a matrix
      if(class(d_graded) == "numeric") {
        d_graded <- t(as.matrix(d_graded))
      }
      if(class(d_gpcm) == "numeric") {
        d_gpcm <- t(as.matrix(d_gpcm))
      }
      dat_graded <- simdata(a = a[-W_ind], d = d_graded,
        itemtype = rep("graded", M - W), Theta = theta)
      dat_gpcm <- simdata(a = a[W_ind], d = d_gpcm,
        itemtype = rep("gpcm", W), Theta = theta)
      cbind(dat_graded, dat_gpcm)
    } else if(W == 0) {
      simdata(a = a, d = D + e, itemtype = rep("graded", M), Theta = theta)
    } else if(W == 10) {
      simdata(a = a, d = cbind(0, D + e), itemtype = rep("gpcm", M), Theta = theta)
    }
    colnames(dat) <- 1:M
    return(dat)
  }

  Analyse <- function(condition, dat, fixed_objects = NULL) {
    ## EM, BFGS, Ramsay, Oakes, 61 quadpts, 1e-4 TOL
    graded <- mirt(dat, 1, "graded", SE = TRUE, technical = list(NCYCLES = 5000))
    gpcm <- mirt(dat, 1, "gpcm", SE = TRUE, technical = list(NCYCLES = 5000))
    return(list(vt = vuongtest(graded, gpcm),
      AIC_graded = extract.mirt(graded, "AIC"),
      AIC_gpcm = extract.mirt(gpcm, "AIC"),
      BIC_graded = extract.mirt(graded, "BIC"),
      BIC_gpcm = extract.mirt(gpcm, "BIC"), M2_graded = M2(graded),
      M2_gpcm = M2(gpcm), graded_check = c(extract.mirt(graded, "converged"),
      extract.mirt(graded, "condnum"), extract.mirt(graded, "secondordertest")),
      gpcm_check = c(extract.mirt(gpcm, "converged"),
      extract.mirt(gpcm, "condnum"), extract.mirt(gpcm, "secondordertest"))))
  }

  Summarise <- function(condition, results, fixed_objects = NULL) {
    return(list(vt_dist = mean(sapply(results, function(x) x$vt$p_omega < 0.05)),
      vt_graded = mean(sapply(results, function(x) x$vt$p_LRT$A < 0.05)),
      vt_gpcm = mean(sapply(results, function(x) x$vt$p_LRT$B < 0.05)),
      AIC = mean(sapply(results, function(x) (x$AIC_graded - x$AIC_gpcm) < 0)),
      BIC = mean(sapply(results, function(x) (x$BIC_graded - x$BIC_gpcm) < 0)),
      M2_graded = mean(sapply(results, function(x) x$M2_graded$p < 0.05)),
      M2_gpcm = mean(sapply(results, function(x) x$M2_gpcm$p < 0.05)),
      graded_cn = mean(sapply(results, function(x) x$graded_check[1])),
      graded_so = mean(sapply(results, function(x) x$graded_check[3])),
      gpcm_cn = mean(sapply(results, function(x) x$gpcm_check[1])),
      gpcm_so = mean(sapply(results, function(x) x$gpcm_check[3]))))
  }

  ### Run
  results <- runSimulation(Design, replications = 1000,
    generate = Generate, analyse = Analyse, summarise = Summarise,
    packages = c("mirt", "nonnest2"), save = TRUE, save_results = TRUE,
    save_seeds = FALSE, parallel = TRUE, progress = TRUE, verbose = FALSE,
    filename = "sim_graded_gpcm",
    save_details = list(save_results_dirname = "sim_graded_gpcm-results"))

  ### clean up
  #rm(list = c("Analyse", "Design", "Generate", "results", "Summarise"))
} else {
  stop("Please set the working directory as described above.\n")
}



###################################################################################################
# Simulation 1.2: Figure 2 and text results
###################################################################################################

res12 <- readRDS("sim_graded_gpcm.rds")
#attributes(res12)

### check convergence
min(res12$graded_cn) # 1
min(res12$gpcm_cn)   # 1
min(res12$graded_so) # 1
min(res12$gpcm_so)   # 1
res12$AIC_graded <- res12$AIC
res12$AIC_gpcm <- 1 - res12$AIC
results12 <- res12[, c(1:14, 21, 22)]

### absolute differences in AIC values
temp <- readRDS("sim_graded_gpcm-results/results-row-1.rds")
abs_aic12_1 <- sapply(temp[[2]], function(x) abs(x$AIC_graded - x$AIC_gpcm))
temp <- readRDS("sim_graded_gpcm-results/results-row-10.rds")
abs_aic12_10 <- sapply(temp[[2]], function(x) abs(x$AIC_graded - x$AIC_gpcm))
temp <- readRDS("sim_graded_gpcm-results/results-row-16.rds")
abs_aic12_16 <- sapply(temp[[2]], function(x) abs(x$AIC_graded - x$AIC_gpcm))

#setEPS()
#postscript("../figure2.eps", width = 4.5, height = 2.4, pointsize = 10)
#pdf("../figure2.pdf", width = 4.5, height = 2.4, pointsize = 10)

par(mfrow = c(1, 3), oma = c(4.0, 3.4, 0, 3.8), mai = c(0, 0, 0.175, 0.05),
  mgp = c(2.7, 0.7, 0), lwd = 1.5)

temp <- results12[results12$N == 500, ]
temp$W <- as.factor(temp$W)
plot(as.numeric(temp$W), temp$vt_dist, type = "l", lty = 4, ylim = c(0, 1),
  axes = FALSE, xlab = "", ylab = "", main = expression(italic(N) == 500))
points(as.numeric(temp$W), temp$vt_graded, type = "l", lty = 1,
  col = "#FF8C00")
points(as.numeric(temp$W), temp$vt_gpcm, type = "l", lty = 1,
  col = "#8FCEF3")
points(as.numeric(temp$W), temp$AIC_graded, type = "l", lty = "92",
  col = "#FF8C00")
points(as.numeric(temp$W), temp$AIC_gpcm, type = "l", lty = "92",
  col = "#8FCEF3")
points(as.numeric(temp$W), temp$M2_graded, type = "l", lty = 3,
  col = "#FF8C00", lwd = 2)
points(as.numeric(temp$W), temp$M2_gpcm, type = "l", lty = 3,
  col = "#8FCEF3", lwd = 2)
axis(1, at = 1:11, labels = c(0, rep(NA, 4), 5, rep(NA, 4), 10))
axis(2, at = c(0, 0.05, 0.25, 0.5, 0.75, 1))
box(lwd = 1)
mtext(expression(italic(P)), 2, 2, FALSE, at = 0.5)

temp <- results12[results12$N == 1000, ]
temp$W <- as.factor(temp$W)
plot(as.numeric(temp$W), temp$vt_dist, type = "l", lty = 4, ylim = c(0, 1),
  axes = FALSE, xlab = "", ylab = "", main = expression(italic(N) == 1000))
points(as.numeric(temp$W), temp$vt_graded, type = "l", lty = 1,
  col = "#FF8C00")
points(as.numeric(temp$W), temp$vt_gpcm, type = "l", lty = 1,
  col = "#8FCEF3")
points(as.numeric(temp$W), temp$AIC_graded, type = "l", lty = "92",
  col = "#FF8C00")
points(as.numeric(temp$W), temp$AIC_gpcm, type = "l", lty = "92",
  col = "#8FCEF3")
points(as.numeric(temp$W), temp$M2_graded, type = "l", lty = 3,
  col = "#FF8C00", lwd = 2)
points(as.numeric(temp$W), temp$M2_gpcm, type = "l", lty = 3,
  col = "#8FCEF3", lwd = 2)
axis(1, at = 1:11, labels = c(0, rep(NA, 4), 5, rep(NA, 4), 10))
box(lwd = 1)

temp <- results12[results12$N == 2000, ]
temp$W <- as.factor(temp$W)
temp$W <- as.factor(temp$W)
plot(as.numeric(temp$W), temp$vt_dist, type = "l", lty = 4, ylim = c(0, 1),
  axes = FALSE, xlab = "", ylab = "", main = expression(italic(N) == 2000))
points(as.numeric(temp$W), temp$vt_graded, type = "l", lty = 1,
  col = "#FF8C00")
points(as.numeric(temp$W), temp$vt_gpcm, type = "l", lty = 1,
  col = "#8FCEF3")
points(as.numeric(temp$W), temp$AIC_graded, type = "l", lty = "92",
  col = "#FF8C00")
points(as.numeric(temp$W), temp$AIC_gpcm, type = "l", lty = "92",
  col = "#8FCEF3")
points(as.numeric(temp$W), temp$M2_graded, type = "l", lty = 3,
  col = "#FF8C00", lwd = 2)
points(as.numeric(temp$W), temp$M2_gpcm, type = "l", lty = 3,
  col = "#8FCEF3", lwd = 2)
axis(1, at = 1:11, labels = c(0, rep(NA, 4), 5, rep(NA, 4), 10))
box(lwd = 1)

mtext(expression(italic(D)), 1, 2, TRUE, 0.495)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend(0.92, 0.971, legend = "Dist", lty = 4, bty = "n", cex = 0.8,
  seg.len = 3)
legend(0.92, 0.941, legend = "AIC", lty = "92", bty = "n", col = "#8FCEF3",
  text.col = "#8FCEF3", cex = 0.8, seg.len = 3)
legend(0.92, 0.911, legend = expression(LRT[v]), lty = 1, bty = "n",
  col = "#8FCEF3", text.col = "#8FCEF3", cex = 0.8, seg.len = 3)
legend(0.92, 0.358, legend = expression(italic(M[2] ^ "*")), lty = 3,
  bty = "n", col = "#8FCEF3", text.col = "#8FCEF3", cex = 0.8, lwd = 2,
  seg.len = 3)
legend(0.92, 0.308, legend = expression(italic(M[2] ^ "*")), lty = 3,
  bty = "n", col = "#FF8C00", text.col = "#FF8C00", cex = 0.8, lwd = 2,
  seg.len = 3)
legend(0.92, 0.258, legend = "AIC", lty = "92", bty = "n", col = "#FF8C00",
  text.col = "#FF8C00", cex = 0.8, seg.len = 3)
legend(0.92, 0.228, legend = expression(LRT[v]), lty = 1, bty = "n",
  col = "#FF8C00", text.col = "#FF8C00", cex = 0.8, seg.len = 3)
legend(0.4, 0.025, c("GRM", "GPCM"), bty = "n", pch = c(19, 19),
  ncol = 2, cex = 0.8, col = c("#FF8C00", "#8FCEF3"))

#dev.off()

### clean up
#rm(list = c("abs_aic12_1", "abs_aic12_10", "abs_aic12_16", "res12", "results12", "temp"))



###################################################################################################
# Simulation 2.1: Code and results or rerun
###################################################################################################

if(length(grep("vuong_mirt_code/orig$", getwd()))) {
  cat("Working with the original simulation result files.\n")
} else if(length(grep("vuong_mirt_code/repl$", getwd()))) {
  if("sim_rm_2pl_mis.rds" %in% dir() | "sim_rm_2pl_mis-results" %in% dir()) {
    stop("New result files exist. You should inspect them prior to a rerun.\n")
  } else {
    cat("Rerunning the simulation. This may take a very long time.\n")
  }

  ### Settings
  N <- c(500, 1000, 2000) # persons
  M <- c(10, 20, 30, 40) # items
  g <- c("RM", "2PL", 0.01, 0.05, 0.25) # guessing parameter
  Design <- expand.grid(N = N, M = M, g = g)
  rm(N, M, g)

  ### SimDesign functions
  Generate <- function(condition, fixed_objects = NULL) {
    Attach(condition)
    a <- if(g != "2PL") {
      g <- if(g == "RM") {
        0
      } else if(g != "RM") {
        as.numeric(as.character(g))
      }
      rep(1, M)
    } else if(g == "2PL") {
      g <- 0
      exp(rnorm(M, 0, 0.25))
    }
    return(simdata(a = a, d = rnorm(M, 0, 1), guess = rep(g, M),
      itemtype = rep("3PL", M), Theta = as.matrix(rnorm(N, 0, 1))))
  }

  Analyse <- function(condition, dat, fixed_objects = NULL) {
    ## EM, BFGS, Ramsay, Oakes, 61 quadpts, 1e-4 TOL
    m1 <- mirt(dat, 1, "Rasch", SE = TRUE, technical = list(NCYCLES = 5000))
    m2 <- mirt(dat, 1, "2PL", SE = TRUE, technical = list(NCYCLES = 5000))
    return(list(vt = vuongtest(m1, m2, nested = TRUE), lrt = anova(m1, m2),
      M2_1 = M2(m1), M2_2 = M2(m2), m1_check = c(extract.mirt(m1, "converged"),
      extract.mirt(m1, "condnum"), extract.mirt(m1, "secondordertest")),
      m2_check = c(extract.mirt(m2, "converged"), extract.mirt(m2, "condnum"),
      extract.mirt(m2, "secondordertest"))))
  }

  Summarise <- function(condition, results, fixed_objects = NULL) {
    return(list(vt_dist = mean(sapply(results, function(x) x$vt$p_omega < 0.05)),
      vt_lrt = mean(sapply(results, function(x) x$vt$p_LRT$A < 0.05)),
      lrt = mean(sapply(results, function(x) x$lrt$p[2] < 0.05)),
      AIC = mean(sapply(results, function(x) (x$lrt$AIC[2] - x$lrt$AIC[1]) < 0)),
      BIC = mean(sapply(results, function(x) (x$lrt$BIC[2] - x$lrt$BIC[1]) < 0)),
      M2_1 = mean(sapply(results, function(x) x$M2_1$p < 0.05)),
      M2_2 = mean(sapply(results, function(x) x$M2_2$p < 0.05)),
      m1_cn = mean(sapply(results, function(x) x$m1_check[1])),
      m1_so = mean(sapply(results, function(x) x$m1_check[3])),
      m2_cn = mean(sapply(results, function(x) x$m2_check[1])),
      m2_so = mean(sapply(results, function(x) x$m2_check[3]))))
  }

  ### Run
  results <- runSimulation(Design, replications = 1000,
    generate = Generate, analyse = Analyse, summarise = Summarise,
    packages = c("mirt", "nonnest2"), save = TRUE, save_results = TRUE,
    save_seeds = FALSE, parallel = TRUE, progress = TRUE, verbose = FALSE,
    filename = "sim_rm_2pl_mis",
    save_details = list(save_results_dirname = "sim_rm_2pl_mis-results"))

  ### clean up
  #rm(list = c("Analyse", "Design", "Generate", "results", "Summarise"))
} else {
  stop("Please set the working directory as described above.\n")
}



###################################################################################################
# Simulation 2.1: Table 2, Figure 3 and text results
###################################################################################################

res21 <- readRDS("sim_rm_2pl_mis.rds")

### check convergence
min(res21$m1_cn) # 1
min(res21$m2_cn) # 1
min(res21$m1_so) # 1
min(res21$m2_so) # 1
results21 <- res21[, 1:14]
results21_1 <- results21[results21$N == 500, ]
results21_2 <- results21[results21$N == 1000, ]
results21_3 <- results21[results21$N == 2000, ]

temp <- readRDS("sim_rm_2pl_mis-results/results-row-3.rds")

#setEPS()
#postscript("../figure3.eps", width = 6, height = 2.4, pointsize = 10)
#pdf("../figure3.pdf", width = 6, height = 2.4, pointsize = 10)

par(mai = c(0.45, 0.45, 0.1, 0.1), mgp = c(2.7, 0.7, 0), lwd = 1)
layout(matrix(c(1, 1, 2, 2), 1, 4))
hist(sapply(temp[[2]], function(x) x$vt$p_omega), main = "", xlab = expression(italic(p)[Dist]),
  col = "darkgrey")
box(bty = "L")
hist(sapply(temp[[2]], function(x) x$vt$p_LRT$A), main = "", xlab = expression(italic(p)[LRT[v]]),
  col = "darkgrey")
box(bty = "L")

#dev.off()

### clean up
#rm(list = c("res21", "results21", "results21_1", "results21_2", "results21_3", "temp"))



###################################################################################################
# Simulation 2.2: Code and results or rerun
###################################################################################################

if(length(grep("vuong_mirt_code/orig$", getwd()))) {
  cat("Working with the original simulation result files.\n")
} else if(length(grep("vuong_mirt_code/repl$", getwd()))) {
  if("sim_2pl_2dim.rds" %in% dir() | "sim_2pl_2dim-results" %in% dir()) {
    stop("New result files exist. You should inspect them prior to a rerun.\n")
  } else {
    cat("Rerunning the simulation. This may take a very long time.\n")
  }

  ### Settings
  N <- 2000 # persons
  M <- c(10, 20, 30, 40) # items
  rho <- c(1, 2 / 3, 1 / 3, 0) # correlation of the dimensions
  Design <- expand.grid(N = N, M = M, rho = rho)
  rm(N, M, rho)

  ### SimDesign functions
  Generate <- function(condition, fixed_objects = NULL) {
    Attach(condition)
    ## rho == 1 unidimensional 2pl
    if(rho == 1) {
      a <- exp(rnorm(M, 0, 0.25))
      theta <- as.matrix(rnorm(N, 0, 1))
    } else {
      a <- cbind(exp(rnorm(M, 0, 0.25)), exp(rnorm(M, 0, 0.25)))
      theta <- mvrnorm(N, c(0, 0), matrix(c(1, rho, rho, 1), 2, 2))
    }
    return(simdata(a = a, d = rnorm(M, 0, 1), itemtype = rep("2PL", M),
      Theta = theta))
  }

  Analyse <- function(condition, dat, fixed_objects = NULL) {
    ## EM, BFGS, Ramsay, Oakes, 61/31 quadpts, 1e-4 TOL, (last item a2 = 0)
    m1 <- mirt(dat, 1, "2PL", SE = TRUE, technical = list(NCYCLES = 5000))
    m2 <- mirt(dat, 2, "2PL", SE = TRUE, technical = list(NCYCLES = 5000))
    return(list(vt = vuongtest(m1, m2, nested = TRUE), lrt = anova(m1, m2),
      M2_1 = M2(m1), M2_2 = M2(m2), m1_check = c(extract.mirt(m1, "converged"),
      extract.mirt(m1, "condnum"), extract.mirt(m1, "secondordertest")),
      m2_check = c(extract.mirt(m2, "converged"), extract.mirt(m2, "condnum"),
      extract.mirt(m2, "secondordertest"))))
  }

  Summarise <- function(condition, results, fixed_objects = NULL) {
    return(list(vt_dist = mean(sapply(results, function(x) x$vt$p_omega < 0.05)),
      vt_lrt = mean(sapply(results, function(x) x$vt$p_LRT$A < 0.05)),
      lrt = mean(sapply(results, function(x) x$lrt$p[2] < 0.05)),
      AIC = mean(sapply(results, function(x) (x$lrt$AIC[2] - x$lrt$AIC[1]) < 0)),
      BIC = mean(sapply(results, function(x) (x$lrt$BIC[2] - x$lrt$BIC[1]) < 0)),
      M2_1 = mean(sapply(results, function(x) x$M2_1$p < 0.05)),
      M2_2 = mean(sapply(results, function(x) x$M2_2$p < 0.05)),
      m1_cn = mean(sapply(results, function(x) x$m1_check[1])),
      m1_so = mean(sapply(results, function(x) x$m1_check[3])),
      m2_cn = mean(sapply(results, function(x) x$m2_check[1])),
      m2_so = mean(sapply(results, function(x) x$m2_check[3]))))
  }

  ### Run
  results <- runSimulation(Design, replications = 1000,
    generate = Generate, analyse = Analyse, summarise = Summarise,
    packages = c("MASS", "mirt", "nonnest2"), save = TRUE, save_results = TRUE,
    save_seeds = FALSE, parallel = TRUE, progress = TRUE, verbose = FALSE,
    filename = "sim_2pl_2dim",
    save_details = list(save_results_dirname = "sim_2pl_2dim-results"))

  ### clean up
  #rm(list = c("Analyse", "Design", "Generate", "results", "Summarise"))
} else {
  stop("Please set the working directory as described above.\n")
}



###################################################################################################
# Simulation 2.2: Figure 4 and text results
###################################################################################################

res22 <- readRDS("sim_2pl_2dim.rds")

### check convergence
min(res22$m1_cn) # 1
min(res22$m2_cn) # 0.959
min(res22$m1_so) # 1
min(res22$m2_so) # 1
results22 <- res22[, 1:14]

### results of replications in which both models converged
res22_oc <- res22
res22_oc$cn_both <- NULL
for(i in 1:dim(res22_oc)[1]) {
  pth <- gsub("-i", paste0("-", i), "sim_2pl_2dim-results/results-row-i.rds")
  temp <- readRDS(pth)
  cn <- which(sapply(temp[[2]], function(x) x$m1_check[1] == 1 & x$m2_check[1] == 1))
  res22_oc$cn_both[i] <- length(cn) / length(temp[[2]])
  ## vt_dist
  res22_oc[i, 4] <- mean(sapply(temp[[2]], function(x) x$vt$p_omega)[cn] < 0.05)
  ## vt_lrt
  res22_oc[i, 5] <- mean(sapply(temp[[2]], function(x) x$vt$p_LRT$A)[cn] < 0.05)
  ## lrt
  res22_oc[i, 6] <- mean(sapply(temp[[2]], function(x) x$lrt$p[2])[cn] < 0.05)
  ## AIC
  res22_oc[i, 7] <- mean(sapply(temp[[2]], function(x) (x$lrt$AIC[2] - x$lrt$AIC[1]) < 0)[cn])
  ## BIC
  res22_oc[i, 8] <- mean(sapply(temp[[2]], function(x) (x$lrt$BIC[2] - x$lrt$BIC[1]) < 0)[cn])
  ## M2_1
  res22_oc[i, 9] <- mean(sapply(temp[[2]], function(x) x$M2_1$p)[cn] < 0.05)
  ## M2_1
  res22_oc[i, 10] <- mean(sapply(temp[[2]], function(x) x$M2_2$p)[cn] < 0.05)
  ## m1_cn
  res22_oc[i, 11] <- mean(sapply(temp[[2]], function(x) x$m1_check[1] == 1)[cn])
  ## m1_so
  res22_oc[i, 12] <- mean(sapply(temp[[2]], function(x) x$m1_check[3] == 1)[cn])
  ## m2_cn
  res22_oc[i, 13] <- mean(sapply(temp[[2]], function(x) x$m2_check[1] == 1)[cn])
  ## m2_so
  res22_oc[i, 14] <- mean(sapply(temp[[2]], function(x) x$m2_check[3] == 1)[cn])
}
results22_oc <- res22_oc[, 1:14]
results22_oc <- cbind(results22_oc, "cn_both" = res22_oc$cn_both)
which.min(results22_oc$cn_both) # 6 N = 2000, M = 20, rho = 2 / 3

#setEPS()
#postscript("../figure4.eps", width = 6, height = 2.4, pointsize = 10)
#pdf("../figure4.pdf", width = 6, height = 2.4, pointsize = 10)

par(mfrow = c(1, 4), oma = c(4.0, 3.4, 0, 3.8), mai = c(0, 0, 0.175, 0.05),
  mgp = c(2.7, 0.7, 0), lwd = 1.5)

temp <- results22_oc[results22_oc$rho == 1, ]
temp$M <- as.factor(temp$M)
plot(as.numeric(temp$M), temp$vt_dist, type = "l", lty = 4, ylim = c(0, 1),
  axes = FALSE, xlab = "", ylab = "", main = expression("2PL"))
points(as.numeric(temp$M), temp$vt_lrt, type = "l", lty = 1)
points(as.numeric(temp$M), temp$lrt, type = "l", lty = 6, col = "#D5D5D5")
points(as.numeric(temp$M), temp$AIC, type = "l", lty = "92", col = "#7D7D7D")
points(as.numeric(temp$M), temp$BIC, type = "l", lty = 2, col = "#7D7D7D")
points(as.numeric(temp$M), temp$M2_1, type = "l", lty = 3,
  col = "#FF8C00", lwd = 2)
points(as.numeric(temp$M), temp$M2_2, type = "l", lty = 3,
  col = "#8FCEF3", lwd = 2)
axis(1, at = 1:4, labels = c(10, 20, 30, 40))
axis(2, at = c(0, 0.05, 0.25, 0.5, 0.75, 1))
box(lwd = 1)
mtext(expression(italic(P)), 2, 2, FALSE, at = 0.5)

temp <- results22_oc[results22_oc$rho == 2 / 3, ]
temp$M <- as.factor(temp$M)
plot(as.numeric(temp$M), temp$vt_dist, type = "l", lty = 4, ylim = c(0, 1),
  axes = FALSE, xlab = "", ylab = "", main = expression(italic(rho) == 2 / 3))
points(as.numeric(temp$M), temp$vt_lrt, type = "l", lty = 1)
points(as.numeric(temp$M), temp$lrt, type = "l", lty = 6, col = "#D5D5D5")
points(as.numeric(temp$M), temp$AIC, type = "l", lty = "92", col = "#7D7D7D")
points(as.numeric(temp$M), temp$BIC, type = "l", lty = 2, col = "#7D7D7D")
points(as.numeric(temp$M), temp$M2_1, type = "l", lty = 3,
  col = "#FF8C00", lwd = 2)
points(as.numeric(temp$M), temp$M2_2, type = "l", lty = 3,
  col = "#8FCEF3", lwd = 2)
axis(1, at = 1:4, labels = c(10, 20, 30, 40))
box(lwd = 1)

temp <- results22_oc[results22_oc$rho == 1 / 3, ]
temp$M <- as.factor(temp$M)
plot(as.numeric(temp$M), temp$vt_dist, type = "l", lty = 4, ylim = c(0, 1),
  axes = FALSE, xlab = "", ylab = "", main = expression(italic(rho) == 1 / 3))
points(as.numeric(temp$M), temp$vt_lrt, type = "l", lty = 1)
points(as.numeric(temp$M), temp$lrt, type = "l", lty = 6, col = "#D5D5D5")
points(as.numeric(temp$M), temp$AIC, type = "l", lty = "92", col = "#7D7D7D")
points(as.numeric(temp$M), temp$BIC, type = "l", lty = 2, col = "#7D7D7D")
points(as.numeric(temp$M), temp$M2_1, type = "l", lty = 3,
  col = "#FF8C00", lwd = 2)
points(as.numeric(temp$M), temp$M2_2, type = "l", lty = 3,
  col = "#8FCEF3", lwd = 2)
axis(1, at = 1:4, labels = c(10, 20, 30, 40))
box(lwd = 1)

temp <- results22_oc[results22_oc$rho == 0, ]
temp$M <- as.factor(temp$M)
plot(as.numeric(temp$M), temp$vt_dist, type = "l", lty = 4, ylim = c(0, 1),
  axes = FALSE, xlab = "", ylab = "", main = expression(italic(rho) == 0))
points(as.numeric(temp$M), temp$vt_lrt, type = "l", lty = 1)
points(as.numeric(temp$M), temp$lrt, type = "l", lty = 6, col = "#D5D5D5")
points(as.numeric(temp$M), temp$AIC, type = "l", lty = "92", col = "#7D7D7D")
points(as.numeric(temp$M), temp$BIC, type = "l", lty = 2, col = "#7D7D7D")
points(as.numeric(temp$M), temp$M2_1, type = "l", lty = 3,
  col = "#FF8C00", lwd = 2)
points(as.numeric(temp$M), temp$M2_2, type = "l", lty = 3,
  col = "#8FCEF3", lwd = 2)
axis(1, at = 1:4, labels = c(10, 20, 30, 40))
box(lwd = 1)

mtext(expression(italic(J)), 1, 2, TRUE, 0.495)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend(0.95, 0.982, legend = expression(LRT[t]), lty = 6, bty = "n",
  col = "#D5D5D5", text.col = "#D5D5D5", cex = 0.8, seg.len = 3)
legend(0.95, 0.950, legend = "AIC", lty = "92", bty = "n", col = "#7D7D7D",
  text.col = "#7D7D7D", cex = 0.8, seg.len = 3)
legend(0.95, 0.920, legend = expression(LRT[v]), lty = 1, bty = "n",
  cex = 0.8, seg.len = 3)
legend(0.95, 0.880, legend = expression(italic(M[2])), lty = 3,
  bty = "n", col = "#FF8C00", text.col = "#FF8C00", cex = 0.8, lwd = 2,
  seg.len = 3)
legend(0.95, 0.810, legend = "Dist", lty = 4, bty = "n",
  cex = 0.8, seg.len = 3)
legend(0.95, 0.283, legend = expression(italic(M[2])), lty = 3,
  bty = "n", col = "#8FCEF3", text.col = "#8FCEF3", cex = 0.8, lwd = 2,
  seg.len = 3)
legend(0.95, 0.243, legend = "BIC", lty = 2, bty = "n", col = "#7D7D7D",
  text.col = "#7D7D7D", cex = 0.8, seg.len = 3)
legend(0.42, 0.025, c("2PL", "2d-2PL"), bty = "n", pch = c(19, 19),
  ncol = 2, cex = 0.8, col = c("#FF8C00", "#8FCEF3"))

#dev.off()

### clean up
#rm(list = c("cn", "i", "pth", "res22", "res22_oc", "results22", "results22_oc", "temp"))



###################################################################################################
# Application: The Nerdy Personality Attributes Scale
###################################################################################################

require(tools) # for checking the md5sum

### data from the Open Source Psychometrics Project
temp <- tempfile()
check <- try(
  download.file("http://openpsychometrics.org/_rawdata/NPAS-data.zip", temp),
  silent = TRUE
)
if(check != 0 | md5sum(temp) != "1f9e96af145de42599477070e8a4fafb") {
  stop("It seems like the dataset has changed. Please open an issue.")
  browseURL("https://github.com/sumny/vuong_mirt_code/issues")
} else {
  raw <- read.csv(unz(temp, "NPAS-data/NPAS-data.csv"))
}
unlink(temp)

### get the three validity check items
temp <- raw[, c(45, 48, 51)]
temp[is.na(temp)] <- 0
### exclude based on the validity check items
check <- which(temp$VCL6 != 0 | temp$VCL9 != 0 | temp$VCL12 != 0)
### subset to the six science related items
dat <- raw[-check, c(1, 2, 6, 13, 22, 23)]
### handle response codes, NA and shift data to zero
dat[dat == 0] <- NA
dat[dat == 6] <- NA
dat[dat == 7] <- NA
full_na <- which(!apply(dat, 1, function(x) any(!(is.na(x)))))
dat <- dat[- full_na, ] - 1

### graded vs. gpcm
graded <- mirt(dat, model = 1, itemtype = "graded", SE = TRUE,
  technical = list(NCYCLES = 5000), verbose = FALSE)
gpcm <- mirt(dat, model = 1, itemtype = "gpcm", SE = TRUE,
  technical = list(NCYCLES = 5000), verbose = FALSE)
M2(graded, na.rm = TRUE, type = "C2")
M2(gpcm, na.rm = TRUE, type = "C2")
extract.mirt(graded, "AIC")
extract.mirt(gpcm, "AIC")
vuongtest(graded, gpcm, nested = FALSE)

### 1d graded vs. 2d graded
graded2 <- mirt(dat, model = 2, itemtype = "graded", SE = TRUE,
  technical = list(NCYCLES = 5000), verbose = FALSE)
M2(graded2, na.rm = TRUE, type = "C2")
extract.mirt(graded2, "AIC")
anova(graded, graded2)
vuongtest(graded, graded2, nested = TRUE)

### clean up
#rm(list = c("check", "dat", "full_na", "gpcm", "graded", "graded2", "raw", "temp"))


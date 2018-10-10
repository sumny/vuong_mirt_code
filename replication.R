###############################################################################
# Packages and versions, original results or rerun the simulations?
###############################################################################

### for exact results use the versions specified below

# R version 3.4.3
require(MASS) # version 7.3.48
require(mirt) # version 1.26.3
require(nonnest2) # version 0.5
require(SimDesign) # version 1.7

### work with the original results or rerun the simulations?
#setwd("~/vuong_mirt_code/orig") # work with the original results
#setwd("~/vuong_mirt_code/repl") # rerun the simulations

# Last mod: October/10/2018, LS (non-content modifications)



###############################################################################
# Simulation 1: Code and results or rerun
###############################################################################

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
    if(W > 0 & W < 10) {
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
      dat <- simdata(a = a[-W_ind], d = d_graded,
        itemtype = rep("graded", M - W), Theta = theta)
      dat <- cbind(dat, simdata(a = a[W_ind], d = d_gpcm,
        itemtype = rep("gpcm", W), Theta = theta))
    } else if (W == 0) {
      dat <- simdata(a = a, d = D + e,
        itemtype = rep("graded", M), Theta = theta)
    } else {
      dat <- simdata(a = a, d = cbind(0, D + e),
        itemtype = rep("gpcm", M), Theta = theta)
    }
    colnames(dat) <- paste("I", 1:M, sep = "")
    return(dat)
  }

  Analyse <- function(condition, dat, fixed_objects = NULL) {
    graded <- mirt(dat, 1, "graded", SE = TRUE)
    gpcm <- mirt(dat, 1, "gpcm", SE = TRUE)
    ret <- list(vt = vuongtest(graded, gpcm),
      aic_graded = extract.mirt(graded, "AIC"),
      aic_gpcm = extract.mirt(gpcm, "AIC"),
      M2_graded = M2(graded), M2_gpcm = M2(gpcm))
    return(ret)
  }

  Summarise <- function(condition, results, fixed_objects = NULL) {
    vt_dist <- mean(sapply(results, function(x) x$vt$p_omega < 0.05))
    vt_graded <- mean(sapply(results, function(x) x$vt$p_LRT$A < 0.05))
    vt_gpcm <- mean(sapply(results, function(x) x$vt$p_LRT$B < 0.05))
    AIC <- mean(sapply(results, function(x) (x$aic_graded - x$aic_gpcm) < 0))
    M2_graded <- mean(sapply(results, function(x) x$M2_graded$p < 0.05))
    M2_gpcm <- mean(sapply(results, function(x) x$M2_gpcm$p < 0.05))
    return(list(vt_dist = vt_dist, vt_graded = vt_graded, vt_gpcm = vt_gpcm,
      AIC = AIC, M2_graded = M2_graded, M2_gpcm = M2_gpcm))
  }

  ### Run
  results <- runSimulation(Design, replications = 1000, generate = Generate,
    analyse = Analyse, summarise = Summarise, packages = c("mirt", "nonnest2"),
    save = TRUE, save_results = TRUE, parallel = TRUE, progress = TRUE,
    verbose = FALSE, filename = "sim_graded_gpcm",
    save_details = list(save_results_dirname = "sim_graded_gpcm-results"))

  ### clean up
  #rm(list = c("Analyse", "Design",  "Generate", "results", "Summarise"))
} else {
  stop("Please set the working directory as described above.\n")
}



###############################################################################
# Simulation 1: Figure 1 and text results
###############################################################################

### the figure will only look nice if exported as .eps or .pdf

res1 <- readRDS("sim_graded_gpcm.rds")
#attributes(res1)

res1$AIC_graded <- res1$AIC
res1$AIC_gpcm <- 1 - res1$AIC

#setEPS()
#postscript("../figure1.eps", width = 4.5, height = 2.4, pointsize = 10)
#pdf("../figure1.pdf", width = 4.5, height = 2.4, pointsize = 10)

par(mfrow = c(1, 3), oma = c(4.0, 3.4, 0, 3.8), mai = c(0, 0, 0.175, 0.05),
  mgp = c(2.7, 0.7, 0), lwd = 1.5)

temp <- res1[res1$N == 500, ]
temp$W <- as.factor(temp$W)
plot(as.numeric(temp$W), temp$vt_dist, type = "l", lty = 4, ylim = c(0, 1),
  axes = FALSE, xlab = "", ylab = "", main = expression(italic(N) == 500))
points(as.numeric(temp$W), temp$vt_graded, type = "l", lty = 1,
  col = "#FF8C00")
points(as.numeric(temp$W), temp$vt_gpcm, type = "l", lty = 1,
  col = "#8FCEF3")
points(as.numeric(temp$W), temp$AIC_graded, type = "l", lty = 2,
  col = "#FF8C00")
points(as.numeric(temp$W), temp$AIC_gpcm, type = "l", lty = 2,
  col = "#8FCEF3")
points(as.numeric(temp$W), temp$M2_graded, type = "l", lty = 3,
  col = "#FF8C00", lwd = 2)
points(as.numeric(temp$W), temp$M2_gpcm, type = "l", lty = 3,
  col = "#8FCEF3", lwd = 2)
axis(1, at = 1:11, labels = c(0, rep(NA, 4), 5, rep(NA, 4), 10))
axis(2, at = c(0, 0.05, 0.25, 0.5, 0.75, 1))
box(lwd = 1)
mtext(expression(italic(P)), 2, 2, FALSE, at = 0.5)

temp <- res1[res1$N == 1000, ]
temp$W <- as.factor(temp$W)
plot(as.numeric(temp$W), temp$vt_dist, type = "l", lty = 4, ylim = c(0, 1),
  axes = FALSE, xlab = "", ylab = "", main = expression(italic(N) == 1000))
points(as.numeric(temp$W), temp$vt_graded, type = "l", lty = 1,
  col = "#FF8C00")
points(as.numeric(temp$W), temp$vt_gpcm, type = "l", lty = 1,
  col = "#8FCEF3")
points(as.numeric(temp$W), temp$AIC_graded, type = "l", lty = 2,
  col = "#FF8C00")
points(as.numeric(temp$W), temp$AIC_gpcm, type = "l", lty = 2,
  col = "#8FCEF3")
points(as.numeric(temp$W), temp$M2_graded, type = "l", lty = 3,
  col = "#FF8C00", lwd = 2)
points(as.numeric(temp$W), temp$M2_gpcm, type = "l", lty = 3,
  col = "#8FCEF3", lwd = 2)
axis(1, at = 1:11, labels = c(0, rep(NA, 4), 5, rep(NA, 4), 10))
box(lwd = 1)

temp <- res1[res1$N == 2000, ]
temp$W <- as.factor(temp$W)
temp$W <- as.factor(temp$W)
plot(as.numeric(temp$W), temp$vt_dist, type = "l", lty = 4, ylim = c(0, 1),
  axes = FALSE, xlab = "", ylab = "", main = expression(italic(N) == 2000))
points(as.numeric(temp$W), temp$vt_graded, type = "l", lty = 1,
  col = "#FF8C00")
points(as.numeric(temp$W), temp$vt_gpcm, type = "l", lty = 1,
  col = "#8FCEF3")
points(as.numeric(temp$W), temp$AIC_graded, type = "l", lty = 2,
  col = "#FF8C00")
points(as.numeric(temp$W), temp$AIC_gpcm, type = "l", lty = 2,
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
  seg.len = 2.5)
legend(0.92, 0.941, legend = "AIC", lty = 2, bty = "n", col = "#8FCEF3",
  text.col = "#8FCEF3", cex = 0.8, seg.len = 2.5)
legend(0.92, 0.911, legend = expression(LRT[v]), lty = 1, bty = "n",
  col = "#8FCEF3", text.col = "#8FCEF3", cex = 0.8, seg.len = 2.5)
legend(0.92, 0.358, legend = expression(italic(M[2] ^ "*")), lty = 3,
  bty = "n", col = "#8FCEF3", text.col = "#8FCEF3", cex = 0.8, lwd = 2,
  seg.len = 2.5)
legend(0.92, 0.308, legend = expression(italic(M[2] ^ "*")), lty = 3,
  bty = "n", col = "#FF8C00", text.col = "#FF8C00", cex = 0.8, lwd = 2,
  seg.len = 2.5)
legend(0.92, 0.258, legend = "AIC", lty = 2, bty = "n", col = "#FF8C00",
  text.col = "#FF8C00", cex = 0.8, seg.len = 2.5)
legend(0.92, 0.228, legend = expression(LRT[v]), lty = 1, bty = "n",
  col = "#FF8C00", text.col = "#FF8C00", cex = 0.8, seg.len = 2.5)
legend(0.4, 0.025, c("GRM", "GPCM"), bty = "n", pch = c(19, 19),
  ncol = 2, cex = 0.8, col = c("#FF8C00", "#8FCEF3"))

#dev.off()

### clean up
#rm(list = c("res1", "temp"))



###############################################################################
# Simulation 2: Code and results or rerun
###############################################################################

if(length(grep("vuong_mirt_code/orig$", getwd()))) {
  cat("Working with the original simulation result files.\n")
} else if(length(grep("vuong_mirt_code/repl$", getwd()))) {
  if("sim_2pl_2dim" %in% dir() | "sim_2pl_2dim-results" %in% dir()) {
    stop("New result files exist. You should inspect them prior to a rerun.\n")
  } else {
    cat("Rerunning the simulation. This may take a very long time.\n")
  }
  ### Settings
  N <- c(2000) # persons
  M <- c(10, 20, 30, 40) # items
  rho <- c(1, 2 / 3, 1 / 3, 0) # correlation of the dimensions
  Design <- expand.grid(N = N, M = M, rho = rho)
  rm(N, M, rho)

  ### SimDesign functions
  Generate <- function(condition, fixed_objects = NULL) {
    Attach(condition)
    mu <- c(0, 0)
    sigma <- matrix(c(1, rho, rho, 1), 2, 2)
    theta <- mvrnorm(N, mu, sigma)
    a <- cbind(exp(rnorm(M, 0, 0.25)), exp(rnorm(M, 0, 0.25)))
    d <- rnorm(M, 0, 1)
    dat <- simdata(a = a, d = d, itemtype = rep("2PL", M), Theta = theta)
    return(dat)
  }

  Analyse <- function(condition, dat, fixed_objects = NULL) {
    m1 <- mirt(dat, 1, "2PL", SE = TRUE)
    m2 <- mirt(dat, 2, "2PL", SE = TRUE)
    ret <- list(vt = vuongtest(m1, m2, nested = TRUE), lrt = anova(m1, m2),
      M2_1 = M2(m1), M2_2 = M2(m2))
    return(ret)
  }

  Summarise <- function(condition, results, fixed_objects = NULL) {
    vt_dist <- mean(sapply(results, function(x) x$vt$p_omega < 0.05))
    vt_lrt <- mean(sapply(results, function(x) x$vt$p_LRT$A < 0.05))
    lrt <- mean(sapply(results, function(x) x$lrt$p[2] < 0.05))
    AIC <- mean(sapply(results, function(x) (x$lrt$AIC[2] - x$lrt$AIC[1]) < 0))
    M2_1 <- mean(sapply(results, function(x) x$M2_1$p < 0.05))
    M2_2 <- mean(sapply(results, function(x) x$M2_2$p < 0.05))
    return(list(vt_dist = vt_dist, vt_lrt = vt_lrt, lrt = lrt, AIC = AIC,
      M2_1 = M2_1, M2_2 = M2_2))
  }

  ### Run
  results <- runSimulation(Design, replications = 1000,
    generate = Generate, analyse = Analyse, summarise = Summarise,
    packages = c("MASS", "mirt", "nonnest2"), save = TRUE, save_results = TRUE,
    parallel = TRUE, progress = TRUE, verbose = FALSE,
    filename = "sim_2pl_2dim",
    save_details = list(save_results_dirname = "sim_2pl_2dim-results"))

  ### clean up
  #rm(list = c("Analyse", "Design",  "Generate", "results", "Summarise"))
} else {
  stop("Please set the working directory as described above.\n")
}



###############################################################################
# Simulation 2: Figure 2 and text results
###############################################################################

### the figure will only look nice if exported as .eps or .pdf

res2 <- readRDS("sim_2pl_2dim.rds")
#attributes(res2)
path <- "sim_2pl_2dim-results"

### get BIC
I <- length(dir(path))
res2$BIC <- numeric(I)
for(i in seq_len(I)) {
  temp <- readRDS(gsub("I", i, paste0(path, "/results-row-I.rds")))
  res2$BIC[i] <- sum(sapply(temp$results, function(x) {
    (x[[2]][1, 4] - x[[2]][2, 4]) > 0
    }) / unique(res2$REPLICATIONS)
  )
}

#setEPS()
#postscript("../figure2.eps", width = 6, height = 2.4, pointsize = 10)
#pdf("../figure2.pdf", width = 6, height = 2.4, pointsize = 10)

par(mfrow = c(1, 4), oma = c(4.0, 3.4, 0, 3.8), mai = c(0, 0, 0.175, 0.05),
  mgp = c(2.7, 0.7, 0), lwd = 1.5)

temp <- res2[res2$N == 2000 & res2$rho == 1, ]
temp$M <- as.factor(temp$M)
plot(as.numeric(temp$M), temp$vt_dist, type = "l", lty = 4, ylim = c(0, 1),
  axes = FALSE, xlab = "", ylab = "",
  main = expression(italic(rho) == 1~("2PL")))
points(as.numeric(temp$M), temp$vt_lrt, type = "l", lty = 1)
points(as.numeric(temp$M), temp$lrt, type = "l", lty = 6, col = "#D5D5D5")
points(as.numeric(temp$M), temp$AIC, type = "l", lty = 2, col = "#7D7D7D")
points(as.numeric(temp$M), temp$BIC, type = "l", lty = 2, col = "#7D7D7D")
points(as.numeric(temp$M), temp$M2_1, type = "l", lty = 3,
  col = "#FF8C00", lwd = 2)
points(as.numeric(temp$M), temp$M2_2, type = "l", lty = 3,
  col = "#8FCEF3", lwd = 2)
axis(1, at = 1:4, labels = c(10, 20, 30, 40))
axis(2, at = c(0, 0.05, 0.25, 0.5, 0.75, 1))
box(lwd = 1)
mtext(expression(italic(P)), 2, 2, FALSE, at = 0.5)

temp <- res2[res2$N == 2000 & res2$rho == 2 / 3, ]
temp$M <- as.factor(temp$M)
plot(as.numeric(temp$M), temp$vt_dist, type = "l", lty = 4, ylim = c(0, 1),
  axes = FALSE, xlab = "", ylab = "",
  main = expression(italic(rho) == 2 / 3))
points(as.numeric(temp$M), temp$vt_lrt, type = "l", lty = 1)
points(as.numeric(temp$M), temp$lrt, type = "l", lty = 6, col = "#D5D5D5")
points(as.numeric(temp$M), temp$AIC, type = "l", lty = 2, col = "#7D7D7D")
points(as.numeric(temp$M), temp$BIC, type = "l", lty = 2, col = "#7D7D7D")
points(as.numeric(temp$M), temp$M2_1, type = "l", lty = 3,
  col = "#FF8C00", lwd = 2)
points(as.numeric(temp$M), temp$M2_2, type = "l", lty = 3,
  col = "#8FCEF3", lwd = 2)
axis(1, at = 1:4, labels = c(10, 20, 30, 40))
box(lwd = 1)

temp <- res2[res2$N == 2000 & res2$rho == 1 / 3, ]
temp$M <- as.factor(temp$M)
plot(as.numeric(temp$M), temp$vt_dist, type = "l", lty = 4, ylim = c(0, 1),
  axes = FALSE, xlab = "", ylab = "",
  main = expression(italic(rho) == 1 / 3))
points(as.numeric(temp$M), temp$vt_lrt, type = "l", lty = 1)
points(as.numeric(temp$M), temp$lrt, type = "l", lty = 6, col = "#D5D5D5")
points(as.numeric(temp$M), temp$AIC, type = "l", lty = 2, col = "#7D7D7D")
points(as.numeric(temp$M), temp$BIC, type = "l", lty = 2, col = "#7D7D7D")
points(as.numeric(temp$M), temp$M2_1, type = "l", lty = 3,
  col = "#FF8C00", lwd = 2)
points(as.numeric(temp$M), temp$M2_2, type = "l", lty = 3,
  col = "#8FCEF3", lwd = 2)
axis(1, at = 1:4, labels = c(10, 20, 30, 40))
box(lwd = 1)

temp <- res2[res2$N == 2000 & res2$rho == 0, ]
temp$M <- as.factor(temp$M)
plot(as.numeric(temp$M), temp$vt_dist, type = "l", lty = 4, ylim = c(0, 1),
  axes = FALSE, xlab = "", ylab = "",
  main = expression(italic(rho) == 0))
points(as.numeric(temp$M), temp$vt_lrt, type = "l", lty = 1)
points(as.numeric(temp$M), temp$lrt, type = "l", lty = 6, col = "#D5D5D5")
points(as.numeric(temp$M), temp$AIC, type = "l", lty = 2, col = "#7D7D7D")
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
legend(0.95, 0.971, legend = expression(LRT[t]), lty = 6, bty = "n",
  col = "#D5D5D5", cex = 0.8, seg.len = 2.5)
legend(0.95, 0.941, legend = "AIC", lty = 2, bty = "n", col = "#7D7D7D",
  text.col = "#7D7D7D", cex = 0.8, seg.len = 2.5)
legend(0.95, 0.911, legend = expression(LRT[v]), lty = 1, bty = "n",
  cex = 0.8, seg.len = 2.5)
legend(0.95, 0.869, legend = expression(italic(M[2])), lty = 3,
  bty = "n", col = "#FF8C00", text.col = "#FF8C00", cex = 0.8, lwd = 2,
  seg.len = 2.5)
legend(0.95, 0.806, legend = "Dist", lty = 4, bty = "n",
  cex = 0.8, seg.len = 2.5)
legend(0.95, 0.293, legend = expression(italic(M[2])), lty = 3,
  bty = "n", col = "#8FCEF3", text.col = "#8FCEF3", cex = 0.8, lwd = 2,
  seg.len = 2.5)
legend(0.95, 0.253, legend = "BIC", lty = 2, bty = "n", col = "#7D7D7D",
  text.col = "#7D7D7D", cex = 0.8, seg.len = 2.5)
legend(0.42, 0.025, c("2PL", "2d-2PL"), bty = "n", pch = c(19, 19),
  ncol = 2, cex = 0.8, col = c("#FF8C00", "#8FCEF3"))

#dev.off()

### clean up
#rm(list = c("i", "I", "path", "res2", "temp"))



###############################################################################
# Application: The Nerdy Personality Attributes Scale
###############################################################################

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

### graded vs. gpcm
graded <- mirt(dat, model = 1, itemtype = "graded", SE = TRUE,
  technical = list(NCYCLES = 20000), verbose = FALSE)
gpcm <- mirt(dat, model = 1, itemtype = "gpcm", SE = TRUE,
  technical = list(NCYCLES = 20000), verbose = FALSE)
extract.mirt(graded, "AIC")
extract.mirt(gpcm, "AIC")
vuongtest(graded, gpcm, nested = FALSE)

### 1d-graded vs. 2d-graded
graded2 <- mirt(dat, model = 2, itemtype = "graded", SE = TRUE,
  technical = list(NCYCLES = 20000), verbose = FALSE)
extract.mirt(graded2, "AIC")
anova(graded, graded2)
vuongtest(graded, graded2, nested = TRUE)

### clean up
#rm(list = c("check", "dat", "gpcm", "graded", "graded2", "raw", "temp"))


# Universal thresholds for interval mapping

rm(list = ls())

setwd("C:/Users/msmg/Desktop/BMWpop/Interval_mapping/")

library(mpMap)
library(lme4)
library(aods3)
library(qtl)
library(VPdtw)

source("Functions/mpIM.R")

# Load function into global environment to use the sourced mpIM function --------------------
sim.sigthr <- function (mpcross, nsim = 100, alpha = 0.05, pindex = 1, step = 0,
                        ncov = 0, ...)
{
  output <- list()
  if (is.null(mpcross$map))
    stop("Must have marker map to generate null simulations")
  vare <- 1
  if (!is.null(mpcross$pheno))
    vare <- var(mpcross$pheno[, pindex], na.rm = T)
  if (!(inherits(mpcross, "mpprob") && attr(mpcross$prob, "step") ==
        step))
    mpp <- mpprob(mpcross, program = "qtl", step = step)
  else mpp <- mpcross
  minp <- vector(length = nsim)
  for (i in 1:nsim) {
    ph <- rnorm(nrow(mpcross$finals), sd = vare^0.5)
    mpp$pheno[, 1] <- ph
    res <- mpIM(object = mpp, ncov = ncov, step = step, responsename = names(mpp$pheno)[1],
                mrkpos = T)
    minp[i] <- min(unlist(res$QTLresults$pvalue))
  }
  output$call <- match.call()
  output$minp <- minp
  output$thr <- sort(minp)[floor(alpha * nsim)]
  return(output)
}

# Perform permutations --------------------

load("BMWpop.Rdata")
set.seed(2020)
output <- sim.sigthr(mpcross = BMWpop, nsim = 10000, step = 1)
-log10(output$thr) # 4.256106

## ultimate calculation of div rates.


library(ape)
library(phangorn)
library(dplyr)
library(TESS)
library(geiger)

## scaling the true species tree's branch lengths

#dataset1Tree <- read.tree("/Users/eren/Documents/GitHub/Chapter3/UltrametricTrees/Dataset1/s_tree.trees")

#dataset1Tree <- rescale(dataset1Tree, "depth", 1)

#write.tree(dataset1Tree, file = "/Users/eren/Documents/GitHub/Chapter3/UltrametricTrees/Dataset1/s_tree.trees")

#dataset2Tree <- read.tree("/Users/eren/Documents/GitHub/Chapter3/UltrametricTrees/Dataset2/s_tree.trees")

#dataset2Tree <- rescale(dataset2Tree, "depth", 1)

#write.tree(dataset2Tree, file = "/Users/eren/Documents/GitHub/Chapter3/UltrametricTrees/Dataset1/s_tree.trees")

#dataset3Tree <- read.tree("/Users/eren/Documents/GitHub/Chapter3/UltrametricTrees/Dataset3/s_tree.trees")

#dataset3Tree <- rescale(dataset3Tree, "depth", 1)

#write.tree(dataset3Tree, file = "/Users/eren/Documents/GitHub/Chapter3/UltrametricTrees/Dataset3/s_tree.trees")

#dataset4Tree <- read.tree("/Users/eren/Documents/GitHub/Chapter3/UltrametricTrees/Dataset4/s_tree.trees")

#dataset4Tree <- rescale(dataset4Tree, "depth", 1)

#write.tree(dataset4Tree, file = "/Users/eren/Documents/GitHub/Chapter3/UltrametricTrees/Dataset4/s_tree.trees")


## # # # # # # # # # # # # # # # # # # # # # # # # # # #

## Birth-death processes with constant rates

Sys.setenv("DISPLAY"=":0.0")

input_dir <- "/data/schwartzlab/eren/Chapter3/UltrametricTrees"

#test_input_dir <- "/Users/eren/Documents/GitHub/Chapter3/UltrametricTrees/Dataset1/"

out_dir <- "/data/schwartzlab/eren/Chapter3/outdir"

#tree_object <- read.tree("../UltrametricTrees/Dataset1/inference_fong_best600.txt.treefile")


for(dataset in list.dirs(input_dir, recursive = F,full.names = F)){
  for(tree in list.files(paste(input_dir,dataset, sep = "/"))){

    tree_object <- read.tree(paste(input_dir,"/",dataset,"/",tree, sep = ""))

    #variables:
    #tree time

    times <- as.numeric(branching.times(tree_object))

    #defining priors for constant birth-death model

    prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
    prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }

    priorsConstBD <- c("diversification"=prior_delta,
                       "turnover"=prior_tau)

    # likelihood function
    likelihoodConstBD <- function(params) {
      speciation <- params[1] + params[2]
      extinction <- params[2]

      lnl <- tess.likelihood(times,
                             lambda = speciation,
                             mu = extinction,
                             samplingProbability = 1.0,
                             log = TRUE)
      return(lnl)
    }


    #samples for constBD

    samplesConstBD <- tess.mcmc(likelihoodFunction = likelihoodConstBD,
                                priors = priorsConstBD,
                                parameters = runif(2,0,1),
                                logTransforms = c(TRUE,TRUE),
                                delta = c(1,1),
                                iterations = 10000,
                                burnin = 1000,
                                thinning = 10,
                                adaptive = TRUE,
                                verbose = TRUE)

    #write.(summary(samplesConstBD), file = paste(out_dir,"/",dataset,"/","sum_samples_",tree,"_samplesConstBD.csv", sep=""))

    write.csv(samplesConstBD, file = paste(out_dir,"/",dataset,"/","samples_",tree,"_samplesConstBD.csv", sep=""))

    pdf(paste(out_dir,"/",dataset,"/","plot_",tree,"_samplesConstBD.pdf", sep=""),width=10, height=10)

    plot(samplesConstBD)

    dev.off()

    ## Birth-death processes with continuously varying rates

    prior_delta <- function(x) { dexp(x,rate=0.1,log=TRUE) }
    prior_lambda <- function(x) { dexp(x,rate=10.0,log=TRUE) }
    prior_alpha <- function(x) { dexp(x,rate=0.1,log=TRUE) }
    priorsDecrBD <- c("turnover"=prior_delta,
                      "initial speciation"=prior_lambda,
                      "speciation decay"=prior_alpha)

    likelihoodDecrBD <- function(params) {

      speciation <- function(t) params[1] + params[2] * exp(-params[3]*t)
      extinction <- function(t) params[1]

      lnl <- tess.likelihood(times,
                             lambda = speciation,
                             mu = extinction,
                             samplingProbability = 1.0,
                             log = TRUE)
      return (lnl)
    }

    samplesDecrBD <- tess.mcmc(likelihoodFunction = likelihoodDecrBD,
                               priors = priorsDecrBD,
                               parameters = runif(3,0,1),
                               logTransforms = c(TRUE,TRUE,TRUE),
                               delta = c(1,1,1),
                               iterations = 10000,
                               burnin = 1000,
                               thinning = 10,
                               adaptive = TRUE,
                               verbose = TRUE)

    #write.csv(summary(samplesDecrBD), file = paste(out_dir,"/",dataset,"/","sum_samples_",tree,"_samplesDecrBD", sep=""))

    write.csv(samplesDecrBD, file = paste(out_dir,"/",dataset,"/","samples_",tree,"_samplesDecrBD", sep=""))

    pdf(paste(out_dir,"/",dataset,"/","plot_",tree,"_samplesDecrBD.pdf", sep=""),height=10, width=10)

    plot(samplesDecrBD)

    dev.off()

    ## Birth-death processes with episodically varying rates

    rateChangeTime <- max( times ) / 2

    prior_delta_before <- function(x) { dexp(x,rate=10.0,log=TRUE) }
    prior_tau_before <- function(x) { dexp(x,rate=10.0,log=TRUE) }
    prior_delta_after <- function(x) { dexp(x,rate=10.0,log=TRUE) }
    prior_tau_after <- function(x) { dexp(x,rate=10.0,log=TRUE) }

    priorsEpisodicBD <- c("diversification before"=prior_delta_before,
                          "turnover before"=prior_tau_before,
                          "diversification after"=prior_delta_after,
                          "turnover after"=prior_tau_after)

    likelihoodEpisodicBD <- function(params) {

      speciation <- c(params[1]+params[2],params[3]+params[4])
      extinction <- c(params[2],params[4])


      lnl <- tess.likelihood.rateshift(times,
                                       lambda = speciation,
                                       mu = extinction,
                                       rateChangeTimesLambda = rateChangeTime,
                                       rateChangeTimesMu = rateChangeTime,
                                       samplingProbability = 1.0,
                                       log = TRUE)
      return(lnl)
    }

    samplesEpisodicBD <- tess.mcmc(likelihoodFunction = likelihoodEpisodicBD,
                                   priors = priorsEpisodicBD,
                                   parameters = runif(4,0,1),
                                   logTransforms = c(TRUE,TRUE,TRUE,TRUE),
                                   delta = c(1,1,1,1),
                                   iterations = 10000,
                                   burnin = 1000,
                                   thinning = 10,
                                   adaptive = TRUE,
                                   verbose = TRUE)

    summary(samplesEpisodicBD)


    #write.csv(summary(samplesEpisodicBD), file = paste(out_dir,"/",dataset,"/","sum_samples_",tree,"_samplesEpisodicBD", sep=""))

    write.csv(samplesEpisodicBD, file = paste(out_dir,"/",dataset,"/","samples_",tree,"_samplesEpisodicBD", sep=""))

    pdf(paste(out_dir,"/",dataset,"/","plot_",tree,"_samplesEpisodicBD.pdf", sep=""),height=10,width=10)

    plot(samplesEpisodicBD)

    dev.off()


    ### 2.4.4 Birth-death processes with explicit mass-extinction events

    survivalProbability <- 0.1

    prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
    prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }
    prior_time <- function(x) { dunif(x,min=max(times)/2,max=max(times),log=TRUE)}
    priorsMassExtinctionBD <- c("diversification"=prior_delta,
                                "turnover"=prior_tau,
                                "mass-extinction time"=prior_time)


    likelihoodMassExtinctionBD <- function(params) {
      speciation <- params[1]+params[2]
      extinction <- params[2]
      time <- params[3]

      lnl <- tess.likelihood(times,
                             lambda = speciation,
                             mu = extinction,
                             massExtinctionTimes = time,
                             massExtinctionSurvivalProbabilities = survivalProbability,
                             samplingProbability = 1.0,
                             log = TRUE)

      return (lnl)

    }

    samplesMassExtinctionBD <- tess.mcmc(likelihoodFunction = likelihoodMassExtinctionBD,
                                         priors = priorsMassExtinctionBD,
                                         parameters = c(runif(2,0,1),max(times)*3/4),
                                         logTransforms = c(TRUE,TRUE,FALSE),
                                         delta = c(1,1,1),
                                         iterations = 10000,
                                         burnin = 1000,
                                         thinning = 10,
                                         adaptive = TRUE,
                                         verbose = TRUE)

    summary(samplesMassExtinctionBD)

    #write.csv(summary(samplesMassExtinctionBD), file = paste(out_dir,"/",dataset,"/","sum_samples_",tree,"_samplesMassExtinctionBD", sep=""))

    write.csv(samplesEpisodicBD, file = paste(out_dir,"/",dataset,"/","samples_",tree,"_samplesMassExtinctionBD", sep=""))

    pdf(paste(out_dir,"/",dataset,"/","plot_",tree,"_samplesMassExtinctionBD.pdf", sep=""),height=10,width=10)

    plot(samplesEpisodicBD)

    dev.off()

  }
}

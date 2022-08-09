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

    ConstBD_prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
    ConstBD_prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }

    priorsConstBD <- c("diversification"=ConstBD_prior_delta,
                       "turnover"=ConstBD_prior_tau)

    # likelihood function
    likelihoodConstBD <- function(params) {
      ConstBD_speciation <- params[1] + params[2]
      ConstBD_extinction <- params[2]

      lnl <- tess.likelihood(times,
                             lambda = ConstBD_speciation,
                             mu = ConstBD_extinction,
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

    DecrBD_prior_delta <- function(x) { dexp(x,rate=0.1,log=TRUE) }
    DecrBD_prior_lambda <- function(x) { dexp(x,rate=10.0,log=TRUE) }
    DecrBD_prior_alpha <- function(x) { dexp(x,rate=0.1,log=TRUE) }
    priorsDecrBD <- c("turnover"=DecrBD_prior_delta,
                      "initial speciation"=DecrBD_prior_lambda,
                      "speciation decay"=DecrBD_prior_alpha)

    likelihoodDecrBD <- function(params) {

      DecrBD_speciation <- function(t) params[1] + params[2] * exp(-params[3]*t)
      DecrBD_extinction <- function(t) params[1]

      lnl <- tess.likelihood(times,
                             lambda = DecrBD_speciation,
                             mu = DecrBD_extinction,
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

    EpisodicBD_prior_delta_before <- function(x) { dexp(x,rate=10.0,log=TRUE) }
    EpisodicBD_prior_tau_before <- function(x) { dexp(x,rate=10.0,log=TRUE) }
    EpisodicBD_prior_delta_after <- function(x) { dexp(x,rate=10.0,log=TRUE) }
    EpisodicBD_prior_tau_after <- function(x) { dexp(x,rate=10.0,log=TRUE) }

    priorsEpisodicBD <- c("diversification before"=EpisodicBD_prior_delta_before,
                          "turnover before"=EpisodicBD_prior_tau_before,
                          "diversification after"=EpisodicBD_prior_delta_after,
                          "turnover after"=EpisodicBD_prior_tau_after)

    likelihoodEpisodicBD <- function(params) {

      EpisodicBD_speciation <- c(params[1]+params[2],params[3]+params[4])
      EpisodicBD_extinction <- c(params[2],params[4])


      lnl <- tess.likelihood.rateshift(times,
                                       lambda = EpisodicBD_speciation,
                                       mu = EpisodicBD_extinction,
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

    MassExtinctionBD_prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
    MassExtinctionBD_prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }
    MassExtinctionBD_prior_time <- function(x) { dunif(x,min=max(times)/2,max=max(times),log=TRUE)}
    priorsMassExtinctionBD <- c("diversification"=MassExtinctionBD_prior_delta,
                                "turnover"=MassExtinctionBD_prior_tau,
                                "mass-extinction time"=MassExtinctionBD_prior_time)


    likelihoodMassExtinctionBD <- function(params) {
      MassExtinctionBD_speciation <- params[1]+params[2]
      MassExtinctionBD_extinction <- params[2]
      MassExtinctionBD_time <- params[3]

      lnl <- tess.likelihood(times,
                             lambda = MassExtinctionBD_speciation,
                             mu = MassExtinctionBD_extinction,
                             massExtinctionTimes = MassExtinctionBD_time,
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

    write.csv(samplesMassExtinctionBD, file = paste(out_dir,"/",dataset,"/","samples_",tree,"_samplesMassExtinctionBD", sep=""))

    pdf(paste(out_dir,"/",dataset,"/","plot_",tree,"_samplesMassExtinctionBD.pdf", sep=""),height=10,width=10)

    plot(samplesMassExtinctionBD)

    dev.off()

    ##model evaluation with Bayesian Factore comparison

    marginalLikelihoodConstBD <- tess.steppingStoneSampling(
      likelihoodFunction = likelihoodConstBD,
      priors = priorsConstBD,
      parameters = runif(2,0,1),
      logTransforms = c(TRUE,TRUE),
      iterations = 1000,
      burnin = 100,
      K = 50)

    marginalLikelihoodDecrBD <- tess.steppingStoneSampling(
      likelihoodFunction = likelihoodDecrBD,
      priors = priorsDecrBD,
      parameters = runif(3,0,1),
      logTransforms = c(TRUE,TRUE,TRUE),
      iterations = 1000,
      burnin = 100,
      K = 50)

      marginalLikelihoodEpisodicBD <- tess.steppingStoneSampling(
      likelihoodFunction = likelihoodEpisodicBD,
      priors = priorsEpisodicBD,
      parameters = runif(4,0,1),
      logTransforms = c(TRUE,TRUE,TRUE,TRUE),
      iterations = 1000,
      burnin = 100,
      K = 50)


    marginalLikelihoodMassExtinctionBD <- tess.steppingStoneSampling(
      likelihoodFunction = likelihoodMassExtinctionBD,
      priors = priorsMassExtinctionBD,
      parameters = c(runif(2,0,1),max(times)*3/4),
      logTransforms = c(TRUE,TRUE,FALSE),
      iterations = 1000,
      burnin = 100,
      K = 50)

    candidateModels <- c("ConstBD"=marginalLikelihoodConstBD,
                         "DecrBD"=marginalLikelihoodDecrBD,
                         "EpisodicBD"=marginalLikelihoodEpisodicBD,
                         "MassExtinctionBD"=marginalLikelihoodMassExtinctionBD)

    marginalLikelihoodGrid <- expand.grid(M0=names(candidateModels),
                                          M1=names(candidateModels))

    marginalLikelihoodGrid$BF <- 2 * (candidateModels[marginalLikelihoodGrid$M0] -
                                        candidateModels[marginalLikelihoodGrid$M1])

    marginalLikelihoodGrid <- marginalLikelihoodGrid[order(marginalLikelihoodGrid$BF,
                                                           decreasing=TRUE),]

    write.csv(marginalLikelihoodGrid, file = paste(out_dir,"/",dataset,"/","BF_",tree,"ML_GRID",sep=""),row.names = F,sep = "\t")



    ## model evaluation with posteior-predictive simulations # ConstBD

    tmrca <- max( times )

    simConstBD <- function(params) {

      ConstBD_speciation <- params[1] + params[2]
      ConstBD_extinction <- params[2]

      repeat {
        ConstBD_tree <- tess.sim.age(n = 1,
                                     age = tmrca,
                                     lambda = ConstBD_speciation,
                                     mu = ConstBD_extinction,
                                     samplingProbability = 1.0,
                                     MRCA = TRUE)[[1]]
        if (ConstBD_tree$Nnode > 1) break }
      return (ConstBD_tree) }

    treesConstBD <- tess.PosteriorPrediction(simConstBD,samplesConstBD)


    # compute the number of species in each simulate ConstBD_tree
    numTaxaConstBD <- c()
    for (i in 1:length(treesConstBD)){
      numTaxaConstBD[i] <- treesConstBD[[i]]$Nnode + 1 }

    numTaxaPPDI_ConstBD <- quantile(numTaxaConstBD,prob=c(0.025,0.975))

    observedGamma <- gammaStat(tree_object)

    ConstBD_ppt <- tess.PosteriorPredictiveTest(treesConstBD,tree_object,
                                                gammaStat)
    ConstBD_gammaPPDI <- quantile(ConstBD_ppt[[1]],prob=c(0.025,0.975))


    pdf(paste(out_dir,"/",dataset,"/","plot_PPD_",tree,"_PPD_samplesConstBD.pdf", sep=""),height=6,width=18)

    par(mfrow=c(1,3))


    plot(density(numTaxaConstBD),main="Number of taxa",xlab="",
         ylab="Posterior Predictive Density",lwd=2) + abline(v=numTaxaPPDI_ConstBD,lty=2,col="gray",lwd=2) +
      points(tree_object$Nnode+1,0,pch="x") + title(sub = "(A)", xlab = "Frequency")


    ltt.plot(treesConstBD[[1]],backward=FALSE,col="gray",log="y",
             ylim=c(1,max(200)),main="LTT-plot") + ltt.lines(tree_object,backward=FALSE,lwd=3) +
      for (i in 2:min(100,length(treesConstBD))) ltt.lines(treesConstBD[[i]], backward=FALSE, col="gray") + title(sub = "(B)")


    plot(density(ppt[[1]]),main="Gamma Statistic",xlab="",
         ylab="Posterior Predictive Density",lwd=2) + abline(v=ConstBD_gammaPPDI,lty=2,col="gray",lwd=2) +
      points(observedGamma,0,pch="x") + title(sub = "(C)", xlab = "Gamma Value")


    dev.off()

    ## PPD DecrBD

    simDecrBD <- function(params) {

      DecrBD_speciation <- function(t) params[1] + params[2] * exp(-params[3]*t)
      DecrBD_extinction <- function(t) params[1]

      repeat {
        DecrBD_tree <- tess.sim.age(n = 1,
                                    age = tmrca,
                                    lambda = DecrBD_speciation,
                                    mu = DecrBD_extinction,
                                    samplingProbability = 1.0,
                                    MRCA = TRUE)[[1]]
        if (DecrBD_tree$Nnode > 1) break }
      return (DecrBD_tree) }

    treesDecrBD <- tess.PosteriorPrediction(simDecrBD,samplesDecrBD)


    # compute the number of species in each simulate DecrBD_tree
    numTaxaDecrBD <- c()
    for (i in 1:length(treesDecrBD)){
      numTaxaDecrBD[i] <- treesDecrBD[[i]]$Nnode + 1 }

    numTaxaPPDI_DecrBD <- quantile(numTaxaDecrBD,prob=c(0.025,0.975))

    observedGamma <- gammaStat(tree_object)

    DecrBD_ppt <- tess.PosteriorPredictiveTest(treesDecrBD,tree_object,
                                               gammaStat)
    DecrBD_gammaPPDI <- quantile(DecrBD_ppt[[1]],prob=c(0.025,0.975))


    pdf(paste(out_dir,"/",dataset,"/","plot_PPD_",tree,"_PPD_samplesDecrBD.pdf", sep=""),height=6,width=18)

    par(mfrow=c(1,3))


    plot(density(numTaxaDecrBD),main="Number of taxa",xlab="",
         ylab="Posterior Predictive Density",lwd=2) + abline(v=numTaxaPPDI_DecrBD,lty=2,col="gray",lwd=2) +
      points(tree_object$Nnode+1,0,pch="x") + title(sub = "(A)", xlab = "Frequency")


    ltt.plot(treesDecrBD[[1]],backward=FALSE,col="gray",log="y",
             ylim=c(1,max(200)),main="LTT-plot") + ltt.lines(tree_object,backward=FALSE,lwd=3) +
      for (i in 2:min(100,length(treesDecrBD))) ltt.lines(treesDecrBD[[i]], backward=FALSE, col="gray") + title(sub = "(B)")


    plot(density(ppt[[1]]),main="Gamma Statistic",xlab="",
         ylab="Posterior Predictive Density",lwd=2) + abline(v=DecrBD_gammaPPDI,lty=2,col="gray",lwd=2) +
      points(observedGamma,0,pch="x") + title(sub = "(C)", xlab = "Gamma Value")


    dev.off()

    ## PDD for Episodic

    simEpisodicBD <- function(params) {

      EpisodicBD_speciation <- c(params[1]+params[2],params[3]+params[4])
      EpisodicBD_extinction <- c(params[2],params[4])

      repeat {
        EpisodicBD_tree <- tess.sim.age(n = 1,
                                        age = tmrca,
                                        lambda = EpisodicBD_speciation,
                                        mu = EpisodicBD_extinction,
                                        samplingProbability = 1.0,
                                        MRCA = TRUE)[[1]]
        if (EpisodicBD_tree$Nnode > 1) break }
      return (EpisodicBD_tree) }

    treesEpisodicBD <- tess.PosteriorPrediction(simEpisodicBD,samplesEpisodicBD)


    # compute the number of species in each simulate EpisodicBD_tree
    numTaxaEpisodicBD <- c()
    for (i in 1:length(treesEpisodicBD)){
      numTaxaEpisodicBD[i] <- treesEpisodicBD[[i]]$Nnode + 1 }

    numTaxaPPDI_EpisodicBD <- quantile(numTaxaEpisodicBD,prob=c(0.025,0.975))

    observedGamma <- gammaStat(tree_object)

    EpisodicBD_ppt <- tess.PosteriorPredictiveTest(treesEpisodicBD,tree_object,
                                                   gammaStat)
    EpisodicBD_gammaPPDI <- quantile(EpisodicBD_ppt[[1]],prob=c(0.025,0.975))


    pdf(paste(out_dir,"/",dataset,"/","plot_PPD_",tree,"_PPD_samplesEpisodicBD.pdf", sep=""),height=6,width=18)

    par(mfrow=c(1,3))


    plot(density(numTaxaEpisodicBD),main="Number of taxa",xlab="",
         ylab="Posterior Predictive Density",lwd=2) + abline(v=numTaxaPPDI_EpisodicBD,lty=2,col="gray",lwd=2) +
      points(tree_object$Nnode+1,0,pch="x") + title(sub = "(A)", xlab = "Frequency")


    ltt.plot(treesEpisodicBD[[1]],backward=FALSE,col="gray",log="y",
             ylim=c(1,max(200)),main="LTT-plot") + ltt.lines(tree_object,backward=FALSE,lwd=3) +
      for (i in 2:min(100,length(treesEpisodicBD))) ltt.lines(treesEpisodicBD[[i]], backward=FALSE, col="gray") + title(sub = "(B)")


    plot(density(ppt[[1]]),main="Gamma Statistic",xlab="",
         ylab="Posterior Predictive Density",lwd=2) + abline(v=EpisodicBD_gammaPPDI,lty=2,col="gray",lwd=2) +
      points(observedGamma,0,pch="x") + title(sub = "(C)", xlab = "Gamma Value")


    dev.off()


    ## PPD for MassExtinction

    simMassExtinctionBD <- function(params) {

      MassExtinctionBD_speciation <- params[1]+params[2]
      MassExtinctionBD_extinction <- params[2]
      MassExtinctionBD_time <- params[3]

      repeat {
        MassExtinctionBD_tree <- tess.sim.age(n = 1,
                                              age = tmrca,
                                              lambda = MassExtinctionBD_speciation,
                                              mu = MassExtinctionBD_extinction,
                                              samplingProbability = 1.0,
                                              MRCA = TRUE)[[1]]
        if (MassExtinctionBD_tree$Nnode > 1) break }
      return (MassExtinctionBD_tree) }

    treesMassExtinctionBD <- tess.PosteriorPrediction(simMassExtinctionBD,samplesMassExtinctionBD)


    # compute the number of species in each simulate MassExtinctionBD_tree
    numTaxaMassExtinctionBD <- c()
    for (i in 1:length(treesMassExtinctionBD)){
      numTaxaMassExtinctionBD[i] <- treesMassExtinctionBD[[i]]$Nnode + 1 }

    numTaxaPPDI_MassExtinctionBD <- quantile(numTaxaMassExtinctionBD,prob=c(0.025,0.975))

    observedGamma <- gammaStat(tree_object)

    MassExtinctionBD_ppt <- tess.PosteriorPredictiveTest(treesMassExtinctionBD,tree_object,
                                                         gammaStat)
    MassExtinctionBD_gammaPPDI <- quantile(MassExtinctionBD_ppt[[1]],prob=c(0.025,0.975))


    pdf(paste(out_dir,"/",dataset,"/","plot_PPD_",tree,"_PPD_samplesMassExtinctionBD.pdf", sep=""),height=6,width=18)

    par(mfrow=c(1,3))


    plot(density(numTaxaMassExtinctionBD),main="Number of taxa",xlab="",
         ylab="Posterior Predictive Density",lwd=2) + abline(v=numTaxaPPDI_MassExtinctionBD,lty=2,col="gray",lwd=2) +
      points(tree_object$Nnode+1,0,pch="x") + title(sub = "(A)", xlab = "Frequency")


    ltt.plot(treesMassExtinctionBD[[1]],backward=FALSE,col="gray",log="y",
             ylim=c(1,max(200)),main="LTT-plot") + ltt.lines(tree_object,backward=FALSE,lwd=3) +
      for (i in 2:min(100,length(treesMassExtinctionBD))) ltt.lines(treesMassExtinctionBD[[i]], backward=FALSE, col="gray") + title(sub = "(B)")


    plot(density(ppt[[1]]),main="Gamma Statistic",xlab="",
         ylab="Posterior Predictive Density",lwd=2) + abline(v=MassExtinctionBD_gammaPPDI,lty=2,col="gray",lwd=2) +
      points(observedGamma,0,pch="x") + title(sub = "(C)", xlab = "Gamma Value")


    dev.off()


  }
}

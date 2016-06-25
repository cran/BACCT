#' @title Bayesian Augmented Control for Binary Responses
#' @description Calling JAGS to implement BAC for binary responses
#' @param yh,nh  Vector of the numbers of events (subjects) in the historical
#'   trial(s). Must be of equal length.
#' @param n1,n2  Number of subjects in the control or treatment arm of the current
#'   trial.
#' @param y1.range,y2.range  Number of events in control or treatment arm of the
#'   current trial. See "Details".
#' @param n.chain  Controls the number of posterior samples. Each chain contains
#'   20,000 samples.
#' @param tau.alpha,tau.beta  Hyperparameters of the inverse gamma distribution
#'   controling the extent of borrowing.
#' @param prior.type  Type of prior on control groups. Currenly, only the
#'   inverse-gamma prior is implemented.
#' @param criterion.type  Type of posterior quantities to be monitored. See
#'   "Details."
#' @param prob.threshold  For \code{criterion.type="prob"} only. See "Details".
#' @param sim.mode  Simulation duration reduces greatly in \code{"express"}
#'   mode, if treatment and control arms are independent. See "Details".
#' @details There are two types of posterior quantities for
#'   \code{criterion.type} argument. With \code{"diff"} option, the quantity
#'   computed is \eqn{p_{T} - p_{C}}; with \code{"prob,"} such quantity is
#'   \eqn{pr(p_{T} - p_{C}>\Delta)}, where \eqn{\Delta} is specified by
#'   \code{prob.threshold} argument.
#'
#'   By default, \code{y1.range} and \code{y2.range} cover all possible outcomes
#'   and should be left unspecified in most cases. However, when \code{n1}
#'   and/or \code{n2} is fairly large, it is acceptable to use a reduced range
#'   that covers the outcomes that are most likely (e.g., within 95\% CI) to be
#'   observed. This may help shorten the time to run MCMC.
#'
#'   Another way that can greatly shorten the MCMC running time is to specify
#'   \code{"express"} mode in \code{sim.mode} argument. Express mode reduces the
#'   number of simulations from \code{length(y1.range)*length(y2.range)} to
#'   \code{length(y1.range)+length(y2.range)}. Express mode is proper when the
#'   treatment arm rate is independent of control arm rate.
#' @return An object of class "BAC".
#' @examples
#' \dontrun{
#' library(BACCT)
#' #borrow from 3 historical trials#
#' yh = c(11,300,52);nh = c(45,877,128)
#' #specify current trial sample sizes#
#' n1 = 20;n2 = 30
#'
#' #Difference criterion type in full simulation mode#
#' obj1 = BAC_binom(yh=yh,nh=nh,n1=n1,n2=n2,n.chain=5,
#' criterion.type="diff",sim.mode="full")
#'
#' #Probability criterion type in express simulation mode#
#' obj2 = BAC_binom(yh=yh,nh=nh,n1=n1,n2=n2,n.chain=5,
#' criterion.type="prob",prob.threshold=0.1,sim.mode="express")
#'
#' #S3 method for class "BAC"
#' summary(obj1)
#' }
#' @author Hongtao Zhang
#' @import rjags
#' @importFrom stats update
#' @export

BAC_binom = function(yh,nh,n1,n2,y1.range = 0:n1,y2.range = 0:n2,n.chain =
                       5,tau.alpha = 0.001, tau.beta = 0.001, prior.type =
                       "nonmixture",criterion.type = c("diff","prob"),prob.threshold,sim.mode =
                       c("full","express"))
{
  posterior.mat = matrix(NA,length(y1.range),length(y2.range))
  cat(JAGS.binom.nonmix,file = "JAGS.model.txt",sep = "\n")

  if (criterion.type == "diff")
  {
    if (sim.mode == "full")
    {
      for (i in 1:length(y1.range))
      {
        for (j in 1:length(y2.range))
        {
          jags.data = list(
            'yh' = yh, 'nh' = nh,'y1' = y1.range[i], 'n1' = n1,'y2' = y2.range[j],'n2' = n2, 'alpha' =
              tau.alpha,'beta' = tau.beta
          )
          jags.obj = jags.model(
            "JAGS.model.txt",data = jags.data,n.chains = n.chain,n.adapt = 5000,quiet =
              T
          )
          update(jags.obj,5000,progress.bar = "none")
          tt = coda.samples(
            jags.obj,variable.names = "delta",n.iter = 20000,thin = 1,progress.bar =
              "none"
          )
          mcmc.sample = unlist(tt)
          tt1 = mean(mcmc.sample)
          posterior.mat[i,j] = tt1
          progress.info = paste("y1=",y1.range[i],";y2=",y2.range[j],sep =
                                  "")
          print(progress.info)
        }
      }
    }

    if (sim.mode == "express")
    {
      p1 = NULL
      for (i in 1:length(y1.range))
      {
        jags.data = list(
          'yh' = yh, 'nh' = nh,'y1' = y1.range[i], 'n1' = n1,'y2' = y2.range[1],'n2' = n2, 'alpha' =
            tau.alpha,'beta' = tau.beta
        )
        jags.obj = jags.model(
          "JAGS.model.txt",data = jags.data,n.chains = n.chain,n.adapt = 5000,quiet =
            T
        )
        update(jags.obj,5000,progress.bar = "none")
        tt = coda.samples(
          jags.obj,variable.names = "p1",n.iter = 20000,thin = 1,progress.bar =
            "none"
        )
        mcmc.sample = unlist(tt)
        p1 = cbind(p1,mcmc.sample)
        progress.info = paste("y1=",y1.range[i],sep =
                                "")
        print(progress.info)
      }

      p2 = NULL
      for (j in 1:length(y2.range))
      {
        jags.data = list(
          'yh' = yh, 'nh' = nh,'y1' = y1.range[1], 'n1' = n1,'y2' = y2.range[j],'n2' = n2,'alpha' =
            tau.alpha,'beta' = tau.beta
        )
        jags.obj = jags.model(
          "JAGS.model.txt",data = jags.data,n.chains = n.chain,n.adapt = 5000,quiet =
            T
        )
        update(jags.obj,5000,progress.bar = "none")
        tt = coda.samples(
          jags.obj,variable.names = "p2",n.iter = 20000,thin = 1,progress.bar =
            "none"
        )
        mcmc.sample = unlist(tt)
        p2 = cbind(p2,mcmc.sample)
        progress.info = paste("y2=",y2.range[j],sep =
                                "")
        print(progress.info)
      }

      tt1 = colMeans(p1) %o% rep(1,length(y2.range))
      tt2 = rep(1,length(y1.range)) %o% colMeans(p2)
      posterior.mat = tt2 - tt1
    }

    rownames(posterior.mat) = y1.range
    colnames(posterior.mat) = y2.range
    file.remove("JAGS.model.txt")
    rlist = list(
      response = "Binary",posterior.mat = posterior.mat,y1.range = y1.range,y2.range =
        y2.range,n1 = n1,n2 = n2,yh = yh,nh = nh,nsample = 20000 * n.chain,criterion.type =
        criterion.type, sim.mode = sim.mode, tau.alpha = tau.alpha, tau.beta = tau.beta
    )
    class(rlist) = "BAC"
    return(rlist)
  }

  if (criterion.type == "prob")
  {
    if (sim.mode == "full")
    {
      for (i in 1:length(y1.range))
      {
        for (j in 1:length(y2.range))
        {
          jags.data = list(
            'yh' = yh, 'nh' = nh,'y1' = y1.range[i], 'n1' = n1,'y2' = y2.range[j],'n2' = n2, 'alpha' =
              tau.alpha,'beta' = tau.beta
          )
          jags.obj = jags.model(
            "JAGS.model.txt",data = jags.data,n.chains = n.chain,n.adapt = 5000,quiet =
              T
          )
          update(jags.obj,5000,progress.bar = "none")
          tt = coda.samples(
            jags.obj,variable.names = "delta",n.iter = 20000,thin = 1,progress.bar =
              "none"
          )
          mcmc.sample = unlist(tt)
          tt1 = mean(mcmc.sample > prob.threshold)
          posterior.mat[i,j] = tt1
          progress.info = paste("y1=",y1.range[i],";y2=",y2.range[j],sep =
                                  "")
          print(progress.info)
        }
      }
    }

    if (sim.mode == "express")
    {
      p1 = NULL
      for (i in 1:length(y1.range))
      {
        jags.data = list(
          'yh' = yh, 'nh' = nh,'y1' = y1.range[i], 'n1' = n1,'y2' = y2.range[1],'n2' = n2, 'alpha' =
            tau.alpha,'beta' = tau.beta
        )
        jags.obj = jags.model(
          "JAGS.model.txt",data = jags.data,n.chains = n.chain,n.adapt = 5000,quiet =
            T
        )
        update(jags.obj,5000,progress.bar = "none")
        tt = coda.samples(
          jags.obj,variable.names = "p1",n.iter = 20000,thin = 1,progress.bar =
            "none"
        )
        mcmc.sample = unlist(tt)
        p1 = cbind(p1,mcmc.sample)
        progress.info = paste("y1=",y1.range[i],sep =
                                "")
        print(progress.info)
      }

      p2 = NULL
      for (j in 1:length(y2.range))
      {
        jags.data = list(
          'yh' = yh, 'nh' = nh,'y1' = y1.range[1], 'n1' = n1,'y2' = y2.range[j],'n2' = n2,'alpha' =
            tau.alpha,'beta' = tau.beta
        )
        jags.obj = jags.model(
          "JAGS.model.txt",data = jags.data,n.chains = n.chain,n.adapt = 5000,quiet =
            T
        )
        update(jags.obj,5000,progress.bar = "none")
        tt = coda.samples(
          jags.obj,variable.names = "p2",n.iter = 20000,thin = 1,progress.bar =
            "none"
        )
        mcmc.sample = unlist(tt)
        p2 = cbind(p2,mcmc.sample)
        progress.info = paste("y2=",y2.range[j],sep =
                                "")
        print(progress.info)
      }

      for (i in 1:length(y1.range))
      {
        for (j in 1:length(y2.range))
        {
          tt1 = p2[,j] - p1[,i]
          posterior.mat[i,j] = mean(tt1 > prob.threshold)
        }
      }
    }

    rownames(posterior.mat) = y1.range
    colnames(posterior.mat) = y2.range
    file.remove("JAGS.model.txt")
    rlist = list(
      response = "Binary",posterior.mat = posterior.mat,y1.range = y1.range,y2.range =
        y2.range,n1 = n1,n2 = n2,yh = yh,nh = nh,nsample = 20000 * n.chain,criterion.type =
        criterion.type, sim.mode = sim.mode, tau.alpha = tau.alpha, tau.beta = tau.beta
    )
    class(rlist) = "BAC"
    return(rlist)
  }
}

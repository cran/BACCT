#' @title Summarizing BAC
#' @description \code{summary} method for class "BAC"
#' @param object  An object of class "BAC"
#' @param ... Argument to be passed to or from other methods
#' @author Hongtao Zhang
#' @export
#'
summary.BAC=function(object,...)
{
  cat("\n****************************")
  cat("\nDATA")
  cat("\n****************************")
  cat("\nHistorical control data:\n")
  cat(paste(object$yh,collapse=","),"out of",paste(object$nh,collapse=","))
  y1.min = min(object$y1.range);y1.max = max(object$y1.range)
  y2.min = min(object$y2.range);y2.max = max(object$y2.range)
  cat("\nRange of control arm in current trial:\n")
  cat(paste(y1.min," to ",y1.max," out of ",object$n1,sep=""))
  cat("\nRange of treatment arm in current trial:\n")
  cat(paste(y2.min," to ",y2.max," out of ",object$n2,sep=""))
  cat("\n****************************")
  cat("\nMCMC PARAMETERS")
  cat("\n****************************")
  cat("\nPosterior quantity monitored:\n")
  cat(object$criterion.type)
  cat("\n# of posterior samples:\n")
  cat(object$nsample)
}

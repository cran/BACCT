#'@title Evaluating a Decision Rule
#'@description Applies a decision rule to a "BAC" class object and provides rule
#'  evaluation
#'@param object  An object of class "BAC".
#'@param decision.rule  A vector of \code{c(a,b)} specifying the thresholds for
#'  claiming significance (or probabilities of making correct go/no-go decisions
#'  at interim look). See "Details".
#'@param control.range  A vector of control rates at which the decision rule is
#'  evaluated.
#'@param es  A vector of treatment arm effect sizes, compared to control arm.
#'@param csv.name  If a name is specified, the output data set is exported in
#'  CSV format.
#'@details The decision rules specified in \code{c(a,b)} may be in the context
#'  of either interim or final analysis. At the interim, a "go" decision is made
#'  if the criterion in the "BAC" object exceeds \code{b} and a "no go" decision
#'  if such criterion is below \code{a}. Otherwise, the decision falls in the
#'  gray zone.
#'
#'  For the final analysis, the decision rule should satisfy \code{a}=\code{b}.
#'  Significance is claimed if the criterion in the "BAC" object exceeds
#'  \code{a}. Specifying an \code{a} larger than \code{b} will lead to an error.
#'
#'  For interim analysis, specified decision rule is evaluated by the
#'  probability of making a correct go or no go decision. For final analysis,
#'  power or type-I error is computed.
#'
#'  Negative \code{es} values are allowed if a lower rate is desirable.
#'@return An object of class "BACdecision".
#'@author Hongtao Zhang
#' @examples
#' \dontrun{
#' #borrow from 3 historical trials#
#' yh = c(11,300,52);nh = c(45,877,128)
#' #specify current trial sample sizes#
#' n1 = 20;n2 = 30
#' obj = BAC_binom(yh=yh,nh=nh,n1=n1,n2=n2,n.chain=5,
#' criterion.type="prob",prob.threshold=0.1,sim.mode="express")
#'
#' rule = decision_eval(obj,decision.rule=c(0.05,0.15),
#' control.range=seq(0.3,0.5,0.01),es=c(0,0.1,0.2),csv.name="result.csv")
#'
#' #S3 method for class "BACdecision"
#' plot(rule,interim=T)
#' }
#' @export
#' @importFrom stats dbinom
#' @importFrom utils write.csv

decision_eval = function(object,decision.rule,control.range,es,csv.name =
                           NULL)
{
  if (object$response == "Binary")
  {
    n1 = object$n1;n2 = object$n2
    y1.range = object$y1.range;y2.range = object$y2.range
    nogo.cut = decision.rule[1];go.cut = decision.rule[2]
    #if it makes no sense...
    if (nogo.cut > go.cut)
    {
      stop("Invalid decision rule: ")
    }
    else
    {
      outdt = data.frame()
      for (k in 1:length(es))
      {
        tt = matrix(nogo.cut,length(y1.range),length(y2.range))
        nogo.decision.mat = -1 * (object$posterior.mat < tt)
        tt = matrix(go.cut,length(y1.range),length(y2.range))
        go.decision.mat = 1 * (object$posterior.mat > tt)
        nogo.probvec = go.probvec = rep(0,length(control.range))
        for (i in 1:length(control.range))
        {
          p1 = control.range[i]
          p2 = p1 + es[k]
          tt1 = dbinom(y1.range,n1,p1)
          tt2 = dbinom(y2.range,n2,p2)
          jointmass.mat = tt1 %o% tt2
          nogo.probvec[i] = sum(-nogo.decision.mat * jointmass.mat)
          go.probvec[i] = sum(go.decision.mat * jointmass.mat)
        }
        #generate (part of) output dataset
        ttt = data.frame(control.range,es[k],nogo.probvec,go.probvec)
        outdt = rbind(outdt,ttt)
      }
      names(outdt) = c("control.rate","ES","nogo.prob","go.prob")
      rlist = list(
        response = object$response,decision.rule = decision.rule,control.range =
          control.range,data = outdt
      )
      class(rlist) = "BACdecision"
      if (!is.null(csv.name))
      {
        write.csv(outdt,csv.name,row.names = F)
      }
      return(rlist)
    }
  }

  if (object$response == "Continuous")
  {
    cat("needs work...")
  }

}

#' @details JAGS software can be downloaded from \url{http://mcmc-jags.sourceforge.net/}.
#' @examples \dontrun{
#' library(BACCT)
#'
#' #############################
#' #Example for binary response#
#' #############################
#'
#' #specify historical data
#' yh = c(11,305,52);nh = c(45,874,120)
#' #specify subjects
#' n1 = 20;n2 = 30
#' #implement BAC and wait patiently
#' post = BAC_binom(yh=yh,nh=nh,n1=n1,n2=n2,n.chain=5,
#'        criterion.type="diff",sim.mode="express")
#' #evaluate the decision
#' rule1 = decision_eval(object=post,decision.rule=c(0.05,0.05),
#'         control.range=seq(0.3,0.5,0.01),es=c(0,0.1,0.15),csv.name="rule1.csv")
#' #plot the decision evaluation
#' (fig1 = plot(rule1))
#' #continue polishing the figure
#' #add data points
#' fig1 + geom_point(size=4)
#' #replace the title
#' fig1 + ggtitle("replace title")
#' #add reference lines
#' fig1 + geom_hline(aes(yintercept=0.05)) +
#'        geom_vline(aes(xintercept=0.42),color="black",type="dashed")
#' }
#'
#' @references Viele, et al., "Use of historical control data for assessing
#'   treatment effects in clinical trials." Pharmaceutical statistics 13(1)
#'   (2014): 41-54.

"_PACKAGE"



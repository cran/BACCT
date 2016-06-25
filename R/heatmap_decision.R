#' @title Heatmap for Decision Rules
#' @description Visualizing a decision rule for binary endpoint using heatmap
#'   plots
#' @param object  An object of "BAC" class.
#' @param decision.rule A vector of \code{c(a,b)} specifying the decision rule.
#'   See help for \code{decision_eval} function.
#' @param y1.display,y2.display  A subset of control/treatment number of events
#'   to be displayed.
#' @examples
#' \dontrun{
#' #borrow from 3 historical trials#
#' yh = c(11,300,52);nh = c(45,877,128)
#' #specify current trial sample sizes#
#' n1 = 20;n2 = 30
#' obj = BAC_binom(yh=yh,nh=nh,n1=n1,n2=n2,n.chain=5,
#' criterion.type="prob",prob.threshold=0.1,sim.mode="express")
#'
#' #generate full heatmap
#' heatmap_decision(obj,decision.rule=c(0.05,0.15))
#' #generate partial heatmap
#' heatmap_decision(obj,decision.rule=c(0.05,0.15),y1.display=5:15,y2.display=10:25)
#'
#' }
#' @author Hongtao Zhang
#' @import ggplot2
#' @import reshape2
#' @export

heatmap_decision=function(object,decision.rule,y1.display=NA,y2.display=NA)
{
  n1 = object$n1;n2 = object$n2
  y1.range = object$y1.range;y2.range = object$y2.range
  nogo.cut = decision.rule[1];go.cut = decision.rule[2]
  tt = matrix(nogo.cut,length(y1.range),length(y2.range))
  nogo.decision.mat = -1 * (object$posterior.mat < tt)
  tt = matrix(go.cut,length(y1.range),length(y2.range))
  go.decision.mat = 1 * (object$posterior.mat >= tt)
  decision.mat = nogo.decision.mat + go.decision.mat

  if(all(is.na(y1.display))){y1.display = y1.range}
  if(all(is.na(y2.display))){y2.display = y2.range}
  decision.long = reshape2::melt(decision.mat[y1.display+1,y2.display+1])
  decision.long$valuecat = factor(decision.long$value)

  if(decision.rule[1]==decision.rule[2])
  {
    cat("needs work")
  }
  if(decision.rule[1] < decision.rule[2])
  {
    #Basic elements
    ggplot(decision.long,aes_string(x="Var2",y="Var1",fill="valuecat")) + geom_tile() + geom_tile(color="black",linetype="longdash",show.legend=F) +
      #Remove gray background
      theme_bw() +
      #Axes
      scale_x_continuous(breaks=y2.display,name="Treatment") + scale_y_continuous(breaks=y1.display,name="Control") +
      #Legend
      scale_fill_manual(breaks=c("-1","0","1"),labels=c("No Go","Gray","Go"),values=c("red","gray","darkgreen"),name="Decision") +
      #Fine tunes
      theme(axis.title.x = element_text(face="bold",size=15),
            axis.text.x  = element_text(size=12,face="bold"),
            axis.title.y = element_text(face="bold",size=15),
            axis.text.y  = element_text(size=12,face="bold"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank()
      )


  }
}

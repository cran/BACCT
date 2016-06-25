#' @title Generateing Plot(s) Used for Decision Rule Evaluation
#' @description \code{plot} method for class "BACdecision"
#' @param x  An object of "BACdecision" class.
#' @param es.null  Effect size under the null hypothesis. Default is 0.
#' @param es.null.side  "=" is the only option now.
#' @param interim  Logical indicator of interim analysis (versus final
#'   analysis). Figures will differ in legends and titles.
#' @param ... Argument to be passed to or from other methods
#' @details If \code{interim=F}, only one power/type I error figure will be
#'   generated. Otherwise, two figures will be generated correponding to "No Go"
#'   and "Go" decisions respectively.
#' @return An object of "ggplot" class. Certain further edits are still allowed,
#'   such as changing title and adding reference lines.
#' @author Hongtao Zhang
#' @import ggplot2
#' @export

plot.BACdecision=function(x,es.null=0,es.null.side,interim=F,...)
{
  evaldata = x$data
  evaldata$EScat = factor(evaldata$ES)
  if(interim)
  {
    evaldata$nogo.status = ifelse(evaldata$ES==es.null,"Correct","Incorrect")
    evaldata$go.status = ifelse(evaldata$ES==es.null,"Incorrect","Correct")
  }
  if(!interim)
  {
    evaldata$nogo.status = ifelse(evaldata$ES==es.null,"Power","Type I error")
    evaldata$go.status = ifelse(evaldata$ES==es.null,"Type I error","Power")
  }

  if(x$response=="Binary")
  {
    nogo.plot =
      #Basic elements
      ggplot(evaldata,aes_string(x="control.rate",y="nogo.prob",group="EScat",linetype="EScat",color="nogo.status")) + geom_line(size=1.5) +
      #Axes
      scale_x_continuous(breaks=x$control.range,name="Control Rate") + scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.05),name="Probability") +
      #Legend box names
      scale_linetype_discrete(name="Effect Size") +
      #Fine tunes of title fonts, legend fonts
      theme(legend.title=element_text(size=15,face="bold"),
            legend.text=element_text(size=12),
            legend.key.size=unit(2,"cm"),
            plot.title = element_text(size=20,face="bold"),
            axis.title.x = element_text(face="bold",size=15),
            axis.text.x  = element_text(angle=90,vjust=0.5, size=12,face="bold"),
            axis.title.y = element_text(face="bold",size=15),
            axis.text.y  = element_text(size=12,face="bold")
            )
    if(interim)
    {
      nogo.plot = nogo.plot +
      ggtitle("No Go Decision Evaluation") +
      scale_colour_manual(name="Decision",values=c("Correct"="green4","Incorrect"="red"))
    }
    if(!interim)
    {
      nogo.plot = nogo.plot +
        ggtitle("Type I Error and Power") +
        scale_colour_manual(name=" ",values=c("Power"="green4","Type I error"="red"))
    }


    go.plot =
      #Basic elements
      ggplot(evaldata,aes_string(x="control.rate",y="go.prob",group="EScat",linetype="EScat",color="go.status")) + geom_line(size=1.5) +
      #Axes
      scale_x_continuous(breaks=x$control.range,name="Control Rate") + scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.05),name="Probability") +
      #Legend box names
      scale_linetype_discrete(name="Effect Size") +
      #Fine tunes of title fonts, legend fonts
      theme(legend.title=element_text(size=15,face="bold"),
            legend.text=element_text(size=12),
            legend.key.size=unit(2,"cm"),
            plot.title = element_text(size=20,face="bold"),
            axis.title.x = element_text(face="bold",size=15),
            axis.text.x  = element_text(angle=90,vjust=0.5, size=12,face="bold"),
            axis.title.y = element_text(face="bold",size=15),
            axis.text.y  = element_text(size=12,face="bold")
            )
    if(interim)
    {
      go.plot = go.plot +
        ggtitle("Go Decision Evaluation") +
        scale_colour_manual(name="Decision",values=c("Correct"="green4","Incorrect"="red"))
      rlist = list(nogo.plot,go.plot)
    }
    if(!interim)
    {
      go.plot = go.plot +
        ggtitle("Type I Error and Power") +
        scale_colour_manual(name=" ",values=c("Power"="green4","Type I error"="red"))
      rlist = list(go.plot)
    }
    return(rlist)
  }
}

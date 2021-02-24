#' Plot KM curve
#'
#'This function will plot a KM curve with possible stratification. You can specifyif you want
#'a legend or confidence bands as well as the units of time used.
#'
#' @param data dataframe containing your data
#' @param response character vector with names of columns to use for response
#' @param group string specifiying the column name of stratification variable
#' @param pos what position you want the legend to be. Current option are bottomleft and topright
#' @param units string specifying what the unit of time is use lower case and plural
#' @param CI boolean to specify if you want confidence intervals
#' @param legend boolean to specify if you want a legend
#' @param title title of plot
#' @importFrom graphics axis legend lines mtext par plot
#' @keywords plot
#' @export
#' @examples
#' require(survival)
#' lung$sex<-factor(lung$sex)
#' plotkm(lung,c("time","status"))
#' plotkm(lung,c("time","status"),"sex")
plotkm<-function(data,response,group=1,pos="bottomleft",units="months",CI=F,legend=T, title=""){
  if(class(group)=="numeric"){
    kfit<-survfit(as.formula(paste("Surv(",response[1],",",response[2],")~1",sep="")),data=data)
    sk<-summary(kfit)$table
    levelnames<-paste("N=",sk[1], ", Events=",sk[4]," (",round(sk[4]/sk[1],2)*100,"%)",sep="")
    if(title=="")  title<-paste("KM-Curve for ",nicename(response[2]),sep="")

  }else if(length(group)>1){
    return("Currently you can only stratify by 1 variable")
  }else{
    if(class(data[,group])!="factor")
      stop("group must be a vactor variable. (Or leave unspecified for no group)")
    lr<-survdiff(as.formula(paste("Surv(",response[1],",",response[2],")~", paste(group,collapse="+"),sep="")),data=data)
    lrpv<-1-pchisq(lr$chisq, length(lr$n)- 1)
    levelnames<-levels(data[,group])
    kfit<-survfit(as.formula(paste("Surv(",response[1],",",response[2],")~", paste(group,collapse="+"),sep="")),data=data)
    if(title=="") title<-paste("KM-Curve for ",nicename(response[2])," stratified by ", nicename(group),sep="")
    levelnames<-sapply(1:length(levelnames), function(x){paste(levelnames[x]," n=",lr$n[x],sep="")})

  }


  plot(kfit,mark.time=T, lty=1:length(levelnames),xlab=paste("Time (",cap(units),")",sep=""),
       ylab="Suvival Probability ",cex=1.1, conf.int=CI,
       main=title)


  if(legend){
    if(class(group)=="numeric"){legend(pos,levelnames,lty=1:length(levelnames),bty="n")
    }else{ legend(pos,c(levelnames,paste("p-value=",pvalue(lrpv)," (Log Rank)",sep="")),
                  col=c(rep(1,length(levelnames)),"white"),lty=1:(length(levelnames)+1),bty="n")}
  }
}

#'Get event time summary dataframe
#'
#'This function will output a dataframe with usefull summary statistics from a coxph model
#'
#'@param data dataframe containing data
#'@param response character vector with names of columns to use for response
#'@param group string specifiying the column name of stratification variable
#'@param times numeric vector of times you want survival time provbabilities for.
#'@keywords dataframe
#'@export
#'@examples
#'require(survival)
#'lung$sex<-factor(lung$sex)
#'etsum(lung,c("time","status"),"sex")
#'etsum(lung,c("time","status"))
#'etsum(lung,c("time","status"),"sex",c(1,2,3))
etsum<- function(data,response,group=1,times=c(12,24)){
  if(class(group)=="numeric"){
    kfit<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",group,sep=""))  ,data=data))
    maxtime=max(kfit$time)
    times[times>maxtime]=maxtime
    kfit2<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",group,sep="")) ,data=data),times=times)
    tab<-as.data.frame(cbind(strata=as.character(kfit2$strata),times=kfit2$time,SR=paste(round(kfit2$surv*100,0)," (",round(kfit2$lower*100,0),"-",round(kfit2$upper*100,0),")",sep="")))
    tbl<-kfit2$table
  }else{
    if(class(data[,group])!="factor")
      stop("group variable must be factor or leave unspecified for no group")
    tab<-lapply(levels(data[,group]),function(level){
      subdata<-subset(data,data[,group]==level)
      kfit<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",1,sep=""))  ,data=subdata))
      maxtime=max(kfit$time)
      times[times>maxtime]=maxtime
      kfit2<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",1,sep="")) ,data=subdata),times=times)
      list(cbind(strata=paste0(group,"=",level),times=kfit2$time,SR=paste(round(kfit2$surv*100,0)," (",round(kfit2$lower*100,0),"-",round(kfit2$upper*100,0),")",sep="")),kfit2$table)})
    tbl=t(sapply(tab,"[[",2))
    rownames(tbl)=sapply(levels(data[,group]),function(level)paste0(group,"=",level))
    tab=do.call(rbind.data.frame,lapply(tab,"[[",1))
  }

  if(class(group)!="numeric"){
    kfit<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",group,sep=""))  ,data=data))
    med=by(data,data[,group],function(x) median(x[,response[1]],na.rm=T))
    min=by(data,data[,group],function(x) min(x[,response[1]],na.rm=T))
    max=by(data,data[,group],function(x) max(x[,response[1]],na.rm=T))
    survtimes<-data.frame(strata=as.character(kfit$strata),kfit$time)
    minst<-round(as.numeric(by(survtimes,survtimes$strata,function(x) min (x[,2]))),1)
    maxst<-round(as.numeric(by(survtimes,survtimes$strata,function(x) max (x[,2]))),1)
    tab<-reshape::cast(tab, strata ~ times)
    names<-names(tab)
    tab<-data.frame(tab)
    names(tab)<-names
    tab[,1]<-levels(data[,group])
    if(length(times)>1){
      indx<-c(0,sapply(sort(as.numeric(names(tab)[-1])),function(x){which(as.numeric(names(tab)[-1])==x)}))+1
      tab<-tab[,indx]
      tab<-tab[c(2:length(tab),1)]
    }else{
      tab<-tab[c(2:length(tab),1)]
    }
    noeventsindx<-ifelse(length(which(tbl[,4]==0))!=0,
                         which(tbl[,4]==0),NA)
    if(!is.na(noeventsindx)){
      for(i in noeventsindx){
        if(i==1){
          minst<-c(0,minst)
          maxst<-c(0,maxst)
        }else if(i>length(minst)){
          minst<-c(minst,0)
          maxst<-c(maxst,0)
        }else{
          minst<-c(minst[1:i-1],0,minst[i:length(minst)])
          maxst<-c(maxst[1:i-1],0,maxst[i:length(maxst)])
        }}}


    tab<-cbind("n"=tbl[,1],"Events"=tbl[,4], "MedKM"=round(tbl[,5],1),
               "LCI"=round(tbl[,6],1), "UCI"=round(tbl[,7],1),
               "MedFU"=round(as.numeric(med),1),
               "MinFU"=round(as.numeric(min),1),"MaxFU"=round(as.numeric(max),1),
               "MinET"=minst,"MaxET"=maxst,tab)
    rownames(tab)<-NULL
  }else{
    med=median(data[,response[1]],na.rm=T)
    min=min(data[,response[1]],na.rm=T)
    max=max(data[,response[1]],na.rm=T)
    if(length(times)>1){
      tab<-data.frame(t(tab))
      rownames(tab)<-NULL
      names(tab)<-as.numeric(as.matrix(tab[1,]))
      tab<-tab[-1,]
    }else{
      rownames(tab)<-NULL
      names(tab)[2]<-times
      tab<-tab[-1]
    }
    tab<-cbind("n"=tbl[1],"Events"=tbl[4],"MedKM"=round(tbl[5],1),"LCI"=round(tbl[6],1), "UCI"=round(tbl[7],1),
               "MedFU"=round(as.numeric(med),1),"MinFU"=round(as.numeric(min),1),"MaxFU"=round(as.numeric(max),1),
               "MinET"=round(min(kfit$time),1),"MaxET"=round(max(kfit$time),1),tab)
    rownames(tab)<-NULL
  }
  return(tab)
}

#'Print LaTeX event time summary
#'
#'Wrapper for the etsum function that prints paragraphs of text in LaTeX
#'
#'@param data dataframe containing data
#'@param response character vector with names of columns to use for response
#'@param group string specifiying the column name of stratification variable
#'@param times numeric vector of times you want survival time provbabilities for.
#'@param units string indicating the unit of time. Use lower case and plural.
#'@keywords print
#'@export
#'@examples
#'require(survival)
#'lung$sex<-factor(lung$sex)
#'petsum(lung,c("time","status"),"sex")
#'petsum(lung,c("time","status"))
#'petsum(lung,c("time","status"),"sex",c(1,2,3),"months")
petsum<-function(data,response,group=1,times=c(12,14),units="months"){
  t<-etsum(data,response,group,times)

  #plotkm(nona,response,group)

  names<-names(t)
  if("strata"%in% names){
    strta<-sapply(t[,"strata"], function(x) paste(x,": ",sep=""))
    offset<-2
    ofst<-1
  }else{
    strta=matrix(c("",""))
    offset<-1
    ofst<-0
  }


  out<-sapply(seq_len(nrow(t)),function(i){

    if(is.na(t[i,3])) {km<-paste("The KM median event time has not been achieved due to lack of events.",sep="")
    }else if (!is.na(t[i,5])){km<-paste("The KM median event time is ",t[i,3]," with 95",sanitizestr("%")," confidence Interval (",t[i,4],",",t[i,5],").",sep="")
    }else{km<-paste("The KM median event time is ",t[i,3], " ",units, " with 95",sanitizestr("%")," confidence Interval (",t[i,4],",",t[i,10],").",sep="")}

    # if at least one event
    if(t[i,2]!=0){
      flet<-paste(" The first and last event times occurred at ",t[i,9],
                  " and ",t[i,10]," ",units," respectively. ",sep="")

      psindex=11:(ncol(t)-ofst)
      psindex=psindex[which(!is.na(t[i,psindex]))]
      if(length(psindex)>1){
        lastindex=psindex[length(psindex)]
        firstindex=psindex[-length(psindex)]
        ps<-paste("The ",paste(names[firstindex], collapse=", ")," and ",names[lastindex], " " , substring(units,1,nchar(units)-1),
                  " probabilities of 'survival' and their 95",sanitizestr("%")," confidence intervals are ",
                  paste(sapply(t[i,firstindex],function(x) paste(x)),collapse=", ")," and ", t[i,lastindex], " percent.",sep="")

      }else{
        ps<-paste("The ",names[psindex]," ", substring(units,1,nchar(units)-1),
                  " probability of 'survival' and 95",sanitizestr("%")," confidence interval is ",
                  t[i,psindex]," percent.",sep="")
      }
      #if no events
    }else{
      km=""
      ps=""
      flet=""
    }


    out<-paste(lbld(sanitizestr(nicename(strta[i])))," There are ",t[i,1]," patients. There were ",t[i,2],
               " (",round(100*t[i,2]/t[i,1],0),sanitizestr("%"),") events. The median and range of the follow-up times is ",
               t[i,6]," (",t[i,7],"-",t[i,8],") ",units,". ", km, flet,ps,sep="")
    cat("\n",out,"\n")
  })
}



#'Returns a dataframe corresponding to a descriptive table
#'
#'@param data dataframe containing data
#'@param covs character vector with the names of columns to include in table
#'@param maincov covariate to stratify table by
#'@param numobs named list overriding the number of people you expect to have the covariate
#'@param markup boolean indicating if you want latex markup
#'@param sanitize boolean indicating if you want to sanitize all strings to not break LaTeX
#'@param nicenames booling indicating if you want to replace . and _ in strings with a space
#'@param IQR boolean indicating if you want to display the inter quantile range (Q1,Q3) as opposed to (min,max) in the summary for continuous variables
#'@param digits.cat number of digits for the proportions when summarizing categorical data (default: 0)
#'@param testcont test of choice for continuous variables, one of \emph{rank-sum} (default) or \emph{ANOVA}
#'@param testcat test of choice for categorical variables, one of \emph{Chi-squared} (default) or \emph{Fisher}
#'@param excludeLevels a named list of levels to exclude in the form list(varname =c('level1','level2')) these levels will be excluded from association tests
#'@keywords dataframe
#'@return A formatted table displaying a summary of the covariates stratified by maincov
#'@export
#'@seealso \code{\link{fisher.test}}, \code{\link{chisq.test}}, \code{\link{wilcox.test}}, \code{\link{kruskal.test}}, and \code{\link{anova}}
covsum<-function(data,covs,maincov=NULL,numobs=NULL,markup=TRUE,sanitize=TRUE,nicenames=TRUE, IQR = FALSE, digits.cat = 0,
                 testcont = c('rank-sum test','ANOVA'), testcat = c('Chi-squared','Fisher'),excludeLevels=NULL){
  # New LA 18 Feb, test for presence of variables in data and convert character to factor
  missing_vars = setdiff(covs,names(data))
  if (length(missing_vars)>0){
    stop(paste('These covariates are not in the data:',missing_vars))
  }
  for (v in c(maincov,covs)) if (class(data[[v]])[1]=='character') data[[v]] <- factor(data[[v]])
  testcont <- match.arg(testcont)
  testcat <- match.arg(testcat)
  if(!markup){
    lbld<-identity
    addspace<-identity
    lpvalue<-identity
  }
  digits.cat <- as.integer(digits.cat)
  if( digits.cat<0 ) stop("parameter 'digits.cat' cannot be negative!")
  if(!sanitize) sanitizestr<-identity
  if(!nicenames) nicename<-identity
  if(!is.null(maincov)){
    levels<-names(table(data[[maincov]]))
    levels<-c(list(levels),as.list(levels))
  }else{
    levels<-"NOMAINCOVNULLNA"
  }
  N=nrow(data)
  if(!is.null(maincov)){
    nmaincov<-c(sum(table(data[[maincov]])),table(data[,maincov]))
  }else{
    nmaincov<-N
    p<-NULL
  }
  out<-lapply(covs,function(cov){
    ismiss=F
    n<-sum(table(data[[cov]]))
    # Exclude specified levels
    if (!is.null(excludeLevels[[cov]])){
      excludeLevel = excludeLevels[[cov]]
    } else excludeLevel = ''

    #Set up the first coulmn
    factornames<-NULL
    if(is.null(numobs[[cov]]))  numobs[[cov]]<-nmaincov
    if(numobs[[cov]][1]-n>0) {ismiss=T
    factornames<-c(factornames,"Missing")
    }
    #if the covariate is a factor
    if(is.factor(data[[cov]])){
      factornames<-c(levels(data[[cov]]),factornames)
      if (!is.null(maincov)) {
        pdata = data[!(data[[cov]] %in% excludeLevel),]

        p <- if( testcat=='Fisher'){ try(fisher.test(pdata[[maincov]], pdata[[cov]])$p.value)
        } else try(chisq.test(pdata[[maincov]], pdata[[cov]])$p.value)
        if (class(p) == "try-error")
          p <- NA
        p <- lpvalue(p)
      }


      #set up the main columns
      onetbl<-mapply(function(sublevel,N){
        missing<-NULL
        if(sublevel[1]!="NOMAINCOVNULLNA"){
          subdata<-subset(data,subset=data[[maincov]]%in%sublevel)
        }else{
          subdata<-data
        }
        table<-table(subdata[[cov]])
        tbl<-table(subdata[[cov]])
        n<-sum(tbl)
        prop <- round(tbl/n,2+digits.cat)*100
        prop <- sapply(prop,function(x){if(!is.nan(x)){x} else{0}})
        prop.fmt <- sprintf(paste0("%.",digits.cat,"f"),prop)
        tbl<-mapply(function(num,prop){paste(num," (", prop,")",sep="")},tbl,prop.fmt)
        if(ismiss) missing<-N-n
        tbl<-c(tbl,lbld(missing))
        return(tbl)
      },levels,numobs[[cov]])

      #if the covariate is not a factor
    }else{
      #setup the first column
      factornames <- c("Mean (sd)", ifelse(IQR, "Median (Q1,Q3)", "Median (Min,Max)"), factornames)
      if (!is.null(maincov)) {
        p <- if( testcont=='rank-sum test'){
          if( length(unique(data[[maincov]]))==2 ){
            try( wilcox.test(data[[cov]] ~ data[[maincov]])$p.value )
          } else try( kruskal.test(data[[cov]] ~ data[[maincov]])$p.value )
        } else try(anova(lm(data[[cov]] ~ data[[maincov]]))[5][[1]][1])
        if (class(p) == "try-error")
          p <- NA
        p <- lpvalue(p)
      }
      #set up the main columns
      onetbl <- mapply(function(sublevel,N){
        missing <- NULL
        if(sublevel[1]!="NOMAINCOVNULLNA"){
          subdata<-subset(data,subset=data[[maincov]]%in%sublevel)
        }else{subdata<-data}
        #if there is a missing in the whole data
        if(ismiss){
          n<-sum(table(subdata[[cov]]))
          missing<-N-n
        }
        # Updated LA to remove NaN from tables
        sumCov <-gsub(" ","",format(round(summary(subdata[[cov]]), 1),1))
        if (sumCov[4]=="NaN"){
          meansd <-''
          mmm <-''} else {meansd <- paste(sumCov[4], " (", round(sd(subdata[[cov]], na.rm = T), 1), ")", sep = "")
          mmm <- if (IQR) {
            paste(sumCov[3], " (", sumCov[2], ",", sumCov[5],
                  ")", sep = "")
          } else paste(sumCov[3], " (", sumCov[1], ",",
                       sumCov[6], ")", sep = "")}
        tbl <- c(meansd, mmm, lbld(missing))

        return(tbl)}
        ,levels,numobs[[cov]])}

    #Add the first column to the main columns and get the matrix ready for later
    factornames<-addspace(sanitizestr(nicename(factornames)))
    # LA Added 20 Jan 2021 to deal with one-level factors
    if (is.null(nrow(onetbl))){onetbl <- matrix(data=onetbl,ncol=length(onetbl),nrow=1) }

    onetbl<-cbind(factornames,onetbl)

    if(!is.null(maincov)){
      onetbl<-rbind(c(lbld(sanitizestr(nicename(cov))),rep("",length(levels[[1]])+1)),onetbl)
      p_NA = rep("", nrow(onetbl) - 1)
      p_NA[factornames %in% excludeLevel] <-'excl'
      onetbl <- cbind(onetbl, c(p,p_NA))

    }else{
      onetbl<-rbind(c(lbld(sanitizestr(nicename(cov))),""),onetbl)
    }
    rownames(onetbl)<-NULL
    colnames(onetbl)<-NULL
    return(onetbl)})
  table <- do.call("rbind", lapply(out, data.frame, stringsAsFactors = FALSE))
  ### unlist each column of the table
  table <- data.frame(apply(table,2,unlist), stringsAsFactors = FALSE)
  rownames(table)<-NULL
  if(!is.null(maincov)){
    colnames(table)<-c("Covariate",paste("Full Sample (n=",N,")",sep=""),
                       mapply(function(x,y){paste(x," (n=",y,")",sep="")},
                              names(table(data[[maincov]])),table(data[[maincov]])),"p-value")
  }else{
    colnames(table)<-c("Covariate",paste("n=",N,sep=""))

  }
  colnames(table)<-sanitizestr(colnames(table))
  return(table)
}


# # Problem:
# # when spss data are imported with haven and converted to factors
# # the code is.factor(test[,x_var]) returns FALSE
# # however, is.factor(test[[x_var]]) returns TRUE
# # Note: I was having trouble debugging with the cov variable, because
# # cov is a function in the Kendall package so I changed the name of the
# # cov2 and cov variables to x_var_str and x_var
#
# # Example of problem
# library(haven)
# path <- system.file("examples", "iris.sav", package = "haven")
# test <- read_sav(path)
# test$new_var <- haven::as_factor(test$Species,ordered=FALSE)
# x_var = 'new_var'
# is.factor(test[,x_var])  # returns FALSE
# is.factor(test[[x_var]]) # returns TRUE as expected
# 5 Jan 2021 - added the option of adjustment covariates (adj_cov). If these are specified, they are adjusted for, ie, included in each of the model fits.

#'Get univariate summary dataframe
#'
#'Returns a dataframe corresponding to a univariate table
#'
#'@param response string vector with name of response
#'@param covs character vector with the names of columns to fit univariate models to
#'@param data dataframe containing data
#'@param adj_cov character vector containing covariates to adjust for in each (no longer univariate) model
#'@param type string indicating he type of univariate model to fit. The function will try and guess what type you want based on your response. If you want to override this you can manually specify the type.
#'Options in clude "linear", "logistic", "coxph", "crr", "boxcox","logistic"
#'@param strata character vector of covariates to stratify by. Only used for coxph and crr
#'@param markup boolean indicating if you want latex markup
#'@param sanitize boolean indicating if you want to sanitize all strings to not break LaTeX
#'@param nicenames booling indicating if you want to replace . and _ in strings with a space
#'@param testing boolean to indicate if you want to print out the covariates before the model fits.
#'@param showN boolean indicating if you want to show sample sizes
#'@param CIwidth width of confidence interval, default is 0.95
#'This will allow you to see which model is not fitting if the function throws an error
#'@keywords dataframe
#'@export
uvsum <- function (response, covs, data, adj_cov=NULL,type = NULL, strata = 1, markup = T,
                      sanitize = T, nicenames = T, testing = F,showN=T,CIwidth=0.95)
{
  # New LA 24 Feb, test for presence of variables in data and convert character to factor
  missing_vars = setdiff(c(response,covs),names(data))
  if (length(missing_vars)>0){
    stop(paste('These covariates are not in the data:',missing_vars))
  }
  for (v in c(response,covs)) if (class(data[[v]])[1]=='character') data[[v]] <- factor(data[[v]])
  if (!markup) {
    lbld <- identity
    addspace <- identity
    lpvalue <- identity
  }
  if (!sanitize)
    sanitizestr <- identity
  if (!nicenames)
    nicename <- identity
  if (class(strata) != "numeric") {
    strataVar=strata
    strata <- sapply(strata, function(stra) {
      paste("strata(", stra, ")", sep = "")
    })
  } else {
    strataVar <-""
    strata <- ""
  }
  if (!is.null(type)) {
    if (type == "logistic") {
      beta <- "OR"
    } else if (type == "linear" | type == "boxcox") {
      beta <- "Estimate"
    } else if (type == "coxph" | type == "crr") {
      beta <- "HR"
    } else {
      stop("type must be either coxph, logisitc, linear, coxbox, crr (or NULL)")
    }
  } else {
    if (length(response) == 2) {
      if (length(unique(data[[response[2]]])) < 3) {
        type <- "coxph"
      } else {
        type <- "crr"
      }
      beta <- "HR"
    } else if (length(unique(data[[response]])) == 2) {
      type <- "logistic"
      beta <- "OR"
    } else {
      type <- "linear"
      beta <- "Estimate"
    }
  }
  beta = betaWithCI(beta,CIwidth)

  if (strata != "" & type != "coxph") {
    stop("strata can only be used with coxph")
  }
  if (is.null(adj_cov)){
    adj_cov=''
  } else{
    adj_cov = paste('+',paste(adj_cov,collapse='+'))
  }
  out <- lapply(covs, function(x_var) {
    testData <- dplyr::select(data,any_of(c(response,x_var,strataVar,adj_cov)))
    testData <-na.omit(testData)
    x_var_str <- x_var
    if (testing)
      print(x_var)
    if (is.factor(data[[x_var]])) {

      testData[[x_var]] <- factor(testData[[x_var]],ordered = FALSE)
      levelnames = sapply(sapply(sapply(levels(testData[[x_var]]),nicename),sanitizestr),addspace)
      #levelnames <- sapply(sapply(sapply(levels(factor(testData[[x_var]])), nicename), sanitizestr), addspace)
      x_var_str <- lbld(sanitizestr(nicename(x_var)))
      title <- NULL
      body <- NULL
      if (type == "coxph") {
        if (adj_cov!='') stop('Adjustment covariates should be specified as strata for coxph models.')
        m2 <- coxph(as.formula(paste(paste("Surv(", response[1],
                                           ",", response[2], ")", sep = ""), "~", x_var,
                                     ifelse(strata == "", "", "+"), paste(strata,
                                                                          collapse = "+"), sep = "")), data = testData)
        hazardratio <- c("Reference", apply(matrix(summary(m2,conf.int=CIwidth)$conf.int[,
                                                                        c(1, 3, 4)], ncol = 3), 1, psthr))
        pvalue <- c("", sapply(summary(m2,conf.int=CIwidth)$coef[, 5],
                               lpvalue))
        title <- c(x_var_str, "", "", lpvalue(summary(m2,conf.int=CIwidth)$waldtest[3]))
      } else if (type == "crr") {
        if (adj_cov!='') stop('Adjustment covariates have not been implemented for competing risk models.')
        m2 <- crrRx(as.formula(paste(paste(response,
                                           collapse = "+"), "~", x_var, sep = "")), data = testData)
        hazardratio <- c("Reference", apply(matrix(summary(m2,conf.int=CIwidth)$conf.int[,
                                                                        c(1, 3, 4)], ncol = 3), 1, psthr))
        pvalue <- c("", sapply(summary(m2)$coef[, 5],
                               lpvalue))
        globalpvalue <- try(aod::wald.test(b = m2$coef, Sigma = m2$var,
                                      Terms = seq_len(length(m2$coef)))$result$chi2[3])
        if (class(globalpvalue) == "try-error")
          globalpvalue <- "NA"
        title <- c(x_var_str, "", "", lpvalue(globalpvalue))
      } else if (type == "logistic") {
        m2 <- glm(as.formula(paste(response, "~", x_var, ifelse(adj_cov=='','',adj_cov),
                                   sep = "")), family = "binomial", data = testData)
        globalpvalue <- try(aod::wald.test(b = m2$coefficients[-1],
                                      Sigma = vcov(m2)[-1, -1], Terms = seq_len(length(m2$coefficients[-1])))$result$chi2[3])
        if (class(globalpvalue) == "try-error")
          globalpvalue <- "NA"
        m <- summary(m2)$coefficients
        Z_mult = qnorm(1-(1-CIwidth)/2)
        hazardratio <- c("Reference", apply(cbind(exp(m[-1,
                                                        1]), exp(m[-1, 1] - Z_mult * m[-1, 2]), exp(m[-1,
                                                                                                      1] + Z_mult * m[-1, 2])), 1, psthr))
        pvalue <- c("", sapply(m[-1, 4], lpvalue))
        title <- c(x_var_str, "", "", lpvalue(globalpvalue))
      } else if (type == "linear" | type == "boxcox") {
        if (type == "linear") {
          m2 <- lm(as.formula(paste(response, "~", x_var,ifelse(adj_cov=='','',adj_cov),
                                    sep = "")), data = testData)
        } else {
          m2 <- boxcoxfitRx(as.formula(paste(response,
                                             "~", x_var, ifelse(adj_cov=='','',adj_cov),
                                             sep = "")), data = testData)
        }
        m <- summary(m2)$coefficients
        globalpvalue <- try(aod::wald.test(b = m2$coefficients[-1],
                                      Sigma = vcov(m2)[-1, -1], Terms = seq_len(length(m2$coefficients[-1])))$result$chi2[3])
        if (class(globalpvalue) == "try-error")
          globalpvalue <- "NA"
        T_mult = qt(1-(1-CIwidth)/2,m2$df.residual)
        hazardratio <- c("Reference", apply(cbind(m[-1,
                                                    1], m[-1, 1] - T_mult * m[-1, 2], m[-1, 1] +
                                                    T_mult * m[-1, 2]), 1, psthr))
        pvalue <- c("", sapply(m[-1, 4], lpvalue))
        title <- c(x_var_str, "", "", lpvalue(globalpvalue))
      }
      if (length(levelnames) == 2) {
        body <- cbind(levelnames, hazardratio, c("",
                                                 ""), c("", ""))
      } else {
        body <- cbind(levelnames, hazardratio, pvalue,
                      rep("", length(levelnames)))
      }
      out <- rbind(title, body)
      if (showN){
        n_by_level = c(nrow(testData),
                       as.vector(table(testData[[x_var]])))
        out <- cbind(out,n_by_level)
      }
      rownames(out) <- NULL
      colnames(out) <- NULL
      return(list(out, nrow(out)))
    } else {
      x_var_str <- lbld(sanitizestr(nicename(x_var)))
      if (type == "coxph") {
        m2 <- coxph(as.formula(paste(paste("Surv(", response[1],
                                           ",", response[2], ")", sep = ""), "~", x_var,
                                     ifelse(strata == "", "", "+"), paste(strata,
                                                                          collapse = "+"), sep = "")), data = data)
        out <- matrix(c(x_var_str, psthr(summary(m2,conf.int=CIwidth)$conf.int[,
                                                              c(1, 3, 4)]), "", lpvalue(summary(m2)$waldtest[3])),
                      ncol = 4)

      } else if (type == "crr") {
        m2 <- crrRx(as.formula(paste(paste(response,
                                           collapse = "+"), "~", x_var, sep = "")), data = data)
        globalpvalue <- try(aod::wald.test(b = m2$coef, Sigma = m2$var,
                                      Terms = seq_len(length(m2$coef)))$result$chi2[3])
        if (class(globalpvalue) == "try-error")
          globalpvalue <- "NA"
        out <- matrix(c(x_var_str, psthr(summary(m2,conf.int=CIwidth)$conf.int[,
                                                              c(1, 3, 4)]), "", lpvalue(globalpvalue)), ncol = 4)
      } else if (type == "logistic") {
        m2 <- glm(as.formula(paste(response, "~", x_var,
                                   sep = "")), family = "binomial", data = data)
        m <- summary(m2)$coefficients
        globalpvalue <- try(aod::wald.test(b = m2$coefficients[-1],
                                      Sigma = vcov(m2)[-1, -1], Terms = seq_len(length(m2$coefficients[-1])))$result$chi2[3])
        if (class(globalpvalue) == "try-error")
          globalpvalue <- "NA"
        Z_mult = qnorm(1-(1-CIwidth)/2)
        out <- matrix(c(x_var_str, psthr(c(exp(m[-1, 1]), exp(m[-1,
                                                                1] - Z_mult * m[-1, 2]), exp(m[-1, 1] + Z_mult *
                                                                                               m[-1, 2]))), "", lpvalue(globalpvalue)), ncol = 4)
      } else if (type == "linear" | type == "boxcox") {
        if (type == "linear") {
          m2 <- lm(as.formula(paste(response, "~", x_var,
                                    sep = "")), data = data)
        } else {
          m2 <- boxcoxfitRx(as.formula(paste(response,
                                             "~", x_var, sep = "")), data = data)
        }
        globalpvalue <- try(aod::wald.test(b = m2$coefficients[-1],
                                      Sigma = vcov(m2)[-1, -1], Terms = seq_len(length(m2$coefficients[-1])))$result$chi2[3])
        if (class(globalpvalue) == "try-error")
          globalpvalue <- "NA"
        m <- summary(m2)$coefficients
        T_mult = qt(1-(1-CIwidth)/2,m2$df.residual)
        out <- matrix(c(x_var_str, psthr(c(m[-1, 1], m[-1,
                                                       1] - T_mult * m[-1, 2], m[-1, 1] + T_mult * m[-1,
                                                                                                     2])), "", lpvalue(globalpvalue)), ncol = 4)
      }
      if (showN){
        out<- cbind(out,nrow(testData))
      }

      return(list(out, nrow(out)))
    }
  })
  table <- lapply(out, function(x) {
    return(x[[1]])
  })
  table <- do.call("rbind", lapply(table, data.frame, stringsAsFactors = FALSE))
  if (showN){
    colnames(table) <- sapply(c("Covariate", sanitizestr(beta),
                                "p-value", "Global p-value","N"), lbld)

  } else{
    colnames(table) <- sapply(c("Covariate", sanitizestr(beta),
                                "p-value", "Global p-value"), lbld)

  }
  return(table)
}

# This function is unchanged, but will call matchcovariate instead of matchcovariate
# matchcovariate fixes a minor bug, code is in helper_functions.R
# change all ordered factors to factors
# expnt is only checked if type = 'glm'. This can be set to FALSE to show raw, instead of exponentiated estimates
# 5 Jan 2021 - replaced normal distribution with t distribution for CI from lm,lme
# TODO: Add support for svyglm and svycoxph functions. May need to check the feasibility of the global p-value here
#'Get multivariate summary dataframe
#'
#'Returns a dataframe corresponding to a univariate table
#'
#'@param model fitted model object
#'@param data dataframe containing data
#'@param showN boolean indicating whether sample sizes should be shown
#'@param markup boolean indicating if you want latex markup
#'@param sanitize boolean indicating if you want to sanitize all strings to not break LaTeX
#'@param nicenames boolean indicating if you want to replace . and _ in strings with a space
#'@param CIwidth level of significance for computing the confidence intervals, default is 0.95
#'@param expnt defaults to null but can be set to false if you want unexponentiated estimates from a log or logit model
#'@keywords dataframe
#'@export
mvsum <-function(model, data, showN = T, markup = T, sanitize = T, nicenames = T,CIwidth=0.95,expnt=NULL)
{
  if (!markup) {
    lbld <- identity
    addspace <- identity
    lpvalue <- identity
  }
  if (!sanitize)
    sanitizestr <- identity
  if (!nicenames)
    nicename <- identity
  if (showN)
    ss_data <- model.frame(model$call,data=data,na.action = 'na.omit')
  call <- paste(deparse(summary(model)$call), collapse = "")
  call <- unlist(strsplit(call, "~", fixed = T))[2]
  call <- unlist(strsplit(call, ",", fixed = T))[1]
  if (substr(call, nchar(call), nchar(call)) == "\"")
    call <- substr(call, 1, nchar(call) - 1)
  call <- unlist(strsplit(call, "\"", fixed = T))[1]
  call <- unlist(strsplit(call, "+", fixed = T))
  call <- unlist(strsplit(call, "*", fixed = T))
  call <- unique(call)
  call <- call[which(is.na(sapply(call, function(cov) {
    charmatch("strata(", cov)
  })) == T)]
  call <- gsub("\\s", "", call)
  type <- class(model)[1]
  if (type=='polr'){
    betanames <- names(model$coefficients)
  } else  if (type != "lme") {
    betanames <- attributes(summary(model)$coef)$dimnames[[1]]
  } else {
    betanames <- names(model$coef$fixed)
  }

  if (type == "glm" ) {
    if (is.null(expnt)) expnt = TRUE
    if (expnt){
      beta <- "OR"
    } else{
      beta <- "Estimate"
    }
    betanames <- betanames[-1]
  } else if (type == "polr" ) {
    if (is.null(expnt)) expnt = TRUE
    beta <- "OR"
  } else if (type == "lm" | type == "lm" | type == "lme") {
    beta <- "Estimate"
    betanames <- betanames[-1]
    expnt = FALSE
  } else if (type == "coxph" | type == "crr") {
    beta <- "HR"
    expnt = FALSE
  } else {
    stop("type must be either polr, coxph, logistic, lm, crr, lme (or NULL)")
  }
  beta = betaWithCI(beta,CIwidth)

  ucall = unique(call)
  indx = matchcovariate(betanames, ucall)
  # My additions to enable tibbles to work
  data = as.data.frame(data)
  if (min(indx) == -1)
    stop("Factor name + level name is the same as another factor name. Please change. Will fix this issue later")
  y <- betaindx(indx)
  if (type %in% c("lm", "glm", "lm", "lme")) {
    y <- lapply(y, function(x) {
      x + 1
    })
    betanames <- c("intercept", betanames)
  }
  out <- lapply(y, function(covariateindex) {
    betaname <- betanames[covariateindex]
    betaname <- strsplit(betaname, ":", fixed = T)
    betaname <- gsub(' - ','-',betaname)  # Added LA, 14 Dec 2020
    betaname <- gsub(' + ','-',betaname)  # Added LA, 14 Dec 2020
    oldcovname <- covnm(betaname[[1]], call)
    oldcovname <- getvarname(oldcovname)
    levelnames <- unlist(lapply(betaname, function(level) {
      paste(mapply(function(lvl, cn) {
        # Changed LA , Dec 15 for level names that also include the varname
        #  result <- unlist(strsplit(lvl, cn, fixed = T))[2]
        result <- unlist(sub(cn,'',lvl))
        out <- ifelse(is.na(result), cn, result)
      }, level, oldcovname), collapse = ":")
    }))
    levelnames <- addspace(sanitizestr(nicename(levelnames)))
    covariatename <- lbld(sanitizestr(nicename(paste(oldcovname,
                                                     collapse = ":"))))
    reference = NULL
    title = NULL
    body = NULL
    if (type == "lme") {
      globalpvalue <- try(aod::wald.test(b = model$coef$fixed[covariateindex],
                                    Sigma = vcov(model)[covariateindex, covariateindex],
                                    Terms = seq_along(covariateindex))$result$chi2[3])
    } else if (type == "polr") {
      globalpvalue <- try(aod::wald.test(b = model$coefficients[covariateindex],
                                    Sigma = vcov(model)[covariateindex, covariateindex],
                                    Terms = seq_along(covariateindex))$result$chi2[3])
    } else if (type != "crr") {
      globalpvalue <- try(aod::wald.test(b = coef(model)[covariateindex],
                                    Sigma = vcov(model)[covariateindex, covariateindex],
                                    Terms = seq_along(covariateindex))$result$chi2[3])
    } else {
      globalpvalue <- try(aod::wald.test(b = model$coef[covariateindex],
                                    Sigma = model$var[covariateindex, covariateindex],
                                    Terms = seq_along(covariateindex))$result$chi2[3])
    }
    if (class(globalpvalue) == "try-error")
      globalpvalue <- "NA"
    globalpvalue <- lpvalue(globalpvalue)
    if (type == "coxph" | type == "crr") {
      hazardratio <- c(apply(matrix(summary(model,conf.int=CIwidth)$conf.int[covariateindex,
                                                            c(1, 3, 4)], ncol = 3), 1, psthr))
      pvalues <- c(sapply(summary(model)$coef[covariateindex,
                                              5], lpvalue))
    } else if (type == "glm" & expnt) {
      m <- summary(model,conf.int=CIwidth)$coefficients
      Z_mult = qnorm(1-(1-CIwidth)/2)
      hazardratio <- apply(cbind(exp(m[covariateindex,1]), exp(m[covariateindex, 1] - Z_mult * m[covariateindex,2]), exp(m[covariateindex, 1] + Z_mult * m[covariateindex,2])), 1, psthr)
      pvalues <- c(sapply(m[covariateindex, 4], lpvalue))
    } else if (type == "polr" & expnt) {
      m <- summary(model)$coefficients
      Z_mult = qnorm(1-(1-CIwidth)/2)
      hazardratio <- apply(cbind(exp(m[covariateindex,1]), exp(m[covariateindex, 1] - Z_mult * m[covariateindex,2]), exp(m[covariateindex, 1] + Z_mult * m[covariateindex,2])), 1, psthr)
      pvalues =  pnorm(abs(m[covariateindex, "Value"]/m[covariateindex, "Std. Error"]),lower.tail = FALSE) * 2
      pvalues <- c(sapply(pvalues, lpvalue))
    } else if (type == "lm" | type == "glm" & !expnt) {
      T_mult = abs(qt((1-CIwidth)/2,model$df.residual))
      m <- summary(model,conf.int=CIwidth)$coefficients
      hazardratio <- apply(cbind(m[covariateindex, 1],m[covariateindex, 1] - T_mult * m[covariateindex,2], m[covariateindex, 1] + T_mult * m[covariateindex,2]), 1, psthr)
      pvalues <- sapply(m[covariateindex, 4], lpvalue)
    } else if (type == "lme") {
      T_mult = abs(qt((1-CIwidth)/2,model$df.residual))
      m <- summary(model,conf.int=CIwidth)$tTable
      hazardratio <- apply(cbind(m[covariateindex, 1],m[covariateindex, 1] - T_mult * m[covariateindex,2], m[covariateindex, 1] + T_mult * m[covariateindex,2]), 1, psthr)
      pvalues <- c(sapply(m[covariateindex, 5], lpvalue))
    }
    if (length(betaname[[1]]) == 1) {
      if (!is.factor(data[, oldcovname])) {
        title <- c(nicename(covariatename), hazardratio,
                   "", globalpvalue)
      } else if (length(levelnames) == 1) {
        title <- c(covariatename, "", "", globalpvalue)
        if (!is.null(data))
          reference <- c(addspace(sanitizestr(names(table(data[,
                                                               which(names(data) == oldcovname)]))[1])),
                         "reference", "", "")
        body <- c(levelnames, hazardratio, "", "")
      } else {
        if (!is.null(data))
          # LisaA 4 Feb, replaced sanitizestr with nicename so that sample sizes could be calculated
          reference <- c(addspace(nicename(names(table(data[,
                                                            which(names(data) == oldcovname)]))[1])),
                         "reference", "", "")
        # reference <- c(addspace(sanitizestr(names(table(data[,
        #                                                        which(names(data) == oldcovname)]))[1])),
        #                  "reference", "", "")
        title <- c(covariatename, "", "", globalpvalue)
        body <- cbind(levelnames, hazardratio, pvalues,
                      rep("", length(levelnames)))
      }
    } else {
      if (length(levelnames) != 1) {
        title <- c(covariatename, "", "", globalpvalue)
        body <- cbind(levelnames, hazardratio, pvalues,
                      rep("", length(levelnames)))
      } else {
        title <- c(covariatename, hazardratio, "", globalpvalue)
      }
    }
    out <- rbind(title, reference, body)
    rownames(out) <- NULL
    colnames(out) <- NULL
    return(list(out, nrow(out)))
  })
  table <- lapply(out, function(x) {
    return(x[[1]])
  })
  index <- unlist(lapply(out, function(x) {
    return(x[[2]])
  }))
  # My new code for sample size calculation
  if (showN ){
    for (i in seq_along(table)){
      tbl <- table[[i]]
      varname <- getvarname(ucall[i])
      level_lookup =  data.frame(nameinDB = unique(ss_data[[varname]]),
                                 nicename = nicename(unique(ss_data[[varname]])))
      N <- nrow(ss_data)
      if (nrow(tbl) > 1){
        for (j in 2:nrow(tbl)){
          dbname = level_lookup$nameinDB[which(level_lookup$nicename==tbl[j,1])]
          if (length(dbname)==0) dbname=tbl[j,1]
          N <- c(N,sum(ss_data[[varname]]==dbname))
        }
      }
      table[[i]] <- cbind(tbl,N)
    }
  }

  table <- do.call("rbind", lapply(table, data.frame, stringsAsFactors = FALSE))
  colnames(table) <- sapply(c("Covariate", sanitizestr(beta),
                              "p-value", "Global p-value"), lbld)
  if (showN) colnames(table)[length(colnames(table))] <- 'N'
  return(table)
}


#' Fit and format an ordinal logistic regression using polr from the {MASS} package. The parallel regression assumption can
#' be tested using the Brant test in the Brant package and visually. Only logistic ordinal regression is supported currently.
#' # Options to perform univariate or multivariate analysis.
#'@param data dataframe containing data [REQUIRED]
#'@param covs character vector with the names of columns to include in table [REQUIRED]
#'@param response ordinal outcome variable [REQUIRED]
#'@param reflevel manual specification of the reference level, must match level exactly
#'@param markup boolean indicating if you want latex markup
#'@param sanitize boolean indicating if you want to sanitize all strings to not break LaTeX
#'@param nicenames booling indicating if you want to replace . and _ in strings with a space
#'@param mv logical indicating whether you want to summarise univariate analyses or a single model
#'@param excludeLevels a named list of levels to exclude from factor variables. Currently, this has only been implemented for the response variable.
#'@param testPO logical, should the proportional odds (parallel regression) assumption be tested with the Brant test, defaults to TRUE
#'@param showN logical, should sample sizes be shown for each lvel, defaults to TRUE
#'@param digits number of digits to display, defaults to
#'@param CIwidth level of significance for computing the confidence intervals, default is 0.95
#'#'@return A formatted table displaying the odds ratio associated with each covariate
#'@example
#'ordsum(data,covs=c('Age','Sex'),maincov='SES Level')
#'@keywords ordinal regression, Brant test
#' @importFrom xtable xtable
#' @importFrom MASS polr
#'@export
#'
ordsum    <- function(data, covs, response,reflevel='NULL',markup=FALSE,sanitize=TRUE,nicenames=TRUE,
                      mv=FALSE,excludeLevels=NULL,testPO=TRUE,showN=TRUE,digits=2,CIwidth=0.95){

  if (!markup) {
    lbld <- identity # not yet used
    addspace <- identity  # not yet used
    lpvalue <- identity
  }
  if (!sanitize)
    sanitizestr <- identity
  if (!nicenames)
    nicename <- identity

  missing_covs = setdiff(covs,names(data))
  if (length(missing_covs)>0) {
    stop(paste('Check the covarariates, the following variables are not in the data:',missing_covs))
  }
  if (!class(data[[response]])[1] %in% c('factor','ordered')) {
    warning('Response variable is not a factor, will be converted to an ordered factor')
    data[[response]] <- factor(data[[response]],ordered=T)
  }
  if (!markup) {
    lbld <- identity
    addspace <- identity
    lpvalue <- identity
  }
  if (!sanitize)
    sanitizestr <- identity
  if (!nicenames)
    nicename <- identity
  if (!is.null(excludeLevels)){
    if (response %in% names(excludeLevels)){
      to_remove = sapply(data[[response]],function (x) {x %in% excludeLevels[[response]]})
      data = data[!to_remove,]
    }
  }
  if (!(mv)){
    out <- lapply(covs, function(x_var) {
      if (!is.null(excludeLevels[[x_var]])){
        excluded_xvar = excludeLevels[[x_var]]
      } else {
        excluded_xvar = NULL
      }
      excluded_xvar <- c(excluded_xvar,names(table(data[[x_var]]))[table(data[[x_var]])==0])
      if (length(excluded_xvar)>0){
        testData = data[!(data[[x_var]] %in% excluded_xvar),]
      } else{
        testData = data
      }

      # if (class(testData[[x_var]])[1]=='ordered'){
      #   testData[[x_var]] <- factor(testData[[x_var]],ordered = FALSE)
      # }
      polr.form = as.formula(paste(response,'~',x_var))
      fit = MASS::polr(data=testData,polr.form,method='logistic',Hess = TRUE)
      brant_test = try(modified_brant(fit,by.var=T),silent = T)
      if (class(brant_test)[1] == "try-error") {
        po_test_omni = data.frame(Covariate = x_var,"PO Test" = 'Not Tested')
        zero_counts <- NA
      } else {
        po_test_omni = data.frame(Covariate=rownames(brant_test$result),brant_test$result)
        po_test_omni$"PO Test" = formatp(po_test_omni$probability)
        zero_counts <-brant_test$zero_count_cells
      }
      coef_tbl <- data.frame(summary(fit)$coef)
      coef_tbl <- coef_tbl[grep(x_var,rownames(summary(fit)$coef)),]
      coef_tbl$p_value <- pnorm(abs(coef_tbl$"t.value"), lower.tail = FALSE) * 2
      coef_tbl$"p-value" = sapply(coef_tbl$p_value,lpvalue)
      coef_tbl$OR <- exp(coef_tbl$"Value")
      qstdN <- qnorm(p=1-(1-CIwidth)/2)
      coef_tbl$LB <- exp(coef_tbl$"Value"-qstdN*coef_tbl$"Std..Error")
      coef_tbl$UB <- exp(coef_tbl$"Value"+qstdN*coef_tbl$"Std..Error")
      coef_tbl$"OR_CI" = paste0(niceNum(coef_tbl$OR,digits = digits),
                                ' (',niceNum(coef_tbl$LB,digits = digits),
                                ',',niceNum(coef_tbl$UB,digits = digits),')')
      tbl <- cbind(Covariate=rownames(coef_tbl),data.frame(coef_tbl[,c(9,5)]))

      if (class(testData[[x_var]])[1] %in% c("ordered", "factor" )){
        brant_test_level = try(modified_brant(model=fit,by.var=F),silent = T)
        if (class(brant_test_level)[1] == "try-error") {
          po_test_level = data.frame(Covariate = tbl,"PO Test" = 'NT')
        } else {
          po_test_level = data.frame(cbind(brant_test_level$result,Covariate=rownames(brant_test_level$result)))
          po_test_level$"PO Test" <- formatp(po_test_level$probability)
          #          po_test_level$"PO Test"[po_test_level$"PO Test"=='1'] <- 'NT' # WHY DID I ADD THIS?
        }

        tbl$"PO Test" = sapply(tbl$Covariate, function (x) po_test_level$"PO Test"[po_test_level$Covariate==x])
        tbl$Covariate = sub(x_var,'',tbl$Covariate)
        reflevel=setdiff(levels(testData[[x_var]]),c(excluded_xvar,tbl$Covariate))
        tbl <- rbind(c(reflevel,"Reference","",""),tbl)
        nterms=length(fit$coefficients)
        globalp <- try(aod::wald.test(Sigma=vcov(fit)[1:nterms,1:nterms],
                                      b=fit$coefficients,Terms=1:nterms)$result$chi2[3],silent = T)
        if (class(globalp)[1]=='try-error') globalp <- NA
        tbl$"globalPval" = ''
        tbl <- rbind(c(x_var,"","",
                       po_test_omni$"PO Test"[po_test_omni$Covariate==x_var],
                       lpvalue(globalp)),tbl)

        n_by_level <- as.vector(sapply(tbl$Covariate[-1],FUN = function(x){ sum(fit$model[[x_var]]==x)}))
        n_by_level <- c(sum(n_by_level),n_by_level)
        tbl <- cbind(tbl,N=n_by_level)


      } else {
        tbl$"PO Test" = po_test_omni$"PO Test"[po_test_omni$Covariate==x_var]
        tbl$"globalPval" = tbl$p.value
        tbl$p.value <- ''

        tbl <- cbind(tbl,N=nrow(fit$fitted.values))

      }
      tbl <- cbind(tbl,ZeroCount=c(zero_counts,rep(NA,nrow(tbl)-1)))
      tbl <- tbl[,c(1,6,2,5,3,4,7)]
      beta = betaWithCI('OR',CIwidth)
      names(tbl) <- c("Covariate","N",sanitizestr(beta),"Global p-value","p-value","PO Test","ZeroCount")
      return(tbl)
    })
    onetbl = do.call("rbind",out)
    onetbl$Covariate <- nicename(onetbl$Covariate)

    if (sum(onetbl$ZeroCount>0,na.rm=T)==0){
      outtext <- "The proportional odds assumption is reasonable for all variables."
    } else {
      caution_star = rep('',nrow(onetbl))
      caution_star[onetbl$ZeroCount>0] <- '*'
      onetbl$`PO Test` <-  paste0(onetbl$`PO Test`,caution_star)
      outtext <- "*Caution, proportional odds test performed with empty cells."
    }
  } # this is the end of the univariate analysis

  onetbl <- onetbl[,-7]
  if (sum(onetbl$`p-value`=='')==nrow(onetbl)) {
    onetbl <- onetbl[,-4]
  }
  if (!showN) {
    onetbl <- onetbl[,-2]
  }

  return(list(onetbl,outtext))
}




#' Plot univariate relationships all on one plot
#' Need a hierarchy of how to display plots sensibly
#' If response is continuous
#'   For a numeric predictor -> scatterplot
#'   For a categorical predictor -> boxplot
#' If response is a factor
#'   For a numeric predictor -> boxplot
#'   For a categorical predictor -> barplot
#' @param response character vector with names of columns to use for response
#' @param covs character vector with names of columns to use for covariates
#' @param data dataframe containing your data
#' @param showN boolean indicating whether sample sizes should be shown on the plots
#' @param na.rm boolean indicating whether na values should be shown or removed
#' @param response_title character value with title of the plot
#' @keywords plot
#' @import ggplot2
#' @importFrom ggpubr ggarrange
#' @export
plot_univariate <- function(response,covs,data,showN=FALSE,na.rm=TRUE,response_title=NULL){
  # if (!class(data[[response]])[1] %in% c('factor','ordered','numeric')) {
  #   stop('Response variable must be numeric or factor.')
  # }
  # bad_covs = sapply(covs,function(x) !class(data[[x]])[1] %in% c('factor','ordered','numeric'))
  # if (sum(bad_covs)>0){
  #   stop(paste('The following variables are neither numeric nor factors and can not be covariates:',c(covs[bad_covs])))
  # }

  for (v in c(response,covs)){
    if (class(data[[v]])=='character') data[[v]] <- factor(data[[v]])
  }

  if (is.null(response_title)) response_title = response
  response_title = niceStr(response_title)
  plist <- NULL
  if (class(data[[response]])[1] %in% c('factor','ordered')){
    levels(data[[response]]) = niceStr(levels(data[[response]]))
    for (x_var in covs){
      # remove missing data, if requested
      if (na.rm) pdata = na.omit(data[,c(response,x_var)]) else pdata = data[,c(response,x_var)]

      if (class(pdata[[x_var]])[1] =='numeric' ){
        p <- ggplot(data=pdata, aes_string(y=response,x=x_var,fill=response)) +
          geom_boxplot()
        if (showN){
          p=  p+
            stat_summary(geom='text',fun.data = lbl_count,vjust=-0.5,hjust=1)
        }
      } else {
        p <- ggplot(data=pdata, aes_string(x=x_var,fill=response)) +
          geom_bar(position='fill') +
          scale_x_discrete(labels= function(x) wrp_lbl(x))
        if (showN){
          p <- p +
            geom_text(aes(label=stat(count)),stat='count',position='fill',vjust=1)
        }
        if (length(unique(pdata[[x_var]]))>8){
          p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
      }
      plist[[x_var]] <- p  +
        theme(axis.text.y=element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = 'bottom',
              plot.title = element_text(size=10),
              plot.margin = unit(c(0,1,0,1), "lines")) +
        labs(title=niceStr(x_var),x='',y='',fill=response_title)
    }
  } else{
    for (x_var in covs){
      # remove missing data, if requested
      if (na.rm) pdata = na.omit(data[,c(response,x_var)]) else pdata = data[,c(response,x_var)]

      if (class(pdata[[x_var]])[1] =='numeric' ){
        p <- ggplot(data=pdata, aes_string(y=response,x=x_var)) +
          geom_point()
      } else
        p <- ggplot(data=pdata, aes_string(y=response,x=x_var,fill=response)) +
          geom_boxplot() +
          scale_x_discrete(labels= function(x) wrp_lbl(x))
      if (showN){
        p=  p+
          stat_summary(geom='text',fun.data = lbl_count,vjust=-0.5,hjust=1)
      }
      if (length(unique(pdata[[x_var]]))>8){
        p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
      plist[[x_var]] <- p  +
        theme(
          legend.position = 'bottom',
          plot.title = element_text(size=10),
          plot.margin = unit(c(0,1,0,1), "lines")) +
        labs(title=niceStr(x_var),x='',y='',fill=response_title)
    }

  }
  ggarrange(plotlist=plist,
            common.legend = T,
            ncol=2,
            nrow=ceiling(length(plist)/2))
}


#'Outputs a ggplot to test whether numeric predictors are linear in logit
#'
#'@param model A glm model object to test
#'@export
plot_lr_check <- function(model){

  if ( class(model)[1]!='glm' | model$family$link !='logit') stop('This function requires a glm logit link model.')
  data = model$model

  predictors <- colnames(data)[-1]
  numeric_vars = data[,predictors]
  numeric_vars = dplyr::select_if(numeric_vars,is.numeric)

  if (ncol(numeric_vars)==0) stop('No numeric predictors to check.')

  numeric_predictors <- colnames(numeric_vars)
  response = names(data)[1]
  probs = predict(model,type='response')
  data$logits = log(probs/(1-probs))

  data <- data[,c(response,'logits',numeric_predictors)]

  pd = reshape::melt(data,id.vars = 1:2)

  ggplot(pd, aes(logits,value)) +
    geom_point(aes_string(col=response)) +
    geom_smooth(method='lm') +
    geom_smooth(col='red',se=F) +
    facet_wrap(~variable,scales='free_y')
}

#' This will provide a graph to check the PO assumption, according to
#' Harrell, RMS, pg 337
#'@param polr.fit A polr model object to test
#'@export
plot_po_check <- function(polr.fit){
  fit=polr.fit

  data = fit$model
  response = names(data)[1]
  responseOptions = levels(data[[response]])
  data$y = as.numeric(data[[response]])-1

  # Build the summary function for any response
  sfStr <- NULL
  for (i in 2:length(responseOptions)){
    entry = paste0("'Y>=",responseOptions[i],"'",'=qlogis(mean(y>=',i-1,'))')
    sfStr  = c(sfStr,entry)
  }
  sfStr  = paste0(sfStr,collapse=',')
  sfStr  = paste0('sf <- function(y) {c(',sfStr,')}',collapse='')
  eval(parse(text=sfStr))

  f = paste('y~',paste(attr(fit$terms,'term.labels'),collapse = '+'))
  s <- summary(as.formula(f),data=data,fun=sf)
  # remove any infinite values - users
  if (sum(is.infinite(s))>0) subtitle = '***Infinite estimates removed***' else subtitle=''
  s[is.infinite(s)] <- NA
  plot(s,which=1:2,pch=1:2,xlab='logit',vnames='names',main=response, width.factor=1.5)
  mtext(subtitle)
}

# # To Do: Improve flow so that this only needs to be called once per knit
# This functionality does not work once the functions are inside a package library
# getOutputFormat <- function() {
#   # This function, if run while a document is being knit
#   # will return the output format, for example:
#   # "bookdown::word_document2"
#
#   # First check for format override
#   if (exists('outputFormatOverride')){
#     return(outputFormatOverride)
#   }
#   output <- try(opts_knit$get("rmarkdown.pandoc.to"),silent=T)
#   print(output)
#   if (class(output)[1]=='try-error' | is.null(output)){
#     # This happens when running from a chunk
#     return('pdf')
#   }
#   if (length(grep('doc',output,ignore.case = T))>0) {
#     return('docx')
#   } else {
#     return('pdf')
#   }
# }

#' Format any generic table like the reportRx tables.
#' Table output defaults to kable, but the kableExtra package doesn't work well with Word.
#' To export nice tables to Word use options('doc_type'='doc')
#' @param tab a table to format, it assumes the first column contains labels to tidy
#' @param rmstr characters to replace with spaces in the output, defaults to . and _
#' @param digits number of digits to round numeric columns to
#' @param to_indent numeric vector the length of nrow(tab) indicating which rows to indent
#' @param to_bold numeric vector the length of nrow(tab) indicating which rows to bold
#' @param caption table caption
#' @param chunk_label only used if out_fmt = doc to allow cross-referencing
#' @param ... other variables passed to covsum and the table output function
#' @export
niceTable <- function(tab,rmstr = "[._]",digits=2,to_indent=numeric(0),to_bold=numeric(0),caption=NULL,chunk_label,...){

  out_fmt = ifelse(is.null(getOption('doc_type')),'pdf',getOption('doc_type'))
  out_fmt =ifelse(out_fmt%in%c('doc','docx'),'doc','pdf')

# Replace '.' and '_' in column names with spaces
  colnames(tab) = sapply(colnames(tab), function(x) gsub(rmstr," ",x),simplify = T)
  tab[,1] = sapply(tab[,1], function(x) gsub(rmstr," ",x),simplify = T)

    # Check for optional arguments
  chunk_label = ifelse(missing(chunk_label),'tbl',chunk_label)

  # set NA to empty in pander and kable
  options(knitr.kable.NA = '')

  if (is.null(to_indent)) to_indent = numeric(0)
  to_indent = as.vector(to_indent)
  rownames(tab) <- NULL


  if (out_fmt=='doc'){
    pander::panderOptions('round', digits)
    caption = paste0('(\\#tab:',chunk_label,')',caption)
    tab[is.na(tab)] <-'&nbsp;' # This is necessary to assign the 'Compact' style to empty cells
    tab[tab==''] <-'&nbsp;'

    tab[[1]][to_indent] <- sapply(tab[[1]][to_indent],function(x) paste('&nbsp;&nbsp;',x))
    if ( length(to_bold)>0) {
      pander::pander(tab,
                     caption=caption,
                     emphasize.strong.rows=to_bold,
                     split.table=Inf, split.cells=15,
                     justify = paste(c('l',rep('r',ncol(tab)-1)),collapse = '',sep=''))

    } else {
      pander::pander(tab,
                     caption=caption,
                     split.table=Inf, split.cells=15,
                     justify = paste(c('l',rep('r',ncol(tab)-1)),collapse = '',sep=''))
    }
  } else {
    if (nrow(tab)>30){
      kout <- knitr::kable(tab, booktabs=TRUE,
                           longtable=TRUE,
                           linesep='',
                           digits = digits,
                           caption=caption,
                           align = paste(c('l',rep('r',ncol(tab)-1)),collapse = '',sep=''))
      if (ncol(tab)>4) {
        kout <- kableExtra::kable_styling(kout,full_width = T,latex_options = c('repeat_header'))
      } else {
        kout <- kableExtra::kable_styling(kout,latex_options = c('repeat_header'))
      }
    } else {
      kout <- knitr::kable(tab, booktabs=TRUE,
                           longtable=FALSE,
                           linesep='',
                           digits = digits,
                           caption=caption,
                           align = paste(c('l',rep('r',ncol(tab)-1)),collapse = '',sep=''))
      if (ncol(tab)>4) kout <- kableExtra::kable_styling(kout,full_width = T)
    }
    kout <- kableExtra::add_indent(kout,positions = to_indent)
    if (length(to_bold)>0){
      kout<- kableExtra::row_spec(kout,to_bold,bold=TRUE)
    }
    kout
  }

}

#' The output function for the print methods
#' Table output defaults to kable, but the kableExtra package doesn't work well with Word.
#' To export nice tables to Word use options('doc_type'='doc')
#' @param tab a table to format
#' @param to_indent numeric vector the length of nrow(tab) indicating which rows to indent
#' @param to_bold numeric vector the length of nrow(tab) indicating which rows to bold
#' @param caption table caption
#' @param chunk_label only used if out_fmt = doc to allow cross-referencing
#' @param ... other variables passed to covsum and the table output function
#' @export
outTable <- function(tab,to_indent=numeric(0),to_bold=numeric(0),caption=NULL,chunk_label,...){

  out_fmt = ifelse(is.null(getOption('doc_type')),'pdf',getOption('doc_type'))
  out_fmt =ifelse(out_fmt%in%c('doc','docx'),'doc','pdf')
  chunk_label = ifelse(missing(chunk_label),'NOLABELTOADD',chunk_label)

  # set NA to empty in pander and kable
  options(knitr.kable.NA = '')

  if (is.null(to_indent)) to_indent = numeric(0)
  to_indent = as.vector(to_indent)
  rownames(tab) <- NULL


  if (out_fmt=='doc'){
    caption = ifelse(chunk_label=='NOLABELTOADD',caption,paste0('(\\#tab:',chunk_label,')',caption))
    tab[is.na(tab)] <-'&nbsp;' # This is necessary to assign the 'Compact' style to empty cells
    tab[tab==''] <-'&nbsp;'

    tab[[1]][to_indent] <- sapply(tab[[1]][to_indent],function(x) paste('&nbsp;&nbsp;',x))
    if (length(to_bold)>0) {
      pander::pander(tab,
             caption=caption,
             emphasize.strong.rows=to_bold,
             split.table=Inf, split.cells=15,
             justify = paste(c('l',rep('r',ncol(tab)-1)),collapse = '',sep=''))

    } else {
      pander::pander(tab,
             caption=caption,
             split.table=Inf, split.cells=15,
             justify = paste(c('l',rep('r',ncol(tab)-1)),collapse = '',sep=''))
    }
  } else {
    if (nrow(tab)>30){
      kout <- knitr::kable(tab, booktabs=TRUE,
                           longtable=TRUE,
                           linesep='',
                           caption=caption,
                           align = paste(c('l',rep('r',ncol(tab)-1)),collapse = '',sep=''))
      if (ncol(tab)>4) {
        kout <- kableExtra::kable_styling(kout,full_width = T,latex_options = c('repeat_header'))
      } else {
        kout <- kableExtra::kable_styling(kout,latex_options = c('repeat_header'))
      }
    } else {
      kout <- knitr::kable(tab, booktabs=TRUE,
                           longtable=FALSE,
                           linesep='',
                           caption=caption,
                           align = paste(c('l',rep('r',ncol(tab)-1)),collapse = '',sep=''))
      if (ncol(tab)>4) kout <- kableExtra::kable_styling(kout,full_width = T)
    }
    kout <- kableExtra::add_indent(kout,positions = to_indent)
    if (length(to_bold)>0){
      kout<- kableExtra::row_spec(kout,to_bold,bold=TRUE)
    }
    kout
  }

}

#' Nicely display a confusion matrix from the caret package
#' @param confMtrx a confusion matrix from the caret package
#' @param caption table caption
#' @export
printConfMatrix <- function(confMtrx,caption=NULL){

  if (class(confMtrx) != "confusionMatrix") stop('This function tidies output from a confusionMatrix from the caret package.')
  t <- confMtrx$table
  tbl <- data.frame(Prediction = dimnames(t)$Prediction,
                    col2 = t[,1],
                    col3 = t[,2])
  names(tbl)[2:3] = paste('Reference',dimnames(t)$Reference)
  if (is.null(caption)) caption = 'Agreement between observed and predicted values.'
  outTable(tbl,caption=caption)

}

#' Nicely display statistics from a confusion matrix from the caret package
#' @param confMtrx a confusion matrix from the caret package
#' @param digits number of digitis to display
#' @param what which stats to show, see caret package for details
#' @param caption table caption
#' @param tblOnly should a dataframe or a formatted object be returned
#' @export
printConfStats <- function(confMtrx,digits=2,what=c(1,2,3,4,12),caption=NULL,tblOnly=FALSE){

  if (class(confMtrx) != "confusionMatrix") stop('This function tidies output from a confusionMatrix from the caret package.')
  if (max(what)>18) stop('The what argument needs to specify indices of possible performance statistics. Set what=NULL for full list.')
  t <- confMtrx$byClass
  tbl1 = tibble(Statistic = names(t),
                Value = round(t,digits))
  t = confMtrx$overall
  tbl2 = tibble(Statistic = names(t),
                Value = round(t,digits))
  if (is.null(what)){
    tab <- rbind(tbl1,tbl2)
  } else {
    tab <- rbind(tbl1,tbl2)[what,]
  }
  if (tblOnly) return(tab)
  if (is.null(caption)) caption = 'Performance statistics for the classification model.'
  outTable(tab,caption=caption)

}



#' Fit and format an ordinal logistic regression using polr from the {MASS} package. The parallel regression assumption can
#' be tested using the Brant test in the Brant package and visually. Only logistic ordinal regression is supported currently.
#'@param data dataframe containing data [REQUIRED]
#'@param covs character vector with the names of columns to include in table [REQUIRED]
#'@param response ordinal outcome variable [REQUIRED]
#'@param reflevel manual specification of the reference level, must match level exactly
#'@param caption Table caption
#'@param showN logical, should sample sizes be shown for each lvel, defaults to TRUE
#'@param mv logical indicating whether you want to summarise univariate analyses or a single model only mv=FALSE currently supported
#'@param excludeLevels a named list of levels to exclude from factor variables. Currently, this has only been implemented for the response variable.
#'@param testPO logical, should the proportional odds (parallel regression) assumption be tested with the Brant test, defaults to TRUE
#'@param digits number of digits to display, defaults to
#'@param CIwidth level of significance for computing the confidence intervals, default is 0.95
#'@return A formatted table displaying the odds ratio associated with each covariate
#'@example
#'printOrdsum(data,covs=c('Age','Sex'),maincov='SES Level')
#'@keywords ordinal regression, Brant test
#'@export
#'
printOrdsum <- function(data, covs, response, reflevel='NULL', caption = NULL, showN=T,
                        mv=FALSE,excludeLevels=NULL,testPO=TRUE,digits=2,CIwidth=0.95){


  rtn <- ordsum(data=data,
                covs=covs,
                response=response,
                reflevel=reflevel,
                markup=FALSE,
                sanitize=FALSE,
                nicenames=T,
                mv=FALSE,
                excludeLevels=excludeLevels,
                testPO=testPO,
                showN=showN,
                digits = digits,
                CIwidth=CIwidth)
  tab <- rtn[[1]]
  if (mv){
    type ='Multivariable'
  } else type='Univariate'

  if (is.null(caption))
    caption = paste(type,'ordinal logistic regression analysis of predictors of',nicename(response),'.')

  # format p-values nicely
  tab$`Global p-value` <- sapply(tab$`Global p-value`,formatp)
  if (length(which(names(tab)=='p-value'))>0)
    tab[,which(names(tab)=='p-value')] <- sapply(tab[[which(names(tab)=='p-value')]],formatp)

  nice_var_names = gsub('_',' ',covs)
  to_indent <- which(!tab$Covariate %in% nice_var_names )
  to_bold <- which(as.numeric(tab[["Global p-value"]])<(1-CIwidth))

  outTable(tab=tab,to_indent=to_indent,to_bold=to_bold,
           caption=caption)

}

#'
#'Returns a dataframe corresponding to a descriptive table
#' The default output is a kable table for use in pdfs or html, but pander tables can be produced
#' for Word documents by specifying options('doc_type'='doc') in the setup chunk of the markdown document.
#'
#'@param data dataframe containing data
#'@param covs character vector with the names of columns to include in table
#'@param maincov covariate to stratify table by
#'@param caption character containing table caption
#'@param excludeLevels a named list of levels to exclude in the form list(varname =c('level1','level2')) these levels will be excluded from association tests
#'@param showP boolean indicating if p values should be displayed (may only want corrected p-values)
#'@param showIQR boolean indicating if you want to display the inter quantile range (Q1,Q3) as opposed to (min,max) in the summary for continuous variables
#'@param tableOnly should a dataframe or a formatted print object be returned
#'@param covTitle character with the names of the covariate column
#'@param chunk_label only used if out_fmt = doc to allow cross-referencing
#'@param ... other variables passed to covsum and the table output function
#'@keywords dataframe
#'@return A formatted table displaying a summary of the covariates stratified by maincov
#'@export
#'@seealso \code{\link{fisher.test}}, \code{\link{chisq.test}}, \code{\link{wilcox.test}}, \code{\link{kruskal.test}}, and \code{\link{anova}}
printCovsum <- function(data,covs,maincov=NULL,caption=NULL,excludeLevels=NULL,showP=TRUE,showIQR=FALSE,tableOnly=FALSE,covTitle='Covariate',
                        chunk_label,...){

  tab <- covsum(data,covs,maincov,markup=FALSE,excludeLevels=excludeLevels,IQR=showIQR,sanitize=FALSE,...)
  if (is.null(maincov))
    showP=FALSE
  to_bold = numeric(0)
  if (showP) {
    # format p-values nicely
    p_vals <- tab[,ncol(tab)]
    new_p <- sapply(p_vals,formatp)
    tab[,ncol(tab)] <- new_p
    to_bold <- which(as.numeric(new_p)<0.05)
  } else {
    tab <- tab[,!names(tab) %in%'p-value']
  }
  nice_var_names = gsub('[_.]',' ',covs)
  to_indent <- which(!tab$Covariate %in% nice_var_names )
  if (showIQR) {
    tab$Covariate <- gsub('[(]Q1,Q3[)]','(IQR)',tab$Covariate)
  }
  if (covTitle !='Covariate') names(tab[1]) <-covTitle
  if (tableOnly){
    return(tab)
  }
  if (is.null(caption)) {
    if (!is.null(maincov)){
      caption = paste0('Summary sample statistics by ',nicename(maincov),'.')
    } else
      caption = 'Summary sample statistics.'
  }

  outTable(tab=tab,to_indent=to_indent,to_bold=to_bold,
           caption=caption,
           chunk_label=ifelse(missing(chunk_label),'NOLABELTOADD',chunk_label))

}

#' Output several univariate models nicely.
#' The default output is a kable table for use in pdfs or html, but pander tables can be produced
#' for Word documents by specifying options('doc_type'='doc') in the setup chunk of the markdown document.
#'
#' @param response string vector with name of response
#' @param covs character vector with the names of columns to fit univariate models to
#' @param data dataframe containing data
#' @param adj_cov character vector containing covariates to adjust for in each (no longer univariate) model
#' @param caption table caption
#' @param showP boolean indicating if p-values should be shown (may only want to show holm-corrected values)
#' @param showN boolean indicating if sample sizes should be shown
#' @param tableOnly boolean indicating if unformatted table should be returned
#' @param removeInf boolean indicating if infinite estimates should be removed from the table
#' @param HolmGlobalp boolean indicting if a Holm-corrected p-value should be presented
#' @param ... other variables passed to uvsum and the table output functions
#' @export
printUnivariateTable <- function(response, covs , data ,adj_cov=NULL, caption=NULL,showP=T,showN=T,tableOnly=FALSE,removeInf=T,HolmGlobalp=FALSE,...){

  # get the table
  tab <- uvsum(response,covs,data,adj_cov=NULL,markup = FALSE,showN=showN,sanitize=FALSE,...)

  cap_warn <- character(0)
  if (removeInf){
    # Do not display unstable estimates
    if (showN) {
      if (any(tab$N<=1)){
        tab[tab$N<=1,2] <- NA
        cap_warn = paste0('Covariates not analysed: ',paste(tab$Covariate[tab$N<=1],collapse=','),'. ')
      }}
    inf_values =  grep('Inf',tab[,2])
    if (length(inf_values)>0){
      tab[inf_values,2:4] <-NA
      cap_warn <- paste0(cap_warn,'Covariates with unstable estimates:',paste(tab$Covariate[inf_values],collapse=','),'.')

    }
  }

  # If HolmGlobalp = T then report an extra column with the adjusted p and only bold these values
  if (HolmGlobalp){
    p_sig <- p.adjust(tab$`Global p-value`,method='holm')
    tab$"Holm Adj p" = p_sig
  } else {
    p_sig <- tab$`Global p-value`
  }

  to_bold <- which(as.numeric(p_sig)<0.05)
  nice_var_names = gsub('[_.]',' ',covs)
  to_indent <- which(!tab$Covariate %in% nice_var_names )

  tab[["p-value"]] <- formatp(tab[["p-value"]])
  tab[["Global p-value"]] <- formatp(tab[["Global p-value"]])
  if (HolmGlobalp){
    tab[["Holm Adj p"]] <- formatp(tab[["Holm Adj p"]])
  }

  # Reorder if the samples sizes are shown
  if (showN) {
    tab <- tab[,c("Covariate","N",setdiff(colnames(tab),c("Covariate","N")))]
  }

  # If all outcomes are continunous (and so all p-values are NA), remove this column & rename Global p-value to p-value
  if (sum(is.na(tab[["p-value"]]))==nrow(tab)) {
    tab <- tab[,-which(names(tab)=="p-value")]
    names(tab) <- gsub('Global p-value','p-value',names(tab))
  }

  if (tableOnly){
    if (nchar(cap_warn)>0) warning(cap_warn)
    return(tab)
  }
  if(is.null(caption)){
    caption = paste0('Univariate analysis of predictors of ',niceStr(response),'.',cap_warn)
  } else if (caption=='none' & identical(cap_warn,character(0))) {
    caption=NULL
  } else caption = paste0(caption,cap_warn)

  outTable(tab=tab,to_indent=to_indent,to_bold=to_bold,
           caption=caption)
}

# # TO FINISH - Look at Hogan study
# printMissingDataComparison <- function(model,data,id_variable,status_var=NULL,caption=NULL){
#
#   data$excluded = if_else(data[[id_variable]] %in% names(model$residuals),"Included","Missing")
#   if (!is.null(status_var)){
#     data[[status_var]] = factor(if_else(data[[status_var]]==1,"1","0")  )
#   }
#
#   printCovsum(data=data,
#               covs=c('Follow_Up','Deaths',surv_predictors),
#               maincov='os_missing',
#               showP= TRUE)
# }

#' Output a multivariable model nicely
#' The default output is a kable table for use in pdfs or html, but pander tables can be produced
#' for Word documents by specifying options('doc_type'='doc') in the setup chunk of the markdown document.
#'
#' @param model model fit
#' @param data data that model was fit on
#' @param caption table caption
#' @param showN boolean indicating if sample sizes should be shown
#' @param tableOnly boolean indicating if unformatted table should be returned
#' @param HolmGlobalp boolean indicting if a Holm-corrected p-value should be presented
#' @param expnt defaults to NULL can be set to FALSE for log and logit link models
#' @export
printMultivariableTable <- function(model , data ,caption=NULL,showN=T,tableOnly=FALSE,HolmGlobalp=FALSE,expnt=NULL){

  # get the table
  tab <- mvsum(model=model,data=data,showN=showN,markup = FALSE, sanitize = FALSE, nicenames = T,expnt=expnt)


  # Reduce the number of significant digits in p-values
  p_val <-  formatp(tab$`p-value`)
  g_p_val = formatp(tab$`Global p-value`)
  # If HolmGlobalp = T then report an extra column with the adjusted p and only bold these values
  if (HolmGlobalp){
    gp <- p.adjust(tab$`Global p-value`,method='holm')
  } else {
    gp <- tab$`Global p-value`
  }
  to_bold <- which(as.numeric(gp)<0.05)
  to_indent <- which(is.na(g_p_val))

  tab$`p-value` <- p_val
  tab$`Global p-value` <- g_p_val

  # Reorder if the samples sizes are shown
  if (showN) {
    tab <- tab[,c(1,5,2:4)]
  }

  # Reorder
  if (HolmGlobalp){
    tab$`Holm Adj p` <- formatp(gp)
  }
  # TO DO: possibly automate this... need to extract response from mvsum
  # if(is.null(caption)){
  #   caption = paste0('Multivariable analysis of predictors of ',niceStr(response),'.')
  # } else if (caption=='none') {
  #   caption=NULL
  # }


  # If all outcomes are contunous (and so all p-values are NA), remove this column
  if (sum(is.na(tab[["p-value"]]))==nrow(tab)) tab <- tab[,-which(names(tab)=="p-value")]

  if (tableOnly){
    return(tab)
  }
  outTable(tab=tab,to_indent=to_indent,to_bold=to_bold,
           caption=caption)

}

#' Function to nicely display median survival time by group
#' @param data datafrmae
#' @param time survival time
#' @param status indication of censored or observed data
#' @param group the variable to group observations by
#' @param grpLabel character value describing the grouping variable
#' @param digits the number of digitis to display
#' @param caption table caption
#' @param tblOnly should a dataframe or a formatted object be returned
#' @importFrom tibble as_tibble
#' @export
survmedian_sum <- function(data,time,status,group,grpLabel='Group',digits=1,caption=NULL,tblOnly=FALSE){
  if (group==1) stop('This function requires a grouping variable.')
  out_fmt <- knitr::opts_knit$get("rmarkdown.pandoc.to")
  if (is.null(out_fmt))  out_fmt='pdf'

  overall_fit <- survfit(as.formula(paste0("Surv(",time,',',status,') ~1')),
                         data = data)

  otbl <- tibble::as_tibble(t(summary(overall_fit)$table))
  otbl[[grpLabel]] <- 'Overall'
  otbl$Expected <- NA

  fit <- survfit(as.formula(paste0("Surv(",time,',',status,') ~',group)),
                 data = data)
  lr_test = survdiff(as.formula(paste0("Surv(",time,',',status,') ~',group)),
                     data = data)
  grptbl <- tibble::as_tibble(summary(fit)$table,rownames=grpLabel)
  grptbl$Expected = lr_test$exp

  xtbl <- dplyr::bind_rows(otbl,grptbl)
  xtbl[[grpLabel]] <- gsub(paste0(group,"="),"",xtbl[[grpLabel]])

  miss_vals = xtbl$records-xtbl$n.start
  xtbl$CI = paste0(niceNum(xtbl[["0.95LCL"]],digits),', ',niceNum(xtbl[["0.95UCL"]],digits))
  xtbl$CI <- gsub('NA','Inf',xtbl$CI)
  xtbl$CI[is.na(xtbl$median)] <- NA
  xtbl$median = niceNum(xtbl$median,digits)
  xtbl$Expected =niceNum(xtbl$Expected,digits)
  xtbl <- dplyr::select(xtbl,c(grpLabel,"n.start","events","Expected","median","CI"))

  lr_p = niceNum(1-pchisq(lr_test$chisq,length(lr_test$n)-1),3)
  if (lr_p=="0.000") lr_p = "<0.001"
  xtbl <- rbind(xtbl,
                c('LR p-value',rep(NA,ncol(xtbl)-2),lr_p))
  colnames(xtbl) <- c(grpLabel,'n','Events',"Expected",'Median','95% CI')

  if (tblOnly){
    return(xtbl)
  }
  to_indent  <- which(xtbl[[grpLabel]] %in% levels(data[[group]]))
  outTable(xtbl,to_indent,caption)

}

#' Function to nicely display survival time summaries at discrete times
#' @param data datafrmae
#' @param time survival time
#' @param status indication of censored or observed data
#' @param group the variable to group observations by
#' @param survtimes numeric vector specifying the times to calculate survival at
#' @param survtimeunit the unit of time that should be indicated in the tables (days, years, etc)
#' @param grpLabel character value describing the grouping variable
#' @param survtimesLbls if supplied, a vector the same length as survtimes with descriptions (useful for displaying years with data provided in months)
#' @param digits the number of digitis to display
#' @param caption table caption
#' @param tblOnly should a dataframe or a formatted object be returned
#' @param showN should the sample sizes be displayed
#' @export
survtime_sum <- function(data,time,status,group,survtimes,survtimeunit,grpLabel='Group',survtimesLbls=NULL,digits=2,caption=NULL,tblOnly=FALSE,showN=FALSE){
  if (group==1) stop('This function requires a grouping variable.')
  if (showN) {
    n_cols <- c('n.risk','n.event','n.censor')
  } else n_cols=NULL
  if (length(survtimesLbls)!=length(survtimes)) stop('If supplied, the survtimesLbls vector must be the same length as survtime')
  if (is.null(survtimesLbls)) survtimesLbls = survtimes
  survtimeLookup <- tibble(time=survtimes,
                           lbldtime = survtimesLbls)
  timelbl <- paste0('Time (',survtimeunit,')')
  names(survtimeLookup)[2] = timelbl
  out_fmt <- knitr::opts_knit$get("rmarkdown.pandoc.to")
  if (is.null(out_fmt)) out_fmt='pdf'

  colsToExtract <- c(2:9,14,15)
  survtime_est <- NULL

  ofit<- summary(survfit(as.formula(paste0("Surv(",time,',',status,') ~1')),
                         data = data),times=survtimes)
  otbl <-NULL
  for (i in colsToExtract) {
    otbl <- cbind(otbl,ofit[[i]])
  }
  colnames(otbl) <- names(ofit)[colsToExtract]
  otbl <- as_tibble(otbl)
  survtime_est[['Overall']] <- otbl

  if (class(data[[group]])[1] %in% c('ordered','factor')){
    levelnames <- levels(data[[group]])
    levelnames <- levelnames[levelnames %in% unique(data[[group]])]
  } else  levelnames <- unique(data[[group]])
  for (g in levelnames){
    gfit <- summary(survfit(as.formula(paste0("Surv(",time,',',status,') ~1')),
                            data = data[data[[group]]==g,]),times=survtimes,extend=T)
    gtbl <-NULL
    for (i in colsToExtract) {
      gtbl <- cbind(gtbl,gfit[[i]])
    }
    colnames(gtbl) <- names(gfit)[colsToExtract]
    gtbl <- as_tibble(gtbl)
    survtime_est[[g]] <- gtbl
  }
  xtbl <- bind_rows(survtime_est,.id=grpLabel)
  xtbl$SR = niceNum(xtbl[["surv"]])
  xtbl$CI = paste0(niceNum(xtbl[["lower"]],digits),', ',niceNum(xtbl[["upper"]],digits))


  xtbl <- dplyr::left_join(xtbl,survtimeLookup)
  tab <- dplyr::select(xtbl,c(grpLabel,timelbl,all_of(n_cols),"SR","CI"))
  tab <-  dplyr::filter(tab,CI!="NA, NA")
    tab <-  dplyr::filter(tab,n.risk+n.event+n.censor>0)


  if (length(n_cols)>0){
    colnames(tab) <- c(grpLabel,timelbl,'At Risk','Events','Censored','Survival Rate','95% CI')
  } else  colnames(tab) <- c(grpLabel,timelbl,'Survival Rate','95% CI')
  tab <- kable_stk_hdr(tab,grpLabel,timelbl,tblOnly = TRUE)

  if (tblOnly){
    return(tab)
  }
  to_indent  <- which(!tab[[timelbl]] %in% c("Overall",levelnames))
  outTable(tab,to_indent,caption)

}

# To Do: update to allow a column to bold rows by
#' Stack columns in a table for clearer viewing
#' @param data dataframe
#' @param head_col character value specifying the column name with the headers
#' @param to_col character value specifying the column name to add the headers into
#' @param caption table caption
#' @param indent should the original values in the to_col be indented
#' @param hdr_prefix character value that will prefix headers
#' @param hdr_suffix character value that will suffix headers
#' @param tblOnly boolean indicating if the table should be formatted for prining or returned as a data frame
kable_stk_hdr <- function(data,head_col,to_col,caption=NULL,indent=TRUE,hdr_prefix='',hdr_suffix='',tblOnly=FALSE){
  data[[to_col]] <- as.character(data[[to_col]])
  new_row = data[1,]
  for (i in 1:ncol(new_row)) new_row[1,i] <- NA
  new_headers = unique(data[[head_col]])
  repeat{
    header_index = which(!duplicated(data[[head_col]]) & !is.na(data[[head_col]]))[1]
    new_row[[to_col]] <- data[[head_col]][header_index]

    data <- tibble::add_row(data,new_row, .before = header_index)
    data[[head_col]][data[[head_col]]==new_row[[to_col]]] <- NA
    if (sum(is.na(data[[head_col]]))==nrow(data)) break
  }
  header_rows <- which(data[[to_col]] %in% new_headers)
  to_indent <- which(!(data[[to_col]] %in% new_headers) )
  data[[to_col]][header_rows] <- paste0(hdr_prefix,data[[to_col]][header_rows],hdr_suffix)
  data <- dplyr::select(data,-all_of(head_col))

  if (tblOnly){
    return(data)
  }

  outTable(tab=data,to_indent=to_indent,caption=caption)
}

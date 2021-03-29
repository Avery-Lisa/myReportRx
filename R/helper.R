
#' Tidy rounding
#' @param x figure to be rounded, can be a vector
#' @param digits number of decimal places to display to, defaults to 2
#' @export
niceNum <- function(x,digits=2){

  rndx = sapply(x, function(x) {format(round(as.numeric(x),digits),nsmall=digits)})
  return(gsub(" ","",rndx))
}

pstprn<-function(x){paste(x[1]," (",paste(x[-1],collapse=","),")",sep="")}

# LA updated to always return a formatted string
psthr<- function (x, y = 2)
{
  x <- sapply(x, function(x) {
    ifelse(abs(x) < 0.01 | abs(x) > 1000, format(x, scientific = TRUE,
                                                 digits = y), format(round(x, y),nsmall = y))
  })
  pstprn(x)
}

formatp<- function(pvalues,sigdigits=2){
  p_out <- sapply(pvalues, function(x){
    x <- signif(as.numeric(x),sigdigits)
    x <- ifelse(is.na(x),NA_character_,ifelse(x<0.001,"<0.001",format(x)))})
  return(p_out)
}

betaWithCI <-function(betaname,CIwidth=0.95){
  paste0(betaname,"(",100*CIwidth,"%CI)")
}

lbl_count <- function(y){
  q75 <- summary(y)[5]
  return(data.frame(y=max(y),  label=paste('n =',length(y))))
}

niceStr <- function (strings)
{
  out <- sapply(strings, function(x) {
    x <- chartr('/',' ',x)
    x <- chartr(".", " ", x)
    x <- chartr("_", " ", x)
    return(x)
  })
  return(out)
}

wrp_lbl <- function(x,width = 10){
  x <- niceStr(x)
  str_wrap(x,width = width)
}

label_wrap_reportRx <- function (width = 25, multi_line = TRUE) {
  fun <- function(labels) {
    labels <- label_value(labels, multi_line = multi_line)
    lapply(labels, function(x) {
      x <- niceStr(x)
      x <- strwrap(x, width = width, simplify = FALSE)
      vapply(x, paste, character(1), collapse = "\n")
    })
  }
  structure(fun, class = "labeller")
}
# LA, 14 Dec 2020
# Fixed a bug in matchcovariate that can occur when variables are centred and spaces are added to variable names
# stored this in helper functions, in the r_code folder
matchcovariate <- function(betanames,ucall){
  out = as.vector(sapply(betanames, function(betaname) {
    splitbetaname = unlist(strsplit(betaname, ":",
                                    fixed = T))
    out = sapply(splitbetaname, function(bname) {
      bname=gsub(" ","",bname)
      indx = which(sapply(ucall, function(cov) grepl(cov,bname, fixed = TRUE)))
      if (length(indx) == 1)
        return(indx)
      indx2 <- which.max(sapply(ucall[indx], nchar))
      if (length(indx2) == 1)
        return(indx[indx2])
      indx3 <- which(sapply(ucall[indx2], function(c) {
        substr(betaname, 1, nchar(c)) == c
      }))
      if (length(indx3) == 1)
        return(ucall[indx[indx2[indx3]]])
      return(-1)
    })
    if (-1 %in% out)
      return(-1)
    result = 0
    n = length(out)
    for (i in 1:length(out)) {
      result = result + out[i] * 100^(n - 1)
      n = n - 1
    }
    return(result)
  }))
  if (-1 %in% out)
    return(-1)
  return(out)
}



betaindx <- function(x){
  i = 1
  out <- 1
  result <- NULL
  while (TRUE) {
    if (i + 1 > length(x)) {
      result <- c(result, list(out))
      return(result)
    }
    else if (alleql(x[[i + 1]], x[[i]])) {
      out <- c(out, i + 1)
    }
    else {
      result <- c(result, list(out))
      out <- i + 1
    }
    i = i + 1
  }
}

cap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

alleql <- function (x, y){
  !any((x == y) == F)
}

covnm <- function(betanames,call){
  sapply(betanames, function(betaname) {
    indx = which(sapply(call, function(cov) grepl(cov, betaname,
                                                  fixed = TRUE)))
    if (length(indx) == 1)
      return(call[indx])
    indx2 <- which.max(sapply(call[indx], nchar))
    if (length(indx2) == 1)
      return(call[indx[indx2]])
    indx3 <- which(sapply(call[indx2], function(c) {
      substr(betaname, 1, nchar(c)) == c
    }))
    if (length(indx3) == 1)
      return(call[indx[indx2[indx3]]])
  })
}


# New function to strip centering from a covariate
getvarname = function(betaname){
  sapply(betaname,function(x){
    x = gsub('I[(]','',x)
    x = gsub('[-+].*','',x)
    x = trimws(x)
    return(x)
  })
}


#'Clean strings for printing
#'
#' Returns strings with . and _ replaced by a space. This is nice when printing column names of your dataframe in a report
#' @param strings vector of strings to give a nice name
#' @keywords helper
#' @export
nicename<-function(strings){
  out<-sapply(strings,function(x){
    x<-chartr(".", " ",x)
    x<-chartr("_", " ",x)
    return(x)})
  return(out)
}

#' Formats p-values
#'
#' Returns <0.001 if pvalue is <0.001. Else rounds the pvalue to 2 significant digits
#'
#' @param x an integer
#' @export
pvalue<-function(x){
  if(is.na(x)|class(x)=="character") return(x)
  else if (x<=0.001) return("<0.001")
  else return(signif(x,2))
}

sanitize <- function(str) {
  result <- str
  result <- gsub("\\\\", "SANITIZE.BACKSLASH", result)
  result <- gsub("$", "\\$", result, fixed = TRUE)
  result <- gsub(">", "$>$", result, fixed = TRUE)
  result <- gsub("<", "$<$", result, fixed = TRUE)
  result <- gsub("|", "$|$", result, fixed = TRUE)
  result <- gsub("{", "\\{", result, fixed = TRUE)
  result <- gsub("}", "\\}", result, fixed = TRUE)
  result <- gsub("%", "\\%", result, fixed = TRUE)
  result <- gsub("&", "\\&", result, fixed = TRUE)
  result <- gsub("_", "\\_", result, fixed = TRUE)
  result <- gsub("#", "\\#", result, fixed = TRUE)
  result <- gsub("^", "\\verb|^|", result, fixed = TRUE)
  result <- gsub("~", "\\~{}", result, fixed = TRUE)
  result <- gsub("SANITIZE.BACKSLASH", "$\\backslash$",
                 result, fixed = TRUE)
  return(result)
}

#' Sanitizes strings to not break LaTeX
#'
#' Strings with special charaters will break LaTeX if returned 'asis' by knitr. This happens every time we use one of the main reportRx functions. We first sanitize our strings with this function to stop LaTeX from breaking.
#'
#'@param str a vector of strings to sanitize
#'@export
sanitizestr<-function(str){
  as.vector(sapply(str,function(char){sanitize(char)}))
}

#'Bold strings in LaTeX
#'
#'Bold strings in LaTeX.
#'
#'@param strings A vector of strings to bold.
#'@export
lbld<-function(strings){sapply(strings,function(x){
  if(is.null(x)) return(x)
  if(is.na(x)) return(x)
  return(paste("\\textbf{",x,"}",sep=""))})}

#'Add spaces to strings in LaTeX
#'
#'Add spaces to strings in LaTeX. Returns appends ~~~ before the string
#'
#'@param x string
#'@export
addspace<-function(x){
  paste("~~~",x,sep="")
}
#' Formats p-values for LaTeX
#'
#' Returns <0.001 if pvalue is <0.001. Else rounds the pvalue to 2 significant digits. Will bold the p-value if it is <= 0.05
#' @param x an integer
#' @export
lpvalue<-function(x){
  if(is.na(x)|class(x)=="character") return(x)
  else if (x<=0.001) return("\\textbf{$<$0.001}")
  else x=signif(x,2)
  if(x<=0.05) return(paste("\\textbf{",x,"}",sep=""))
  else return(x)
}



removedollar<-function(x){
  colnms<-strsplit(x,":")
  indx<-unlist(lapply(colnms,function(colnm) sapply(colnm, function(coln) regexpr("$",coln,fixed=T)[1]+1)))
  if(length(unique(indx))==1){
    if(unique(indx)!=0) x<-unlist(lapply(colnms,function(colnm) paste(substring(colnm,indx[1]),collapse=":")))
  }
  return(x)
}

modelmatrix<-function(f,data=NULL){
  k<-as.character(f)
  y<-NULL
  if(!length(k)%in%c(2,3)) stop("formula not properly formed")
  if(length(k)==3) {
    f<-as.formula(paste("~",k[2],"+",k[3],sep=""))
    y<-model.matrix(as.formula(paste("~",k[2],sep="")),data)[,-1,drop=F]}
  x<-model.matrix(f,data)[,-1,drop=F]
  colnames(x)<-removedollar(colnames(x))
  if(!is.null(y)){
    return(list(x[,1:ncol(y),drop=F],x[,(ncol(y)+1):ncol(x),drop=F]))
  }else{
    return(x)
  }}


# (ggsurv) ---------------------------------------------------------

pstprn0 <- function (x)
{
  paste0(x[1], "(", paste0(x[-1], collapse = ","),
         ")", sep = "")
}

psthr0 <- function (x, y = 2)
{
  x <- sapply(x, function(x) {
    ifelse(abs(x) < 0.01 | abs(x) > 1000, format(x, scientific = TRUE,
                                                 digits = y), round(x, y))
  })
  pstprn0(x)
}
break_function <- function(xmax){

  xmax_length <- ifelse(xmax>1,nchar(round(xmax)),round(abs(log10(xmax))))

  byx <- if(xmax>1) {round(xmax/10,digits = 2-xmax_length)
  }else round(xmax/10,digits = xmax_length+1)

  breaks <- seq(0,xmax,by=byx)
  if(max(breaks)<byx) breaks <- c(breaks,max(breaks)+byx)
  return(breaks)
}

lpvalue2 <- function (x)
{
  if (is.na(x) | class(x) == "character")
    return(x)
  else if (x <= 0.001)
    return("<0.001")
  else x = signif(x, 2)
}

.extract_ggplot_colors <- function(p, grp.levels){
  g <- ggplot_build(p)
  .cols <- unlist(unique(g$data[[1]]["colour"]))
  if(!is.null(grp.levels)){
    if(length(.cols)==1) .cols <- rep(.cols, length(grp.levels))
    names(.cols) <- grp.levels
  }
  .cols
}

.set_large_dash_as_ytext <- function(ggp){
  ggp + theme(axis.text.y = element_text(size = 50, vjust = 0.35),
              axis.ticks.y = element_blank())
}

##This function is used by the survfit package
survfit_confint <- function(p, se, logse=TRUE, conf.type, conf.int=0.95,
                            selow, ulimit=TRUE) {
  zval <- qnorm(1- (1-conf.int)/2, 0,1)
  if (missing(selow)) scale <- 1.0
  else scale <- ifelse(selow==0, 1.0, selow/se)  # avoid 0/0 at the origin
  if (!logse) se <- ifelse(se==0, 0, se/p)   # se of log(survival) = log(p)

  if (conf.type=='plain') {
    se2 <- se* p * zval  # matches equation 4.3.1 in Klein & Moeschberger
    if (ulimit) list(lower= pmax(p -se2*scale, 0), upper = pmin(p + se2, 1))
    else  list(lower= pmax(p -se2*scale, 0), upper = p + se2)
  }
  else if (conf.type=='log') {
    #avoid some "log(0)" messages
    xx <- ifelse(p==0, NA, p)
    se2 <- zval* se
    temp1 <- exp(log(xx) - se2*scale)
    temp2 <- exp(log(xx) + se2)
    if (ulimit) list(lower= temp1, upper= pmin(temp2, 1))
    else  list(lower= temp1, upper= temp2)
  }
  else if (conf.type=='log-log') {
    xx <- ifelse(p==0 | p==1, NA, p)
    se2 <- zval * se/log(xx)
    temp1 <- exp(-exp(log(-log(xx)) - se2*scale))
    temp2 <- exp(-exp(log(-log(xx)) + se2))
    list(lower = temp1 , upper = temp2)
  }
  else if (conf.type=='logit') {
    xx <- ifelse(p==0, NA, p)  # avoid log(0) messages
    se2 <- zval * se *(1 + xx/(1-xx))

    temp1 <- 1- 1/(1+exp(log(p/(1-p)) - se2*scale))
    temp2 <- 1- 1/(1+exp(log(p/(1-p)) + se2))
    list(lower = temp1, upper=temp2)
  }
  else if (conf.type=="arcsin") {
    xx <- ifelse(p==0, NA, p)
    se2 <- .5 *zval*se * sqrt(xx/(1-xx))
    list(lower= (sin(pmax(0, asin(sqrt(xx)) - se2*scale)))^2,
         upper= (sin(pmin(pi/2, asin(sqrt(xx)) + se2)))^2)
  }
  else stop("invalid conf.int type")
}

# (ordinal regression) ---------------------------------------------------------

# This is taken from the brant package
getCombiCoefs <- function (model)
{
  classes = attr(model$terms, "dataClasses")
  factors = ifelse(classes[2:length(classes)] != "numeric",
                   T, F)
  f = i = var = 1
  result = data.frame(i = 1:length(coef(model)), var = NA)
  for (factor in factors) {
    if (factor) {
      n = length(unlist(model$xlevels[f]))
      for (j in 1:(n - 1)) {
        result[i, "var"] = var
        i = i + 1
      }
      var = var + 1
      f = f + 1
    }
    else {
      result[i, "var"] = var
      var = var + 1
      i = i + 1
    }
  }
  return(result)
}
# WE NEED TO BE ATTRITBUTE TO THE ORIGINAL AUTHOR IF THIS IS POSTED ON CRAN!
# Benjamin Schlegel, Marco Steenbergen
# https://cran.r-project.org/web/packages/brant/brant.pdf
# This is code from the brant package, with  this line of code:
# result.matrix = print.testresult(model, X2, df.v, by.var)
# replaced by the contents of brant:::print.testresult) to not print to screen
# Changed: LA 16 December 2020 return the result matrix, and the number of empty cells separately in a list
# the number of empty cells
modified_brant <- function (model, by.var = F)
{
  m_model <- model$call
  if (is.matrix(eval.parent(m_model$data)))
    m_model$data <- as.data.frame(data)
  m_model[-which(names(m_model) %in% c("", "formula", "data"))] <- NULL
  m_model[[1L]] <- quote(stats::model.frame)
  m_model <- eval.parent(m_model)
  Terms <- attr(m_model, "terms")
  x <- model.matrix(Terms, m_model)
  xint <- match("(Intercept)", colnames(x), nomatch = 0L)
  x <- x[, -xint, drop = FALSE]
  y <- as.numeric(model.response(m_model))
  x.variables = names(m_model)[-1]
  temp.data = data.frame(m_model, y)
  if (grepl(":", paste0(colnames(x), collapse = "")) & by.var) {
    by.var = FALSE
    warning("by.var = TRUE currently not supported for interactions, setting by.var to FALSE")
  }
  x.factors = c()
  for (name in x.variables) {
    if (!is.numeric(m_model[, name])) {
      x.factors = c(x.factors, name)
    }
  }
  if (length(x.factors) > 0) {
    tab = table(data.frame(temp.data$y, m_model[, x.factors]))
    count0 = sum(tab==0)
    # count0 = colSums(tab==0)
  }
  else {
    count0 = 0
  }
  J = max(y, na.rm = T)
  K = length(coef(model))
  for (m in 1:(J - 1)) {
    temp.data[[paste0("z", m)]] = ifelse(y > m, 1, 0)
  }
  binary.models = list()
  beta.hat = matrix(NA, nrow = J - 1, ncol = K + 1, byrow = T)
  var.hat = list()
  for (m in 1:(J - 1)) {
    mod = glm(paste0("z", m, " ~ ", as.character(formula(model)[3])),
              data = temp.data, family = "binomial")
    binary.models[[paste0("model", m)]] = mod
    beta.hat[m, ] = coef(mod)
    var.hat[[m]] = vcov(mod)
  }
  X = cbind(1, x)
  tau = matrix(model$zeta, nrow = 1, ncol = J - 1, byrow = T)
  pi.hat = matrix(NA, nrow = length(model$model[, 1]), ncol = J -
                    1, byrow = T)
  for (m in 1:(J - 1)) {
    pi.hat[, m] = binary.models[[m]]$fitted.values
  }
  varBeta = matrix(NA, nrow = (J - 1) * K, ncol = (J - 1) *
                     K)
  for (m in 1:(J - 2)) {
    for (l in (m + 1):(J - 1)) {
      Wml = Matrix::Diagonal(x = pi.hat[, l] - pi.hat[,
                                                      m] * pi.hat[, l])
      Wm = Matrix::Diagonal(x = pi.hat[, m] - pi.hat[,
                                                     m] * pi.hat[, m])
      Wl = Matrix::Diagonal(x = pi.hat[, l] - pi.hat[,
                                                     l] * pi.hat[, l])
      Xt = t(X)
      varBeta[((m - 1) * K + 1):(m * K), ((l - 1) * K +
                                            1):(l * K)] = as.matrix((solve(Xt %*% Wm %*%
                                                                             X) %*% (Xt %*% Wml %*% X) %*% solve(Xt %*% Wl %*%
                                                                                                                   X))[-1, -1])
      varBeta[((l - 1) * K + 1):(l * K), ((m - 1) * K +
                                            1):(m * K)] = varBeta[((m - 1) * K + 1):(m *
                                                                                       K), ((l - 1) * K + 1):(l * K)]
    }
  }
  betaStar = c()
  for (m in 1:(J - 1)) {
    betaStar = c(betaStar, beta.hat[m, -1])
  }
  for (m in 1:(J - 1)) {
    varBeta[((m - 1) * K + 1):(m * K), ((m - 1) * K + 1):(m *
                                                            K)] = var.hat[[m]][-1, -1]
  }
  I = diag(1, K)
  E0 = diag(0, K)
  for (i in 1:(J - 2)) {
    for (j in 1:(J - 1)) {
      if (j == 1) {
        temp = I
      }
      else if (j == i + 1) {
        temp = cbind(temp, -I)
      }
      else {
        temp = cbind(temp, E0)
      }
    }
    if (i == 1) {
      D = temp
    }
    else {
      D = rbind(D, temp)
    }
  }
  X2 = t(D %*% betaStar) %*% solve(D %*% varBeta %*% t(D)) %*%
    (D %*% betaStar)
  df.v = (J - 2) * K
  if (by.var) {
    combinations = getCombiCoefs(model)
    for (v in unique(combinations$var)) {
      k = subset(combinations, var == v)$i
      s = c()
      df.v.temp = 0
      for (e in k) {
        s = c(s, seq(from = e, to = K * (J - 1), by = K))
        df.v.temp = df.v.temp + J - 2
      }
      s = sort(s)
      Ds = D[, s]
      if (!is.null(dim(Ds))) {
        Ds = Ds[which(!apply(Ds == 0, 1, all)), ]
      }
      if (!is.null(dim(Ds)))
        X2 = c(X2, t(Ds %*% betaStar[s]) %*% solve(Ds %*%
                                                     varBeta[s, s] %*% t(Ds)) %*% (Ds %*% betaStar[s]))
      else X2 = c(X2, t(Ds %*% betaStar[s]) %*% solve(Ds %*%
                                                        varBeta[s, s] %*% t(t(Ds))) %*% (Ds %*% betaStar[s]))
      df.v = c(df.v, df.v.temp)
    }
  }
  else {
    for (k in 1:K) {
      s = seq(from = k, to = K * (J - 1), by = K)
      Ds = D[, s]
      if (!is.null(dim(Ds))) {
        Ds = Ds[which(!apply(Ds == 0, 1, all)), ]
      }
      if (!is.null(dim(Ds)))
        X2 = c(X2, t(Ds %*% betaStar[s]) %*% solve(Ds %*%
                                                     varBeta[s, s] %*% t(Ds)) %*% (Ds %*% betaStar[s]))
      else X2 = c(X2, t(Ds %*% betaStar[s]) %*% solve(Ds %*%
                                                        varBeta[s, s] %*% t(t(Ds))) %*% (Ds %*% betaStar[s]))
      df.v = c(df.v, J - 2)
    }
  }
  # result.matrix = print.testresult(model, X2, df.v, by.var)
  p.values = pchisq(X2, df.v, lower.tail = FALSE)
  if (by.var) {
    var.names = unlist(strsplit(as.character(formula(model))[3],
                                split = " \\+ "))
  }
  else {
    var.names = names(coef(model))
  }
  longest.char = max(nchar(var.names))
  n.tabs = ceiling(longest.char/7)
  n.tabs = ifelse(n.tabs < 2, 2, n.tabs)
  result.matrix = matrix(c(X2, df.v, p.values), ncol = 3)
  rownames(result.matrix) = c("Omnibus", var.names)
  colnames(result.matrix) = c("X2", "df", "probability")
  # invisible(result.matrix)
  # Changed: LA 16 December 2020 return the result matrix, and the number of empty cells separately in a list
  list(result = result.matrix,zero_count_cells=count0)
}

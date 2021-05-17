#' get a column number (or vector) of Excel columns specified as unquoted letters
#' @param ... unquoted excel column headers ie excelCol(A,CG,AA)
#' @export
excelCol<- function(...){
  args <- as.list(match.call())[-1]
  args <-unname(unlist(lapply(args,function(x) {rlang::as_string(x)})))
  rtn<-sapply(args, function(x){
    colHead <- toupper(x)
    if (nchar(colHead)>1){
      l1 = substr(colHead,1,1)
      l2 = substr(colHead,2,2)
      rtn <- 26*which(LETTERS==l1)+which(LETTERS==l2)
    } else {
      rtn <- which(LETTERS==colHead)
    }
  })
  names(rtn) <- toupper(names(rtn))
  return(rtn)
}

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

bib_ReadGatherTidy <- function(file){
  # Note that this is an amalgamation of the bib2df read, gather and tidy functions which are included here because they are not exported from the bib2df namespace
  bib <- readLines(file)
  bib <- stringr::str_replace_all(bib, "[^[:graph:]]", " ")

  from <- which( stringr::str_extract(bib, "[:graph:]") == "@")
  to <- c(from[-1] - 1, length(bib))
  if (!length(from)) {
    return(empty)
  }
  itemslist <- mapply(function(x, y) return(bib[x:y]), x = from,
                      y = to - 1, SIMPLIFY = FALSE)
  keys <- lapply(itemslist, function(x) {
     stringr::str_extract(x[1], "(?<=\\{)[^,]+")
  })
  fields <- lapply(itemslist, function(x) {
     stringr::str_extract(x[1], "(?<=@)[^\\{]+")
  })
  fields <- lapply(fields, toupper)
  categories <- lapply(itemslist, function(x) {
     stringr::str_extract(x, "[:graph:]+")
  })
  dupl <- sum(unlist(lapply(categories, function(x) sum(duplicated(x[!is.na(x)])))))
  if (dupl > 0) {
    message("Some BibTeX entries may have been dropped.\n            The result could be malformed.\n            Review the .bib file and make sure every single entry starts\n            with a '@'.")
  }
  values <- lapply(itemslist, function(x) {
     stringr::str_extract(x, "(?<==).*")
  })
  values <- lapply(values, function(x) {
     stringr::str_extract(x, "(?![\"\\{\\s]).*")
  })
  values <- lapply(values, function(x) {
    gsub("?(^[\\{\"])", "", x)
  })
  values <- lapply(values, function(x) {
    gsub("?([\\}\"]\\,$)", "", x)
  })
  values <- lapply(values, function(x) {
    gsub("?([\\}\"]$)", "", x)
  })
  values <- lapply(values, function(x) {
    gsub("?(\\,$)", "", x)
  })
  values <- lapply(values, trimws)
  items <- mapply(cbind, categories, values, SIMPLIFY = FALSE)
  items <- lapply(items, function(x) {
    x <- cbind(toupper(x[, 1]), x[, 2])
  })
  items <- lapply(items, function(x) {
    x[complete.cases(x), ]
  })
  items <- mapply(function(x, y) {
    rbind(x, c("CATEGORY", y))
  }, x = items, y = fields, SIMPLIFY = FALSE)
  items <- lapply(items, t)
  items <- lapply(items, function(x) {
    colnames(x) <- x[1, ]
    x <- x[-1, ]
    return(x)
  })
  items <- lapply(items, function(x) {
    x <- t(x)
    x <- data.frame(x, stringsAsFactors = FALSE)
    return(x)
  })
  empty = data.frame(CATEGORY=character(), BIBTEXKEY=character(), ADDRESS=character(), ANNOTE=character(), AUTHOR=character(), BOOKTITLE=character(), CHAPTER=character(), CROSSREF=character(), EDITION=character(), EDITOR=character(), HOWPUBLISHED=character(), INSTITUTION=character(), JOURNAL=character(), KEY=character(), MONTH=character(), NOTE=character(), NUMBER=character(), ORGANIZATION=character(), PAGES=character(), PUBLISHER=character(), SCHOOL=character(), SERIES=character(), TITLE=character(), TYPE=character(), VOLUME=character(), YEAR=character())
  dat <- dplyr::bind_rows(c(list(empty), items))
  dat <- as.data.frame(dat)
  dat$BIBTEXKEY <- unlist(keys)
  bib <- dat
  if (dim(bib)[1] == 0) {
    return(bib)
  }
  AUTHOR <- EDITOR <- YEAR <- CATEGORY <- NULL
  if ("AUTHOR" %in% colnames(bib)) {
    bib$AUTHOR = strsplit(bib$AUTHOR, " and ",fixed = TRUE)
  }
  if ("EDITOR" %in% colnames(bib)) {
    bib$EDITOR = strsplit(bib$EDITOR, " and ",fixed = TRUE)
  }
  if ("YEAR" %in% colnames(bib)) {
    if (sum(is.na(as.numeric(bib$YEAR))) == 0) {
      bib$YEAR = as.numeric(bib$YEAR)
    } else { warning('Check YEAR in bibfile, may be some missing or character strings.')}
  }
  return(bib)
}

hmisc_summaryF <- function (formula, data = NULL, subset = NULL, na.action = NULL,
          fun = NULL, method = c("response", "reverse",
                                 "cross"), overall = method == "response" |
            method == "cross", continuous = 10, na.rm = TRUE,
          na.include = method != "reverse", g = 4, quant = c(0.025,
                                                             0.05, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 0.95,
                                                             0.975), nmin = if (method == "reverse") 100 else 0,
          test = FALSE, conTest = conTestkw, catTest = catTestchisq,
          ordTest = ordTestpo, ...)
{
  # NOTE: This is the Hmisc:::summary.formula  function copied directly here for testing of the po assumption in ordinal regression
  call <- match.call()
  missmethod <- missing(method)
  method <- match.arg(method)
  if (grepl(".*\\+.*~", paste(deparse(formula), collapse = "")))
    return(summaryM(formula, data = data, subset = subset,
                    na.action = na.action, overall = overall, continuous = continuous,
                    na.include = na.include, quant = quant, nmin = nmin,
                    test = test, conTest = conTest, catTest = catTest,
                    ordTest = ordTest))
  X <- match.call(expand.dots = FALSE)
  X$fun <- X$method <- X$na.rm <- X$na.include <- X$g <- X$overall <- X$continuous <- X$quant <- X$nmin <- X$test <- X$conTest <- X$catTest <- X$... <- NULL
  if (missing(na.action))
    X$na.action <- na.retain
  Terms <- if (missing(data))
    terms(formula, "stratify")
  else terms(formula, "stratify", data = data)
  X$formula <- Terms
  X[[1]] <- as.name("model.frame")
  X <- eval(X, sys.parent())
  Terms <- attr(X, "terms")
  resp <- attr(Terms, "response")
  if (resp == 0 && missmethod)
    method <- "reverse"
  if (test && method != "reverse")
    stop("test=TRUE only allowed for method=\"reverse\"")
  if (method != "reverse" && resp != 1)
    stop("must have a variable on the left hand side of the formula")
  nact <- attr(X, "na.action")
  nvar <- ncol(X) - 1
  strat <- attr(Terms, "specials")$stratify
  getlab <- function(x, default) {
    lab <- attr(x, "label")
    if (!length(lab) || lab == "")
      default
    else lab
  }
  if (length(strat)) {
    if (method != "response")
      stop("stratify only allowed for method=\"response\"")
    temp <- untangle.specials(Terms, "stratify")
    strat.name <- var.inner(Terms)[temp$terms]
    strat <- if (length(temp$vars) == 1)
      as.factor(X[[temp$vars]])
    else stratify(X[, temp$vars])
    strat.label <- getlab(X[, temp$vars[1]], strat.name)
    X[[temp$vars]] <- NULL
  }
  else {
    strat <- factor(rep("", nrow(X)))
    strat.name <- strat.label <- ""
  }
  nstrat <- length(levels(strat))
  if (resp > 0) {
    Y <- X[[resp]]
    yname <- as.character(attr(Terms, "variables"))[2]
    ylabel <- getlab(Y, yname)
    if (!is.matrix(Y))
      Y <- matrix(Y, dimnames = list(names(Y), yname))
  }
  else {
    yname <- ylabel <- NULL
  }
  if (method != "reverse") {
    if (!length(fun)) {
      fun <- function(y) apply(y, 2, mean)
      uy <- unique(Y[!is.na(Y)])
      r <- range(uy, na.rm = TRUE)
      funlab <- if (length(uy) == 2 && r[1] == 0 & r[2] ==
                    1)
        "Fraction"
      else "Mean"
      funlab <- paste(funlab, "of", yname)
    }
    else if (is.character(fun) && fun == "%") {
      fun <- function(y) {
        stats <- 100 * apply(y, 2, mean)
        names(stats) <- paste(dimnames(y)[[2]], "%")
        stats
      }
      funlab <- paste("% of", yname)
    }
    s <- if (inherits(Y, "Surv"))
      as.vector((1 * is.na(unclass(Y))) %*% rep(1, ncol(Y)) >
                  0)
    else ((if (is.character(Y))
      Y == "" | Y == "NA"
      else is.na(Y)) %*% rep(1, ncol(Y))) > 0
    stats <- if (length(dim(Y)))
      fun(Y[!s, , drop = FALSE])
    else fun(Y[!s])
    nstats <- length(stats)
    name.stats <- if (length(dn <- dimnames(stats)) == 2)
      as.vector(outer(dn[[1]], dn[[2]], FUN = function(a,
                                                       b) paste(b, a)))
    else names(stats)
    if (length(fun)) {
      if (length(de <- deparse(fun)) == 2) {
        de <- as.list(fun)
        de <- as.character(de[[length(de)]])
        funlab <- if (de[1] == "apply")
          de[length(de)]
        else de[1]
      }
      else funlab <- as.character(substitute(fun))
    }
    if (funlab[1] == "")
      funlab <- yname
    if (length(name.stats) == 0) {
      name.stats <- if (nstats == 1)
        yname
      else paste0(yname, 1:nstats)
    }
  }
  if (method == "response") {
    X[[resp]] <- NULL
    s <- if (!na.rm)
      FALSE
    else if (inherits(Y, "Surv"))
      as.vector((1 * is.na(unclass(Y))) %*% rep(1, ncol(Y)) >
                  0)
    else ((if (is.character(Y))
      Y == "" | Y == "NA"
      else is.na(Y)) %*% rep(1, ncol(Y))) > 0
    nmissy <- sum(s)
    if (nmissy) {
      X <- X[!s, , drop = FALSE]
      Y <- Y[!s, , drop = FALSE]
      strat <- strat[!s]
    }
    nc <- nstrat * (1 + nstats)
    colname <- rep(c("N", name.stats), nstrat)
    rowname <- vname <- vlabel <- vunits <- res <- NULL
    dm <- dim(X)
    nx <- dm[2]
    n <- dm[1]
    nlevels <- integer(nx)
    labels <- character(nx)
    units <- labels
    i <- 0
    nams <- c(names(X), if (overall) "Overall")
    for (v in nams) {
      i <- i + 1
      x <- if (v == "Overall")
        factor(rep("", n))
      else X[[v]]
      if (inherits(x, "mChoice"))
        x <- as.numeric(x)
      labels[i] <- getlab(x, nams[i])
      units[i] <- if (length(l <- attr(x, "units")))
        l
      else ""
      if (!(ismc <- is.matrix(x))) {
        s <- is.na(x)
        if (!is.factor(x)) {
          xu <- unique(x[!s])
          lu <- length(xu)
          x <- if (lu < continuous) {
            r <- range(xu)
            if (lu == 2 && r[1] == 0 && r[2] == 1)
              factor(x, labels = c("No", "Yes"))
            else factor(x)
          }
          else cut2(x, g = g, ...)
        }
        if (na.include && any(s)) {
          x <- na.include(x)
          levels(x)[is.na(levels(x))] <- "NA"
        }
        xlev <- levels(x)
        if (nmin > 0) {
          nn <- table(x)
          xlev <- names(nn)[nn >= nmin]
        }
      }
      else {
        xlev <- dimnames(x)[[2]]
        if (!length(xlev))
          stop("matrix variables must have column dimnames")
        if (!is.logical(x)) {
          if (is.numeric(x))
            x <- x == 1
          else {
            x <- structure(casefold(x), dim = dim(x))
            x <- x == "present" | x == "yes"
          }
        }
        if (nmin > 0) {
          nn <- apply(x, 2, sum, na.rm = TRUE)
          xlev <- xlev[nn >= nmin]
        }
      }
      nlevels[i] <- length(xlev)
      for (lx in xlev) {
        r <- NULL
        for (js in levels(strat)) {
          j <- if (ismc)
            strat == js & x[, lx]
          else strat == js & x == lx
          if (!na.include)
            j[is.na(j)] <- FALSE
          nj <- sum(j)
          f <- if (nj) {
            statz <- unlist(fun(Y[j, , drop = FALSE]))
            if (length(statz) != nstats)
              stop(paste("fun for stratum", lx,
                         js, "did not return", nstats, "statistics"))
            matrix(statz, ncol = nstats, byrow = TRUE)
          }
          else rep(NA, nstats)
          r <- c(r, nj, f)
        }
        res <- rbind(res, r)
      }
      rowname <- c(rowname, xlev)
      bl <- rep("", length(xlev) - 1)
      vname <- c(vname, v, bl)
      vlabel <- c(vlabel, labels[i], bl)
      vunits <- c(vunits, units[i], bl)
    }
    rowname[rowname == "NA"] <- "Missing"
    dimnames(res) <- list(rowname, colname)
    at <- list(formula = formula, call = call, n = n, nmiss = nmissy,
               yname = yname, ylabel = ylabel, ycolname = if (length(d <- dimnames(Y)[[2]])) d else yname,
               funlab = funlab, vname = vname, vlabel = vlabel,
               nlevels = nlevels, labels = labels, units = units,
               vunits = vunits, strat.name = strat.name, strat.label = strat.label,
               strat.levels = levels(strat))
    attributes(res) <- c(attributes(res), at)
    attr(res, "class") <- "summary.formula.response"
    return(res)
  }
  if (method == "reverse") {
    quants <- unique(c(quant, 0.025, 0.05, 0.125, 0.25, 0.375,
                       0.5, 0.625, 0.75, 0.875, 0.95, 0.975))
    if (resp) {
      group <- as.factor(X[[resp]])[, drop = TRUE]
      group.freq <- table(group)
      group.freq <- group.freq[group.freq > 0]
      if (overall)
        group.freq <- c(group.freq, Combined = sum(group.freq))
    }
    else {
      group <- rep(0, nrow(X))
      group.freq <- NULL
    }
    nv <- ncol(X) - resp
    n <- integer(nv)
    type <- n
    nams <- names(X)
    comp <- dat <- vector("list", nv)
    names(comp) <- names(dat) <- if (resp)
      nams[-1]
    else nams
    labels <- Units <- vector("character", nv)
    if (test) {
      testresults <- vector("list", nv)
      names(testresults) <- names(comp)
    }
    for (i in 1:nv) {
      w <- X[[resp + i]]
      if (length(attr(w, "label")))
        labels[i] <- attr(w, "label")
      if (length(attr(w, "units")))
        Units[i] <- attr(w, "units")
      if (!inherits(w, "mChoice")) {
        if (!is.factor(w) && !is.logical(w) && length(unique(w[!is.na(w)])) <
            continuous)
          w <- as.factor(w)
        s <- !is.na(w)
        if (na.include && !all(s) && length(levels(w))) {
          w <- na.include(w)
          levels(w)[is.na(levels(w))] <- "NA"
          s <- rep(TRUE, length(s))
        }
        n[i] <- sum(s)
        w <- w[s]
        g <- group[s, drop = TRUE]
        if (is.factor(w) || is.logical(w)) {
          tab <- table(w, g)
          if (test) {
            if (is.ordered(w))
              testresults[[i]] <- ordTest(g, w)
            else testresults[[i]] <- catTest(tab)
          }
          if (nrow(tab) == 1) {
            b <- casefold(dimnames(tab)[[1]], upper = TRUE)
            pres <- c("1", "Y", "YES",
                      "PRESENT")
            abse <- c("0", "N", "NO",
                      "ABSENT")
            jj <- match(b, pres, nomatch = 0)
            if (jj > 0)
              bc <- abse[jj]
            else {
              jj <- match(b, abse, nomatch = 0)
              if (jj > 0)
                bc <- pres[jj]
            }
            if (jj) {
              tab <- rbind(tab, rep(0, ncol(tab)))
              dimnames(tab)[[1]][2] <- bc
            }
          }
          if (overall)
            tab <- cbind(tab, Combined = apply(tab, 1,
                                               sum))
          comp[[i]] <- tab
          type[i] <- 1
        }
        else {
          sfn <- function(x, quant) {
            o <- options("digits")
            options(digits = 15)
            on.exit(options(o))
            c(quantile(x, quant), Mean = mean(x), SD = sqrt(var(x)))
          }
          qu <- tapply(w, g, sfn, simplify = TRUE, quants)
          if (test)
            testresults[[i]] <- conTest(g, w)
          if (overall)
            qu$Combined <- sfn(w, quants)
          comp[[i]] <- matrix(unlist(qu), ncol = length(quants) +
                                2, byrow = TRUE, dimnames = list(names(qu),
                                                                 c(format(quants), "Mean", "SD")))
          if (any(group.freq <= nmin))
            dat[[i]] <- lapply(split(w, g), nmin = nmin,
                               function(x, nmin) if (length(x) <= nmin)
                                 x
                               else NULL)
          type[i] <- 2
        }
      }
      else {
        w <- as.numeric(w) == 1
        n[i] <- nrow(w)
        g <- as.factor(group)
        ncat <- ncol(w)
        tab <- matrix(NA, nrow = ncat, ncol = length(levels(g)),
                      dimnames = list(dimnames(w)[[2]], levels(g)))
        if (test) {
          pval <- numeric(ncat)
          names(pval) <- dimnames(w)[[2]]
          d.f. <- stat <- pval
        }
        for (j in 1:ncat) {
          tab[j, ] <- tapply(w[, j], g, sum, simplify = TRUE,
                             na.rm = TRUE)
          if (test) {
            tabj <- rbind(table(g) - tab[j, ], tab[j,
            ])
            st <- catTest(tabj)
            pval[j] <- st$P
            stat[j] <- st$stat
            d.f.[j] <- st$df
          }
        }
        if (test)
          testresults[[i]] <- list(P = pval, stat = stat,
                                   df = d.f., testname = st$testname, statname = st$statname,
                                   namefun = st$namefun, latexstat = st$latexstat,
                                   plotmathstat = st$plotmathstat)
        if (overall)
          tab <- cbind(tab, Combined = apply(tab, 1,
                                             sum))
        comp[[i]] <- tab
        type[i] <- 3
      }
    }
    labels <- ifelse(nchar(labels), labels, names(comp))
    return(structure(list(stats = comp, type = type, group.name = if (resp) nams[1] else NULL,
                          group.label = ylabel, group.freq = group.freq, labels = labels,
                          units = Units, quant = quant, data = dat, N = sum(!is.na(group)),
                          n = n, testresults = if (test) testresults else NULL,
                          call = call, formula = formula), class = "summary.formula.reverse"))
  }
  if (method == "cross") {
    X[[resp]] <- NULL
    Levels <- vector("list", nvar)
    nams <- names(X)
    names(Levels) <- names(X)
    labels <- character(nvar)
    for (i in 1:nvar) {
      xi <- X[[i]]
      if (inherits(xi, "mChoice"))
        xi <- factor(format(xi))
      else if (is.matrix(xi) && ncol(xi) > 1)
        stop("matrix variables not allowed for method=\"cross\"")
      labels[i] <- getlab(xi, nams[i])
      if (is.factor(xi))
        xi <- xi[, drop = TRUE]
      if (!is.factor(xi) && length(unique(xi[!is.na(xi)])) >=
          continuous)
        xi <- cut2(xi, g = g, ...)
      X[[i]] <- na.include(as.factor(xi))
      levels(X[[i]])[is.na(levels(X[[i]]))] <- "NA"
      Levels[[i]] <- c(levels(X[[i]]), if (overall) "ALL")
    }
    df <- expand.grid(Levels)
    nl <- nrow(df)
    N <- Missing <- integer(nl)
    na <- is.na(Y %*% rep(1, ncol(Y)))
    S <- matrix(NA, nrow = nl, ncol = nstats, dimnames = list(NULL,
                                                              name.stats))
    chk <- function(z, nstats) {
      if (length(z) != nstats)
        stop(paste("fun did not return", nstats,
                   "statistics for a stratum"))
      z
    }
    if (nvar == 1) {
      df1 <- as.character(df[[1]])
      x1 <- X[[1]]
      for (i in 1:nl) {
        s <- df1[i] == "ALL" | x1 == df1[i]
        w <- if (na.rm)
          s & !na
        else s
        N[i] <- sum(w)
        Missing[i] <- sum(na[s])
        S[i, ] <- if (any(w))
          chk(fun(Y[w, , drop = FALSE]), nstats)
        else rep(NA, nstats)
      }
    }
    else if (nvar == 2) {
      df1 <- as.character(df[[1]])
      df2 <- as.character(df[[2]])
      x1 <- X[[1]]
      x2 <- X[[2]]
      for (i in 1:nl) {
        s <- (df1[i] == "ALL" | x1 == df1[i]) &
          (df2[i] == "ALL" | x2 == df2[i])
        w <- if (na.rm)
          s & !na
        else s
        N[i] <- sum(w)
        Missing[i] <- sum(na[s])
        S[i, ] <- if (any(w))
          chk(fun(Y[w, , drop = FALSE]), nstats)
        else rep(NA, nstats)
      }
    }
    else if (nvar == 3) {
      df1 <- as.character(df[[1]])
      df2 <- as.character(df[[2]])
      df3 <- as.character(df[[3]])
      x1 <- X[[1]]
      x2 <- X[[2]]
      x3 <- X[[3]]
      for (i in 1:nl) {
        s <- (df1[i] == "ALL" | x1 == df1[i]) &
          (df2[i] == "ALL" | x2 == df2[i]) & (df3[i] ==
                                                "ALL" | x3 == df3[i])
        w <- if (na.rm)
          s & !na
        else s
        N[i] <- sum(w)
        Missing[i] <- sum(na[s])
        S[i, ] <- if (any(w))
          chk(fun(Y[w, , drop = FALSE]), nstats)
        else rep(NA, nstats)
      }
    }
    else stop("no more than 3 independent variables allowed")
    lab <- names(df)
    lab2 <- if (length(lab) > 1)
      paste(lab, collapse = ", ")
    else lab
    heading <- paste(funlab, "by", lab2)
    S <- S[, , drop = TRUE]
    attr(S, "label") <- yname
    df$S <- S
    df$N <- N
    df$Missing <- Missing
    a <- list(heading = heading, byvarnames = lab2, Levels = Levels,
              labels = labels, na.action = nact, formula = formula,
              call = call, yname = yname, ylabel = ylabel, class = c("summary.formula.cross",
                                                                     "data.frame"))
    attributes(df) <- c(attributes(df), a)
    df
  }
}

# Training-only full indicator coding; unseen levels have all-zero indicators.
fold_prepare <- function(train, test) {
  med <- median(train$age_months,na.rm=TRUE)
  if(!is.finite(med)) med <- 0
  train$age_numeric <- ifelse(is.na(train$age_months),med,train$age_months)
  test$age_numeric <- ifelse(is.na(test$age_months),med,test$age_months)
  list(train=train,test=test,median=med)
}
fold_matrix <- function(train,test,formula) {
  vars <- all.vars(formula)
  cats <- intersect(vars,c('coarse_cell_type','intervention','species','dataset_id','time_bin','age_bin'))
  z <- rbindlist(list(train,test),fill=TRUE)
  blocks <- list()
  for(v in vars) {
    if(v %in% cats) {
      lev <- sort(unique(as.character(train[[v]])))
      stopifnot(length(lev)>0L,!anyNA(lev))
      b <- outer(as.character(z[[v]]),lev,FUN=function(a,b) a==b)*1
      colnames(b) <- paste0(v,"[",lev,"]")
    } else b <- matrix(z[[v]],ncol=1,dimnames=list(NULL,v))
    blocks[[v]] <- b
  }
  tm <- terms(formula)
  parts <- if(attr(tm,"intercept")) list(matrix(1,nrow(z),1,dimnames=list(NULL,"(Intercept)"))) else list()
  for(term in attr(tm,"term.labels")) {
    vs <- strsplit(term,":",fixed=TRUE)[[1]]
    b <- blocks[[vs[1]]]
    if(length(vs)>1L) for(v in vs[-1]) {
      a <- b; c <- blocks[[v]]
      b <- do.call(cbind,lapply(seq_len(ncol(a)),function(j) sweep(c,1,a[,j],"*")))
      colnames(b) <- unlist(lapply(colnames(a),function(n) paste(n,colnames(c),sep=":")))
    }
    parts[[length(parts)+1L]] <- b
  }
  x <- do.call(cbind,parts)
  stopifnot(nrow(x)==nrow(train)+nrow(test),all(is.finite(x)))
  x
}

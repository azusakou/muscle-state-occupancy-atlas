
library(data.table)
source("code/train_only_encoding.R")
tr <- data.table(age_months=c(6,NA,24,6),baseline_logit=c(-2,-1,0,1),
 coarse_cell_type=c("a","b","a","b"),intervention=c("load","rest","load","rest"),species=rep("mouse",4))
te <- data.table(age_months=c(NA,60),baseline_logit=c(1,2),coarse_cell_type=c("a","b"),
 intervention=c("novel","rest"),species=c("human","mouse"))
f <- fold_prepare(copy(tr),copy(te))
form <- ~baseline_logit*coarse_cell_type+intervention*coarse_cell_type+age_numeric+species
x <- fold_matrix(f$train,f$test,form)
te2 <- copy(te);te2[,age_months:=999];te2[,intervention:="another_novel"]
g <- fold_prepare(copy(tr),te2);y <- fold_matrix(g$train,g$test,form)
stopifnot(identical(x[1:4,],y[1:4,]),identical(colnames(x),colnames(y)),
 f$median==6,f$test$age_numeric[1]==6,
 all(x[5,grepl("^intervention",colnames(x))]==0),
 all(x[5,grepl("^species",colnames(x))]==0),
 !identical(x[5,grepl("^intervention",colnames(x))],x[1,grepl("^intervention",colnames(x))]))
cat("PASS: training invariance, unseen all-zero coding, fixed column vocabulary, training age imputation\n")

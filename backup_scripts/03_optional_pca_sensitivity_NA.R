# Optional sensitivity check: exclude the 'NA' class from PCA.
# This is NOT the published model. Do not substitute these results
# in the manuscript without updating all PCA tables and figures.
# Requires packages installed by scripts/01_reproduce_article.R.
stopifnot(requireNamespace("readxl",quietly=TRUE),
          requireNamespace("janitor",quietly=TRUE))
root <- getwd()
if (!file.exists(file.path(root,"data","raw","Bacia_percentual.xlsx")))
  stop("Open the repository root first")
inputs <- c(Basin="Bacia_percentual.xlsx",DOF="fod_percentual.xlsx",
            SSF="fes_percentual.xlsx",DSF="fed_percentual.xlsx")
run_one <- function(filename,omit_na=FALSE) {
  df <- janitor::clean_names(readxl::read_xlsx(file.path(root,"data","raw",filename)))
  df <- df[order(df$ano),,drop=FALSE]
  x <- df[,setdiff(names(df), c("ano",if(omit_na) "na")),drop=FALSE]
  x <- x[,vapply(x,is.numeric,logical(1)),drop=FALSE]
  x <- x[,vapply(x,function(v) sum(is.finite(v))>=2 && sd(v,na.rm=TRUE)>0,
                  logical(1)),drop=FALSE]
  if(anyNA(x)) stop("Missing values detected")
  p <- prcomp(x,center=TRUE,scale.=TRUE)
  v <- 100*p$sdev^2/sum(p$sdev^2)
  c(PC1=v[1],PC2=v[2],Cum=v[1]+v[2],n_class=ncol(x))
}
result <- do.call(rbind,lapply(names(inputs),function(k) {
  a <- run_one(inputs[[k]],FALSE);b <- run_one(inputs[[k]],TRUE)
  data.frame(Unit=k,PC1_published=a[1],PC2_published=a[2],
             Cumulative_published=a[3],Cumulative_no_NA=b[3],
             Difference_pp=b[3]-a[3],
             Classes_published=a[4],Classes_no_NA=b[4])
}))
print(result,row.names=FALSE,digits=5)
out<-file.path(root,"outputs");dir.create(out,showWarnings=FALSE)
write.csv(result,file.path(out,"PCA_sensitivity_NA.csv"),row.names=FALSE)

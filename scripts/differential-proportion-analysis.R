#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(edgeR))
parser <- ArgumentParser(description='Differential Analysis for Relative Abundance of Two Counts')
parser$add_argument('--matrix-1', type='character', required=TRUE,
    help='first count matrix (numerator). Rows are genes. Columns are samples.')
parser$add_argument('--matrix-2', type='character', required=TRUE,
    help='second count matrix. sum of two matrix is denominator. Rows are genes. Columns are samples.')
parser$add_argument('--metadata', type='character', required=TRUE,
    help = 'input metadata, first line should be a header, the first column should be sample id')
parser$add_argument('--label-field', type='character', required=TRUE,
    help = 'label field for differential analysis, should present metadata table')
parser$add_argument('--covariate-fields', type='character', required=FALSE,
    help = 'Covariates to adjust in differential analysis.')
parser$add_argument('--case-label', type='character', required=TRUE,
    help = 'label of experiment group / case group.')
parser$add_argument('--control-label', type='character', required=TRUE,
    help = 'label of control group.')
parser$add_argument('--test', type='character', default="binomial",
    choices=c('binomial', 'quasibinomial','betabinomial'),
    help='differential expression method to use')
parser$add_argument('--correction', type='character', default="BH",
    choices=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY"),
    help='method for mutiple test correction')
parser$add_argument('--output', type='character', required=TRUE,
    help='output differential expression table')
args <- parser$parse_args()

library(aod)

message(paste0("[",Sys.time(),"] ", "Load first counts matrix ..."))
count.matrix.one <- read.table(args$matrix_1, header = TRUE, row.names=1, check.names=FALSE, sep='\t')
message(paste0("[",Sys.time(),"] ", "Load second counts matrix ..."))
count.matrix.two <- read.table(args$matrix_2, header = TRUE, row.names=1, check.names=FALSE, sep='\t')

count.matrix.two <- count.matrix.two[rownames(count.matrix.one),colnames(count.matrix.one)]

message(paste0("[",Sys.time(),"] ", "Load metadata ..."))
metadata <- read.table(args$metadata, header = TRUE, row.names=1, check.names=TRUE, sep='\t', stringsAsFactors = T)

sample.ids <- intersect(rownames(metadata),colnames(count.matrix.one))
message(paste0("[",Sys.time(),"] ", length(sample.ids), " samples appear in both matrix ."))

count.matrix.one <- count.matrix.one[,sample.ids]
count.matrix.two <- count.matrix.two[,sample.ids]
metadata <- metadata[sample.ids,,drop=FALSE]
metadata[[args$label_field]] <- relevel(metadata[[args$label_field]],ref=args$control_label) 

message(paste0("[",Sys.time(),"] ", "Filter features by total coverage ..."))
count.matrix.total <- count.matrix.one + count.matrix.two
y <- edgeR::DGEList(counts=count.matrix.total)
keep <- edgeR::filterByExpr(y, group = metadata[[args$label_field]])
message(paste0("[",Sys.time(),"] ", "Among ",nrow(count.matrix.total)," features, ",sum(keep)," passed the filtered"))

count.matrix.one <- count.matrix.one[keep,]
count.matrix.two <- count.matrix.two[keep,]
rm(count.matrix.total)
rm(y)

message(paste0("[",Sys.time(),"] ", "Test for difference ..."))

proportion.test <- function(counts.one,counts.two,metadata,test.method,regressors,label_field,case_label){
  metadata[["counts.one"]] <- as.numeric(counts.one[1,])
  metadata[["counts.two"]] <- as.numeric(counts.two[1,])
  formula <- as.formula(paste("cbind(counts.one,counts.two)~", paste(regressors, collapse="+")))
  coefficient <- paste0(label_field,case_label)
  if(test.method %in% c('binomial', 'quasibinomial')){
   res <- glm(formula, family=test.method, data=metadata)  
   res.sum <- summary(res)
   if(coefficient %in% rownames(res.sum[["coefficients"]])){
     return(res.sum[["coefficients"]][coefficient,])
   }else{
    return(NULL)
   }
  }else if(test.method=="betabinomial"){
   #random.formula <- as.formula(paste("~",paste(regressors, collapse="+")))
   random.formula <- as.formula("~ 1")
   state <-  try(capture<-capture.output(res<-betabin(formula,random=random.formula, data=metadata)), silent=TRUE)
   if("try-error" %in% class(state)){return(NULL);}
   res.sum <- summary(res)
   if(coefficient %in% rownames(res.sum@Coef)){
     res.coef <- t(res.sum@Coef[coefficient,])
     colnames(res.coef) <- NULL
     return(res.coef)
   }else{
    return(NULL)
   }
  }
}

regressors <- c(args$label_field)
if("covariate.fields" %in% names(args)){
    covariate.fields <- unlist(strsplit(args$covariate_fields,","))
    regressors <- c(regressors,covariate.fields)
}

gene.ids <- rownames(count.matrix.one)
results <- lapply(gene.ids,
                  function(gene.id){
                           proportion.test(count.matrix.one[gene.id,],
                                           count.matrix.two[gene.id,],
                                           metadata,
                                           args$test,
                                           regressors,
                                           args$label_field,args$case_label)})
names(results) <- gene.ids
results <- results[unlist(lapply(results,function(x){!is.null(x)}))]
case.vs.control.table <-  t(as.data.frame(results,check.names=F))
# p.adjust(p, method = p.adjust.methods, n = length(p))
case.vs.control.table[["p.adjust"]] <-  p.adjust(case.vs.control.table[["Pr(> |z|)"]],method=args$correction)
message(paste0("[",Sys.time(),"] ", "Save results to ", args$output, " ..."))
write.table(case.vs.control.table, args$output, sep='\t', quote=FALSE, row.names=TRUE)
message(paste0("[",Sys.time(),"] ", "All done ."))

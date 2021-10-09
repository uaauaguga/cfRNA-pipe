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
parser$add_argument('--test', type='character', default="betabinomial",
    choices=c('binomial', 'quasibinomial','betabinomial'),
    help='differential expression method to use')
parser$add_argument('--correction', type='character', default="BH",
    choices=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY"),
    help='method for mutiple test correction')
parser$add_argument('--pseudocount', type='integer', default=1,
    help='pseudocount to add')
parser$add_argument('--cores', type='integer', default=1,
    help='number of cores to used')
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

count.matrix.one <- count.matrix.one[,sample.ids] + args$pseudocount
count.matrix.two <- count.matrix.two[,sample.ids] + args$pseudocount
metadata <- metadata[sample.ids,,drop=FALSE]
metadata[[args$label_field]] <- relevel(metadata[[args$label_field]],ref=args$control_label) 

message(paste0("[",Sys.time(),"] ", "Filter features by total coverage ..."))
count.matrix.total <- count.matrix.one + count.matrix.two
y <- edgeR::DGEList(counts=count.matrix.total)
keep <- edgeR::filterByExpr(y, group = metadata[[args$label_field]], min.count = 20)
message(paste0("[",Sys.time(),"] ", "Among ",nrow(count.matrix.total)," features, ",sum(keep)," passed the filtered"))

count.matrix.one <- count.matrix.one[keep,]
count.matrix.two <- count.matrix.two[keep,]
rm(count.matrix.total)
rm(y)


proportion.test <- function(x){
  counts.one <- x[[1]]
  counts.two <- x[[2]]
  if(x[[3]]%%1000==0){message(paste0("[",Sys.time(),"] ", x[[3]], " genes processed."))}
  counts.data <- metadata
  counts.data[["counts.one"]] <- counts.one
  counts.data[["counts.two"]] <- counts.two
  formula <- as.formula(paste("cbind(counts.one,counts.two)~", paste(regressors, collapse="+")))
   coefficient <- paste0(args$label_field,args$case_label)
  if(args$test %in% c('binomial', 'quasibinomial')){
   res <- glm(formula, family=args$test, data=counts.data)  
   res.sum <- summary(res)
   if(coefficient %in% rownames(res.sum[["coefficients"]])){
     return(res.sum[["coefficients"]][coefficient,])
   }else{
    return(NULL)
   }
  }else if(args$test=="betabinomial"){
   random.formula <- as.formula(paste("~",paste(regressors, collapse="+")))
   #random.formula <- as.formula("~1")
   state <-  try(capture<-capture.output(res<-betabin(formula,random=random.formula, data=counts.data)), silent=TRUE)
   if("try-error" %in% class(state)){print(capture);return(NULL);}
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
if("covariate_fields" %in% names(args)){
    covariate.fields <- unlist(strsplit(args$covariate_fields,","))
    regressors <- c(regressors,covariate.fields)
}


gene.ids <- rownames(count.matrix.one)

message(paste0("[",Sys.time(),"] ", "Prepare data ..."))
counts.list <- list()

count.matrix.one <- as.matrix(count.matrix.one)
count.matrix.two <- as.matrix(count.matrix.two)
gene.ids <- rownames(count.matrix.one)
for(i in 1:nrow(count.matrix.one)){
  if(i%%10000==0){message(paste0("[",Sys.time(),"] ", i, " genes loaded."))}
  counts.list[[i]] <- list(as.numeric(count.matrix.one[i,]),as.numeric(count.matrix.two[i,]),i)
}
names(counts.list) <- gene.ids


library(parallel)
message(paste0("[",Sys.time(),"] ", "Use ",args$cores, " cores for differential testing ..."))
results <- mclapply(counts.list,proportion.test,mc.cores=args$cores)
# results <- lapply(counts.list,proportion.test)
# for(i in 1:length(results)){if("try-error" %in% class(results[[i]])){print(gene.ids[i]);print(results[[i]])}}
names(results) <- gene.ids
results <- results[unlist(lapply(results,function(x){(!is.null(x))&&(!"try-error" %in% class(x))}))]
case.vs.control.table <-  as.data.frame(t(as.data.frame(results,check.names=F)))
colnames(case.vs.control.table) <- c("log.odds.ratio","std.err","z.score","p.value")
p.adjusted <- p.adjust(case.vs.control.table[,"p.value"],method=args$correction)
names(p.adjusted) <- NULL
case.vs.control.table$p.adjusted <- p.adjusted
message(paste0("[",Sys.time(),"] ", "Save results to ", args$output, " ..."))
write.table(case.vs.control.table, args$output, sep='\t', quote=FALSE, row.names=TRUE)
message(paste0("[",Sys.time(),"] ", "All done ."))

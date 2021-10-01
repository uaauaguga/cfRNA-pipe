#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description='Differential Expression Analysis')
parser$add_argument('--matrix', type='character', required=TRUE,
    help='input count matrix. Rows are genes. Columns are samples.')
parser$add_argument('--metadata', type='character', required=TRUE,
    help = 'input metadata, first line should be a header, the first column should be sample id')
parser$add_argument('--label-field', type='character', required=TRUE,
    help = 'label field for differential analysis, should present in metadata table')
parser$add_argument('--covariate-fields', type='character', required=FALSE,
    help = 'Covariates to adjust in differential analysis.')
parser$add_argument('--case-label', type='character', required=TRUE,
    help = 'label of experiment group / case group.')
parser$add_argument('--control-label', type='character', required=TRUE,
    help = 'label of control group.')
parser$add_argument('--test', type='character', default="edger-glmqlf",
    choices=c('edger-glmqlf', 'edger-glmlrt', 'limma-voom', 'limma-trend'),
    help='differential expression method to use')
parser$add_argument('--adjust', type='character', default="BH",
    choices=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY"),
    help='method for mutiple test correction')
# https://www.rdocumentation.org/packages/stats/versions/3.5.0/topics/p.adjust
parser$add_argument('--normalize', type='character', default='TMM',
    choices=c('RLE', 'CPM', 'TMM', 'upperquartile'))
parser$add_argument('--output', type='character', required=TRUE,
    help='output differential expression table')
args <- parser$parse_args()


suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))

message(paste0("[",Sys.time(),"] ", "Load counts matrix ..."))
count.matrix <- read.table(args$matrix, header = TRUE, row.names=1, check.names=FALSE, sep='\t')

message(paste0("[",Sys.time(),"] ", "Load metadata ..."))
metadata <- read.table(args$metadata, header = TRUE, row.names=1, check.names=TRUE, sep='\t', stringsAsFactors = T)

sample.ids <- intersect(rownames(metadata),colnames(count.matrix))
message(paste0("[",Sys.time(),"] ", length(sample.ids), " samples appear in both matrix ."))

count.matrix <- count.matrix[,sample.ids]
metadata <- metadata[sample.ids,]

message(paste0("[",Sys.time(),"] ", "Prepare design matrix ..."))
metadata[[args$label_field]] <- relevel(metadata[[args$label_field]],ref=args$control_label) 
regressors <- c(args$label_field)
if("covariate.fields" %in% names(args)){
    covariate.fields <- unlist(strsplit(args$covariate_fields,","))
    regressors <- c(regressors,covariate.fields)
}
design <- model.matrix(as.formula(paste("~", paste(regressors, collapse="+"))),data=metadata)

message(paste0("[",Sys.time(),"] ", "Filter gene by expression value ... " ))
y <- edgeR::DGEList(counts=count.matrix)
keep <- edgeR::filterByExpr(y,group = metadata[[args$label_field]])
y <- y[keep, , keep.lib.sizes=FALSE]
message(paste0("[",Sys.time(),"] ", "Calculate size factor with ",args$normalize," ..."))
y <- edgeR::calcNormFactors(y, method=args$normalize)


coef <- paste0(args$label_field,args$case_label)

message(paste0("[",Sys.time(),"] ", "Test for differential expression with ", args$test ," ..."))

if(grepl("^edger",args$test)){
  y <- edgeR::estimateDisp(y, design)
  if(args$test=="edger-glmqlf"){
    fit <- edgeR::glmQLFit(y, design)
    case.vs.control <- edgeR::glmQLFTest(fit, coef=coef)
  }else if (args$test=="edger-glmlrt"){
    fit <- edgeR::glmFit(y, design)
    case.vs.control <- edgeR::glmLRT(fit, coef=coef)
  }
    case.vs.control.table <- edgeR::topTags(case.vs.control, n=Inf, adjust.method=args$adjust)
}else{
  if(args$test=="limma-trend"){
    log.cpm <- edgeR::cpm(y, log=TRUE, prior.count=3)
    fit <- limma::lmFit(log.cpm, design)
    fit <- limma::eBayes(fit, trend=TRUE)
  }else{
    v <- limma::voom(y, design, plot=FALSE)
    fit <- lmFit(v, design)
    fit <- limma::eBayes(fit, trend=FALSE)
  }
    case.vs.control.table <- limma::topTable(fit, coef=coef,number=Inf,adjust.method=args$adjust)    
}


message(paste0("[",Sys.time(),"] ", "Save results to ", args$output, " ..."))
write.table(case.vs.control.table, args$output, sep='\t', quote=FALSE, row.names=TRUE)
message(paste0("[",Sys.time(),"] ", "All done ."))



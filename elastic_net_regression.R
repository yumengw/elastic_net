#############################################################################
## This is the code to learn pathway network from molecular data through 
## elastic net regression
## -Input data: pathway score, mutation, SCNA, gene expression, DNA methy
## -Output data: weights of selected features in each cancer type
## (c) Copyright The University of Texas MD Anderson Cancer Center,
## (c) Baylor College of Medicine,
## (c) Dr. Han Liang Lab (http://odin.mdacc.tmc.edu/~hliang1/)
## @author Yumeng Wang (yumengw@bcm.edu) 09-10-2017
############################################################################

library(glmnet)
source("preselection.R")
op <- options(warn = (-1)) # suppress warnings 

## cancers
cancer_list = scan("cancer_list_9125", what=character(), sep="\n")

## sample file
sample_list = read.table("/rsrch2/bcb/ywang50/pathway/9125_samples.tsv", sep="\t", header=T, check.names=F, stringsAsFactors=F)

## pathway genes
pathway_gene_list <- scan("/rsrch2/bcb/ywang50/pathway/pathway_gene_list", what=character(), sep="\n")

## generate results dir if not exist
results_dir <- "results/"
if (!dir.exists(results_dir)) 
	dir.create(results_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


############################################################################
## main 
############################################################################

for (cancer in cancer_list) {
	
	## x: matrix of all predictors
    sample_list_cancer <- sample_list[sample_list$DISEASE == cancer, "PATIENT_BARCODE"]

	## y: pathway score
	target_score = read.table(paste("/rsrch2/bcb/ywang50/pathway/elastic_net/targe_score/", cancer, "_pathway_score_9125.tsv", sep=""), header = T, sep = '\t', stringsAsFactors=FALSE, check.names=F)
	target_score <- target_score[match(sample_list_cancer, rownames(target_score)), , drop=F]
	stopifnot((identical(sample_list_cancer, rownames(target_score)))
	y <- as.vector(as.numeric(target_score$TargetZscore))

	## mRNA
	x_mRNA = read.table(paste("/rsrch2/bcb/ywang50/pathway/elastic_net/mRNA/", cancer, "_mRNA_9125.tsv", sep=""), header = T, sep = '\t', stringsAsFactors=FALSE, check.names=F)
	x_mRNA <- x_mRNA[, !colnames(x_mRNA) %in% c("AMOTL2", "LATS2"),drop=F] ## remove target gene
	x_mRNA <- x_mRNA[match(sample_list_cancer, rownames(x_mRNA)), ,drop=F]
	stopifnot(identical(sample_list_cancer, rownames(x_mRNA)))
	keep_link <- preselection(x_mRNA, target_score, "mRNA", cancer)
	x_mRNA <- x_mRNA[, keep_link, drop=F]
	colnames(x_mRNA) <- paste(colnames(x_mRNA), "(mRNA)", sep="")

	## CNV
	x_CNV = read.table(paste("/rsrch2/bcb/ywang50/pathway/elastic_net/copy_number_value_hippo/", cancer, "_cnv_value_9125.tsv", sep=""), header = T, sep = '\t', stringsAsFactors=FALSE, check.names=F)
    x_CNV <- x_CNV[match(sample_list_cancer, rownames(x_CNV)), ,drop=F]
    stopifnot(identical(sample_list_cancer, rownames(x_CNV)))
    keep_link <- preselection(x_CNV, target_score, "cnv")
    x_CNV <- x_CNV[, keep_link, drop=F]
    colnames(x_CNV) <- paste(colnames(x_CNV), "(cnv)", sep="")

	## Mutation
	x_mut = read.table(paste("/rsrch2/bcb/ywang50/pathway/elastic_net/mutation_hippo/", cancer, "_mut_matrix_9125.tsv", sep=""), header = T, sep = '\t', stringsAsFactors=FALSE, check.names=F)
	colnames(x_mut) <- paste(colnames(x_mut), "(mut)", sep="")
    x_mut <- x_mut[match(sample_list_cancer, rownames(x_mut)), ,drop=F]
	x_mut <- apply(x_mut, 2, as.factor) ## imporant, has to be factor
    stopifnot(identical(sample_list_cancer, rownames(x_mut)))

	## Methylation
	x_methy = read.table(paste("/rsrch2/bcb/ywang50/pathway/elastic_net/methylation/process_450K_data/", cancer, "_methy_9125.tsv", sep=""), header = T, sep = '\t', stringsAsFactors=FALSE, check.names=F)
    x_methy <- x_methy[, colnames(x_methy) %in% pathway_gene_list, drop=F]
    x_methy <- x_methy[match(sample_list_cancer, rownames(x_methy)), ,drop=F]
    stopifnot(identical(sample_list_cancer, rownames(x_methy)))
    keep_link <- preselection(x_methy, target_score, "methy")
    x_methy <- x_methy[, keep_link, drop=F]
    colnames(x_methy) <- paste(colnames(x_methy), "(methy)", sep="")

	## normalize binary mutation data into factor
	x_mut_factor <- model.matrix(y ~ ., data=data.frame(x_mut))[, -1, drop=F]
	x <- scale(as.matrix(data.frame(x_mRNA, x_CNV, x_methy, check.names=F))) 
	x <- as.matrix(data.frame(x, x_mut_factor, check.names=F))

    ## fit model
	rep_time <- 100
	mse <- rep(NA, rep_time) # mse
	coef <- matrix(NA, ncol(x)+1, rep_time, dimnames=list(c(), c(1:rep_time)))
	for (ii in c(1:rep_time)) { # do this rep_time times

		train_rows <- sample(1:length(y), 0.8*length(y))
		x.train <- x[train_rows, ]
		y.train <- y[train_rows]

		x.test <- x[-train_rows, ]
		y.test <- y[-train_rows]

		fit <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=0.5, standardize=FALSE, family="gaussian") ## standardize=FALSE
		yhat <- predict(fit, s=fit$lambda.1se, newx=x.test)
		mse[ii] <- mean((y.test - yhat)^2)

		model <- glmnet(x, y, alpha=0.5, family="gaussian", standardize=FALSE)
		coef.apprx <- coef(model, s=fit$lambda.1se, exact=F)
		coef[,ii] <- as.vector(coef.apprx[, 1])
		rownames(coef) <- rownames(coef.apprx)
	}
	write.table(coef, file=paste(results_dir, cancer, "_coeff.tsv", sep=""), sep="\t", quote=F)
	write.table(mse, file=paste(results_dir, cancer, "_mse.tsv", sep=""), sep="\t", quote=F)
	
	## find top predictors being selected in more than 80% of the simulation
	coef <- coef[-1, ]
	count <- apply(coef, 1, function(x) sum(x!=0))
	top_predictor <- names(count[count >= 0.8*rep_time])
	top_coef <- coef[rownames(coef) %in% top_predictor, ,drop=F]
	top_coef[top_coef==0] <- NA
	average_weight <- apply(top_coef, 1, mean, na.rm=T)
	average_weight <- average_weight[order(-average_weight)]
	write.table(data.frame(average_weight), file=paste(results_dir, cancer, "_top_predictor.tsv", sep=""), sep="\t", quote=F)

}

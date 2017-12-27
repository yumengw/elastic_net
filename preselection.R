# given data, calculate correlations with pathway target score and 
# decide if it is at the right direction
# 2017-10-15

preselection <- function(x, y, type) {
	
	keep_colnames <- c()
	pos_corr_cutoff <- 0.3
	neg_corr_cutoff <- -0.3

	stopifnot(identical(rownames(y), rownames(x)))

	y <- as.vector(as.numeric(y$TargetZscore))

	pathway = scan("/rsrch2/bcb/ywang50/pathway/pathway_gene_list", sep="\n", what=character())
	pathway_onco = scan("/rsrch2/bcb/ywang50/pathway/pathway_gene_list_oncogene", sep="\n", what=character())
	pathway_suppressor <- setdiff(pathway, pathway_onco)
	if (type == "mRNA" | type == "cnv") {
		for (colindex in c(1:ncol(x))) {
			tmpx <- as.vector(as.numeric(x[, colindex]))
			corr <- cor.test(tmpx, y, method="spearman")
			#print(colnames(x)[colindex])
			#print(corr$estimate[[1]])
			if (!colnames(x)[colindex] %in% pathway)
				keep_colnames <- c(keep_colnames, colnames(x)[colindex])
			else if (corr$estimate[[1]] > neg_corr_cutoff & 
					colnames(x)[colindex] %in% pathway_onco) 
				keep_colnames <- c(keep_colnames, colnames(x)[colindex])
			else if (corr$estimate[[1]] < pos_corr_cutoff & 
					colnames(x)[colindex] %in% pathway_suppressor)
				keep_colnames <- c(keep_colnames, colnames(x)[colindex])
		}
		return(keep_colnames)
	}

	else if (type == "methy") {
        for (colindex in c(1:ncol(x))) {
            tmpx <- as.vector(as.numeric(x[, colindex]))
            corr <- cor.test(tmpx, y, method="spearman")
			#print(colnames(x)[colindex])
			#print(corr$estimate[[1]])
            if (corr$estimate[[1]] < pos_corr_cutoff & 
					colnames(x)[colindex] %in% pathway_onco) 
                keep_colnames <- c(keep_colnames, colnames(x)[colindex])
            if (corr$estimate[[1]] > neg_corr_cutoff & 
					colnames(x)[colindex] %in% pathway_suppressor)
                keep_colnames <- c(keep_colnames, colnames(x)[colindex])
        }
        return(keep_colnames)
    }

	else if (type == "RPPA") {
		pathway_onco_rppa <- c("YAP", "TAZ", "YAP_pS127")
		pathway_suppressor_rppa <- c("ECADHERIN", "NF2", "AMPKALPHA_pT172", "LKB1")
		for (colindex in c(1:ncol(x))) {
            tmpx <- as.vector(as.numeric(x[, colindex]))
            corr <- cor.test(tmpx, y, method="spearman")
            if (corr$estimate[[1]] > neg_corr_cutoff & 
					colnames(x)[colindex] %in% pathway_onco_rppa)
                keep_colnames <- c(keep_colnames, colnames(x)[colindex])
            if (corr$estimate[[1]] < pos_corr_cutoff & 
					colnames(x)[colindex] %in% pathway_suppressor_rppa)
                keep_colnames <- c(keep_colnames, colnames(x)[colindex])
        }
        return(keep_colnames)
	}
}

## Class declaration
###############################################################################

setClass("TranslatomeDataset", 
	representation(expr.matrix="matrix", 
								 cond.a="character", cond.b="character", 
								 cond.c="character", cond.d="character", 
								 label.level="character", label.condition="character", 
								 data.type="character", DEGs="DEGs")
)


## Generics declaration (getters, setters and methods)
###############################################################################

setGeneric("getExprMatrix",
	function(object) standardGeneric("getExprMatrix"))
	
setGeneric("getConditionA",
	function(object) standardGeneric("getConditionA"))
	
setGeneric("getConditionB",
	function(object) standardGeneric("getConditionB"))
	
setGeneric("getConditionC",
	function(object) standardGeneric("getConditionC"))
	
setGeneric("getConditionD",
	function(object) standardGeneric("getConditionD"))
	
setGeneric("getDataType",
	function(object) standardGeneric("getDataType"))
	
setGeneric("getConditionLabels", 
	function(object) standardGeneric("getConditionLabels"))
	
setGeneric("getLevelLabels",
	function(object) standardGeneric("getLevelLabels"))
	
setGeneric("getDEGs",
	function(object) standardGeneric("getDEGs"))

setGeneric("computeDEGs", signature="object",
	function(object, method="limma", significance.threshold= 0.05, 
					 FC.threshold= 0,	log.transformed = FALSE, mult.cor=TRUE) 
	standardGeneric("computeDEGs"))


## Getters implementation
###############################################################################

setMethod("getExprMatrix", "TranslatomeDataset",
	function(object) { object@expr.matrix }
)

setMethod("getConditionA", "TranslatomeDataset",
	function(object) { object@cond.a }
)

setMethod("getConditionB", "TranslatomeDataset",
	function(object) { object@cond.b }
)

setMethod("getConditionC", "TranslatomeDataset",
	function(object) { object@cond.c }
)

setMethod("getConditionD", "TranslatomeDataset",
	function(object) { object@cond.d }
)

setMethod("getDataType", "TranslatomeDataset",
	function(object) { object@data.type }
)

setMethod("getConditionLabels", "TranslatomeDataset",
	function(object) { object@label.condition }
)

setMethod("getLevelLabels", "TranslatomeDataset",
	function(object) { object@label.level }
)

setMethod("getDEGs", "TranslatomeDataset",
	function(object) { object@DEGs }
)


## Methods implementation
###############################################################################

setMethod("computeDEGs", "TranslatomeDataset", 
	function(object, method="limma", significance.threshold= 0.05, 
					 FC.threshold= 0,	log.transformed = FALSE, mult.cor=TRUE) {
	
		# check input parameters for correctness
		if (method != "none" & 
				(length(object@cond.a)<2 | length(object@cond.b)<2 | 
				 length(object@cond.c)<2 | length(object@cond.d)<2)) 
			stop('You need at least two replicates for each condition 
			to calculate DEGs with statistical methods!')
		
		if (!(method %in% 
					c("limma", "SAM", "t-test", "TE", "RP", "ANOTA", "DESeq", "edgeR", "none"))) 
			stop('This method is not recognized!')
		
		# conditions for the two levels (first is 1,2 and second is 3,4)
		cond1 <- object@expr.matrix[,object@cond.a]
		cond2 <- object@expr.matrix[,object@cond.b]
		cond3 <- object@expr.matrix[,object@cond.c]
		cond4 <- object@expr.matrix[,object@cond.d]

		if (!(object@data.type == "ngs") && !log.transformed) {
			cond1 <- log(cond1, base=2)
			cond2 <- log(cond2, base=2)  
			cond3 <- log(cond3, base=2)
			cond4 <- log(cond4, base=2)
		}	
	
		cond <- cbind(cond1, cond2)
		cond.vector <- c(rep(0, ncol(cond1)), rep(1, ncol(cond2)))

		cond.2 <- cbind(cond3,cond4)
		cond.2.vector <- c(rep(0, ncol(cond3)), rep(1, ncol(cond4)))
		
		# if the chosen method is translational efficiency, put together samples
		# as first & second level case and first & second level control
		# meaning: Pol Case vs Sub/Tot Case and Pol Ctrl vs Sub/Tot Ctrl
		if (method == "TE") {
			cond <- cbind(cond1, cond3)
			cond.vector <- c(rep(0, ncol(cond1)), rep(1, ncol(cond3)))

			cond.2 <- cbind(cond2,cond4)
			cond.2.vector <- c(rep(0, ncol(cond2)), rep(1, ncol(cond4)))
		}
		
		#Calculation of FC, avg, and sd for the first level
		FC <- apply(cond, 1,
					function(x) mean(x[which(cond.vector == 1)], na.rm=TRUE) - 
											mean(x[which(cond.vector == 0)], na.rm=TRUE))
		if (object@data.type  ==  "ngs") 
        FC <- apply(cond, 1,
					function(x) log(mean(x[which(cond.vector == 1)], na.rm=TRUE) / 
													mean(x[which(cond.vector == 0)], na.rm=TRUE), base=2))
		avg.trt <- apply(cond, 1,
								function(x) mean(x[which(cond.vector == 1)], na.rm=TRUE))
		sd.trt <- apply(cond, 1,
								function(x) sd(x[which(cond.vector == 1)], na.rm=TRUE))
		avg.ctr <- apply(cond, 1,
								function(x) mean(x[which(cond.vector == 0)], na.rm=TRUE))
		sd.ctr <- apply(cond, 1,
								function(x) sd(x[which(cond.vector == 0)], na.rm=TRUE))

		#Calculation of FC, avg, and sd for the second level
		FC2 <- apply(cond.2, 1,
						function(x) mean(x[which(cond.2.vector == 1)],na.rm=TRUE) - 
												mean(x[which(cond.2.vector == 0)],na.rm=TRUE))
		if (object@data.type  ==  "ngs") 
      FC2 <- apply(cond.2, 1,
						function(x) log(mean(x[which(cond.2.vector == 1)],na.rm=TRUE) / 
														mean(x[which(cond.2.vector == 0)],na.rm=TRUE),base=2))
    avg.trt2 <- apply(cond.2, 1,
									function(x) mean(x[which(cond.2.vector == 1)], na.rm=TRUE))
		sd.trt2 <- apply(cond.2, 1,
									function(x) sd(x[which(cond.2.vector == 1)], na.rm=TRUE))
		avg.ctr2 <- apply(cond.2, 1,
									function(x) mean(x[which(cond.2.vector == 0)], na.rm=TRUE))
		sd.ctr2 <- apply(cond.2, 1,
									function(x) sd(x[which(cond.2.vector == 0)], na.rm=TRUE))

		#execution of the selected significance computation method 
		#(sig.matrix is set with the default values for "none" method)
		sig.matrix <- matrix(0,nrow=nrow(cond),ncol=4)
		
		if (method == "RP") 
			sig.matrix <- methodRP(cond, cond.2, cond.vector, cond.2.vector, mult.cor)
		if (method == "t-test") 
			sig.matrix <- methodTTest(cond, cond.2, cond.vector, cond.2.vector)
		if (method == "TE")
			# compute translational efficiency p-values with Limma as if it was
			# the normal condition, but cond and cond.2 have been built in a
			# different way (tot/sub + pol ctrl) and (tot/sub + pol case)
			sig.matrix <- methodLimma(cond, cond.2, cond.vector, cond.2.vector)
		if (method == "SAM") 
			sig.matrix <- methodSAM(cond, cond.2, cond.vector, cond.2.vector)
		if (method == "limma") 
			sig.matrix <- methodLimma(cond, cond.2, cond.vector, cond.2.vector)
		if (method == "ANOTA") 
			sig.matrix <- methodANOTA(cond, cond.2, cond.vector, cond.2.vector)
		if (method == "DESeq") 
			sig.matrix <- methodDESeq(cond, cond.2, cond.vector, cond.2.vector)
		if (method == "edgeR") 
			sig.matrix <- methodEdgeR(cond, cond.2, cond.vector, cond.2.vector)

		#generation of the final matrix
		
		col1 <- ifelse(mult.cor, 2, 1)
		col2 <- ifelse(mult.cor, 4, 3)

		level1Up <- sig.matrix[,col1] < significance.threshold & FC > FC.threshold
		level1Down <- sig.matrix[,col1] < significance.threshold & FC < -FC.threshold
		level2Up <- sig.matrix[,col2] < significance.threshold & FC2 > FC.threshold
		level2Down <- sig.matrix[,col2] < significance.threshold & FC2 < -FC.threshold

		UpUp <- level1Up == 1 & level2Up == 1			 
		DownUp <- level1Down == 1 & level2Up == 1
		UpDown <- level1Up == 1 & level2Down == 1
		DownDown <- level1Down == 1 & level2Down == 1

		final.matrix <- cbind(FC, avg.ctr, sd.ctr, 
													avg.trt, sd.trt, sig.matrix[,1:2],
													FC2, avg.ctr2, sd.ctr2, 
													avg.trt2, sd.trt2, sig.matrix[,3:4],
													level1Up, level1Down, level2Up, level2Down,
													UpUp, DownUp, UpDown, DownDown)
													
		colnames(final.matrix) <- c(
			LOG2.FC1=paste("log2 FC (", 
										 object@label.level[1], ")", sep=""), 
			AVG11=paste("avg(", object@label.level[1], ", ", 
									object@label.condition[1], ")", sep=""),  
			SD11=paste("sd(", object@label.level[1], ", ",
								 object@label.condition[1], ")", sep=""), 
			AVG12=paste("avg(", object@label.level[1], ", ",
									object@label.condition[2], ")", sep=""),  
			SD12=paste("sd(", object@label.level[1], ", ", 
								 object@label.condition[2], ")", sep=""), 
			PVAL.1=paste(method, ".pval(", 
									 object@label.level[1], ")", sep=""),  
			PVAL.MTC.1=paste(method, ".pval.mtc(", 
											 object@label.level[1], ")", sep=""), 
			LOG2.FC2=paste("log2 FC (", 
										 object@label.level[2], ")", sep=""), 
			AVG21=paste("avg(", object@label.level[2], ", ", 
									object@label.condition[1], ")", sep=""),  
			SD21=paste("sd(", object@label.level[2], ", ", 
								 object@label.condition[1], ")", sep=""), 
			AVG22=paste("avg(", object@label.level[2], ", ", 
									object@label.condition[2], ")", sep=""),  
			SD22=paste("sd(", object@label.level[2], ", ", 
								 object@label.condition[2], ")", sep=""), 
			PVAL.2=paste(method, ".pval(", 
									 object@label.level[2], ")", sep=""),  
			PVAL.MTC.2=paste(method, ".pval.mtc(", 
											 object@label.level[2], ")", sep=""), 
			"level1Up",  "level1Down",  "level2Up",  "level2Down",  
			"UpUp",  "DownUp",  "UpDown",  "DownDown")

		# store and return the obtained DEGs
		object@DEGs <- new("DEGs", DEGs.table=final.matrix, method=method, 
											 significance.threshold=significance.threshold, 
											 FC.threshold=FC.threshold, 
											 label.level.DEGs=object@label.level, 
											 label.condition=object@label.condition)
		return(object@DEGs)
	}
)


setMethod("show", "TranslatomeDataset",
	function(object) {
		print(head(object@expr.matrix))
		print(object@cond.a)
		print(object@cond.b)
		print(object@cond.c)
		print(object@cond.d)
		print(object@label.condition)
		print(object@label.level)
		print(object@data.type)
		print(object@DEGs)
	}
)


## Helper functions
###############################################################################

# Constructor method for users
newTranslatomeDataset <- function(expr.matrix, cond.a, cond.b, cond.c, cond.d,
																	 data.type="ngs", 
																	 label.level=c("1st level","2nd level"), 
																	 label.condition=c("control","treated")) {

	# check input parameters for completeness and correctness
  if (missing(expr.matrix) | missing(cond.a) | missing(cond.b) | 
			missing(cond.c) | missing(cond.d)) 
		stop('Some of the mandatory arguments are missing!')
  
  # if the input dataset is a Biobase ExpressionSet, 
  # extract the expression matrix from it
  finalMatrix = expr.matrix
  if (class(expr.matrix) == "ExpressionSet") finalMatrix = exprs(expr.matrix)
  
	return(new(Class="TranslatomeDataset",expr.matrix = finalMatrix, 
						 cond.a=cond.a, cond.b=cond.b, cond.c=cond.c, cond.d=cond.d, 
             data.type=data.type, label.level=label.level, 
             label.condition=label.condition))
}

# Implementation of the RP helper function

methodRP <- function(cond, cond.2, cond.vector, cond.2.vector, mult.cor) {
	
	num.perm = ifelse(mult.cor, 1000, 10)

	RP <- RP(cond, cond.vector, num.perm=num.perm, 
					 logged=TRUE, gene.names=rownames(cond))
	rp.pval.1 <- apply(RP$pval, 1, min)
	rp.pfp.1 <- apply(RP$pfp, 1, min)
	
	RP2 <- RP(cond.2, cond.2.vector, num.perm=num.perm, 
						logged=TRUE, gene.names=rownames(cond.2))
	rp.pval.2 <- apply(RP2$pval, 1, min)
	rp.pfp.2 <-  apply(RP2$pfp, 1, min)
	
	# build the significance p-values matrix and return it	
	return(cbind(rp.pval.1, rp.pfp.1, rp.pval.2, rp.pfp.2))
}

# Implementation of the t-test helper function
methodTTest <- function(cond, cond.2, cond.vector, cond.2.vector) {
	
	# t.test.pval for first level
	t.test.pval <- calcTStatFast(cond, cond.vector)$pval
	t.test.pval.adj <- p.adjust(t.test.pval, 
															method="BH", n=length(t.test.pval))

	# t.test.pval for second level
	t.test.pval2 <- calcTStatFast(cond.2, cond.2.vector)$pval
	t.test.pval.adj2 <- p.adjust(t.test.pval2, 
															 method="BH", n=length(t.test.pval2))

	# build the significance p-values matrix and return it	
	return(cbind(t.test.pval, t.test.pval.adj, 
							 t.test.pval2, t.test.pval.adj2))	
}


# Implementation of the SAM helper function
methodSAM <- function(cond, cond.2, cond.vector, cond.2.vector) {
	
	data.tot <- list(x=cond, y=(cond.vector+1), logged2=TRUE)
	sam <- samr(data.tot, resp.type="Two class unpaired", nperms=1000)
	pv.sam <- samr.pvalues.from.perms(sam$tt, sam$ttstar)
	pv.sam.adj <- p.adjust(pv.sam, method="BH", n=length(pv.sam))

	data.tot2 <- list(x=cond.2, y=(cond.2.vector+1), logged2=TRUE)
	sam2 <- samr(data.tot2, resp.type="Two class unpaired", nperms=1000)
	pv.sam2 <- samr.pvalues.from.perms(sam2$tt, sam2$ttstar)
	pv.sam.adj2 <- p.adjust(pv.sam2, method="BH", n=length(pv.sam2))

	# build the significance p-values matrix and return it	
	return(cbind(pv.sam, pv.sam.adj, pv.sam2, pv.sam.adj2))
}

# Implementation of the limma helper function
methodLimma <- function(cond, cond.2, cond.vector, cond.2.vector) {
	
	eset <- data.frame(cond.a=cond[,cond.vector==0], 
										 cond.b=cond[,cond.vector==1], 
										 cond.c=cond.2[,cond.2.vector==0], 
										 cond.d=cond.2[,cond.2.vector==1])
										  	
	Target <- c(rep("cond.a", sum(cond.vector==0)), 
							rep("cond.b", sum(cond.vector==1)), 
							rep("cond.c", sum(cond.2.vector==0)), 
							rep("cond.d", sum(cond.2.vector==1)))
							
	colnames(eset) <- Target
	FileName <- colnames(eset)
	targets <- data.frame(FileName=FileName, Target=Target)

	lev <- unique(Target)
	f <- factor(targets$Target, levels = lev)
	design <- model.matrix(~0 + f)
	colnames(design) <- lev
	fit <- lmFit(eset, design)

	cont.firstlevel <- makeContrasts("cond.b - cond.a", levels = design)
	fitfirstlevel <- contrasts.fit(fit, cont.firstlevel)
	fitfirstlevel <- eBayes(fitfirstlevel)
	BH.firstlevel <- p.adjust(fitfirstlevel$F.p.value, method = "BH") 

	cont.secondlevel <- makeContrasts("cond.d - cond.c", levels = design)
	fitsecondlevel <- contrasts.fit(fit, cont.secondlevel)
	fitsecondlevel <- eBayes(fitsecondlevel)
	BH.secondlevel <- p.adjust(fitsecondlevel$F.p.value, method = "BH") 
		
	# build the significance p-values matrix and return it	
	return(cbind(fitfirstlevel$F.p.value, BH.firstlevel, 
							 fitsecondlevel$F.p.value, BH.secondlevel))
}

# Implementation of the ANOTA helper function
methodANOTA <- function(cond, cond.2, cond.vector, cond.2.vector) {
	#library(anota)
	
	if (length(cond.vector) != length(cond.2.vector) | 
			!all(cond.vector == cond.2.vector))
		stop('The ANOTA method cannot be applied if the two levels 
		have a different experimental design!')
	
	#first level
	anotaQcOut <- anotaPerformQc(dataT = cond, dataP = cond.2, 
															 phenoVec = cond.vector)
	anotaSigGeneOut <- anotaGetSigGenes(dataT = cond, dataP = cond.2, 
															 phenoVec = cond.vector, anotaQcObj = anotaQcOut)

	pvalues.1st <- anotaSigGeneOut$apvStats[[1]][,"apvP"]
	pvalues.1st.BH <- p.adjust(pvalues.1st, method = "BH") 

	#second level
	anotaQcOut.2nd <- anotaPerformQc(dataT = cond.2, dataP = cond, 
																	 phenoVec = cond.vector)
	anotaSigGeneOut.2nd <- anotaGetSigGenes(dataT = cond.2, dataP = cond, 
																	 phenoVec = cond.vector, 
																	 anotaQcObj = anotaQcOut.2nd)

	pvalues.2nd <- anotaSigGeneOut.2nd$apvStats[[1]][,"apvP"]
	pvalues.2nd.BH <- p.adjust(pvalues.2nd, method = "BH") 

	# build the significance p-values matrix and return it	
	return(cbind(pvalues.2nd, pvalues.2nd.BH, pvalues.1st, pvalues.1st.BH))
}

# Implementation of the DESeq helper function (for NGS data)
methodDESeq <- function(cond, cond.2, cond.vector, cond.2.vector) {
	#require(DESeq)
	
	cond.1.deseq <- newCountDataSet(cond, cond.vector)
	cond.1.deseq <- estimateSizeFactors(cond.1.deseq)
  # if GLM fit fails, switch to local fitting
	tmp <- tryCatch(cond.1.deseq <- estimateDispersions(cond.1.deseq),
								 error=function(ex){ return(-1) })
	if (is.numeric(tmp))
		cond.1.deseq <- estimateDispersions(cond.1.deseq, 
																			fitType="local",
																			sharingMode="fit-only")	
	res.1 <- nbinomTest(cond.1.deseq, "0", "1")
	
	cond.2.deseq <- newCountDataSet(cond.2, cond.2.vector)
	cond.2.deseq <- estimateSizeFactors(cond.2.deseq)
  # if GLM fit fails, switch to local fitting
  tmp <- tryCatch(cond.2.deseq <-estimateDispersions(cond.2.deseq), 
								 error=function(ex){ return(-1) })
	if (is.numeric(tmp))
		cond.2.deseq <- estimateDispersions(cond.2.deseq, 
																			fitType="local",
																			sharingMode="fit-only")	
	res.2 <- nbinomTest(cond.2.deseq, "0", "1")
	
	# build the significance p-values matrix and return it	
	return(cbind(res.1$pval, res.1$padj, res.2$pval, res.2$padj))
}

# Implementation of the edgeR helper function (for NGS data)
methodEdgeR <- function(cond, cond.2, cond.vector, cond.2.vector) {
	#require(edgeR)
	
	cond.1.edger <- DGEList(counts=cond, group=cond.vector)
	cond.1.edger <- calcNormFactors(cond.1.edger)
	cond.1.edger <- estimateCommonDisp(cond.1.edger)
	cond.1.edger <- estimateTagwiseDisp(cond.1.edger)
	res.1 <- exactTest(cond.1.edger)
	
	cond.2.edger <- DGEList(counts=cond.2, group=cond.2.vector)
	cond.2.edger <- calcNormFactors(cond.2.edger)
	cond.2.edger <- estimateCommonDisp(cond.2.edger)
	cond.2.edger <- estimateTagwiseDisp(cond.2.edger)
	res.2 <- exactTest(cond.2.edger)
	
	# build the significance p-values matrix and return it	
	return(cbind(res.1$table$PValue, 
							 p.adjust(res.1$table$PValue,
												method="BH",n=length(res.1$table$PValue)), 
               res.2$table$PValue, 
               p.adjust(res.2$table$PValue,
												method="BH",n=length(res.2$table$PValue))))
}
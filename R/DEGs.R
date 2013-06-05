## Class declaration
###############################################################################

setClass("DEGs",
	representation(
		DEGs.table="matrix", method="character", 
		significance.threshold="numeric", FC.threshold="numeric",
		label.level.DEGs="character", label.condition="character"
	)
)


## Generics declaration (getters, setters and methods)
###############################################################################

setGeneric("DEGs.table", signature="object", 
					function(object) standardGeneric("DEGs.table"))
					
setGeneric("getDEGsMethod", signature="object", 
					function(object) standardGeneric("getDEGsMethod"))

setGeneric("significance.threshold", signature="object", 
					function(object) standardGeneric("significance.threshold"))
					
setGeneric("FC.threshold", signature="object", 
					function(object) standardGeneric("FC.threshold"))
					
setGeneric("label.level.DEGs", signature="object", 
					function(object) standardGeneric("label.level.DEGs"))
					
setGeneric("label.condition", signature="object", 
					function(object) standardGeneric("label.condition"))

setGeneric("MAplot", signature="object", 
					function(object, outputformat="on screen",track="") 
					standardGeneric("MAplot"))
					
setGeneric("CVplot", signature="object", 
					function(object, outputformat="on screen",track="") 
					standardGeneric("CVplot"))
					
setGeneric("SDplot", signature="object", 
					function(object, outputformat="on screen",track="") 
					standardGeneric("SDplot"))
					
setGeneric("Histogram", signature="object", 
					function(object, plottype="summary", outputformat="on screen") 
					standardGeneric("Histogram"))
					
setGeneric("Scatterplot", signature="object",
					function(object, outputformat="on screen",track="") 
					standardGeneric("Scatterplot"))
					
setGeneric("GOEnrichment", signature="object", 
					function(object, ontology="all", 
									classOfDEGs="both", test.method="classic",
									test.threshold = 0.05, mult.cor=TRUE) 
					standardGeneric("GOEnrichment"))
		
setGeneric("RegulatoryEnrichment", signature="object", 
					function(object, classOfDEGs="both", 
									 significance.threshold = 0.05, mult.cor=TRUE)
					standardGeneric("RegulatoryEnrichment"))


## Getters implementation
###############################################################################

setMethod("DEGs.table", "DEGs",
	function(object) {
		object@DEGs.table
	}
)

setMethod("getDEGsMethod", "DEGs",
	function(object) {
		object@method
	}
)

setMethod("significance.threshold", "DEGs",
	function(object) {
		object@significance.threshold
	}
)

setMethod("FC.threshold", "DEGs",
	function(object) {
		object@FC.threshold
	}
)

setMethod("label.level.DEGs", "DEGs",
	function(object) {
		object@label.level.DEGs
	}
)

setMethod("label.condition", "DEGs",
	function(object) {
		object@label.condition
	}
)


## Methods implementation
###############################################################################

## Implementation of the show method
setMethod("show", "DEGs",
	function(object) {
		print(head(object@DEGs.table))
		print(object@method)
		print(object@significance.threshold)
		print(object@FC.threshold)
		print(object@label.condition)
		print(object@label.level.DEGs)
	}
)


## Implementation of the CVplot method
setMethod("CVplot", "DEGs",
	function(object,outputformat="on screen",track="") {
		resultMatrix <- DEGs.table(object)
		a <- grep(c("FC"), colnames(resultMatrix),value=TRUE)
		c <- grep(c("avg"), colnames(resultMatrix),value=TRUE)
		b <- grep(c("sd"), colnames(resultMatrix),value=TRUE)

		vector1 <- rowMeans(resultMatrix[,b[1:2]]) / 
												abs(rowMeans(resultMatrix[,c[1:2]]))
		vector2 <- resultMatrix[,a[1]] 

		vector3 <- rowMeans(resultMatrix[,b[3:4]]) /
												abs(rowMeans(resultMatrix[,c[3:4]]))
		vector4 <- resultMatrix[,a[2]] 

		xlabel <- "CV"
		ylabel <- "log2 FC"

		if (outputformat == "pdf") pdf(file="CVplot.pdf")
		if (outputformat == "postscript") postscript(file="CVplot.ps")
		if (outputformat == "jpeg") jpeg(filename="CVplot.jpeg")

		par(fig=c(0,1,0,1), mar=c(4,4,1,2), mgp=c(2, 0.75, 0))	
		layout(matrix(c(1,2), 2, 1, byrow = TRUE))
		plot(vector1, vector2, xlab=xlabel, ylab=ylabel, 
					col="grey", pch=20, frame.plot=TRUE,cex=0.8, 
					main=substr(a[1], 10, nchar(a[1])-1))

		list1 <- rownames(resultMatrix)[which(resultMatrix[,"level1Up"]==1 | 
																					resultMatrix[,"level1Down"]==1)]
		variable <- which(rownames(resultMatrix) %in% list1)
		color <- "steelblue3"
		points(vector1[variable],vector2[variable],col=color,pch=20,cex=0.8)
		abline(h=0,lty=2,col="darkgray")
		legend("topright", c("DEGs"),fill=c(color), bty="n", cex=0.8)

		listm <- rownames(resultMatrix)[rownames(resultMatrix)%in%track]
		if (length(listm)>0) {
			variable <- which(rownames(resultMatrix) %in% track)
			points(vector1[variable],vector2[variable],col="white",pch=4,cex=0.5)
			text(vector1[variable],vector2[variable],labels=listm,pos=3,cex=0.7)
		}

		plot(vector3, vector4, xlab=xlabel, ylab=ylabel, 
					col="grey", pch=20, frame.plot=TRUE,cex=0.8,
					main=substr(a[2], 10, nchar(a[2])-1))
		abline(h=0,lty=2,col="darkgray")

		list2 <- rownames(resultMatrix)[which(resultMatrix[,"level2Up"]==1 | 
																					resultMatrix[,"level2Down"]==1)]
		variable <- which(rownames(resultMatrix) %in% list2)
		color <- "gold1"
		points(vector3[variable],vector4[variable],col=color,pch=20,cex=0.8)
		abline(h=0,lty=2,col="darkgray")
		legend("topright", c("DEGs"), fill=c(color), bty="n", cex=0.8)

		listm <- rownames(resultMatrix)[rownames(resultMatrix) %in% track]
		if (length(listm) > 0) {
			variable <- which(rownames(resultMatrix) %in% track)
			points(vector3[variable],vector4[variable],col="white",pch=4,cex=0.5)
			text(vector3[variable],vector4[variable],labels=listm,pos=3,cex=0.7)
		}

		# if a file device is open, close it.
		if (!(outputformat == "on screen")) dev.off()
	}
)


## Implementation of the MAplot method
setMethod("MAplot", "DEGs",
	function(object,outputformat="on screen",track="") {
		resultMatrix <- DEGs.table(object)
		a <- grep(c("FC"),colnames(resultMatrix),value=TRUE)
		b <- grep(c("avg"),colnames(resultMatrix),value=TRUE)

		vector1 <- rowMeans(resultMatrix[,b[1:2]])
		vector2 <- resultMatrix[,a[1]] 

		vector3 <- rowMeans(resultMatrix[,b[3:4]])
		vector4 <- resultMatrix[,a[2]] 

		xlabel="log2 avg signal (A)"
		ylabel="log2 FC (M)"

		if (outputformat == "pdf") pdf(file="MAplot.DEGs.pdf")
		if (outputformat == "postscript") postscript(file="MAplot.DEGs.ps")
		if (outputformat == "jpeg") jpeg(filename="MAplot.DEGs.jpeg")

		par(fig=c(0,1,0,1), mar=c(4,4,1,2), mgp=c(2, 0.75, 0))	
		layout(matrix(c(1,2), 2, 1, byrow = TRUE))
		plot(vector1, vector2, xlab=xlabel, ylab=ylabel, 
					col="grey", pch=20, frame.plot=TRUE, cex=0.8, 
					main=substr(a[1], 10, nchar(a[1])-1))

		list1 <- rownames(resultMatrix)[which(resultMatrix[,"level1Up"]==1 | 
																					resultMatrix[,"level1Down"]==1)]
		variable <- which(rownames(resultMatrix) %in% list1)
		color <- "steelblue3"
		points(vector1[variable],vector2[variable],col=color,pch=20,cex=0.8)
		abline(h=0,lty=2,col="darkgray")
		legend("topright", c("DEGs"), fill=c(color), bty="n", cex=0.8)

		listm <- rownames(resultMatrix)[rownames(resultMatrix) %in% track]
		if (length(listm)>0) {
			variable <- which(rownames(resultMatrix) %in% track)
			points(vector1[variable],vector2[variable],col="white",pch=4,cex=0.5)
			text(vector1[variable],vector2[variable],labels=listm,pos=3,cex=0.7)
		}

		plot(vector3, vector4, xlab=xlabel, ylab=ylabel, 
					col="grey", pch=20, frame.plot=TRUE,cex=0.8,
					main=substr(a[2], 10, nchar(a[2])-1))
		abline(h=0,lty=2,col="darkgray")

		list2 <- rownames(resultMatrix)[which(resultMatrix[,"level2Up"]==1 | 
																					resultMatrix[,"level2Down"]==1)]
		variable <- which(rownames(resultMatrix) %in% list2)
		color <- "gold1"
		points(vector3[variable],vector4[variable],col=color,pch=20,cex=0.8)
		abline(h=0,lty=2,col="darkgray")
		legend("topright", c("DEGs"), fill=c(color), bty="n", cex=0.8)

		listm <- rownames(resultMatrix)[rownames(resultMatrix) %in% track]
		if (length(listm)>0) {
			variable <- which(rownames(resultMatrix) %in% track)
			points(vector3[variable], vector4[variable], col="white", pch=4, cex=0.5)
			text(vector3[variable], vector4[variable], labels=listm, pos=3, cex=0.7)
		}

		# if a file device is open, close it.
		if (!(outputformat == "on screen")) dev.off()
	}
)


## Implementation of the Histogram method
setMethod("Histogram", "DEGs",
	function(object, plottype="summary", outputformat="on screen") {

		resultMatrix <- DEGs.table(object)
		a <- grep(c("FC"), colnames(resultMatrix), value=TRUE)
		list1stlevelUP <- rownames(resultMatrix)[
															which(resultMatrix[,"level1Up"]==1)]
		list1stlevelDOWN <- rownames(resultMatrix)[
															which(resultMatrix[,"level1Down"]==1)]

		list2ndlevelUP <- rownames(resultMatrix)[which(resultMatrix[,"level2Up"]==1)]
		list2ndlevelDOWN <- rownames(resultMatrix)[which(resultMatrix[,"level2Down"]==1)]

		listUPUP <- rownames(resultMatrix)[which(resultMatrix[,"UpUp"]==1)]
		listUPDOWN <- rownames(resultMatrix)[which(resultMatrix[,"UpDown"]==1)]
		listDOWNUP <- rownames(resultMatrix)[which(resultMatrix[,"DownUp"]==1)]
		listDOWNDOWN <- rownames(resultMatrix)[which(resultMatrix[,"DownDown"]==1)]

		if (outputformat == "pdf") pdf(file="Histogram.DEGs.pdf")
		if (outputformat == "postscript") postscript(file="Histogram.DEGs.ps")
		if (outputformat == "jpeg") jpeg(filename="Histogram.DEGs.jpeg")
			
		if (plottype=="summary") {
			inputforsummaryplot <- matrix(c(length(list1stlevelUP),
																			length(list1stlevelDOWN), 
																			length(list2ndlevelUP), 
																			length(list2ndlevelDOWN),
																			0, 0), nrow=2, byrow=FALSE)

			par(fig=c(0, 1, 0, 1), mar=c(4, 5, 4, 2), mgp=c(2, 0.75, 0))
			barplot(as.table(inputforsummaryplot), xlab="number of DEGs",
							names.arg=c(substr(a[1], 10, nchar(a[1])-1), 
													substr(a[2], 10, nchar(a[2])-1),""), 
							col=c("grey25","grey75"), las=1, 
							border=NA, horiz=TRUE, beside=FALSE) 
			
			legend("topright", c("up", "down"), fill=c("grey25", "grey75"), bty="n")
		}

		if (plottype == "detailed") {

			inputfordetailedplot <- matrix(
				c(length(list1stlevelUP) - (length(listUPUP) + length(listUPDOWN)),
				length(list1stlevelDOWN) - (length(listDOWNDOWN) + length(listDOWNUP)),
				length(list2ndlevelUP) - (length(listUPUP) + length(listDOWNUP)),
				length(list2ndlevelDOWN) - (length(listDOWNDOWN) + length(listUPDOWN)),
				length(listUPUP), length(listDOWNDOWN), 
				length(listUPDOWN), length(listDOWNUP)),
				nrow=1, byrow=FALSE)
				
			par(fig=c(0,1,0,1), mar=c(4,7,4,2), mgp=c(2, 0.75, 0))
			barplot(as.table(inputfordetailedplot),
				names.arg= c("up / -","down / -","- / up","- / down", 
									   "up / up","down / down","up / down","down / up"),
				las=1,border=NA,horiz=TRUE,
				col=c("steelblue","steelblue4","gold","orange", 
							"chartreuse","chartreuse4","red","red3"),
				beside=TRUE,xlab="number of DEGs")
				
			mtext(paste(substr(a[1], 10, nchar(a[1])-1),
									substr(a[2], 10, nchar(a[2])-1), sep="/"),
						adj=0, line=0, cex=1)
		}

		# if a file device is open, close it.
		if (!(outputformat == "on screen")) dev.off()
	}
)
  
  
## Implementation of the Scatterplot method  
setMethod("Scatterplot", "DEGs",
	function(object, outputformat="on screen", track="") {
	
		resultMatrix <- DEGs.table(object)
		a <- grep(c("FC"), colnames(resultMatrix), value=TRUE)
		 
		vector1 <- resultMatrix[,a[1]]
		vector2 <- resultMatrix[,a[2]]

		originalFCvalX <- 2 ^ vector1
		originalFCvalY <- 2 ^ vector2

		xlabel=substr(a[1], 6, nchar(a[1]))
		ylabel=substr(a[2], 6, nchar(a[2]))

		if (outputformat == "pdf") pdf(file="Scatterplot.DEGs.pdf")
		if (outputformat == "postscript") postscript(file="Scatterplot.DEGs.ps")
		if (outputformat == "jpeg") jpeg(filename="Scatterplot.DEGs.jpeg")

		plot(originalFCvalX, originalFCvalY, xlab=xlabel, ylab=ylabel, 
				 col="grey", pch=20, frame.plot=TRUE, cex=0.8, log="xy")

		mtext(paste("Spearman all:", 
					round(cor.test(vector1, vector2, method="spearman")$estimate, 2),
					" "),
					side=1, adj=1, line=-2.5, cex=0.8)

		list0 <- rownames(resultMatrix)[c(which(resultMatrix[,"level1Up"]==1), 
																			which(resultMatrix[,"level1Down"]==1),
																			which(resultMatrix[,"level2Up"]==1), 
																			which(resultMatrix[,"level2Down"]==1))]
		variable <- which(rownames(resultMatrix)%in%list0)

		if (length(list0)>1)
			mtext(paste("Spearman DEGs:",
						round(cor.test(vector1[list0], 
													 vector2[list0],
													 method="spearman")$estimate, 2),
						" "), 
						side=1, adj=1, line=-1.5, cex=0.8)

		leg <- NULL
		leg.col <- NULL
		list1 <- rownames(resultMatrix)[which(resultMatrix[,"level1Up"]==1 | 
																					resultMatrix[,"level1Down"]==1)]
		variable <- which(rownames(resultMatrix) %in% list1)
		color <- "steelblue3"
		points(originalFCvalX[variable], originalFCvalY[variable],
					 col=color, pch=20, cex=0.8)
		leg <- c(leg, paste(substr(a[1], 10, nchar(a[1])-1), "only"))
		leg.col <- c(leg.col, color)

		list2 <- rownames(resultMatrix)[which(resultMatrix[,"level2Up"]==1 | 
																					resultMatrix[,"level2Down"]==1)]
		variable <- which(rownames(resultMatrix) %in% list2)
		color <- "gold1"
		points(originalFCvalX[variable], originalFCvalY[variable], 
					 col=color, pch=20, cex=0.8)
		leg <- c(leg, paste(substr(a[2], 10, nchar(a[2])-1), "only"))
		leg.col <- c(leg.col, color)

		list3 <- rownames(resultMatrix)[which(resultMatrix[,"UpUp"]==1 |
																					resultMatrix[,"DownDown"]==1)]
		variable <- which(rownames(resultMatrix) %in% list3)
		color <- "chartreuse3"
		points(originalFCvalX[variable], originalFCvalY[variable], 
					 col=color, pch=20, cex=0.8)
		leg <- c(leg, "homodirectional")
		leg.col <- c(leg.col, color)

		list4 <- rownames(resultMatrix)[which(resultMatrix[,"UpDown"]==1 | 
																					resultMatrix[,"DownUp"]==1)]
		variable <- which(rownames(resultMatrix) %in% list4)
		color <- "red2"
		points(originalFCvalX[variable], originalFCvalY[variable],
					 col=color, pch=20, cex=0.8)
		leg <- c(leg, "opposite change")
		leg.col <- c(leg.col, color)

		abline(h=1, lty=2, col="darkgray")
		abline(v=1, lty=2, col="darkgray")
		legend("topleft", leg, fill=leg.col, bty="n")

		listm <- rownames(resultMatrix)[rownames(resultMatrix) %in% track]
		if (length(listm) > 0) {
			variable <- which(rownames(resultMatrix) %in% track)
			points(originalFCvalX[variable], originalFCvalY[variable], 
						 col="white", pch=4, cex=0.5)
			text(originalFCvalX[variable], originalFCvalY[variable], 
					 labels=listm, pos=3, cex=0.7)
		}

		# if a file device is open, close it.
		if (!(outputformat == "on screen")) dev.off()
	}
)


## Implementation of the SDplot method
setMethod("SDplot", "DEGs",
	function(object, outputformat="on screen", track="") {
	
		resultMatrix <- DEGs.table(object)
		a <- grep(c("FC"), colnames(resultMatrix), value=TRUE)
		b <- grep(c("sd"), colnames(resultMatrix), value=TRUE)

		vector1 <- rowMeans(resultMatrix[,b[1:2]])
		vector2 <- resultMatrix[,a[1]] 

		vector3 <- rowMeans(resultMatrix[,b[3:4]])
		vector4 <- resultMatrix[,a[2]] 

		xlabel="SD"
		ylabel="log2 FC"

		if (outputformat == "pdf") pdf(file="SDplot.DEGs.pdf")
		if (outputformat == "ps") postscript(file="SDplot.DEGs.ps")
		if (outputformat == "jpeg") jpeg(filename="SDplot.DEGs.jpeg")
		
		par(fig=c(0,1,0,1), mar=c(4,4,1,2), mgp=c(2, 0.75, 0))	
		layout(matrix(c(1,2), 2, 1, byrow = TRUE))
		plot(vector1, vector2, xlab=xlabel, ylab=ylabel, 
				 col="grey", pch=20, frame.plot=TRUE,cex=0.8,
				 main=substr(a[1], 10, nchar(a[1])-1))

		list1 <- rownames(resultMatrix)[which(resultMatrix[,"level1Up"]==1 | 
																					resultMatrix[,"level1Down"]==1)]
		variable <- which(rownames(resultMatrix) %in% list1)
		color <- "steelblue3"
		points(vector1[variable], vector2[variable], 
					 col=color, pch=20, cex=0.8)
		abline(h=0, lty=2, col="darkgray")
		legend("topright", c("DEGs"), fill=c(color), bty="n", cex=0.8)

		listm <- rownames(resultMatrix)[rownames(resultMatrix) %in% track]
		if (length(listm) > 0) {
			variable <- which(rownames(resultMatrix) %in% track)
			points(vector1[variable], vector2[variable], col="white", pch=4, cex=0.5)
			text(vector1[variable], vector2[variable], labels=listm, pos=3, cex=0.7)
		}

		plot(vector3, vector4, xlab=xlabel, ylab=ylabel, 
				 col="grey", pch=20, frame.plot=TRUE, cex=0.8,
				 main=substr(a[2], 10, nchar(a[2])-1))
		abline(h=0, lty=2, col="darkgray")

		list2 <- rownames(resultMatrix)[which(resultMatrix[,"level2Up"]==1 | 
																					resultMatrix[,"level2Down"]==1)]
		variable <- which(rownames(resultMatrix) %in% list2)
		color <- "gold1"
		points(vector3[variable], vector4[variable],
					 col=color, pch=20, cex=0.8)
		abline(h=0, lty=2, col="darkgray")
		legend("topright", c("DEGs"), fill=c(color), bty="n", cex=0.8)

		listm <- rownames(resultMatrix)[rownames(resultMatrix) %in% track]
		if (length(listm) > 0) {
			variable <- which(rownames(resultMatrix) %in% track)
			points(vector3[variable], vector4[variable], 
						 col="white", pch=4, cex=0.5)
			text(vector3[variable], vector4[variable], labels=listm, pos=3, cex=0.7)
		}

		# if a file device is open, close it.
		if (!(outputformat=="on screen")) dev.off()
	}
)


## Implementation of the GOenrichment method
setMethod("GOEnrichment", "DEGs",
	function(object, ontology="all", classOfDEGs="both", test.method="classic", 
					 test.threshold = 0.05, mult.cor=TRUE) {
		
		resultMatrix <- DEGs.table(object)
		a <- grep(c("FC"), colnames(resultMatrix), value=TRUE)

		myInterestedGenes1stlevel <- rownames(resultMatrix)[
																	which(resultMatrix[,"level1Up"]==1 | 
																				resultMatrix[,"level1Down"]==1)]
		myInterestedGenes2ndlevel <- rownames(resultMatrix)[
																	which(resultMatrix[,"level2Up"]==1 | 
																				resultMatrix[,"level2Down"]==1)]
		
		if (classOfDEGs == "up") {
			myInterestedGenes1stlevel <- rownames(resultMatrix)[
																		which(resultMatrix[,"level1Up"]==1)]
			myInterestedGenes2ndlevel <- rownames(resultMatrix)[
																		which(resultMatrix[,"level2Up"]==1)]
		}
			
		if (classOfDEGs == "down") {
			myInterestedGenes1stlevel <- rownames(resultMatrix)[
																		which(resultMatrix[,"level1Down"]==1)]
			myInterestedGenes2ndlevel <- rownames(resultMatrix)[
																		which(resultMatrix[,"level2Down"]==1)]
		}

		if (length(myInterestedGenes1stlevel) == 0 && 
				length(myInterestedGenes2ndlevel) == 0) {
			print("There are no genes selected for enrichment analysis.")
			print("Please try to use a less stringent gene filter.")
		}
		else {
			background <- unique((toTable(org.Hs.egSYMBOL)[,"symbol"]))
			label.level.DEGs <- c(substr(a[1], 10, nchar(a[1])-1),
														substr(a[2], 10, nchar(a[2])-1))

			# Choice of the ontology
			
			if (ontology !=  "all") {
				onto.first <- createspecifictable(background, 
											myInterestedGenes1stlevel, ontology, 
											test.method, test.threshold, mult.cor, 
											label.level.DEGs[1])
											
				onto.second <- createspecifictable(background,
											 myInterestedGenes2ndlevel, ontology, 
											 test.method, test.threshold, mult.cor, 
											 label.level.DEGs[2])
											 
				final.table <- rbind(onto.first, onto.second)
						
				allRes <- new("GOsets", enriched.table = final.table, 
											label.level.enriched = label.level.DEGs)
			} 
			else {		
				BP.first <- createspecifictable(background, myInterestedGenes1stlevel, 
																				"BP", test.method, test.threshold, 
																				mult.cor, label.level.DEGs[1])
				MF.first <- createspecifictable(background, myInterestedGenes1stlevel, 
																				"MF", test.method, test.threshold,
																				mult.cor, label.level.DEGs[1])
				CC.first <- createspecifictable(background, myInterestedGenes1stlevel, 
																				"CC", test.method, test.threshold,
																				mult.cor, label.level.DEGs[1])
				
				BP.second <- createspecifictable(background, myInterestedGenes2ndlevel, 
																				"BP", test.method, test.threshold,
																				mult.cor, label.level.DEGs[2])
				MF.second <- createspecifictable(background, myInterestedGenes2ndlevel, 
																				"MF", test.method, test.threshold,
																				mult.cor, label.level.DEGs[2])
				CC.second <- createspecifictable(background, myInterestedGenes2ndlevel,
																				"CC", test.method, test.threshold,
																				mult.cor, label.level.DEGs[2])

				final.table <- rbind(BP.first, MF.first, CC.first, 
														 BP.second, MF.second, CC.second)
				allRes <- new("GOsets", enriched.table=final.table, 
											label.level.enriched=label.level.DEGs)
			}
		}
		return(allRes)
	}
)


## Implementation of the RegulatoryEnrichment method
setMethod("RegulatoryEnrichment", "DEGs", 
          function(object, classOfDEGs="both", 
									 significance.threshold= 0.05, mult.cor=TRUE) {
  
	resultMatrix <- DEGs.table(object)
	label.fc <- grep(c("FC"), colnames(resultMatrix), value=TRUE)
	label.level.DEGs <- c(substr(label.fc[1], 10, nchar(label.fc[1])-1),
												substr(label.fc[2], 10, nchar(label.fc[2])-1))

	
	genes.1stlevel <- rownames(resultMatrix)[
															which(resultMatrix[,"level1Up"]==1 | 
																		resultMatrix[,"level1Down"]==1)]
	genes.2ndlevel <- rownames(resultMatrix)[
															which(resultMatrix[,"level2Up"]==1 | 
																		resultMatrix[,"level2Down"]==1)]
		
	if (classOfDEGs == "up") {
		genes.1stlevel <- rownames(resultMatrix)[
															which(resultMatrix[,"level1Up"]==1)]
		genes.2ndlevel <- rownames(resultMatrix)[
															which(resultMatrix[,"level2Up"]==1)]
	}
			
	if (classOfDEGs == "down") {
		genes.1stlevel <- rownames(resultMatrix)[
															which(resultMatrix[,"level1Down"]==1)]
		genes.2ndlevel <- rownames(resultMatrix)[
															which(resultMatrix[,"level2Down"]==1)]
	}

	if (length(genes.1stlevel) == 0 && length(genes.2ndlevel) == 0) {
		print("There are no genes selected for enrichment analysis.")
		stop("Please try to use a less stringent gene filter.")
	}
	else {
		enriched.1 <- computeGeneListEnrichment(genes.1stlevel,
											label.level.DEGs[1], significance.threshold, mult.cor)
		enriched.2 <- computeGeneListEnrichment(genes.2ndlevel,
											label.level.DEGs[2], significance.threshold, mult.cor)
		
		return(new("EnrichedSets", enriched.table=rbind(enriched.1, enriched.2), 
												 label.level.enriched=label.level.DEGs))	
	}
}
)


## Helper functions
###############################################################################

createspecifictable <- function(background,myInterestedGenes,ontology,
	test.method,test.threshold,mult.cor,level.label) {
	
	xx <- annFUN.org(ontology, feasibleGenes=background, 
									mapping="org.Hs.eg.db", ID="symbol")
	xxxx <- inverseList(xx)
	geneNames <- names(xxxx)
	geneList <- factor(as.integer(geneNames %in% myInterestedGenes))
	names(geneList) <- geneNames
	str(geneList)
	
	GOdata <- new("topGOdata", ontology = ontology, allGenes = geneList,
								annot = annFUN.gene2GO, gene2GO = xxxx, nodeSize= 5)
		
	if (test.method == "classic") 
		resultFisher <- runTest(GOdata, algorithm = "classic", 
														statistic = "fisher")
	if (test.method == "elim") 
		resultFisher <- runTest(GOdata, algorithm = "elim", 
														statistic = "fisher")
	if (test.method == "weight") 
		resultFisher <- runTest(GOdata, algorithm = "weight", 
														statistic = "fisher")
	if (test.method == "weight01") 
		resultFisher <- runTest(GOdata, algorithm = "weight01",
														statistic = "fisher")
	if (test.method == "parentchild") 
		resultFisher <- runTest(GOdata, algorithm = "parentchild", 
														statistic = "fisher")

	s <- score(resultFisher)
	
	enrich.table <- GenTable(GOdata, Fisher.result = resultFisher, 
													 orderBy = "Fisher.result", 
													 ranksOf = "Fisher.result", topNodes = length(s))
	enrich.table$Fisher.classic.BH <- p.adjust(
																		as.numeric(enrich.table[,"Fisher.result"]),
																		"BH", 
																		n=length(enrich.table[,"Fisher.result"]))
	
	if (!mult.cor) {
		enrich.table <- enrich.table[
											which(
												enrich.table[,"Fisher.result"] <= test.threshold),]
	}
	else {
		enrich.table <- enrich.table[
											which(
												enrich.table[,"Fisher.classic.BH"] <= test.threshold),]
	}
		
	enrich.table <- enrich.table[order(enrich.table[,"Fisher.result"], 
															 na.last = TRUE, decreasing = FALSE),]
	
	OntologyLevel <- matrix(,nrow=dim(enrich.table)[1], ncol=2)
	colnames(OntologyLevel) <- c("ontology", "level")
	OntologyLevel[,1] <- ontology
	OntologyLevel[,2] <- level.label
	final.table <-  cbind (OntologyLevel, enrich.table)
	colnames(final.table) <- c("ontology", "level","GO.ID", "term", "annotated",
														 "significant", "expected",
														 "pv.fisher", "pv.fisher.BH")
	
	return (final.table)
}


computeGeneListEnrichment <- 
		function(DEGs.genes, level.label, significance.threshold, mult.cor) {
	
	# explicitly define the two tables (contained in tRanslatomeDataset)
	# to avoid Rcmd check complaining as it does not see these	
	regulatory.elements.regulated <- NULL
	regulatory.elements.counts <- NULL
	
	# we need data contained in our tRanslatome dataset, so load it
	data(tRanslatomeSampleData, envir=environment())
	
		enrichments <- c()
		
		for (i in c(1:nrow(regulatory.elements.regulated))) {
			regulated.list <- strsplit(
								as.character(regulatory.elements.regulated[i,2]), 
								split=",")[[1]]
			
			regel.name <- as.character(regulatory.elements.regulated[i,1])
			regelIdx <- which(regulatory.elements.counts[,1] == regel.name)
			regel.regulated <- as.integer(regulatory.elements.counts[regelIdx,2])
			regel.nonregulated <- as.integer(regulatory.elements.counts[regelIdx,3])
			
			# get genes from input list regulated by this regulator
			regel.input.regulated <- DEGs.genes[
																		which(DEGs.genes %in% regulated.list)] 
			
			# build the contingency table for fisher test
			cont.table <- matrix(nrow=2, ncol=2)
			# number of input genes regulated by this element
			cont.table[1,1] <- length(regel.input.regulated)
			# number of input genes NOT regulated by this element
			cont.table[1,2] <- length(DEGs.genes) - cont.table[1,1]
			# number of non-input genes regulated by this element
			cont.table[2,1] <- as.integer(regel.regulated) - cont.table[1,1]
			# number of non-input genes NOT regulated by this element
			cont.table[2,2] <- (regel.regulated + regel.nonregulated) - 
												 cont.table[2,1]
			
			# compute the fisher test p-value and report the regulator only if
			# the enrichment p-value is below the given significance threshold
			pvalue <- fisher.test(cont.table)$p.value
			if (pvalue <= significance.threshold) 
				enrichments <- rbind(enrichments, c(regel.name, level.label, 
														 cont.table[1,1],
														 paste(regel.input.regulated, collapse=","),
														 pvalue))
		} 
	
	colnames(enrichments) <- c("ID", "level", "number","list","pv.fisher")															 	
	# if requested, adjust the p-value for multiple testing and add
	# the corrected p-value to the enrichments table
	if (mult.cor)
		enrichments <- cbind(enrichments, 
												 "pv.fisher.BH"=p.adjust(enrichments[,5], 
												 method="BH",n=nrow(regulatory.elements.counts)))
												
	return(as.data.frame(enrichments[order(as.numeric(enrichments[,5])),],stringsAsFactors=F))
}
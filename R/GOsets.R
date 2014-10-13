## Class declaration
###############################################################################

setClass("GOsets", contains="EnrichedSets")


## Generics declaration (getters, setters and methods)
###############################################################################

setGeneric("GOComparison", signature="object", 
          function(object) standardGeneric("GOComparison"))


## Implementation of getters
###############################################################################


## Methods implementation
###############################################################################

setMethod("Radar", "GOsets", 
	function(object, outputformat="on screen", 
           n.nodes.1stlevel="5", n.nodes.2ndlevel="5", mult.cor=TRUE,
           ontology="MF") {
	
		pvalList1stlevel <- NULL
		pvalList2ndlevel <- NULL

		if (outputformat == "pdf") pdf(file="Radar.pdf")
		if (outputformat == "postscript") postscript(file="Radar.ps")
		if (outputformat == "jpeg") 
      jpeg(filename="Radar.jpeg", width = 650, height = 500, units = "px", 
					 pointsize = 12, bg = "transparent")
		if (outputformat == "png") 
			png(filename="Radar.png", width = 650, height = 500, units = "px", 
					pointsize = 12, bg = "transparent")

		if (length(unique(enriched.table(object)[,1])) == 1 ) {
			Atable <- enriched.table(object)[,c("level", "GO.ID", 
																					"term", "pv.fisher.BH")]
			if (!mult.cor) 
				Atable <- enriched.table(object)[,c("level", "GO.ID", 
																					"term", "pv.fisher") ] 	 
		}
		else {
			# take the rows relevant only to the selected ontology
			index <- apply(enriched.table(object), 1, 
										function(row) all(ontology %in% row))
			Atable <- enriched.table(object)[index,
																c("level", "GO.ID", "term", "pv.fisher.BH")]
			if (!mult.cor) 
				Atable <- enriched.table(object)[index,
																c("level", "GO.ID", "term", "pv.fisher")] 
		}

		indexfirst <- apply(Atable, 1, 
									function(row) all(label.level.enriched(object)[1] %in% row))
		indexsecond <- apply(Atable, 1, 
									function(row) all(label.level.enriched(object)[2] %in% row))
		List1stlevel <- Atable[indexfirst, -c(1) ] 
		List2ndlevel <- Atable[indexsecond, -c(1)]

		uniquelist <- unique(c(List1stlevel[1:n.nodes.1stlevel,"GO.ID"], 
												   List2ndlevel[1:n.nodes.2ndlevel,"GO.ID"]))
		uniquelistTerms <- unique(c(List1stlevel[1:n.nodes.1stlevel,"term"], 
																List2ndlevel[1:n.nodes.2ndlevel,"term"]))

		uniquelist <- uniquelist[!is.na(uniquelist)]
		uniquelistTerms <- uniquelistTerms[!is.na(uniquelistTerms)]

		if (length(List1stlevel[,1]) == 0 | length(List2ndlevel[,1]) == 0) {
			print("One of the lists is empty!")
		}
		else {
			for(i in 1:length(uniquelist)) {
				print(uniquelist[i])
				print(uniquelistTerms[i])
				
				if ((uniquelist[i] %in% List1stlevel[,"GO.ID"]) == TRUE) {
					pvalList1stlevel <- c(pvalList1stlevel, List1stlevel[
																which(List1stlevel[,"GO.ID"] == uniquelist[i]),
																3])
				}
				else {
					pvalList1stlevel <- c(pvalList1stlevel, 1)
				}
				
				if ((uniquelist[i] %in% List2ndlevel[,"GO.ID"]) == TRUE) {
					pvalList2ndlevel  <- c(pvalList2ndlevel, List2ndlevel[
																 which(List2ndlevel[,"GO.ID"]==uniquelist[i]), 
																 3])
				}
				else {
					pvalList2ndlevel <- c(pvalList2ndlevel, 1)
				}
			}
		}

		matrix <- matrix(c(-log(as.numeric(pvalList1stlevel), base=10), 
											 -log(as.numeric(pvalList2ndlevel), base=10)), 
										 nrow=2, byrow=TRUE)
		print(matrix)

		if ((sum(matrix == 0)) == (nrow(matrix) * ncol(matrix))) {
			print("There are no significant enriched terms in both levels.")
			print("Try to set mult.cor to FALSE.")
		}
		else {
			par(fig=c(0, 1, 0, 1), mar=c(8, 8, 8, 8), mgp=c(2, 0.75, 0))
			radial.plot(matrix, labels=uniquelistTerms, rp.type="p", 
									grid.unit="-log10 p-value", line.col=c("steelblue3","gold1"),
									show.grid.labels=3, lwd=3, start=1,
									point.symbols=18, point.col="black")
									
			legend("topleft", c(label.level.enriched(object)[1], 
													label.level.enriched(object)[2]), 
						 fill=c("steelblue3", "gold1"), bty="n")
		}
		
		# if a file device is open, close it.
		if (!(outputformat == "on screen")) dev.off()
	}
)


setMethod("Heatmap", "GOsets",
	function(object, outputformat="on screen", 
					 n.nodes.1stlevel="5", n.nodes.2ndlevel="5", mult.cor=TRUE, 
					 ontology="MF") {
	
		pvalList1stlevel <- NULL
		pvalList2ndlevel <- NULL

		if (outputformat == "pdf") pdf(file="GO.Heatmap.pdf")
		if (outputformat =="postscript") postscript(file="GO.Heatmap.ps")
		if (outputformat == "jpeg") 
			jpeg(filename="GO.Heatmap.jpeg", width = 650, height = 500, units = "px",
					 pointsize = 12, bg = "transparent")
		if (outputformat == "png") 
			png(filename="GO.Heatmap.png", width = 650, height = 500, units = "px", 
					pointsize = 12, bg = "transparent")

		if (length(unique(enriched.table(object)[,1])) == 1 ) {
			Atable <- enriched.table(object)[,c("level", "GO.ID", 
																					"term", "pv.fisher.BH")]
			if (!mult.cor) 
				Atable <- enriched.table(object)[,c("level", "GO.ID", 
																						"term", "pv.fisher") ] 			 
		}
		else {
			# take the rows relevant only to the selected ontology
			index <- apply(enriched.table(object), 1, 
										 function(row) all(ontology %in% row))
			Atable <- enriched.table(object)[index, 
																 c("level","GO.ID","term","pv.fisher.BH")]
			if (!mult.cor) 
				Atable <- enriched.table(object)[index,
																	 c("level","GO.ID","term","pv.fisher")]  
		}

		indexfirst <- apply(Atable, 1, 
									function(row) all(label.level.enriched(object)[1] %in% row))
		indexsecond <- apply(Atable, 1, 
									function(row) all(label.level.enriched(object)[2] %in% row))
		List1stlevel <- Atable[indexfirst, -c(1) ] 
		List2ndlevel <- Atable[indexsecond, -c(1)]

		uniquelist <- unique(c(List1stlevel[1:n.nodes.1stlevel, "GO.ID"], 
													 List2ndlevel[1:n.nodes.2ndlevel, "GO.ID"]))
		uniquelistTerms <- unique(c(List1stlevel[1:n.nodes.1stlevel, "term"], 
																List2ndlevel[1:n.nodes.2ndlevel, "term"]))

		uniquelist <- uniquelist[!is.na(uniquelist)]
		uniquelistTerms <- uniquelistTerms[!is.na(uniquelistTerms)]

		if (length(List1stlevel[,1]) == 0 | length(List2ndlevel[,1]) == 0) {
			print("One of the lists is empty!")
		}
		else {
			for(i in 1:length(uniquelist)) {

				print(uniquelist[i])
				print(uniquelistTerms[i])
				
				if ((uniquelist[i] %in% List1stlevel[,"GO.ID"]) == TRUE) {
					pvalList1stlevel  <-  c(pvalList1stlevel, List1stlevel[
																 which(List1stlevel[,"GO.ID"]==uniquelist[i]),
																 3])
				}
				else {
					pvalList1stlevel <- c(pvalList1stlevel, 1)
				}
					
				if ((uniquelist[i] %in% List2ndlevel[,"GO.ID"]) == TRUE) {
					pvalList2ndlevel  <- c(pvalList2ndlevel, List2ndlevel[
																 which(List2ndlevel[,"GO.ID"]==uniquelist[i]),
																 3])
				}
				else {
					pvalList2ndlevel <- c(pvalList2ndlevel, 1)
				}
			}
		}

		matrix <- matrix(c(-log(as.numeric(pvalList1stlevel), base=10), 
											 -log(as.numeric(pvalList2ndlevel), base=10)), 
										 ncol=2)
		print(matrix)

		print(dim(matrix))
		print(dim(uniquelistTerms))

		if ((sum(matrix == 0)) == (nrow(matrix) * ncol(matrix))) {
			print("There are no significant enriched terms in both levels.")
			print("Try to set mult.cor to FALSE.")
		}
		else {
			rownames(matrix) <- uniquelistTerms
			colnames(matrix) <- c(label.level.enriched(object)[1], 
														label.level.enriched(object)[2])
			print(matrix)
			heatmap.2(matrix, col = RGBColVec(50)[c(26:50)], key=TRUE, 
								margins = c(7,20), keysize=1.5, denscol="white", na.rm = TRUE,
							  scale="none", dendrogram="both", trace="none", 
							  Colv=TRUE, Rowv=TRUE, labRow=NULL, labCol=NULL, 
							  cexRow=1, cexCol=1)
		}
		
		# if a file device is open, close it.
		if (!(outputformat == "on screen")) dev.off()
	}
)


setMethod("GOComparison", "GOsets",
	function(object) {
		resulttable <- NULL

		# identity - stays the same for both of the cases
		# first find the duplicates, then when encountering that duplicate 
		# check wether it is already covered, 
		# 	if not, leave that row as it is, put the label "both levels" 
		# and put it in "coveredDuplicates", 
		# 	if it is already covered mark it for deletion
		
		Atable <- enriched.table(object)[,c("ontology", "level", "GO.ID", "term")]
		duplicates <- Atable[duplicated(Atable[,"GO.ID"]), "GO.ID"]
		Atable$level <- as.character(Atable$level)
		coveredDuplicates <- c()
		rowsForDeletion <- c()

		for (i in 1:nrow(Atable)) {
			if (any(duplicates == Atable[i, "GO.ID"])) {
				if (all(coveredDuplicates != Atable[i, "GO.ID"])) {
					Atable[i,"level"] <- "both levels"
					coveredDuplicates <- c(coveredDuplicates, Atable[i,"GO.ID"])
				}
				else {
					rowsForDeletion <- c(rowsForDeletion, i)
				}
			}
			else {
				Atable[i,"level"] <- paste(Atable[i,"level"],"only")
			}
		}

		Atable <- Atable[-c(rowsForDeletion),]

		if (length(unique(enriched.table(object)[,1])) == 1 ) {
			Btable <- enriched.table(object)[,c("level","GO.ID","term")]
			ontology <- as.character(unique(enriched.table(object)[,1]))

			indexfirst <- apply(Btable, 1, 
									function(row) all(label.level.enriched(object)[1] %in% row))
			indexsecond <- apply(Btable, 1, 
									function(row) all(label.level.enriched(object)[2] %in% row))
			List1stlevel <- Btable[indexfirst,-c(1)] 
			List2ndlevel <- Btable[indexsecond,-c(1)]

			if (length(List1stlevel[,1]) == 0 | length(List2ndlevel[,1]) == 0) {
				print("One of the lists is empty!")
			}
			else {
				resulttable1to2 <- comparisonBetweenTwoLists(List1stlevel, 
													List2ndlevel, ontology, 
													paste(label.level.enriched(object)[1], " to ", 
																label.level.enriched(object)[2]))
																
				resulttable2to1 <- comparisonBetweenTwoLists(List2ndlevel, 
													List1stlevel, ontology, 
													paste(label.level.enriched(object)[2], " to ",
																label.level.enriched(object)[1]))
				
				#bind the two result tables
				Similarity.Matrix <- rbind(resulttable1to2, resulttable2to1)

				totalsimilarity <- c(mean(c(as.numeric(resulttable1to2[,7]), 
																		as.numeric(resulttable2to1[,7]))))
				names(totalsimilarity) <- c(
									as.character(unique(enriched.table(object)[,1])))
			}
		}
		else {	
			ontologyCC <- specificOntologyResult(object, "CC")
			ontologyMF <- specificOntologyResult(object, "MF")
			ontologyBP <- specificOntologyResult(object, "BP")
			
			totalsimilarity <- c(mean(c(as.numeric(ontologyCC$result1[,7]), 
																	as.numeric(ontologyCC$result2[,7]))),
                           mean(c(as.numeric(ontologyMF$result1[,7]), 
																	as.numeric(ontologyMF$result2[,7]))),
                           mean(c(as.numeric(ontologyBP$result1[,7]), 
																	as.numeric(ontologyBP$result2[,7]))))
			names(totalsimilarity) <- c("CC","MF","BP")

			Similarity.Matrix <- rbind(ontologyCC$result1, ontologyCC$result2, 
																 ontologyMF$result1, ontologyMF$result2, 
																 ontologyBP$result1, ontologyBP$result2)
		}
		
		# return the obtained results table
		return(new("GOsims", 
							similarity.matrix = Similarity.Matrix, 
							identity.matrix = as.matrix(Atable), 
							average.similarity.scores = totalsimilarity))
	}
)


## Helper functions
###############################################################################

specificOntologyResult <- function(object, ontology) {

	index <- apply(enriched.table(object), 1, 
								 function(row) all(ontology %in% row))
	ontologyTable <- enriched.table(object)[index,c("level", "GO.ID", "term") ] 
	
	indexfirst <- apply(ontologyTable, 1, 
								function(row) all(label.level.enriched(object)[1] %in% row))
	indexsecond <- apply(ontologyTable, 1, 
								function(row) all(label.level.enriched(object)[2] %in% row))
	List1stlevel <- ontologyTable[indexfirst, -c(1)] 
	List2ndlevel <- ontologyTable[indexsecond, -c(1)]

	if (length(List1stlevel[,1]) == 0 | length(List2ndlevel[,1]) == 0) {
		print("There aren't enough BP GO enriched terms in one of the levels!")
	}
	else {
		resulttable1to2 <- comparisonBetweenTwoLists(List1stlevel, List2ndlevel, 
											 ontology, paste(label.level.enriched(object)[1], " to ",
																			label.level.enriched(object)[2]))
		resulttable2to1 <- comparisonBetweenTwoLists(List2ndlevel, List1stlevel, 
											 ontology, paste(label.level.enriched(object)[2], " to ", 
																			 label.level.enriched(object)[1]))
		
		return(list(result1 = resulttable1to2, result2 = resulttable1to2))
	}
}

comparisonBetweenTwoLists <- function(go1,go2,ontology,direction) {

	finalmat <- NULL
	for (i in 1:nrow(go1)) {
		mat <- NULL
		for (j in 1:nrow(go2)) {
			sim <- goSim(go1[i, 1], go2[j, 1], 
									 organism="human", measure="Wang", ont = ontology)
			mat <- rbind(mat, c(go1[i, 1],go1[i, 2],go2[j, 1],go2[j, 2],sim))	
		}
		
		max1.index <- which.max(mat[,5])
		max1 <- mat[max1.index, 5]	
		
		firstgotermID <- mat[max1.index, 3]
		firstgotermTERM <- mat[max1.index, 4]	
		
		finalmat <- rbind(finalmat, c(go1[i, 1], go1[i, 2], 
																	firstgotermID, firstgotermTERM, max1))
		colnames(finalmat) <- c("start.GO.ID", "start.term", 
														"end.GO.ID", "end.term", 
														"similarity.score")			
	}
	
	OntologyAndDirection <- matrix(,nrow=dim(finalmat)[1], ncol=2)
	colnames(OntologyAndDirection) <- c("ontology", "direction")
	OntologyAndDirection[,1] <- ontology
	OntologyAndDirection[,2] <- direction
	
	return(cbind(OntologyAndDirection, finalmat))
}

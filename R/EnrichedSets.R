## Class declaration
###############################################################################

setClass("EnrichedSets", representation(enriched.table="data.frame", 
																				label.level.enriched="character"))


## Generics declaration (getters, setters and methods)
###############################################################################

setGeneric("enriched.table", signature="object", 
          function(object) standardGeneric("enriched.table"))
          
setGeneric("label.level.enriched", signature="object", 
          function(object) standardGeneric("label.level.enriched"))

setGeneric("Radar", signature="object", 
          function(object, outputformat="on screen",
                   n.nodes.1stlevel="5", n.nodes.2ndlevel="5", mult.cor=TRUE, ...) 
          standardGeneric("Radar"))
          
setGeneric("Heatmap", signature="object", 
          function(object, outputformat="on screen", 
                   n.nodes.1stlevel="5", n.nodes.2ndlevel="5", mult.cor=TRUE, ...)
          standardGeneric("Heatmap"))
 

## Implementation of getters
###############################################################################

setMethod("enriched.table", "EnrichedSets",
	function(object) {
		object@enriched.table
	}
)

setMethod("label.level.enriched", "EnrichedSets",
	function(object) {
		object@label.level.enriched
	}
)


## Methods implementation
###############################################################################

setMethod("show", "EnrichedSets",
	function(object) {
		print(apply(object@enriched.table, c(1,2), function(x){
							if(nchar(x) > 30) 
								return(paste(substr(x,1,30),"...",sep=""))
							else
								return(x)
							})
				 )
		print(object@label.level.enriched)
	}
)


setMethod("Radar", "EnrichedSets", 
	function(object, outputformat="on screen", 
           n.nodes.1stlevel="5", n.nodes.2ndlevel="5", mult.cor=TRUE) {
	
		pvalList1stlevel <- NULL
		pvalList2ndlevel <- NULL

		if (outputformat == "pdf") pdf(file="Enrichment.Radar.pdf")
		if (outputformat == "postscript") postscript(file="Enrichment.Radar.ps")
		if (outputformat == "jpeg") 
      jpeg(filename="Enrichment.Radar.jpeg", width = 650, height = 500, 
					 units = "px", pointsize = 12, bg = "transparent")
		if (outputformat == "png") 
			png(filename="Enrichment.Radar.png", width = 650, height = 500, 
					units = "px", pointsize = 12, bg = "transparent")

		Atable <- enriched.table(object)[,c("level","ID","list","pv.fisher") ]
		if (mult.cor) 
			Atable <- enriched.table(object)[,c("level","ID","list", "pv.fisher.BH") ] 	 

		indexfirst <- apply(Atable, 1, 
												function(row) all(label.level.enriched(object)[1] %in% row))
		indexsecond <- apply(Atable, 1, 
												function(row) all(label.level.enriched(object)[2] %in% row))
		List1stlevel <- Atable[indexfirst, -c(1)] 
		List2ndlevel <- Atable[indexsecond, -c(1)]

		uniquelist <- unique(c(List1stlevel[1:n.nodes.1stlevel,"ID"], 
												   List2ndlevel[1:n.nodes.2ndlevel,"ID"]))
		uniquelist <- uniquelist[!is.na(uniquelist)]

		if (length(List1stlevel[,1]) == 0 | length(List2ndlevel[,1]) == 0) {
			print("One of the lists is empty!")
		}
		else {
			for(i in 1:length(uniquelist)) {
				print(uniquelist[i])
				
				if ((uniquelist[i] %in% List1stlevel[,"ID"]) == TRUE) {
					pvalList1stlevel <- c(pvalList1stlevel, List1stlevel[
																which(List1stlevel[,"ID"] == uniquelist[i]),
																3])
				}
				else {
					pvalList1stlevel <- c(pvalList1stlevel, 1)
				}
				
				if ((uniquelist[i] %in% List2ndlevel[,"ID"]) == TRUE) {
					pvalList2ndlevel  <- c(pvalList2ndlevel, List2ndlevel[
																 which(List2ndlevel[,"ID"]==uniquelist[i]), 
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
			print("There are no significant enriched regulators in both levels.")
			print("Try to set mult.cor to FALSE.")
		}
		else {
			par(fig=c(0, 1, 0, 1), mar=c(8, 8, 8, 8), mgp=c(2, 0.75, 0))
			radial.plot(matrix, labels=uniquelist, rp.type="p", 
									grid.unit="-log10 p-value", line.col=c("steelblue3","gold1"),
									show.grid.labels=3, lwd=3, start=1,
									point.symbols=18, point.col="black")
									
			legend("topleft", c(label.level.enriched(object)[1], label.level.enriched(object)[2]), 
						 fill=c("steelblue3", "gold1"), bty="n")
		}
		
		# if a file device is open, close it.
		if (!(outputformat == "on screen")) dev.off()
	}
)


setMethod("Heatmap", "EnrichedSets",
	function(object, outputformat="on screen", 
					 n.nodes.1stlevel="5", n.nodes.2ndlevel="5", mult.cor=TRUE) {
	
		pvalList1stlevel <- NULL
		pvalList2ndlevel <- NULL

		if (outputformat == "pdf") pdf(file="Enrichment.Heatmap.pdf")
		if (outputformat =="postscript") postscript(file="Enrichment.Heatmap.ps")
		if (outputformat == "jpeg") 
			jpeg(filename="Enrichment.Heatmap.jpeg", width = 650, height = 500, 
					 units = "px", pointsize = 12, bg = "transparent")
		if (outputformat == "png") 
			png(filename="Enrichment.Heatmap.png", width = 650, height = 500, 
					units = "px", pointsize = 12, bg = "transparent")

		Atable <- enriched.table(object)[,c("level", "ID", "list", "pv.fisher")]
		if (mult.cor) 
			Atable <- enriched.table(object)[,c("level", "ID", "list", "pv.fisher.BH")] 

		indexfirst <- apply(Atable, 1, 
												function(row) all(label.level.enriched(object)[1] %in% row))
		indexsecond <- apply(Atable, 1, 
												function(row) all(label.level.enriched(object)[2] %in% row))
		List1stlevel <- Atable[indexfirst, -c(1) ] 
		List2ndlevel <- Atable[indexsecond, -c(1)]

		uniquelist <- unique(c(List1stlevel[1:n.nodes.1stlevel, "ID"], 
													 List2ndlevel[1:n.nodes.2ndlevel, "ID"]))
		uniquelist <- uniquelist[!is.na(uniquelist)]
		
		if (length(List1stlevel[,1]) == 0 | length(List2ndlevel[,1]) == 0) {
			print("One of the lists is empty!")
		}
		else {
			for(i in 1:length(uniquelist)) {

				print(uniquelist[i])
				
				if ((uniquelist[i] %in% List1stlevel[,"ID"]) == TRUE) {
					pvalList1stlevel  <-  c(pvalList1stlevel, List1stlevel[
																 which(List1stlevel[,"ID"]==uniquelist[i]),
																 3])
				}
				else {
					pvalList1stlevel <- c(pvalList1stlevel, 1)
				}
					
				if ((uniquelist[i] %in% List2ndlevel[,"ID"]) == TRUE) {
					pvalList2ndlevel  <- c(pvalList2ndlevel, List2ndlevel[
																 which(List2ndlevel[,"ID"]==uniquelist[i]),
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
		print(length(uniquelist))

		if ((sum(matrix == 0)) == (nrow(matrix) * ncol(matrix))) {
			print("There are no significant enriched regulators in both levels.")
			print("Try to set mult.cor to FALSE.")
		}
		else {
			rownames(matrix) <- uniquelist
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


## Helper functions
###############################################################################

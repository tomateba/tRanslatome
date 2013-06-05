## Class declaration
###############################################################################

setClass("GOsims", representation(similarity.matrix="matrix", 
																	identity.matrix="matrix", 
																	average.similarity.scores = "vector"))


## Generics declaration (getters, setters and methods)
###############################################################################

setGeneric("similarity.matrix", 
	function(object) standardGeneric("similarity.matrix"))
	
setGeneric("identity.matrix", 
	function(object) standardGeneric("identity.matrix"))
	
setGeneric("average.similarity.scores", 
	function(object) standardGeneric("average.similarity.scores"))
	
setGeneric("IdentityPlot", 
	function(object, outputformat="on screen") standardGeneric("IdentityPlot"))
	
setGeneric("SimilarityPlot", 
	function(object, outputformat="on screen") standardGeneric("SimilarityPlot"))


## Implementation of getters
###############################################################################

setMethod("similarity.matrix", "GOsims",
	function(object) {
		object@similarity.matrix
	}
)

setMethod("identity.matrix", "GOsims",
	function(object) {
		object@identity.matrix
	}
)

setMethod("average.similarity.scores", "GOsims",
	function(object) {
		object@average.similarity.scores
	}
)


## Methods implementation
###############################################################################

setMethod("show", "GOsims",
	function(object) {
		print(head(identity.matrix(object)))
		print(head(similarity.matrix(object)))
		print(head(average.similarity.scores(object)))
	}
)


setMethod("IdentityPlot", "GOsims",
	function(object,outputformat="on screen") {
		
		if (length(average.similarity.scores(object)) != 3) {
			print("Similarity scores for some of the ontologies are missing!")
		}
		else {
			if (outputformat == "pdf") pdf(file="IdentityPlot.pdf")
			if (outputformat == "postscript") postscript(file="IdentityPlot.ps")
			if (outputformat == "jpeg") jpeg(filename="IdentityPlot.jpeg")
				
			barplot(c(average.similarity.scores(object)[3], 
								average.similarity.scores(object)[2], 
								average.similarity.scores(object)[1]), 
							names.arg=c("BP","MF","CC"), col=c("grey25","grey75","grey"), 
							las=1, border=NA, beside=FALSE, ylab="semantic similarity")	

			# if a file device is open, close it.
			if (!(outputformat == "on screen")) dev.off()		
		}
	}
)


setMethod("SimilarityPlot", "GOsims",
	function(object,outputformat="on screen")	{
	
		if (length(average.similarity.scores(object)) != 3) 
			print("Similarity scores for some of the ontologies may be missing!")
		
		if (outputformat == "pdf") pdf(file="SimilarityPlot.pdf")
		if (outputformat == "postscript") postscript(file="SimilarityPlot.ps")
		if (outputformat == "jpeg") jpeg(filename="SimilarityPlot.jpeg")
				
		namesOfLevels <- unique(identity.matrix(object)[,"level"])	
		namesOfLevels <- namesOfLevels[-c(namesOfLevels[]=="both levels")]
			
		#CC	
		CC <- identity.matrix(object)[identity.matrix(object) 
																		[,"ontology"] == "CC",]
		CC1 <- dim(CC[CC[,"level"] == namesOfLevels[1],])[1]
		CC2 <- dim(CC[CC[,"level"] == namesOfLevels[2],])[1]
		CC3 <- dim(CC[CC[,"level"] == "both levels",])[1]
			
		#MF
		MF <- identity.matrix(object)[identity.matrix(object) 
																		[,"ontology"] == "MF",]
		MF1 <- dim(MF[MF[,"level"] == namesOfLevels[1],])[1]
		MF2 <- dim(MF[MF[,"level"] == namesOfLevels[2],])[1]
		MF3 <- dim(MF[MF[,"level"] == "both levels",])[1]
			
		#BP
		BP <- identity.matrix(object)[identity.matrix(object) 
																		[,"ontology"] == "BP",]
		BP1 <- dim(BP[BP[,"level"] == namesOfLevels[1],])[1]
		BP2 <- dim(BP[BP[,"level"] == namesOfLevels[2],])[1]
		BP3 <- dim(BP[BP[,"level"] == "both levels",])[1]
			
		#plot the similarity between levels and ontologies
		inputmatrix <- matrix(c(CC1, CC2, CC3, MF1, MF2, MF3, BP1, BP2, BP3), 
													nrow=3, ncol=3)		
		barplot(inputmatrix, beside=TRUE, las=1, border=NA, xlab="GO ontology",
						ylab="number of GO enriched terms", names.arg=c("CC","MF","BP"), 
						col=c("steelblue3","chartreuse3","gold1"))	
						
		legend("topleft", c(namesOfLevels[1],"both levels",namesOfLevels[2]), 
					 bty="n", fill=c("steelblue3","chartreuse3","gold1"))
		
		# if a file device is open, close it.
		if (!(outputformat == "on screen")) dev.off()			
	}
)
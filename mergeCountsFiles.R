# Input arguments:
# 	- dirWithCountsForAnAccessNum: full path of the dir that contains the folders 
#	with the processed.count files (./1/proc.cnts, ./2/proc.cnts, ect...) 
#
#	- numOfCountsFiles: number of processed.count files (or alternatively ./1/processed.counts, 
#	./2/processed.counts, ect folders) to merge (by adding up elements from the 2nd column 
#	that have the same row name at the 1st column

## Read input
args <- commandArgs(trailingOnly = TRUE);
dirWithCountsForAnAccessNum <- args[1];
numOfCountsFiles <- args[2];
pathToSaveTheMergedCountsFile <- args[3];
curAccessNum <- args[4];

##### Functions #####
mergeBaseDfWithCountsFile <- function(baseDf, countFileToMergePath){

#df1 <- read.table(textConnection(skip_FirstAndLast1), header = FALSE, stringsAsFactors = FALSE, sep ="\t", skip=1, col.names=c("x", "y"))
        colnames(baseDf) = c("x", "y")


        allNewCountsToMergeFile <- file(countFileToMergePath)
	allNewCountsToMerge <- readLines(allNewCountsToMergeFile)
        allNewCountsToMergeExceptLastLine <- head(allNewCountsToMerge, -1)
        allNewCountsToMergeExceptLastLineTCon <- textConnection(allNewCountsToMergeExceptLastLine) 
	
	dfToMerge <- read.table(allNewCountsToMergeExceptLastLineTCon, header = FALSE, stringsAsFactors = FALSE, sep ="\t", skip=1, col.names=c("x", "z"))
	close(allNewCountsToMergeExceptLastLineTCon)	
	close(allNewCountsToMergeFile)

        dfm <- merge(baseDf, dfToMerge, all=TRUE)

        mergedDf <- data.frame(x = dfm$x, y = rowSums(dfm[, -1], na.rm = TRUE))
        return(mergedDf)
}
#####################

# Vector containing the full paths for all the produced processed.counts files for an accession number
countsFilesVector <- paste(dirWithCountsForAnAccessNum, 1:numOfCountsFiles, "processed.counts", sep="/")

# Create the base data frame from the processed.counts file of the 1/ folder
allCountsFile <- file(countsFilesVector[1])
allCounts <- readLines(allCountsFile)
allCountsExceptLastRow <- head(allCounts, -1)
allCountsExceptLastRowTCon <- textConnection(allCountsExceptLastRow)


baseDf <- read.table(allCountsExceptLastRowTCon, header = FALSE, stringsAsFactors = FALSE, sep ="\t", skip=1, col.names=c("x", "y"))

close(allCountsExceptLastRowTCon)
close(allCountsFile)

# merge base data frame with the rest of the processed.counts files
if(numOfCountsFiles > 1){
	for (fId in 2:numOfCountsFiles){
		newDfBase <- mergeBaseDfWithCountsFile(baseDf, countsFilesVector[fId])
		baseDf <- newDfBase
	}
}

colnames(baseDf) <- c("miRNA_id", "depth")

# sort by 'depth' column (descending)
baseDf <- baseDf[with(baseDf, order(-depth)), ]
row.names(baseDf) <- 1:nrow(baseDf)

# final merged processed.counts 'data frame'
#baseDf

#countsFinalFileName <- paste(cuAccessionNum, "processed.counts", sep="_")
processedCountsFname <- paste(curAccessNum, "_processed.counts", sep="")
mergedCountsFileFullName <- paste(pathToSaveTheMergedCountsFile, processedCountsFname, sep="/")

#write.table(baseDf, file=mergedCountsFileFullName, sep=",");
write.csv(baseDf, file=mergedCountsFileFullName)


library(dada2)
#library(dada2);package_version("dada2")
path <- "C:\\Users\\mclar\\Downloads\\OneDrive-2020-01-08 (1)"
list.files(path)

#Now we read in the names of the fastq files, and perform some string manipulation to get matched lists of the forward and reverse fastq files.
#Forward and Reverse fastq filename have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001_fastq
fnFs <- sort(list.files(path,pattern = '_R1.fastq', full.names = TRUE))
fnRs <- sort(list.files(path,pattern = '_R2.fastq', full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspect read quality profiles
#Now we start by visualizing the quality profiles of the forward reads:
plotQualityProfile(fnFs)
#Now we visualize the quality profile of the reverse reads:
plotQualityProfile(fnRs)

#FILTER AND TRIM
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,240),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

# LEARNING THE ERROR RATES
#every amplicon dataset has a different set of error rates. The learnErrors method learns this error model from the data, 
#by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution.

errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)
head(errF)

#visualize the estimated error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Points are the observed error rates for each consensus quality score. The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. Here the estimated error rates (black line) are a good fit to the observed rates (points), 
#and the error rates drop with increased quality as expected. Everything looks reasonable and we proceed with confidence.


#DEREPLICATE READS INTO UNIQUE SEQUENCES
#The next step is to remove  sequences that share 100% similarity,  since most of the reads in our dataset should be exactly the same 
derepFwd <- derepFastq(filtFs) # dereplicating fwd reads 
derepRev <- derepFastq(filtRs) # dereplicating rev reads 


#To see how many unreplicated sequences we had after doing the unreplication, we run the script above:
countUniqueAndTotalFromDerepObjcList <- function(deRepObj) {
  
  ### count the total no. of unique seqs in the dereplication dada2 object
  
  listOut <- list()
  
  uniq <- sapply(X = deRepObj, FUN = summary)[1,]
  uniqOut <- as.numeric(uniq)
  names(uniqOut) <- names(deRepObj)
  
  total <- sapply(X = deRepObj, FUN = summary)[3,]
  totalOut <- as.numeric(total)
  names(totalOut) <- names(deRepObj)
  
  listOut[["unique"]] <- uniqOut
  listOut[["total"]] <- totalOut
  
  return(listOut)
  
}

seqDerepFwd <- countUniqueAndTotalFromDerepObjcList(deRepObj = derepFwd)
uniqDerepFwd <- seqDerepFwd$unique
total <- seqDerepFwd$total 
seqDerepRev <- countUniqueAndTotalFromDerepObjcList(derepRev) 
uniqDerepRev <- seqDerepRev$unique 
derepMtxAbund <- rbind(total, uniqDerepFwd, uniqDerepRev)
barplot(derepMtxAbund, main = "Unique/dereplicated sequences (absolute abundance)", xlab = "Samples", ylab = "Unique sequences (absolute abundance)", col = c("green","blue", "lightblue"), legend = rownames(x = derepMtxAbund), beside=TRUE)


Unique_Fwd_Perc <- uniqDerepFwd / total * 100
Unique_Rev_Perc <- uniqDerepRev / total * 100
derepMtxPerc <- rbind(Unique_Fwd_Perc, Unique_Rev_Perc) 
barplot(derepMtxPerc, main = "Unique/dereplicated sequences (%)", xlab = "Samples", ylab = "Unique sequences (%)", col = c("gold", "magenta"), legend = rownames(x = derepMtxPerc), beside=TRUE)
#It was then observed that the number of total strings decreases approximately 80% for fowards and 50% for reverse

#DENOISE UNIQUE SEQUENCES (AND EXCLUDE SINGLETONS)
#Using the error rates obtained previously, we used the function dada to do a denoising of the amplicons, 
#that is, remove the sequences with errors

dadaFwd <- dada(derepFwd, errF, multithread = TRUE) # denoise fwd seqs

dadaRev <- dada(derepRev, errR, multithread = TRUE) # denoise rev seqs

dadaFwd

dadaRev

#MERGE DENOISED FORWARD AND REVERSE READS
#After denoising, we will merge our forward reads with the respective reverse-complement of reverse reads.
#By default, if forward and the reverse-complement of reverse reads does not overlap, reads are discarded as well forward and reverse reads that does not have their pair.
mergePE <- mergePairs(dadaF = dadaFwd, derepF = derepFwd, dadaR = dadaRev, derepR = derepRev, verbose = TRUE)

countMergedSeqsFromDadaObjcList <- function(mergedObj) {
  
  ### count the total no. of merged seqs in the dada2 object
  
  merged <- sapply((sapply(X = mergedObj, `[`, 2)), sum)
  mergedOut <- as.numeric(merged)
  names(mergedOut) <- names(mergedObj)
  
  return(mergedOut)
  
}
## minimum
min(countMergedSeqsFromDadaObjcList(mergePE))

## maximum
max(countMergedSeqsFromDadaObjcList(mergePE))

## Average of percentage of sequences lost 

mean((out[,1] - countMergedSeqsFromDadaObjcList(mergePE)) / out[,1] * 100)

#to plot it instead, you can also do the following (in this case we represented the sequences kept):
mergedPE <- countMergedSeqsFromDadaObjcList(mergePE) # get the number of merged reads
countTable <- cbind(out, mergedPE)
colnames(countTable) <- c("Initial", "Filtered", "Merged")
rownames(countTable) <- sample.names
countTable <- t(countTable)
barplot(countTable, col = c("#EFF3FF", "#BDD7E7", "#6BAED6"), main = "Comparison of no. of sequences kept across pipeline steps: initial, filtered and merged", ylab = "Absolute no. of sequences", xlab = "Samples", legend = rownames(countTable), beside = TRUE)

# Doing the same plot but in percentage 
countTablePerc <- apply(X = countTable, MARGIN = 1, FUN = function(x) x / countTable[1,] * 100) 
countTablePerc <- t(countTablePerc)
barplot(countTablePerc, col = c("#EFF3FF", "#BDD7E7", "#6BAED6"), main = "Comparison of proportion of sequences kept across pipeline steps: initial, filtered and merged", ylab = "Percentage of sequences (%)", xlab = "Samples", legend = rownames(countTablePerc), beside = TRUE)
abline(h = min(countTablePerc[3,]))

#Then, proceed to the construction of an ASV table using the makeSequenceTable () function to tabulate the ASVs using the assembly sequence and the getSequences () function to organize the table

### Make a ASV table

asvTbl <- makeSequenceTable(samples = mergePE) # tabulate ASVs

histSeqLen <- table(nchar(getSequences(asvTbl)))

dim(asvTbl)

#REMOVE CHIMERAS
# Remove chimeras from the ASV table
asvTblNoChim <- removeBimeraDenovo(unqs = asvTbl, method = "consensus", multithread = TRUE, verbose = TRUE) 

#Summarize the no. of sequences kept in each pipeline step
getN <- function(x) sum(getUniques(x)) # function that sums `sum(getUniques(x)` the no. of unique sequences `getUniques(x)`

## build a matrix with all the sequences kept in each pipeline step
summaryTblSeq <- cbind(out, # initial reads and filtered/trimmed reads
                       sapply(dadaFwd, getN), sapply(dadaRev, getN), # denoised sequences 
                       sapply(mergePE, getN), # merged PE sequences
                       rowSums(asvTblNoChim)) # non-chimeric sequences

## rename the column and row names 
colnames(summaryTblSeq) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(summaryTblSeq) <- sample.names

## create a second summary table seq with one column for the samples 
summaryTblSeq2 <- cbind("Samples" = sample.names, summaryTblSeq)

dir.create("output2") # let's create 'output' folder
write.table(x = summaryTblSeq2, file = "./output/summaryTblSeq.tsv", sep = "\t", row.names = FALSE)

knitr::kable(summaryTblSeq)

# Barplot with the abs. abundance of sequences
summaryTblSeqTrans <- t(summaryTblSeq) # transpose 
barplot(summaryTblSeqTrans, main = "Absolute no. of sequences kept through the pipeline", ylab = "Absolute no. of sequences", xlab = "Samples", col = c("gray", "#EFF3FF", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C"), legend = rownames(summaryTblSeqTrans), beside = TRUE)

# Barplot in percentage
summaryTblSeqPerc <- apply(X = summaryTblSeq, MARGIN = 2, function(x) x / summaryTblSeq[,1] * 100) # get the correspondent percentage table
summaryTblSeqPercTrans <- t(summaryTblSeqPerc) # transpose
barplot(summaryTblSeqPercTrans, main = "Percentage of sequences kept through the pipeline", ylab = "Percentage of sequences (%)", xlab = "Samples", col = c("gray", "#EFF3FF", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C"), legend = rownames(summaryTblSeqPercTrans), beside = TRUE) # plot it 

#ASSINING OUR TAXONOMY
#In this step we will use the SILVA training reference database to assign a taxonomy to our samples.
#The first file assigns a taxonomy to the table and the second adds the species name if there is a match in our ASVs
taxTbl <- assignTaxonomy(seqs = asvTblNoChim, refFasta = "C:\\Users\\mclar\\Downloads\\silva_nr_v132_train_set.fa.gz", multithread = FALSE)
taxTbl <- addSpecies(taxtab = taxTbl, refFasta = "C:\\Users\\mclar\\Downloads\\silva_species_assignment_v132.fa.gz")
taxTbl.print <- taxTbl 
rownames(taxTbl.print) <- NULL
head(taxTbl.print)


# Save ASV and taxonomy tables in R format 
saveRDS(object = asvTblNoChim, file = "./output/asvTblNoChim.rds") # save the ASV table
saveRDS(object = taxTbl, file = "./output/taxTbl.rds") # save the ASV taxonomy

tax2biom <- function(taxTable) {
  
  lenNrow = nrow(taxTable)
  lenNcol = ncol(taxTable)
  id = c()
  fullTax = c()
  
  
  for (irow in 1:lenNrow) {
    rowID = as.character(rownames(taxTable)[irow])
    id = append(id, rowID)
    tax = ""
    
    for (icol in 1:lenNcol) {
      colName <- colnames(taxTable)[icol]
      colSign <- paste0(tolower(strtrim(colName,1)), "__")
      
      if (tax == "") {
        tax <- paste0(colSign, taxTable[irow,icol])
        #tax <- taxTable[irow,icol]
      }
      else {
        tax <- paste(tax, paste0(colSign, taxTable[irow,icol]), sep = "; ")
        #tax <- paste(tax, taxTable[irow,icol], sep = "; ")
      }
    }
    fullTax <- append(fullTax, tax)
  }
  
  taxTable2Biom <- data.frame(id, fullTax)
  colnames(taxTable2Biom) <- c("ASV_ID", "taxonomy")
  
  return(taxTable2Biom)
}
taxTbl2 <- cbind(taxTbl, "ASV" = paste0("ASV_", 1:nrow(taxTbl)))
rownames(taxTbl2) <- taxTbl2[,8]
asvTblNoChim2 <- asvTblNoChim 
colnames(asvTblNoChim2) <- taxTbl2[,8] 
asvTblNoChim2 <- t(asvTblNoChim2) 
asvTblNoChim2 <- as.data.frame(asvTblNoChim2)
asvTblNoChim2[,"ASV_ID"] <- rownames(asvTblNoChim2)
write.table(x = taxTbl2, file = "./output/taxTbl.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(x = asvTblNoChim2, file = "./output/asvTblNoChim.txt", sep = "\t", row.names = FALSE, quote = FALSE)


taxTbl2 <- tax2biom(taxTbl2) 
asvTaxTbl <- cbind(asvTblNoChim2, "taxonomy" = taxTbl2[,-1]) 
write.table(x = asvTaxTbl, file = "./output/asvTaxTbl.txt", sep = "\t", row.names = FALSE, quote = FALSE) 
summary(asvTaxTbl)

source("C:\\Users\\mclar\\Downloads\\Rscript.R") # import R script with built-in functions
taxTbl2 <- tax2biom(taxTbl2) 

BiocManager::install('biomformat')
#install.packages('command -v biom')
library('biom')
convertTab2Biom <- function(inFile, outFile) {
  
  if (system("command - v biom", ignore.stdout = TRUE, ignore.stderr = TRUE) !=0)  {
    
    stop("biom program is not installed or it is not accessible!\n  Exiting...")
    
  }
  
  system(paste("biom convert", "-i", inFile, "-o", outFile, "--to-hdf5", 
               '--table-type="OTU table"', "--process-obs-metadata taxonomy"))
  
}
## Join ASV and Taxonomy tables into one
asvTaxTbll <- cbind(asvTblNoChim2, "taxonomy" = taxTbl2[,-1]) # exclude the "ID" first column from "taxTbl2" because "asvTblNoChim2" has already this information
write.table(x = asvTaxTbll, file = "./output/asvTaxTbl.txt", sep = "\t", row.names = FALSE, quote = FALSE) # save ASV-taxonomy tables

convertTab2Biom(inFile = "./output/asvTaxTbl.txt", outFile = "./output/asvTable.biom")

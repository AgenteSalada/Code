# All graphics plotted here are showed in 16S rRNA gene amplicon using dada2pipeline Rplots file
# Installing and importing packages 


library(dada2)
path <- "C:\\Users\\mclar\\Downloads\\OneDrive-2020-01-08 (1)"
list.files(path)
#[1] "B1IV_R1.fastq"  "B1IV_R2.fastq"  "B1V_R1.fastq"   "B1V_R2.fastq"   
#"B2II_R1.fastq"  "B2II_R2.fastq"  "B2III_R1.fastq" "B2III_R2.fastq"
#[9] "B3I_R1.fastq"   "B3I_R2.fastq"   "B3III_R1.fastq" "B3III_R2.fastq" 


#Now we read in the names of the fastq files, and perform some string manipulation to get matched lists of the forward and reverse fastq files.
#Forward and Reverse fastq filename have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001_fastq
fnFs <- sort(list.files(path,pattern = '_R1.fastq', full.names = TRUE))
fnRs <- sort(list.files(path,pattern = '_R2.fastq', full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#INSPECT READ QUALITY PROFILE 
#Now we start by visualizing the quality profiles of the forward reads:
plotQualityProfile(fnFs)
#Rplot1

#Now we visualize the quality profile of the reverse reads:
plotQualityProfile(fnRs)
#Rplot2

#FILTER AND TRIM
#In order to improve the quality of reads, a purification step was carried out
#Looking at the fastQC graphics of the read foward and reverse, we chose to make cuts from 240 bp for both reads
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,240),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows we set multithread=FALSE
head(out)

##                      reads.in  reads.out
## B1IV_R1.fastq       	 128413	    100201
## B1V_R1.fastq	         119773	     95330
## B2II_R1.fastq	       118935	     92027
## B2III_R1.fastq      	 119098	     95522
##B3I_R1.fastq	         131019      93694
## B3III_R1.fastq	       136112	    100492



# LEARNING THE ERROR RATES
#every amplicon dataset has a different set of error rates. The learnErrors method learns this error model from the data, 
#by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution.

errF <- learnErrors(filtFs, multithread=FALSE)
## 114425760 total bases in 476774 reads from 5 samples will be used for learning the error rates.

errR <- learnErrors(filtRs, multithread=FALSE)
## 114425760 total bases in 476774 reads from 5 samples will be used for learning the error rates.

#visualize the estimated error rates
plotErrors(errF, nominalQ=TRUE)
## Rplot3

plotErrors(errR, nominalQ=TRUE)
## Rplot4

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

#Ploting abundance
seqDerepFwd <- countUniqueAndTotalFromDerepObjcList(deRepObj = derepFwd)
uniqDerepFwd <- seqDerepFwd$unique
total <- seqDerepFwd$total 
seqDerepRev <- countUniqueAndTotalFromDerepObjcList(derepRev) 
uniqDerepRev <- seqDerepRev$unique 
derepMtxAbund <- rbind(total, uniqDerepFwd, uniqDerepRev)
barplot(derepMtxAbund, main = "Unique/dereplicated sequences (absolute abundance)", xlab = "Samples", ylab = "Unique sequences (absolute abundance)", col = c("green","blue", "lightblue"), legend = rownames(x = derepMtxAbund), beside=TRUE)
## Rplot5

#Ploting percentage
Unique_Fwd_Perc <- uniqDerepFwd / total * 100
Unique_Rev_Perc <- uniqDerepRev / total * 100
derepMtxPerc <- rbind(Unique_Fwd_Perc, Unique_Rev_Perc) 
barplot(derepMtxPerc, main = "Unique/dereplicated sequences (%)", xlab = "Samples", ylab = "Unique sequences (%)", col = c("gold", "magenta"), legend = rownames(x = derepMtxPerc), beside=TRUE)
## Rplot6

#It was then observed that the number of total strings decreases approximately 80% for fowards and 50% for reverse

#DENOISE UNIQUE SEQUENCES (AND EXCLUDE SINGLETONS)
#Using the error rates obtained previously, we used the function dada to do a denoising of the amplicons, 
#that is, remove the sequences with errors

dadaFwd <- dada(derepFwd, errF, multithread = TRUE) # denoise fwd seqs
## Sample 1 - 100201 reads in 28062 unique sequences.
## Sample 2 - 95330 reads in 19427 unique sequences.
## Sample 3 - 92027 reads in 19221 unique sequences.
## Sample 4 - 95522 reads in 15457 unique sequences.
## Sample 5 - 93694 reads in 17859 unique sequences.
## Sample 6 -100492 reads in 17223 unique sequences

dadaRev <- dada(derepRev, errR, multithread = TRUE) # denoise rev seqs
## Sample 1 - 100201 reads in 48359 unique sequences.
## Sample 2 - 95330 reads in 40975 unique sequences.
## Sample 3 - 92027 reads in 46015 unique sequences.
## Sample 4 - 95522 reads in 41347 unique sequences.
## Sample 5 - 93694 reads in 45736 unique sequences.
## Sample 6 - 100492 reads in 45484 unique sequences.


#MERGE DENOISED FORWARD AND REVERSE READS
#After denoising, we will merge our forward reads with the respective reverse-complement of reverse reads.
#By default, if forward and the reverse-complement of reverse reads does not overlap, reads are discarded as well forward and reverse reads that does not have their pair.
mergePE <- mergePairs(dadaF = dadaFwd, derepF = derepFwd, dadaR = dadaRev, derepR = derepRev, verbose = TRUE)
## 86882 paired-reads (in 4533 unique pairings) successfully merged out of 92628 (in 8745 pairings) input.
## 87495 paired-reads (in 2598 unique pairings) successfully merged out of 89687 (in 3865 pairings) input.
## 85521 paired-reads (in 2496 unique pairings) successfully merged out of 87217 (in 3319 pairings) input.
## 90245 paired-reads (in 1506 unique pairings) successfully merged out of 91769 (in 2194 pairings) input.
## 87013 paired-reads (in 2242 unique pairings) successfully merged out of 89390 (in 3437 pairings) input.
## 94577 paired-reads (in 1667 unique pairings) successfully merged out of 96520 (in 2670 pairings) input.

countMergedSeqsFromDadaObjcList <- function(mergedObj) {
  
  ### count the total no. of merged seqs in the dada2 object
  
  merged <- sapply((sapply(X = mergedObj, `[`, 2)), sum)
  mergedOut <- as.numeric(merged)
  names(mergedOut) <- names(mergedObj)
  
  return(mergedOut)
  
}

## minimum
min(countMergedSeqsFromDadaObjcList(mergePE))
## [1] 85521
## maximum

max(countMergedSeqsFromDadaObjcList(mergePE))
## [1] 9457

## Average of percentage of sequences lost 

mean((out[,1] - countMergedSeqsFromDadaObjcList(mergePE)) / out[,1] * 100)
## [1] 29.28574

#to plot it instead, we can also do the following (in this case we represented the sequences kept):

#Ploting abundance
mergedPE <- countMergedSeqsFromDadaObjcList(mergePE) # get the number of merged reads
countTable <- cbind(out, mergedPE)
colnames(countTable) <- c("Initial", "Filtered", "Merged")
rownames(countTable) <- sample.names
countTable <- t(countTable)
barplot(countTable, col = c("#EFF3FF", "#BDD7E7", "#6BAED6"), main = "Comparison of no. of sequences kept across pipeline steps: initial, filtered and merged",
         ylab = "Absolute no. of sequences", xlab = "Samples", legend = rownames(countTable), beside = TRUE)
## Rplot7

# Doing the same plot but in percentage 
countTablePerc <- apply(X = countTable, MARGIN = 1, FUN = function(x) x / countTable[1,] * 100) 
countTablePerc <- t(countTablePerc)
barplot(countTablePerc, col = c("#EFF3FF", "#BDD7E7", "#6BAED6"), main = "Comparison of proportion of sequences kept across pipeline steps: initial, filtered and merged", 
        ylab = "Percentage of sequences (%)", xlab = "Samples", legend = rownames(countTablePerc), beside = TRUE)
abline(h = min(countTablePerc[3,]))
## Rplot8

#Then, proceed to the construction of an ASV table using the makeSequenceTable () function to tabulate the ASVs using the assembly sequence and the getSequences () function to organize the table

### CONSTRUCT AN ASV TABLE 
# Then, proceed to the construction of an ASV table using the 'makeSequenceTable' function to tabulate the ASVs, 
# using the assembly sequence and the 'getSequences' function to organize the table
asvTbl <- makeSequenceTable(samples = mergePE) # tabulate ASVs

histSeqLen <- table(nchar(getSequences(asvTbl)))

dim(asvTbl)
## 6 11651
#The table is composed of six samples and has a total of 11651 ASVs

#REMOVE CHIMERAS
# Chimeras are formed during PCR amplification when a sequence is incompletely amplified
# Chimeras must be removed separately from the ASV table using the removeBimeraDenovo () function
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


##       input   filtered  denoisedF  denoisedR  merged  nonchim
## B1IV	 128413	   100201	     95837	    94914		86882    35032
## B1V	 119773	    95330	     92254	    90954	  87495	   47313
## B2II	 118935	    92027	     89227	    88537	  85521	   35555
## B2III 119098	    95522	     93315	    92732	  90245	   40582
## B3I	 131019	    93694	     91658	    90565	  87013	   34550
## B3III 136112	    100492	   98519	    97500	  94577	   38301

## create a second summary table seq with one column for the samples 
summaryTblSeq2 <- cbind("Samples" = sample.names, summaryTblSeq)

dir.create("output2") # let's create 'output' folder
write.table(x = summaryTblSeq2, file = "./output/summaryTblSeq.tsv", sep = "\t", row.names = FALSE)

knitr::kable(summaryTblSeq)


# Barplot with the abs. abundance of sequences
summaryTblSeqTrans <- t(summaryTblSeq) # transpose 
barplot(summaryTblSeqTrans, main = "Absolute no. of sequences kept through the pipeline", ylab = "Absolute no. of sequences", xlab = "Samples", col = c("gray", "#EFF3FF", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C"),
        legend = rownames(summaryTblSeqTrans), beside = TRUE)
## Rplot9

# Barplot in percentage
summaryTblSeqPerc <- apply(X = summaryTblSeq, MARGIN = 2, function(x) x / summaryTblSeq[,1] * 100) # get the correspondent percentage table
summaryTblSeqPercTrans <- t(summaryTblSeqPerc) # transpose
barplot(summaryTblSeqPercTrans, main = "Percentage of sequences kept through the pipeline", ylab = "Percentage of sequences (%)", xlab = "Samples", col = c("gray", "#EFF3FF", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C"),
        legend = rownames(summaryTblSeqPercTrans), beside = TRUE) # plot it 
## Rplot10

#ASSINING OUR TAXONOMY
#In this step we will use the SILVA training reference database to assign a taxonomy to our samples.
#The first file assigns a taxonomy to the table and the second adds the species name if there is a match in our ASVs
taxTbl <- assignTaxonomy(seqs = asvTblNoChim, refFasta = "C:\\Users\\mclar\\Downloads\\silva_nr_v132_train_set.fa.gz", multithread = FALSE)
taxTbl <- addSpecies(taxtab = taxTbl, refFasta = "C:\\Users\\mclar\\Downloads\\silva_species_assignment_v132.fa.gz")
taxTbl.print <- taxTbl 
rownames(taxTbl.print) <- NULL
head(taxTbl.print)
##       Kingdom    Phylum           Class                 Order             Family           Genus                  Species
## [1,] "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Pseudomonadales" "Moraxellaceae"  "Acinetobacter"        NA     
## [2,] "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Pseudomonadales" "Moraxellaceae"  "Acinetobacter"        NA     
## [3,] "Bacteria" "Cyanobacteria"  "Oxyphotobacteria"    "Chloroplast"     NA               NA                     NA     
## [4,] "Bacteria" "Cyanobacteria"  "Oxyphotobacteria"    "Nostocales"      "Microcystaceae" "Microcystis_PCC-7914" NA     
## [5,] "Bacteria" "Cyanobacteria"  "Oxyphotobacteria"    "Chloroplast"     NA               NA                     NA     
## [6,] "Bacteria" "Cyanobacteria"  "Oxyphotobacteria"    "Nostocales"      "Microcystaceae" "Microcystis_PCC-7914" NA   

# Save ASV and taxonomy tables in R format 
saveRDS(object = asvTblNoChim, file = "./output/asvTblNoChim.rds") # save the ASV table
saveRDS(object = taxTbl, file = "./output/taxTbl.rds") # save the ASV taxonomy

#Taxonomy - Construction of the ASV table
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


source("C:\\Users\\mclar\\Downloads\\Rscript.R") # import R script with built-in functions
taxTbl2 <- tax2biom(taxTbl2) 

## Join ASV and Taxonomy tables into one
asvTaxTbll <- cbind(asvTblNoChim2, "taxonomy" = taxTbl2[,-1]) # exclude the "ID" first column from "taxTbl2" because "asvTblNoChim2" has already this information
write.table(x = asvTaxTbll, file = "./output/asvTaxTbl.txt", sep = "\t", row.names = FALSE, quote = FALSE) # save ASV-taxonomy tables
##        B1IV  B1V  B2II B2III    B3I   B3III   ASV_ID   taxonomy
## ASV_1	2070	3345	38	  47	     0	     0	   ASV_1	 k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Acinetobacter; s__NA; a__ASV_1
## ASV_2	1636	3610	56	  63	     0	     0	   ASV_2	 k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Acinetobacter; s__NA; a__ASV_2
## ASV_3	  42    37 483	4442     218	    62	   ASV_3	 k__Bacteria; p__Cyanobacteria; c__Oxyphotobacteria; o__Chloroplast; f__NA; g__NA; s__NA; a__ASV_3
## ASV_4	   0	   0	20	   0    1529	  3551	   ASV_4	 k__Bacteria; p__Cyanobacteria; c__Oxyphotobacteria; o__Nostocales; f__Microcystaceae; g__Microcystis_PCC-7914; s__NA; a__ASV_4
## ASV_5	  40	  42 562	4337	    87	    31	   ASV_5	 k__Bacteria; p__Cyanobacteria; c__Oxyphotobacteria; o__Chloroplast; f__NA; g__NA; s__NA; a__ASV_5
## ASV_6	   0	   0	 0	   0	  2303	  2715	   ASV_6	 k__Bacteria; p__Cyanobacteria; c__Oxyphotobacteria; o__Nostocales; f__Microcystaceae; g__Microcystis_PCC-7914; s__NA; a__ASV_6
## ASV_7	   0	   0	47	   0	   950	  3767	   ASV_7	 k__Bacteria; p__Cyanobacteria; c__Oxyphotobacteria; o__Nostocales; f__Microcystaceae; g__Microcystis_PCC-7914; s__NA; a__ASV_7
## ASV_8	  43	  58 482	3735	   205	    87	   ASV_8	 k__Bacteria; p__Cyanobacteria; c__Oxyphotobacteria; o__Chlor


### Import metadata and put it in a biom format too
# After obtaining the ASV tables, the biomformat software would be used to convert the files into .biom format for analysis using Biome-Shiny

metadata <- read.table("C:\\Users\\mclar\\Downloads\\miseqsopdata (3)\\MiSeq_SOP\\mouse.time.design", header = TRUE)
rownames(metadata) <- metadata[,1]
colnames(metadata) <- c("SampleID", "Condition")
write.csv(x = metadata, file = "./output/metadata.csv", quote = FALSE, row.names = TRUE)


BiocManager::install('biomformat')
#install.packages('command -v biom')
install.packages ('biomformat')



library('biomformat')
convertTab2Biom <- function(inFile, outFile) {
  
  if (system("command - v biomformat", ignore.stdout = TRUE, ignore.stderr = TRUE) !=0)  {
    
    stop("biom program is not installed or it is not accessible!\n  Exiting...")
    
  }
  
  system(paste("biom convert", "-i", inFile, "-o", outFile, "--to-hdf5", 
               '--table-type="OTU table"', "--process-obs-metadata taxonomy"))
  
}


convertTab2Biom(inFile = ".\\output\\asvTaxTbl.txt", outFile = ".\\output\\asvTable.biom")
## this procedure did not go as expected because we had an error 
#Error in convertTab2Biom(inFile = ".\\output\\asvTaxTbl.txt", outFile = ".\\output\\asvTable.biom") : 
  #biom program is not installed or it is not accessible!
  #Exiting...
#In addition: Warning message:
  #In system("command - v biomformat", ignore.stdout = TRUE, ignore.stderr = TRUE) :
  #'command' not found
# this procedure did not go as expected due to R being unable to locate this package and not having found the access command to execute.




library(dada2); packageVersion("dada2") 

#Change to the directory containing the fastq files
path <- "$Sequence_File_Directory" list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

plotQualityProfile(fnFs[1:5]) # visualizing the quality profiles of the forward reads
#plotQualityProfile(fnFs)
plotQualityProfile(fnRs[1:5]) # visualizing the quality profiles of the reverse reads
#plotQualityProfile(fnRs)
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#Choose points where they will both over lap don't choose numbers in the bad zone
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,180),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
out

errF <- learnErrors(filtFs, multithread=TRUE) # estimate error rates in forward sequences
errR <- learnErrors(filtRs, multithread=TRUE) # estimate error rates in reverse sequences

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE) # visualize estimated error rates

# dereplicate filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)


# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE) #Infer the sequence variants in each sample

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]

#merge denoised forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


#Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


###make PICRUSt2 required files after running DADA2###
#PICRUST2 needs: seqfile.fna, biomfile.biom 
library(biomformat);packageVersion('biomformat')
library(ShortRead);packageVersion('ShortRead')

#biom file
zymo.biom <- make_biom(t(seqtab.nochim))
write_biom(zymo.biom, "$File_directory/biom_file.biom")

#fasta files 
sq <- getSequences(seqtab.nochim)
id <- paste0("Abundance=", colSums(seqtab.nochim))
names(sq) <- id
writeFasta(sq, file="$File_directory/fasta_file.fna")

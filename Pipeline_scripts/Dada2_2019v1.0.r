#!/usr/bin/env Rscript

# RUN with Rscript Dada2_2019v1.0.r Dir_path
args = commandArgs(trailingOnly=TRUE)
#test there was input
if (length(args)==0) {
  stop("At least one argument must be supplied (input directory).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

path <- args[1] #Directory to run Dada2 over
print(path)
truncF <- as.numeric(args[2]) #trunc len for forward reads from param file
truncR <- as.numeric(args[3]) #trunc len for reverse reads from param file

list.files(path=path, pattern='*trimmed', recursive=TRUE)

#Run dada2

library(dada2); packageVersion("dada2")

#Set Dada2 Paramaters to increase sensitivity
#setDadaOpt(OMEGA_A = 1e-4)
#setDadaOpt(OMEGA_C = 1e-4)
#setDadaOpt(DETECT_SINGLETONS = TRUE)

#Based off these tutorials:
#https://benjjneb.github.io/dada2/tutorial.html
#https://astrobiomike.github.io/amplicon/dada2_workflow_ex

# Forward and reverse fastq filenames have format: "${C_LIB}"_R1_trimmed.fq and "${C_LIB}"_R2_trimmed.fq
fnFs <- sort(list.files(path=path, pattern="_R1_trimmed.fq", full.names = TRUE, recursive=TRUE))
fnRs <- sort(list.files(path=path, pattern="_R2_trimmed.fq", full.names = TRUE, recursive=TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
print(sample.names)

#We start by visualizing the quality profiles of the forward reads:
pdf(file = "Qual_profiles.pdf")
plotQualityProfile(fnFs[1:2])

#Now we visualize the quality profile of the reverse reads:
plotQualityProfile(fnRs[1:2])

dev.off()
#Assign the filenames for the filtered fastq.gz files.
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#COI
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(truncF,truncR),
                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=TRUE)

#learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

pdf(file = "errF_errR.pdf")
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

print('Run Sample Inference')
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool=TRUE)

#merge paired reads
print('Merge Paired Reads')
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
print('Dimensions of Sequence Table:')
dim(seqtab)

# Inspect distribution of sequence lengths
print('Distribution of Sequence Lengths: ')
table(nchar(getSequences(seqtab)))

#remove chimeras
print('Remove Chimeras')
#Consensus
print('removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)')
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Pooled
#print('removeBimeraDenovo(seqtab, method="pooled", multithread=TRUE, verbose=TRUE)')
#seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

print('Percent Non-Chimera Reads')
sum(seqtab.nochim)/sum(seqtab)

#Track reads through pipeline
print('Track Read counts through Dada2')
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
#head(track)
#export Table
write.table(track, "ASVs_track.tsv", sep="\t", quote=F, col.names=NA)
track

#Export ASV_table and fasta file

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

print('Finished with Dada2 processing.')

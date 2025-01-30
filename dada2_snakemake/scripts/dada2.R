# pRun snakemake
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(dada2)
library(stringr)


path <- snakemake@input[[1]]
fnFs <- sort(list.files(path, pattern="R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#sample_some_seq <- sample(6,1:length(fnFs),replace=F)

#dada2::plotQualityProfile(fnFs[sample_some_seq])
#dada2::plotQualityProfile(fnRs[sample_some_seq])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


left2trunk <- 20
  
right2trunk <- snakemake@params[[1]]

out <- filterAndTrim(fnFs, 
                    filtFs,
                    fnRs,
                    filtRs,
                    trimLeft=left2trunk,
                    trimRight=right2trunk,
                    maxN=0,
                    maxEE=c(2,2),
                    truncQ=2,
                    rm.phix=TRUE,
                    compress=TRUE,
                    multithread=snakemake@params[[2]]) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=snakemake@params[[2]])
errR <- learnErrors(filtRs, multithread=snakemake@params[[2]])

dadaFs <- dada(filtFs, err=errF, multithread=snakemake@params[[2]])
dadaRs <- dada(filtRs, err=errR, multithread=snakemake@params[[2]])

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

tabseq <- makeSequenceTable(mergers)

# Inspect distribution of sequence lengths
table(nchar(getSequences(tabseq)))

seqtab.nochim <- removeBimeraDenovo(tabseq, method="consensus", multithread=snakemake@params[[2]], verbose=TRUE)

sum(seqtab.nochim)/sum(tabseq)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(tabseq,
                      snakemake@params[[3]],
                      multithread=snakemake@params[[2]],
                      outputBootstraps = T)

taxa$tax <-cbind(sequences=rownames(taxa$tax),taxa$tax)
colnames(taxa$boot) <- paste0(colnames(taxa$boot),"_boots")
taxa <- cbind(taxa$tax,taxa$boot)
# seqtab.nochim$samples <- rownames(seqtab.nochim)

write.csv2(tabseq,snakemake@output[[3]])
write.csv2(seqtab.nochim, snakemake@output[[2]])
write.csv2(taxa, snakemake@output[[1]],row.names = F)
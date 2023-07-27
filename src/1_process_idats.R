

#!/usr/bin/env Rscript

## Library
suppressMessages(library(minfi))
suppressMessages(library(conumee))

args = commandArgs( trailingOnly=TRUE )

if (length(args)!= 1) {
  stop("At least one argument must be supplied (input dir).n", call.=FALSE)
}

idat_path= args[1]
stopifnot(dir.exists(idat_path))

## load list of bad probes
nonspec = read.csv("../../Datasets/Methylation/illumina450k_filtering/48639-non-specific-probes-Illumina450k.csv")
multimap = read.csv("../../Datasets/Methylation/illumina450k_filtering/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt", header = FALSE)


process_methylation = function(file_paths, out_name){
for ( i in file_paths){
  makeMethSheet(dir = i, group = "cancer", outname = paste0(i,"RefSheet.csv" ) )
  # Record platform
  targets = read.metharray.sheet(i)
  RGset = read.metharray.exp(targets=targets, force=TRUE)
  
  # Flag probes which failed in more than 50% samples
  detP = detectionP(RGset )
  failed = detP > 0.01
  colMeans(failed) # Fraction of failed positions per sample
  sum(rowMeans(failed)>0.5)
  failed.probes = rownames(detP[rowMeans(failed) > 0.5,])
  
  bad.probes = unique(c(failed.probes, nonspec$TargetID, multimap$V1))
  
  Mset = preprocessNoob(RGset, dyeMethod = "single")
  Rset = ratioConvert(Mset, what = "beta", keepCN = FALSE)
  GRset = mapToGenome(Rset)
  GRset = dropLociWithSnps(GRset, snps=c("SBE","CpG"), maf=0 )
  GRset = dropMethylationLoci(GRset, dropRS = TRUE, dropCH = TRUE)
  pheno = pData(GRset)
  beta = getBeta(GRset )
  
  ## Remove probes for X and Y / bad probes
  dropvec = as.character(seqnames(GRset@rowRanges)) %in% c("chrX", "chrY")
  xyprobes = names(GRset@rowRanges)[dropvec]
  beta = beta[!( rownames(beta) %in% c( bad.probes, xyprobes) ), ]
  
  outname = paste0( out_name, ".rds")
  saveRDS(beta,  outname)
  }
}  


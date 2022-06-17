setwd("K:/")

#Since peak calling setting varies from ChIP-Seq to ChIP-Seq, we recommand that the user follows the "csaw" tutorial to generate
#the differential binding sites using the settings that best suits its data.

#Alternatively, any bed file can serve as the input at this step.

#Modify the reference genome as needed.

#####rGADEM
library(rGADEM)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(rtracklayer)
library(IRanges)
library(BiocGenerics)
library(GenomicAlignments)
library(GenomicRanges)
library(regioneR)
library(dplyr)


#List the name of your bed files (or converted csaw output) placed in the "bed"
mlist = c( "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8") 

for(i in mlist){
  gr = toGRanges(paste0("bed/",i,".bed"))
  Seqs_see=Biostrings::getSeq(Hsapiens, GRanges(gr) )
  time.0 <- Sys.time()
  gadem<-GADEM(Sequences = Seqs_see,verbose=1,genome=Hsapiens)
  print(Sys.time()-time.0)
  save(gadem, file = paste0("GADEM/GADEM_",i,".RData")) 
  
}



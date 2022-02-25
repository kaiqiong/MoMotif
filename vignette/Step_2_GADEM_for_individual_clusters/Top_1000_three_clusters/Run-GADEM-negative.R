# Run the analysis using GADEM


.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")

# Object out.ranges stores all the genomic Ranges -- a total of 78601 ranges
out.ranges<- readRDS("/scratch/greenwood/kaiqiong.zhao/kaiqiong.zhao/Projects/WitcherKI/Results/out_ranges_June_28.rds")
#o <- order(out.ranges$PValue) 
# The ordered results are in 
#out.ranges[o]


# the best is the index of the window that is the most significant in each cluster,
#bets.logFC is the log-fold change of that 'best' window


library(rGADEM)
library(BSgenome.Hsapiens.UCSC.hg19)
#library(rtracklayer)
#--------------------------------------------------------------------
# Step 1: Generate a 'XString' Object from the positively significant clusters from csaw package
#--------------------------------------------------------------------
#Object 'Sequences' should be of type 'XStringViews', 'DNAStringSet' or 'RangedData'

out_order <- out.ranges[order(out.ranges$PValue),]

# The top 100 down/up sequences

out_down <-  out_order[out_order$direction =="down",]
out_up <-  out_order[out_order$direction =="up",]

gr_up <- out_up[1:1000,]
gr_down <- out_down[1:1000, ]
gr_null <- out_order[(length(out_order)-999):length(out_order),]
# convert GRanges back to a RangedData


gr <- gr_down


#see<-RangedData(gr@ranges,space = gr@seqnames)

#see <- GRanges(seqnames=gr@seqnames, ranges = gr@ranges, seqinfo = gr@seqinfo, strand = gr@strand)

Seqs_see=Biostrings::getSeq(Hsapiens, GRanges(gr) )

# Run GADEM

time.0 <- Sys.time()
gadem<-GADEM(Sequences = Seqs_see,verbose=1,genome=Hsapiens)
print(Sys.time()-time.0)
setwd("/scratch/greenwood/kaiqiong.zhao/kaiqiong.zhao/Projects/WitcherKI/GADEM_top_1000")
save(gadem, file = "GADEM_down_1000.RData")


load("/scratch/greenwood/kaiqiong.zhao/kaiqiong.zhao/Projects/WitcherKI/2-GADEM_top_1000/GADEM_down_1000.RData")

gadem@motifList
length(gadem@motifList)  # the algorithm found 4 motif?


consensus(gadem)

startPos(gadem)

motifs <- getPWM(gadem)

library(MotIV)

# motifMatch -- search for motifs matches (from the JASPER dataset) corresponding to the motifs found from GADEM

summary(motifMatch(motifs))

MotIV.obj <- motifMatch(motifs)
# plot function provides a summary of each identified TF (from JASPER dataset) associated to the input motifs, 
#the sequence logo, the name of the motif match and the p-value of the alignment
setwd("/scratch/greenwood/kaiqiong.zhao/kaiqiong.zhao/Projects/WitcherKI/2-GADEM_top_1000")
pdf("Negative.pdf", width = 10, height=12)

plot(MotIV.obj, top=5, main ="Logo (top matching JASPAR)", bysim=T) # the default reference database if JASPAR

# identification of primary motifs from the JASPAR database

dev.off()
# ALIGNMENT E-values

# a visualization of the repartition of TF found

plot(motifMatch(motifs), gadem, ncol = 2)

#viewAlignments(motifMatch(motifs))

plot(motifMatch(motifs), gadem, ncol = 2, type ='distance')








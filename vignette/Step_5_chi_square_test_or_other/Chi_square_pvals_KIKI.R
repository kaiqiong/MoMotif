library(magrittr)
library(rGADEM)
library(ggplot2)
library(reshape2)
#library(MotIV)
library(BSgenome.Hsapiens.UCSC.hg19)

# 1- Summary without any extensions

#load("~/Documents/Packages/diffMotif/vignette/Step_2_GADEM_for_individual_clusters/Top_1000_three_clusters/GADEM_up_1000.RData")
load("/Users/kaiqiongzhao/Documents/Projects_since_2021/WitcherKIKI/GADEM_top_1000/GADEM_up_1000.RData")
gadem_up <- gadem

load("/Users/kaiqiongzhao/Documents/Projects_since_2021/WitcherKIKI/GADEM_top_1000/GADEM_down_1000.RData")
#load("~/Documents/Packages/diffMotif/vignette/Step_2_GADEM_for_individual_clusters/Top_1000_three_clusters/GADEM_down_1000.RData")
gadem_down <- gadem
#load("~/Documents/Packages/diffMotif/vignette/Step_2_GADEM_for_individual_clusters/Top_1000_three_clusters/GADEM_null_1000.RData")
load("/Users/kaiqiongzhao/Documents/Projects_since_2021/WitcherKIKI/GADEM_top_1000/GADEM_null_1000.RData")
gadem_null <- gadem

motif11 <- gadem_up@motifList[[1]]
motif21 <- gadem_down@motifList[[1]]
motif31 <- gadem_null@motifList[[1]]

summary_motif <- function(motif1){
  seqID_vec <- rep(NA, length(motif1@alignList))
  seq_vec <- strand_vec <- pval <- seqID_vec
  
  for ( i in 1:length(seqID_vec)){
    seqID_vec[i] <- motif1@alignList[[i]]@seqID
    pval[i] <- motif1@alignList[[i]]@pval
    strand_vec[i] <- motif1@alignList[[i]]@strand
    seq_vec[i] <- motif1@alignList[[i]]@seq
  }
  
  out = list(seq_vec=seq_vec, strand_vec = strand_vec, seqID_vec = seqID_vec)
  return(out)
}

rr1 = summary_motif(motif11)
rr2 = summary_motif(motif21)
rr3 = summary_motif(motif31)

# 1-5 --- check the identified motif by GADEM -- before extending




motifs <- getPWM(gadem_up)
seqLogo(motifs[[1]])
seqLogo(consensusMatrix( (DNAStringSet(rr1$seq_vec )), as.prob = TRUE)[1:4,])   #  --- direction-- reverse

seqLogo(consensusMatrix( reverseComplement(DNAStringSet(rr1$seq_vec )), as.prob = TRUE)[1:4,])   #  --- direction-- reverse

seqLogo(getPWM(gadem_down)[[1]])  # ---- direction --- not reverse
seqLogo(consensusMatrix( reverseComplement(DNAStringSet(rr2$seq_vec )), as.prob = TRUE)[1:4,])

seqLogo(getPWM(gadem_null)[[1]]) # ---- direction --- not reverse
seqLogo(consensusMatrix( reverseComplement(DNAStringSet(rr3$seq_vec )), as.prob = TRUE)[1:4,])





# 2- Load the raw sequences

out.ranges<- readRDS("/Users/kaiqiongzhao/Documents/Projects_since_2021/WitcherKIKI/Results/out_ranges_July_8.rds")
out_order <- out.ranges[order(out.ranges$PValue),]

# The top 100 down/up sequences

out_down <-  out_order[out_order$direction =="down",]
out_up <-  out_order[out_order$direction =="up",]

gr_up <- out_up[1:1000,]
gr_down <- out_down[1:1000, ]
gr_null <- out_order[(length(out_order)-999):length(out_order),]

# 3- Extension based on Biostrings::getSeq

Seqs_up=Biostrings::getSeq(Hsapiens, GRanges(gr_up) + 60 )
Seqs_down=Biostrings::getSeq(Hsapiens, GRanges(gr_down) + 60 )
Seqs_null=Biostrings::getSeq(Hsapiens, GRanges(gr_null) +60 )  


# 4- The main step: Use the extended sequences, mapped on the rr1$seq_vec, i.e. the core motif identified from the original rGADEM

seq_set_extract <- function(rr1, motif11, Seqs_up,bp.use = 21, reverse = FALSE, rrWidth=13, start = 2, end = 14){
  ext = (bp.use-1)/2 - (rrWidth-1)/2
  #if(width(rr1$seq_vec) %>% unique()==13){  # up sequences identified by rGADEM
  
  #}
  #if(width(rr1$seq_vec) %>% unique() ==18){ # down sequences
  #  start = 1; end = 11; center = 6
  #}
  #if(width(rr1$seq_vec) %>% unique() ==14){ # null sequences
  #  start = 3; end = 13; center = 8
  #}
  seqs_set <- rep(NA, length(rr1$strand_vec))
  
  for ( i in 1:length(rr1$strand_vec)){
    if(rr1$strand_vec[i] == "+"){
      seq_nn <- Seqs_up[motif11@alignList[[i]]@seqID] #the raw input long sequence
      #rr1$seq_vec[i]
      oo = gregexpr(DNAString( rr1$seq_vec[i]), seq_nn)
      
      seqs_set[i]<- subseq(seq_nn , start= oo[[1]][1] + (start-1) - ext, end = oo[[1]][1]+(end-1) + ext )
      
      # id_now = rr1$seqID_vec[i]
      #seqs_set[i]<- subseq(Seqs_up[id_now],
      ##   start= motif11@alignList[[i]]@pos+(start-1) - ext ,
      #  motif11@alignList[[i]]@pos+(end-1) + ext  )
    }
    if(rr1$strand_vec[i] == "-"){
      seq_nn <- Seqs_up[motif11@alignList[[i]]@seqID] #the raw input long sequence -- this is in the reverse strand
      
      #oo = gregexpr(DNAString( rr1$seq_vec[i])%>% reverseComplement(), seq_nn)
      #if(oo[[1]][1] != -1){
      
      #seqs_set[i]<-subseq(seq_nn %>% reverseComplement(), start= width(seq_nn)-oo[[1]][1]+1-10 -ext-(width(rr1$seq_vec[i])-end) , #width(seq_nn)-oo[[1]][1]+1 + ext -(width(rr1$seq_vec[i])-end) )
      
      
      # Modify the reverse strand extraction-- July 1 reverse the original input strand instead --> advantage: less tuning for start and end for the - strand
      
      oo = gregexpr(DNAString(rr1$seq_vec[i]), seq_nn %>% reverseComplement())
      if(oo[[1]][1] != -1){
        
        seqs_set[i]<- subseq(seq_nn %>% reverseComplement() , start= oo[[1]][1] + (start-1) - ext, end = oo[[1]][1]+(end-1) + ext )
        
        
        
        # if(width(rr1$seq_vec) %>% unique()==18){
        
        #  subseq(seq_nn , start= oo[[1]][1] + (start-1) - ext, end = oo[[1]][1]+(end-1) + ext ) %>% reverseComplement()
        #seqs_set[i]<-subseq(seq_nn %>% reverseComplement(), start= width(seq_nn)-oo[[1]][1], width(seq_nn)-oo[[1]][1]+1+11 )
        #}
        
      }
    }
    
  }
  seqs_set <- seqs_set[!is.na(seqs_set)]
  
  if(reverse){
    res =consensusMatrix( reverseComplement(DNAStringSet(seqs_set )), as.prob = TRUE)[1:4,]
    
    res_counts = consensusMatrix( reverseComplement(DNAStringSet(seqs_set )), as.prob = FALSE)[1:4,]
  }else{
    res =consensusMatrix(DNAStringSet(seqs_set ), as.prob = TRUE)[1:4,]
    res_counts =consensusMatrix(DNAStringSet(seqs_set ), as.prob = FALSE)[1:4,]
  }
  res_melt1<-melt(res)
  colnames(res_melt1)<- c("Nucleotide", "Position", "Freq.")
  return(list(res = res,res_counts = res_counts, res_melt = res_melt1, seqs_set=seqs_set)) 
  # April 5, 2021, extracting the raw sequence matrix as well
}



res_ext1 <- seq_set_extract(rr1, motif11, Seqs_up, bp.use = 61, reverse = FALSE, rrWidth = 11, start = 3, end = 13)
#-------
# for the up cluster -- the orignal rGADEM identified sequences of length 1-14, center is on position 8
#----------------

# start: start in the original subsequences rr1$seq_vec[i]
# end: end in the original subsequences rr1$seq_vec[i]
# make everything symetric; facilitate any follow-up manipulation for the reverse strand

# For the reverse = TRUE sequences, carefully select the start and end

res_ext2 <- seq_set_extract(rr2, motif21, Seqs_down, bp.use = 61, reverse = FALSE, rrWidth = 11, start = 2, end = 12)

# Down cluster -- 16 bp, center is on position 8

res_ext3<- seq_set_extract(rr3, motif31, Seqs_null, bp.use = 61, reverse = FALSE,
                           rrWidth = 11, start = 3, end = 13)
# 5- Store the observed counts 



apply(res_ext1$res_counts,2, sum)


apply(res_ext2$res_counts,2, sum)


apply(res_ext3$res_counts,2, sum)


#  Lost vs. Stable
pval_vec <- rep(NA, 61)

for(posnow in 1:61){
  temp_tab <- data.frame(lost = res_ext2$res_counts[,posnow], stable= res_ext3$res_counts[,posnow])

  tryCatch({pval_vec[posnow] <- chisq.test(temp_tab)$p.value}, 
           warning=function(w) print(posnow))
  
}

pval_vec_plot <- pval_vec
pval_vec_plot[which(pval_vec_plot==0)] <- 10^(-300)

setwd("~/Documents/Packages/diffMotif/vignette/Step_5_chi_square_test_or_other")
pdf("Chi_square_p_vals_KIKI_Lost_Stable_correct.pdf", height = 5, width = 8)
plot(-log10(pval_vec_plot), pch = 19, xlab = "Position",
     ylab = "-log10(p-value)", 
     main = "Lost cluster v.s Stable Cluster",
     ylim = c(0, 50))
abline(v = 31, col = "blue", lty = 4)
abline(v = 25, lty = 2, col = 1)
abline(v = 48, lty = 2, col = 1)
dev.off()


#  Gained vs. Stable

pval_vec1 <- rep(NA, 61)

for(posnow in 1:61){
  temp_tab <- data.frame(gained = res_ext1$res_counts[,posnow], stable= res_ext3$res_counts[,posnow])
  tryCatch({pval_vec1[posnow] <- chisq.test(temp_tab)$p.value}, 
           warning=function(w) print(posnow))
}


pval_vec_plot <- pval_vec1
pval_vec_plot[which(pval_vec_plot==0)] <- 10^(-300)

#pval_vec1[which(pval_vec1==0)] <- 10^(-300)

pdf("Chi_square_p_vals_KIKI_Gained_Stable_correct.pdf", height = 5, width = 8)
plot(-log10(pval_vec_plot), pch = 19, xlab = "Position",
     ylab = "-log10(p-value)", main = "Gained cluster v.s Stable Cluster",
     ylim = c(0, 50))
abline(v = 31, col = "blue", lty = 4)
abline(v = 25, lty = 2, col = 1)
abline(v = 48, lty = 2, col = 1)
dev.off()



p_vals_KI <- data.frame("Lost.stable" = pval_vec, "Gained.stable" = pval_vec1)

write.csv(p_vals_KI, file="Chi_square_p_vals_KIKI_correct.csv")


plot(-log10(pval_vec_plot), pch = 19)

plot(-log10(p.adjust(pval_vec_plot)), pch = 19)

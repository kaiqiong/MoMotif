setwd("K:/")

library(magrittr)
library(rGADEM)
library(ggplot2)
library(reshape2)
library(BSgenome.Hsapiens.UCSC.hg19)

##### 1- Summary without any extensions #####
##EDIT AND RUN## 
#Vector for the "GADEM_XXX.RData" files in the "GADEM" folder to analyze
gad_file_vec = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")
n_gf = length(gad_file_vec)

#Vector for the gadem_xxx suffix you want to give the files
gad_var_vec = c("c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8")
n_gv = length(gad_var_vec)

##Write the name of all the .bed files (1st column: chr?, 2nd: start, 3rd: end ; no header; in "bed" folder) 
#associated with each gadem file, in the same order
bed_vec = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8" )

#Name of the Frequency Table
outname1 = "Freq_C_All"


##

##VERIFY## 
#Make sure n_gf and n_gv are the same size

##

##RUN## 
#Load the gadem files
for(i in 1:n_gv){
  load(paste0("GADEM/GADEM_",gad_file_vec[i],".RData"))
  assign(paste0("gadem_",gad_var_vec[i]), gadem)
  rm(gadem)
  
}

#Make the list of all gadem files
gad_list = mget(paste0("gadem_",gad_var_vec[1:n_gv] ))

#Define the function
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

#Extract the list of primary motif of every site for each file
for(i in 1:n_gv){
  gad = gad_list[[i]]
  motif1l = gad@motifList[[1]]
  assign(paste0("motif1_",gad_var_vec[i]), motif1l)
  rm(gad)
  rm(motif1l)
  
}

#Make the list of all motif1 files
mot_list = mget(paste0("motif1_",gad_var_vec[1:n_gv] ))

#Make the summary file for each file
for(i in 1:n_gv){
  rr = summary_motif(mot_list[[i]])
  assign(paste0("rr_",gad_var_vec[i]), rr)
  rm(rr)
}

#Make list of all rr files
rr_list = mget(paste0("rr_",gad_var_vec[1:n_gv] ))

#Get the logo of forward and reverse consensus motif of all groups

for(i in 1:n_gv){
  motifs = getPWM(gad_list[[i]])
  seqLogo(motifs[[1]])
  rr = rr_list[[i]]
  seqLogo(consensusMatrix( reverseComplement(DNAStringSet(rr$seq_vec )), as.prob = TRUE)[1:4,])
  rm(motifs)
  rm(rr)
}


##

##Verify##

#Look at the motif in each group and make sure that they are all similar. 
#If not, look at other called motif to find an appropriate one
#You can use the following code to verify this

#EDIT# Suffix of the file to look at, its rank in the given lists and level of motif to look at
suf = "down"
rank = 2
lev = 3
#RUN# all the following lines
gadl=mget(paste0("gadem_",suf))
motifs = getPWM(gadl[[1]])
seqLogo(motifs[[lev]])
#If the good motif is not number 1, save the right one in mot list
gad = gad_list[[rank]]
motif1l = gad@motifList[[lev]]
assign(paste0("motif1_",suf), motif1l)
rr = summary_motif(motif1l)
assign(paste0("rr_",suf), rr)

rr_list = mget(paste0("rr_",gad_var_vec[1:n_gv] ))
mot_list = mget(paste0("motif1_",gad_var_vec[1:n_gv] ))

rm(gadl)
rm(motifs)
rm(suf)
rm(lev)
#

##

##### 2- Load the raw sequences

for(i in 1:n_gv){
  gr = read.table(file = paste0("bed/",bed_vec[i],".bed"), header = F)
  names(gr) = c("chr", "start", "end")
  assign(paste0("gr_", gad_var_vec[i]), gr)
  rm(gr)
}

#Make list of bed gr files
gr_list = mget(paste0("gr_",gad_var_vec[1:n_gv] ))

##### 3- Extension based on Biostrings::getSeq
for(i in 1:n_gv){
  Seqs = Biostrings::getSeq(Hsapiens, GRanges(gr_list[[i]]) + 60 )
  assign(paste0("Seqs_", gad_var_vec[i]), Seqs)
  rm(Seqs)
}

#Make list of Seqs
Seqs_list = mget(paste0("Seqs_",gad_var_vec[1:n_gv] ))

# 4- The main step: Use the extended sequences, mapped on the rr1$seq_vec, i.e. the core motif identified from the original rGADEM
#How big is the region to look at and how wide is the motif
bp.use = 61
rrWidth = 12
#Vectors of the reverse status (TRUE or FALSE), start postion and end position of the motif in the seqLogo output previously used

rev_vec = c(T,F,F,T,T,F,T,F,F,T)
st_vec = c(5,2,2,2,2,2,1,2,2,1)
end_vec = c(16,13,13,13,13,13,12,13,13,12)

#rev_vec = c(T,T,F)
#st_vec = c(2,2,3)
#end_vec = c(13,13,14)

for(k in 1:n_gv){
  rr = rr_list[[k]]
  motif = mot_list[[k]]
  Seqs = Seqs_list[[k]]
  reverse = rev_vec[k]
  start = st_vec[k]
  end = end_vec[k]
  
  
  ext = (bp.use-1)/2 - (rrWidth-1)/2
  seqs_set <- rep(NA, length(rr$strand_vec))
  for ( i in 1:length(rr$strand_vec)){
    if(rr$strand_vec[i] == "+"){
      seq_nn <- Seqs[motif@alignList[[i]]@seqID] #the raw input long sequence
      oo = gregexpr(DNAString( rr$seq_vec[i]), seq_nn)
      ss = oo[[1]][1] + (start-1) - ext
      ee = oo[[1]][1]+(end-1) + ext
      if(ss >= 0 & ee >= 0){
        seqs_set[i]<- subseq(seq_nn , start= ss, end = ee )
      }
      
    }
    if(rr$strand_vec[i] == "-"){
      seq_nn <- Seqs[motif@alignList[[i]]@seqID] #the raw input long sequence -- this is in the reverse strand
      oo = gregexpr(DNAString(rr$seq_vec[i]), seq_nn %>% reverseComplement())
      ss = oo[[1]][1] + (start-1) - ext
      ee = oo[[1]][1]+(end-1) + ext
      if(ss >= 0 & ee >= 0){
        if(oo[[1]][1] != -1){
          seqs_set[i]<- subseq(seq_nn %>% reverseComplement() , start= ss, end = ee )
          
        }
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
  res_ext = list(res = res,res_counts = res_counts, res_melt = res_melt1, seqs_set=seqs_set)
  
  assign(paste0("res_ext_",gad_var_vec[k]), res_ext)
  
  
}

#Make list of res_ext
res_list = mget(paste0("res_ext_",gad_var_vec[1:n_gv] ))



#Make Res table
for(i in 1:n_gv){
  res = res_list[[i]]
  res1 = res$res_melt
  names(res1) = c("Nucleotide", "Position", gad_var_vec[i] )
  if(i == 1){
    restab = res1
  }else{
    restab = cbind(restab,res1[,3])
  }
  
  
}

names(restab) = c("Nucleotide", "Position", gad_var_vec)


write.csv(restab, paste0(outname1,".csv"), row.names = F)


#Make folder for p-value
if (file.exists("pvalue")){
} else {
  dir.create("pvalue")
}
#Make p-value table and graph
for(i in 1:n_gv){
  for(k in 1:n_gv){
    if(i != k){
      pval_vec <- rep(NA, bp.use) 
      for(posnow in 1:bp.use){
        res_ext2 = res_list[[i]]
        res_ext3 = res_list[[k]]
        temp_tab <- data.frame(first = res_ext2$res_counts[,posnow], second= res_ext3$res_counts[,posnow])
        tryCatch({pval_vec[posnow] <- chisq.test(temp_tab)$p.value}, 
                 warning=function(w) print(posnow))
      }
      mlog = -log10(pval_vec)
      pos = as.vector(1:bp.use)
      pv_df = as.data.frame(cbind(pos,pval_vec, mlog))
      names(pv_df) = c("Position", "p-value", "negLog10")
      write.csv(pv_df, paste0("pvalue/pvalue_",gad_var_vec[i],"_vs_",gad_var_vec[k],".csv"), row.names = F)
      pval_vec_plot <- pval_vec
      pval_vec_plot[which(pval_vec_plot==0)] <- 10^(-300)
      plot(-log10(pval_vec_plot), pch = 19, xlab = "Position",
           ylab = "-log10(p-value)", 
           main = paste0(gad_var_vec[i],"_vs_",gad_var_vec[k]),
           ylim = c(0, 100))
    }
  }
}

##VERIFY##
#Look at the frequency graph and identify the region you want to extract according to p-value and frequency
sel_reg_st = 10
sel_reg_end = 50

##

##RUN##
#Consensus in the identified region
for(i in 1:n_gv){
  if(isTRUE(rev_vec[i])){
   res_ext = res_list[[i]]
   seq_set_sub =  DNAStringSet(res_ext$seqs_set) %>% reverseComplement()
   assign(paste0("sss_",gad_var_vec[i]), seq_set_sub)
   rm(seq_set_sub)
   rm(res_ext)
  }else{
    res_ext = res_list[[i]]
    seq_set_sub =  res_ext$seqs_set
    assign(paste0("sss_",gad_var_vec[i]), seq_set_sub)
    rm(seq_set_sub)
    rm(res_ext)
  }
}
#Make the list of seq_set_sub
sss_list = mget(paste0("sss_",gad_var_vec[1:n_gv] ))

for(k in 1:n_gv){
  sss = sss_list[[k]]
  ssb = sss
  for(i in 1:length(sss)){
    ssb[i] <- subseq(sss[i],sel_reg_st , sel_reg_end)
  }
  assign(paste0("ssb_",gad_var_vec[k]), ssb)
  rm(sss)
  rm(ssb)
}

#Make the list of ssb
ssb_list = mget(paste0("ssb_",gad_var_vec[1:n_gv] ))

#Plot the consensus motif of these regions
for(i in 1:n_gv){
  ssb = ssb_list[[i]]
  seqLogo( consensusMatrix(DNAStringSet(ssb) , as.prob = TRUE)[1:4,])
}

#Extract the matrching status
if (file.exists("Status_Table")){
} else {
  dir.create("Status_Table")
}
for(i in 1:n_gv){
  rr = rr_list[[i]]
  gr = gr_list[[i]]
  ssb = ssb_list[[i]]
  matchID = rr$seqID_vec
  df = as.data.frame(gr)
  names(gr) = c("chr", "start", "end")
  df$motif_like = 0
  df$motif_like[matchID] <- 1
  df$sequence = NA
  unq_seqID <- matchID %>% unique()
  use_id <- match( unq_seqID, matchID )
  df$sequence[unique(matchID)] <- ssb[use_id]
  assign(paste0("df_",gad_var_vec[i]), df)  
  write.csv(df, paste0("Status_Table/",gad_var_vec[i],".csv"), row.names = F)  
}



##

##VERIFY##
#Make sure that the core motif are align on all plots.
#If not, change the start and end in st_vec and end_vec and 
#rerun the whole program until every single one is aligned

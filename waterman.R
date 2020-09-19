install.packages("seqinr", repos="http://cran.us.r-project.org")
library(seqinr)

waterman <- function(fasta){
  # sequences
  fasta <- "test.fasta"
  sequences = read.fasta(file = fasta, seqtype = "DNA", forceDNAtolower = FALSE)
  seq1 <- c(NA,getSequence(sequences)[[1]])
  seq2 <- c(NA,getSequence(sequences)[[2]])
 
  # ---------------------------- build data frame
  aln1 <- matrix(, nrow = length(seq1), ncol = length(seq2))
  rownames(aln1) <- seq1
  colnames(aln1) <- seq2
  aln1.df <- as.data.frame(aln1)
  #  aln1.df
  
  # ---------------------------- copy data frame
  traceback.df <- as.data.frame(aln1.df)
  #  traceback.df
  
  # ---------------------------- assign variables
  gap <- -1
  match <- 2
  mismatch <- -2
  b <- 0
  r <- 0
  d <- 0 
  
  # ---------------------------- initialize
  aln1.df[1,] <- 0
  aln1.df[,1] <- 0
  
  # ---------------------------- rest of matrix
  for(i in 2:length(seq1)){
    for(j in 2:length(seq2)){
      b <- aln1.df[i-1, j] + gap 
      r <- aln1.df[i, j-1] + gap 
      if(seq1[i] == seq2[j]){
        d <- aln1.df[i-1, j-1] + match
      } else {
        d <- aln1.df[i-1, j-1] + mismatch
      }
      aln1.df[i,j] <- max(b,r,d,0)
    }
  }
  
  print(aln1.df)
  # create empty vectors for alignments
  seqfinal1 <- c()
  seqfinal2 <- c()
  
  # find max value in data frame 
  q <- which(aln1.df==max(aln1.df), arr.ind=T)
  s <- 1 # this sets index value for while loop
  
  while(s <= nrow(q)){
    seqfinal1 <- c()
    seqfinal2 <- c()
    i <- q[s,1]
    j <- q[s,2]
    while(i > 1 && j > 1){
      if(aln1.df[i,j] == 0){
        i <- 1
        j <- 1
        # this is written in to esentially stop the sequence from running when it hits a 0 in the alignment data frame
      } else if((aln1.df[i,j] == aln1.df[i-1,j-1] + match) && (seq1[i] == seq2[j])){ # checking to see if it's a match and moves diagonally
        seqfinal1 <- c(seq1[i], seqfinal1) 
        seqfinal2 <- c(seq2[j], seqfinal2)
        j <- j - 1
        i <- i - 1
      } else if(aln1.df[i,j] == aln1.df[i-1,j] + gap){ # 
        seqfinal1 <- c(seq1[i], seqfinal1)
        seqfinal2 <- c("-", seqfinal2)
        i <- i - 1
      } else if(aln1.df[i,j] == aln1.df[i,j-1] + gap){
        seqfinal1 <- c("-", seqfinal1)
        seqfinal2 <- c(seq2[j], seqfinal2)
        j <- j - 1
      } else if((aln1.df[i,j] == aln1.df[i-1,j-1] + mismatch) && (seq1[i] != seq2[j])){
        seqfinal1 <- c(seq1[i], seqfinal1) 
        seqfinal2 <- c(seq2[j], seqfinal2)
        j <- j - 1
        i <- i - 1
      }
      # how can I make this next part output 4 different fasta files instead of these 2
      # 
      write.fasta(sequences = list(seqfinal1,seqfinal2), names = names(sequences), file.out = paste0("aligned_",basename(fasta)))
      s <- s+1
    }
}
}

if(!interactive()) {
  ARGS = commandArgs(trailingOnly = TRUE)
  waterman(ARGS[1])
}
# install.packages("seqinr", repos="http://cran.us.r-project.org")
library(seqinr)

wunsch <- function(fasta){
  # sequences
  # fasta <- "test.fasta"
  sequences = read.fasta(file = fasta, seqtype = "DNA", forceDNAtolower = FALSE)
  seq1 <- c(NA,getSequence(sequences)[[1]])
  seq2 <- c(NA,getSequence(sequences)[[2]])
  
  #seq1 <- c(NA, strsplit(seq1.st, split = "", useBytes = TRUE)[[1]])
  #seq2 <- c(NA, strsplit(seq2.st, split = "", useBytes = TRUE)[[1]])
  
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
  
  
  # ---------------------------- int rows
  for(i in 1:length(seq2)){
    aln1.df[1,i] <- gap * (i-1)
    traceback.df[1,i] <- "R"
  }
  
  # ---------------------------- int columns
  for(i in 1:length(seq1)){
    aln1.df[i,1] <- gap * (i-1)
    traceback.df[i,1] <- "B"
  }
  
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
      aln1.df[i,j] <- max(b,r,d)
    }
  }
  
    print(aln1.df)
    # create empty vectors for alignments
    seqfinal1 <- c()
    seqfinal2 <- c()
    
    # redefine i and j
    i <- nrow(aln1.df)
    j <- ncol(aln1.df)
    
    while(i > 1 && j > 1){
      if((aln1.df[i,j] == aln1.df[i-1,j-1] + match) && (seq1[i] == seq2[j])){ # checking to see if it's a match and moves diagonally
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
    }
    
    # Save alignment in fasta format
    write.fasta(sequences = list(seqfinal1,seqfinal2), names = names(sequences), file.out = paste0("aligned_",basename(fasta)))
    # print(seqfinal1)
    # print(seqfinal2)
}

if(!interactive()) {
  ARGS = commandArgs(trailingOnly = TRUE)
  wunsch(ARGS[1])
}
needleman <- function(){
  
  # string to vector
  seq1.st <- "ACCGGTAT"
  seq2.st <- "ACCTATC"
  
  seq1 <- c(NA, strsplit(seq1.st, split = "", useBytes = TRUE)[[1]])
  seq2 <- c(NA, strsplit(seq2.st, split = "", useBytes = TRUE)[[1]])
  
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
#  aln1.df
#  traceback.df
  
  # ---------------------------- int columns
  for(i in 1:length(seq1)){
    aln1.df[i,1] <- gap * (i-1)
    traceback.df[i,1] <- "B"
  }
#  aln1.df
#  traceback.df
  
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
    m <- which.max(c(b,r,d))
    traceback.df[i,j] <- ifelse(m == 1, "B", ifelse(m == 2, "R", "D") )
    #traceback.df[i,j] <- ifelse(max(b,r,d) == b, "B", ifelse(max(b,r,d) == r, "R", "D") )
    }
  }
#  aln1.df
  return(traceback.df)
}
# ---------------------------- traceback

?which.max # explore this 

# Create df for values to generate into 
tb.aln <- matrix (nrow = 2, ncol = length(seq1) - 1)
final.df <- as.data.frame(tb.aln)
rownames(final.df) <- c("Sequence 1", "Sequence 2")
colnames(final.df) <- c(1:(length(seq1)-1))
final.df

# ----------------------------- traceback part 2 
# use while loop? 

backtrace = function(traceback.df){

  # empty alignments
  finalseq1 <- c()
  finalseq2 <- c()

  #reset i and j
  i <- nrow(traceback.df)
  j <- ncol(traceback.df)

  while(i > 1 && j > 1){
   if(traceback.df[i,j] == "D"){
      finalseq1 <- c(rownames(traceback.df)[i], finalseq1)
      finalseq2 <- c(colnames(traceback.df)[j], finalseq2)
    j <- j - 1
    i <- i - 1
    } else if(traceback.df[i,j] == "B"){
      finalseq1 <- c(rownames(traceback.df)[i], finalseq1)
      finalseq2 <- c("-", finalseq2)
      i <- i - 1
    } else if(traceback.df[i,j] == "R"){
      finalseq1 <- c("-", finalseq1)
      finalseq2 <- c(colnames(traceback.df)[j], finalseq2)
      j <- j - 1
    } 
  }
  finalalignment <- rbind(finalseq1, finalseq2)
  return(finalalignment)
}

alignment <- backtrace(traceback.df)
alignment

add <- function(a,b) {
  return(a+b)
}

result <- add(2,3)

traceback.df


# next steps
# figure out other traceback 
# run script from terminal 
# read sequences from FASTA 

# old traceback 
#for(i in length(seq1):2){
#  for(j in length(seq2):2){
#    x <- aln1.df[i,j]
#    y <- aln1.df[i+1, j]
#    z <- aln1.df[i, j+1]
#    if(x > y){
      
#    }
#  }}


############################## OLD CODE  - delete when final 

# SEQUENCES TO ALIGN 
# seq1 <- c(NA, "A", "C", "T", "G", "A", "T", "T", "C", "A")
# seq2 <- c(NA, "A", "C", "G", "C", "A", "T", "C", "A")

#  --  C  TC
#  AC   -  AT
  

  # for(i in aln1.df){
  #   if(seq1[i] == seq2[i]){
  #    x <- match
  #    } else {
  #      if(seq1[i] == seq2[i-1]){
#        x <- 1
#      } else { 
#        if(seq1[i-1] == seq2[i]){ 
#          x <- 1}
#        }
#    }
#  if(x == 0){
#    aln1.df[i+1,i+1] <- -2
#  }
#  x <- 0
#}

# for(i in aln1.df){
#   print("hello")
# }



# string to vector
seq1.st <- "ACCGGTAT"
seq2.st <- "ACCTATC"

seq1 <- c(NA, strsplit(seq1.st, split = "", useBytes = TRUE)[[1]])
seq2 <- c(NA, strsplit(seq2.st, split = "", useBytes = TRUE)[[1]])

# generating matrices 
P <- matrix(, nrow = length(seq1), ncol = length(seq2))
D <- matrix(, nrow = length(seq1), ncol = length(seq2))
Q <- matrix(, nrow = length(seq1), ncol = length(seq2))

for(i in P){
  
}
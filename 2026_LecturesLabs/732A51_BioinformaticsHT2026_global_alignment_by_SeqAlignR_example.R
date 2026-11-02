## This file accompanies the course 732A51 Bioinformatics 
## This code example is due to Leonard Persson Norblad (author of SeqAlignR)

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


library(SeqAlignR)
seq1 <- "GCATGCG"
seq2 <- "GATTACA"
example_alignment <- align_sequences(seq1, seq2, d = -1, mismatch = -1, match = 1, method="needleman")
print(example_alignment)
png("NW_seqalign_ex.png");plot(example_alignment);dev.off()


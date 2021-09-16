install.packages("seqinr")
library(seqinr)



getncbiseq1 <- function(accession)
{
  require("seqinr") # this function requires the SeqinR R package
  # first find which ACNUC database the accession is stored in:
  dbs <- c("genbank","refseq","refseqViruses","bacterial")
  numdbs <- length(dbs)
  for (i in 1:numdbs)
  {
    db <- dbs[i]
    choosebank(db)
    # check if the sequence is in ACNUC database 'db':
    resquery <- try(query(".tmpquery", paste("AC=", accession)), silent = TRUE)
    
    if (!(inherits(resquery, "try-error"))) {
      queryname <- "query2"
      thequery <- paste("AC=", accession, sep="")
      query2 <- query(queryname, thequery)
      # see if a sequence was retrieved:
      seq <- getSequence(query2$req[[1]])
      closebank()
      return(seq)
    }
    closebank()
  }
  print(paste("ERROR: accession",accession,"was not found"))
}
#Getting seqences from NCBI #####
dengueseq <- getncbiseq1("NC_001477")

dengueseq[1:24]
write.fasta(names = "Den_1", sequences = dengueseq, file.out = "den1.fasta")

#Reading in Fasta files from Computer####
dengue <- read.fasta(file = "/home/ubuntu/R/Seqdata/Dungue.fasta")  
head(dengue, 20)

#Length of a DNA sequence####
length(dengueseq)
#Base composition of a DNA sequence####
table(dengueseq)

#GC Content of DNA####
(2240+2770)*100/(3426+2240+2770+2299)
GC(dengueseq)




#DNA words (counting the different patterens ####
count(dengueseq,1)  #same as table(dengueseq)
count(dengueseq, 2)
count(dengueseq, 3)
count(dengueseq, 10)

denguetable <- count(dengueseq,1)
denguetable[[3]]  #this is telling you that there are 2770 g in this file
#or you can use
denguetable[["g"]]

tail(dengue)


#############DNA seq stats########################

x <- 100
log10(x)
myvector <- c(30,16,303,99,11,111)
mean(myvector)
myvector[3]
seq(1, 100, by = 1) # creating a seq of 1 - 100 by 1, you could do it by 2,3,4 ect by changing by = 
for (i in 1:10) { print (i*i) }
myvector1 <- c(10, 15, 22, 35, 43)
myvector2 <- c(3, 3.2, 3.9, 4.1, 5.2)
plot(myvector1, myvector2, xlab="myvector1", ylab="myvector2")
plot(myvector1, myvector2, xlab="myvector1", ylab="myvector2", type="b")  #type=b adds a  line

dengueseq[452:535]

###### A sliding window analysis of GC content
GC(dengueseq)  #whole file
GC(dengueseq[1:2000])   
GC(dengueseq[2001:4000])
GC(dengueseq[4001:6000])
GC(dengueseq[6001:8000])
GC(dengueseq[8001:10000])
GC(dengueseq[10001:10735])


#####using a loop to get GC as above
starts <- seq(1, length(dengueseq)-2000, by = 2000)
starts
n <- length(starts)  
for (i in 1:n) {
  chunk <- dengueseq[starts[i]:(starts[i]+1999)]
  chunkGC <- GC(chunk)
  print (chunkGC)
}

###plotting gc

chunkGCs <- numeric(n)
for (i in 1:n) {
  chunk <- dengueseq[starts[i]:(starts[i]+1999)]
  chunkGC <- GC(chunk)
  print(chunkGC)
  chunkGCs[i] <- chunkGC
}
plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content")

######using a function to help with GC sliding window plot
slidingwindowplot <- function(windowsize, inputseq)
{
  starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
  n <- length(starts)    # Find the length of the vector "starts"
  chunkGCs <- numeric(n) # Make a vector of the same length as vector "starts", but just containing zeroes
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
    chunkGC <- GC(chunk)
    print(chunkGC)
    chunkGCs[i] <- chunkGC
  }
  plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content")
}
  

slidingwindowplot(3000, dengueseq)  #using number like in by= 3000
slidingwindowplot(300, dengueseq)    #by = 300

##############Over-represented and under-represented DNA words #########

count(dengueseq, 1) # Get the number of occurrences of 1-nucleotide DNA words
2770/(3426+2240+2770+2299) # Get fG  ρ(GC) = fGC/(fG * fC)
2240/(3426+2240+2770+2299) # Get fC   ρ(GC) = fGC/(fG * fC)
count(dengueseq, 2) # Get the number of occurrences of 2-nucleotide DNA words
500/(1108+720+890+708+901+523+261+555+976+500+787+507+440+497+832+529) # Get fGC
0.04658096/(0.2580345*0.2086633) # Get rho(GC)


####################Sequence Databases#######################

#mysequence1
#ACATGAGACAGACAGACCCCCAGAGACAGACCCCTAGACACAGAGAGAG
#TATGCAGGACAGGGTTTTTGCCCAGGGTGGCAGTATG

# List all the sub-databases in ACNUC
choosebank()
choosebank("genbank")  # Specify that we want to search the 'genbank' ACNUC sub-database
choosebank("refseq") # Specify that we want to search the 'refseq' ACNUC sub-database

query("RefSeqBact", "SP=Bacteria")
closebank()

choosebank("genbank")
query("SchistosomamRNA", "SP=Schistosoma mansoni AND M=mrna")
closebank()
choosebank("refseqViruses")
Dengue1 <- query("Dengue1", "AC=NC_001477")
help("query")
attributes(Dengue1)
Dengue1$nelem
attr(Dengue1, "names")
dengueseq <- getSequence(Dengue1$req[[1]])  # “getSequence()” to retrieve the sequence data for the DEN-1 Dengue virus genome, and puts the sequence into a variable dengueseq:
annots <- getAnnot(Dengue1$req[[1]])  # you can retrieved its annotations by using the “getAnnot()” function.
annots[1:20]
closebank()



#finding the sequences published in Nature 460:352-358###

choosebank("genbank") # Specify that we want to search the 'genbank' ACNUC sub-database
naturepaper <- query('naturepaper', 'R=Nature/460/352')
naturepaper$nelem






########################Pairwise Sequence Alignment########################################

choosebank()
leprae <- read.fasta(file = "/home/ubuntu/R/Seqdata/Q9CD83.fasta")
lepraeseq <- leprae[[1]]
ulcerans <- read.fasta(file = "/home/ubuntu/R/Seqdata/A0PQ23.fasta")
ulceransseq <- ulcerans[[1]]
ulceransseq

### dotplot of the sequences for the chorismate lyase proteins 
dotPlot(lepraeseq, ulceransseq)  #the M. leprae sequence is plotted along the x-axis; M. ulcerans sequence is plotted along the y-axis:
help(dotplot)


##################Pairwise global alignment of DNA sequences using the Needleman-Wunsch algorithm#######
library(Biostrings)
sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
sigma
help(Biostrings)


s1 <- "GAATTC"
s2 <- "GATTA"
globalAligns1s2 <- pairwiseAlignment(s1, s2, substitutionMatrix = sigma, gapOpening = -2,
                                     gapExtension = -8, scoreOnly = FALSE)
globalAligns1s2 # Print out the optimal alignment and its score



#####Pairwise global alignment of protein sequences using the Needleman-Wunsch algorithm#############

data("BLOSUM50")
BLOSUM50
data(package="Biostrings")
s3 <- "PAWHEAE"

s4 <- "HEAGAWGHEE"
globalAligns3s4 <- pairwiseAlignment(s3, s4, substitutionMatrix = "BLOSUM50", gapOpening = -2,gapExtension = -8, scoreOnly = FALSE)
globalAligns3s4 # Print out the optimal global alignment and its score

#Aligning UniProt sequences###########

lepraeseqstring <- c2s(lepraeseq)  # Make a string that contains the sequence in "lepraeseq"  convert the vectors lepraeeq and ulceransseq into strings using c2s.
ulceransseqstring <- c2s(ulceransseq) # Make a string that contains the sequence in "ulceransseq"
lepraeseq
lepraeseqstring

lepraeseqstring <- toupper(lepraeseqstring)   ##convert from lower case to upper case
ulceransseqstring <- toupper(ulceransseqstring)

globalAlignLepraeUlcerans <- pairwiseAlignment(lepraeseqstring, ulceransseqstring,substitutionMatrix = BLOSUM50, gapOpening = -2, gapExtension = -8, scoreOnly = FALSE)
globalAlignLepraeUlcerans # Print out the optimal global alignment and its score



##########Viewing a long pairwise alignment###################

printPairwiseAlignment <- function(alignment, chunksize=60, returnlist=FALSE)# function “printPairwiseAlignment()” below will do this for you:
{
  require(Biostrings)           # This function requires the Biostrings package
  seq1aln <- pattern(alignment) # Get the alignment for the first sequence
  seq2aln <- subject(alignment) # Get the alignment for the second sequence
  alnlen  <- nchar(seq1aln)     # Find the number of columns in the alignment
  starts  <- seq(1, alnlen, by=chunksize)
  n       <- length(starts)
  seq1alnresidues <- 0
  seq2alnresidues <- 0
  for (i in 1:n) {
    chunkseq1aln <- substring(seq1aln, starts[i], starts[i]+chunksize-1)
    chunkseq2aln <- substring(seq2aln, starts[i], starts[i]+chunksize-1)
    # Find out how many gaps there are in chunkseq1aln:
    gaps1 <- countPattern("-",chunkseq1aln) # countPattern() is from Biostrings package
    # Find out how many gaps there are in chunkseq2aln:
    gaps2 <- countPattern("-",chunkseq2aln) # countPattern() is from Biostrings package
    # Calculate how many residues of the first sequence we have printed so far in the alignment:
    seq1alnresidues <- seq1alnresidues + chunksize - gaps1
    # Calculate how many residues of the second sequence we have printed so far in the alignment:
    seq2alnresidues <- seq2alnresidues + chunksize - gaps2
    if (returnlist == 'FALSE')
    {
      print(paste(chunkseq1aln,seq1alnresidues))
      print(paste(chunkseq2aln,seq2alnresidues))
      print(paste(' '))
    }
  }
  if (returnlist == 'TRUE')
  {
    vector1 <- s2c(substring(seq1aln, 1, nchar(seq1aln)))
    vector2 <- s2c(substring(seq2aln, 1, nchar(seq2aln)))
    mylist <- list(vector1, vector2)
    return(mylist)
  }
}



printPairwiseAlignment(globalAlignLepraeUlcerans, 60)


########Pairwise local alignment of protein sequences using the Smith-Waterman algorithm#######


localAlignLepraeUlcerans <- pairwiseAlignment(lepraeseqstring, ulceransseqstring,
                                              substitutionMatrix = BLOSUM50, gapOpening = -2, gapExtension = -8, scoreOnly = FALSE, type="local")

localAlignLepraeUlcerans # Print out the optimal local alignment and its score

printPairwiseAlignment(localAlignLepraeUlcerans, 60)

################Calculating the statistical significance of a pairwise global alignment#########

generateSeqsWithMultinomialModel <- function(inputsequence, X)
{
  # Change the input sequence into a vector of letters
  require("seqinr") # This function requires the SeqinR package.
  inputsequencevector <- s2c(inputsequence)
  # Find the frequencies of the letters in the input sequence "inputsequencevector":
  mylength <- length(inputsequencevector)
  mytable <- table(inputsequencevector)
  # Find the names of the letters in the sequence
  letters <- rownames(mytable)
  numletters <- length(letters)
  probabilities <- numeric() # Make a vector to store the probabilities of letters
  for (i in 1:numletters)
  {
    letter <- letters[i]
    count <- mytable[[i]]
    probabilities[i] <- count/mylength
  }
  # Make X random sequences using the multinomial model with probabilities "probabilities"
  seqs <- numeric(X)
  for (j in 1:X)
  {
    seq <- sample(letters, mylength, rep=TRUE, prob=probabilities) # Sample with replacement
    seq <- c2s(seq)
    seqs[j] <- seq
  }
  # Return the vector of random sequences
  return(seqs)
}

##We can use this function to generate 1000 7-letter amino acid sequences using a multinomial model in which the probabilities of the letters are set equal to their frequencies in ‘PAWHEAE’ ##


randomseqs <- generateSeqsWithMultinomialModel('PAWHEAE',1000)
randomseqs[1:10] # Print out the first 10 random sequences

s4 <- "HEAGAWGHEE"
pairwiseAlignment(s4, randomseqs[1], substitutionMatrix = "BLOSUM50", gapOpening = -2,
                  gapExtension = -8, scoreOnly = FALSE)
pairwiseAlignment(s4, randomseqs[1], substitutionMatrix = "BLOSUM50", gapOpening = -2,
                  gapExtension = -8, scoreOnly = TRUE)

#We can then compare the actual score for aligning ‘PAWHEAE’ to ‘HEAGAWGHEE’ (ie. -5) to the distribution of scores for aligning ‘HEAGAWGHEE’ to the random sequences.


randomscores <- double(1000) # Create a numeric vector with 1000 elements
for (i in 1:1000)
{
  score <- pairwiseAlignment(s4, randomseqs[i], substitutionMatrix = "BLOSUM50",
                             gapOpening = -2, gapExtension = -8, scoreOnly = TRUE)
  randomscores[i] <- score
}



hist(randomscores, col="red") # Draw a red histogram
sum(randomscores >= -5)



########Multiple Alignment and Phylogenetic trees############

#retrieve several sequences from UniProt
retrieveseqs <- function(seqnames,acnucdb)
{
  myseqs <- list()   # Make a list to store the sequences
  require("seqinr")  # This function requires the SeqinR R package
  choosebank(acnucdb)
  for (i in 1:length(seqnames))
  {
    seqname <- seqnames[i]
    print(paste("Retrieving sequence",seqname,"..."))
    queryname <- "query2"
    query <- paste("AC=",seqname,sep="")
    query2 <- query(`queryname`,`query`)  ###had to fix this i used the code in getncbiseq1
    seq <- getSequence(query2$req[[1]]) # Makes a vector "seq" containing the sequence
    myseqs[[i]] <- seq
  }
  closebank()
  return(myseqs)
}

closebank()
seqnames <- c("P06747", "P0C569", "O56773", "Q5VKP1")  # Make a vector containing the names of the sequences
seqs <- retrieveseqs(seqnames,"swissprot")             # Retrieve the sequences and store them in list variable "seqs"
length(seqs)    # Print out the number of sequences retrieved
seq1 <- seqs[[1]]    # Get the first sequence
seq1[1:20]     # Print out the first 20 letters of the first sequence
seq2 <- seqs[[2]]   # Get the second sequence
seq3 <-seqs[[3]]
seq4 <- seqs[[4]]
seq4[1:20]   
write.fasta(seqs, seqnames, file="/home/ubuntu/R/Seqdata/phosphoproteins.fasta")

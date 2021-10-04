library(seqinr)
library(msa)  ##this library is the alignment software
library(ape)  # this aids in creating a tree
install.packages("ips")
install.packages("reactable")

## make sure to set the working directory to the folder with all the sequences
setwd("~/programs/R/merging_fasta/outbreak1")
############# works by reading in single files but you have to create each vari###############
fa1 = read.fasta(file = "gisaid_hcov-19_2021_08_10_20.fasta", as.string = TRUE)
fa2 = read.fasta(file = "gisaid_hcov-19_2021_08_10_20_a.fasta", as.string = TRUE)


##################### TRYING TO AUTOMATE IT WITH LIST OF FILES ############
files <- list.files(pattern = "\\.fasta$")

DF <-  read.fasta(files[1], as.string = TRUE)

for (f in files[-1]){
  df <- read.fasta((f), as.string = TRUE)      # read the file
  DF <- rbind(DF, df)    # append the current file
  write.fasta(DF, names = getName(DF), file.out = 'merged.fasta')
}

## Reading the document back into R and Checking the names of the sequences
merged_seq = read.fasta(file = "merged.fasta")
getName(merged_seq)


## this step preps the file with all the sequences for the alignment

my_seq = readDNAStringSet("merged.fasta")

#this step creates the alignment.  This step with take a bit of time.
# the step is using the default setting within MSA (clustW)
COValign = msa(my_seq)  

#####################TEST _ DELETE NEEDING TO FIGURE OUT HOW TO GET PRETTY PRINT TO WORK #####################################
msaPrettyPrint(x=COValign, file="COVIDAlignOutbreak.pdf",
               paperWidth=15.3, paperHeight=4.5,
               shadingMode="functional",
               shadingModeArg="structure",
               shadingColors="greens",
               logoColors="rasmol",
               showLogoScale="right",
               showLegend=FALSE,
               askForOverwrite=FALSE)


Sys.which("pdflatex")
Sys.getenv("PATH")
Sys.setenv(PATH=paste(Sys.getenv("PATH"),"/usr/texbin",sep=":"))

### this produces a file output of your aligned sequences #

alignment2Fasta <- function(alignment, filename) {
  sink(filename)
  
  n <- length(rownames(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    the.sequence <- toString(unmasked(alignment)[[i]])
    cat(the.sequence)
    cat('\n')  C
  }
  
  sink(NULL)
}


alignment2Fasta(COValign, 'outalign.fasta')  


###### Creating trees by hand ####
#this aids in prepping the document for creating the tree (creates the matrix)
COValign2 = msaConvert(COValign, type = "seqinr::alignment")
?msaConvert
#dist.alignment(COValign2, "identity")
d = as.matrix(dist.alignment(COValign2, "identity"))
View(d)
COVIDalign3 = as.DNAbin(COValign)
COVIDalign_dist = dist.dna(COVIDalign3)
COVIDalign_dist









####### DIFFERENT TREES####################################################
### this is a neighbor joining tree ####
COVtree = read.tree(d)
COVtree_upgma = 
COVtree_nj = nj(d)
COVtree2 = nj(COVIDalign_dist)

tree.Ult.MPL <- chronoMPL(COVtree_nj)
tree.Ult.S <- chronos(COVtree_nj)
plot(tree.Ult.MPL, main="Smaller Outbreak Neighbor Joining Tree", edge.width = 2,
     cex = 0.5, cex.main = 2); foo()
plot(tree.Ult.S, main="Smaller Outbreak Neighbor Joining Tree", edge.width = 2,
     cex = 0.5, cex.main = 2); foo()

plot(COVtree_nj, main="Smaller Outbreak Neighbor Joining Tree", edge.width = 2,
     cex = 0.5, cex.main = 2); foo()
plot(COVtree2, main="Smaller Outbreak Neighbor Joining Tree", edge.width = 2,
     cex = 0.5, cex.main = 2); foo()
plot(COVtree_nj, "u", main="Smaller Outbreak unrooted Tree", edge.width = 2,
     cex = 2, cex.main = 2); foo()

plot(COVtree_nj, "c", main="Smaller Outbreak Cladogram Tree"); foo()
plot(COVtree_nj, "f", main="Smaller Outbreak Cladogram Tree"); foo()
box("outer")
add.scale.bar(cex = 0.7, font = 1, col = "red")
write.tree(COVtree_nj, file = "COVTREE.nwk")


### Minimum Spanning Tree #####
COVtree_mst = mst(d)
plot(COVtree_mst, main = "Smaller outbreak Minimum Spanning Tree")
write.tree(COVtree_mst, file = "COVTREE_MST.nwk")




###############  playing around ######
add.scale.bar(cex = 0.7, font = 2, col = "red")
layout(1)
plot(COVtree_nj, main="Outbreak #2", type = "unrooted")
plot(COVtree, main="Outbreak #2", type = "fan")
plot(COVtree, main="Outbreak #2", type = "cladogram")
COVtree$tip.label
COVtree$edge
COVtree$Nnode
####
#function to make the tree pretty
mytr <- read.tree(text = "((Pan:5, Homo:5'):2, Gorilla:7);")
foo <- function() {
  col <- "green"
  for (i in 1:2)
    axis(i, col = col, col.ticks = col, col.axis = col, las = 1)
  box(lty = "19")
}


############# other code to try ####
#https://rpubs.com/vuongyenxuan/project
COValign$seq



#########################Function to produce a rooted tree ######
#https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter5.html

rootedNJtree <- function(alignment, theoutgroup, type)
{
  # load the ape and seqinR packages:
  require("ape")
  require("seqinr")
  # define a function for making a tree:
  makemytree <- function(alignmentmat, outgroup=`theoutgroup`)
  {
    alignment <- ape::as.alignment(alignmentmat)
    if      (type == "protein")
    {
      mydist <- dist.alignment(alignment)
    }
    else if (type == "DNA")
    {
      alignmentbin <- as.DNAbin(alignment)
      mydist <- dist.dna(alignmentbin)
    }
    mytree <- nj(mydist)
    mytree <- makeLabel(mytree, space="") # get rid of spaces in tip names.
    myrootedtree <- root(mytree, outgroup, r=TRUE)
    return(myrootedtree)
  }
  # infer a tree
  mymat  <- as.matrix.alignment(alignment)
  myrootedtree <- makemytree(mymat, outgroup=theoutgroup)
  # bootstrap the tree
  myboot <- boot.phylo(myrootedtree, mymat, makemytree)
  # plot the tree:
  plot.phylo(myrootedtree, type="p")  # plot the rooted phylogenetic tree
  nodelabels(myboot,cex=0.7)          # plot the bootstrap values
  mytree$node.label <- myboot   # make the bootstrap values be the node labels
  return(mytree)
}

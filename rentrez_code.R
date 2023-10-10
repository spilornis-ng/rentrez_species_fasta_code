library(rentrez)

setwd("your/working/directory") #set you working directory where you want to download the sequences

sp<-read.csv("species_list.csv", header = T) #Read in a species list

g<-"ND2" #assign the gene you would like to download (check NCBI for most common gene symbols) - for example ND2

for (k in 1:nrow(sp)){
  f<-sp$SCINAME[k]
  i<-paste0(f, "[ORGN]", " ", "AND", " ", "ND2")  #Here we are defining the search term - Species name, and the gene
  x<-entrez_search(db = "nucleotide", term = i) #This performs the search
  l<-length(x$ids) #here we see how many entries are there for each species for a given gene
  if(l==0) next #if there is no entry it will skip to the next species
  if(l>10){ #if there are more than 10 entries then I have limited myself to only first 10 search hits
    l<-10
  } 
  for (j in 1:l){
    q<-entrez_summary(db = "nucleotide", id = x$ids[j]) #We gather the summary for the Genbank IDs that we got in the search results
    if(q$moltype!="dna") next #sometimes it will also document rna, which I was not interested in so i skip the entry if its not dna
    if(q$slen>2000) next #most of the genes that we need are going to be shorter than 2000bp
    #sometimes a mitogenome will be one of the first hits, we do not want to downlaod a full mitogenome, so we avoid by skipping anything longer than 2000bp
    p<-entrez_fetch(db = "nucleotide", id = x$ids[j], rettype = "fasta") #here we fetch the fasta files for the filtered Genbank IDs
    r<-paste0(f,"_",g, "_", j, ".fasta") #now you can define a name - for example here its going to be Species_ND2_genbankID.fasta
    try(write(p, file = r)) #finally you write the output
  }
}
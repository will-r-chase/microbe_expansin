library(readxl)
library(plyr)

# Function read fasta file into dataframe
fasta2df<-function(file) {
  fasta<-readLines(file)
  ind<-grep(">", fasta)
  s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasta)))
  seqs<-rep(NA, length(ind))
  for(i in 1:length(ind)) {
    seqs[i]<-paste(fasta[s$from[i]:s$to[i]], collapse="")
  }
  DF<-data.frame(name=gsub(">", "", fasta[ind]), sequence=seqs, stringsAsFactors = F)
  return(DF)
}

# Function write data frame as fasta file
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"sequence"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# Function extract sequence for domain based on HMMER CAZYzyme annotation
# Input the domain name you want to extract, the sequence file with fusion proteins and the HMMER results in excel file
extractDomain<-function(domain.name, sequence.file, dbCan.file){
  fusion_seqs<-fasta2df(sequence.file)
  for(i in 1:length(fusion_seqs$name)){
    ID<-strsplit(fusion_seqs$name[i], "\\.")[[1]][1]
    fusion_seqs$id[i]<-ID
  }
  
  CAZY <- read_excel(dbCan.file)
  
  for(i in 1:length(CAZY$`Gene ID`)){
    ID2<-strsplit(CAZY$`Gene ID`[i], "\\.")[[1]][1]
    CAZY$id[i]<-ID2
  }
  
  CAZY<-data.frame(id=CAZY$id, domains=CAZY$HMMER, stringsAsFactors = F)
  fusion_dat<-join(fusion_seqs, CAZY, by="id")
  fusion_dat2<-fusion_dat[grep(domain.name, fusion_dat$domains),]
  
  for(i in 1:length(fusion_dat2$name)){
    domains_split<-unlist(strsplit(fusion_dat2$domains[i], "\\+"))
    fusion_dat2$domain.name[i]<-domains_split[grep(domain.name, domains_split)]
  }
  
  fusion_dat2$coords<-gsub(".*\\((.*)\\).*", "\\1", fusion_dat2$domain.name)
  
  for(i in 1:length(fusion_dat2$name)){
    coords_split<-strsplit(fusion_dat2$coords[i], "\\-")
    fusion_dat2$start[i]<-as.numeric(coords_split[[1]][1])
  }
  for(i in 1:length(fusion_dat2$name)){
    coords_split<-strsplit(fusion_dat2$coords[i], "\\-")
    fusion_dat2$end[i]<-as.numeric(coords_split[[1]][2])
  }
  
  fusion_dat2$start<-as.numeric(fusion_dat2$start)
  fusion_dat2$end<-as.numeric(fusion_dat2$end)
  
  fusion_dat2$domain.name_seq<-substr(fusion_dat2$sequence, start=fusion_dat2$start, stop=fusion_dat2$end)
  
  domain.name_seqs<-data.frame(name=fusion_dat2$name, sequence=fusion_dat2$domain.name_seq, stringsAsFactors = F)
  writeFasta(domain.name_seqs, paste(domain.name, ".txt", sep=""))
  
}
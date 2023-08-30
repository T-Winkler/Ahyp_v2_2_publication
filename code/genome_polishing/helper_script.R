# set working directory
setwd("/home/tom/Documents/projects/Ahyp_v2_2_publication/")

# read in list of all headers with the correct order
headers <- read.table("data/NextPolish/input/out.headers.txt")
headers <- gsub(">","",headers$V1)
headers <- gsub("quiver_","quiver",headers)

# read in prefiltered fasta index
prefilter <- read.table("data/NextPolish/input/out.prefiltered.renamed.txt.fai")

# use the 
# every sequence that is in 
no_seq <- headers[!headers %in% prefilter[,1]]
no_seq <- sub("",">",no_seq)

write.table(no_seq, file="data/NextPolish/processed/header_without_sequence.fa", 
            quote=F,
            row.names = F,
            col.names = F)

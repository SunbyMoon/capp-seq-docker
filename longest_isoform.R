## get the longest isoform per gene
# code from https://gist.github.com/jwaageSnippets/5133941
session <- browserSession("UCSC")
genome(session)<-"hg19"
query <- ucscTableQuery(session, "refGene")
tableName(query) <- "refGene"
getTable(query) -> refseq
refseq[,c(2,3,4,5,6,7,8,13)] -> refseq
refseq$"width" <- refseq$"txEnd"-refseq$"txStart"

as.character(refseq[,1]) -> refseq[,1]
as.character(refseq[,8]) -> refseq[,8]

split(refseq, refseq$"name2") -> refseqSplit

getGene <- function(x)
{
  x[which.max(x$"width"),]
}

lapply(refseqSplit, FUN=getGene) -> result
do.call("rbind", result) -> result

write.table(result, 'longest_resfseq.txt', row.names=F)

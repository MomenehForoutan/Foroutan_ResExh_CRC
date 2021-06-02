sigPath <- "../data/Signatures_literature/"

immuneStroma<- read.table(paste0(sigPath, "GeneSig.Stromal.Immune.txt"),
                           header = T, sep = "\t")

immune <- immuneStroma$Gene[immuneStroma$Set == "Immune141_UP"]
stroma <- immuneStroma$Gene[immuneStroma$Set == "Stromal141_UP"]

cib <- lapply(Sys.glob(paste0(sigPath, "CIBERSORT/*.tsv")), read.delim)
names(cib) <- list.files(paste0(sigPath, "CIBERSORT/"))

names(cib) <- gsub(".tsv", "", names(cib))

cib <- lapply(cib, function(x) {
  x <- x$X[x$Dir == "Up"]

})


cibImmStr <- c(cib, 
               list(Immune = immune), 
               list(Stroma = stroma))

ImmSigGeneSet <- lapply(cibImmStr, GSEABase::GeneSet)

ImmSigsGeneSetName <- lapply(names(ImmSigGeneSet), function(x) {
  currentSig <- ImmSigGeneSet[[x]]
  GSEABase::setName(currentSig) <- x
  return(currentSig)
})

names(ImmSigsGeneSetName) <- names(ImmSigGeneSet)

ImmSigsCollection <- GSEABase::GeneSetCollection(ImmSigsGeneSetName)

cibeNames <- names(ImmSigsCollection)



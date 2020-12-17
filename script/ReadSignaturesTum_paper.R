##--------- Read Marker signatures

library(singscore)

outPathCD8 <- "/Users/mfor0011/Documents/projects/Monash/ResExh_CRC/output/CD8/"
outPathCD4 <- "/Users/mfor0011/Documents/projects/Monash/ResExh_CRC/output/CD4/"
outPathNK <- "/Users/mfor0011/Documents/projects/Monash/ResExh_CRC/output/NK/"
sigPath <- "/Users/mfor0011/Documents/data/signatures/"
KalliesPath <- "/Users/mfor0011/Documents/projects/Monash/ResExh_CRC/output/Zhang/"
CollinsPath <- "/Users/mfor0011/Documents/projects/Monash/NK_Bulk/output/Collins/"
tcgaPath0 <- "/Users/mfor0011/Documents/projects/Monash/ResExh_CRC/output/TCGA/"


cd8Sigs <-
  read.table(
    paste0(outPathCD8, "TrmExhMarkers_CD8_CLs_LCM.txt"),
    header = T,
    sep = "\t",
    stringsAsFactors = F,
    check.names = F
  )


cd4Sigs <-
  read.table(
    paste0(outPathCD4, "TrmExhMarkers_CD4_CLs_LCM.txt"),
    header = T,
    sep = "\t",
    stringsAsFactors = F,
    check.names = F
  )


nkSigs <-
  read.table(
    paste0(outPathNK, "TrmExhMarkers_NK_CLs_LCM.txt"),
    header = T,
    sep = "\t",
    stringsAsFactors = F,
    check.names = F
  )


nkSigs$Cell.Type <- "NK"
cd8Sigs$Cell.Type <- "CD8"
cd4Sigs$Cell.Type <- "CD4"

cd8Sigs$Marker <- gsub("Trm", "Res", cd8Sigs$Marker)
cd4Sigs$Marker <- gsub("Trm", "Res", cd4Sigs$Marker)
nkSigs$Marker <- gsub("Trm", "Res", nkSigs$Marker)

allMarkers <- rbind(nkSigs, cd8Sigs, cd4Sigs)

exhData <- allMarkers[allMarkers$Marker == "Exh", ]
trmData <- allMarkers[allMarkers$Marker == "Res", ]

exhData_CanPass <- exhData[exhData$CancerPass, ]
trmData_CanPass <- trmData[trmData$CancerPass, ]




All_Exh <- unique(exhData$gene[! exhData$gene %in% c("IL12RB2", "MARCH3")])
All_Trm <- unique(trmData$gene[! trmData$gene %in% c("IL12RB2", "MARCH3")])

All_Exh_CancerPass <-  unique(exhData_CanPass$gene[! exhData_CanPass$gene %in% c("IL12RB2")])
All_Trm_CancerPass <-  unique(trmData_CanPass$gene[! trmData_CanPass$gene %in% c("IL12RB2")])

CD8NK_Exh <- unique(exhData_CanPass$gene[exhData_CanPass$Cell.Type  %in% c("NK", "CD8")])
CD8NK_Trm <- unique(trmData_CanPass$gene[trmData_CanPass$Cell.Type  %in% c("NK", "CD8")])


##------------ Marker genes
NK_ExhMarker <- nkSigs$gene[nkSigs$Marker == "Exh"]
NK_TrmMarker <- nkSigs$gene[nkSigs$Marker == "Res"]

CD4_ExhMarker <- cd4Sigs$gene[cd4Sigs$Marker == "Exh"]
CD4_TrmMarker <- cd4Sigs$gene[cd4Sigs$Marker == "Res"]

CD8_ExhMarker <- cd8Sigs$gene[cd8Sigs$Marker == "Exh"]
CD8_TrmMarker <- cd8Sigs$gene[cd8Sigs$Marker == "Res"]


NK_ExhMarkerCanPass <- nkSigs$gene[nkSigs$Marker == "Exh" & nkSigs$CancerPass]
NK_TrmMarkerCanPass <- nkSigs$gene[nkSigs$Marker == "Res" & nkSigs$CancerPass]

CD4_ExhMarkerCanPass <- cd4Sigs$gene[cd4Sigs$Marker == "Exh" & cd4Sigs$CancerPass]
CD4_TrmMarkerCanPass <- cd4Sigs$gene[cd4Sigs$Marker == "Res" & cd4Sigs$CancerPass]

CD8_ExhMarkerCanPass <- cd8Sigs$gene[cd8Sigs$Marker == "Exh" & cd8Sigs$CancerPass]
CD8_TrmMarkerCanPass <- cd8Sigs$gene[cd8Sigs$Marker == "Res" & cd8Sigs$CancerPass]



#--------- Read all other signatures ------------

###----------------- Sigs from immune and tum
sigsImmune <- readRDS(paste0(sigPath, "ProcessedSigs/SigLists_Immune.RDS"))
sigsTum <- readRDS(paste0(sigPath, "ProcessedSigs/SigLists_Tumour.RDS"))


##---- some cancer related signatures:
# Suppl table 5: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2816644/#sup1
# 6 genes overlapping between the two Hypoxia signatures
# [1] "BNIP3"  "SLC2A1" "PGK1"   "NDRG1"  "P4HA1" 
# [6] "ADM"  

# write.table(data.frame(Genes = Hypoxia_Buffa), "/Users/mfor0011/Documents/data/signatures/HPX/Hypoxia_Buffa.txt", sep = "\t", row.names = F)

Hypoxia_Buffa <- as.character(quote(c(GPI,
                                      UTP11L,
                                      ENO1,
                                      P4HA1,
                                      PFKP,
                                      CA9,
                                      PGAM1,
                                      NP,
                                      HIG2,
                                      C20orf20,
                                      AK3L1,
                                      ALDOA,
                                      ADM,
                                      SEC61G,
                                      ACOT7,
                                      PSMA7,
                                      LDHA,
                                      HK2,
                                      NDRG1,
                                      TPI1,
                                      SLC2A1,
                                      MRPL15,
                                      SLC25A32,
                                      TUBB6,
                                      DDIT4,
                                      CDKN3,
                                      VEGFA,
                                      MRPS17,
                                      PGK1,
                                      BNIP3,
                                      CORO1C,
                                      ANKRD37,
                                      MAP7D1,
                                      MIF,
                                      MCTS1,
                                      MAD2L2,
                                      MRPL13,
                                      SHCBP1,
                                      GAPDH,
                                      SLC16A1,
                                      YKT6,
                                      RBM35A,
                                      KIF20A,
                                      TUBA1B,
                                      TUBA1C,
                                      CHCHD2,
                                      ANLN,
                                      PSRC1,
                                      KIF4A,
                                      CTSL2,
                                      LRRC42)))[-1]



## curated list of antigen presentation genes
AgPresentation <-
  c("IFNGR1",
    "IFNGR2",
    "TAP1",
    "TAP2",
    "TAPBP",
    "JAK1",
    "JAK2",
    "STAT1",
    "CASP8",
    "B2M")  

## also read a few other signatures from MSigDB according to this paper:
# https://www-nature-com.ezproxy.lib.monash.edu.au/articles/s41590-020-0725-2#Sec46

hum <- msigdf::msigdf.human

checkSigName <- c(
  "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_PEROXISOME",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_HYPOXIA", 
  "HALLMARK_ANGIOGENESIS",
  "PID_NFAT_TFPATHWAY",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)

checkSigID <- hum[hum$geneset %in% checkSigName, ]
# glyco <- as.character(hum$entrez[hum$geneset == "HALLMARK_GLYCOLYSIS"])

library(org.Hs.eg.db)
checkSigID$Symbol <- annotate::getSYMBOL(as.character(checkSigID$entrez), data='org.Hs.eg')



selImmuneSigs <- c(
  list(

    ## add some tum signatures
    TGFbEMT = sigsTum$TGFb_EMT$tgfb_Up,
    Epi_Thiery = sigsTum$EpiT,
    Epi_Byers = sigsTum$EpiB,
    Mes_Thiery = sigsTum$MesT,
    Mes_Byers = sigsTum$MesB,
    HPX_EGF_EMT = sigsTum$HpxEGF_EMT$hpxEGF_Up,
    
    MYC_neg = sigsTum$MYC_negative,
    MYC_pos = sigsTum$MYC_positive,
    TP53_neg = sigsTum$TP53_negative,
    TP53_pos = sigsTum$TP53_positive,
    DDR = sigsTum$DDR_All,
    
    Antigen_Presentation = AgPresentation,
    Hypoxia20 = sigsTum$Hypoxia,
    Hypoxia_Buffa = Hypoxia_Buffa,
    Hypoxia = checkSigID$Symbol[checkSigID$geneset == "HALLMARK_HYPOXIA"],
    ReactiveOx = checkSigID$Symbol[checkSigID$geneset == "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY"],
    OXPHOS = checkSigID$Symbol[checkSigID$geneset == "HALLMARK_OXIDATIVE_PHOSPHORYLATION"],
    Peroxisome = checkSigID$Symbol[checkSigID$geneset == "HALLMARK_PEROXISOME"],
    KRAS_signaling = checkSigID$Symbol[checkSigID$geneset == "HALLMARK_KRAS_SIGNALING_UP"],
    TP53_signaling = checkSigID$Symbol[checkSigID$geneset == "HALLMARK_P53_PATHWAY"],
    WNT_signaling = checkSigID$Symbol[checkSigID$geneset == "HALLMARK_WNT_BETA_CATENIN_SIGNALING"],
    TGFb_Signaling = checkSigID$Symbol[checkSigID$geneset == "HALLMARK_TGF_BETA_SIGNALING"],
    EMT = checkSigID$Symbol[checkSigID$geneset == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"],
    Glycolysis = checkSigID$Symbol[checkSigID$geneset == "HALLMARK_GLYCOLYSIS"],
    
    IFNa_Response = checkSigID$Symbol[checkSigID$geneset == "HALLMARK_INTERFERON_ALPHA_RESPONSE"],
    IFNg_Response = checkSigID$Symbol[checkSigID$geneset == "HALLMARK_INTERFERON_GAMMA_RESPONSE"],
    PID_NFAT_TF = checkSigID$Symbol[checkSigID$geneset == "PID_NFAT_TFPATHWAY"],
    IL6_JAK_STAT3 = checkSigID$Symbol[checkSigID$geneset == "HALLMARK_IL6_JAK_STAT3_SIGNALING"],
    DNA_Repair = checkSigID$Symbol[checkSigID$geneset == "HALLMARK_DNA_REPAIR"],
    

    NK_Exh_Bulk = NK_ExhMarkerCanPass,
    NK_Res_Bulk = NK_TrmMarkerCanPass,

    CD4_Exh_Bulk = c(CD4_ExhMarkerCanPass),
    CD4_Res_Bulk = c(CD4_TrmMarkerCanPass),

    CD8_Exh_Bulk = c(CD8_ExhMarkerCanPass),
    CD8_Res_Bulk = c(CD8_TrmMarkerCanPass),
    
    All_Exh_Bulk = All_Exh_CancerPass,
    All_Res_Bulk = All_Trm_CancerPass,
    
    CD8NK_Exh_Bulk = CD8NK_Exh,
    CD8NK_Res_Bulk = CD8NK_Trm

  )
)



AllSigs <- selImmuneSigs


AllSigsGeneSet <- lapply(AllSigs, GSEABase::GeneSet)

AllSigsGeneSetName <- lapply(names(AllSigsGeneSet), function(x) {
  currentSig <- AllSigsGeneSet[[x]]
  GSEABase::setName(currentSig) <- x
  return(currentSig)
})

names(AllSigsGeneSetName) <- names(AllSigsGeneSet)

AllSigsCollection <- GSEABase::GeneSetCollection(AllSigsGeneSetName)

scoreNames <- names(selImmuneSigs)[order(names(selImmuneSigs))]















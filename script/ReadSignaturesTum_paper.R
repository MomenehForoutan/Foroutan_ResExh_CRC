##--------- load libraries and set the paths

library(singscore)
library(org.Hs.eg.db)


# outPathCD8 <- "../output/CD8/"
# outPathCD4 <- "../output/CD4/"
# outPathNK <- "../output/NK/"
outCompSigs <- "../output/CompareSigs/"
sigPath <- "../data/Signatures_literature/"

##--------- Read Marker signatures
# 
# cd8Sigs <-
#   read.table(
#     paste0(outPathCD8, "ResExhMarkers_CD8_CLs_LCM.txt"),
#     header = T,
#     sep = "\t",
#     stringsAsFactors = F,
#     check.names = F
#   )
# 
# 
# cd4Sigs <-
#   read.table(
#     paste0(outPathCD4, "ResExhMarkers_CD4_CLs_LCM.txt"),
#     header = T,
#     sep = "\t",
#     stringsAsFactors = F,
#     check.names = F
#   )
# 
# 
# nkSigs <-
#   read.table(
#     paste0(outPathNK, "ResExhMarkers_NK_CLs_LCM.txt"),
#     header = T,
#     sep = "\t",
#     stringsAsFactors = F,
#     check.names = F
#   )
# 
# 
# nkSigs$Cell.Type <- "NK"
# cd8Sigs$Cell.Type <- "CD8"
# cd4Sigs$Cell.Type <- "CD4"
# 
# # cd8Sigs$Marker <- gsub("Trm", "Res", cd8Sigs$Marker)
# # cd4Sigs$Marker <- gsub("Trm", "Res", cd4Sigs$Marker)
# # nkSigs$Marker <- gsub("Trm", "Res", nkSigs$Marker)
# 
# allMarkers <- rbind(nkSigs, cd8Sigs, cd4Sigs)
# 
# exhData <- allMarkers[allMarkers$Marker == "Exh", ]
# trmData <- allMarkers[allMarkers$Marker == "Res", ]
# 
# exhData_CanPass <- exhData[exhData$CancerPass, ]
# trmData_CanPass <- trmData[trmData$CancerPass, ]
# 
# 
# # All_Exh <- unique(exhData$gene[! exhData$gene %in% c("IL12RB2", "MARCH3")])
# # All_Res <- unique(trmData$gene[! trmData$gene %in% c("IL12RB2", "MARCH3")])
# 
# All_Exh_CancerPass <-  unique(exhData_CanPass$gene[! exhData_CanPass$gene %in% c("IL12RB2")])
# All_Res_CancerPass <-  unique(trmData_CanPass$gene[! trmData_CanPass$gene %in% c("IL12RB2")])
# 
# CD8NK_Exh <- unique(exhData_CanPass$gene[exhData_CanPass$Cell.Type  %in% c("NK", "CD8")])
# CD8NK_Res <- unique(trmData_CanPass$gene[trmData_CanPass$Cell.Type  %in% c("NK", "CD8")])


##------------ Marker genes
# NK_ExhMarker <- nkSigs$gene[nkSigs$Marker == "Exh"]
# NK_ResMarker <- nkSigs$gene[nkSigs$Marker == "Res"]
# 
# CD4_ExhMarker <- cd4Sigs$gene[cd4Sigs$Marker == "Exh"]
# CD4_TrmMarker <- cd4Sigs$gene[cd4Sigs$Marker == "Res"]
# 
# CD8_ExhMarker <- cd8Sigs$gene[cd8Sigs$Marker == "Exh"]
# CD8_ResMarker <- cd8Sigs$gene[cd8Sigs$Marker == "Res"]
# 
# 
# NK_ExhMarkerCanPass <- nkSigs$gene[nkSigs$Marker == "Exh" & nkSigs$CancerPass]
# NK_ResMarkerCanPass <- nkSigs$gene[nkSigs$Marker == "Res" & nkSigs$CancerPass]
# 
# CD4_ExhMarkerCanPass <- cd4Sigs$gene[cd4Sigs$Marker == "Exh" & cd4Sigs$CancerPass]
# CD4_ResMarkerCanPass <- cd4Sigs$gene[cd4Sigs$Marker == "Res" & cd4Sigs$CancerPass]
# 
# CD8_ExhMarkerCanPass <- cd8Sigs$gene[cd8Sigs$Marker == "Exh" & cd8Sigs$CancerPass]
# CD8_ResMarkerCanPass <- cd8Sigs$gene[cd8Sigs$Marker == "Res" & cd8Sigs$CancerPass]

##---- read signatures with stroma information - here we remove the genes that 
## are high in stroma and define anew version of signatures without them
## we then use these to score TCGA data and see if there is any difference?

markersCanPass0 <-
  read.table(
    paste0(outCompSigs, "Markers_ResExh_Bulk_Stroma.txt"),
    header = T,
    sep = "\t"
  )

NK_ExhMarkerCanPass_rmStroma <-
  markersCanPass0$gene[markersCanPass0$Cell_Marker == "NK_Exh" &
                         !markersCanPass0$HighInStroma]

NK_ResMarkerCanPass_rmStroma <-
  markersCanPass0$gene[markersCanPass0$Cell_Marker == "NK_Res" &
                         !markersCanPass0$HighInStroma]

CD4_ExhMarkerCanPass_rmStroma <-
  markersCanPass0$gene[markersCanPass0$Cell_Marker == "CD4_Exh" &
                         !markersCanPass0$HighInStroma]

CD4_ResMarkerCanPass_rmStroma <-
  markersCanPass0$gene[markersCanPass0$Cell_Marker == "CD4_Res" &
                         !markersCanPass0$HighInStroma]

CD8_ExhMarkerCanPass_rmStroma <-
  markersCanPass0$gene[markersCanPass0$Cell_Marker == "CD8_Exh" &
                         !markersCanPass0$HighInStroma]

CD8_ResMarkerCanPass_rmStroma <-
  markersCanPass0$gene[markersCanPass0$Cell_Marker == "CD8_Res" &
                         !markersCanPass0$HighInStroma]

CD8NK_Exh_rmStroma <-
  unique(c(NK_ExhMarkerCanPass_rmStroma, CD8_ExhMarkerCanPass_rmStroma))
CD8NK_Res_rmStroma <-
  unique(c(NK_ResMarkerCanPass_rmStroma, CD8_ResMarkerCanPass_rmStroma))

All_Exh_rmStroma <- unique(
  c(
    NK_ExhMarkerCanPass_rmStroma,
    CD8_ExhMarkerCanPass_rmStroma,
    CD4_ExhMarkerCanPass_rmStroma
  )
)
All_Res_rmStroma <- unique(
  c(
    NK_ResMarkerCanPass_rmStroma,
    CD8_ResMarkerCanPass_rmStroma,
    CD4_ResMarkerCanPass_rmStroma
  )
)




##-------- Read other immune signatures
# sigsImm <- readRDS(paste0(sigPath, "SigLists_Immune.RDS")) 

#--------- Read all other signatures ------------

sigsTum <- readRDS(paste0(sigPath, "TumourSigList_CRC_ResExh.RDS"))

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

# library(org.Hs.eg.db)
checkSigID$Symbol <- annotate::getSYMBOL(as.character(checkSigID$entrez), data = 'org.Hs.eg')



selSigs <- c(
  list(
    
    ## add some tum signatures
    TGFbEMT = sigsTum$TGFb_EMT$tgfb_Up,
    Epi_Thiery = sigsTum$EpiT,
    # Epi_Byers = sigsTum$EpiB,
    Mes_Thiery = sigsTum$MesT,
    # Mes_Byers = sigsTum$MesB,
    HPX_EGF_EMT = sigsTum$HpxEGF_EMT$hpxEGF_Up,
    
    # MYC_neg = sigsTum$MYC_negative,
    # MYC_pos = sigsTum$MYC_positive,
    # TP53_neg = sigsTum$TP53_negative,
    # TP53_pos = sigsTum$TP53_positive,
    DDR = sigsTum$DDR_All,
    
    Antigen_Presentation = sigsTum$AgPresentation,
    # Hypoxia20 = sigsTum$Hypoxia,
    # Hypoxia_Buffa = sigsTum$Hypoxia_Buffa,
    
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
    
    # NK_Exh_Bulk_old = NK_ExhMarkerCanPass,
    # NK_Res_Bulk_old = NK_ResMarkerCanPass,
    # 
    # CD4_Exh_Bulk_old = c(CD4_ExhMarkerCanPass),
    # CD4_Res_Bulk_old = c(CD4_ResMarkerCanPass),
    # 
    # CD8_Exh_Bulk_old = c(CD8_ExhMarkerCanPass),
    # CD8_Res_Bulk_old = c(CD8_ResMarkerCanPass),
    # 
    # All_Exh_Bulk_old = All_Exh_CancerPass,
    # All_Res_Bulk_old = All_Res_CancerPass,
    # 
    # CD8NK_Exh_Bulk_old = CD8NK_Exh,
    # CD8NK_Res_Bulk_old = CD8NK_Res,
    
    ## add signatures after removal of bad genes (high in stroma)
    
    NK_Exh_Bulk = NK_ExhMarkerCanPass_rmStroma,
    NK_Res_Bulk = NK_ResMarkerCanPass_rmStroma,
    
    CD4_Exh_Bulk = CD4_ExhMarkerCanPass_rmStroma,
    CD4_Res_Bulk = CD4_ResMarkerCanPass_rmStroma,
    
    CD8_Exh_Bulk = CD8_ExhMarkerCanPass_rmStroma,
    CD8_Res_Bulk = CD8_ResMarkerCanPass_rmStroma,
    
    CD8NK_Exh_Bulk = CD8NK_Exh_rmStroma,
    CD8NK_Res_Bulk = CD8NK_Res_rmStroma, 

    All_Exh_Bulk = All_Exh_rmStroma,
    All_Res_Bulk = All_Res_rmStroma
    
  )
)



AllSigsGeneSet <- lapply(selSigs, GSEABase::GeneSet)

AllSigsGeneSetName <- lapply(names(AllSigsGeneSet), function(x) {
  currentSig <- AllSigsGeneSet[[x]]
  GSEABase::setName(currentSig) <- x
  return(currentSig)
})

names(AllSigsGeneSetName) <- names(AllSigsGeneSet)

AllSigsCollection <- GSEABase::GeneSetCollection(AllSigsGeneSetName)

scoreNames <- names(selSigs)[order(names(selSigs))]















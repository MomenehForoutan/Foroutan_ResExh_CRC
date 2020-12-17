##--------- Read Marker signatures

outPathCD8 <- "/Users/mfor0011/Documents/projects/Monash/ResExh_CRC/output/CD8/"
outPathCD4 <- "/Users/mfor0011/Documents/projects/Monash/ResExh_CRC/output/CD4/"
outPathNK <- "/Users/mfor0011/Documents/projects/Monash/ResExh_CRC/output/NK/"
sigPath <- "/Users/mfor0011/Documents/data/signatures/ProcessedSigs/"
KalliesPath <- "/Users/mfor0011/Documents/projects/Monash/ResExh_CRC/output/Zhang/"

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

allMarkers <- rbind(nkSigs, cd8Sigs, cd4Sigs)

exhData <- allMarkers[allMarkers$Marker == "Exh", ]
trmData <- allMarkers[allMarkers$Marker == "Trm", ]

exhData_CanPass <- exhData[exhData$CancerPass, ]
trmData_CanPass <- trmData[trmData$CancerPass, ]

# write.table(allMarkers, paste0(outPath, "All_TrmExh_Markers_TNK.txt"),
# sep = "\t", row.names = F)

# table(allMarkers$Cell.Type)
# 
# CD4 CD8  NK 
# 114  70  80

##------------ Marker genes
NK_ExhMarker <- nkSigs$gene[nkSigs$Marker == "Exh"]
NK_TrmMarker <- nkSigs$gene[nkSigs$Marker == "Trm"]

CD4_ExhMarker <- cd4Sigs$gene[cd4Sigs$Marker == "Exh"]
CD4_TrmMarker <- cd4Sigs$gene[cd4Sigs$Marker == "Trm"]

CD8_ExhMarker <- cd8Sigs$gene[cd8Sigs$Marker == "Exh"]
CD8_TrmMarker <- cd8Sigs$gene[cd8Sigs$Marker == "Trm"]


NK_ExhMarkerCanPass <- nkSigs$gene[nkSigs$Marker == "Exh" & nkSigs$CancerPass]
NK_TrmMarkerCanPass <- nkSigs$gene[nkSigs$Marker == "Trm" & nkSigs$CancerPass]

CD4_ExhMarkerCanPass <- cd4Sigs$gene[cd4Sigs$Marker == "Exh" & cd4Sigs$CancerPass]
CD4_TrmMarkerCanPass <- cd4Sigs$gene[cd4Sigs$Marker == "Trm" & cd4Sigs$CancerPass]

CD8_ExhMarkerCanPass <- cd8Sigs$gene[cd8Sigs$Marker == "Exh" & cd8Sigs$CancerPass]
CD8_TrmMarkerCanPass <- cd8Sigs$gene[cd8Sigs$Marker == "Trm" & cd8Sigs$CancerPass]


#--------- Read all other signatures ---

# Read signatures ontained from Collins et al. as well as the processed immune signatures I have collected.
###------------------------- Sigs from Collins

library(singscore)
CollinsPath <- "/Users/mfor0011/Documents/projects/Monash/NK_Bulk/output/Collins/"

collinsSig <-
  read.table(
    paste0(
      CollinsPath,
      "Collins_Signatures_ILCs_CD56bright_CD56dim_GPCRs.txt"
    ),
    header = T,
    sep = "\t"
  )

## subset to up sets
collinsSig <- collinsSig[collinsSig$Direction == "Up", ]

# "TCF7"    "BCL6"    "LEF1"    "TBX21"   "CD226"   "SLAMF1"  "TNFSF14"  "IL7R"    "IL2RA"   "IL18R1"  "CCR7"    "CXCR6"   "XCL1"
# 
# TBX21 -- more in CX3CR1
# LEF1, TCF7, CCR7, SLAMF1, IL7R

###----------------- Sigs from immune and tum and TRMs

sigsImmune <- readRDS(paste0(sigPath, "SigLists_Immune.RDS"))

## This includes TRM sigs in several tissues (lung, liver, skin, brain, gut) 
## as well as TEM and TCM
sigsTRM <- readRDS(paste0(sigPath, "TRM_signatures_List.RDS"))
## subset to Up genes only; Gut_Up has 124 genes
sigsTRM_tissue <- sigsTRM$TRM_Tissues[grepl("Up", sigsTRM$TRM_Tissues$Signature), ]
sigsTRM_tissue <- sigsTRM_tissue[!sigsTRM_tissue$Signature %in% c("TCM_Up", "TEM_Up"), ]
sigsTRM_tissue <- sigsTRM_tissue[complete.cases(sigsTRM_tissue$HGNC.symbol), ]

sigsTRM_gut <- sigsTRM_tissue[sigsTRM_tissue$Signature == "Gut_Up", ]

sigsTRM_gut$HGNC.symbol[sigsTRM_gut$MGI.symbol == "Cd244"] <- "CD244"
sigsTRM_gut$HGNC.symbol[sigsTRM_gut$MGI.symbol == "Gpr132"] <- "GPR132"
sigsTRM_gut$HGNC.symbol[sigsTRM_gut$MGI.symbol == "Gp49a"] <- "LILRA5"
sigsTRM_gut$HGNC.symbol[sigsTRM_gut$MGI.symbol == "Lilrb4"] <- "LILRB4"

Gut_TRM_genes <- sigsTRM_gut$HGNC.symbol
Gut_TRM_genes <- Gut_TRM_genes[! duplicated(Gut_TRM_genes)]
Gut_TRM_genes <- Gut_TRM_genes[complete.cases(Gut_TRM_genes)]
Gut_TRM_genes <- c(Gut_TRM_genes, "XCL1")

sigsImmune <- c(sigsImmune, list(TRM_Gut = Gut_TRM_genes))

sigsTum <- readRDS(paste0(sigPath, "SigLists_Tumour.RDS"))


CD8_TEX_JA <- sigsImmune$T.CD8.EXHAUSTED
CD4_TEX_JA <- sigsImmune$T.CD4.EXHAUSTED
TEX_ImmuneSig <- sigsImmune$ExhaustedT   ## check where it comes from

###----------------------- Sigs from Guo (NSCLC)

CD8_TEX_Guo <- read.csv(paste0(sigPath, "../TRM/Guo_NSCLC/CD8_TEX_Guo.csv"),
                        stringsAsFactors = F)
CD4_TEX_Guo <- read.csv(paste0(sigPath, "../TRM/Guo_NSCLC/CD4_TEX_Guo.csv"),
                        stringsAsFactors = F)

##----------------------- Sigs from Li et al Melanoma

CD8_Dysfunc_Li <- read.csv( "/Users/mfor0011/Documents/data/signatures/TRM/Li_Melanoma/Li_Melanoma_CD8Dysfunctional.csv",
                            stringsAsFactors = F)
CD4_TregsTum_Li <- read.csv( "/Users/mfor0011/Documents/data/signatures/TRM/Li_Melanoma/Li_Melanoma_Tregs_vs_naive.csv",
                             stringsAsFactors = F)

###------------------------- Sigs from Corridoni (Colitis)
# We subset to Trm, IL26 (similar to exhausted cells), Double Positive(DP) cells (CD8 Tregs!), and IEL cells (TYROBP+/-), as well as GZMK-eff(2).


corridoni <- read.csv(paste0(sigPath, "../TRM/Corridoni_Colitis_ClusterMarkers.csv"), stringsAsFactors = F)


corri <- corridoni[corridoni$Cluster %in% c("Trm", 
                                            "CD8+ CD4+ FOXP3+", 
                                            "GZMK+ Effectors (2)",
                                            "IL26+", 
                                            "TYROBP- IELs", 
                                            "TYROBP+ IELs"), ]

## subset to those with higher logFC and lower FDR
# plot(corri$avg_logFC, -log10(corri$FDR))

corri <- corri %>% filter(FDR < 0.05 & avg_logFC > 0.3)

## 583
corriGenes <- unique(corri$Gene)

###------------------------- Sigs from From Nicole's paper
## DEG shared across all cell types and tissues

resident3Genes <- c("RGS1", "CXCR6", "ITGAE")

## DEGs shared between trNK cells in lung and BM, but not with CD8+ TRM
residentNK_lungBM_notTRM <- c(
  "IFNG",
  "STK16",
  "KIAA1671",
  "ATXN1",
  "CCL5",
  "MTFP1",
  "AGTRAP",
  "IKZF3",
  "RGS2",
  "PLA2G16",
  "CCDC167",
  "GZMA"
)


residentNK_lungOnly <- c(
  "LGALS3",
  "HLA-DQA2",
  "HLA-DRB5",
  "HLA-DQA1",
  "HLA-DQB1",
  "LGALS1",
  "HLA-DRB1",
  "HLA-DMA",
  "ANXA2",
  "S100A6"
)
## positive at least in two cpmparisons in the figure:
# - NK resident (lung/BM) and TRM (lung and spleen)

residentGenes <-  as.character(quote(
  c(
    RGS1,
    CXCR6,
    ITGAE,
    ALOX5AP,
    SLC6A6,
    DAPK2,
    GLIPR1,
    TENT5C,
    CD96,
    METTL9,
    SPRY2,
    ZNF331,
    IFNG,
    STK16,
    KIAA1671,
    ATXN1,
    CCL5,
    MTFP1,
    AGTRAP,
    IKZF3,
    RGS2,
    PLA2G16,
    CCDC167,
    GZMA,
    ITGA1,
    JAML,
    KRT81,
    LAG3,
    HSPA1A,
    KRT86,
    PERP,
    GLUL,
    CITED2,
    TBCD,
    SAMSN1,
    LDLRAD4,
    CYTIP,
    PTPN22,
    PLP2,
    IRF4,
    CAMK4,
    RBPJ,
    ITM2C,
    FRMD4B,
    CD82,
    DUSP4,
    ZNF683
  )
))[-1]



##---- only in Trm:
lungSpleenTRM <-  as.character(quote(
  c(
    ITGA1,
    JAML,
    KRT81,
    LAG3,
    HSPA1A,
    KRT86,
    PERP,
    GLUL,
    CITED2,
    TBCD,
    SAMSN1,
    LDLRAD4,
    CYTIP,
    PTPN22,
    PLP2,
    IRF4,
    CAMK4,
    RBPJ,
    ITM2C,
    FRMD4B,
    CD82,
    DUSP4,
    ZNF683,
    RGS1,
    CXCR6,
    ITGAE,
    ALOX5AP,
    SLC6A6,
    DAPK2,
    GLIPR1,
    TENT5C,
    CD96,
    METTL9,
    SPRY2,
    ZNF331,
    ## In spleen
    EGR1,
    Aff3,
    IFNG
  )
))[-1]

## read all signatures from Nicole Maq
resNichole <- read.csv(paste0(sigPath, "../TRM/TRNK_Nichole_LungSpleenBM.csv"),
                       stringsAsFactors = F)

## get logFC> 0
resNichole <- resNichole[resNichole$log2FoldChange > 0, ]

# They all are  NKG2A+CD16- NK cells , so we remove these:

resNichole$Comparison <- gsub(" NKG2A\\+CD16\\- NK cells", "", resNichole$Comparison )
resNichole$Comparison <- gsub(" NKG2A\\+CD16\\-", "", resNichole$Comparison )
resNichole$Comparison <- gsub(" NK cells", "", resNichole$Comparison )

table(resNichole$Comparison)
## there is only one genes for CD69+CD49a+CD103+ versus CD69+CD49a+CD103-, which we remove: HSPA1A
resNichole[resNichole$Comparison == "CD69+CD49a+CD103+ versus CD69+CD49a+CD103-", ]
resNichole <- resNichole[! resNichole$Comparison == "CD69+CD49a+CD103+ versus CD69+CD49a+CD103-", ]

resNicholeSigs <- sapply(unique(resNichole$Comparison), function(x){
  d0 <- resNichole[resNichole$Comparison == x, ]
  return(as.character(d0$Gene.name))
})

names(resNicholeSigs) <- c(
  "CD69pCD49apCD103p_vs_CD69n",
  "CD69pCD49apCD103n_vs_CD69n",
  "CD69sp_vs_CD69n",
  "CD69pCD49apCD103p_vs_CD69sp",
  "CD69pCD49apCD103n_vs_CD69sp"
)


###-------------------- Signatures from Kallies

KalliesSigs <- readRDS(paste0(KalliesPath, "Kallies_Sigs_MouseHumanSymbols.RDS"))

Tex <- KalliesSigs$Tex
Tpex <- KalliesSigs$Tpex
TRM_IDs_up <- KalliesSigs$TRM_IDs_up
ExhCore <- KalliesSigs$ExhCore


##---- NK genes from review paper:
nkGenes <- as.character(quote(
  c(
    ##------- NK receptor genes:
    ## PB CD56bright
    KLRD1,
    ## KLRD1  is the same as CD94
    KLRF1,
    ## KLRF1 is the same as NKp80
    NCAM1,
    ## same as CD56 also in resident NK, not in Kidney
    IL2RA,
    IL7RA,
    KIT,
    ## CD117
    SELL,
    ## CD62L
    NKG2A,
    ## in most resident NKs but not in uterine, and not in CD56 dim
    ITGA5,
    ## CD49e
    CXCR3,
    CCR7,
    NCR1,
    ## NKp46
    
    ## PB CD56dim
    FCGR3A,
    ## CD16
    ##KIRs  also in lung and uterine resident NK
    KIR3DL3,
    KIR2L3,
    KIR2DP1,
    KIR2DL1,
    KIR3DP1,
    KIR2DL4,
    KIR3DL1,
    KIR2DS4,
    KIR3DL2,
    KIR2DS1,
    KIR2DS2,
    KIR2DS3,
    KIR2DS5,
    KIR2DL2,
    KIR2DL5,
    KIR3DS1,
    S1PR5,
    ## S1P5
    CX3CR1,
    CXCR2,
    CXCR1,
    ## PB adaptive NK:
    CD2,
    KLRC2,
    ## NKG2C
    LILRB1,
    ## CD85j
    B3GAT1,
    ## CD57
    ## Bone marrow , SLT,Gut:
    ITGA1,
    ## CD49a,
    ITGAE,
    ## CD103,
    CD69,
    CCR5,
    ## not in lung, uterine and kidney
    CXCR6,
    ## not in lung, uterine and kidney
    ## Uterine:
    CD9,
    ##----- NK secretory proteins
    IFNG,
    ## In PB CD56bright and PB adaptive NK, and not resident NK
    PRF1,
    ## In PB CD56dim
    GZMA,
    GZMB,
    GZMH,
    GZMK,
    GZMM ## In PB CD56dim)
  )
)) [-1]


## from de Andrade:
# NK cells are: CD45+CD56+CD3–CD14–CD15–CD163–
"/Users/mfor0011/Documents/data/signatures/ProcessedSigs/../TRM/ILCs_NK/"
# From Björklund et al https://www-nature-com.ezproxy.lib.monash.edu.au/articles/ni.3368

# mle: log2 fold expression difference value, reported as maximum likelihood estimate (mle)
# lb, ub: 95% lower bound (lb) and upper bound (ub) of the log2 fold expression difference value
# ce: conservative estiamte of the log2 fold expression difference value (upper or lower bound, whichever is closer to 0, or 0 if the 95% confidence interval includes 0)
# Z: statistical significance of the expression differences, expressed as a signed Z score (can be directly converted to P values; roughly a p-value of 1e-2 would correspond to a Z score around 2.3, and a p-value of 1e-3 to 3.1; One can convert Z scores back to p-values using pnorm(Z,lower.tail=F) for upregulated genes and pnorm(Z,lower.tail=T) for downregulated genes; Z scores are signed, so that the positive values correspond to up-regulated genes, negative to downregulated)
# cZ: Z score adjusted for multiple hypothesis testing

nk_Bjorklund <- read.csv(paste0(sigPath, "../TRM/ILCs_NK/NK_vs_ILCs_Bjorklund.csv"),
                         stringsAsFactors = F)
ILC1_Bjorklund <- read.csv(paste0(sigPath, "../TRM/ILCs_NK/ILC1_vs_ILC2.3_Bjorklund.csv"),
                           stringsAsFactors = F)
ILC2_Bjorklund <- read.csv(paste0(sigPath, "../TRM/ILCs_NK/ILC2_vs_ILC1.3_Bjorklund.csv"),
                           stringsAsFactors = F)
ILC3_Bjorklund <- read.csv(paste0(sigPath, "../TRM/ILCs_NK/ILC3_vs_ILC1.2_Bjorklund.csv"),
                           stringsAsFactors = F)

NK_vs_ILCs <- nk_Bjorklund[nk_Bjorklund$pvalue2sided.adj < 0.05 & 
                             nk_Bjorklund$mle > 2, ]

ILC1_vs_ILC2.3 <- ILC1_Bjorklund[ILC1_Bjorklund$pvalue2sided.adj < 0.05 & 
                                   ILC1_Bjorklund$mle > 2, ]

ILC2_vs_ILC1.3 <- ILC2_Bjorklund[ILC2_Bjorklund$pvalue2sided.adj < 0.05 & 
                                   ILC2_Bjorklund$mle > 2, ]

ILC3_vs_ILC1.2 <- ILC3_Bjorklund[ILC3_Bjorklund$pvalue2sided.adj < 0.01 & 
                                   ILC3_Bjorklund$mle > 3, ]

dim(NK_vs_ILCs)
dim(ILC1_vs_ILC2.3)
dim(ILC2_vs_ILC1.3)
dim(ILC3_vs_ILC1.2)


##---- some cancer related signatures:
# Suppl table 5: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2816644/#sup1
# 6 genes overlapping between the two Hypoxia signatures
# [1] "BNIP3"  "SLC2A1" "PGK1"   "NDRG1"  "P4HA1" 
# [6] "ADM"  
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




CanEvasion <-
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
# 
# The glycolysis score for each cell was calculated as the average z-score over all genes in the HALLMARK_GLYCOLYSIS gene set55. ROS score for each cell was calculated as the average z-score of all genes in the GO_RESPONSE_TO_OXIDATIVE_STRESS gene set56. NFAT score for each cell was calculated as the average z-score over all genes in the PID_NFAT_TFPATHWAY gene set

hum <- msigdf::msigdf.human

checkSigName <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "CHUANG_OXIDATIVE_STRESS_RESPONSE_UP",
  "MIKHAYLOVA_OXIDATIVE_STRESS_RESPONSE_VIA_VHL_UP",
  "GO_RESPONSE_TO_OXIDATIVE_STRESS",
  "HALLMARK_GLYCOLYSIS",
  "PID_NFAT_TFPATHWAY"
)

checkSigID <- hum[hum$geneset %in% checkSigName, ]
# glyco <- as.character(hum$entrez[hum$geneset == "HALLMARK_GLYCOLYSIS"])

library(org.Hs.eg.db)
checkSigID$Symbol <- annotate::getSYMBOL(as.character(checkSigID$entrez), data='org.Hs.eg')



selImmuneSigs <- c(sigsImmune[c(
  "TRM_Review",
  "TRM_TGFb_Unstim_Up",
  "TRM_TGFbIL2_Unstim_Up",
  "TRM_TGFbIL2_IL2_Up",
  "TRM_Nizard_Up", 
  "TRM_Gut",
  "ExhaustedT", 
  "Stem_T",
  "TerminallyDiff_T",
  "ILC1", 
  "Interferon_imsig", 
  "Proliferation_imsig", 
  "Translation_imsig", 
  "TEM_Up",
  "TCM_Up",
  "NK_cur",
  "NK",
  "NK_Xiong",
  "NK_ImmGene",
  "NK_imsig",
  "T_ImmGene",
  "T_imsig",
  "T_Xiong",
  "T_gd1",
  "T_gd2",
  "T.CD4",
  "T.CD8" ,
  "T.CELL", 
  "T.CD8.EXHAUSTED",
  "T.CD4.EXHAUSTED",
  "T.CD4.NAIVE",
  "T.CD8.NAIVE",
  "T.CD4.TREG",
  "T.CD8.CYTOTOXIC"
)], 
resNicholeSigs, 
list(resident3Genes = resident3Genes, 
     residentGenes = residentGenes,
     residentNK_lungOnly = residentNK_lungOnly,
     residentNK_lungBM_notTRM = residentNK_lungBM_notTRM, 
     TRM_lungSpleen = lungSpleenTRM, 
     TGFbEMT = sigsTum$TGFb_EMT$tgfb_Up, 
     CD8_TEX_Kallies = Tex$HGNC.symbol,
     CD8_TPEX_Kallies = Tpex$HGNC.symbol,
     CD8_ExhCore_Kallies = ExhCore$HGNC.symbol,
     CD8_TRM_Kallies = TRM_IDs_up$HGNC.symbol,
     CD8_Dysfunc_Li = CD8_Dysfunc_Li$X,
     CD4_TEX_Guo = CD4_TEX_Guo$Gene.Symbol,
     CD4_TEX_JA = CD4_TEX_JA,
     CD4_TregsTum_Li = CD4_TregsTum_Li$X,
     NK_vs_ILCs = NK_vs_ILCs$Gene,
     ILC1_vs_ILC2.3 = ILC1_vs_ILC2.3$Gene,
     ILC2_vs_ILC1.3 = ILC2_vs_ILC1.3$Gene,
     ILC3_vs_ILC1.2 = ILC2_vs_ILC1.3$Gene,
     
     ## Also add some tum signatures
     DDR = sigsTum$DDR_All,
     TP53_neg = sigsTum$TP53_negative,
     Hypoxia20 = sigsTum$Hypoxia, 
     Hypoxia_Buffa = Hypoxia_Buffa, 
     Antigen_Presentation = CanEvasion, 
     
     OX_StressResponse = checkSigID$Symbol[checkSigID$geneset == "CHUANG_OXIDATIVE_STRESS_RESPONSE_UP"], 
     OX_StressResponseGO = checkSigID$Symbol[checkSigID$geneset == "GO_RESPONSE_TO_OXIDATIVE_STRESS"], 
     Glycolysis = checkSigID$Symbol[checkSigID$geneset == "HALLMARK_GLYCOLYSIS"],
     OXPHOS = checkSigID$Symbol[checkSigID$geneset == "HALLMARK_OXIDATIVE_PHOSPHORYLATION"], 
     OX_StressResponseVHL = checkSigID$Symbol[checkSigID$geneset == "MIKHAYLOVA_OXIDATIVE_STRESS_RESPONSE_VIA_VHL_UP"], 
     PID_NFAT_TFPATHWAY = checkSigID$Symbol[checkSigID$geneset == "PID_NFAT_TFPATHWAY"], 
     
     NK_Exh = NK_ExhMarker,
     NK_Exh_CancerPass = NK_ExhMarkerCanPass,
     NK_Trm = NK_TrmMarker,
     NK_Trm_CancerPass = NK_TrmMarkerCanPass,
     
     CD4_Exh = c(CD4_ExhMarker),
     CD4_Exh_CancerPass = c(CD4_ExhMarkerCanPass),
     CD4_Trm = c(CD4_TrmMarker),
     CD4_Trm_CancerPass = c(CD4_TrmMarkerCanPass),
     
     CD8_Exh = c(CD8_ExhMarker),
     CD8_Exh_CancerPass = c(CD8_ExhMarkerCanPass),
     CD8_Trm = c(CD8_TrmMarker),
     CD8_Trm_CancerPass = c(CD8_TrmMarkerCanPass),
     
     All_Exh = unique(exhData$gene),
     All_Trm = unique(trmData$gene)
     # All_Exh = rownames(exhData),
     # All_Trm = rownames(trmData)
))




## change the name of ILC1 in the selImmuneSig as we have another ILC1 signature from Collins data below.
names(selImmuneSigs)[names(selImmuneSigs) == "ILC1"] <- "ILC1_Fer"

AllSigs <- c(list(
  CD56bright = collinsSig$Genes[collinsSig$Signature == "CD56bright"], 
  CD56dim = collinsSig$Genes[collinsSig$Signature == "CD56dim"], 
  CD56dimCD57neg = collinsSig$Genes[collinsSig$Signature == "CD56dimCD57neg"], 
  CD56dimCD57pos = collinsSig$Genes[collinsSig$Signature == "CD56dimCD57pos"], 
  ILC1 = collinsSig$Genes[collinsSig$Signature == "ILC1"], 
  ILC3 = collinsSig$Genes[collinsSig$Signature == "ILC3"]
  ),
  selImmuneSigs)
  
AllSigsGeneSet <- lapply(AllSigs, GSEABase::GeneSet)

AllSigsGeneSetName <- lapply(names(AllSigsGeneSet), function(x) {
  currentSig <- AllSigsGeneSet[[x]]
  GSEABase::setName(currentSig) <- x
  return(currentSig)
})

names(AllSigsGeneSetName) <- names(AllSigsGeneSet)

AllSigsCollection <- GSEABase::GeneSetCollection(AllSigsGeneSetName)


kpGenes <-
  unique(c(unique(unlist(selImmuneSigs)), 
           unique(collinsSig$Genes)))

# saveRDS(kpGenes, "kpGenes.RDS")















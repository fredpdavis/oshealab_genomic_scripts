# Author: Alejandro Villarino
# 181218_1419

# load packages
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(msigdbr)
library(org.Mm.eg.db)
library(plyr)
library(dplyr)
library(data.table)

# name variables
WorkDir_Path = "~/Google Drive/STAT5 & Metabolism - Metabolon/STAT5 & Metabolism - RNAseq/111117 S5 Metabolism RNAseq - Updated Analysis/S5 Metabolism RNAseq - Pathway Analysis/S5 Metabolism RNAseq - clusterProfiler/111818_Th0_Rest_IL2_IL21_TCR_select_clusterProfiler/WT_Rest_IL2_v_WT_Rest_clusterProfiler"

DEG_Master_Path = "~/Google Drive/STAT5 & Metabolism - Metabolon/STAT5 & Metabolism - RNAseq/111117 S5 Metabolism RNAseq - Updated Analysis/S5 Metabolism RNAseq - DEG/RNAseq DEGs - EdgeR_glm/111818_Th0_Rest_IL2_IL21_TCR_select_EDgeR_glmTest/WT_Rest_IL2-WT_Rest_eRglmQLFTest_DEG_logFC1_BHp05.txt"

RL_Master_Path = "~/Google Drive/STAT5 & Metabolism - Metabolon/STAT5 & Metabolism - RNAseq/111117 S5 Metabolism RNAseq - Updated Analysis/S5 Metabolism RNAseq - DEG/RNAseq DEGs - EdgeR_glm/111818_Th0_Rest_IL2_IL21_TCR_select_EDgeR_glmTest/WT_Th0_Rest_IL2_IL21_TCR_select_eRglmQLFTest_results_table_111818.txt"

Output_file_prefix = "WT_Rest_IL2_v_WT_Rest_CTL_eRTest_DEG_logFC1_BHp05"
GSEA_output_file_prefix = "WT_Rest_IL2_v_WT_Rest_CTL_eRTest"

# set working directory
setwd(WorkDir_Path)

# read edgeR DEG table
DEG_Master = read.delim(DEG_Master_Path, row.names=NULL)

# read edgeR results table_includes GSEArank
RL_Master = read.delim(RL_Master_Path, row.names=NULL)

# prepare DEG dataframe for hypergeometric testing
DEG1_RL <- subset(DEG_Master, select=c("SYMBOL","GSEArank"))

# prepare Ranked List dataframe for GSEA
RL1 <- subset(RL_Master, select=c("Table1.SYMBOL","Table1.GSEArank"))
RL1 <- rename(RL1, c("Table1.SYMBOL" = "SYMBOL", "Table1.GSEArank" = "GSEArank"))

# Add ENTREZ and KEGG ID to DEG dataframe
ENTREZ_anno1 <- bitr(DEG1_RL$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
KEGG_anno1 <- bitr_kegg(ENTREZ_anno1$ENTREZID, fromType='ncbi-geneid', toType='kegg', organism='mmu')
KEGG_anno1 <- rename(KEGG_anno1, c("kegg" = "KEGG", "ncbi-geneid" = "ENTREZID"))
DEG1_RL <- merge(DEG1_RL, ENTREZ_anno1, by="SYMBOL")
DEG1_RL <- merge(DEG1_RL, KEGG_anno1, by="ENTREZID")
DEG1_RL <- DEG1_RL[c("SYMBOL","GSEArank","ENTREZID","KEGG")]

# Add ENTREZ and KEGG ID to Ranked List dataframe
ENTREZ_anno2 <- bitr(RL1$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
KEGG_anno2 <- bitr_kegg(ENTREZ_anno2$ENTREZID, fromType='ncbi-geneid', toType='kegg', organism='mmu')
KEGG_anno2 <- rename(KEGG_anno2, c("kegg" = "KEGG", "ncbi-geneid" = "ENTREZID"))
RL1 <- merge(RL1, ENTREZ_anno2, by="SYMBOL")
RL1 <- merge(RL1, KEGG_anno2, by="ENTREZID")
RL1 <- RL1[c("ENTREZID","GSEArank","SYMBOL","KEGG")]

# # Prepare GSEA input list_must be sorted by rank in decreasing order
RL1_GR = RL1[,"GSEArank"]
names(RL1_GR) = as.character(RL1[,"ENTREZID"])
RL1_GSEA = sort(RL1_GR, decreasing = TRUE)

# Create dataframes from selected MsigDB collections
mouse_msigDB = msigdbr(species = "Mus musculus")
mouse_msigDB_summary <- count(mouse_msigDB, c("gs_cat", "gs_subcat"))
mouse_msigDB_H = msigdbr(species = "Mus musculus", category = "H")
mouse_msigDB_H = subset(mouse_msigDB_H, select=c("gs_name","entrez_gene"))
mouse_msigDB_C2_CGP = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
mouse_msigDB_C2_CGP = subset(mouse_msigDB_C2_CGP, select=c("gs_name","entrez_gene")) 
mouse_msigDB_C2_CP = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP")
mouse_msigDB_C2_CP = subset(mouse_msigDB_C2_CP, select=c("gs_name","entrez_gene"))
mouse_msigDB_C2_CPb = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:BIOCARTA")
mouse_msigDB_C2_CPb = subset(mouse_msigDB_C2_CPb, select=c("gs_name","entrez_gene"))
mouse_msigDB_C3_TFT = msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT")
mouse_msigDB_C3_TFT = subset(mouse_msigDB_C3_TFT, select=c("gs_name","entrez_gene"))
mouse_msigDB_C4 = msigdbr(species = "Mus musculus", category = "C4")
mouse_msigDB_C4 = subset(mouse_msigDB_C4, select=c("gs_name","entrez_gene"))
mouse_msigDB_C6 = msigdbr(species = "Mus musculus", category = "C6")
mouse_msigDB_C6 = subset(mouse_msigDB_C6, select=c("gs_name","entrez_gene")) 
mouse_msigDB_C7 = msigdbr(species = "Mus musculus", category = "C7")
mouse_msigDB_C7 = subset(mouse_msigDB_C7, select=c("gs_name","entrez_gene"))

#Tests KEGG, GO, Reactome and MsigDB
# Run Hypergeometric test for pathway enrichment among DEGs
KEGG1_BHp20 <- enrichKEGG(gene = DEG1_RL$KEGG,
                          organism = 'mmu',
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.2)
KEGG1_BHp20 <- setReadable(KEGG1_BHp20, 'org.Mm.eg.db', 'ENTREZID')
KEGG1_BHp20_DF <- as.data.frame(KEGG1_BHp20)

KEGG1_all <- enrichKEGG(gene = DEG1_RL$KEGG,
                        organism = 'mmu',
                        minGSSize    = 1,
                        maxGSSize    = 10000,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1)
KEGG1_all <- setReadable(KEGG1_all, 'org.Mm.eg.db', 'ENTREZID')
KEGG1_all_DF <- as.data.frame(KEGG1_all)

KEGG1_M_all <- enrichMKEGG(gene = DEG1_RL$KEGG, 
                           organism = 'mmu',
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 1,
                           qvalueCutoff  = 1)
KEGG1_M_all <- setReadable(KEGG1_M_all, 'org.Mm.eg.db', 'ENTREZID')
KEGG1_M_all_DF <- as.data.frame(KEGG1_M_all)

GO1_BP_BHp001 <- enrichGO(gene = DEG1_RL$SYMBOL,
                          OrgDb = org.Mm.eg.db,
                          keyType = 'SYMBOL',
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.001)
GO1_BP_BHp001_DF <- as.data.frame(GO1_BP_BHp001)

GO1_MF_BHp05 <- enrichGO(gene = DEG1_RL$SYMBOL,
                         OrgDb = org.Mm.eg.db,
                         keyType = 'SYMBOL',
                         ont = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05)
GO1_MF_BHp05_DF <- as.data.frame(GO1_MF_BHp05)

Reactome1_BHp05 <- enrichPathway(DEG1_RL$ENTREZID, 
                                 organism = "mouse",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05,
                                 readable = TRUE)
Reactome1_BHp05_DF <- as.data.frame(Reactome1_BHp05)

Reactome1_all <- enrichPathway(DEG1_RL$ENTREZID, 
                               organism = "mouse",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 1,
                               qvalueCutoff = 1,
                               minGSSize = 1,
                               maxGSSize = 10000,
                               readable = TRUE)
Reactome1_all_DF <- as.data.frame(Reactome1_all)

MsigDB_H_BHp20 <- enricher(gene = DEG1_RL$ENTREZID, 
                           pAdjustMethod = "BH",                               	
                           pvalueCutoff = 0.2,
                           TERM2GENE = mouse_msigDB_H)
MsigDB_H_BHp20 <- setReadable(MsigDB_H_BHp20, 'org.Mm.eg.db', 'ENTREZID')
MsigDB_H_BHp20_DF <- as.data.frame(MsigDB_H_BHp20)

MsigDB_H_all <- enricher(gene = DEG1_RL$ENTREZID, 
                         pAdjustMethod = "BH",                               	
                         pvalueCutoff = 1,
                         qvalueCutoff = 1,
                         TERM2GENE = mouse_msigDB_H)
MsigDB_H_all <- setReadable(MsigDB_H_all, 'org.Mm.eg.db', 'ENTREZID')
MsigDB_H_all_DF <- as.data.frame(MsigDB_H_all)

MsigDB_C2_CGP_BHp0001 <- enricher(gene = DEG1_RL$ENTREZID, 
                                  pAdjustMethod = "BH",                               	
                                  pvalueCutoff = 0.0001,
                                  TERM2GENE = mouse_msigDB_C2_CGP)
MsigDB_C2_CGP_BHp0001 <- setReadable(MsigDB_C2_CGP_BHp0001, 'org.Mm.eg.db', 'ENTREZID')
MsigDB_C2_CGP_BHp0001_DF <- as.data.frame(MsigDB_C2_CGP_BHp0001)

MsigDB_C2_CP_BHp20 <- enricher(gene = DEG1_RL$ENTREZID, 
                               pAdjustMethod = "BH",                               	
                               pvalueCutoff = 0.2,
                               TERM2GENE = mouse_msigDB_C2_CP)
MsigDB_C2_CP_BHp20 <- setReadable(MsigDB_C2_CP_BHp20, 'org.Mm.eg.db', 'ENTREZID')
MsigDB_C2_CP_BHp20_DF <- as.data.frame(MsigDB_C2_CP_BHp20)

MsigDB_C2_CPb_BHp20 <- enricher(gene = DEG1_RL$ENTREZID, 
                                pAdjustMethod = "BH",                               	
                                pvalueCutoff = 0.2,
                                TERM2GENE = mouse_msigDB_C2_CPb)
MsigDB_C2_CPb_BHp20 <- setReadable(MsigDB_C2_CPb_BHp20, 'org.Mm.eg.db', 'ENTREZID')
MsigDB_C2_CPb_BHp20_DF <- as.data.frame(MsigDB_C2_CPb_BHp20)

MsigDB_C3_TFT_BHp20 <- enricher(gene = DEG1_RL$ENTREZID, 
                                pAdjustMethod = "BH",                               	
                                pvalueCutoff = 0.2,
                                TERM2GENE = mouse_msigDB_C3_TFT)
MsigDB_C3_TFT_BHp20 <- setReadable(MsigDB_C3_TFT_BHp20, 'org.Mm.eg.db', 'ENTREZID')
MsigDB_C3_TFT_BHp20_DF <- as.data.frame(MsigDB_C3_TFT_BHp20)

MsigDB_C6_BHp20 <- enricher(gene = DEG1_RL$ENTREZID, 
                            pAdjustMethod = "BH",                               	
                            pvalueCutoff = 0.2,
                            TERM2GENE = mouse_msigDB_C6)
MsigDB_C6_BHp20 <- setReadable(MsigDB_C6_BHp20, 'org.Mm.eg.db', 'ENTREZID')
MsigDB_C6_BHp20_DF <- as.data.frame(MsigDB_C6_BHp20)

MsigDB_C7_BHp0001 <- enricher(gene = DEG1_RL$ENTREZID, 
                              pAdjustMethod = "BH",                               	
                              pvalueCutoff = 0.0001,
                              TERM2GENE = mouse_msigDB_C7)
MsigDB_C7_BHp0001 <- setReadable(MsigDB_C7_BHp0001, 'org.Mm.eg.db', 'ENTREZID')
MsigDB_C7_BHp0001_DF <- as.data.frame(MsigDB_C7_BHp0001)

# Count pathways meeting specified PValues and output statistics to dataframe
hSum1 <- as.data.frame(nrow(KEGG1_BHp20_DF))
hSum2 <- as.data.frame(nrow(KEGG1_all_DF))
hSum3 <- as.data.frame(nrow(KEGG1_M_all_DF))
hSum4 <- as.data.frame(nrow(GO1_BP_BHp001_DF))
hSum5 <- as.data.frame(nrow(GO1_MF_BHp05_DF))
hSum6 <- as.data.frame(nrow(Reactome1_BHp05_DF))
hSum7 <- as.data.frame(nrow(MsigDB_H_BHp20_DF))
hSum8 <- as.data.frame(nrow(MsigDB_H_all_DF))
hSum9 <- as.data.frame(nrow(MsigDB_C2_CGP_BHp0001_DF))
hSum10 <- as.data.frame(nrow(MsigDB_C2_CP_BHp20_DF))
hSum11 <- as.data.frame(nrow(MsigDB_C2_CPb_BHp20_DF))
hSum12 <- as.data.frame(nrow(MsigDB_C3_TFT_BHp20_DF))
hSum13 <- as.data.frame(nrow(MsigDB_C6_BHp20_DF))
hSum14 <- as.data.frame(nrow(MsigDB_C7_BHp0001_DF))
hSum15 <- as.data.frame(nrow(Reactome1_all_DF))
hSummary <- as.data.frame(cbind(mget(ls(pattern = "hSum\\d+"))))

# Rank calculated by multiplying LogFC and negative log10 FDR for all genes passing EdgeR filering
# Run GSEA against transcriptome wide ranked list
gKEGG1_BHp10 <- gseKEGG(geneList = RL1_GSEA,
                        organism = 'mmu',
                        nPerm = 10000,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.1,
                        verbose = FALSE)
gKEGG1_BHp10 <- setReadable(gKEGG1_BHp10, 'org.Mm.eg.db', 'ENTREZID')
gKEGG1_BHp10_DF <- as.data.frame(gKEGG1_BHp10)

gKEGG1_all <- gseKEGG(geneList = RL1_GSEA,
                      organism = 'mmu',
                      nPerm = 10000,
                      minGSSize    = 1,
                      maxGSSize    = 10000,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1,
                      verbose = FALSE)
gKEGG1_all <- setReadable(gKEGG1_all, 'org.Mm.eg.db', 'ENTREZID')
gKEGG1_all_DF <- as.data.frame(gKEGG1_all)

gGO1_BP_BHp05 <- gseGO(geneList = RL1_GSEA,
                       OrgDb = org.Mm.eg.db,
                       ont = "BP",
                       nPerm = 10000,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       minGSSize = 10,
                       verbose = FALSE)
gGO1_BP_BHp05 <- setReadable(gGO1_BP_BHp05, 'org.Mm.eg.db', 'ENTREZID')
gGO1_BP_BHp05_DF <- as.data.frame(gGO1_BP_BHp05)
gGO1_BP_BHp05_light_DF <- as.data.frame(simplify(gGO1_BP_BHp05, cutoff=0.7, by="p.adjust", select_fun=min))

gGO1_MF_BHp10 <- gseGO(geneList = RL1_GSEA,
                       OrgDb = org.Mm.eg.db,
                       ont = "MF",
                       nPerm = 10000,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.1,
                       minGSSize = 10,
                       verbose = FALSE)
gGO1_MF_BHp10 <- setReadable(gGO1_MF_BHp10, 'org.Mm.eg.db', 'ENTREZID')
gGO1_MF_BHp10_DF <- as.data.frame(gGO1_MF_BHp10)
gGO1_MF_BHp10_light_DF <- as.data.frame(simplify(gGO1_MF_BHp10, cutoff=0.7, by="p.adjust", select_fun=min))

gReactome1_BHp05 <- gsePathway(RL1_GSEA, 
                               organism = "mouse",
                               nPerm = 10000,
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,
                               verbose=FALSE)
gReactome1_BHp05 <- setReadable(gReactome1_BHp05, 'org.Mm.eg.db', 'ENTREZID')
gReactome1_BHp05_DF <- as.data.frame(gReactome1_BHp05)

gReactome1_all <- gsePathway(RL1_GSEA, 
                             organism = "mouse",
                             nPerm = 10000,
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 1,
                             minGSSize = 1,
                             maxGSSize = 10000,
                             verbose=FALSE)
gReactome1_all <- setReadable(gReactome1_all, 'org.Mm.eg.db', 'ENTREZID')
gReactome1_all_DF <- as.data.frame(gReactome1_all)

gMsigDB_H_BHp20 <- GSEA(gene = RL1_GSEA, 
                        nPerm = 10000,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.2,
                        TERM2GENE = mouse_msigDB_H,
                        verbose=FALSE)
gMsigDB_H_BHp20 <- setReadable(gMsigDB_H_BHp20, 'org.Mm.eg.db', 'ENTREZID')
gMsigDB_H_BHp20_DF <- as.data.frame(gMsigDB_H_BHp20)

gMsigDB_H_all <- GSEA(gene = RL1_GSEA, 
                      nPerm = 10000,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1,
                      TERM2GENE = mouse_msigDB_H,
                      verbose=FALSE)
gMsigDB_H_all <- setReadable(gMsigDB_H_all, 'org.Mm.eg.db', 'ENTREZID')
gMsigDB_H_all_DF <- as.data.frame(gMsigDB_H_all)

gMsigDB_C2_CGP_BHp01 <- GSEA(gene = RL1_GSEA, 
                             nPerm = 10000,
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.01,
                             TERM2GENE = mouse_msigDB_C2_CGP,
                             verbose=FALSE)
gMsigDB_C2_CGP_BHp01 <- setReadable(gMsigDB_C2_CGP_BHp01, 'org.Mm.eg.db', 'ENTREZID')
gMsigDB_C2_CGP_BHp01_DF <- as.data.frame(gMsigDB_C2_CGP_BHp01)

gMsigDB_C2_CP_BHp20 <- GSEA(gene = RL1_GSEA, 
                            nPerm = 10000,
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.2,
                            TERM2GENE = mouse_msigDB_C2_CP,
                            verbose=FALSE)
gMsigDB_C2_CP_BHp20 <- setReadable(gMsigDB_C2_CP_BHp20, 'org.Mm.eg.db', 'ENTREZID')
gMsigDB_C2_CP_BHp20_DF <- as.data.frame(gMsigDB_C2_CP_BHp20)

gMsigDB_C2_CPb_BHp20 <- GSEA(gene = RL1_GSEA, 
                             nPerm = 10000,
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.2,
                             TERM2GENE = mouse_msigDB_C2_CPb,
                             verbose=FALSE)
gMsigDB_C2_CPb_BHp20 <- setReadable(gMsigDB_C2_CPb_BHp20, 'org.Mm.eg.db', 'ENTREZID')
gMsigDB_C2_CPb_BHp20_DF <- as.data.frame(gMsigDB_C2_CPb_BHp20)

gMsigDB_C3_TFT_BHp20 <- GSEA(gene = RL1_GSEA, 
                             nPerm = 10000,
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.2,
                             TERM2GENE = mouse_msigDB_C3_TFT,
                             verbose=FALSE)
gMsigDB_C3_TFT_BHp20 <- setReadable(gMsigDB_C3_TFT_BHp20, 'org.Mm.eg.db', 'ENTREZID')
gMsigDB_C3_TFT_BHp20_DF <- as.data.frame(gMsigDB_C3_TFT_BHp20)

gMsigDB_C6_BHp20 <- GSEA(gene = RL1_GSEA, 
                         nPerm = 10000,
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.2,
                         TERM2GENE = mouse_msigDB_C6,
                         verbose=FALSE)
gMsigDB_C6_BHp20 <- setReadable(gMsigDB_C6_BHp20, 'org.Mm.eg.db', 'ENTREZID')
gMsigDB_C6_BHp20_DF <- as.data.frame(gMsigDB_C6_BHp20)

gMsigDB_C7_BHp005 <- GSEA(gene = RL1_GSEA, 
                          nPerm = 10000,
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.005,
                          TERM2GENE = mouse_msigDB_C7,
                          verbose=FALSE)
gMsigDB_C7_BHp005 <- setReadable(gMsigDB_C7_BHp005, 'org.Mm.eg.db', 'ENTREZID')
gMsigDB_C7_BHp005_DF <- as.data.frame(gMsigDB_C7_BHp005)

# Count pathways meeting specified PValues and output statistics to dataframe
gSum1 <- as.data.frame(nrow(gKEGG1_BHp10_DF))
gSum2 <- as.data.frame(nrow(gKEGG1_all_DF))
gSum3 <- as.data.frame(nrow(gGO1_BP_BHp05_light_DF))
gSum4 <- as.data.frame(nrow(gGO1_MF_BHp10_light_DF))
gSum5 <- as.data.frame(nrow(gReactome1_BHp05_DF))
gSum6 <- as.data.frame(nrow(gMsigDB_H_BHp20_DF))
gSum7 <- as.data.frame(nrow(gMsigDB_H_all_DF))
gSum8 <- as.data.frame(nrow(gMsigDB_C2_CGP_BHp01_DF))
gSum9 <- as.data.frame(nrow(gMsigDB_C2_CP_BHp20_DF))
gSum10 <- as.data.frame(nrow(gMsigDB_C2_CPb_BHp20_DF))
gSum11 <- as.data.frame(nrow(gMsigDB_C3_TFT_BHp20_DF))
gSum12 <- as.data.frame(nrow(gMsigDB_C6_BHp20_DF))
gSum13 <- as.data.frame(nrow(gMsigDB_C7_BHp005_DF))
gSum14 <- as.data.frame(nrow(gReactome1_all))
gSummary <- as.data.frame(cbind(mget(ls(pattern = "gSum\\d+"))))

# Write DEG hypergeometric test results to tab delimited text
write.table(KEGG1_BHp20_DF, file=paste(Output_file_prefix, "KEGG1_BHp20.txt", sep="_"), sep='\t', quote=FALSE, row.names=FALSE)
write.table(KEGG1_all_DF, file=paste(Output_file_prefix, "KEGG1_all.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(KEGG1_M_all_DF, file=paste(Output_file_prefix, "KEGG1_M_all.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(GO1_BP_BHp001_DF, file=paste(Output_file_prefix, "GO1_BP_BHp001.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(GO1_MF_BHp05_DF, file=paste(Output_file_prefix, "GO1_MF_BHp05.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(Reactome1_BHp05_DF, file=paste(Output_file_prefix, "Reactome1_BHp05.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(Reactome1_all_DF, file=paste(Output_file_prefix, "Reactome1_BHp05.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(MsigDB_H_BHp20_DF, file=paste(Output_file_prefix, "MsigDB_H_BHp20.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(MsigDB_H_all_DF, file=paste(Output_file_prefix, "MsigDB_H_all.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(MsigDB_C2_CGP_BHp0001_DF, file=paste(Output_file_prefix, "MsigDB_C2_CGP_BHp0001.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(MsigDB_C2_CP_BHp20_DF, file=paste(Output_file_prefix, "MsigDB_C2_CP_BHp20.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(MsigDB_C2_CPb_BHp20_DF, file=paste(Output_file_prefix, "MsigDB_C2_CPb_BHp20.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(MsigDB_C3_TFT_BHp20_DF, file=paste(Output_file_prefix, "MsigDB_C3_TFT_BHp20.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(MsigDB_C6_BHp20_DF, file=paste(Output_file_prefix, "MsigDB_C6_BHp20.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(MsigDB_C7_BHp0001_DF, file=paste(Output_file_prefix, "MsigDB_C7_BHp0001.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(hSummary$V1, file=paste(Output_file_prefix, "clusterProfiler_Summary.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)

# Write transcriptome GSEA results to tab delimited text
write.table(gKEGG1_BHp10_DF, file=paste(GSEA_output_file_prefix, "gKEGG1_BHp10.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(gKEGG1_all_DF, file=paste(GSEA_output_file_prefix, "gKEGG1_all.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(gGO1_BP_BHp05_light_DF, file=paste(GSEA_output_file_prefix, "gGO1_BP_BHp05_light.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(gGO1_MF_BHp10_light_DF, file=paste(GSEA_output_file_prefix, "gGO1_MF_BHp10_light.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(gReactome1_BHp05_DF, file=paste(GSEA_output_file_prefix, "gReactome1_BHp05.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(gReactome1_all_DF, file=paste(GSEA_output_file_prefix, "gReactome1_BHp05.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(gMsigDB_H_all_DF, file=paste(GSEA_output_file_prefix, "gMsigDB_H_all.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(gMsigDB_C2_CGP_BHp01_DF, file=paste(GSEA_output_file_prefix, "gMsigDB_C2_CGP_BHp01.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(gMsigDB_C2_CP_BHp20_DF, file=paste(GSEA_output_file_prefix, "gMsigDB_C2_CP_BHp20.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(gMsigDB_C2_CPb_BHp20_DF, file=paste(GSEA_output_file_prefix, "gMsigDB_C2_CPb_BHp20.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(gMsigDB_C3_TFT_BHp20_DF, file=paste(GSEA_output_file_prefix, "gMsigDB_C3_TFT_BHp20.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(gMsigDB_C6_BHp20_DF, file=paste(GSEA_output_file_prefix, "gMsigDB_C6_BHp20.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(gMsigDB_C7_BHp005_DF, file=paste(GSEA_output_file_prefix, "gMsigDB_C7_BHp005.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)
write.table(gSummary$V1, file=paste(GSEA_output_file_prefix, "clusterProfiler_GSEA_Summary.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)

# Plot DGE hypergeometric test data and export to single file
pdf(file=paste(Output_file_prefix, "clusterProfiler_plots.pdf", sep="_"), width=8, heigh=8, title="clusterProfiler_plots")
dotplot(KEGG1_all, showCategory=20,font.size = 8, title = "KEGG1_all")
dotplot(KEGG1_M_all, showCategory=20,font.size = 8, title = "KEGG1_M_all")
dotplot(GO1_BP_BHp001, showCategory=20,font.size = 8, title = "GO1_BP_BHp001")
dotplot(GO1_MF_BHp05, showCategory=20,font.size = 8, title = "GO1_MF_BHp05")
dotplot(Reactome1_all, showCategory=20,font.size = 8, title = "Reactome1_BHp05")
dotplot(MsigDB_H_all, showCategory=20,font.size = 8, title = "gMsigDB_H_all")
dotplot(MsigDB_C2_CGP_BHp0001, showCategory=20,font.size = 8, title = "MsigDB_C2_CGP_BHp0001")
dotplot(MsigDB_C2_CP_BHp20, showCategory=20,font.size = 8, title = "MsigDB_C2_CP_BHp20")
dotplot(MsigDB_C2_CPb_BHp20, showCategory=20,font.size = 8, title = "MsigDB_C2_CPb_BHp20")
dotplot(MsigDB_C3_TFT_BHp20, showCategory=20,font.size = 8, title = "MsigDB_C3_TFT_BHp20")
dotplot(MsigDB_C6_BHp20, showCategory=20,font.size = 8, title = "MsigDB_C6_BHp20")
dotplot(MsigDB_C7_BHp0001, showCategory=20,font.size = 8, title = "MsigDB_C7_BHp0001")
emapplot(KEGG1_all, showCategory=20)+ggplot2::ggtitle('KEGG1_all')
emapplot(KEGG1_M_all, showCategory=20)+ggplot2::ggtitle('KEGG1_M_all')
emapplot(GO1_BP_BHp001, showCategory=20)+ggplot2::ggtitle('GO1_BP_BHp001')
emapplot(GO1_MF_BHp05, showCategory=20)+ggplot2::ggtitle('GO1_MF_BHp05')
emapplot(Reactome1_all, showCategory=20)+ggplot2::ggtitle('Reactome1_BHp05')
emapplot(MsigDB_H_all, showCategory=20)+ggplot2::ggtitle('MsigDB_H_all')
emapplot(MsigDB_C2_CGP_BHp0001, showCategory=20)+ggplot2::ggtitle('MsigDB_C2_CGP_BHp0001')
emapplot(MsigDB_C2_CP_BHp20, showCategory=20)+ggplot2::ggtitle('MsigDB_C2_CP_BHp20')
emapplot(MsigDB_C2_CPb_BHp20, showCategory=20)+ggplot2::ggtitle('MsigDB_C2_CPb_BHp20')
emapplot(MsigDB_C3_TFT_BHp20, showCategory=20)+ggplot2::ggtitle('MsigDB_C3_TFT_BHp20')
emapplot(MsigDB_C6_BHp20, showCategory=20)+ggplot2::ggtitle('MsigDB_C6_BHp20')
emapplot(MsigDB_C7_BHp0001, showCategory=20)+ggplot2::ggtitle('MsigDB_C7_BHp0001')
heatplot(KEGG1_all, showCategory=25)+ggplot2::ggtitle('KEGG1_all')
cnetplot(KEGG1_all, colorEdge = TRUE, node_label = FALSE, showCategory=5)+ggplot2::ggtitle('KEGG1_all')
upsetplot(KEGG1_all)
pmcplot(KEGG1_all$Description[1:10], 2010:2018, proportion=FALSE)+
  ggplot2::theme_classic() + ggtitle('gKEGG1_all')
dev.off()

# Plot transcriptome GSEA data and export to single file
pdf(file=paste(GSEA_output_file_prefix, "clusterProfiler_GSEA_plots.pdf", sep="_"), width=8, heigh=8, title="clusterProfiler_plots")
dotplot(gKEGG1_all, showCategory=20,font.size = 8, title = "gKEGG1_all")
dotplot(gGO1_BP_BHp05, showCategory=20,font.size = 8, title = "gGO1_BP_BHp05")
dotplot(gGO1_MF_BHp10, showCategory=20,font.size = 8, title = "gGO1_MF_BHp10")
dotplot(gReactome1_all, showCategory=20,font.size = 8, title = "gReactome1_BHp05")
dotplot(gMsigDB_H_BHp20, showCategory=20,font.size = 8, title = "gMsigDB_H_BHp20")
dotplot(gMsigDB_H_all, showCategory=20,font.size = 8, title = "gMsigDB_H_all")
dotplot(gMsigDB_C2_CGP_BHp01, showCategory=20,font.size = 8, title = "gMsigDB_C2_CGP_BHp01")
dotplot(gMsigDB_C2_CP_BHp20, showCategory=20,font.size = 8, title = "gMsigDB_C2_CP_BHp20")
dotplot(gMsigDB_C2_CPb_BHp20, showCategory=20,font.size = 8, title = "gMsigDB_C2_CPb_BHp20")
dotplot(gMsigDB_C3_TFT_BHp20, showCategory=20,font.size = 8, title = "gMsigDB_C3_TFT_BHp20")
dotplot(gMsigDB_C6_BHp20, showCategory=20,font.size = 8, title = "gMsigDB_C6_BHp20")
dotplot(gMsigDB_C7_BHp005, showCategory=20,font.size = 8, title = "gMsigDB_C7_BHp005")            
emapplot(gKEGG1_all, showCategory=20)+ggplot2::ggtitle('gKEGG1_all')
emapplot(gGO1_BP_BHp05, showCategory=20)+ggplot2::ggtitle('gGO1_BP_BHp05')
emapplot(gGO1_MF_BHp10, showCategory=20)+ggplot2::ggtitle('gGO1_MF_BHp10')
emapplot(gReactome1_all, showCategory=20)+ggplot2::ggtitle('gReactome1_BHp05')
emapplot(gMsigDB_H_all, showCategory=20)+ggplot2::ggtitle('gMsigDB_H_all')
emapplot(gMsigDB_C2_CGP_BHp01, showCategory=20)+ggplot2::ggtitle('gMsigDB_C2_CGP_BHp01')
emapplot(gMsigDB_C2_CP_BHp20, showCategory=20)+ggplot2::ggtitle('gMsigDB_C2_CP_BHp20')
emapplot(gMsigDB_C2_CPb_BHp20, showCategory=20)+ggplot2::ggtitle('gMsigDB_C2_CPb_BHp20')
emapplot(gMsigDB_C3_TFT_BHp20, showCategory=20)+ggplot2::ggtitle('gMsigDB_C3_TFT_BHp20')
emapplot(gMsigDB_C6_BHp20, showCategory=20)+ggplot2::ggtitle('gMsigDB_C6_BHp20')
emapplot(gMsigDB_C7_BHp005, showCategory=20)+ggplot2::ggtitle('gMsigDB_C7_BHp005')
gseaplot2(gKEGG1_all, geneSetID = 1:5, pvalue_table = TRUE, color = c("blue", "red", "green", "gold", "violet"), ES_geom = "line", base_size = 11, subplots = 1:3, rel_heights = c(2, 0.5, 1), title = "gKEGG1_all")
gseaplot2(gMsigDB_H_all, geneSetID = 1:5, pvalue_table = TRUE, color = c("blue", "red", "green", "gold", "violet"), ES_geom = "line", base_size = 11, subplots = 1:3, rel_heights = c(2, 0.5, 1), title = "gMsigDB_H_all")
gseaplot(gMsigDB_H_all, geneSetID = 1, title = gMsigDB_H_all$Description[1])
gseaplot(gMsigDB_H_all, geneSetID = 2, title = gMsigDB_H_all$Description[2])
gseaplot(gMsigDB_H_all, geneSetID = 3, title = gMsigDB_H_all$Description[3])
gseaplot(gMsigDB_H_all, geneSetID = 4, title = gMsigDB_H_all$Description[4])
gseaplot(gMsigDB_H_all, geneSetID = 5, title = gMsigDB_H_all$Description[5])
gsearank(gMsigDB_H_all, 1, title = gMsigDB_H_all[1, "Description"])
gsearank(gMsigDB_H_all, 2, title = gMsigDB_H_all[2, "Description"])
gsearank(gMsigDB_H_all, 3, title = gMsigDB_H_all[3, "Description"])
gsearank(gMsigDB_H_all, 4, title = gMsigDB_H_all[4, "Description"])
gsearank(gMsigDB_H_all, 5, title = gMsigDB_H_all[5, "Description"])
heatplot(gMsigDB_H_all)+ggplot2::ggtitle('gMsigDB_H_all')
cnetplot(gMsigDB_H_all, colorEdge = TRUE, node_label = FALSE, showCategory=5)+ggplot2::ggtitle('gMsigDB_H_all')
pmcplot(gKEGG1_all$Description[1:10], 2010:2018, proportion=FALSE)+
  ggplot2::theme_classic() + ggtitle('gKEGG1_all')
dev.off()

writeLines(capture.output(sessionInfo(package = NULL)), paste(GSEA_output_file_prefix, "clusterProfiler_Rinfo.txt",sep="_"))




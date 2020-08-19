params <-
list(isSlides = "no")

## ----setup_varManS2, include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo = TRUE,cache = TRUE,cache.lazy = FALSE, tidy = T)
# AsSlides <- TRUE
#
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressPackageStartupMessages(library(SNPlocs.Hsapiens.dbSNP144.GRCh37))
suppressPackageStartupMessages(library(ggplot2))
vcf <- readVcf("data/SAMN01882168_filt.vcf.gz","hg19")
#vcf <- readVcf("../data/SAMN01882168_filt.vcf.gz","hg19")
rd <- rowRanges(vcf)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# Genomic Variants (part 2)

---
"    
  )
  
}



## ----annoRS_varMan------------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)


## ----chrSub,eval=TRUE,tidy=FALSE,echo=TRUE------------------------------------
names(vcf)[1:2]
grepl(names(vcf),pattern = "chr1:")[1:2]
vcf_chr1 <- vcf[grepl(names(vcf),pattern = "chr1:")]
rd_chr1 <- rowRanges(vcf_chr1)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Annotating Variants

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Annotating Variants

---
"    
  )
  
}



## ----dbSNPv_varMan------------------------------------------------------------
all_snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
all_snps


## ----dbSNPv_varMan_load_fail,eval=FALSE---------------------------------------
## tar_chr <- as.vector(seqnames(rd_chr1)@values)
## my_snps <- snpsBySeqname(all_snps,c(tar_chr))

## ----dbSNPv_varMan_load_fail2,eval=TRUE,echo=FALSE,fig.align="center",out.width="75%"----
knitr::include_graphics("imgs/vcfMan_fig5r.png")


## ----check_seqLvl-------------------------------------------------------------
# seqlevels ~ rd
seqlevels(rd_chr1)
# seqlevels ~ dbSNP
seqlevels(all_snps)


## ----check_seqStyle-----------------------------------------------------------
seqlevelsStyle(rd_chr1)
seqlevelsStyle(all_snps)


## ----check_genome-------------------------------------------------------------
genome(rd_chr1)
genome(all_snps)


## ----dbSNPv_varMan_load-------------------------------------------------------
tar_chr <- as.vector(seqnames(rd_chr1)@values)
tar_chr <- gsub("chr","",tar_chr)
tar_chr[grepl(tar_chr,pattern = "M")] <- "MT"
my_snps <- snpsBySeqname(all_snps,c(tar_chr))
my_snps[1:2]


## ----change_seqlvl------------------------------------------------------------
# change seqlevelsStyle
seqlevelsStyle(my_snps) <- "UCSC"
# change genome
genome(my_snps) <- "hg19"


## ----makeTab_varMan-----------------------------------------------------------
snp_ID <- data.frame(
  posIDX=paste0(seqnames(my_snps),":",pos(my_snps)),
  rsID=my_snps$RefSNP_id,stringsAsFactors = FALSE)
head(snp_ID)


## ----dbSNPV_varMan_tbl_1------------------------------------------------------
matV1 <- data.frame(Variant=names(rd_chr1),stringsAsFactors = FALSE)
matV1[1:2,]


## ----dbSNPV_varMan_tbl_2------------------------------------------------------
matV1$chromosome <- gsub("(.*):(.*)_(.*)/(.*)","\\1",matV1$Variant)
matV1$start <- gsub("(.*):(.*)_(.*)/(.*)","\\2",matV1$Variant)
matV1$end <- gsub("(.*):(.*)_(.*)/(.*)","\\2",matV1$Variant)
matV1$ref_allele <- gsub("(.*):(.*)_(.*)/(.*)","\\3",matV1$Variant)
matV1$alt_allele <- gsub("(.*):(.*)_(.*)/(.*)","\\4",matV1$Variant)
matV1$posIDX <- gsub("(.*)_(.*)","\\1",matV1$Variant)
matV1[1:2,]


## ----dbSNPV_varMan_tbl_3,tidy=FALSE-------------------------------------------
matS <- merge(matV1,snp_ID,all.x=TRUE,by="posIDX")
matS <- dplyr::select(matS,-posIDX)
matS[1:2,]


## ----dbSNPV_varMan_muCt,tidy=FALSE--------------------------------------------
taC2 <- table(!is.na(matS$rsID))
taC2_dat <- as.data.frame(taC2)
taC2


## ----dbSNPv_varMan_muCt_disp1,tidy=FALSE,eval=FALSE,echo=TRUE-----------------
## ggplot(taC2_dat,aes(x=Var1,y=Freq,fill=Var1))+
##   geom_bar(stat='Identity')+
##   labs(x="",y="Counts",fill="in_dbSNP")+
##   theme(legend.position = "none")


## ----dbSNPv_varMan_muCt_disp2,tidy=FALSE,echo=FALSE,eval=TRUE,fig.align="center"----
ggplot(taC2_dat,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat='Identity')+
  labs(x="",y="Counts",fill="in_dbSNP")+
  theme(legend.position = "none")


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Variants and Amino Acid Changes

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Variants and Amino Acid Changes

---
"    
  )
  
}



## ----aaCh_varMan--------------------------------------------------------------
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


## ----aaCh_varMan_txdb---------------------------------------------------------
txdb


## ----aaCh_varMan_pred---------------------------------------------------------
coding <- predictCoding(vcf_chr1, txdb, seqSource=Hsapiens)


## ----aaCh_varMan_pres2--------------------------------------------------------
coding[1]


## ----aaCh_varMan_frame,eval=TRUE,tidy=FALSE,echo=TRUE-------------------------
matA <- data.frame(Variant=names(coding),
                   chromosome=seqnames(coding),
                   start=start(coding),end=end(coding),
                   ref_allele=as.character(coding$REF),
                   alt_allele=unlist(lapply(lapply(
                     coding$ALT,`[[`,1),as.character)),
                   GeneID=coding$GENEID,
                   TxID=coding$TXID,
                   Protein_posi=unlist(lapply(lapply(
                     coding$PROTEINLOC,`[[`,1),as.integer)),
                   ref_AA=as.character(coding$REFAA),
                   alt_AA=as.character(coding$VARAA),
                   Type=coding$CONSEQUENCE)
matA$aaChange <- paste0("p.",matA$ref_AA,matA$Protein_posi,matA$alt_AA)
matA <- dplyr::select(matA,-Protein_posi,-ref_AA,-alt_AA)


## ----aaCh_varMan_tbl----------------------------------------------------------
matA[1:2,]


## ----aaCh_varMan_muCt---------------------------------------------------------
var_in_coding <- data.frame(varName=names(vcf_chr1),
                            in_coding=names(vcf_chr1)
                            %in% matA$Variant,
                            stringsAsFactors = FALSE)
table(var_in_coding$in_coding)


## ----aaCh_varMan_muType-------------------------------------------------------
taC <- table(matA$Type)
taC_dat <- as.data.frame(taC)
taC


## ----aaCh_varMan_muType_disp1,tidy=FALSE,echo=TRUE,eval=FALSE-----------------
## ggplot(taC_dat,aes(x=Var1,y=Freq,fill=Var1))+
##   geom_bar(stat='Identity')+
##   labs(x="",y="Counts",fill="")+
##   theme(legend.position = "none")


## ----aaCh_varMan_muType_disp2,tidy=FALSE,echo=FALSE,eval=TRUE,fig.align="center"----
ggplot(taC_dat,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat='Identity')+
  labs(x="",y="Counts",fill="")+
  theme(legend.position = "none")


## ----comb_varMan--------------------------------------------------------------
matS$GeneID <- matA$GeneID[match(matS$Variant,matA$Variant)]
matS$AAChange <- matA$GeneID[match(matS$Variant,matA$Variant)]
matS[1:2,]


params <-
list(isSlides = "no")

## ----setup_varMan, include=FALSE----------------------------------------------
knitr::opts_chunk$set(echo = TRUE,cache = TRUE,cache.lazy = FALSE, tidy = T)
# AsSlides <- TRUE
#
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressPackageStartupMessages(library(SNPlocs.Hsapiens.dbSNP144.GRCh37))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GenomicFeatures))


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# Genomic Variants (part 1)

---
"    
  )
  
}



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# VCF files

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# VCF files

---
"    
  )
  
}



## ----varLoad_intro,echo=FALSE,out.width = "75%",fig.align="center"------------
knitr::include_graphics("imgs/vcfMan_fig2.png")


## ----varLoad_varMan-----------------------------------------------------------
library(VariantAnnotation)
vcf <- readVcf("data/SAMN01882168_filt.vcf.gz","hg19")
vcf


## ----gInfo_varMan-------------------------------------------------------------
header(vcf)


## ----sample_varMan------------------------------------------------------------
sampleID <- samples(header(vcf))
sampleID


## ----metaOV_META_varMan-------------------------------------------------------
meta(header(vcf))


## ----meta_META_varMan---------------------------------------------------------
# File format
meta(header(vcf))$fileformat
# Source used for variant calling
meta(header(vcf))$source


## ----meta_contig_varMan-------------------------------------------------------
meta(header(vcf))$contig


## ----range_varMan-------------------------------------------------------------
rd <- rowRanges(vcf)
rd[1:2]


## ----range_varMan_posi--------------------------------------------------------
as.vector(seqnames(rd)[1:2]) # Chromosome
start(rd)[1:2] # Start position
end(rd)[1:2] # End position


## ----range_varMan_baseInfo_Ref1-----------------------------------------------
refBase <- ref(vcf)
refBase[1:2]


## ----range_varMan_baseInfo_Ref3-----------------------------------------------
refBase <- as.character(refBase)
refBase[1:2]


## ----range_varMan_baseInfo_Alt1-----------------------------------------------
altBase <- alt(vcf)
alt(vcf)[1:2] 


## ----range_varMan_baseInfo_Alt3-----------------------------------------------
# get the 1st vector of the list
altBase <- lapply(altBase,`[[`,1) 
altBase[1:2]


## ----range_varMan_baseInfo_Alt4-----------------------------------------------
altBase <- unlist(lapply(altBase,as.character)) 
altBase[1:2]


## ----info_varMan--------------------------------------------------------------
info(header(vcf))[1:2,]


## ----info_varMan_disp---------------------------------------------------------
info(vcf)[1:2,]


## ----geno_varMan--------------------------------------------------------------
geno(header(vcf))[1:2,]


## ----genoGT_varMan------------------------------------------------------------
paste0("GT: ",geno(header(vcf))[1,3])
matGT <- geno(vcf)$GT
matGT[1:2,]


## ----genoGT_varMan_tbl--------------------------------------------------------
tbl <- table(geno(vcf)$GT)
tbl_dat <- as.data.frame(tbl)
tbl


## ----genoGT_varMan_disp1,echo=TRUE,tidy=FALSE,eval=FALSE----------------------
## ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
##   geom_bar(stat='Identity')+
##   labs(x="",y="Counts",fill="")+
##   theme_classic()


## ----genoGT_varMan_disp2,echo=FALSE,tidy=FALSE,eval=TRUE,fig.align="center"----
ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat='Identity')+
  labs(x="",y="Counts",fill="")+theme_classic()


## ----genoDP_varMan------------------------------------------------------------
paste0("DP: ",geno(header(vcf))[3,3])
matDP <- geno(vcf)$DP
matDP[1:2,]


## ----genoDP_varMan_dist1,eval=FALSE,echo=TRUE,tidy=FALSE----------------------
## summary(as.vector(matDP))


## ----genoDP_varMan_dist2,eval=TRUE,echo=FALSE,results=TRUE--------------------
summary(as.vector(matDP))


## ----genoDP_varMan_distPres1,fig.align="center",eval=FALSE,echo=TRUE,tidy=FALSE----
## ggplot(as.data.frame(matDP),aes(x=SAMN01882168))+geom_histogram()+
##   labs(x="",y="Counts")+
##   scale_x_log10()+
##   theme_classic()


## ----genoDP_varMan_distPres2,fig.align="center",eval=TRUE,echo=FALSE----------
ggplot(as.data.frame(matDP),aes(x=SAMN01882168))+geom_histogram()+
  labs(x="",y="Counts")+
  scale_x_log10()+
  theme_classic()


## ----genoGQ_varMan------------------------------------------------------------
paste0("GQ: ",geno(header(vcf))[4,3])
matGQ <- geno(vcf)$GQ
matGQ[1:2,]


## ----genoGQ_varMan_dist1,eval=FALSE,echo=TRUE,tidy=FALSE----------------------
## summary(as.vector(matGQ))


## ----genoGQ_varMan_dist2,eval=TRUE,echo=FALSE,tidy=FALSE,results=TRUE---------
summary(as.vector(matGQ))


## ----genoGQ_varMan_distPres1,fig.align="center",eval=FALSE,echo=TRUE,tidy=FALSE----
## ggplot(as.data.frame(matGQ),aes(x=SAMN01882168))+geom_histogram()+
##   labs(x="",y="Counts")+
##   scale_x_log10()+
##   theme_classic()


## ----genoGQ_varMan_distPres2,fig.align="center",eval=TRUE,echo=FALSE,warning=FALSE----
ggplot(as.data.frame(matGQ),aes(x=SAMN01882168))+geom_histogram()+
  labs(x="",y="Counts")+scale_x_log10()+theme_classic()


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Gathering variant information

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Gathering variant information

---
"    
  )
  
}



## ----gatGT_info1,echo=TRUE,tidy=FALSE-----------------------------------------
var_1 <- rownames(geno(vcf)$GT)[
  geno(vcf)$GT=="0/1" | 
    geno(vcf)$GT=="1/1"]


## ----gatGT_info2,echo=TRUE,tidy=FALSE-----------------------------------------
varTab1 <- data.frame(variant=names(rd)[names(rd) %in% var_1],
                      chr=as.vector(seqnames(rd)[names(rd) %in% var_1]),
                      start=start(rd)[names(rd) %in% var_1],
                      end=end(rd)[names(rd) %in% var_1],
                      stringsAsFactors = FALSE)


## ----gatGT_info3,echo=TRUE,tidy=FALSE-----------------------------------------
ref_base <- ref(vcf)[rownames(vcf) %in% var_1]
ref_base[1:2]
varTab1$refBase <- as.character(ref_base)


## ----gatGT_info4,echo=TRUE,tidy=FALSE-----------------------------------------
alt_base <- lapply(alt(vcf)[rownames(vcf) %in% var_1],`[[`,1)
alt_base[1]
alt_base <- lapply(alt_base,as.character)
alt_base[1]
varTab1$altBase <- unlist(alt_base)


## ----gatGT_info5,echo=TRUE,tidy=FALSE-----------------------------------------
adCount <- geno(vcf)$AD[rownames(geno(vcf)$AD) %in% var_1]
adCount[1]
varTab1$refCount <- unlist(lapply(adCount,`[[`,1))
varTab1$altCount <- unlist(lapply(adCount,`[[`,2))


## ----gatGT_info6,echo=TRUE,tidy=FALSE-----------------------------------------
varTab1$genoType <- geno(vcf)$GT[rownames(geno(vcf)$GT) %in% var_1]
varTab1$gtQuality <- geno(vcf)$GQ[rownames(geno(vcf)$GQ) %in% var_1]


## ----gatGT_info7,echo=TRUE,tidy=FALSE-----------------------------------------
var_2 <- rownames(geno(vcf)$GT)[geno(vcf)$GT=="1/2"]
varTab2 <- data.frame(variant=names(rd)[names(rd) %in% var_2],
                      chr=as.vector(seqnames(rd)[names(rd) %in% var_2]),
                      start=start(rd)[names(rd) %in% var_2],
                      end=end(rd)[names(rd) %in% var_2],
                      refBase=unlist(lapply(lapply(
                        alt(vcf)[rownames(vcf) %in% var_2],`[[`,1),as.character)),
                      altBase=unlist(lapply(lapply(
                        alt(vcf)[rownames(vcf) %in% var_2],`[[`,2),as.character)),
                      refCount=unlist(lapply(
                        geno(vcf)$AD[rownames(geno(vcf)$AD) %in% var_2],`[[`,2)),
                      altCount=unlist(
                        lapply(geno(vcf)$AD[rownames(geno(vcf)$AD) %in% var_2],`[[`,3)),
                      genoType=geno(vcf)$GT[rownames(geno(vcf)$GT) %in% var_2],
                      gtQuality=geno(vcf)$GQ[rownames(geno(vcf)$GQ) %in% var_2],
                      stringsAsFactors = FALSE)


## ----gatGT_merge--------------------------------------------------------------
varTab <- rbind(varTab1,varTab2)
varTab[1:2,]


## ----mutType_gen--------------------------------------------------------------
# differentiate SNP/INS/DEL/Others
for(k in 1:length(varTab$variant)){
  if(width(varTab$refBase[k]) < width(varTab$altBase[k])){
    varTab$mutType[k] <- "INS"
  }else if(width(varTab$refBase[k]) > width(varTab$altBase[k])){
    varTab$mutType[k] <- "DEL"
  }else if(width(varTab$refBase[k])==1&width(varTab$altBase[k])==1){
    varTab$mutType[k] <- "SNP"
  }else{
    varTab$mutType[k] <- "Others"}}
#
tbl <- table(varTab$mutType)
tbl_dat <- as.data.frame(tbl)
tbl


## ----mutType_pres1,echo=TRUE,eval=FALSE,tidy=FALSE,fig.align="center"---------
## ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
##   geom_bar(stat = 'identity')+
##   labs(x="",y="Mutations",fill="")+
##   theme_classic()


## ----mutType_pres2,echo=FALSE,eval=TRUE,fig.align="center"--------------------
ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+geom_bar(stat = 'identity')+
  labs(x="",y="Mutations",fill="")+theme_classic()


## ----TiTb_gen,tidy=FALSE------------------------------------------------------
# Transition (Ti)
ti <- c("A>G","G>A","C>T","T>C")
# Transveersion (Tv)
tv <- c("A>T","A>C","G>T","G>C","C>A","C>G","T>A","T>G")
varTab$nuSub <- paste0(varTab$refBase,">",varTab$altBase)
varTab$TiTv[varTab$nuSub %in% ti] <- "Ti"
varTab$TiTv[varTab$nuSub %in% tv] <- "Tv"
varTab[1:2,]


## ----TiTv_pres_nuSub----------------------------------------------------------
varX <- varTab[varTab$mutType=="SNP",]
tbl <- table(varX$nuSub)
tbl_dat <- as.data.frame(tbl)
tbl


## ----TiTv_pres_nuSubPres1,eval=FALSE,tidy=FALSE,echo=TRUE,fig.align="center"----
## ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
##   geom_bar(stat = 'identity')+
##   labs(x="",y="Mutations",fill="")+
##   theme(legend.position = "none")


## ----TiTv_pres_nuSubPres2,eval=TRUE,echo=FALSE,tidy=FALSE,fig.align="center"----
ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+geom_bar(stat = 'identity')+
  labs(x="",y="Mutations",fill="")+theme(legend.position = "none")


## ----TiTv_pres_TiTv1----------------------------------------------------------
tbl <- table(varX$TiTv)
tbl_dat <- as.data.frame(tbl)
tbl


## ----TiTv_pres_TiTv2,echo=TRUE,eval=FALSE,tidy=FALSE--------------------------
## ggplot(as.data.frame(table(varX$TiTv)),aes(x=Var1,y=Freq,fill=Var1))+
##   geom_bar(stat = 'identity')+labs(x="",y="Mutations",fill="")+
##   theme(legend.position = "none")


## ----TiTv_pres_TiTv=3,echo=FALSE,eval=TRUE,tidy=FALSE,fig.align="center"------
ggplot(as.data.frame(table(varX$TiTv)),aes(x=Var1,y=Freq,fill=Var1))+geom_bar(stat = 'identity')+labs(x="",y="Mutations",fill="")+theme(legend.position = "none")


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Trinucleotide motif analysis

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Trinucleotide motif analysis

---
"    
  )
  
}



## ----motif_load_advAn---------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicFeatures)
library(stringr)
#
rd_idx<-str_split(names(rd),"_", simplify=T)
rd_idx[1:2,]
rd_sub <- rd[rd_idx[,2]=="C/T"]
rd_sub[1:2,]


## ----motif_seqExt_advAn,tidy=FALSE--------------------------------------------
rd_sub$triNu <- getSeq(Hsapiens,
                 seqnames(rd_sub),
                 start=start(rd_sub)-1,
                 end=end(rd_sub)+1)
rd_sub[1:2]


## ----motif_seqPat_advan-------------------------------------------------------
tbl <- table(rd_sub$triNu)
tbl_dat <- as.data.frame(tbl)
tbl


## ----motif_seqPat_advan2,eval=FALSE,echo=TRUE,tidy=FALSE----------------------
## ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
##   geom_bar(stat='identity')+
##   labs(x="",y="Variants",fill="")+
##   theme(legend.position = "none")


## ----motif_seqPat_advan3,eval=TRUE,echo=FALSE,tidy=FALSE,fig.align="center"----
ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat='identity')+
  labs(x="",y="Variants",fill="")+
  theme(legend.position = "none")


## ----motif_ApoTar_advan1------------------------------------------------------
# TCW: TCA/TCT
tbl_dat$APOBEC_target <- tbl_dat$Var1 %in% c("TCA","TCT")
tbl_dat[1:2]
# collapse by APOBEC_target
apobec_dat <- aggregate(Freq ~ APOBEC_target,tbl_dat,FUN=sum,na.rm=TRUE)
apobec_dat


## ----motif_ApoTar_advan2,eval=FALSE,echo=TRUE,tidy=FALSE,fig.align="center"----
## ggplot(apobec_dat,aes(x=APOBEC_target,y=Freq,fill=APOBEC_target))+
##   geom_bar(stat='identity')+
##   labs(x="",y="Variants",fill="")+
##   theme(legend.position = "none")


## ----motif_ApoTar_advan3,eval=TRUE,echo=FALSE,tidy=FALSE,fig.align="center"----
ggplot(apobec_dat,aes(x=APOBEC_target,y=Freq,fill=APOBEC_target))+
  geom_bar(stat='identity')+
  labs(x="",y="Variants",fill="")+
  theme(legend.position = "none")


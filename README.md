DPA: Differential presentation analysis to prioritize etiologically relevant T cell epitopes across autoimmune diseases
===============================================

The 'DPA' package provides a streamlined framework of differential presentation analysis for a given set of HLA-I- and HLA-II-restricted T cell epitope sequences.

Concept
------------------------
The underlying concept of the DPA framework is simple. Several genetic associations between HLA alleles and various autoimmune diseases have been identified from previous penome-wide association studies. Given the biological function of HLA molecules, here it is assumed that HLA molecules whose alleles are genetically associated with disease predisposition presumably present either pathogenic epitopes more or protective epitopes less, and likewise, HLA molecules whose alleles are genetically associated with disease protection presumably present either pathogenic epitopes less or protective epitopes more. 

Method
------------------------
All four-digit HLA alleles significantly associated with a specific autoimmune disease X were used for HLA binding prediction of an epitope Y that has been studied in the context of X. NetMHCpan 4.0 and NetMHCIIpan 3.2 were utilized for HLA binding prediction with default parameter sets. Predicted percentile rank was chosen as a metric for the strength of epitope binding because this metric is not affected by inherent bias of specific HLA molecules towards higher or lower mean predicted affinities and thus allows a direct comparison between different HLA molecules. The highest values of sign-inverted log10-transformed percentile ranks, corresponding to the lowest percentile ranks, among predisposing and protective HLA molecules were adopted for differential presentation analysis. Differential presentation index (DPI) was defined as the transformed value of predisposing alleles subtracted by that of protective alleles. DPI is disease-dependent because the sets of predisposing and protective alleles vary between diseases. Epitopes were then categorized in a binary fashion; epitopes with DPI of higher than 0.5 and lower than -0.5 were considered putatively disease-predisposing and putatively disease-protective. Note that epitopes predicted not to bind to any of the disease-associated alleles, with the thresholds being 2% and 10% for HLA-I and HLA-II binding prediction, respectively, were excluded from this binary categorization. 

Installation
------------------------
Install the latest version from [GitHub](https://github.com/masato-ogishi/DPA) as follows:
``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("masato-ogishi/DPA")
```
-   You might be prompted to install some packages before installling DPA. Follow the message(s).

Usage
------------------
0. Working environment
``` r
library(tidyverse)
library(data.table)
library(DPA)
```
1. Datasets
-   The following datasets are provided along with the package.
``` r
SummaryDF_Disease_HLA          ## A summary of associations between autoimune diseases and HLA alleles
SummaryDF_Disease_HLAString    ## A summary of HLA allele strings associated with autoimmune diseases (used for NetMHC binding prediction)
SummaryDF_Disease_Peptide      ## A summary of T cell epitope sequences collected from IEDB studied in the context of autoimmune diseases
autoimmunityDT                 ## A summary of autoimmunity-associated T cell epitope sequences with full analysis results
```
2. NetMHC prediction
-   Scripts and source FASTA files for NetMHC prediction can be automatically generated as follows.
``` r
peptide_hlastring_pairs <- merge(SummaryDF_Disease_Peptide[,.(Disease,Peptide,HLA_Class)], SummaryDF_Disease_HLAString)
NetMHCScript_HLAI(
  peptide_hlastring_pairs[HLA_Class=="ClassI",],
  outDir="./NetMHC/",
  outputFileName="NetMHCScript_HLA-I.txt"   ## For NetMHCpan 4.0
)
NetMHCScript_HLAII(
  peptide_hlastring_pairs[HLA_Class=="ClassII",],
  outDir="./NetMHC/",
  outputFileName="NetMHCScript_HLA-II.txt"  ## For NetMHCIIpan 3.2
)
```
3. Differential presentation analysis
-   Predicted binding to multiple disease-associated HLA molecules is summarized into a unidimentional scale, termed "differential presentation index (DPI)". Peptides are then categorized into "putatively disease-predisposing" and "putatively disease-protective" based on their DPI values accordingly. 
``` r
dt_HLAI <- differentialPresentationAnalysis(
  list.files(pattern="HLAI_.+.xls", path="./NetMHC/", full.names=T), 
  hlaClass="ClassI", rankPerThr=2, DPIThr=0.5
) ## %Rank threshold is set 2%
dt_HLAI <- merge(dt_HLAI, SummaryDF_Disease_Peptide[HLA_Class=="ClassI",], all.x=T, all.y=F)
dt_HLAII <- differentialPresentationAnalysis(
  list.files(pattern="HLAII_.+.xls", path="./NetMHC/", full.names=T), 
  hlaClass="ClassII", rankPerThr=10, DPIThr=0.5
) ## %Rank threshold is set 10%
dt_HLAII <- merge(dt_HLAII, SummaryDF_Disease_Peptide[HLA_Class=="ClassII"], all.x=T, all.y=F)
dt_autoimmunity <- rbind(dt_HLAI, dt_HLAII)
```
4. Similarity-to-self analysis
-   Similarity-to-self is determined through pairwise sequence alignment of the epitope and the entire human proteome.  
``` r
proteome_Human <- UniProt_Proteome_Import("UniProt_HomoSapiens_UP000005640.fasta")
UPMatch <- UniProt_Proteome_ExactMatch(dt_autoimmunity$"Peptide", proteome_Human)
dt_autoimmunity[,UniProtMatch:=UPMatch]
dt_autoimmunity[,Origin:="S"]                   ## Self
dt_autoimmunity[UniProtMatch=="",Origin:="NS"]  ## Non-self
dt_autoimmunity[,Origin:=factor(Origin, levels=c("S","NS"))]
UPAlScores <- UniProt_Proteome_AlignmentScores(dt_autoimmunity$"Peptide", proteome_Human)
dt_autoimmunity[,UniProtAlScore_Max:=UPAlScores$"Highest"]  ## Reflecting the "best-match" in the human proteome
dt_autoimmunity[,UniProtAlScore_Ave:=UPAlScores$"Average"]
dt_autoimmunity[,UniProtAlScore_Min:=UPAlScores$"Lowest"]

View(dt_autoimmunity)
```

Reference
------------------------
Ogishi, M. (2019) "Prioritizing putatively etiological T cell epitopes across autoimmune diseases." bioRxiv.


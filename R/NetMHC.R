#' NetMHC analysis utility.
#'
#' @param dt A data.table of peptides, associated diseases and HLA strings. Note: NetMHCpan and NetMHCIIpan only accept peptides of 8~14aa and 9~32aa, respectively. Our script-generating functions remove peptides of other lengths intetrnally.
#' @param hlaDigit The digit of HLA alleles from PheWAS to be included in the NetMHC prediction. Either "Two", "Four", or "TwoAndFour". Note: all possible matched four-digit HLA alleles are used in place of two-digits HLA alleles.
#' @param outDir A path to the output directory for NetMHC analysis.
#' @param outputFileName A file name for the generated NetMHC script.
#' @param fileName A NetMHC result file to be parsed. It should be named like <Class>_<Disease>_<Length>. Eg. HLAI_AS_9mer.xls
#' @export
#' @rdname NetMHC
#' @name NetMHC
NetMHCScript_HLAI <- function(dt, hlaDigit="Four", outDir="./NetMHC/", outputFileName="NetMHCScript_HLA-I.txt"){
  dir.create(outDir, showWarnings=F, recursive=T)
  script <- function(peptideSet, diseaseName, hlaString){
    peptideLength <- unique(nchar(peptideSet))
    seqinr::write.fasta(as.list(peptideSet), peptideSet, paste0(outDir, "HLAI_", diseaseName, "_", peptideLength, "mer.fasta"), as.string=T)
    if(nchar(hlaString)>1000){
      hla <- unlist(strsplit(hlaString, ","))
      hla <- BBmisc::chunk(hla, n.chunks=ceiling(nchar(hlaString)/1000))
      sapply(1:length(hla), function(i){
        paste0("./netMHCpan -BA -a ", paste0(hla[[i]], collapse=","), " -l ", peptideLength, " -xls -xlsfile ./Peptides/HLAI_", diseaseName, "_", peptideLength, "mer_", i, ".xls -f ./Peptides/HLAI_", diseaseName, "_", peptideLength, "mer.fasta")
      })
    }else{
      paste0("./netMHCpan -BA -a ", hlaString, " -l ", peptideLength, " -xls -xlsfile ./Peptides/HLAI_", diseaseName, "_", peptideLength, "mer.xls -f ./Peptides/HLAI_", diseaseName, "_", peptideLength, "mer.fasta")
    }
  }
  dt[,Peptide_Length:=nchar(Peptide)]
  dt <- dt[Peptide_Length %in% 8:14,]
  switch(hlaDigit,
         "Two"=dt[,HLA_String_NetMHC:=HLA_String_TwoDigits_AllFourDigits],
         "Four"=dt[,HLA_String_NetMHC:=HLA_String_FourDigits],
         "TwoAndFour"=dt[,HLA_String_NetMHC:=HLA_String_TwoAndFourDigits],
         return("hlaDigit should be either \"Two\", \"Four\", or \"TwoAndFour\"!")
  )
  res <- pbapply::pblapply(
    split(dt, by=c("Peptide_Length", "Disease")),
    function(d){script(d$"Peptide", d$"Disease"[1], d$"HLA_String_NetMHC"[1])}
  ) %>% unlist() %>% sort() %>% as.data.frame()
  write.table(res, paste0(outDir, outputFileName), row.names=F, col.names=F, quote=F)
  return(paste0(outDir, outputFileName))
}

#' @export
#' @rdname NetMHC
#' @name NetMHC
NetMHCScript_HLAII <- function(dt, hlaDigit="Four", outDir="./NetMHC/", outputFileName="NetMHCScript_HLA-II.txt"){
  dir.create(outDir, showWarnings=F, recursive=T)
  script <- function(peptideSet, diseaseName, hlaString){
    peptideLength <- unique(nchar(peptideSet))
    seqinr::write.fasta(as.list(peptideSet), peptideSet, paste0(outDir, "HLAII_", diseaseName, "_", peptideLength, "mer.fasta"), as.string=T)
    paste0("./netMHCIIpan -a ", hlaString, " -length ", peptideLength, " -xls -xlsfile ./Peptides/HLAII_", diseaseName, "_", peptideLength, "mer.xls -f ./Peptides/HLAII_", diseaseName, "_", peptideLength, "mer.fasta")
  }
  dt[,Peptide_Length:=nchar(Peptide)]
  dt <- dt[Peptide_Length %in% 9:32,]
  switch(hlaDigit,
         "Two"=dt[,HLA_String_NetMHC:=HLA_String_TwoDigits_AllFourDigits],
         "Four"=dt[,HLA_String_NetMHC:=HLA_String_FourDigits],
         "TwoAndFour"=dt[,HLA_String_NetMHC:=HLA_String_TwoAndFourDigits],
         return("hlaDigit should be either \"Two\", \"Four\", or \"TwoAndFour\"!")
  )
  res <- pbapply::pbsapply(
    split(dt, by=c("Peptide_Length", "Disease")),
    function(d){script(d$"Peptide", d$"Disease"[1], d$"HLA_String_NetMHC"[1])}
  ) %>% sort() %>% as.data.frame()
  write.table(res, paste0(outDir, outputFileName), row.names=F, col.names=F, quote=F)
  return(paste0(outDir, outputFileName))
}

#' @export
#' @rdname NetMHC
#' @name NetMHC
NetMHC_Import_HLAI <- function(fileName){
  diseaseName <- unlist(stringr::str_split(basename(fileName), "_"))[2]
  hlaList <- unlist(data.table::fread(fileName, nrows=1, header=F)[1,])
  hlaList <- hlaList[!is.na(hlaList)]
  hlaList <- stringr::str_replace_all(hlaList, "HLA-", "")
  hlaList <- stringr::str_replace_all(hlaList, ":", "")
  stringr::str_sub(hlaList, 2, 1) <- "_"
  res <- data.table::fread(fileName, skip=1)[Pos==0,]
  colnames(res)[seq(4, 3+5*length(hlaList))] <- paste0(rep(hlaList, each=5), "&", colnames(res)[seq(4, 3+5*length(hlaList))])
  res <- res %>%
    dplyr::select(Peptide, dplyr::matches("nM"), dplyr::matches("Rank")) %>%
    tidyr::gather(Variable, Value, -Peptide) %>%
    tidyr::separate(col="Variable", into=c("HLA", "Variable"), sep="&") %>%
    dplyr::mutate(Variable=dplyr::if_else(Variable=="nM", "Affinity", "Rank")) %>%
    dplyr::mutate(Disease=diseaseName, HLA_Class="ClassI") %>%
    tidyr::spread(Variable, Value)
  return(res)
}

#' @export
#' @rdname NetMHC
#' @name NetMHC
NetMHC_Import_HLAII <- function(fileName){
  diseaseName <- unlist(stringr::str_split(basename(fileName), "_"))[2]
  hlaList <- unlist(data.table::fread(fileName, nrows=1, header=F)[1,])
  hlaList <- hlaList[!is.na(hlaList)]
  hlaList <- stringr::str_replace_all(hlaList, "HLA-", "")
  hlaList <- stringr::str_replace_all(hlaList, ":", "")
  res <- data.table::fread(fileName, skip=1)[Pos==1,]
  colnames(res)[seq(4, 3+3*length(hlaList))] <- paste0(rep(hlaList, each=3), "&", colnames(res)[seq(4, 3+3*length(hlaList))])
  res <- res %>%
    dplyr::select(Peptide, dplyr::matches("nM"), dplyr::matches("Rank")) %>%
    tidyr::gather(Variable, Value, -Peptide) %>%
    tidyr::separate(col="Variable", into=c("HLA", "Variable"), sep="&") %>%
    dplyr::mutate(Variable=dplyr::if_else(Variable=="nM", "Affinity", "Rank")) %>%
    dplyr::mutate(Disease=diseaseName, HLA_Class="ClassII") %>%
    tidyr::spread(Variable, Value) %>%
    data.table::as.data.table()
  res.DR <- res[stringr::str_detect(HLA, "DR"),] %>%
    dplyr::mutate(HLA_Type="DR")
  res.DPA <- res[stringr::str_detect(HLA, "DP"),] %>%
    tidyr::separate(col="HLA", into=c("HLA", "HLA_B"), sep="-") %>%
    dplyr::group_by(Peptide, HLA, Disease, HLA_Class) %>%
    dplyr::summarise(Affinity=median(Affinity), Rank=median(Rank)) %>%
    dplyr::mutate(HLA_Type="DP")
  stringr::str_sub(res.DPA$"HLA", 5, 4) <- "_"
  res.DPB <- res[stringr::str_detect(HLA, "DP"),] %>%
    tidyr::separate(col="HLA", into=c("HLA_A", "HLA"), sep="-") %>%
    dplyr::group_by(Peptide, HLA, Disease, HLA_Class) %>%
    dplyr::summarise(Affinity=median(Affinity), Rank=median(Rank)) %>%
    dplyr::mutate(HLA_Type="DP")
  stringr::str_sub(res.DPB$"HLA", 5, 4) <- "_"
  res.DQA <- res[stringr::str_detect(HLA, "DQ"),] %>%
    tidyr::separate(col="HLA", into=c("HLA", "HLA_B"), sep="-") %>%
    dplyr::group_by(Peptide, HLA, Disease, HLA_Class) %>%
    dplyr::summarise(Affinity=median(Affinity), Rank=median(Rank)) %>%
    dplyr::mutate(HLA_Type="DQ")
  stringr::str_sub(res.DQA$"HLA", 5, 4) <- "_"
  res.DQB <- res[stringr::str_detect(HLA, "DQ"),] %>%
    tidyr::separate(col="HLA", into=c("HLA_A", "HLA"), sep="-") %>%
    dplyr::group_by(Peptide, HLA, Disease, HLA_Class) %>%
    dplyr::summarise(Affinity=median(Affinity), Rank=median(Rank)) %>%
    dplyr::mutate(HLA_Type="DQ")
  stringr::str_sub(res.DQB$"HLA", 5, 4) <- "_"
  return(data.table::rbindlist(list(res.DR, res.DPA, res.DPB, res.DQA, res.DQB)))
}

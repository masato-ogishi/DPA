#' Differential presentation analysis.
#'
#' NetMHC prediction results are parsed and the mean, maximum, and minimum rank percentiles among disease-predisposing and -protective HLA molecules are computed. Note: currently, only four-digit HLA alleles are used to determine HLA-disease association, and thus prediction results from two-digit HLAs are discarded.
#'
#' @param fileNames The set of NetMHC result files to be parsed. See the document of \code{NetMHC}.
#' @param hlaClass A string, either "ClassI" or "ClassII".
#' @param rankPerThr A threshold for rank percentile. For HLA class I and II, 2 and 10 are recommended, respectively.
#' @param DPIThr A threshold of differential presentation index, or DPI. We recommend 0.5. In this case, peptides with DPIs higher than 0.5 and lower than -0.5 are considered "predisposing" and "protective" epitopes, respectively.
#' @export
#' @rdname DifferentialPresentationAnalysis
#' @name DifferentialPresentationAnalysis
differentialPresentationAnalysis <- function(fileNames, hlaClass, rankPerThr=2, DPIThr=0.5){
  cat("Parsing NetMHC results...\n")
  dt_disease_hla <- data.table::copy(SummaryDF_Disease_HLA)[HLA_Digit==4,][,HLA:=gsub("HLA_", "", HLA)]
  if(hlaClass=="ClassI"){
    dt <- data.table::rbindlist(pbapply::pblapply(fileNames, NetMHC_Import_HLAI))
    dt <- merge(dt, dt_disease_hla[HLA_Class=="ClassI",])
  }else if(hlaClass=="ClassII"){
    dt <- data.table::rbindlist(pbapply::pblapply(fileNames, NetMHC_Import_HLAII))
    dt <- merge(dt, dt_disease_hla[HLA_Class=="ClassII",])
  }else{
    return("hlaClass must be either \"ClassI\" or \"ClassII\"!")
  }

  cat("Computing differential presentation index and classifying peptides...\n")
  rankPerThr <- -log10(rankPerThr/100)  ## Converting to the inverted log-scale
  dt <- dt %>%
    dplyr::group_by(Peptide, Disease, Association) %>%
    dplyr::summarise(
      Rank_Ave=mean(-log10(Rank/100)),
      Rank_Lo=min(-log10(Rank/100)), ## Weakest binding among HLAs tested
      Rank_Up=max(-log10(Rank/100))  ## Strongest binding among HLAs tested
    ) %>%
    dplyr::ungroup() %>%
    tidyr::gather(Stat, Value, -Peptide, -Disease, -Association) %>%
    dplyr::transmute(Peptide, Disease, Stat=paste0(Stat, "_", Association), Value) %>%
    tidyr::spread(Stat, Value) %>%
    data.table::as.data.table()
  dt[, DPI:=Rank_Up_Predisposing-Rank_Up_Protective]
  dt[, Category:="Others"]
  dt[DPI>DPIThr, Category:="Predisposing"]
  dt[-DPI>DPIThr, Category:="Protective"]
  dt[Rank_Up_Predisposing<rankPerThr&Rank_Up_Protective<rankPerThr, Category:="Others"]
  dt[,Category:=factor(Category, levels=c("Predisposing","Protective","Others"))]
  return(dt)
}

#' Sequence alignment-based homology scoring to UniProt proteomes.
#'
#' @param fastaFileName A UniProt proteome fasta file name to be parsed.
#' @param peptideSet A set of peptides.
#' @param proteome A set of protein sequences parsed by \code{UniProt_Proteome_Import}.
#' @param coreN The number of cores to be used for parallelization.
#' @export
#' @rdname UniProtAlignment
#' @name UniProtAlignment
UniProt_Proteome_Import <- function(fastaFileName){
  prot <- seqinr::read.fasta(fastaFileName, seqtype="AA", as.string=T)
  names(prot) <- data.table::transpose(data.table::as.data.table(strsplit(names(prot), "|", fixed=T)))[[2]]
  prot <- unlist(prot)
  sequenceFilter <- function(sequenceSet){
    s <- sequenceSet[!is.na(sequenceSet)]
    s <- toupper(s)
    letters <- unique(unlist(stringr::str_split(s, "")))
    letters.exclude <- setdiff(letters, Biostrings::AA_STANDARD)
    for (l in letters.exclude) {
      s <- s[!stringr::str_detect(s, stringr::fixed(l))]
    }
    return(s)
  }
  prot <- sequenceFilter(prot)
}

#' @export
#' @rdname UniProtAlignment
#' @name UniProtAlignment
UniProt_Proteome_ExactMatch <- function(peptideSet, proteome,
                                        coreN=parallel::detectCores(logical=F)){
  cl <- parallel::makeCluster(coreN, type="PSOCK")
  IDs <- names(proteome)
  seqs <- unlist(proteome)
  names(seqs) <- NULL
  res <- pbapply::pbsapply(
    peptideSet,
    function(p){
      paste0(IDs[which(stringr::str_detect(seqs, p))], collapse="|")
    },
    cl=cl
  )
  parallel::stopCluster(cl)
  gc();gc()
  return(res)
}

#' @export
#' @rdname UniProtAlignment
#' @name UniProtAlignment
UniProt_Proteome_AlignmentScores <- function(peptideSet, proteome,
                                             coreN=parallel::detectCores(logical=F)){
  cl <- parallel::makeCluster(coreN, type="PSOCK")
  proteome_AA <- Biostrings::AAStringSet(proteome)
  scoreMat <- pbapply::pbsapply(
    peptideSet,
    function(p){
      al <- Biostrings::pairwiseAlignment(
        pattern=proteome_AA,
        subject=p,
        type="local",
        substitutionMatrix="PAM30", gapOpening=9, gapExtension=1, ## Same as the parameters used in the blastp-short program
        scoreOnly=T
      )
      return(c(max(al),mean(al),min(al)))
    },
    cl=cl
  )
  parallel::stopCluster(cl)
  gc();gc()
  return(list("Highest"=scoreMat[1,], "Average"=scoreMat[2,], "Lowest"=scoreMat[3,]))
}


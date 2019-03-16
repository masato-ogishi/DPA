#' Selfness score.
#'
#' @param peptideSet A set of peptide sequences.
#' @param featureDT A peptide feature data.table.
#' @param sampling A sampling method passed to caret::trainControl.
#' @param seedSet A set of random seeds for ERT computation.
#' @param mtry A parameter passed to extraTrees::extraTrees.
#' @param numRandomCuts A parameter passed to extraTrees::extraTrees.
#' @param modelList A list of trained ERT models.
#' @export
#' @rdname Selfness
#' @name Selfness
Selfness_Features_HLAI <- function(peptideSet){
  peptidesToDescriptors <- function(peptideSet){
    dt.1 <- data.table::transpose(data.table::as.data.table(sapply(peptideSet, Peptides::crucianiProperties)))
    dt.2 <- data.table::transpose(data.table::as.data.table(sapply(peptideSet, Peptides::kideraFactors)))
    dt.3 <- data.table::transpose(data.table::as.data.table(sapply(peptideSet, Peptides::fasgaiVectors)))
    dt.4 <- data.table::transpose(data.table::as.data.table(sapply(peptideSet, Peptides::protFP)))
    dt.5 <- data.table::transpose(data.table::as.data.table(sapply(peptideSet, Peptides::blosumIndices)))
    colnames(dt.1) <- paste0("PP", 1:3)
    colnames(dt.2) <- paste0("KF", 1:10)
    colnames(dt.3) <- paste0("F", 1:6)
    colnames(dt.4) <- paste0("ProtFP", 1:8)
    colnames(dt.5) <- paste0("BLOSUM", 1:10)
    dt <- cbind(dt.1, dt.2, dt.3, dt.4, dt.5)
    dt[,Peptide:=peptideSet]
    return(dt)
  }
  peptidesToDummyFeatures <- function(peptideSet, maxLength=14){
    dt <- data.table::transpose(data.table::as.data.table(strsplit(peptideSet, "|")))
    for(i in 1:ncol(dt)) dt[[i]] <- factor(dt[[i]], levels=Biostrings::AA_STANDARD)
    colnames(dt) <- paste0("P", formatC(1:ncol(dt), flag="0", width=2))
    df <- mlr::createDummyFeatures(as.data.frame(dt))
    cols <- apply(data.table::CJ(paste0("P", formatC(1:maxLength, flag="0", width=2)), Biostrings::AA_STANDARD), 1, function(v){paste0(v[1], ".", v[2])})
    cols <- setdiff(cols, colnames(df))
    if(length(cols)>=1){
      for(i in cols){
        df[[i]] <- 0
      }
    }
    df$"Peptide" <- peptideSet
    return(df)
  }
  replaceNAToZero <- function(x){ for(j in names(x)) set(x,which(is.na(x[[j]])),j,0) }
  dt_feat <- peptidesToDescriptors(peptideSet)
  dt_feat[, Peptide_Length:=nchar(Peptide)]
  dt_feat_dummy <- lapply(split(dt_feat, by="Peptide_Length", keep.by=T, sorted=T), function(dt){
    peptidesToDummyFeatures(dt$"Peptide")
  }) %>% data.table::rbindlist(fill=T)
  replaceNAToZero(dt_feat_dummy)
  data.table::setcolorder(dt_feat_dummy, sort(colnames(dt_feat_dummy)))
  dt_feat <- merge(dt_feat, dt_feat_dummy, by="Peptide", sort=T)
  data.table::setcolorder(dt_feat, c("Peptide","Peptide_Length", setdiff(colnames(dt_feat), c("Peptide","Peptide_Length"))))
  return(dt_feat)
}

#' @export
#' @rdname Selfness
#' @name Selfness
Selfness_Features_HLAII <- function(peptideSet){
  peptidesToDescriptors <- function(peptideSet){
    dt.1 <- data.table::transpose(data.table::as.data.table(sapply(peptideSet, Peptides::crucianiProperties)))
    dt.2 <- data.table::transpose(data.table::as.data.table(sapply(peptideSet, Peptides::kideraFactors)))
    dt.3 <- data.table::transpose(data.table::as.data.table(sapply(peptideSet, Peptides::fasgaiVectors)))
    dt.4 <- data.table::transpose(data.table::as.data.table(sapply(peptideSet, Peptides::protFP)))
    dt.5 <- data.table::transpose(data.table::as.data.table(sapply(peptideSet, Peptides::blosumIndices)))
    colnames(dt.1) <- paste0("PP", 1:3)
    colnames(dt.2) <- paste0("KF", 1:10)
    colnames(dt.3) <- paste0("F", 1:6)
    colnames(dt.4) <- paste0("ProtFP", 1:8)
    colnames(dt.5) <- paste0("BLOSUM", 1:10)
    dt <- cbind(dt.1, dt.2, dt.3, dt.4, dt.5)
    dt[,Peptide:=peptideSet]
    return(dt)
  }
  dt_feat <- peptidesToDescriptors(peptideSet)
  dt_feat[, Peptide_Length:=nchar(Peptide)]
  data.table::setcolorder(dt_feat, c("Peptide","Peptide_Length", setdiff(colnames(dt_feat), c("Peptide","Peptide_Length"))))
  return(dt_feat)
}

#' @export
#' @rdname Selfness
#' @name Selfness
Selfness_Score <- function(featureDT, sampling=NULL, seedSet=1:5, mtry=16, numRandomCuts=16){
  # Random splitting
  trainIDList <- lapply(seedSet, function(s){set.seed(s);lapply(BBmisc::chunk(1:nrow(featureDT), n.chunks=5, shuffle=T), function(v){setdiff(1:nrow(featureDT), v)})})

  # The main workflow
  main <- function(df, trainIDs){
    df_train <- df[trainIDs, ]
    df_test <- df[-trainIDs, ]
    trgt <- df_train$"Origin"
    tab <- as.numeric(table(trgt))
    w <- 1/tab[trgt]
    modelCtrl <- caret::trainControl(method="cv", number=10, classProbs=T, sampling=sampling)
    model <- caret::train(Origin~., data=dplyr::select(df_train, -Peptide),
                          trControl=modelCtrl, metric="Kappa", weights=w,
                          method="extraTrees", numThreads=4,
                          tuneGrid=data.frame("mtry"=mtry, "numRandomCuts"=numRandomCuts))
    predDT <- data.table::data.table(
      "Peptide"=df_test$"Peptide",
      "Selfness"=predict(model, df_test, type="prob")[["Self"]]
    )
    return(list(
      "ertModel"=model, "predDT"=predDT
    ))
  }

  # Compute peptide "Selfness"
  resList <- foreach::foreach(i=1:length(seedSet))%do%{
    cat("Random seed = ", seedSet[i], "\n", sep="")
    set.seed(seedSet[i])
    pbapply::pblapply(trainIDList[[i]], function(trainIDs){main(featureDT, trainIDs)})
  }
  ertModelList <- resList %>%
    lapply(function(res){lapply(res, function(r){r$"ertModel"})}) %>%
    purrr::flatten()
  dt_self <- resList %>%
    lapply(function(res){lapply(res, function(r){r$"predDT"})}) %>%
    purrr::flatten() %>%
    data.table::rbindlist()
  dt_self <- dt_self[, .(Selfness=mean(Selfness), Selfness.sd=sd(Selfness)), by=Peptide]
  dt_self[, Selfness.cv:=Selfness.sd/Selfness]
  dt_self <- dt_self[, c("Peptide", "Selfness", "Selfness.cv"), with=F]
  return(list("ERTModels"=ertModelList, "SelfnessDT"=dt_self))
}

#' @export
#' @rdname Selfness
#' @name Selfness
Selfness_Score_Extrapolation <- function(featureDT, modelList, seedSet=1:5){
  seedSet <- rep(seedSet, each=length(modelList)/length(seedSet))
  dt_self <- data.table::rbindlist(pbapply::pblapply(1:length(modelList), function(i){
    set.seed(seedSet[i])
    data.table::data.table("Peptide"=featureDT$"Peptide", "Selfness"=predict(modelList[[i]], featureDT, type="prob")[["Self"]])
  }))
  dt_self <- dt_self[, .(Selfness=mean(Selfness), Selfness.sd=sd(Selfness)), by=Peptide]
  dt_self[, Selfness.cv:=Selfness.sd/Selfness]
  dt_self <- dt_self[, c("Peptide", "Selfness", "Selfness.cv"), with=F]
  return(dt_self)
}

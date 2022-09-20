library(parallel)
source('utilities.R')
source('data_utils.R')
source('constants.R')
source('prediction.R')

run.intend.using.tcga.lvo.preiction.model.customized <- function(exp.met.data, lvo.subtype, custom.subtype, verbose)
{
  expression.methylation <- exp.met.data
  regression.output <- load.expression.methylation.regression.output(output.file.path = sprintf(
    REGRESSION.OUTPUT.TRAIN.SUBTYPES.PATH.CV.LVO.FORMAT, lvo.subtype, UPSTREAM.MARGIN, DOWNSTREAM.MARGIN))
  
  predicted.expression <- load.precalculated.predicted.expression(
    sprintf(CV.LVO.SUBTYPE.NAME.FORMAT, lvo.subtype))
  
  stopifnot(all(colnames(expression.methylation$methylation) %in% colnames(predicted.expression)))
  predicted.expression <- predicted.expression[, colnames(expression.methylation$methylation)]
  
  projected.expression.and.methylation <- get.projected.expression.methylation.based.on.predicted.expression(
    expression.test = expression.methylation$expression,
    predicted.expression.test = predicted.expression,
    regression.output =  regression.output)
  
  projected.data.similarity <- projected.expression.and.methylation$W
  
  if (verbose) {
    utils.log(sprintf(
      "The weights matrix score for subtype %s is: %#.2f", 
      custom.subtype,
      get.cross.manifold.similarity.matrix.score(projected.data.similarity) * 100))
  }
  
  projected.expression <- projected.expression.and.methylation$X
  projected.methylation <- projected.expression.and.methylation$Y
  
  results.dir <- file.path(INTEND.RESULTS.DIR, custom.subtype)
  dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)
  
  write.table(projected.expression,
              file = file.path(results.dir, expression.projected.file.name()), sep = '\t')
  write.table(projected.methylation,
              file = file.path(results.dir, methylation.projected.file.name()), sep = '\t')
}

run.intend.using.tcga.leave.one.subtype.out.preiction.model <- function(subtype, verbose = FALSE)
{
  utils.log(sprintf("Running INTEND for subtype: %s", subtype))
  
  expression.methylation <- get.expression.methylation.data.for.expression.prediction(subtype)
  run.intend.using.tcga.lvo.preiction.model.customized(
    exp.met.data = expression.methylation,
    lvo.subtype = subtype,
    custom.subtype = subtype,
    verbose = verbose
  )
}

run.intend.using.preiction.model.for.multi.subtypes <- function(verbose = FALSE)
{
  train.types <- setdiff(TCGA.SUBTYPES, SUBTYPES.FOR.MULTI.SUBTYPE.ANALYSIS)
  test.types <- SUBTYPES.FOR.MULTI.SUBTYPE.ANALYSIS
  
  test.multi.type <- get.multi.subtypes.name(test.types)
  train.multi.type <- get.multi.subtypes.name(train.types)
  utils.log(sprintf("Running INTEND for subtype: %s", test.multi.type))
  
  expression.methylation <- get.expression.methylation.data.for.expression.prediction(test.multi.type)
  
  regression.output <- load.expression.methylation.regression.output(output.file.path = sprintf(
    REGRESSION.OUTPUT.TRAIN.SUBTYPES.PATH.FORMAT, UPSTREAM.MARGIN, DOWNSTREAM.MARGIN, train.multi.type))
  
  predicted.expression <- do.call(
    cbind,
    mclapply(test.types, function(test.type) load.precalculated.predicted.expression(
      sprintf(SUBTYPE.BASED.ON.TRAIN.SUBTYPES.FORMAT, test.type, train.multi.type)))
  )
  
  projected.expression.and.methylation <- get.projected.expression.methylation.based.on.predicted.expression(
    expression.test = expression.methylation$expression,
    predicted.expression.test = predicted.expression,
    regression.output =  regression.output)
  
  projected.data.similarity <- projected.expression.and.methylation$W
  
  if (verbose) {
    utils.log(sprintf(
      "The weights matrix score for subtype %s is: %#.2f", 
      test.multi.type,
      get.cross.manifold.similarity.matrix.score(projected.data.similarity) * 100))
  }
  
  projected.expression <- projected.expression.and.methylation$X
  projected.methylation <- projected.expression.and.methylation$Y
  
  results.dir <- file.path(INTEND.RESULTS.DIR, test.multi.type)
  dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)
  
  write.table(projected.expression,
              file = file.path(results.dir, expression.projected.file.name()), sep = '\t')
  write.table(projected.methylation,
              file = file.path(results.dir, methylation.projected.file.name()), sep = '\t')
}

# Integrate LUAD datasets: DNA methylation from TCGA and gene expression from LCCS
run.intend.for.luad.tcga.lccs.integration <- function()
{
  luad.lccs.expression <- get.lccs.luad.expression.preprocessed.and.normalized()
  
  luad.tcga.regression.output <- load.expression.methylation.regression.output(output.file.path = sprintf(
    REGRESSION.OUTPUT.TRAIN.SUBTYPES.PATH.CV.LVO.FORMAT, "LUAD", UPSTREAM.MARGIN, DOWNSTREAM.MARGIN))
  
  # Predicted from methylation data, using the leave LUAD subtype out model with 10 other cancer subtypes 
  # from TCGA as the regression train set
  luad.tcga.predicted.expression <- load.precalculated.predicted.expression(
    sprintf(CV.LVO.SUBTYPE.NAME.FORMAT, "LUAD"))
  
  projected.expression.and.methylation <- get.projected.expression.methylation.based.on.predicted.expression(
    expression.test = luad.lccs.expression,
    predicted.expression.test = luad.tcga.predicted.expression,
    regression.output =  luad.tcga.regression.output)
  
  projected.expression <- projected.expression.and.methylation$X
  projected.methylation <- projected.expression.and.methylation$Y
  
  results.dir <- file.path(INTEND.RESULTS.DIR, LUAD.TCGA.LCCS.INTEGRATION.SUBTYPE)
  dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)
  
  write.table(projected.expression,
              file = file.path(results.dir, expression.projected.file.name()), sep = '\t')
  write.table(projected.methylation,
              file = file.path(results.dir, methylation.projected.file.name()), sep = '\t')
}

run.cca <- function(data1, data2, required.dimension, standardise = TRUE)
{
  reduced.results <- run.cca.all.dimensions(data1, data2, standardise)
  
  return(list(
    reduced.1 = reduced.results$reduced.1[1:required.dimension, ],
    reduced.2 = reduced.results$reduced.2[1:required.dimension, ]
  ))
}

run.cca.all.dimensions <- function(data1, data2, standardise = TRUE)
{
  if (standardise) {
    data1 <- scale(data1)
    data2 <- scale(data2)
  }
  
  cross.product <- crossprod(data1, data2)
  cca.svd <- svd(cross.product)
  
  reduced.data.1 <- t(cca.svd$u)
  reduced.data.2 <- t(cca.svd$v)
  
  rownames(reduced.data.1) <- paste0("CC", 1:nrow(reduced.data.1))
  rownames(reduced.data.2) <- paste0("CC", 1:nrow(reduced.data.2))
  
  colnames(reduced.data.1) <- colnames(data1)
  colnames(reduced.data.2) <- colnames(data2)
  
  return(list(reduced.1 = reduced.data.1, reduced.2 = reduced.data.2))
}

get.intend.projections <- function(subtype, dimension, results.dir = INTEND.RESULTS.DIR)
{
  results.dir <- file.path(results.dir, subtype)
  
  projected.expression <- as.matrix(read.table(file = file.path(results.dir, expression.projected.file.name())))
  projected.methylation <- as.matrix(read.table(file = file.path(results.dir, methylation.projected.file.name())))
  
  projected.results.reduced <- run.cca(projected.expression,
                                       projected.methylation,
                                       required.dimension = dimension)
  
  return(list(
    expression = projected.results.reduced$reduced.1,
    methylation =  projected.results.reduced$reduced.2))
}

create.and.run.intend.for.skcm.split.data <- function() {
  
  skcm.exp.met.orig <- get.expression.methylation.data.without.reverse.methylation("SKCM")
  
  seeds <- c(1506848755, -408267596, 204949605, -2099428689, -143742477,
             1569065296, 39514193, -10279358, -1384063582, 618835332,
             -1631594997, -903834144, -1922619212, -1735708480, -864213837,
             617094112, -1171008603, -118026804, 26347911, -1269484645,
             690673515, 500715500, 241174634, 139965732, 828805704,
             -600683108, 152995369, -618807196, 987457920, 638808615)
  
  for (i in 1:length(seeds)) {
    set.seed(seeds[i])
    custom.subtype <- sprintf("SKCM_SPLIT_%d", i)
    print(custom.subtype)
    
    samples <- colnames(skcm.exp.met.orig$gene.expression.data)
    ge.samples <- sample(samples, size = as.integer(length(samples) * 0.5))
    dm.samples <- setdiff(samples, ge.samples)
    
    skcm.exp.met <- list()
    skcm.exp.met$gene.expression.data <- skcm.exp.met.orig$gene.expression.data[, ge.samples]
    skcm.exp.met$gene.expression.data.aggregated <- skcm.exp.met.orig$gene.expression.data.aggregated[, ge.samples]
    skcm.exp.met$methy.data <- skcm.exp.met.orig$methy.data[, dm.samples]
    skcm.exp.met$gene.methylation.data <- skcm.exp.met.orig$gene.methylation.data[, dm.samples]
    
    expression.methylation.data.file.path <- sprintf(PREPROCESSED.EXPRESSION.METHYLATION.DATA.FORMAT,
                                                     "without_reverse_methylation",
                                                     custom.subtype)
    
    saveRDS(skcm.exp.met, file = expression.methylation.data.file.path)
    
    skcm.exp.met <- get.expression.methylation.data.for.expression.prediction(custom.subtype)
    
    
    run.intend.using.tcga.lvo.preiction.model.customized(
      exp.met.data = skcm.exp.met, 
      lvo.subtype = "SKCM",
      custom.subtype = custom.subtype,
      verbose = FALSE # important as this will fail the FOSCTTM tests 
    )
    
    skcm.umap <- run.umap(custom.subtype, skcm.exp.met)
    
  }
}

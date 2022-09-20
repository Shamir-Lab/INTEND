library(glmnet)
library(dplyr)
library(parallel)
source('utilities.R')
source('gene_mapping_utils.R')
source('constants.R')
source('data_utils.R')

get.lasso.regression <- function(X.train, Y.train, n.fold = 10)
{
  # X.train is an n*m matrix containing n samples of the m independent variables
  # Y.train is an n-length vector containing the dependent variable values in each sample in the training set
  # n.fold is regarding choosing the best model parameters on the training data
  return(cv.glmnet(x = X.train, y = Y.train, nfold = n.fold))
}

load.precalculated.predicted.expression <- function(new.subtype)
{
  PREDICTED.EXPRESSION.FOR.NEW.SUBTYPE.PATH.FORMAT <-
    "CachedData/expression_methylation_regression/predicted_expression_%s.rda"
  
  predicted.expression.file.path <- 
    sprintf(PREDICTED.EXPRESSION.FOR.NEW.SUBTYPE.PATH.FORMAT, new.subtype)
  
  stopifnot(file.exists(predicted.expression.file.path))
  return(readRDS(file = predicted.expression.file.path))
}

get.predicted.expression <- function(
  methylation,
  gene.to.cpg.site.mapping,
  regression.output,
  new.subtype,
  force.update = FALSE)
{
  PREDICTED.EXPRESSION.FOR.NEW.SUBTYPE.PATH.FORMAT <-
    "CachedData/expression_methylation_regression/predicted_expression_%s.rda"
  
  predicted.expression.file.path <- 
    sprintf(PREDICTED.EXPRESSION.FOR.NEW.SUBTYPE.PATH.FORMAT, new.subtype)
  
  if (file.exists(predicted.expression.file.path) & !force.update) {
    return(readRDS(file = predicted.expression.file.path))
  }
  
  gc()
  predicted.expression <- do.call(rbind, mclapply(1:nrow(regression.output), function(i) {
    
    if (i %% 500 == 1) {
      utils.log(sprintf("Done calculating predicted expression for gene number %d", i - 1))
    }
    
    gene <- rownames(regression.output)[i]
    cg.sites <- gene.to.cpg.site.mapping[gene, "cg_sites"][[1]]
    
    X <- t(methylation[cg.sites, , drop = FALSE])
    fit <- regression.output$fit.parameters[[i]]
    
    return(as.vector(predict(fit, newx = X, s = "lambda.min")))
  },
  mc.cores = LIMITED.NUMBER.OF.CORES)) # Set to low number to avoid memory exhaustion
  
  rownames(predicted.expression) <- rownames(regression.output)
  colnames(predicted.expression) <- colnames(methylation)
  
  saveRDS(predicted.expression, file = predicted.expression.file.path)
  return(predicted.expression)
}

load.expression.methylation.regression.output <- function(output.file.path)
{
  stopifnot(file.exists(output.file.path))
  return(readRDS(file = output.file.path))
}

get.expression.methylation.regression.output <- 
  function(expression.methylation.train.test, gene.to.cpg.site.mapping, output.file.path, force.update = FALSE)
{
    if (file.exists(output.file.path) & !force.update) {
      return(readRDS(file = output.file.path))
    }
    
    expression.train <- expression.methylation.train.test$expression.train
    methylation.train <- expression.methylation.train.test$methylation.train
    
    regression.output <- data.frame(row.names = gene.to.cpg.site.mapping$gene)
    regression.output$fit.parameters <- list(NA)
    regression.output$predicted.expression.train <- list(NA)
    regression.output$R.squared.train <- rep(NA)
    gc()
    output.for.test.genes <- mclapply(1:nrow(regression.output), function(i) {
      
      test.gene <- gene.to.cpg.site.mapping$gene[i]
      test.cg.sites <- gene.to.cpg.site.mapping$cg_sites[[i]]
      
      
      output.for.test.gene <- tryCatch({
        
        X.train <- t(methylation.train[test.cg.sites, , drop = FALSE])
        Y.train <- expression.train[test.gene,]
        
        fit <- get.lasso.regression(X.train, Y.train)
        
        predicted.expression.train <- predict(fit, newx = X.train, s = "lambda.min")
        residuals.train <- Y.train - predicted.expression.train
        R.squared.train <- (var(as.vector(Y.train)) - var(as.vector(residuals.train))) / var(as.vector(Y.train))
        list(fit.parameters = fit,
             predicted.expression.train = predicted.expression.train,
             R.squared.train = R.squared.train
        )
      }, error = function(err) {
        list(fit.parameters = NA,
             predicted.expression.train = NA,
             R.squared.train = NA
        )
        
      }, finally = function() {})
      
      return(output.for.test.gene)
    },
    mc.cores = LIMITED.NUMBER.OF.CORES) # Set to low number to avoid memory exhaustion
    
    for (i in 1:nrow(regression.output))  {
      
      regression.output$fit.parameters[[i]] <- output.for.test.genes[[i]]$fit.parameters
      regression.output$predicted.expression.train[[i]] <- output.for.test.genes[[i]]$predicted.expression.train
      regression.output$R.squared.train[i] <- output.for.test.genes[[i]]$R.squared.train
    }
    
    # Remove NA values
    regression.output <- regression.output %>% 
      filter(!is.na(fit.parameters) & !is.na(predicted.expression.train) & !is.na(R.squared.train))
    
    saveRDS(regression.output, file = output.file.path)
    return(regression.output)
}

get.expression.methylation.for.train.subtypes <- function(
  train.subtypes,
  expression.features.all.subtypes,
  methylation.features.all.subtypes,
  force.update = FALSE)
{
  subtypes.str <- get.multi.subtypes.name(train.subtypes)
  cached.data.for.train.subtypes.file.path <-
    sprintf("CachedData/expression_methylation_regression/expression_methylation_train_%s.rda", subtypes.str)
  
  if (file.exists(cached.data.for.train.subtypes.file.path) & !force.update) {
    return(readRDS(file = cached.data.for.train.subtypes.file.path))
  }
  
  expression.all.train <- matrix(
    nrow = length(expression.features.all.subtypes),
    ncol = 0,
    dimnames = list(expression.features.all.subtypes))
  
  methylation.all.train <- matrix(
    nrow = length(methylation.features.all.subtypes),
    ncol = 0,
    dimnames = list(methylation.features.all.subtypes))
  
  samples.indices.by.subtype <- list()

  
  for (subtype in train.subtypes) {
    utils.log(sprintf("Adding %s", subtype))
    expression.methylation <- get.expression.methylation.data.for.expression.prediction(subtype)
    
    subtype.samples.start.index <- ncol(expression.all.train) + 1
    subtype.samples.end.index <- subtype.samples.start.index + ncol(expression.methylation$expression) - 1
    
    samples.indices.by.subtype[[subtype]] <- subtype.samples.start.index:subtype.samples.end.index
    
    expression.all.train <- cbind(expression.all.train,
                                  expression.methylation$expression[expression.features.all.subtypes, ])
    methylation.all.train <- cbind(methylation.all.train,
                                   expression.methylation$methylation[methylation.features.all.subtypes, ])
    
  }
  
  # We set test data to NA, and save the regression parameters learned on the train subtypes data
  # for a later use of predicting the expression from the expression in each of the test subtypes data.
  expression.methylation.for.train.subtypes <- list(
    expression.train = expression.all.train,
    methylation.train = methylation.all.train,
    samples.indices.by.subtype = samples.indices.by.subtype)
  
  saveRDS(expression.methylation.for.train.subtypes, file = cached.data.for.train.subtypes.file.path)
  return(expression.methylation.for.train.subtypes)
}

get.projected.expression.methylation.based.on.predicted.expression <- function(
  expression.test, 
  predicted.expression.test,
  regression.output)
{
  NUMBER.OF.GENES.FOR.SIMILARITY.COMPUTATION <- 2000
  
  # We're not using here R.squared.test as it's based on Y.test 
  genes.with.best.R.squared <- rownames(regression.output)[order(
    regression.output$R.squared.train, decreasing = TRUE)[1:NUMBER.OF.GENES.FOR.SIMILARITY.COMPUTATION]]
  
  
  expression.variation.by.gene <- apply(expression.test, 1, var)
  real.high.var.expressed.genes <- names(head(sort(expression.variation.by.gene, decreasing = TRUE),
                                              n = NUMBER.OF.GENES.FOR.SIMILARITY.COMPUTATION))
  
  predicted.expression.variavtion.by.gene <- apply(predicted.expression.test, 1, var)
  predicted.high.var.expressed.genes <- names(
    head(sort(predicted.expression.variavtion.by.gene, decreasing = TRUE),
         n = NUMBER.OF.GENES.FOR.SIMILARITY.COMPUTATION))
  
  high.var.expressed.genes <- intersect(real.high.var.expressed.genes, predicted.high.var.expressed.genes)
  
  genes.for.similarity.computation <- intersect(genes.with.best.R.squared, high.var.expressed.genes)
  
  utils.log(sprintf("Proceeding with %d genes for similarity computation",
                    length(genes.for.similarity.computation)))
  
  real.expression.for.selected.genes <- expression.test[genes.for.similarity.computation, ]
  predicted.expression.for.selected.genes <- predicted.expression.test[genes.for.similarity.computation, ]
  
  real.expression.for.selected.genes <- scale.data(real.expression.for.selected.genes)
  predicted.expression.for.selected.genes <- scale.data(predicted.expression.for.selected.genes)
  
  # Compute expression methylation similarity matrix. The values used for the methylation are the predicted
  # expression according to the methylation
  
  expression.methylation.similarity <- get.inter.similarity.matrix(
    real.expression.for.selected.genes,
    predicted.expression.for.selected.genes)
  
  data <- list(
    X = real.expression.for.selected.genes,
    Y = predicted.expression.for.selected.genes,
    W = expression.methylation.similarity)
  
  return(data)
}

train.expression.prediction.for.all.subtypes.using.constant.train.set <- function()
{
  train.subtypes <- setdiff(TCGA.SUBTYPES, SUBTYPES.FOR.MULTI.SUBTYPE.ANALYSIS)
  test.subtypes <- SUBTYPES.FOR.MULTI.SUBTYPE.ANALYSIS
  
  common.features.for.all.subtypes <- get.common.features.for.all.subtypes()
  
  mapping.all.subtypes <- get.filtered.gene.to.cg.sites.mapping.all.subtypes(
    common.features.for.all.subtypes$expression,
    common.features.for.all.subtypes$methylation,
    UPSTREAM.MARGIN,
    DOWNSTREAM.MARGIN)
  
  rownames(mapping.all.subtypes) = mapping.all.subtypes$gene

  utils.log(sprintf("Getting expression and methylation for train subtypes: %s",
                    paste(train.subtypes, collapse = ", ")))
  
  expression.methylation.for.train.subtypes <- get.expression.methylation.for.train.subtypes(
    train.subtypes = train.subtypes,
    expression.features.all.subtypes = common.features.for.all.subtypes$expression,
    methylation.features.all.subtypes = common.features.for.all.subtypes$methylation)
  
  train.subtypes.str <- get.multi.subtypes.name(train.subtypes)
  
  utils.log(sprintf("Learning regression model from train model..."))
  
  regression.output.from.train.subtypes <- get.expression.methylation.regression.output(
    expression.methylation.train.test = expression.methylation.for.train.subtypes,
    gene.to.cpg.site.mapping = mapping.all.subtypes,
    output.file.path = 
      sprintf(REGRESSION.OUTPUT.TRAIN.SUBTYPES.PATH.FORMAT,
              UPSTREAM.MARGIN,
              DOWNSTREAM.MARGIN,
              train.subtypes.str))
  
  utils.log(sprintf("Done learning"))
  
  for (test.subtype in test.subtypes) {
    
    utils.log(sprintf("Running for test subtype: %s", test.subtype))
    
    expression.methylation.for.test.subtype <-
      get.expression.methylation.data.for.expression.prediction(test.subtype)
    
    expression <- expression.methylation.for.test.subtype$expression[
      common.features.for.all.subtypes$expression, ]
    methylation <- expression.methylation.for.test.subtype$methylation[
      common.features.for.all.subtypes$methylation, ]
    
    new.subtype <- sprintf(SUBTYPE.BASED.ON.TRAIN.SUBTYPES.FORMAT, test.subtype, train.subtypes.str)
    
    predicted.expression <- get.predicted.expression(
      methylation = methylation,
      gene.to.cpg.site.mapping = mapping.all.subtypes,
      regression.output = regression.output.from.train.subtypes,
      new.subtype = new.subtype)
    
    projected.expression.methylation <- get.projected.expression.methylation.based.on.predicted.expression(
      expression.test = expression,
      predicted.expression.test = predicted.expression,
      regression.output = regression.output.from.train.subtypes)

    projected.data.similarity.matrix <- projected.expression.methylation$W
    utils.log(sprintf("The weights matrix score for test subtype %s is: %#.2f", 
                      test.subtype,
                      get.cross.manifold.similarity.matrix.score(projected.data.similarity.matrix) * 100))
    
    rm(expression.methylation.for.test.subtype,
       expression, 
       methylation,
       predicted.expression,
       projected.expression.methylation)
    gc()
  }
}

train.expression.prediction.for.all.subtypes.leaving.one.subtype.out <- function()
{
  common.features.for.all.subtypes <- get.common.features.for.all.subtypes()
  
  mapping.all.subtypes <- get.filtered.gene.to.cg.sites.mapping.all.subtypes(
    common.features.for.all.subtypes$expression,
    common.features.for.all.subtypes$methylation,
    UPSTREAM.MARGIN,
    DOWNSTREAM.MARGIN)
  
  rownames(mapping.all.subtypes) = mapping.all.subtypes$gene
  
  cross.validation.subtypes <- TCGA.SUBTYPES

  expression.methylation.for.subtypes <- get.expression.methylation.for.train.subtypes(
    train.subtypes = cross.validation.subtypes,
    expression.features.all.subtypes = common.features.for.all.subtypes$expression,
    methylation.features.all.subtypes = common.features.for.all.subtypes$methylation)
  
  for (test.subtype in cross.validation.subtypes) {
    
    utils.log(sprintf("Running for test subtype: %s", test.subtype))
    
    test.subtype.samples.indices <- 
      expression.methylation.for.subtypes$samples.indices.by.subtype[[test.subtype]]
    
    expression.methylation.for.train.subtypes <- list(
      expression.train = 
        expression.methylation.for.subtypes$expression.train[, -test.subtype.samples.indices],
      methylation.train = 
        expression.methylation.for.subtypes$methylation.train[, -test.subtype.samples.indices]
    )
  
    regression.output.from.train.subtypes <- get.expression.methylation.regression.output(
      expression.methylation.train.test = expression.methylation.for.train.subtypes,
      gene.to.cpg.site.mapping = mapping.all.subtypes,
      output.file.path = 
        sprintf(REGRESSION.OUTPUT.TRAIN.SUBTYPES.PATH.CV.LVO.FORMAT, 
                test.subtype, 
                UPSTREAM.MARGIN,
                DOWNSTREAM.MARGIN))
    
    expression <- expression.methylation.for.subtypes$expression.train[, test.subtype.samples.indices]
    methylation <- expression.methylation.for.subtypes$methylation.train[, test.subtype.samples.indices]
    
    new.subtype <- sprintf(CV.LVO.SUBTYPE.NAME.FORMAT, test.subtype)
    
    predicted.expression <- get.predicted.expression(
      methylation = methylation,
      gene.to.cpg.site.mapping = mapping.all.subtypes,
      regression.output = regression.output.from.train.subtypes,
      new.subtype = new.subtype)

    projected.expression.methylation <- get.projected.expression.methylation.based.on.predicted.expression(
      expression.test = expression,
      predicted.expression.test = predicted.expression,
      regression.output = regression.output.from.train.subtypes)
    
    projected.data.similarity.matrix <- projected.expression.methylation$W
    utils.log(sprintf("The weights matrix score for test subtype %s is: %#.2f", 
                      test.subtype,
                      get.cross.manifold.similarity.matrix.score(projected.data.similarity.matrix) * 100))
    
    rm(expression.methylation.for.train.subtypes,
       regression.output.from.train.subtypes,
       expression, 
       methylation,
       predicted.expression,
       projected.expression.methylation)
    gc()
  }
}

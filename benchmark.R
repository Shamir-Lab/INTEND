source('utilities.R')
source('data_utils.R')
source('constants.R')
source('projection_functions.R')
source('intend.R')
library(parallel)
library(rliger)
library(Seurat)
library(glue)

# Assuming correspondence between gene expression and gene methylation features
get.gene.expression.methylation.data <- function(subtype, number.of.feature.to.keep)
{
  expression.methylation.data <- get.expression.methylation.data.before.feature.selection.and.scaling(subtype)
  
  gene.expression.data.aggregated <- expression.methylation.data$gene.expression.data.aggregated
  gene.methylation.data <- expression.methylation.data$gene.methylation.data
  
  # In order to compute the expression-methylation similarity matrix, we select the highly variable genes
  # according to gene.expression.data.aggregated and intersect gene.methylation.data with those genes
  top.features <-
    get.highly.variable.features(data = gene.expression.data.aggregated, number.of.feature.to.keep)
  
  gene.expression.data.aggregated <- gene.expression.data.aggregated[top.features, ]
  gene.methylation.data <- gene.methylation.data[top.features, ]
  
  # Scale & center the datasets
  gene.methylation.data <- scale.data(gene.methylation.data)
  gene.expression.data.aggregated <- scale.data(gene.expression.data.aggregated)

  # Compute expression methylation similarity matrix
  expression.methylation.similarity <- get.inter.similarity.matrix(
    gene.expression.data.aggregated,
    gene.methylation.data)
  
  data <- list(
    X = gene.expression.data.aggregated,
    Y = gene.methylation.data,
    W = expression.methylation.similarity)
  
  return(data)
}

get.gene.expression.methylation.data.for.liger <- function(subtype)
{
  expression.methylation.data <- get.expression.methylation.data.before.feature.selection.and.scaling(subtype)
  
  gene.expression.data <- expression.methylation.data$gene.expression.data.aggregated
  gene.methylation.data <- expression.methylation.data$gene.methylation.data
  
  # rename methyl_data colnames as liger does not accept same 'cells' (colnames) in the 2 datasets
  colnames(gene.methylation.data) <- sapply(colnames(gene.methylation.data), function(x) paste(x, "changed", sep = "_"))
  
  data.for.liger <- list(gene.expression.data = gene.expression.data, gene.methylation.data = gene.methylation.data)
  return(data.for.liger)
}

run_manifold_alignment <- function(data, mu, k=0) # k is only for "documentation" in the return value
{
  X <- data$X
  Y <- data$Y
  W <- data$W
  
  d <- min(dim(X)[2], dim(Y)[2])
  
  projection_mappings <- generate_projection_mappings(X, Y, W, d, mu)
  X_projected <- t(projection_mappings$alpha) %*% X
  Y_projected <- t(projection_mappings$beta) %*% Y
  
  output <- list(
    X = X,
    Y = Y,
    W = W,
    mu = mu,
    d = d,
    k = k,
    projection_mappings = projection_mappings,
    X_projected = X_projected,
    Y_projected = Y_projected
  )
  return(output)
}

run.manifold.alignment.for.expression.methylation.integration.with.correspondence <- 
  function(results.dir, subtype, expression.methylation.data)
{
  results.dir <- file.path(results.dir, subtype)
  dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)
  
  write.table(expression.methylation.data$X, file = file.path(results.dir, expression.file.name()), sep = '\t')
  write.table(expression.methylation.data$Y, file = file.path(results.dir, methylation.file.name()), sep = '\t')
  write.table(expression.methylation.data$W,
              file = file.path(results.dir, 'W_cross_similarity.tsv'),
              sep = '\t')

  # Algorithm hyper-parameters
  mu <- 1
  
  output <- run_manifold_alignment(expression.methylation.data, mu)
  
  write.table(output$X_projected, file = file.path(results.dir, expression.projected.file.name()), sep = '\t')
  write.table(output$Y_projected, file = file.path(results.dir, methylation.projected.file.name()), sep = '\t')
  
  write.table(output$projection_mappings$alpha, file = file.path(results.dir, 'alpha.tsv'), sep = '\t')
  write.table(output$projection_mappings$beta, file = file.path(results.dir, 'beta.tsv'), sep = '\t')
}

get.manifold.alignment.projections <- function(subtype, dimension, results.dir)
{
  results.dir <- file.path(results.dir, subtype)
  
  X.projected <- as.matrix(read.table(file = file.path(results.dir, expression.projected.file.name())))
  Y.projected <- as.matrix(read.table(file = file.path(results.dir, methylation.projected.file.name())))
  
  return(list(expression = X.projected[1:dimension, ], methylation =  Y.projected[1:dimension, ]))
}

compute.projection.scores.for.manifold.alignment.projections <- 
  function(results.dir, subtype, d.range, verbose = FALSE)
{
  results.dir <- file.path(results.dir, subtype)
  
  X <- as.matrix(read.table(file = file.path(results.dir, expression.file.name())))
  Y <- as.matrix(read.table(file = file.path(results.dir, methylation.file.name())))
  
  X.projected <- as.matrix(read.table(file = file.path(results.dir, expression.projected.file.name())))
  Y.projected <- as.matrix(read.table(file = file.path(results.dir, methylation.projected.file.name())))
  
  cross.manifold.similarity.scores <-
    get_projection_cross_manifold_similarity_scores_for_d_range(X.projected, Y.projected, d.range)
  best.cross.manifold.similarity.score <- min(cross.manifold.similarity.scores)
  best.d <- which.min(cross.manifold.similarity.scores) + (d.range[1] - 1)
  
  x.similarity.preservation.score <-
    get_projection_manifold_similarity_preservation_score(X, X.projected[1:best.d, ])
  
  y.similarity.preservation.score <-
    get_projection_manifold_similarity_preservation_score(Y, Y.projected[1:best.d, ])
  
  if (verbose)
  {
    utils.log(sprintf("The best projection cross manifold similarity score is: %f (d=%d)", 
                  best.cross.manifold.similarity.score, 
                  best.d))
    
    
    utils.log(sprintf("The projection similarity preservation score for X and d=%d is: %f",
                  best.d,
                  x.similarity.preservation.score))
    
    
    utils.log(sprintf("The projection similarity preservation score for Y and d=%d is: %f",
                  best.d,
                  y.similarity.preservation.score))
  }
  
  return(cross.manifold.similarity.scores)
}

run.liger.for.expression.methylation.integration <- function(
  results.dir,
  subtype,
  dimensions,
  expression.methylation.data.for.liger
  )
{
  gene.expression.data <- expression.methylation.data.for.liger$gene.expression.data
  gene.methylation.data <- expression.methylation.data.for.liger$gene.methylation.data
  
  rna.met <- createLiger(list(rna = gene.expression.data, met = gene.methylation.data))
  rna.met <- normalize(rna.met)
  rna.met <- selectGenes(rna.met, datasets.use = c(1))
  rna.met <- scaleNotCenter(rna.met)
  
  # The methylation is alreday with reversed sign, and also the proportional nature of the gene body 
  # methylation makes it unnecessary to normalize and scale the methylation data. Thus, we simply overwrite
  # the scaled methylation data with the reversed methylation data.
  rna.met@scale.data$met <- t(as.matrix(rna.met@raw.data$met[rna.met@var.genes, ]))
  
  X <- t(rna.met@scale.data$rna)
  Y <- t(rna.met@scale.data$met)
  
  results.dir <- file.path(results.dir, subtype)
  dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)
  
  write.table(X, file = file.path(results.dir, expression.file.name()), sep = '\t')
  write.table(Y, file = file.path(results.dir, methylation.file.name()), sep = '\t')
  
  mclapply(dimensions, function(liger.dimension)
  {
    utils.log(sprintf("Running liger optimization for d=%d...", liger.dimension))
    
    rna.met <- tryCatch(optimizeALS(rna.met, k = liger.dimension), error = function(err) NULL)
    # Due to a possible bug in LIGER in which for some values of k we encounter an error:
    # "Error in if (trim > 0 && n) { : missing value where TRUE/FALSE needed"
    if (is.null(rna.met)) 
    {
      return()
    }
    rna.met <- quantileAlignSNF(rna.met, center = T)
    
    X.projected <- t(rna.met@H.norm[1:ncol(gene.expression.data),])
    Y.projected <- t(rna.met@H.norm[(ncol(gene.expression.data) + 1):
                                      (ncol(gene.expression.data) + ncol(gene.methylation.data)), ])
    
    write.table(
      X.projected,
      file = file.path(results.dir, expression.projected.file.name(liger.dimension)),
      sep = '\t')
    
    write.table(
      Y.projected, 
      file = file.path(results.dir, methylation.projected.file.name(liger.dimension)),
      sep = '\t')
    
  })
}

get.liger.projections <- function(subtype, dimension, results.dir = LIGER.RESULTS.DIR)
{
  results.dir <- file.path(results.dir, subtype)
  
  return(list(
    expression = as.matrix(read.table(file.path(results.dir, expression.projected.file.name(dimension)))),
    methylation =  as.matrix(read.table(file.path(results.dir, methylation.projected.file.name(dimension))))
  ))
}

compute.projection.scores.for.liger.projections <- function(
  subtype,
  dimensions,
  results.dir = LIGER.RESULTS.DIR,
  verbose = FALSE)
{
  cross.manifold.similarity.scores <- mcmapply(function(liger.dimension) {
    
    cross.manifold.similarity.score.for.dimension <- tryCatch({
      
      projections <- get.liger.projections(subtype, liger.dimension, results.dir )

      cross.manifold.similarity.score <- 
        get_projection_cross_manifold_similarity_score(projections$expression, projections$methylation)
      
      if (verbose) {
        utils.log(sprintf("Liger projection cross manifold similarity score is: %f (liger.dimension=%d)", 
                          cross.manifold.similarity.score, 
                          liger.dimension))
      }
      
      cross.manifold.similarity.score
      
    }, error = function(err) {
      
      NA
      
    }, finally = function() {})
    
    return(cross.manifold.similarity.score.for.dimension)
  },
  dimensions)
  
  names(cross.manifold.similarity.scores) <- dimensions

  return(cross.manifold.similarity.scores)
}

write.matrices.to.tsv.file.for.mmd.ma <- function(results.dir, subtype, expression.methylation.data)
{
  X <- expression.methylation.data$X
  Y <- expression.methylation.data$Y
  W_x <- get_similarity_matrix(X)
  diag(W_x) <- 0
  W_y <- get_similarity_matrix(Y)
  diag(W_y) <- 0
  
  write.table(W_x, file = file.path(results.dir, expression.similarities.file.name()), quote = FALSE,
              sep = '\t', col.names = FALSE, row.names = FALSE)
  write.table(W_y, file = file.path(results.dir, methylation.similarities.file.name()), quote = FALSE,
              sep = '\t', col.names = FALSE, row.names = FALSE)

  return(list(W_x = W_x, W_y = W_y))
}

run.mmd.ma.script <- function(
  results.dir,
  subtype, 
  dimensions,
  expression.methylation.data) {
  
  # Currently we ignore the dimensions parameter and it's hardcoded in the MMD-MA python script
  stopifnot(identical(dimensions, dimensions.for.mmd.ma))
  
  # We call here normalizePath as we pass these as an arguments to an external process
  # and the paths are different on Windows and unix-based systems
  results.dir <- normalizePath(file.path(results.dir, subtype))
  expression.matrix <- normalizePath(file.path(results.dir, expression.similarities.file.name()))
  methylation.matrix <- normalizePath(file.path(results.dir, methylation.similarities.file.name()))
  python.mmd.ma.script <- normalizePath(file.path('../MMD-MA',
                                             "manifoldAlignDistortionPen_mmd_multipleStarts.py"))
  
  python.exe <- ifelse(IS.UNIX, "python", "C:\\venv\\Scripts\\python.exe")
  
  dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)
  
  # First write the input matrices for the MMD-MA script
  write.matrices.to.tsv.file.for.mmd.ma(
    results.dir = results.dir,
    subtype = subtype,
    expression.methylation.data = expression.methylation.data)
  
  command <- glue("
                    {python.exe} {python.mmd.ma.script} {expression.matrix} {methylation.matrix} \\
                    --p {max(dimensions)} --results_dir {results.dir}
                    ")
    
  exit.code <- system(command, ignore.stderr = TRUE)
  
  if (exit.code != 0) {
    stop(glue("Failed running MMD-MA (exit code: {exit.code})"))
  }
    
}

get.mmd.ma.projections <- function(subtype, dimension, results.dir)
{
  results.dir <- file.path(results.dir, subtype)
  
  W_x <- read.table(file = file.path(results.dir, expression.similarities.file.name()))
  W_y <- read.table(file = file.path(results.dir, methylation.similarities.file.name()))
  
  X_mmd_projection_mapping <- read.table(file = file.path(results.dir, sprintf('alpha_final_%d_0.txt', dimension)))
  Y_mmd_projection_mapping <- read.table(file = file.path(results.dir, sprintf('beta_final_%d_0.txt', dimension)))
  
  X_mmd_projected <- t(as.matrix(X_mmd_projection_mapping)) %*% as.matrix(W_x)
  Y_mmd_projected <- t(as.matrix(Y_mmd_projection_mapping)) %*% as.matrix(W_y)
  
  return(list(expression = X_mmd_projected, methylation =  Y_mmd_projected))
}

compute.projection.scores.for.mmd.ma.projections <- 
  function(results.dir, subtype, dimensions, verbose = FALSE)
{
  cross.manifold.similarity.scores <- mcmapply(function(d) {
    
    projections <- get.mmd.ma.projections(subtype, d, results.dir)
    
    cross.manifold.similarity.score <- 
      get_projection_cross_manifold_similarity_score(projections$expression, projections$methylation)
    
    if (verbose) {
      utils.log(sprintf("MMD-MA projection cross manifold similarity score is: %f (dimension=%d)",
                        cross.manifold.similarity.score,
                        d))
    }
    
    return(cross.manifold.similarity.score)
  },
  dimensions)
  
  names(cross.manifold.similarity.scores) <- dimensions
  
  return(cross.manifold.similarity.scores)
}

run.seurat.for.expression.methylation.integration <-
  function(results.dir, subtype, dimensions, expression.methylation.data.for.seurat)
{
  rna <- CreateSeuratObject(counts = expression.methylation.data.for.seurat$gene.expression.data,
                            project = "RNA",
                            assay = "RNA")
  met <- CreateSeuratObject(counts = expression.methylation.data.for.seurat$gene.methylation.data,
                            project = "MET",
                            assay = "MET")
  
  # This line is not necessary according to Seurat's documentation,
  # but seems to significantly improve its results
  rna <- NormalizeData(rna, normalization.method = "RC", verbose = FALSE)

  rna <- FindVariableFeatures(object = rna, selection.method = "vst", verbose = FALSE)
  met <- FindVariableFeatures(object = met, selection.method = "vst", verbose = FALSE)
  
  rna <- ScaleData(rna)
  met <- ScaleData(met)
  
  results.dir <- file.path(results.dir, subtype)
  dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)
  
  save(rna, file = file.path(results.dir, expression.file.name()))
  save(met, file = file.path(results.dir, methylation.file.name()))
  
  mclapply(dimensions, function(dimension) {
    
    rna.met.anchors <- FindIntegrationAnchors(list(rna, met), k.filter = 30, dims = 1:dimension)
    rna.met.integrated <- IntegrateData(anchorset = rna.met.anchors, dims = 1:dimension)
    
    DefaultAssay(rna.met.integrated) <- "integrated"
    
    rna.met.integrated <- ScaleData(rna.met.integrated, verbose = FALSE)
    rna.met.integrated <- RunPCA(rna.met.integrated, npcs = dimension, verbose = FALSE)

    # We use the PCA for dimension reduction in order to compare the results of Seurat
    # to the other algorithms' results properly
    expression.integrated.pca.data <- 
      t(rna.met.integrated@reductions$pca@cell.embeddings)[, rna.met.integrated@meta.data$orig.ident == "RNA"]
    
    methylation.integrated.pca.data <- 
      t(rna.met.integrated@reductions$pca@cell.embeddings)[, rna.met.integrated@meta.data$orig.ident == "MET"]
    
    write.table(
      expression.integrated.pca.data,
      file = file.path(results.dir, expression.projected.file.name(dimension)),
      sep = '\t')
    
    write.table(
      methylation.integrated.pca.data, 
      file = file.path(results.dir, methylation.projected.file.name(dimension)),
      sep = '\t')
  })
}

get.seurat.projections <- function(subtype, dimension, results.dir = SEURAT.RESULTS.DIR)
{
  results.dir <- file.path(results.dir, subtype)
  
  return(list(
    expression = as.matrix(read.table(file.path(results.dir, expression.projected.file.name(dimension)))),
    methylation =  as.matrix(read.table(file.path(results.dir, methylation.projected.file.name(dimension))))
    ))
}

compute.projection.scores.for.seurat.projections <- function(
  subtype, 
  dimensions,
  results.dir = SEURAT.RESULTS.DIR,
  verbose = FALSE)
{
  cross.manifold.similarity.scores <- mcmapply(function(dimension) {
    
    cross.manifold.similarity.score.for.dimension <- tryCatch({
      
      projections <- get.seurat.projections(subtype, dimension, results.dir)
      
      cross.manifold.similarity.score <- 
        get_projection_cross_manifold_similarity_score(projections$expression, projections$methylation)
      
      if (verbose) {
        utils.log(sprintf("Seurat projection cross manifold similarity score is: %f (dimension=%d)",
                          cross.manifold.similarity.score,
                          dimension))
      }
    
      cross.manifold.similarity.score
      
    }, error = function(err) {
      
      NA
      
    }, finally = function() {})
    
    return(cross.manifold.similarity.score.for.dimension)
  },
  dimensions)
  
  names(cross.manifold.similarity.scores) <- dimensions
    
  return(cross.manifold.similarity.scores)
}

compute.projection.scores.for.intend <- 
  function(subtype, dimensions, verbose = FALSE)
{
  results.dir <- file.path(INTEND.RESULTS.DIR, subtype)
  
  projected.expression <- as.matrix(read.table(file = file.path(results.dir, expression.projected.file.name())))
  projected.methylation <- as.matrix(read.table(file = file.path(results.dir, methylation.projected.file.name())))
  projected.data.similarity <- get.inter.similarity.matrix(projected.expression, projected.methylation)
  
  if (verbose) {
    utils.log(sprintf(
      "The weights matrix score for subtype %s is: %#.2f", 
      subtype,
      get.cross.manifold.similarity.matrix.score(projected.data.similarity) * 100))
  }
  
  reduced.results <- run.cca.all.dimensions(data1 = projected.expression, data2 = projected.methylation)
  
  gc() # the scores computation uses mclapply which forks multiple child processes
  cross.manifold.similarity.scores <-
    get_projection_cross_manifold_similarity_scores_for_d_range(
      reduced.results$reduced.1,
      reduced.results$reduced.2,
      dimensions)
  
  if (verbose) {
    best.cross.manifold.similarity.score <- min(cross.manifold.similarity.scores)
    best.d <- names(cross.manifold.similarity.scores)[which.min(cross.manifold.similarity.scores)]
    utils.log(sprintf("Best CCA score: %#.3f (%s)", best.cross.manifold.similarity.score * 100, best.d))
    utils.log(sprintf("CCA Score with 30 dimensions: %#.3f",
                      cross.manifold.similarity.scores["30"] * 100))
  }
  
  return(cross.manifold.similarity.scores)
}

run.integration.algorithms.for.cancer.subtype <- function(subtype) {
  
  utils.log(glue("Gathering expression and methylation data for all algorithms for cancer type {subtype}..."))
  
  # Data for JLMA and MMD-MA is going through feature selection and scaling
  expression.methylation.data.for.jlma.and.mmd.ma.wfci <- mclapply(
    possible.num.features,
    function(num.features)
      get.gene.expression.methylation.data(subtype = subtype, number.of.feature.to.keep = num.features))
  utils.log("Gathered data for JLMA/MMD-MA WFCI")

  expression.methylation.data.for.mmd.ma <- get.expression.methylation.raw.data(subtype)
  utils.log("Gathered raw data for MMD-MA")

  expression.methylation.data.for.liger.and.seurat <-
    get.gene.expression.methylation.data.for.liger(subtype = subtype)
  utils.log("Gathered data for LIGER and Seurat v3")

  utils.log("Finished gathering expression and methylation data for all algorithms")


  utils.log("Running JLMA WFCI...")
  mclapply(1:length(possible.num.features), function(i) {
    run.manifold.alignment.for.expression.methylation.integration.with.correspondence(
      results.dir = file.path(JLMA.RESULTS.DIR, possible.num.features[i]),
      subtype = subtype,
      expression.methylation.data = expression.methylation.data.for.jlma.and.mmd.ma.wfci[[i]])
  })

  utils.log("Running LIGER...")
  run.liger.for.expression.methylation.integration(
    results.dir = LIGER.RESULTS.DIR,
    subtype = subtype,
    dimensions = dimensions.for.liger.and.seurat,
    expression.methylation.data.for.liger = expression.methylation.data.for.liger.and.seurat)

  utils.log("Running Seurat v3...")
  run.seurat.for.expression.methylation.integration(
    results.dir = SEURAT.RESULTS.DIR,
    subtype = subtype,
    dimensions = dimensions.for.liger.and.seurat,
    expression.methylation.data.for.seurat = expression.methylation.data.for.liger.and.seurat)

  utils.log("Running MMD-MA...")
  run.mmd.ma.script(
    results.dir = MMDMA.RESULTS.DIR,
    subtype = subtype,
    dimensions = dimensions.for.mmd.ma,
    expression.methylation.data = expression.methylation.data.for.mmd.ma)

  utils.log("Running MMD-MA WFCI...")
  mclapply(1:length(possible.num.features), function(i) {
    run.mmd.ma.script(
      results.dir = file.path(MMDMA.WFCI.RESULTS.DIR, possible.num.features[i]),
      subtype = subtype,
      dimensions = dimensions.for.mmd.ma,
      expression.methylation.data = expression.methylation.data.for.jlma.and.mmd.ma.wfci[[i]])
    })

  utils.log("Running INTEND using TCGA leave one out subtype prediction model...")
  run.intend.using.tcga.leave.one.subtype.out.preiction.model(subtype)
  
  utils.log(glue("Finished expression and methylation integaration with all algorithms for cancer type {subtype}"))
}

compute.projection.scores.for.integration.algorithms <- function(subtype) {
  
  utils.log(glue("Computing scores for subtype {subtype}"))
  
  min.d <- min(min(dimensions.for.mmd.ma), min(dimensions.for.liger.and.seurat))
  max.d <- max(max(dimensions.for.mmd.ma), max(dimensions.for.liger.and.seurat))
  
  utils.log("Computing scores for JLMA WFCI...")
  jlma.scores <- lapply(possible.num.features,  function(num.features)
    compute.projection.scores.for.manifold.alignment.projections(
      results.dir = file.path(JLMA.RESULTS.DIR, num.features),
      subtype = subtype,
      min.d:max.d))
  names(jlma.scores) <- sprintf("JLMA WFCI (%d)", possible.num.features)
  
  utils.log("Computing scores for LIGER...")
  liger.scores <- list(compute.projection.scores.for.liger.projections(
    subtype = subtype,
    dimensions = dimensions.for.liger.and.seurat))
  names(liger.scores) <- "LIGER"
  
  utils.log("Computing scores for Seurat v3...")
  seurat.scores <- list(compute.projection.scores.for.seurat.projections(
    subtype = subtype,
    dimensions = dimensions.for.liger.and.seurat))
  names(seurat.scores) <- "Seurat v3"
  
  utils.log("Computing scores for MMD-MA WFCI...")
  mmd.ma.scores.with.corresponsence <- lapply(possible.num.features,  function(num.features)
    compute.projection.scores.for.mmd.ma.projections(
      results.dir = file.path(MMDMA.WFCI.RESULTS.DIR, num.features),
      subtype = subtype,
      dimensions = dimensions.for.mmd.ma))
  names(mmd.ma.scores.with.corresponsence) <-
    sprintf("MMD-MA WFCI (%d)", possible.num.features)

  utils.log("Computing scores for MMD-MA...")
  mmd.ma.scores <- list(compute.projection.scores.for.mmd.ma.projections(
    results.dir = MMDMA.RESULTS.DIR,
    subtype = subtype,
    dimensions = dimensions.for.mmd.ma))
  names(mmd.ma.scores) <- "MMD-MA"
  
  utils.log("Computing scores for INTEND using trained model...")
  intend.scores <- list(INTEND = compute.projection.scores.for.intend(
    subtype = subtype,
    dimensions = min.d:max.d))
  
  utils.log("Done computing scores")
  projection.scores.list <- c(
    intend.scores,
    liger.scores,
    seurat.scores,
    mmd.ma.scores,
    mmd.ma.scores.with.corresponsence,
    jlma.scores)

  scores <- list(projection.scores.list = projection.scores.list,
              min.d = min.d,
              max.d = max.d,
              subtype = subtype)
  
  return(scores)
}

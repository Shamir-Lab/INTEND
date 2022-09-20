library(parallel)
library(umap)

if (!exists('ASSEMBLER.PATH')) {
  ASSEMBLER.PATH = 'TCGA-Assembler-2-master/TCGA-Assembler'
  source(file.path(ASSEMBLER.PATH, 'Module_A.R'))
  source(file.path(ASSEMBLER.PATH, 'Module_B.R'))
}
source('utilities.R')
source('constants.R')
source('tcga_assembler_preprocessing.R')

get.expression.methylation.data.for.gene.regions <- function(subtype, gene.regions)
{
  # get.(rna|methy).data handles missing values, hence the returned data does not contain any missing values
  gene.expression.data <- get.rna.data(subtype)
  methy.data <- get.methy.data(subtype)
  
  methylation.samples = tcga.ids.to.sample(colnames(methy.data$Data))
  rownames(methy.data$Data) <- methy.data$Des[,"REF"]
  colnames(methy.data$Data) = methylation.samples
  
  # Use TCGA-assembler function to compute per gene methylation
  gene.methylation.data <- CalculateSingleValueMethylationData(
    input = methy.data,
    regionOption = gene.regions,
    DHSOption = "Both",
    outputFileName = "processed_gene_methylation",
    outputFileFolder = file.path(preprocessed_tcga_assembler_data_dir, subtype),
    chipAnnotationFile = file.path(ASSEMBLER.PATH, "SupportingFiles/MethylationChipAnnotation.rda"))
  
  rownames(gene.methylation.data$Data) = gene.methylation.data$Des[,"GeneSymbol"]
  colnames(gene.methylation.data$Data) = methylation.samples
  
  gene.expression.samples = tcga.ids.to.sample(colnames(gene.expression.data$Data))
  rownames(gene.expression.data$Data) <- gene.expression.data$Des[,"GeneSymbol"]
  colnames(gene.expression.data$Data) = gene.expression.samples
  
  # Ignore EntrezID field in the expression data and aggregate rows with the same GeneSymbol, taking the mean
  # value per each sample.
  unique.gene.symbols <- unique(gene.expression.data$Des[, "GeneSymbol"])
  gene.expression.data.aggregated <-
    matrix(
      NA,
      nrow = length(unique.gene.symbols),
      ncol = ncol(gene.expression.data$Data)
    )
  
  for (i in 1:length(unique.gene.symbols)) {
    gene.expression.data.aggregated[i,] <- 
      colMeans(
        gene.expression.data$Data[
          which(gene.expression.data$Des[, "GeneSymbol"] == unique.gene.symbols[i]),, drop = FALSE]
      )
  }
  rownames(gene.expression.data.aggregated) <- unique.gene.symbols
  colnames(gene.expression.data.aggregated) <- colnames(gene.expression.data$Data)
  
  # 'Des' field is is not needed anymore
  gene.methylation.data <- gene.methylation.data$Data
  methy.data <- methy.data$Data
  gene.expression.data <- gene.expression.data$Data
  
  # Leave only samples which appear in both datasets
  common.samples <- intersect(gene.expression.samples, methylation.samples)
  methy.data <- methy.data[, common.samples]
  gene.methylation.data <- gene.methylation.data[, common.samples]
  gene.expression.data <- gene.expression.data[, common.samples]
  gene.expression.data.aggregated <- gene.expression.data.aggregated[, common.samples]
  
  # log transformation
  gene.expression.data <- log(1 + gene.expression.data)
  gene.expression.data.aggregated <- log(1 + gene.expression.data.aggregated)
  
  return(list(
    gene.expression.data = gene.expression.data,
    gene.expression.data.aggregated = gene.expression.data.aggregated,
    methy.data = methy.data,
    gene.methylation.data = gene.methylation.data))
}

get.expression.methylation.data.without.reverse.methylation <- function(subtype, force.update = F)
{
  expression.methylation.data.file.path <- sprintf(PREPROCESSED.EXPRESSION.METHYLATION.DATA.FORMAT,
                                                   "without_reverse_methylation",
                                                   subtype)
  
  if (file.exists(expression.methylation.data.file.path) & !force.update) {
    return(readRDS(file = expression.methylation.data.file.path))
  }
  
  #
  # We consider only the gene regions which may be associated with the gene's promoter (according to our
  # analysis in AML data we use "TSS1500", "TSS200", "5'UTR" and "1stExon").
  # CalculateSingleValueMethylationData computes the mean methylation in those regions for every gene,
  # according to Illumina_HumanMethylation450_BeadChip annotations.
  # 
  
  data <- get.expression.methylation.data.for.gene.regions(
    subtype = subtype,
    gene.regions = c("TSS1500", "TSS200", "5'UTR", "1stExon"))
  
  saveRDS(data, file = expression.methylation.data.file.path)
  
  return(data)
}

get.expression.methylation.data.for.expression.prediction <- function(subtype, force.update = FALSE)
{
  expression.methylation.data <- get.expression.methylation.data.without.reverse.methylation(
    subtype, force.update)
  
  expression.methylation.data <- list(
    subtype = subtype,
    expression = expression.methylation.data$gene.expression.data.aggregated,
    methylation = expression.methylation.data$methy.data,
    methylation.per.gene = expression.methylation.data$gene.methylation.data)
  
  return(expression.methylation.data)
}

get.common.features.for.all.subtypes <- function(force.update = FALSE)
{
  common.features.all.subtypes.file.path <-
    "CachedData/expression_methylation_regression/common_features_all.rda"
  
  if (file.exists(common.features.all.subtypes.file.path) && !force.update) {
    return(readRDS(file = common.features.all.subtypes.file.path))
  }
  
  expression.features.all.subtypes <- c()
  methylation.features.all.subtypes <- c()
  
  first.subtype <- TCGA.SUBTYPES[1]
  expression.methylation <- get.expression.methylation.data.for.expression.prediction(first.subtype)
  expression.features.all.subtypes <- rownames(expression.methylation$expression)
  methylation.features.all.subtypes <- rownames(expression.methylation$methylation)
  
  stopifnot(identical(
    colnames(expression.methylation$expression), colnames(expression.methylation$methylation)))
  
  for (subtype in setdiff(TCGA.SUBTYPES, first.subtype)) {
    expression.methylation <- get.expression.methylation.data.for.expression.prediction(subtype)
    expression.features.all.subtypes <- intersect(expression.features.all.subtypes,
                                                  rownames(expression.methylation$expression))
    methylation.features.all.subtypes <- intersect(methylation.features.all.subtypes,
                                                   rownames(expression.methylation$methylation))
    
    stopifnot(identical(
      colnames(expression.methylation$expression), colnames(expression.methylation$methylation)))
  }
  
  common.features.for.all.subtypes = list(
    expression = expression.features.all.subtypes,
    methylation = methylation.features.all.subtypes)
  
  saveRDS(common.features.for.all.subtypes, file = common.features.all.subtypes.file.path)
  return(common.features.for.all.subtypes)
}

construct.multi.subtype.expression.methylation.data <- function(subtypes)
{
  multi.subtypes.data <- mclapply(subtypes, function(subtype)
    get.expression.methylation.data.without.reverse.methylation(subtype))
  
  common.features <- get.common.features.for.all.subtypes()
  
  multi.subtypes.exp.met <- list()
  
  multi.subtypes.exp.met$gene.expression.data <- do.call(
    cbind,
    mclapply(
      multi.subtypes.data,
      function(subtype.data) subtype.data$gene.expression.data[common.features$expression, ]
    )
  )
  multi.subtypes.exp.met$gene.expression.data.aggregated <- do.call(
    cbind,
    mclapply(
      multi.subtypes.data,
      function(subtype.data) subtype.data$gene.expression.data.aggregated[common.features$expression, ]
    )
  )
  multi.subtypes.exp.met$methy.data <- do.call(
    cbind,
    mclapply(
      multi.subtypes.data,
      function(subtype.data) subtype.data$methy.data[common.features$methylation, ]
    )
  )
  
  gene.methylation.common.features <- Reduce(
    intersect,
    sapply(multi.subtypes.data, function(subtype.data) rownames(subtype.data$gene.methylation.data))
  )
  multi.subtypes.exp.met$gene.methylation.data <- do.call(
    cbind,
    mclapply(
      multi.subtypes.data,
      function(subtype.data) subtype.data$gene.methylation.data[gene.methylation.common.features, ]
    )
  )
  
  multi.subtypes.exp.met.file.path <- sprintf(PREPROCESSED.EXPRESSION.METHYLATION.DATA.FORMAT,
                                              "without_reverse_methylation",
                                              get.multi.subtypes.name(subtypes))
  
  saveRDS(multi.subtypes.exp.met, file = multi.subtypes.exp.met.file.path)
  
  
  multi.subtypes.attibution <- unlist(mclapply(
    1:length(subtypes),
    function(i) rep(subtypes[i], ncol(multi.subtypes.data[[i]]$gene.expression.data))))
  
  multi.subtypes.attribution.file.path <- sprintf(PREPROCESSED.EXPRESSION.METHYLATION.DATA.FORMAT,
                                                  SUBTYPE.ATTR.PREFIX,
                                                  get.multi.subtypes.name(subtypes))
  saveRDS(multi.subtypes.attibution, file = multi.subtypes.attribution.file.path)
}

get.expression.methylation.raw.data <- function(subtype)
{
  expression.methylation.data <- 
    get.expression.methylation.data.without.reverse.methylation(subtype = subtype)
  
  return(list(
    X = expression.methylation.data$gene.expression.data.aggregated,
    Y = expression.methylation.data$methy.data))
}

get.expression.methylation.data.before.feature.selection.and.scaling <- function(subtype) 
{
  expression.methylation.data <- 
    get.expression.methylation.data.without.reverse.methylation(subtype = subtype)
  
  gene.expression.data <- expression.methylation.data$gene.expression.data
  gene.expression.data.aggregated <- expression.methylation.data$gene.expression.data.aggregated
  methy.data <- expression.methylation.data$methy.data
  gene.methylation.data <- expression.methylation.data$gene.methylation.data
  
  # "Reverse sign" of methylation in gene.methylation.data, as we assume an anti-correlation between the gene
  # expression and the methylation level in the cpg sites which are close to the gene's promoter The
  # methylation level are expressed as beta values, hence the maximum value is 1
  gene.methylation.data <- 1 - gene.methylation.data
  methy.data <- 1 - methy.data
  
  # In order to compute the expression-methylation similarity matrix, leave only genes which are both in
  # gene.methylation.data and gene.expression.data.aggregated
  common.features <- intersect(rownames(gene.methylation.data), rownames(gene.expression.data.aggregated))
  gene.methylation.data <- gene.methylation.data[common.features, ]
  gene.expression.data.aggregated <- gene.expression.data.aggregated[common.features, ]
  
  
  return(list(
    gene.expression.data = gene.expression.data,
    gene.expression.data.aggregated = gene.expression.data.aggregated,
    methy.data = methy.data,
    gene.methylation.data = gene.methylation.data))
}

get.lccs.luad.expression.preprocessed.and.normalized <- function(force.update = FALSE)
{
  luad.expression.preprocessed.file.path <-
    'CachedData/preprocessed_expression_methylation/lccs_luad_expression.rds'
  
  if (file.exists(luad.expression.preprocessed.file.path) & !force.update) {
    return(readRDS(file = luad.expression.preprocessed.file.path))
  }
  
  luad.lccs <- read.table(
    'Data/GIS031/GSK_RSEM_rerun_expCounts_172_Tumor.tsv',
    header = TRUE, 
    sep = "\t")
  
  luad.lccs.expression <- as.matrix(subset(luad.lccs, select = -c(Ensembl.ID, Gene.symbol, Gene.type)))
  rownames(luad.lccs.expression) <- luad.lccs$Gene.symbol
  
  # Verify there are no missing values
  stopifnot(!anyNA(luad.lccs.expression))
  
  # Normalization:
  # Compute upper quantile (75th percentile) of every colmumn (sample), after removing zero read counts
  data.quantile.expressed <-
    apply(luad.lccs.expression, 2, function(x){quantile(x[x > 0], 0.75)})
  
  # Divide by the computed quantile and multiply by 1000 to get the normalized read counts
  luad.lccs.expression.norm <- do.call(cbind, lapply(1:ncol(luad.lccs.expression), function(i) {
    return(luad.lccs.expression[,i] / data.quantile.expressed[i] * 1000)
  }))
  # Restore column names
  colnames(luad.lccs.expression.norm) <- colnames(luad.lccs.expression)
  
  common.genes <- intersect(
    get.common.features.for.all.subtypes()$expression,
    rownames(luad.lccs.expression.norm)
  )
  
  luad.lccs.expression.filtered <- luad.lccs.expression.norm[
    rownames(luad.lccs.expression.norm) %in% common.genes, ]
  
  # Aggregate genes with same symbol to support inputs with gene symbols only (and not Ensemble ID)
  luad.lccs.expression.aggregated.df <- aggregate(
    luad.lccs.expression.filtered, 
    list(row.names(luad.lccs.expression.filtered)),
    mean
  )
  
  luad.lccs.expression.aggregated <- 
    as.matrix(subset(luad.lccs.expression.aggregated.df, select = -Group.1))
  rownames(luad.lccs.expression.aggregated) = luad.lccs.expression.aggregated.df$Group.1
  
  luad.lccs.expression.aggregated <- luad.lccs.expression.aggregated[common.genes, ]
  
  # Log transform
  luad.lccs.expression.aggregated <- log(1 + luad.lccs.expression.aggregated)
  
  saveRDS(luad.lccs.expression.aggregated, file = luad.expression.preprocessed.file.path)
  return(luad.lccs.expression.aggregated)
}

run.umap <- function(subtype, expression.methylation, force.update = FALSE)
{
  expression.methylation.umap.file.path <- sprintf(PREPROCESSED.EXPRESSION.METHYLATION.DATA.FORMAT,
                                                   "umap_cosine",
                                                   subtype)
  
  if (file.exists(expression.methylation.umap.file.path) & !force.update) {
    return(readRDS(file = expression.methylation.umap.file.path))
  }
  
  umap.data <- list()
  set.seed(123)
  umap.data$expression <- umap(t(expression.methylation$expression), metric = "cosine")
  set.seed(123)
  umap.data$methylation <- umap(t(expression.methylation$methylation), metric = "cosine")
  set.seed(123)
  umap.data$methylation.per.gene <- umap(t(expression.methylation$methylation.per.gene), metric = "cosine")
  
  umap.data <- lapply(umap.data, function(umap.entry) {
    umap.entry[["data"]] <- NULL
    return(umap.entry)
  })
  
  saveRDS(umap.data, file = expression.methylation.umap.file.path)
  
  return(umap.data)
}

run.umap.for.sutype <- function(subtype, force.update = FALSE)
{
  expression.methylation.umap.file.path <- sprintf(PREPROCESSED.EXPRESSION.METHYLATION.DATA.FORMAT,
                                                   "umap_cosine",
                                                   subtype)
  
  if (file.exists(expression.methylation.umap.file.path) & !force.update) {
    return(readRDS(file = expression.methylation.umap.file.path))
  }
  
  expression.methylation <- get.expression.methylation.data.for.expression.prediction(subtype)
  
  return(run.umap(subtype, expression.methylation, force.update))
}

UPSTREAM.MARGIN = 10000
DOWNSTREAM.MARGIN = 10000

IS.UNIX <- (.Platform$OS.type == "unix")

# This is only for the algorithms which don't have their own feature selection method 
# (currently JLMA and MMD-MA)
possible.num.features <- c(500, 2000)
dimensions.for.mmd.ma <- c(2,10,20,30, 40)
dimensions.for.liger.and.seurat <- c(2:40)

TCGA.SUBTYPES <- c("AML", "BLCA", "COAD", "LGG", "LIHC", "LUAD", "PAAD", "PRAD", "SARC", "SKCM", "THCA")

LIMITED.NUMBER.OF.CORES <- 4

DIMENSION.FOR.INTEGRATION.COMPARISON <- 40

SUBTYPES.FOR.MULTI.SUBTYPE.ANALYSIS <- c("COAD", "LIHC", "SARC", "SKCM")

SUBTYPE.ATTR.PREFIX <- "subtype_attribution"

#File Paths

REGRESSION.OUTPUT.TRAIN.SUBTYPES.PATH.CV.LVO.FORMAT <- 
  "CachedData/expression_methylation_regression/regression_output_train_subtypes_except_%s_%d_%d.rda"

REGRESSION.OUTPUT.TRAIN.SUBTYPES.PATH.FORMAT <- 
  "CachedData/expression_methylation_regression/regression_output_train_subtypes_%d_%d_%s.rda"

PREPROCESSED.EXPRESSION.METHYLATION.DATA.FORMAT <- "CachedData/preprocessed_expression_methylation/%s_%s"

JLMA.RESULTS.DIR <- '../jlma_results'
MMDMA.WFCI.RESULTS.DIR <- '../mmdma_results'
MMDMA.RESULTS.NO.CORRESPONDENCE.DIR <- '../mmdma_results_no_correspondence'
MMDMA.RESULTS.DIR <- '../mmdma_results_raw'
LIGER.RESULTS.DIR <- '../liger_results'
SEURAT.RESULTS.DIR <- '../seurat_results'
INTEND.RESULTS.DIR <- '../intend_results'

# File Names

LUAD.TCGA.LCCS.INTEGRATION.SUBTYPE <- "LUAD_TCGA_LCCS"
CV.LVO.SUBTYPE.NAME.FORMAT <- "CV-LVO.%s"
SUBTYPE.BASED.ON.TRAIN.SUBTYPES.FORMAT <- "%s_with_train_%s"

expression.file.name <- function() {
  return('X_exp.tsv')
}

methylation.file.name <- function() {
  return('Y_met.tsv')
}

expression.similarities.file.name <- function() {
  return('W_x_exp.tsv')
}

methylation.similarities.file.name <- function() {
  return('W_y_met.tsv')
}

expression.projected.file.name <- function(dimension = NA) {
  
  if (is.na(dimension)) {
    
    return('X_exp_projected.tsv')
    
  } else {
    
    return(sprintf('X_exp_projected_%d.tsv', dimension))
    
  }
}

methylation.projected.file.name <- function(dimension = NA) {
  
  if (is.na(dimension)) {
    
    return('Y_met_projected.tsv')
    
  } else {
    
    return(sprintf('Y_met_projected_%d.tsv', dimension))
    
  }
}

WEAK.GRAY <- 'gray85'

EXPRESSION.GROUP <- "Gene Expression"
METHYLATION.GROUP <- "DNA Methylation"
BY.MODALITY <- "Modality"
BY.SUBTYPE <- "Subtype"


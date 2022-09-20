library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)

source('constants.R')
source('gene_mapping_utils.R')
source('data_utils.R')
source('prediction.R')
source('utilities.R')
source('intend.R')
source('benchmark.R')

# library(umap)
# library(ggrepel)
# library(RColorBrewer)
# library(batchelor)
# library(FNN)
# library(class)
# library(hrbrthemes)
# library(survminer)
# library(glue)
# source('tcga_assembler_preprocessing.R')
# source('survival.R')
# source('projection_functions.R')

####################################################################################################

# Common Utils

# The MMA-MA WFCI version and JLMA WFCI are not used in most figures
# Hence, they are referenced as additional algorithms
get.projections <- function(subtype, base.algorithms = T, additional.algorithms = F)
{
  dimension <- DIMENSION.FOR.INTEGRATION.COMPARISON
  
  base <- list()
  if (base.algorithms) {
    base <- list(
      INTEND = get.intend.projections(subtype, dimension),
      LIGER = get.liger.projections(subtype, dimension),
      `Seurat v3` = get.seurat.projections(subtype, dimension),
      `MMD-MA` = get.mmd.ma.projections(subtype, dimension, MMDMA.RESULTS.DIR)
    )
  }
  
  additional <- list()
  if (additional.algorithms) {
    additional <- list(
      `MMD-MA WFCI (500)` = get.mmd.ma.projections(subtype, dimension, file.path(MMDMA.WFCI.RESULTS.DIR, "500")),
      `MMD-MA WFCI (2000)` = get.mmd.ma.projections(subtype, dimension, file.path(MMDMA.WFCI.RESULTS.DIR, "2000")),
      `JLMA WFCI (500)` = get.manifold.alignment.projections(subtype, dimension, file.path(JLMA.RESULTS.DIR, "500")),
      `JLMA WFCI (2000)` = get.manifold.alignment.projections(subtype, dimension, file.path(JLMA.RESULTS.DIR, "2000"))
      )
  }
  
  return(c(base, additional))
}

####################################################################################################

# Supplementary Figures 1 + 2 + 3

plot.model.histograms <- function()
{
  common.features <- get.common.features.for.all.subtypes()
  gene.to.cg.site.mapping <- get.filtered.gene.to.cg.sites.mapping.all.subtypes(
    common.features$expression,
    common.features$methylation,
    UPSTREAM.MARGIN,
    DOWNSTREAM.MARGIN)
  
  number.of.sites.per.gene <- sapply(gene.to.cg.site.mapping$cg_sites, length)
  
  # Supplementary Figure 1 
  plot.hiatogram.of.methylation.sites.per.gene(number.of.sites.per.gene)
  
  regression.output <- load.expression.methylation.regression.output(output.file.path = sprintf(
    REGRESSION.OUTPUT.TRAIN.SUBTYPES.PATH.CV.LVO.FORMAT, "LUAD", UPSTREAM.MARGIN, DOWNSTREAM.MARGIN))
  # The first coefficient is excluded as it is the intercept value. This is why `[c(-1),]` is used.
  number.of.sites.per.gene.after.lasso.shrinkage <-  sapply(
    regression.output$fit.parameters, function(fit) sum(
      coef(fit, s = "lambda.min")[c(-1),] != 0, na.rm = T))
  
  # Supplementary Figure 2 
  plot.hiatogram.of.methylation.sites.per.gene(number.of.sites.per.gene.after.lasso.shrinkage)
  
  r.squared.train <- regression.output$R.squared.train
  r.squared.train.top.2000.genes <- sort(r.squared.train, decreasing = T)[1:2000]
  
  # Supplementary Figure 3 
  plot.histogram.for.prediction.r.squared(r.squared.train)
  plot.histogram.for.prediction.r.squared(r.squared.train.top.2000.genes)
}

plot.hiatogram.of.methylation.sites.per.gene <- function(number.of.sites.per.gene) {
  df <- data.frame(number.of.sites = number.of.sites.per.gene)
  df.lines <- data.frame(values =  c(mean(df$number.of.sites), quantile(df$number.of.sites)[2:4]),
                         labels = c("Average", "Q1", "Median", "Q3"))
  
  dev.new(noRStudioGD = T)
  ggplot() +
    geom_histogram(data = subset(df, number.of.sites <= 100),
                   mapping = aes(x = number.of.sites),
                   binwidth = 1,
                   color = "grey20",
                   fill = "grey95") +
    geom_segment(data = df.lines,
                 aes(x = values, y = 0, xend = values, yend = Inf, color = labels),
                 size = 1.2,
                 show.legend = FALSE,
                 linetype = "dashed") +
    xlab("Number of methylation probes") +
    ylab("Number of genes") +
    geom_text(data = df.lines,
              mapping = aes(x = values, y = -30, label = labels, color = labels),
              size = 5,
              angle = 45,
              show.legend = FALSE) +
    theme_cowplot(font_size = 18)
}

plot.histogram.for.prediction.r.squared <- function(r.squared.train)
{
  df <- data.frame(R.squared.train =r.squared.train)
  df.lines <- data.frame(values =  c(mean(df$R.squared.train)),
                         labels = c("Average"))
  
  dev.new(noRStudioGD = T)
  ggplot() +
    geom_histogram(data = df,
                   mapping = aes(x = R.squared.train),
                   binwidth = 0.02,
                   color = "grey20",
                   fill = "grey95") +
    geom_segment(data = df.lines,
                 aes(x = values, y = 0, xend = values, yend = Inf, color = labels),
                 size = 1.2,
                 show.legend = FALSE,
                 linetype = "dashed") +
    xlab("R squared (training set)") +
    ylab("Number of genes") +
    geom_text(data = df.lines,
              mapping = aes(x = values, y = -12, label = sprintf("%s\n%.02f", labels, values), color = labels),
              size = 4.8,
              show.legend = FALSE) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
    theme_cowplot(font_size = 16)
}

####################################################################################################

# Supplementary Figures 5

# scores.data: the result of compute.projection.scores.for.integration.algorithms (utils.R)
# for a specified subtype

plot.legend.of.score.vs.dimensions <- function(scores.data)
{
  dev.new(width = 1600, height = 500, unit = "px", noRStudioGD = T)
  plot(NULL ,xaxt = 'n',yaxt = 'n',bty = 'n',ylab = '',xlab = '', xlim = 0:1, ylim = 0:1)
  
  projection.scores.list <- scores.data$projection.scores.list
  legend.text <- names(projection.scores.list)
  
  legend(
    "topleft",
    legend = legend.text,
    col = 1:length(projection.scores.list),
    ncol = 4,
    pch = 19,
    cex = 1.2,
    bg = "#f9f9f9")
}


plot.projection.scores.vs.dimensions <- function
(scores.data,
 show.legend = F, 
 show.best.score = F,
 pos = 1,
 new.window = T) 
{
  if (new.window) {
    dev.new(width = 1600, height = 1000, unit = "px", noRStudioGD = T)
  }
  
  projection.scores.list <- scores.data$projection.scores.list
  min.d <- scores.data$min.d
  max.d <- scores.data$max.d
  subtype <- scores.data$subtype
  
  
  if (show.legend) {
    margines = c(5, 5, 4, 15)
    par(mar = margines)
  }
  
  plot(NULL,
       type = "n",
       xlab = "Dimension",
       ylab = "FOSCTTM (%)",
       main = subtype,
       xlim = c(min.d, max.d), 
       ylim = c(0,60),
       cex.lab = 1.3,
       cex.main = 1.5,
       pch = 19)
  
  
  for (i in 1:length(projection.scores.list)) {
    
    scores <- projection.scores.list[[i]]
    dimensions <- as.integer(names(scores))
    
    indices.of.relevant.dimensions <- (dimensions >= min.d & dimensions <= max.d)
    
    points(dimensions[indices.of.relevant.dimensions],
           scores[indices.of.relevant.dimensions] * 100,
           col = i,
           pch = 19,
           cex = 0.8)
  }
  
  par(xpd = FALSE)
  grid()
  
  if (show.legend) {
    par(mar = margines, xpd = TRUE)
    
    if (show.best.score) {
      legend.text <- paste(
        names(projection.scores.list),
        lapply(projection.scores.list, function(x) {
          min.score <- min(x, na.rm = TRUE)
          d.for.min.score <- names(x)[which(x == min.score)]
          sprintf("\n(%#.2f)", min.score*100)
        }))
    }
    else {
      legend.text <- names(projection.scores.list)
    }
    
    legend(
      "right",
      inset = c(-0.25, 0),
      legend = legend.text,
      col = 1:length(projection.scores.list),
      pch = 19,
      cex = 1.2,
      y.intersp = 2,
      bg = "#f9f9f9")
  }
}

####################################################################################################

# Figure 2 + Supplementary Figure 17

plot.projection.score.per.sample.for.all.algorithms <- function(
  subtype,
  dimension,
  show.axis.titles = F,
  show.legend = F,
  additional.algorithms = T
)
{
  projections.list <- get.projections(subtype, dimension, additional.algorithms = additional.algorithms)
  
  scores.per.sample <- lapply(projections.list, function(projections)
    get_projection_cross_manifold_similarity_score_per_sample(
      X_projected = projections$expression,
      Y_projected = projections$methylation
    )
  )
  
  projections.algorithm <- unlist(lapply(1:length(scores.per.sample), function(i) 
    rep(names(scores.per.sample)[i], length(scores.per.sample[[i]]))))
  scores.aggregated <- unlist(scores.per.sample) * 100
  
  projections.algorithm <- factor(projections.algorithm, levels = names(scores.per.sample))
  
  scores.per.sample.df <- data.frame(algorithm = projections.algorithm, scores = scores.aggregated)
  
  boxplots <- ggplot(scores.per.sample.df, aes(x = algorithm, y = scores, color = algorithm, fill = algorithm)
                     ) +
    geom_boxplot() +
    ylab("FOSCTTM (%)") +
    xlab("Algorithm") +
    theme(
      text = element_text(size = 15),
      axis.text.x = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.title = element_blank(),
      legend.position = ifelse(test = show.legend, yes = "right", no = "none")
    ) + scale_fill_brewer(palette = "Pastel2") + scale_color_brewer(palette = "Dark2")
  
  if (!show.axis.titles) {
    boxplots <- boxplots + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  }
  
  return(boxplots)
}

# Figure 2
plot.projection.score.per.sample.for.all.subtypes <- function()
{
  plot.list <- lapply(
    TCGA.SUBTYPES, 
    function(subtype) plot.projection.score.per.sample.for.all.algorithms(
      subtype,
      DIMENSION.FOR.INTEGRATION.COMPARISON
    ) + ggtitle(subtype)
  )
  
  plot <- plot_grid(plotlist = plot.list, nrow = 4)
  
  
  y.grob <- textGrob("FOSCTTM (%)", 
                     gp = gpar(fontsize = 15), rot = 90)
  
  x.grob <- textGrob("Algorithm", 
                     gp = gpar(fontsize = 18))
  
  dev.new(noRStudioGD = T)
  grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob))
}

# Supplementary Figure 17
plot.projection.score.per.sample.for.multi.subtype.example <- function() {
  multi.subtype.name <- get.multi.subtypes.name(SUBTYPES.FOR.MULTI.SUBTYPE.ANALYSIS)
  plot <- plot.projection.score.per.sample.for.all.algorithms(
    multi.subtype.name,
    DIMENSION.FOR.INTEGRATION.COMPARISON,
    show.axis.titles = T,
    show.legend = T,
    additional.algorithms = F
  ) 

  dev.new(noRStudioGD = T)
  plot(plot)
}

####################################################################################################

# Supplementary Figure 4

plot.expression.methylation.correlations.per.genomic.region <- 
  function(correlations.per.region.for.all.subtypes)
  {
    keys <- unique(unlist(lapply(correlations.per.region.for.all.subtypes, names)))
    correlations.per.region <- setNames(
      do.call(mapply,c(FUN = c, lapply(correlations.per.region.for.all.subtypes, `[`, keys))),
      keys
    )
    
    # Remove NA values
    correlations.per.region <- lapply(correlations.per.region, function(cor) cor[!is.na(cor)])
    
    correlations.concatenated <- unlist(correlations.per.region)
    regions.attribution <- rep(names(correlations.per.region), times = sapply(correlations.per.region, length))
    regions.attribution <- factor(regions.attribution, levels = names(correlations.per.region))
    
    df <- data.frame(regions = regions.attribution, correlations = correlations.concatenated)
    
    dev.new(noRStudioGD = T)
    boxplots <- ggplot(df, aes(x = regions, y = correlations, color = regions)) +
      geom_violin(show.legend = FALSE) +
      geom_hline(yintercept = 0, size = 1.2, linetype = "dashed") +
      stat_summary(fun = mean,
                   geom = "crossbar",
                   width = 0.75, 
                   mapping = aes(color = regions),
                   show.legend = FALSE) +
      geom_point(data = subset(df, FALSE), size = 5, shape = 15) + # Hack for legend key control
      guides(color = guide_legend(title = NULL, label.theme = element_text(size = 16))) +
      xlab("Gene Region") + 
      ylab("Expression-Methylation Correlations per Gene") +
      theme(
        text = element_text(size = 18),
        panel.grid.minor = element_line(colour = "white", size = 0.5),
        panel.grid.major = element_line(colour = "white", size = 0.8)
      ) +
      scale_y_continuous(minor_breaks = seq(-1 , 1, 0.05), breaks = seq(-1, 1, 0.1)) +
      scale_fill_brewer(palette = "Set1") 
    
    plot(boxplots)
  }

plot.expression.methylation.correlations.per.genomic.region.for.all.subtypes.in.plot.grid <-
  function(correlations.per.region.for.all.subtypes)
  {
    plot.list <- lapply(1:length(TCGA.SUBTYPES), function(i)
    {
      plot.expression.methylation.correlations.per.genomic.region.for.subtype(
        correlations.per.region.for.all.subtypes[[i]],
        TCGA.SUBTYPES[i])
    }
    )
    
    plot <- plot_grid(plotlist = plot.list)
    
    
    y.grob <- textGrob("Expression-Methylation Correlations per Gene", 
                       gp = gpar(fontsize = 18), rot = 90)
    
    x.grob <- textGrob("Gene Region", 
                       gp = gpar(fontsize = 15))
    
    dev.new(noRStudioGD = T)
    grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob))
  }

plot.expression.methylation.correlations.per.genomic.region.for.subtype <- 
  function(correlations.per.region, subtype)
  {
    # Remove NA values
    correlations.per.region <- lapply(correlations.per.region, function(cor) cor[!is.na(cor)])
    
    correlations.concatenated <- unlist(correlations.per.region)
    regions.attribution <- rep(names(correlations.per.region), times = sapply(correlations.per.region, length))
    regions.attribution <- factor(regions.attribution, levels = names(correlations.per.region))
    
    df <- data.frame(regions = regions.attribution, correlations = correlations.concatenated)
    
    boxplots <- ggplot(df, aes(x = regions, y = correlations, color = regions)) +
      geom_violin(show.legend = FALSE) +
      geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
      stat_summary(fun = mean,
                   geom = "crossbar",
                   width = 0.75, 
                   mapping = aes(color = regions),
                   show.legend = FALSE) +
      ggtitle(subtype) +
      theme(
        text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(size = 15),
        panel.grid.minor = element_line(colour = "white", size = 0.5),
        panel.grid.major = element_line(colour = "white", size = 0.8)
      ) +
      scale_y_continuous(minor_breaks = seq(-1 , 1, 0.1), breaks = seq(-1, 1, 0.5)) +
      scale_fill_brewer(palette = "Set1") 
    
    return(boxplots)
  }


####################################################################################################

create.original.datasets.data.for.plot <- function(
  expression.reduced,
  methylation.reduced,
  reduction.method,
  multi.subtypes.attribution.expression = NULL,
  multi.subtypes.attribution.methylation = multi.subtypes.attribution.expression
)
{
  rownames(expression.reduced) <- paste(rownames(expression.reduced), "1", sep = "_")
  rownames(methylation.reduced) <- paste(rownames(methylation.reduced), "2", sep = "_")
  
  modalities <- c(rep(EXPRESSION.GROUP, nrow(expression.reduced)),
                  rep(METHYLATION.GROUP, nrow(methylation.reduced)))
  
  sample.num <- c(1:nrow(expression.reduced), 1:nrow(methylation.reduced))
  
  data.for.plot <- list(
    reduced.data = rbind(expression.reduced, methylation.reduced), 
    reduction.method = reduction.method,
    sample.num = sample.num
  )
  data.for.plot[[BY.MODALITY]] = modalities
  
  if (!is.null(multi.subtypes.attribution.expression)) {
    
    stopifnot(length(multi.subtypes.attribution.expression) == nrow(expression.reduced) &&
                length(multi.subtypes.attribution.methylation) == nrow(methylation.reduced))
    
    data.for.plot[[BY.SUBTYPE]] = c(
      multi.subtypes.attribution.expression,
      multi.subtypes.attribution.methylation
    )
  }

  return(data.for.plot)
}

get.inetgrated.results.for.plot <- function(
  projections,
  multi.subtypes.attribution.expression = NULL,
  multi.subtypes.attribution.methylation = multi.subtypes.attribution.expression
)
{
  projected.expression <- projections$expression
  projected.methylation <- projections$methylation
  
  colnames(projected.expression) <- paste(colnames(projected.expression), "1", sep = "_")
  colnames(projected.methylation) <- paste(colnames(projected.methylation), "2", sep = "_")
  
  integrated.data <- cbind(projected.expression, projected.methylation)
  modalities <- c(rep(EXPRESSION.GROUP, ncol(projected.expression)),
                  rep(METHYLATION.GROUP, ncol(projected.methylation)))
  
  sample.num <- c(1:ncol(projected.expression), 1:ncol(projected.methylation))
  
  if (!is.null(multi.subtypes.attribution.expression)) {
    
    stopifnot(length(multi.subtypes.attribution.expression) == ncol(projected.expression) &&
                length(multi.subtypes.attribution.methylation) == ncol(projected.methylation))
  }
  
  
  return(list(integrated.data = integrated.data,
              modalities = modalities,
              sample.num = sample.num,
              subtype.attribution = c(
                multi.subtypes.attribution.expression,
                multi.subtypes.attribution.methylation
                )
              ))
}

preprocess.plot.data <- function(plot.data, dims, color.by, shape.by = NA)
{
  df <- data.frame(row.names = rownames(plot.data$reduced.data))
  
  c1 <- glue("{plot.data$reduction.method}_1")
  c2 <- glue("{plot.data$reduction.method}_2")
  
  df[[c1]] <- plot.data$reduced.data[, dims[1]]
  df[[c2]] <- plot.data$reduced.data[, dims[2]]
  df[[color.by]] <- plot.data[[color.by]]
  if (!is.na(shape.by)) {
    df$shape <- as.factor(plot.data[[shape.by]])
  }
  else {
    df$shape <- as.factor(rep(19, times = nrow(df)))
  }
  df$sample.num <- plot.data$sample.num
  
  return(list(df = df, x = c1, y = c2, color.by = color.by))
}

preprocess.plot.data.with.facets <- function(plot.data, dims, color.by, facets.group.by)
{
  data <- preprocess.plot.data(plot.data, dims, color.by)
  data$df[[facets.group.by]] <- plot.data[[facets.group.by]]
  return(data)
}

plot.as.seperate.datasets.2d <- function(
  subtype,
  plot.data,
  dims = c(1,2),
  color.by = BY.MODALITY,
  facets.group.by = BY.MODALITY,
  new.window = FALSE,
  continuous.scale = FALSE) 
{
  if (new.window) {
    dev.new(height = 500, width = 1000, unit = "px",noRStudioGD = TRUE)
  }
  
  data <- preprocess.plot.data.with.facets(plot.data, dims, color.by, facets.group.by)
  
  plot <-  
    ggplot(data = data$df, mapping = aes_string(x = data$x, y = data$y, color = data$color.by)) +
    facet_wrap(as.formula(paste("~", facets.group.by)), scales = "free") +
    geom_point(
      size = 1,
      # shape = ".",
      show.legend = F) +
    theme_cowplot() +
    panel_border() + 
    theme(
      text = element_text(size = 10),
      axis.text = element_text(size = 10),
      strip.text.x = element_text(size = (14))
    )
  
  if (continuous.scale) {
    plot <- plot + scale_colour_gradient2(mid = "yellow", high = "red", na.value = "grey80")
  }
  
  if (new.window) {
    plot(plot)
  } else {
    return(plot)
  }
}

plot.as.seperate.datasets.2d.with.labeled.samples <- function(
  subtype,
  plot.data, 
  dims = c(1,2), 
  color.by = BY.MODALITY,
  facets.group.by = BY.MODALITY,
  new.window = FALSE) 
{
  if (new.window) {
    dev.new(height = 500, width = 1000, unit = "px",noRStudioGD = TRUE)
  }  
  data <- preprocess.plot.data.with.facets(plot.data, dims, color.by, facets.group.by)
  
  set.seed(123456) # For the reproduction of the data points sampling
  sampled.points <- sample(unique(data$df$sample.num), 10)
  
  plot <-  
    ggplot(data = data$df, mapping = aes_string(x = data$x, y = data$y, color = data$color.by)) +
    facet_wrap(as.formula(paste("~", facets.group.by)), scales = "free") +
    geom_point(color = WEAK.GRAY, show.legend = FALSE) +
    geom_point(
      data =  subset(data$df, sample.num %in% sampled.points),
      show.legend = FALSE) +
    geom_text_repel(
      # data = data$df,
      # mapping = aes(label = sample.num),
      # size = 2.5,
      data =  subset(data$df, sample.num %in% sampled.points),
      mapping = aes(label = rep(1:length(sampled.points), 2)),
      size = 5,
      show.legend = FALSE) +
    theme_cowplot() + 
    panel_border() + 
    theme(
      text = element_text(size = 10),
      axis.text = element_text(size = 10),
      strip.text.x = element_text(size = (14))
    )
  
  if (new.window) {
    plot(plot)
  } else {
    return(plot)
  }
}

plot.as.integrated.dataset.2d <- function(
  algorithm,
  plot.data,
  dims = c(1,2),
  color.by = BY.MODALITY,
  shape.by = NA,
  new.window = FALSE,
  show.legend = FALSE,
  show.title = TRUE,
  continuous.scale = FALSE)  
{
  if (new.window) {
    dev.new(height = 500, width = 600, unit = "px",noRStudioGD = TRUE)
  }

  data <- preprocess.plot.data(plot.data = plot.data, dims = dims, color.by = color.by, shape.by = shape.by)
  
  plot <- 
    ggplot(data = data$df, aes_string(x = data$x, y = data$y, color = data$color.by)) +
    geom_point(data = subset(data$df, FALSE), size = 5, shape = 15) + # Hack for legend key control
    guides(color = guide_legend(title = NULL, label.theme = element_text(size = 16))) +
    geom_point(
      size = 1,
      show.legend = FALSE,
      aes(shape = shape)) +
    scale_shape_manual(values = c(19, 0, 3, 15)) +
    theme_cowplot() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      text = element_text(size = 10),
      axis.text = element_text(size = 10),
    )
  
  if (continuous.scale) {
    plot <- plot + scale_colour_gradient2(mid = "yellow", high = "red", na.value = "grey80")
  }
  
  if (show.title) {
    plot <- plot + ggtitle(algorithm)
  }
  if (!show.legend) {
    plot <- plot + theme(legend.position = "none")
  }
  
  if (new.window) {
    plot(plot)
  } else {
    return(plot)
  }
}

plot.as.integrated.dataset.2d.with.labeled.samples <- function(
  algorithm,
  plot.data, 
  dims = c(1,2),
  color.by = BY.MODALITY,
  new.window = FALSE,
  show.legend = FALSE) 
{
  if (new.window) {
    dev.new(height = 500, width = 600, unit = "px",noRStudioGD = TRUE)
  }
  
  data <- preprocess.plot.data(plot.data = plot.data, dims = dims, color.by = color.by)
  
  set.seed(123456) # For the reproduction of the data points sampling
  sampled.points <- sample(unique(data$df$sample.num), 10)
  
  plot <- 
    ggplot(data = data$df, aes_string(x = data$x, y = data$y, color = data$color.by)) +
    geom_point(data = subset(data$df, FALSE), size = 5, shape = 15) + # Hack for legend key control
    guides(color = guide_legend(title = NULL, label.theme = element_text(size = 16))) +
    geom_point(color = WEAK.GRAY, show.legend = FALSE) +
    geom_point(
      data =  subset(data$df, sample.num %in% sampled.points),
      show.legend = FALSE) +
    geom_text_repel(
      # data = data$df,
      # mapping = aes(label = sample.num),
      # size = 2.5,
      data =  subset(data$df, sample.num %in% sampled.points),
      mapping = aes(label = rep(1:length(sampled.points), 2)),
      size = 5,
      show.legend = FALSE) +
    theme_cowplot() + 
    theme(
      text = element_text(size = 10),
      axis.text = element_text(size = 10),
    )
  
  
  if (!show.legend) {
    plot <- plot + theme(legend.position = "none")
  }
  
  if (new.window) {
    plot(plot)
  } else {
    return(plot)
  }
}

get.filtered.indices.with.survival.data <- function(
  survival.data,
  exp.met.data,
  surv.plot.data.orig,
  surv.plot.data.integrated
) 
{
  expression.samples <- colnames(exp.met.data$expression)
  methylation.samples <- colnames(exp.met.data$methylation)
  
  # Filter only primary tumor samples, to avoid more than one sample for the same patient
  filtered.exp.samples <- filter.only.tumor.samples(expression.samples)
  filtered.met.samples <- filter.only.tumor.samples(methylation.samples)
  filtered.exp.patients <- get.short.tcga.samples.names(filtered.exp.samples)
  filtered.met.patients <- get.short.tcga.samples.names(filtered.met.samples)
  
  # Verify that the filtered patients list does not contain duplicates
  stopifnot(anyDuplicated(filtered.exp.patients) == 0)
  stopifnot(anyDuplicated(filtered.met.patients) == 0)
  
  TCGA.PATIENT.LENGTH <- 12
  filtered.exp.indices <- which(rownames(surv.plot.data.orig$reduced.data) %in% filtered.exp.samples)
  filtered.exp.indices.2 <- which(rownames(surv.plot.data.integrated$reduced.data) %in% filtered.exp.samples)
  stopifnot(identical(filtered.exp.indices, filtered.exp.indices.2))
  
  filtered.met.indices <- which(rownames(surv.plot.data.orig$reduced.data) %in% filtered.met.samples)
  filtered.met.indices.2 <- which(rownames(surv.plot.data.integrated$reduced.data) %in% filtered.met.samples)
  stopifnot(identical(filtered.met.indices, filtered.met.indices.2))
  
  filtered.exp.indices <- filtered.exp.indices[substr(rownames(surv.plot.data.orig$reduced.data)[
    filtered.exp.indices], 1, TCGA.PATIENT.LENGTH) %in% names(survival.data)]
  filtered.met.indices <- filtered.met.indices[substr(rownames(surv.plot.data.orig$reduced.data)[
    filtered.met.indices], 1, TCGA.PATIENT.LENGTH) %in% names(survival.data)]
  
  stopifnot(length(intersect(filtered.exp.indices, filtered.met.indices)) == 0)
  
  return(list(
    exp = filtered.exp.indices,
    met = filtered.met.indices - length(expression.samples),
    both = c(filtered.exp.indices, filtered.met.indices)
  ))
}

filter.plot.data <- function(plot.data, indices, additional.1d.array.fields = NULL)
{
  stopifnot(length(plot.data) == 4 + length(additional.1d.array.fields))
  plot.data$reduced.data <- plot.data$reduced.data[indices, ]
  plot.data$sample.num <- plot.data$sample.num[indices]
  plot.data[[BY.MODALITY]] <- plot.data[[BY.MODALITY]][indices]
  for (field in additional.1d.array.fields)
  {
    plot.data[[field]] <- plot.data[[field]][indices]
  }
  return(plot.data)
}

get.subtype.classification.with.survival.analysis.for.skcm <- function()
{
  repetitions <- 30
  p.val <- matrix(nrow = 11, ncol = repetitions)
  p.val.full <- list()
  best.num.clusters <- matrix(nrow = 11, ncol = repetitions)
  num.samples.with.survivial <- matrix(nrow = 2, ncol = repetitions)
  rownames(num.samples.with.survivial) <- c("exp", "met")
  
  for (i in 1:repetitions) {
    subtype <- sprintf("SKCM_SPLIT_%d", i)
    algorithm <- "INTEND"
    
    utils.log(sprintf("Collecting data for subtype: %s", subtype))
    
    exp.met.data <- get.exp.met.data.for.cluster.plot(subtype, expect.same.samples.for.exp.met = FALSE)
    
    original.datasets.for.plot <- get.original.datasets.for.plot.with.umap(subtype)
    intend.projections <- get.intend.projections(subtype, DIMENSION.FOR.INTEGRATION.COMPARISON)
    
    colnames(intend.projections$expression) <-  
      gsub(pattern = "[.]", replacement = "-", x = colnames(intend.projections$expression))
    colnames(intend.projections$methylation) <-  
      gsub(pattern = "[.]", replacement = "-", x = colnames(intend.projections$methylation))
    
    integrated.data.for.plot <- get.integrated.data.for.plot.with.umap(intend.projections)
    
    orig.data.plot <- plot.as.seperate.datasets.2d(subtype, original.datasets.for.plot)
    integrated.plot <- plot.as.integrated.dataset.2d(
      algorithm,
      integrated.data.for.plot, 
      show.legend = T,
      show.title = F
    )
    common.legend <- get_legend(integrated.plot + theme(legend.direction = "horizontal", legend.justification = "center"))
    integrated.plot <- integrated.plot + theme(legend.position = "none")
    
    clinical.data <- get.clinical.data("SKCM")
    survival.data <- get.survival.data(clinical.data)
    surv.plot.data.orig <- original.datasets.for.plot
    surv.plot.data.integrated <- integrated.data.for.plot
    filtered.indices <- get.filtered.indices.with.survival.data(
      survival.data,
      exp.met.data,
      surv.plot.data.orig,
      surv.plot.data.integrated)
    
    num.samples.with.survivial["exp", i] <- length(filtered.indices$exp)
    num.samples.with.survivial["met", i] <- length(filtered.indices$met)
    
    exp.met.data.filtered <- list(expression = exp.met.data$expression[, filtered.indices$exp],
                                  methylation = exp.met.data$methylation[, filtered.indices$met])
    intend.projections.filtered <- list(expression = intend.projections$expression[, filtered.indices$exp],
                                  methylation = intend.projections$methylation[, filtered.indices$met])
    
    surv.plot.data.orig <- filter.plot.data(surv.plot.data.orig, filtered.indices$both)
    surv.plot.data.integrated <- filter.plot.data(surv.plot.data.integrated, filtered.indices$both)
    
    skcm.res <- get.clustering.res.for.subtype(
      subtype,
      surv.plot.data.orig,
      surv.plot.data.integrated,
      survival.data,
      intend.projections.filtered,
      exp.met.data.filtered)
    
    best.num.clusters[, i] <- skcm.res$best.num.of.clusters
    p.val.full[[i]] <- skcm.res$pvalues
    p.val[, i] <- skcm.res$pvalues[cbind(seq_len(nrow(skcm.res$pvalues)), skcm.res$best.num.of.clusters - 1)]
    rownames(best.num.clusters) = rownames(p.val) = rownames(skcm.res$pvalues)
    
    cat(sapply(rownames(skcm.res$pvalues), function(x) sprintf(
      "%s:\t%.03f (%d)\n",
      x,
      skcm.res$pvalues[x, as.character(skcm.res$best.num.of.clusters[x])],
      skcm.res$best.num.of.clusters[x]
    )), sep = "")
  }
  
  num.clusters.exp.orig <- skcm.res$best.num.of.clusters["kmeans.origin.exp"]
  num.clusters.met.orig <- skcm.res$best.num.of.clusters["kmeans.origin.met"]
  num.clusters.met.by.exp <- skcm.res$best.num.of.clusters["met.proj.by.exp.orig.umap.none"]

  clusters.exp.orig <- skcm.res$clustering.results[[toString(num.clusters.exp.orig)]]$kmeans$kmeans.origin.exp$cluster
  clusters.met.orig <- skcm.res$clustering.results[[toString(num.clusters.met.orig)]]$kmeans$kmeans.origin.met$cluster
  clusters.met.by.exp <- skcm.res$clustering.results[[
    toString(num.clusters.met.by.exp)]]$kmeans$met.proj.by.exp.orig.umap.none$cluster
  
  TCGA.PATIENT.LENGTH <- 12
  exp.orig.patient.names <- substr(names(clusters.exp.orig), 1, TCGA.PATIENT.LENGTH)
  met.orig.patient.names <- substr(names(clusters.met.orig), 1, TCGA.PATIENT.LENGTH)
  met.by.exp.patient.names <- substr(names(clusters.met.by.exp), 1, TCGA.PATIENT.LENGTH)
  
  surv.vs.clusters.exp.orig <- data.frame(
    surv = survival.data[exp.orig.patient.names], clusters = clusters.exp.orig)
  surv.vs.clusters.met.orig <- data.frame(
    surv = survival.data[met.orig.patient.names], clusters = clusters.met.orig)
  surv.vs.clusters.met.by.exp <- data.frame(
    surv = survival.data[met.by.exp.patient.names], clusters = clusters.met.by.exp)
  
  surv.analysis.exp.orig <- get.empirical.surv(surv.vs.clusters.exp.orig$surv, surv.vs.clusters.exp.orig$clusters)
  surv.analysis.met.orig <- get.empirical.surv(surv.vs.clusters.met.orig$surv, surv.vs.clusters.met.orig$clusters)
  surv.analysis.met.by.exp <- get.empirical.surv(surv.vs.clusters.met.by.exp$surv, surv.vs.clusters.met.by.exp$clusters)

  surv.analysis.exp.orig$pvalue
  surv.analysis.met.orig$pvalue
  surv.analysis.met.by.exp$pvalue
  
  fit.exp.orig <- survfit(surv ~ clusters, data = surv.vs.clusters.exp.orig)
  fit.met.orig <- survfit(surv ~ clusters, data = surv.vs.clusters.met.orig)
  fit.met.by.exp <- survfit(surv ~ clusters, data = surv.vs.clusters.met.by.exp)
  
  plot.data.clusters.exp.orig <- 
    skcm.res$clustering.results[[toString(num.clusters.exp.orig)]]$plot.data$plot.data.clustered.origin
  exp.indices <- which(plot.data.clusters.exp.orig[[BY.MODALITY]] == EXPRESSION.GROUP)
  plot.data.clusters.exp.orig <- filter.plot.data(plot.data.clusters.exp.orig, exp.indices, "cluster")

  plot.data.clusters.met.orig <- 
    skcm.res$clustering.results[[toString(num.clusters.met.orig)]]$plot.data$plot.data.clustered.origin
  met.indices <- which(plot.data.clusters.met.orig[[BY.MODALITY]] == METHYLATION.GROUP)
  plot.data.clusters.met.orig <- filter.plot.data(plot.data.clusters.met.orig, met.indices, "cluster")

  plot.data.clusters.met.by.exp <- 
    skcm.res$clustering.results[[
      toString(num.clusters.met.by.exp)]]$plot.data$plot.data.met.proj.by.exp.orig.umap.none

  color.by.for.plots <- "cluster"
  shape.by <- NA
  show.legend = FALSE
  
  plot.clustered.exp.orig <- 
    plot.as.integrated.dataset.2d(
      algorithm,
      plot.data.clusters.exp.orig,
      color.by = color.by.for.plots,
      shape.by = shape.by,
      show.legend = show.legend,
      show.title = F)
  plot.clustered.exp.orig <- plot.clustered.exp.orig + scale_color_brewer(palette = "Paired")
  
  plot.clustered.met.orig <- 
    plot.as.integrated.dataset.2d(
      algorithm, 
      plot.data.clusters.met.orig,
      color.by = color.by.for.plots,
      shape.by = shape.by,
      show.legend = show.legend,
      show.title = F)
  plot.clustered.met.orig <- plot.clustered.met.orig + scale_color_brewer(palette = "Dark2")
  
  plot.clustered.met.by.exp <-
    plot.as.integrated.dataset.2d(
      algorithm,
      plot.data.clusters.met.by.exp,
      color.by = color.by.for.plots,
      shape.by = shape.by,
      show.legend = show.legend,
      show.title = F)
  plot.clustered.met.by.exp <- plot.clustered.met.by.exp + scale_color_brewer(palette = "Paired")
  
  
  plot.survival.exp.orig <- ggsurvplot(fit.exp.orig)$plot + 
    xlab("Time (days)") + theme(legend.position = "none") + scale_color_brewer(palette = "Paired")
  
  plot.survival.met.orig <- ggsurvplot(fit.met.orig)$plot + 
    xlab("Time (days)") + theme(legend.position = "none") + scale_color_brewer(palette = "Dark2")
  
  plot.survival.met.by.exp <- ggsurvplot(fit.met.by.exp)$plot + 
    xlab("Time (days)") + theme(legend.position = "none") + scale_color_brewer(palette = "Paired")
  

  tot.withinss.exp.orig.df <- data.frame(tot.withinss = skcm.res$tot.withinss["kmeans.origin.met", ],
                                     num.clusters = as.integer(colnames(skcm.res$tot.withinss)))
  tot.withinss.met.orig.df <- data.frame(tot.withinss = skcm.res$tot.withinss["kmeans.origin.met", ],
                                     num.clusters = as.integer(colnames(skcm.res$tot.withinss)))

  highlighted.num.clusters.exp.orig <- tot.withinss.exp.orig.df %>% filter(num.clusters == num.clusters.exp.orig)
  highlighted.num.clusters.met.orig <- tot.withinss.met.orig.df %>% filter(num.clusters == num.clusters.met.orig)

  tot.withinss.exp.orig.plot <- ggplot(tot.withinss.exp.orig.df, mapping = aes(num.clusters, tot.withinss)) +
    geom_point() +
    geom_point(data = highlighted.num.clusters.exp.orig, color = 'red', size = 3) +
    ylab("Total within-cluster sum of squares") + 
    xlab("Number of clusters") +
    theme_bw()
  
  tot.withinss.met.orig.plot <- ggplot(tot.withinss.met.orig.df, mapping = aes(num.clusters, tot.withinss)) +
    geom_point() +
    geom_point(data = highlighted.num.clusters.met.orig, color = 'red', size = 3) +
    ylab("Total within-cluster sum of squares") + 
    xlab("Number of clusters") +
    theme_bw()
  
  plot.row.1 <- plot_grid(
    orig.data.plot,
    integrated.plot, 
    labels = c('A','B'),
    nrow = 1,
    rel_widths = c(2, 1)
  )
  plot.row.2 <- plot_grid(
    plot.clustered.met.orig,
    plot.clustered.exp.orig,
    plot.clustered.met.by.exp,
    labels = c('C', 'D', 'E'),
    rel_widths = c(1, 1, 1),
    nrow = 1
  )
  
  plot.row.3 <- plot_grid(tot.withinss.met.orig.plot, tot.withinss.exp.orig.plot, labels = c('F', 'G'), nrow = 1)
  
  dev.new(height = 900, width = 900, unit = "px", noRStudioGD = T)
  plot_grid(plot.row.1, plot.row.2, plot.row.3, ncol = 1, 
            rel_heights = c(1, 1, 1.5))
  
  dev.new(height = 300, width = 900, unit = "px", noRStudioGD = T)
  plot_grid(plot.survival.met.orig, plot.survival.exp.orig, plot.survival.met.by.exp, labels = c('A', 'B', 'C'), nrow = 1)

}

get.met.proj.by.exp.orig.results <- function(
  exp.proj.data,
  met.proj.data,
  clustered.origin.data,
  met.indices,
  surv.plot.data.integrated,
  neighbors.to.consider
) {
  met.proj.data.clustering.by.exp.orig.clustering <- knn(
    train = exp.proj.data,
    test = met.proj.data, 
    cl = clustered.origin.data$kmeans.exp$cluster, 
    k = neighbors.to.consider
  )
  names(met.proj.data.clustering.by.exp.orig.clustering) <- rownames(met.proj.data)
  met.proj.data.clustering.by.exp.orig.plot.data <- add.cluster.info.to.plot.data(
    filter.plot.data(surv.plot.data.integrated, met.indices),
    met.proj.data.clustering.by.exp.orig.clustering
  )
  met.proj.clustering.by.exp.orig.fake.kmeans <- list(
    tot.withinss = clustered.origin.data$kmeans.exp$tot.withinss,
    cluster = met.proj.data.clustering.by.exp.orig.clustering)
  
  return(list(
    fake.kmeans.res = met.proj.clustering.by.exp.orig.fake.kmeans,
    plot.data = met.proj.data.clustering.by.exp.orig.plot.data
  ))
}

get.clustering.res.for.subtype <- function(
  subtype,
  surv.plot.data.orig,
  surv.plot.data.integrated,
  survival.data,
  intend.projections.filtered,
  exp.met.data.filtered) {
  
  possible.clusters <- 1:10
  names(possible.clusters) <- possible.clusters
  set.seed(213354734)
  TCGA.PATIENT.LENGTH <- 12
  
  exp.proj.data <- t(intend.projections.filtered$expression)
  met.proj.data <- t(intend.projections.filtered$methylation)
  rownames(exp.proj.data) <- paste(rownames(exp.proj.data), "1", sep = "_")
  rownames(met.proj.data) <- paste(rownames(met.proj.data), "2", sep = "_")
  
  exp.proj.data.umap <-
    surv.plot.data.integrated$reduced.data[surv.plot.data.integrated$Modality == "Gene Expression",]
  met.proj.data.umap <-
    surv.plot.data.integrated$reduced.data[surv.plot.data.integrated$Modality == "DNA Methylation", ]
  
  clustering.results <- (lapply(possible.clusters, function(number.of.clusters) {
    # utils.log(sprintf("Number of clusters: %d", number.of.clusters))
    
    clustered.integrated.data <- get.clustered.integrated.data.for.plot(
      surv.plot.data.integrated,
      number.of.clusters = number.of.clusters,
      cluster.umap = TRUE, #doesn't work without it
      intend.projections = intend.projections.filtered
    )
    
    exp.indices <- which(surv.plot.data.integrated[[BY.MODALITY]] == EXPRESSION.GROUP)
    clustered.projected.exp.data <- get.clustered.integrated.data.for.plot(
      filter.plot.data(surv.plot.data.integrated, exp.indices),
      number.of.clusters = number.of.clusters,
      cluster.umap = TRUE, #doesn't work without it
      intend.projections = list(expression = intend.projections.filtered$expression, methylation = NULL)
    )
    
    met.indices <- which(surv.plot.data.integrated[[BY.MODALITY]] == METHYLATION.GROUP)
    clustered.projected.met.data <- get.clustered.integrated.data.for.plot(
      filter.plot.data(surv.plot.data.integrated, met.indices),
      number.of.clusters = number.of.clusters,
      cluster.umap = TRUE, #doesn't work without it
      intend.projections =  list(expression = NULL, methylation = intend.projections.filtered$methylation)
    )
    
    clustered.origin.data <- get.clustered.origin.data.for.plot(
      surv.plot.data.orig,
      number.of.clusters = number.of.clusters,
      cluster.umap = FALSE,
      exp.met.data = exp.met.data.filtered
    )
    clustered.origin.data.umap <- get.clustered.origin.data.for.plot(
      surv.plot.data.orig,
      number.of.clusters = number.of.clusters,
      cluster.umap = TRUE,
      exp.met.data = exp.met.data.filtered
    )
    
    met.proj.by.exp.orig.results.5.umap.none <- get.met.proj.by.exp.orig.results(
      exp.proj.data,
      met.proj.data, 
      clustered.origin.data,
      met.indices,
      surv.plot.data.integrated,
      neighbors.to.consider = 5)

    met.proj.by.exp.orig.results.5.umap.orig <- get.met.proj.by.exp.orig.results(
      exp.proj.data,
      met.proj.data, 
      clustered.origin.data.umap,
      met.indices,
      surv.plot.data.integrated,
      neighbors.to.consider = 5)
    
    met.proj.by.exp.orig.results.5.umap.proj <- get.met.proj.by.exp.orig.results(
      exp.proj.data.umap,
      met.proj.data.umap, 
      clustered.origin.data,
      met.indices,
      surv.plot.data.integrated,
      neighbors.to.consider = 5)
    
    met.proj.by.exp.orig.results.5.umap.both <- get.met.proj.by.exp.orig.results(
      exp.proj.data.umap,
      met.proj.data.umap, 
      clustered.origin.data.umap,
      met.indices,
      surv.plot.data.integrated,
      neighbors.to.consider = 5)
    
    return(list(
      kmeans = list(kmeans.int = clustered.integrated.data$kmeans.res,
           kmeans.int.proj.exp = clustered.projected.exp.data$kmeans.res,
           kmeans.int.proj.met = clustered.projected.met.data$kmeans.res,
           kmeans.origin.exp = clustered.origin.data$kmeans.exp,
           kmeans.origin.met = clustered.origin.data$kmeans.met,
           kmeans.origin.exp.umap = clustered.origin.data.umap$kmeans.exp,
           kmeans.origin.met.umap = clustered.origin.data.umap$kmeans.met,
           
           met.proj.by.exp.orig.umap.none = met.proj.by.exp.orig.results.5.umap.none$fake.kmeans.res,
           met.proj.by.exp.orig.umap.orig = met.proj.by.exp.orig.results.5.umap.orig$fake.kmeans.res,
           met.proj.by.exp.orig.umap.proj = met.proj.by.exp.orig.results.5.umap.proj$fake.kmeans.res,
           met.proj.by.exp.orig.umap.both = met.proj.by.exp.orig.results.5.umap.both$fake.kmeans.res
      ),
      plot.data = list(plot.data.clustered.int = clustered.integrated.data$plot.data,
           plot.data.clustered.int.proj.exp = clustered.projected.exp.data$plot.data,
           plot.data.clustered.int.proj.met = clustered.projected.met.data$plot.data,
           plot.data.clustered.origin = clustered.origin.data$plot.data,
           plot.data.clustered.origin.umap = clustered.origin.data.umap$plot.data,
           
           plot.data.met.proj.by.exp.orig.umap.none = met.proj.by.exp.orig.results.5.umap.none$plot.data,
           plot.data.met.proj.by.exp.orig.umap.orig = met.proj.by.exp.orig.results.5.umap.orig$plot.data,
           plot.data.met.proj.by.exp.orig.umap.proj = met.proj.by.exp.orig.results.5.umap.proj$plot.data,
           plot.data.met.proj.by.exp.orig.umap.both = met.proj.by.exp.orig.results.5.umap.both$plot.data
      )
    ))
  }))

  clustering.results.kmeans <- lapply(clustering.results, function(res) res$kmeans)
  tot.withinss <- sapply(
    clustering.results.kmeans,
    function(results.for.num.clusters) sapply(
      results.for.num.clusters, 
      function(kmeans.res) kmeans.res$tot.withinss
    )
  )
  
  best.num.of.clusters <- apply(
    tot.withinss,
    1, 
    function(values) 
    {
      second.derivatives <- sapply(
        2:(length(values) - 1),
        function(i) (values[i + 1] + values[i - 1] - 2 * values[i])
      )
      return(which.max(second.derivatives) + 1)
    }
  )
    
  clustering.results.2.and.above <- clustering.results.kmeans[-1]
  pvalues <- sapply(
    seq_along(clustering.results.2.and.above),
    function(clustering.res, num.clusters, i) sapply(
      clustering.res[[i]], 
      function(kmeans.res) 
      {
        surv.vs.clusters <- data.frame(
          surv = survival.data[substr(names(kmeans.res$cluster), 1, TCGA.PATIENT.LENGTH)],
          clusters = kmeans.res$cluster
        )
        pval <- get.p.val(survdiff(surv ~ clusters, data = surv.vs.clusters), as.numeric(num.clusters[[i]]))
        # pval <- get.empirical.surv(
        #   survival.data[substr(names(kmeans.res$cluster), 1, TCGA.PATIENT.LENGTH)],
        #   kmeans.res$cluster
        # )$pvalue
        
        return(pval)
      }
    ),
    clustering.res = clustering.results.2.and.above,
    num.clusters = names(clustering.results.2.and.above)
  )
  colnames(pvalues) <- names(clustering.results.2.and.above)
  
  return(list(
    clustering.results = clustering.results,
    tot.withinss = tot.withinss,
    best.num.of.clusters = best.num.of.clusters,
    pvalues = pvalues
  ))
}

get.exp.met.data.for.cluster.plot <- function(subtype, expect.same.samples.for.exp.met = TRUE)
{
  exp.met.data <- get.expression.methylation.data.for.expression.prediction(subtype)
  
  num.selected.features <- 2000

  exp.met.data$expression <- select.highly.variable.features(exp.met.data$expression, num.selected.features)
  exp.met.data$expression <- scale.data(exp.met.data$expression)

  exp.met.data$methylation <- select.highly.variable.features(exp.met.data$methylation, num.selected.features)
  exp.met.data$methylation <- scale.data(exp.met.data$methylation)
  
  stopifnot(!expect.same.samples.for.exp.met ||
              identical(colnames(exp.met.data$expression), colnames(exp.met.data$methylation)))
  
  colnames(exp.met.data$expression) <- paste(colnames(exp.met.data$expression), "1", sep = "_")
  colnames(exp.met.data$methylation) <- paste(colnames(exp.met.data$methylation), "2", sep = "_")
  
  return(exp.met.data)
}

get.clustered.origin.data.for.plot <- function(
  plot.data,
  number.of.clusters,
  cluster.umap,
  exp.met.data = NULL
)   
{
  max.iterations <- 100
  number.of.random.start.sets <- 100
  
  if (cluster.umap) {
    k.means.exp <- kmeans(
      plot.data$reduced.data[which(plot.data$Modality == EXPRESSION.GROUP), ],
      centers = number.of.clusters, 
      iter.max = max.iterations,
      nstart = number.of.random.start.sets
    )
    
    k.means.met <- kmeans(
      plot.data$reduced.data[which(plot.data$Modality == METHYLATION.GROUP), ],
      centers = number.of.clusters, 
      iter.max = max.iterations,
      nstart = number.of.random.start.sets
    )
  } else {
    k.means.exp <- kmeans(
      t(exp.met.data$expression),
      centers = number.of.clusters, 
      iter.max = max.iterations,
      nstart = number.of.random.start.sets
    )
    k.means.met <- kmeans(
      t(exp.met.data$methylation),
      centers = number.of.clusters, 
      iter.max = max.iterations,
      nstart = number.of.random.start.sets
    )
  }
  
  plot.data$cluster <- rep(NA, length(plot.data$Modality))
  names(plot.data$cluster) <- row.names(plot.data$reduced.data)
  plot.data$cluster[names(k.means.exp$cluster)] <- k.means.exp$cluster
  plot.data$cluster[names(k.means.met$cluster)] <- k.means.met$cluster
  stopifnot(!anyNA(plot.data$cluster))
  
  plot.data$cluster <- as.factor(plot.data$cluster)
  
  return(list(plot.data = plot.data, kmeans.exp = k.means.exp, kmeans.met = k.means.met))
}
  

plot.as.seperate.datasets.2d.by.clusters <- function(
  subtype,
  plot.data,
  number.of.clusters,
  cluster.umap,
  exp.met.data = NULL,
  facets.group.by = BY.MODALITY,
  new.window = FALSE) 
{
  plot.data <- get.clustered.origin.data.for.plot(
    plot.data,
    number.of.clusters,
    cluster.umap,
    exp.met.data
  )$plot.data
  
  return(plot.as.seperate.datasets.2d(
    subtype,
    plot.data,
    facets.group.by = facets.group.by,
    color.by = "cluster",
    new.window = new.window
  ))
}

add.cluster.info.to.plot.data <- function(plot.data, cluster.info) 
{
  plot.data$cluster <- rep(NA, length(plot.data$Modality))
  names(plot.data$cluster) <- row.names(plot.data$reduced.data)
  plot.data$cluster[names(cluster.info)] <- cluster.info
  stopifnot(!anyNA(plot.data$cluster))
  
  plot.data$cluster <- as.factor(plot.data$cluster)
  return(plot.data)
}

get.clustered.integrated.data.for.plot <- function(
  plot.data,
  number.of.clusters,
  cluster.umap,
  intend.projections = NULL
)   
{
  max.iterations <- 100
  number.of.random.start.sets <- 100
  
  if (cluster.umap) {
    k.means.int <- kmeans(
      plot.data$reduced.data,
      centers = number.of.clusters, 
      iter.max = max.iterations,
      nstart = number.of.random.start.sets
    )
  } else {
    if (!is.null(intend.projections$expression)) {
      colnames(intend.projections$expression) <- paste(colnames(intend.projections$expression), "1", sep = "_")
    }
    if (!is.null(intend.projections$methylation)) {
      colnames(intend.projections$methylation) <- paste(colnames(intend.projections$methylation), "2", sep = "_")
    }
    
    k.means.int <- kmeans(
      t(cbind(intend.projections$expression, intend.projections$methylation)),
      centers = number.of.clusters, 
      iter.max = max.iterations,
      nstart = number.of.random.start.sets
    )
  }
  
  plot.data$cluster <- rep(NA, length(plot.data$Modality))
  names(plot.data$cluster) <- row.names(plot.data$reduced.data)
  plot.data$cluster[names(k.means.int$cluster)] <- k.means.int$cluster
  stopifnot(!anyNA(plot.data$cluster))
  
  plot.data$cluster <- as.factor(plot.data$cluster)
  
  return(list(
    plot.data = add.cluster.info.to.plot.data(plot.data, k.means.int$cluster),
    kmeans.res = k.means.int
  ))
}


plot.as.integrated.dataset.2d.by.clusters <- function(
  algorithm,
  plot.data,
  number.of.clusters,
  cluster.umap,
  intend.projections = NULL,
  new.window = FALSE,
  show.legend = FALSE,
  show.title = TRUE
)   
{
  plot.data <- get.clustered.integrated.data.for.plot(
    plot.data,
    number.of.clusters,
    cluster.umap,
    intend.projections
  )$plot.data
    
  return(plot.as.integrated.dataset.2d(
    algorithm,
    plot.data,
    color.by = "cluster",
    new.window = new.window,
    show.legend = show.legend,
    show.title = show.title
  ))
}

plot.tcga.lccs.luad.data <- function()
{
  subtype <- LUAD.TCGA.LCCS.INTEGRATION.SUBTYPE
  dimension <- DIMENSION.FOR.INTEGRATION.COMPARISON
  
  luad.tcga.lccs.umap.data <- run.umap.for.sutype("LUAD")
  luad.tcga.lccs.umap.data$methylation.per.gene <- NULL
  
  luad.lccs.expression <- get.lccs.luad.expression.preprocessed.and.normalized()
  set.seed(123)
  luad.tcga.lccs.umap.data$expression <- umap(t(luad.lccs.expression), metric = "cosine")
  
  original.datasets.for.plot <- create.original.datasets.data.for.plot(
    luad.tcga.lccs.umap.data$expression$layout,
    luad.tcga.lccs.umap.data$methylation$layout,
    "UMAP")
  
  intend.projections <- get.intend.projections(subtype, dimension)
  colnames(intend.projections$methylation) <-  
    gsub(pattern = "[.]", replacement = "-", x = colnames(intend.projections$methylation))
  
  integrated.data.for.plot <- get.integrated.data.for.plot.with.umap(intend.projections)
  
  plot.1 <- plot.as.seperate.datasets.2d(subtype, original.datasets.for.plot)
  plot.2 <- plot.as.integrated.dataset.2d(subtype, integrated.data.for.plot, show.legend = T, show.title = F)
  
  common.legend <- get_legend(plot.2 + theme(legend.direction = "horizontal", legend.justification = "center"))
  plot.2 <- plot.2 + theme(legend.position = "none")
  
  plot.row.1 <- plot_grid(plot.1, plot.2, labels = c('A','B'), nrow = 1, rel_widths = c(0.66, 0.33))
  plot.row.2 <- plot_grid(plotlist = common.legend)
  
  dev.new(height = 400, width = 1200, unit = "px", noRStudioGD = T)
  plot_grid(plot.row.1, plot.row.2, ncol = 1, rel_heights = c(1, 0.1))
}

get.original.datasets.for.plot.with.umap <- function(
  subtype,
  multi.subtypes.attribution.expression = NULL,
  multi.subtypes.attribution.methylation = multi.subtypes.attribution.expression
) 
{
  original.umap.data <- run.umap.for.sutype(subtype)
  
  return(create.original.datasets.data.for.plot(
    original.umap.data$expression$layout,
    original.umap.data$methylation$layout,
    "UMAP",
    multi.subtypes.attribution.expression,
    multi.subtypes.attribution.methylation))
}

get.integrated.data.for.plot.with.umap <- function(
  projections,
  multi.subtypes.attribution.expression = NULL,
  multi.subtypes.attribution.methylation = multi.subtypes.attribution.expression
) 
{
  integrated.results <- get.inetgrated.results.for.plot(
    projections,
    multi.subtypes.attribution.expression,
    multi.subtypes.attribution.methylation)
  
  set.seed(123)
  integrated.umap <- umap(t(integrated.results$integrated.data), metric = "cosine")
  
  integrated.data.for.plot <- list(
    reduced.data = integrated.umap$layout,
    reduction.method = "UMAP",
    sample.num = integrated.results$sample.num)
  integrated.data.for.plot[[BY.MODALITY]] <- integrated.results$modalities
  integrated.data.for.plot[[BY.SUBTYPE]] <- integrated.results$subtype.attribution
  
  return(integrated.data.for.plot)
}

get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(list(legend))
}

plot.original.and.integrated.data.for.subtype <- function(subtype) 
{

  original.datasets.for.plot <- get.original.datasets.for.plot.with.umap(subtype)

  original.plot.list <- c(
    list(plot.as.seperate.datasets.2d(subtype, original.datasets.for.plot)),
    list(plot.as.seperate.datasets.2d.with.labeled.samples(subtype, original.datasets.for.plot)))
  
  projections <- get.projections(subtype)

  plots.of.integrated.data <- sapply(seq_along(projections), function(i) {
    
    integrated.data.for.plot <- get.integrated.data.for.plot.with.umap(projections[[i]])
    
    return(c(
      list(plot.as.integrated.dataset.2d(
        names(projections)[i],
        integrated.data.for.plot,
        show.legend = (i == 1))),
      list(plot.as.integrated.dataset.2d.with.labeled.samples(
        names(projections)[i],
        integrated.data.for.plot))
    ))
  })
  
  common.legend <- get_legend(plots.of.integrated.data[[1]] + 
                                theme(legend.direction = "horizontal", legend.justification = "center"))

  plots.of.integrated.data[[1]] <- plots.of.integrated.data[[1]]  + theme(legend.position = "none")
  
  dev.new(height = 1200, width = 1500, unit = "px", noRStudioGD = T)
  plot.row.1 <- plot_grid(plotlist = original.plot.list, labels = c('A','B'), nrow = 1)
  plot.row.2 <- plot_grid(plotlist = plots.of.integrated.data[c(1,3,5,7)], labels = c('C','D','E','F'), nrow = 1)
  plot.row.3 <- plot_grid(plotlist = plots.of.integrated.data[c(2,4,6,8)], labels = c('G','H','I','J'), nrow = 1)
  plot.row.4 <- plot_grid(plotlist = common.legend)
  
  plot_grid(plot.row.1, plot.row.2, plot.row.3, plot.row.4, ncol = 1, rel_heights = c(1.1, 1.2, 1, 0.15))
}

get.classifiction.confusion.matrices.plot.list <- function(
  integrated.data.for.plot,
  multi.subtypes.attribution
)
{
  plot.list <- lapply(integrated.data.for.plot, function(p)
    {
      NEIGHBORS.TO.CONSIDER <- 5
      
      # reference <- t(p$expression)
      # query <- t(p$methylation)
      reference <- p$reduced.data[p$Modality == "Gene Expression",]
      query <- p$reduced.data[p$Modality == "DNA Methylation", ]
      
      query.class <- knn(
        train = reference,
        test = query, 
        cl = multi.subtypes.attribution, 
        k = NEIGHBORS.TO.CONSIDER
      )
      confusion <- confusionMatrix(
        data = query.class,
        reference = as.factor(multi.subtypes.attribution)
      )
      confusion.matrix <- confusion$table
      confusion.matrix.norm <- as.table(apply(confusion.matrix, 2, function(x) x/sum(x)))
      
      plot <- 
        ggplot(as.data.frame(confusion.matrix.norm), aes(Reference, Prediction, fill =  Freq*100)) +
        geom_tile() + 
        geom_text(aes(label = sprintf("%.02f",Freq*100)), size = 3) +
        scale_fill_gradient(name = "Prediction (% of Reference samples)", low = "lightyellow", high = "darkred")
    }
  )
  return(plot.list)
}

plot.original.and.integrated.data.for.multi.subtypes <- function(subtypes) 
{
  subtype <- get.multi.subtypes.name(subtypes)
  
  multi.subtypes.attribution.file.path <- sprintf(
    PREPROCESSED.EXPRESSION.METHYLATION.DATA.FORMAT,
    SUBTYPE.ATTR.PREFIX,
    subtype
  )
  multi.subtypes.attribution <- readRDS(multi.subtypes.attribution.file.path)
  
  original.datasets.for.plot <- get.original.datasets.for.plot.with.umap(subtype, multi.subtypes.attribution)
  
  original.plot.list <- c(
    list(plot.as.seperate.datasets.2d(subtype, original.datasets.for.plot)),
    list(plot.as.seperate.datasets.2d(subtype, original.datasets.for.plot, color.by = BY.SUBTYPE)))
  
  projections <- get.projections(subtype)
  integrated.data.for.plot <- lapply(
    projections, 
    function(p) get.integrated.data.for.plot.with.umap(p, multi.subtypes.attribution)
  )
  
  plots.of.integrated.data <- sapply(seq_along(projections), function(i) {
    return(c(
      list(plot.as.integrated.dataset.2d(
        names(projections)[i],
        integrated.data.for.plot[[i]],
        show.legend = (i == 1))),
      list(plot.as.integrated.dataset.2d(
        names(projections)[i],
        integrated.data.for.plot[[i]],
        color.by = BY.SUBTYPE,
        show.legend = (i == 1),
        show.title = FALSE))
    ))
  })
  
  common.legend.by.omic <- get_legend(plots.of.integrated.data[[1]] + 
                                theme(legend.direction = "horizontal", legend.justification = "center"))
  common.legend.by.cancer.type <- get_legend(plots.of.integrated.data[[2]] + 
                                        theme(legend.direction = "horizontal", legend.justification = "center"))
  
  plots.of.integrated.data[[1]] <- plots.of.integrated.data[[1]]  + theme(legend.position = "none")
  plots.of.integrated.data[[2]] <- plots.of.integrated.data[[2]]  + theme(legend.position = "none")
  
  plots.of.confusion.matrices <- get.classifiction.confusion.matrices.plot.list(
    integrated.data.for.plot,
    multi.subtypes.attribution
  )
  
  common.legend.for.confusion.matrices <- get_legend(
    plots.of.confusion.matrices[[1]] + 
      theme(legend.direction = "horizontal", legend.justification = "center"))
  
  plots.of.confusion.matrices <- lapply(plots.of.confusion.matrices, function(plot) 
    {
      return(plot + theme(legend.position = "none"))
    }
  )
  
  dev.new(height = 1500, width = 1500, unit = "px", noRStudioGD = T)
  plot.row.1 <- plot_grid(plotlist = c(common.legend.by.omic, common.legend.by.cancer.type))
  plot.row.2 <- plot_grid(plotlist = original.plot.list, labels = c('A','B'), nrow = 1)
  plot.row.3 <- plot_grid(plotlist = plots.of.integrated.data[c(1,3,5,7)], labels = c('C','D','E','F'), nrow = 1)
  plot.row.4 <- plot_grid(plotlist = plots.of.integrated.data[c(2,4,6,8)], labels = c('G','H','I','J'), nrow = 1)
  plot.row.5 <- plot_grid(plotlist = plots.of.confusion.matrices, labels = c("K", "L", "M", "N"), nrow = 1)
  plot.row.6 <- plot_grid(plotlist = common.legend.for.confusion.matrices)
  
  plot_grid(
    plot.row.1,
    plot.row.2, 
    plot.row.3, 
    plot.row.4, 
    plot.row.5,
    plot.row.6,
    ncol = 1,
    rel_heights = c(0.15, 1.1, 1.2, 1, 1.1, 0.25)
  )
}

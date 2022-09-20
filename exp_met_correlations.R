library(batchelor)
library(dplyr)
library(tibble)
library(ggplot2)
library(cowplot)
library(parallel)
source('utilities.R')
source('constants.R')
source('gene_mapping_utils.R')
source('data_utils.R')
source('intend.R')

# Use case of LUAD samples: gene expression from LCCS and DNA methylation from TCGA

get.luad.integrated.data <- function()
{
  integrated.data <- get.intend.projections(
    subtype = LUAD.TCGA.LCCS.INTEGRATION.SUBTYPE, 
    dimension = DIMENSION.FOR.INTEGRATION.COMPARISON
  )
  return(integrated.data)
}

get.luad.original.data <- function()
{
  expression.methylation <- get.expression.methylation.data.for.expression.prediction("LUAD")
  original.data <- list(
    expression = get.lccs.luad.expression.preprocessed.and.normalized(),
    methylation = expression.methylation$methylation
  )

  return(original.data)
}

get.gene.cpg.sites.mappings <- function()
{
  MARGIN.1.MB <- 1000000
  common.features <- get.common.features.for.all.subtypes()
  
  gene.cpg.sites.mapping.1.mb <- get.filtered.gene.to.cg.sites.mapping.all.subtypes(
    common.features$expression,
    common.features$methylation,
    MARGIN.1.MB,
    MARGIN.1.MB)
  
  mapping.1.mb <- gene.cpg.sites.mapping.1.mb$cg_sites
  names(mapping.1.mb) <- gene.cpg.sites.mapping.1.mb$gene
  
  return(mapping.1.mb)
}

get.mutual.neighbors <- function(integrated.data, original.data)
{
  NEIGHBORS.TO.CONSIDER <- 5
  mutual.neigbors <- findMutualNN(
    t(integrated.data$expression),
    t(integrated.data$methylation),
    k1 = NEIGHBORS.TO.CONSIDER
  )
  
  # Verify integrated expression columns
  stopifnot(identical(colnames(original.data$expression), colnames(integrated.data$expression)))
  # Verify integrated methylation columns
  stopifnot(identical(
    colnames(original.data$methylation),
    gsub(pattern = "[.]", replacement = "-", x = colnames(integrated.data$methylation))))

  return(mutual.neigbors)
}

get.expression.methylation.correlations.for.gene <- function(
  relevant.cpg.sites,
  gene.symbol,
  mutual.neighbors, 
  original.data)
{
  correlations <- lapply(
    relevant.cpg.sites,
    function(site) cor.test(original.data$expression[gene.symbol, mutual.neighbors$first],
                       original.data$methylation[site, mutual.neighbors$second])[c("p.value", "estimate")]
  )
  names(correlations) <- relevant.cpg.sites
  return(correlations)
}

get.metylation.sites.info <- function()
{
  methylation.manifest <- read.csv(file = HUMAN.METHYLATION.450.MANIFEST.PATH)
  methylation.sites.info <- methylation.manifest %>% 
    column_to_rownames("IlmnID") %>%
    filter(is.valid.value(CHR) & is.valid.value(MAPINFO) & Genome_Build == 37) %>%
    select(chr = CHR, location = MAPINFO)

  return(methylation.sites.info)
}

compute.exp.met.cor <- function()
{
  integrated.data <- get.luad.integrated.data()
  original.data <- get.luad.original.data()
  original.data.tcga <- get.expression.methylation.data.for.expression.prediction("LUAD")

  methylation.sites.info <- get.metylation.sites.info()
  mappings <- get.gene.cpg.sites.mappings()
  
  genes <- intersect(rownames(original.data$expression), rownames(original.data.tcga$expression))
  genes <- intersect(genes, names(mappings))
  mappings <- mappings[genes]
  mutual.neighbors <- get.mutual.neighbors(
    integrated.data = integrated.data,
    original.data = original.data
  )

  gc()
  
  exp.met.cor <- mclapply(
    genes,
    function(gene.symbol)
    {
      sites <- mappings[[gene.symbol]]

      if (length(sites) < 3)
      {
        return(NA)
      }
      correlations <- get.expression.methylation.correlations.for.gene(
        relevant.cpg.sites = sites,
        gene.symbol = gene.symbol, 
        mutual.neighbors = mutual.neighbors,
        original.data = original.data)
      
      correlation.tcga <- lapply(
        sites,
        function(site) cor.test(original.data.tcga$expression[gene.symbol, ],
                                original.data.tcga$methylation[site, ])[c("p.value", "estimate")]
      )
      names(correlation.tcga) <- sites
      
      
      correlated.sites <- subset(methylation.sites.info, row.names(methylation.sites.info) %in% names(correlations))
      correlated.sites$correlation.pval <- sapply(correlations[rownames(correlated.sites)], function(x) x$p.value)
      correlated.sites$correlation.estimate <- sapply(correlations[rownames(correlated.sites)], function(x) x$estimate)

      correlated.sites$correlation.tcga.pval <- sapply(correlation.tcga[rownames(correlated.sites)], function(x) x$p.value)
      correlated.sites$correlation.tcga.estimate <- sapply(correlation.tcga[rownames(correlated.sites)], function(x) x$estimate)
      
      return(correlated.sites)
    }
  )
  
  names(exp.met.cor) <- genes
  
  saveRDS(exp.met.cor, file = "CachedData/exp_met_cor.rds")
  return(exp.met.cor)
}

get.plot.correlation.analysis.of.tk1 <- function(exp.met.cor) {
  
  correlations.per.gene <- exp.met.cor$TK1
  
  R.squared <- cor(
    correlations.per.gene$correlation.estimate,
    correlations.per.gene$correlation.tcga.estimate) ^ 2
  
  eq <- substitute(~~italic(R)^2~"="~r2, list(r2 = format(R.squared, digits = 3)))
  
  pval.threshold <- 0.00001
  
  correlation.vs.correlation.plot <- ggplot(
    data = correlations.per.gene,
    mapping = aes(x = correlation.tcga.estimate, y = correlation.estimate)) +
    geom_point(show.legend = FALSE) +
    geom_point(
      data = subset(correlations.per.gene, correlations.per.gene$correlation.pval >= pval.threshold |
                      correlations.per.gene$correlation.tcga.pval >= pval.threshold),
      mapping = aes(x = correlation.tcga.estimate, y = correlation.estimate), col = "powderblue") +
    geom_point(
      data = subset(correlations.per.gene, correlations.per.gene$correlation.pval < pval.threshold &
                      correlations.per.gene$correlation.tcga.pval < pval.threshold),
      mapping = aes(x = correlation.tcga.estimate, y = correlation.estimate), col = "royalblue4") +
    xlab("Correlations based on TCGA multi-omic data") +
    ylab("Correlations based on INTEND projections") +
    theme_cowplot(font_size = 12) +
    background_grid(major = "xy", minor = "xy")
  
  TK1.start <- 76170160
  TK1.end <- 76183314
  margins <- 100000
  minor.breaks.res <- 10000
  major.breaks.res <- 50000
  
  start <- TK1.start - margins
  end <- TK1.end + margins
  
  # Extracted from GeneHancer - regions with 4 or more GH sources
  top.enhancers <- data.frame(
    GH17J078120 = c(start = 76116641, end = 76131420),
    GH17J078183 = c(start = 76179781, end = 76184882),
    GH17J078250 = c(start = 76246401, end = 76251482),
    GH17J078255 = c(start = 76251641, end = 76255910),
    GH17J078136 = c(start = 76132687, end = 76138938),
    GH17J078410 = c(start = 76406730, end = 76414888),
    GH17J078170 = c(start = 76167076, end = 76173481),
    GH17J078112 = c(start = 76107683, end = 76112280)
  )
  # Transpose
  top.enhancers <- as.data.frame(t(as.matrix(top.enhancers)))
  
  round_any <- function(x, accuracy, f=round){f(x / accuracy) * accuracy}
  
  correlation.vs.genomic.coordinates.plot <- 
    ggplot(data = correlations.per.gene) +
    coord_cartesian(xlim = c(start, end)) +
    scale_x_continuous(breaks = round_any(seq(start, end, by = major.breaks.res),major.breaks.res),
                       minor_breaks = round_any(seq(start, end, by = minor.breaks.res),minor.breaks.res)) +
    geom_point(
      data = subset(correlations.per.gene, correlations.per.gene$correlation.pval >= pval.threshold),
      mapping = aes(x = location, y = correlation.estimate), col = "powderblue") +
    geom_point(
      data = subset(correlations.per.gene, correlations.per.gene$correlation.pval < pval.threshold),
      mapping = aes(x = location, y = correlation.estimate), col = "royalblue4") +
    geom_vline(xintercept = TK1.start) +
    geom_vline(xintercept = TK1.end) +
    xlab("Genomic location on chromosome 17") +
    ylab("Correlation coefficient") +
    geom_rect(
      data = top.enhancers,
      inherit.aes = FALSE,
      aes(
        xmin = start,
        xmax = end,
        ymin = min(min(correlations.per.gene$correlation.estimate),min(correlations.per.gene$correlation.estimate)),
        ymax = max(max(correlations.per.gene$correlation.estimate),max(correlations.per.gene$correlation.estimate))),
      color = "transparent",
      fill = "burlywood1",
      alpha = 0.3) +
    theme_cowplot(font_size = 12) +
    background_grid(major = "xy", minor = "xy")
  
  return(
    list(
      correlation.vs.correlation.plot = correlation.vs.correlation.plot,
      correlation.vs.genomic.coordinates.plot = correlation.vs.genomic.coordinates.plot
    )
  )
}

analyze.correlation.significance.and.sign.consensus <- function(correlations.flattened) {
  pval.threshold <- 0.01
  indices.of.sites.with.significant.correlation.integration <-
    which(correlations.flattened$intend.integration.pval < pval.threshold)
  indices.of.sites.with.significant.correlation.tcga <-
    which(correlations.flattened$tcga.multi.omic.pval < pval.threshold)
  
  indices.of.sites.with.significant.correlation.both <- intersect(
    indices.of.sites.with.significant.correlation.integration,
    indices.of.sites.with.significant.correlation.tcga
  )
  significant.correlations <- correlations.flattened[indices.of.sites.with.significant.correlation.both,]
  fraction.both.out.of.all.sites <- nrow(significant.correlations) / nrow(correlations.flattened)
  number.of.sites.with.sign.consensus <- 
    sum(sign(significant.correlations$intend.integration) == sign(significant.correlations$tcga.multi.omic))
  fraction.of.sign.consensus <- number.of.sites.with.sign.consensus / nrow(significant.correlations)
  
  utils.log(sprintf("Percentage of sites with significant correlation in both methods: %.02f",
                    fraction.both.out.of.all.sites * 100))
  utils.log(sprintf("Percentage of sites with correlation in both methods (out of the significant sites): %.02f",
                    fraction.of.sign.consensus * 100))
}

analyze.exp.met.cor <- function() {
  
  exp.met.cor <- readRDS("CachedData/exp_met_cor.rds")
  
  common.features <- get.common.features.for.all.subtypes()
  
  mapping <- get.filtered.gene.to.cg.sites.mapping.all.subtypes(
    common.features$expression,
    common.features$methylation,
    UPSTREAM.MARGIN,
    DOWNSTREAM.MARGIN)
  
  exp.met.cor <- exp.met.cor[names(exp.met.cor) %in% mapping$gene]
  exp.met.cor <- exp.met.cor[!is.na(exp.met.cor)]
  
  correlations <- lapply(
    seq_along(exp.met.cor),
    function(i)  
      data.frame(
        intend.integration = exp.met.cor[[i]]$correlation.estimate,
        intend.integration.pval = exp.met.cor[[i]]$correlation.pval,
        tcga.multi.omic = exp.met.cor[[i]]$correlation.tcga.estimate,
        tcga.multi.omic.pval = exp.met.cor[[i]]$correlation.tcga.pval,
        pair.used.in.model = rownames(exp.met.cor[[i]]) %in%
          mapping$cg_sites[mapping$gene == names(exp.met.cor)[i]][[1]]
      )
  )
  names(correlations) <- names(exp.met.cor)
  
  correlations.flattened <- bind_rows(correlations)
  analyze.correlation.significance.and.sign.consensus(correlations.flattened)
  
  correlations.zeroed <- lapply(
    correlations,
    function(df) {
      df$intend.integration[is.na(df$intend.integration)] <- 0
      df$tcga.multi.omic[is.na(df$tcga.multi.omic)] <- 0
      return(df)
    }
  )

  cor.test.of.correlations.per.gene <- lapply(
    correlations.zeroed, 
    function(cor.per.gene) cor.test(
      cor.per.gene$intend.integration,
      cor.per.gene$tcga.multi.omic,
      alternative = "greater")
  )
  p.val.per.gene <- sapply(cor.test.of.correlations.per.gene, function(x) x$p.value)
  cor.of.cors.per.gene <- sapply(cor.test.of.correlations.per.gene, function(x) x$estimate)

  
  correlations.flattened.filtered <- correlations.flattened %>% 
    filter((!is.na(intend.integration) & (!is.na(tcga.multi.omic))))
  
  observed.vs.expected.correlations.plot <- 
    ggplot(data = correlations.flattened.filtered, aes(
      x = tcga.multi.omic, 
      y = intend.integration, 
      color = pair.used.in.model)) +
    geom_point(shape = ".", show.legend = FALSE) +
    theme_cowplot() +
    xlab("Correlations based on TCGA multi-omic data") +
    ylab("Correlations based on INTEND projections") +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      text = element_text(size = 12),
      axis.text = element_text(size = 12),
    ) +
    scale_color_brewer(palette = "Dark2")
  
  R.squared.for.all.correlations <- cor(
    correlations.flattened.filtered$tcga.multi.omic,
    correlations.flattened.filtered$intend.integration) ^ 2

  correlations.between.observed.and.expected.per.gene.histogram = ggplot() +
      geom_histogram(data = data.frame(cor.of.cors.per.gene),
                     mapping = aes(x = cor.of.cors.per.gene),
                     color = "grey20",
                     fill = "grey95",
                     binwidth = 0.1) +
      xlim(c(-1,1)) +
      xlab("Correlation between estimated and
obtained gene-site correlations") +
      ylab("Number of genes") +
      theme_cowplot(font_size = 12)
    
    tk1.plots <- get.plot.correlation.analysis.of.tk1(exp.met.cor)
    
    plot.c <- plot_grid(observed.vs.expected.correlations.plot, NULL, ncol = 1, rel_heights = c(97.4, 2.6))
    plot.d <- plot_grid(correlations.between.observed.and.expected.per.gene.histogram)
    plot.e <- plot_grid(tk1.plots$correlation.vs.correlation.plot, NULL, ncol = 1, rel_heights = c(97.4, 2.6))
    plot.f <- plot_grid(tk1.plots$correlation.vs.genomic.coordinates.plot)
    
    dev.new(height = 400, width = 400, unit = "px", noRStudioGD = T)
    plot_grid(plot.c, plot.d, plot.e, plot.f, nrow = 2, ncol = 2, labels = c("C","D","E","F"))
}





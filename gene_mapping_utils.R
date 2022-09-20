library(rtracklayer)
library(dplyr)
library(parallel)
source('utilities.R')

HG19.REF.GENE.PATH <- "Data/hg19/hg19.refGene.gtf"
GENE.COORDINATES.PATH <- "CachedData/hg19/gene.coordinates.rda"

HUMAN.METHYLATION.450.MANIFEST.PATH <-
  "Data/Methylation/Illumina_HumanMethylation450_BeadChip/HumanMethylation450_15017482_v1-2.csv"

GENE.TO.CG.SITE.MAPPING.PATH.FORMAT <- "CachedData/gene_to_cg_site_mappings/mapping_%d_%d.rda"

get.full.ref.gene.df <- function() 
{
  hg19.ref.gene <- rtracklayer::import(HG19.REF.GENE.PATH)
  hg19.ref.gene.df <- as.data.frame(hg19.ref.gene)
  return(hg19.ref.gene.df)
}

get.gene.coordinates.df <- function(force.update = FALSE) 
{
  if (file.exists(GENE.COORDINATES.PATH) & !force.update) {
    return(readRDS(file = GENE.COORDINATES.PATH))
  }
  
  hg19.ref.gene.df <- get.full.ref.gene.df()
  
  first.chr.index <- which(levels(hg19.ref.gene.df$seqnames) == "chr1")
  last.chr.index <- which(levels(hg19.ref.gene.df$seqnames) == "chrY")
  
  gene.coordinates.df <- hg19.ref.gene.df %>% 
    filter(type == "transcript" &
             as.integer(seqnames) >= first.chr.index & as.integer(seqnames) <= last.chr.index) %>%
    select(gene_id, chr = seqnames, strand, start, end) %>% 
    distinct()
  
  saveRDS(gene.coordinates.df, file = GENE.COORDINATES.PATH)
  return(gene.coordinates.df)
}

create.gene.to.cg.site.mapping <- function(
  gene.coordinates.df, upstream.margin.in.bases, downstream.margin.in.bases, force.update = FALSE)
{
  mapping.file.path <- 
    sprintf(GENE.TO.CG.SITE.MAPPING.PATH.FORMAT, upstream.margin.in.bases, downstream.margin.in.bases)
  
  if (file.exists(mapping.file.path) & !force.update) {
    return(readRDS(file = mapping.file.path))
  }
  
  methylation.manifest <- read.csv(file = HUMAN.METHYLATION.450.MANIFEST.PATH)
  
  methylation.manifest.filtered <- methylation.manifest %>%
    filter(is.valid.value(CHR) & is.valid.value(MAPINFO) & Genome_Build == 37)
  
  gene.to.cg.site.mapping <- data.frame(gene = unique(gene.coordinates.df$gene_id))
  gene.to.cg.site.mapping$cg_sites <- lapply(1:nrow(gene.to.cg.site.mapping), function(x) c())
  
  cg.sites.per.gene <- mclapply(1:nrow(gene.to.cg.site.mapping), function(gene_idx) {
    gene <- gene.to.cg.site.mapping$gene[gene_idx]
    gene.coordinates.filtered <- gene.coordinates.df %>% filter(gene_id == gene)
    
    #convert chromosome to string comparable with the one in methylation.manifest.filtered$CHR
    gene.coordinates.filtered$chr <- sapply(gene.coordinates.filtered$chr, function(chr.str) 
      substr(as.character(chr.str), start = nchar("chr") + 1, stop = nchar(as.character(chr.str))))
    
    gene.cg.sites <- c()
    
    for (i in 1:nrow(gene.coordinates.filtered)) {
      chr <- gene.coordinates.filtered$chr[i]
      
      if (gene.coordinates.filtered$strand[i] == "+") {
        
        start <- gene.coordinates.filtered$start[i] - upstream.margin.in.bases
        end <- gene.coordinates.filtered$end[i] + downstream.margin.in.bases
        
      } else {
        
        start <- gene.coordinates.filtered$start[i] - downstream.margin.in.bases
        end <- gene.coordinates.filtered$end[i] + upstream.margin.in.bases
      }
      
      cg.sites.in.range <- (methylation.manifest.filtered %>% 
                              filter(tolower(CHR) == tolower(chr) & MAPINFO >= start & MAPINFO <= end) %>%
                              select(Name))$Name
      
      gene.cg.sites <- c(gene.cg.sites, cg.sites.in.range)
    }
    
    if (gene_idx %% 1000 == 1) {
      utils.log(sprintf("Done mapping gene number %d", gene_idx - 1))
    }
    
    return(gene.cg.sites)
  })
  
  gene.to.cg.site.mapping$cg_sites <- lapply(cg.sites.per.gene, unique)
  
  saveRDS(gene.to.cg.site.mapping, file = mapping.file.path)
  return(gene.to.cg.site.mapping)
}

get.filtered.gene.to.cg.sites.mapping.for.exp.met.features <- function(
  expression.features,
  methylation.features,
  upstream.margin.in.bases,
  downstream.margin.in.bases,
  force.update = FALSE)
{
  gene.to.cg.site.mapping <- create.gene.to.cg.site.mapping(
    gene.coordinates.df = get.gene.coordinates.df(force.update = force.update),
    upstream.margin.in.bases = upstream.margin.in.bases,
    downstream.margin.in.bases =  downstream.margin.in.bases,
    force.update = force.update)
  
  # Ignore genes which are not on in the input gene expression data
  mapping <- gene.to.cg.site.mapping %>% 
    filter(gene %in% expression.features)
  
  # Ignore CpG sites which are not on in the input methylation data
  mapping$cg_sites <- mclapply(mapping$cg_sites, function(cpg.sites.for.gene) 
    intersect(cpg.sites.for.gene, methylation.features), mc.cores = LIMITED.NUMBER.OF.CORES)
  
  # Remove rows with < 2 CpG sites mapped to the gene
  cpg.sites.lengths <- sapply(mapping$cg_sites, length)
  mapping <- mapping[which(cpg.sites.lengths > 1), ]
  
  return(mapping)
}

get.filtered.gene.to.cg.sites.mapping.all.subtypes <- function(
  expression.features,
  methylation.features,
  upstream.margin.in.bases,
  downstream.margin.in.bases,
  force.update = FALSE)
{
  GENE.TO.CG.SITE.MAPPING.ALL.SUBTYPES.PATH.FORMAT <-
    "CachedData/gene_to_cg_site_mappings/mapping_all_%d_%d.rda"
  
  mapping.file.path <- 
    sprintf(GENE.TO.CG.SITE.MAPPING.ALL.SUBTYPES.PATH.FORMAT, upstream.margin.in.bases, downstream.margin.in.bases)
  
  if (file.exists(mapping.file.path) & !force.update) {
    return(readRDS(file = mapping.file.path))
  }
  
  mapping <- get.filtered.gene.to.cg.sites.mapping.for.exp.met.features(
    expression.features,
    methylation.features,
    upstream.margin.in.bases,
    downstream.margin.in.bases,
    force.update
  )
  
  saveRDS(mapping, file = mapping.file.path)
  return(mapping)
}

library(survival)
library(dplyr)
library(parallel)

get.survival.data <- function(tcga.clinical.data)
{
  # Filter patients with missing both death time and censor time 
  missing.value.strings <- c("[Not Available]", "[Not Applicable]")
  tcga.clinical.data.filtered <- tcga.clinical.data %>% 
    filter(!(tcga.clinical.data$last_contact_days_to %in% missing.value.strings &
             tcga.clinical.data$death_days_to %in% missing.value.strings))
  
  # Filter patients with both days_to_death <= 0 and last_contact_days_to <= 0
  tcga.clinical.data.filtered <- tcga.clinical.data.filtered %>% 
    filter(!(tcga.clinical.data.filtered$last_contact_days_to <= 0 &
               tcga.clinical.data.filtered$death_days_to <= 0))
  
  names <- rownames(tcga.clinical.data.filtered)
  time <- as.integer(ifelse(tcga.clinical.data.filtered$last_contact_days_to == "[Not Available]",
                                     tcga.clinical.data.filtered$death_days_to,
                                     tcga.clinical.data.filtered$last_contact_days_to))
  
  # Supress warnings for this line:
  old.warn <- getOption("warn")
  options(warn = -1)
  
  event <- ifelse(!is.na(as.integer(tcga.clinical.data.filtered$death_days_to)), 1, 0)
  
  # Restore warning state
  options(warn = old.warn)
  
  surv <- Surv(time, event)
  names(surv) <- names

  
  return(surv)
}

get.short.tcga.samples.names <- function(tcga.samples.names)
{
  TCGA.SAMPLE.SHORT.LENGTH <- 12
  TCGA.SAMPLE.OPTIONAL.SUFFIX.START <- 16
  
  return(paste0(
    substring(tcga.samples.names, first = 1, last = TCGA.SAMPLE.SHORT.LENGTH),
    substring(tcga.samples.names, first = TCGA.SAMPLE.OPTIONAL.SUFFIX.START)
  ))
}

filter.only.tumor.samples <- function(tcga.samples.names)
{
  SAMPLE.TYPE.CODE.SUFFIX.START <- 14
  SAMPLE.TYPE.CODE.SUFFIX.END <- 15
  PRIMARY.TUMOR.OR.METASTATIC.SAMPLE.TYPE.CODES <- c("01", "03", "06", "09")
  
  sample.type.codes <- substring(
    tcga.samples.names,
    SAMPLE.TYPE.CODE.SUFFIX.START,
    SAMPLE.TYPE.CODE.SUFFIX.END)
  
  samples <- tcga.samples.names[sample.type.codes %in% PRIMARY.TUMOR.OR.METASTATIC.SAMPLE.TYPE.CODES]
  samples <- sort(samples)
  
  i <- 1
  while (i < (length(samples) - 1)) {
    TCGA.PATIENT.LENGTH <- 12
    if (substr(samples[i], 1, TCGA.PATIENT.LENGTH) == substr(samples[i + 1], 1, TCGA.PATIENT.LENGTH)) {
      samples <- samples[-i]
    } else {
      i <- i + 1
    }
  }
  
  return(samples)
}

check.survival <- function(survival.data, clusters) {
  df <- data.frame(surv = survival.data, clust = clusters)
  return(survdiff(surv ~ clust, data = df))
}

get.p.val <- function(survdiff.res, number.of.clusters)
{
  return(pchisq(survdiff.res$chisq, df = number.of.clusters - 1, lower.tail = FALSE))
}

get.logrank.pvalue <- function(survdiff.res) {
  1 - pchisq(survdiff.res$chisq, length(survdiff.res$n) - 1)  
}

get.empirical.surv <- function(survival.data, clusters) {
  set.seed(42)
  surv.ret = check.survival(survival.data, clusters)
  orig.chisq = surv.ret$chisq
  orig.pvalue = get.logrank.pvalue(surv.ret)
  # The initial number of permutations to run
  num.perms = round(min(max(10 / orig.pvalue, 1e4), 1e6))
  should.continue = T
  
  total.num.perms = 0
  total.num.extreme.chisq = 0
  
  while (should.continue) {
    print('Another iteration in empirical survival calculation')
    print(num.perms)
    perm.chisq = as.numeric(mclapply(1:num.perms, function(i) {
      new.indices <- sample(1:length(clusters))
      cur.clusters = clusters[new.indices]
      cur.chisq = check.survival(survival.data, cur.clusters)$chisq
      return(cur.chisq)
    }, mc.cores = 1))
    
    total.num.perms = total.num.perms + num.perms
    total.num.extreme.chisq = total.num.extreme.chisq + sum(perm.chisq >= orig.chisq)
    
    binom.ret = binom.test(total.num.extreme.chisq, total.num.perms)
    cur.pvalue = binom.ret$estimate
    cur.conf.int = binom.ret$conf.int
    
    print(c(total.num.extreme.chisq, total.num.perms))
    print(cur.pvalue)
    print(cur.conf.int)
    
    sig.threshold = 0.05
    is.conf.small = ((cur.conf.int[2] - cur.pvalue) < min(cur.pvalue / 10, 0.01)) & ((cur.pvalue - cur.conf.int[1]) < min(cur.pvalue / 10, 0.01))
    is.threshold.in.conf = cur.conf.int[1] < sig.threshold & cur.conf.int[2] > sig.threshold
    if ((is.conf.small & !is.threshold.in.conf) | (total.num.perms > 2e7)) {
      should.continue = F
    } else {
      num.perms = 1e5
    }
  }
  
  return(list(pvalue = cur.pvalue, conf.int = cur.conf.int, total.num.perms = total.num.perms, 
              total.num.extreme.chisq = total.num.extreme.chisq))
}
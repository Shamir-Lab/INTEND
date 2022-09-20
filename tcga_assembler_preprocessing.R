preprocessed_tcga_assembler_data_dir <- "Data/TCGA_assembler_preprocessed"

get.omic.data <- function(subtype, file.name) {
  subtype.dataset.dir = file.path(preprocessed_tcga_assembler_data_dir, subtype)
  cur.files = list.files(subtype.dataset.dir)
  if (!any(grepl(file.name, cur.files))) {
    error('no file found')
  } else {
    # Loads Data and Des
    load(file.path(subtype.dataset.dir, paste0(file.name, '.rda')))
    remove.feat = rowSums(is.na(Data)) >= ncol(Data) * 5 / 100
    Data = Data[!remove.feat, ]
    Des = Des[!remove.feat,,drop = F]
    remove.pat = colSums(is.na(Data)) >= nrow(Data) * 5 / 100
    Data = Data[, !remove.pat]
    Data = do.call(rbind, lapply(1:nrow(Data), function(i) {
      cur.row = Data[i,]
      cur.row[is.na(cur.row)] = mean(cur.row, na.rm = T)
      return(cur.row)
    }))
    return(list(Data = Data, Des = Des))
  }
}

get.rna.data <- function(subtype) {return(get.omic.data(subtype, 'processed_rna'))}
get.methy.data <- function(subtype) {return(get.omic.data(subtype, 'processed_methy'))}

# We only treat the project-TSS-paticipant-sample part of the barcode and ignore the
# vial-portion-analyte-plate-center part.
tcga.ids.to.sample <- function(tcga.ids) {
  SAMPLE.PREFIX.LENGTH <- 15
  substr(tcga.ids, 1, SAMPLE.PREFIX.LENGTH)
}

get.clinical.data <- function(subtype)
{
  subtype.dataset.dir <- file.path(preprocessed_tcga_assembler_data_dir, subtype)
  clinical.data.file.path <- file.path(subtype.dataset.dir, 'clinical.txt')
  
  clinical.raw.csv <- read.csv(file = clinical.data.file.path, sep = "\t", header = T)
  first.relevant.row <- 3
  first.relevant.column <- 3
  
  rownames.column <- "bcr_patient_barcode"
  
  clinical.data <- clinical.raw.csv[first.relevant.row:nrow(clinical.raw.csv),
                                      first.relevant.column:ncol(clinical.raw.csv)]
  
  rownames(clinical.data) <- clinical.raw.csv[first.relevant.row:nrow(clinical.raw.csv), rownames.column]
  
  return(clinical.data)
}
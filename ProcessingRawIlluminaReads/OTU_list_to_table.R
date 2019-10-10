#--File to which OTU data is added
raw.fe <- read.csv(paste0(dat.dir,'CultureBased_only/20190206_CB_Taxonomy.csv'),
                   as.is = T)

raw.fe[c('otu.90','otu.95','otu.97','otu.99','esv')] <- NA

#--90% sequence similarity
otu.90 <- read.csv(paste0(data.out,'CultureBased_only/MobyleSnapWB/20190205_',
                          'CB_90%sequencesimilarity.csv'),
                   as.is = T)

for(i in colnames(otu.90)){
  sample <- as.matrix(otu.90[i])
  for (s in sample) {
    raw.fe[raw.fe$sequence_name == s, 'otu.90'] <- colnames(sample)
  }
}

#--95% sequence similarity
otu.95 <- read.csv('data/CultureFree_CultureBased/ITS2_CBandCFfromfiltered_95%.csv',
                   as.is = T)

for(i in colnames(otu.95)){
  sample <- as.matrix(otu.95[i])
  for (s in sample[!sample == ""]) {
    raw.fe[raw.fe$sequence_name == s, 'otu.95'] <- colnames(sample)
  }
}

#--97% sequence similarity
otu.97 <- read.csv(paste0(data.out,'CultureBased_only/MobyleSnapWB/20190205_',
                          'CB_97%sequencesimilarity.csv'),
                   as.is = T)

for(i in colnames(otu.97)){
  sample <- as.matrix(otu.97[i])
  for (s in sample[sample != '',]) {
    raw.fe[raw.fe$sequence_name == s, 'otu.97'] <- colnames(sample)
  }
}

#--99% sequence similarity
otu.99 <- read.csv(paste0(data.out,'CultureBased_only/MobyleSnapWB/20190205_',
                          'CB_99%sequencesimilarity.csv'),
                   as.is = T)

for(i in colnames(otu.99)){
  sample <- as.matrix(otu.99[i])
  for (s in sample) {
    raw.fe[raw.fe$sequence_name == s, 'otu.99'] <- colnames(sample)
  }
}


#--Exact sequence variance
esv <- read.csv(paste0(data.out,'CultureBased_only/MobyleSnapWB/20190205_',
                       'CB_ExactSequenceVariance.csv'),
                as.is = T)

for(i in colnames(otu.99)){
  sample <- as.matrix(otu.99[i])
  for (s in sample) {
    raw.fe[raw.fe$sequence_name == s, 'esv'] <- colnames(sample)
  }
}

write.csv(raw.fe, paste0(data.out,'/CultureBased_only/20190206_CB_SeqSimTaxonomy.csv'),
          row.names = F)

no.otu.designation <- data.frame(raw.fe[is.na(raw.fe$esv),'sequence_name'])

write.csv(no.otu.designation,
          paste0(data.out,'/CultureBased_only/20190206_NoOtuDesignation.csv'),
          row.names = F)


#install.packages('seqinr')
library(seqinr)

cb.cf <- read.fasta(paste0('~/../../Volumes/Cenococcum/Illumina/',
                  'All_SIS_CBandCF_comparison/filtered/cb_ITSx_out.ITS1.fasta'),
                  as.string = T, forceDNAtolower = F)
sample.data <- read.csv('data/CultureBased_only/20190905_CB_Taxonomy.csv')

sample.names <- unique(sample.data$tree)

#--Isolate SCM data
sample.data.scm <- sample.data[sample.data$range == 'SantaCatalina',]
sample.names.scm <- sample.data[sample.data$range == 'SantaCatalina', 'sequence_name']

#--Change names of scm samples from sample numbers to site numbers
for(i in 1:nrow(sample.data.scm)){
  for(t in sample.names.scm){
    if(sample.data.scm[i, 'sequence_name'] == t){
      tree.t <- sample.data.scm[i, 'tree']
      names(cb.cf) <- gsub(t, paste0(tree.t,'z'), names(cb.cf))} 
  }
}

#--Isolate all data except scm
sample.names.other <- unique(sample.data$tree)
sample.names.other <- sample.names.other[1:86]

#--Replace sequence names for all samples except SCM
for(t in sample.names.other){
    names(cb.cf) <- gsub(t, paste0(t,'z'), names(cb.cf))
}

cb.cf.names <- names(cb.cf)

write.fasta(cb.cf, cb.cf.names,paste0('~/../../Volumes/Cenococcum/Illumina/',
            'All_SIS_CBandCF_comparison/filtered/cb_ITSx_out.ITS1_cleaned.fasta'))

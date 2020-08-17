## The influnce of varying purity Estimates:
## First, I will look into CCF estimations:

# First I will start and see how purity estimates influence CCF estimations:
# For this we will use FacetsSuite: ccf_Annotate_MAF

maf = vroom::vroom('../../CPNA_analysis/TCGA/RawData/mc3.v0.2.8.PUBLIC.maf', delim = '\t')
maf = as.data.frame(maf[, c('Hugo_Symbol', 'Tumor_Sample_Barcode', 
                            'Chromosome', 'Start_Position', 
                            'End_Position', 't_ref_count', 't_alt_count', 
                            'n_depth', 'n_ref_count', 'n_alt_count')])

maf$Sample = substr(maf$Tumor_Sample_Barcode, start = 1, stop = 16)
files = list.files('~/Desktop/tmp_TCGA/tmp_TCGA/', full.names = T)
samples_to_keep = substr(files, start = 78, stop = 93)

# subset maf ~ easier working in R
maf_keep = maf[maf$Sample %in% samples_to_keep, ]
rm(maf)

ccf_all_out = data.frame()
for(i in 1:length(files)){
  print(i)
  
  try({
    sample = substr(files[i], start = 78, nchar(files)[i]-56)
    
    data.in = facets::readSnpMatrix(filename = files[i])
    data.processed = facetsSuite::run_facets(data.in, 
                                        cval = 100,
                                        genome = 'hg19')
    data.segs = data.processed$segs
    data.purity = data.processed$purity
    data.purity = ifelse(is.na(data.purity), 0, data.purity)
    
    # maf
    maf_ccf = maf_keep[grep(sample, maf_keep$Tumor_Sample_Barcode), ]
    
    CCF_tmp = facetsSuite::ccf_annotate_maf(maf = maf_ccf,
                                            segs = data.segs,
                                            purity = data.purity)
    
    ccf_all_out = rbind(ccf_all_out, CCF_tmp)
    
    rm(data.in)
    rm(data.processed)
    rm(data.segs)
    rm(data.purity)
    rm(CCF_tmp)
    rm(maf_ccf)
    rm(sample)
    
  })
}

# write.table(ccf_all_out, file = '~/Desktop/ccf_facets_out.txt', sep = '\t', row.names = F, quote = F)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# investigate the CCF output on n = 435 TCGA-PRAD samples;

# how many patient are monoclonal:
ccf_facets = read.csv('~/Desktop/ccf_facets_out.txt', sep = '\t')

ccf_facets_out = data.frame()
for(i in unique(ccf_facets$Tumor_Sample_Barcode)){
  sample.patient = ccf_facets[which(ccf_facets$Tumor_Sample_Barcode == i), ]
  table.out = as.data.frame(table(sample.patient$clonality))
  summary.data = data.frame(sample = i,
                            n.clonal = ifelse(identical(table.out$Freq[which(table.out$Var1 == 'CLONAL')],
                                                        integer(0)), 0, table.out$Freq[which(table.out$Var1 == 'CLONAL')]),
                            n.subclonal = ifelse(identical(table.out$Freq[which(table.out$Var1 == 'SUBCLONAL')],
                                                           integer(0)), 0, table.out$Freq[which(table.out$Var1 == 'SUBCLONAL')]),
                            n.inter = ifelse(identical(table.out$Freq[which(table.out$Var1 == 'INDETERMINATE')],
                                                       integer(0)), 0, table.out$Freq[which(table.out$Var1 == 'INDETERMINATE')]))
  ccf_facets_out = rbind(ccf_facets_out, summary.data)
  rm(sample.patient)
  rm(summary.data)
}

# assign monoclonal VS polyclonal
ccf_facets_out$diversity = ifelse(ccf_facets_out$n.clonal / 
                                 (ccf_facets_out$n.clonal + ccf_facets_out$n.subclonal + ccf_facets_out$n.inter) >= 0.80,
                               'monoclonal', 'polyclonal')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# look into ABSOLUTE calls and estimate 
Absolute.calls = read.csv('~/Documents/MSKCC/CPNA_analysis/TCGA/RawData/TCGA_ABSOLUTE_CPN.txt', sep = '\t')
Absolute.calls = Absolute.calls[Absolute.calls$Sample %in% substr(ccf_facets$Tumor_Sample_Barcode, start = 1, stop = 15), ]

maf <- read.csv('~/Documents/MSKCC/CPNA_analysis/TCGA/RawData/mc3.v0.2.8.PUBLIC.maf', sep = '\t')
maf = as.data.frame(maf[, c('Hugo_Symbol', 'Tumor_Sample_Barcode', 
                            'Chromosome', 'Start_Position', 
                            'End_Position', 't_ref_count', 't_alt_count', 
                            'n_depth', 'n_ref_count', 'n_alt_count')])

maf.absolute = maf[maf$Tumor_Sample_Barcode %in% ccf_facets$Tumor_Sample_Barcode, ]
maf.absolute$Tumor_Sample_Barcode = substr(maf.absolute$Tumor_Sample_Barcode, start = 1, stop = 15)
colnames(maf.absolute)[2] = 'Sample'

TCGA.purity = TCGAbiolinks::Tumor.purity
TCGA.purity = TCGA.purity[substr(TCGA.purity$Sample.ID, start = 1, stop = 15) %in% maf.absolute$Sample, ]
TCGA.purity$Sample.ID = substr(TCGA.purity$Sample.ID, start = 1, stop = 15)
TCGA.purity$CPE = as.numeric(as.character(gsub(pattern = ',', '.', TCGA.purity$CPE)))
TCGA.purity = TCGA.purity[!is.na(TCGA.purity$CPE), ]


tcga.samples = Reduce(intersect, list(Absolute.calls$Sample, TCGA.purity$Sample.ID, maf.absolute$Sample))
tcga.samples = intersect(tcga.samples, Absolute.calls$Sample)
tcga.samples = unique(tcga.samples)

ccf_absolute_out = data.frame()
for(i in 1:length(tcga.samples)){
  try({
    print(i)
    maf.selected = maf.absolute[which(maf.absolute$Sample == tcga.samples[i]), ]
    maf.selected$Chromosome = as.character(maf.selected$Chromosome)
    data.table::setDT(maf.selected, key = c('Chromosome', 'Start_Position', 'End_Position'))
    
    cpn.selected = Absolute.calls[which(Absolute.calls$Sample == tcga.samples[i]), ]
    cpn.selected$Chromosome = as.character(cpn.selected$Chromosome)
    data.table::setDT(cpn.selected, key = c('Chromosome', 'Start', 'End'))
    
    # data merged:
    data.merged = data.table::foverlaps(maf.selected,
                                        cpn.selected,
                                        by.x = c('Chromosome', 'Start_Position', 'End_Position'),
                                        by.y = c('Chromosome', 'Start', 'End'),
                                        type = 'within',
                                        mult = 'first',
                                        nomatch = NA)
    
    purity = TCGA.purity[which(TCGA.purity$Sample.ID == tcga.samples[i]), 'CPE']
    
    # CCF estimation:
    #data.selected = data.merged[,c(1, 2, 9, 21, 23, 24, 25, 26, 27, 28, 29)]
    data.selected = as.data.frame(data.merged)
    
    for(mutation in 1:nrow(data.selected)){
      
      mutant_copies = expected_mutant_copies(t_var_freq = data.selected$t_alt_count[mutation] / (data.selected$t_ref_count[mutation] + data.selected$t_alt_count[mutation]), 
                                             total_copies = data.selected$Modal_Total_CN[mutation],
                                             purity = purity)
      
      out = estimate_ccf(purity = purity,
                         total_copies = data.selected$Modal_Total_CN[mutation],
                         mutant_copies = mutant_copies,
                         t_alt_count = data.selected$t_alt_count[mutation],
                         t_depth = data.selected$t_alt_count[mutation] + data.selected$t_ref_count[mutation])
      
      ccf_mutant = out[[1]]
      ccf_upper = out[[3]]
      clonality = ifelse(ccf_mutant >= 0.8 | (ccf_mutant >= 0.7 & ccf_upper >= 0.9), 'Clonal',
                         ifelse(ccf_mutant < 0.8 & ccf_mutant > 0.001 | (ccf_mutant < 0.7 | ccf_mutant > 0.001 & ccf_upper < 0.9 | ccf_upper > 0.001), 
                                'Subclonal', 'NA'))
        
      report.out = data.frame(sample = tcga.samples[i],
                              Hugo_Symbol = data.selected$Hugo_Symbol[mutation],
                              chrom = data.selected$Chromosome[mutation],
                              Start = data.selected$Start_Position[mutation],
                              purity = purity,
                              totalCPN = data.selected$Modal_Total_CN[mutation],
                              mutant_copies = mutant_copies,
                              CCF = ccf_mutant,
                              cnA1 = data.selected$Cancer_cell_frac_a1[mutation],
                              cnA2 = data.selected$Cancer_cell_frac_a2[mutation],
                              clonality = clonality)
      
      
      ccf_absolute_out = rbind(ccf_absolute_out, report.out)
      
    }
    
  })
  
  rm(maf.selected)
  rm(cpn.selected)
  rm(data.merged)
  rm(data.selected)
  rm(out)
  rm(ccf_mutant)
  rm(ccf_upper)
  rm(report.out)
  rm(clonality)
  rm(mutant_copies)
  rm(purity)
  
}

# write.table(ccf_absolute_out, file = '~/Desktop/ccf_absolute_out.txt', sep = '\t')

ccf_absolute_out = read.csv('~/Desktop/ccf_absolute_out.txt', sep = '\t')

## diversity estimate on ABSOLUTE:
ccf_out = data.frame()
for(i in unique(ccf_absolute_out$sample)){
  sample.patient = ccf_absolute_out[which(ccf_absolute_out$sample == i), ]
  table.out = as.data.frame(table(sample.patient$clonality))
  summary.data = data.frame(sample = i,
                            n.Clonal = ifelse(identical(table.out$Freq[which(table.out$Var1 == 'Clonal')],
                                                        integer(0)), 0, table.out$Freq[which(table.out$Var1 == 'Clonal')]),
                            n.Subclonal = ifelse(identical(table.out$Freq[which(table.out$Var1 == 'Subclonal')],
                                                           integer(0)), 0, table.out$Freq[which(table.out$Var1 == 'Subclonal')]),
                            n.NA = length(sample.patient$clonality[which(is.na(sample.patient$clonality))]))
  ccf_out = rbind(ccf_out, summary.data)
  rm(sample.patient)
  rm(summary.data)
}

ccf_out$diversity = ifelse(ccf_out$n.Clonal / 
                                 (ccf_out$n.Clonal + ccf_out$n.Subclonal + ccf_out$n.NA) >= 0.80,
                               'monoclonal', 'polyclonal')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# comparison between tools;
ccf_out$category = rep('TCGA', nrow(ccf_out))
ccf_facets_out$category = rep('FACETS', nrow(ccf_facets_out))
colnames(ccf_facets_out) = c('sample', 'n.Clonal', 'n.Subclonal', 'n.NA', 'diversity', 'category')
ccf_merged = rbind(ccf_out, ccf_facets_out)

library(ggplot2)
library(ggbeeswarm)

# clonal gene mutations
comparison = ggplot(subset(ccf_merged, n.Clonal <= 150), 
                    aes(x = category, y = n.Clonal, color = category)) +
  geom_quasirandom(dodge.width = 0.8,
                   shape = 16,
                   size = 1.5,
                   alpha = 0.9) +
  theme_std(base_size = 14) +
  theme(aspect.ratio = 1, legend.position = 'none') +
  scale_y_continuous(expand = expansion(mult = 0.01), 
                     limits = c(0, 150),
                     breaks = seq(0, 150, 50)) +
  scale_color_manual(values = c('mediumaquamarine', 'springgreen4')) +
  labs(x = '', y = 'predicted clonally mutated genes')

comparison

# subclonal
comparison.subclonal = ggplot(subset(ccf_merged, n.Subclonal < 300), 
                    aes(x = category, y = n.Subclonal, color = category)) +
  geom_quasirandom(dodge.width = 0.8,
                   shape = 16,
                   size = 1.5,
                   alpha = 0.9) +
  theme_std(base_size = 14) +
  theme(aspect.ratio = 1, legend.position = 'none') +
  scale_y_continuous(expand = expansion(mult = 0.01), 
                     limits = c(0, 150),
                     breaks = seq(0, 150, 50)) +
  scale_color_manual(values = c('mediumaquamarine', 'springgreen4')) +
  labs(x = '', y = 'predicted subclonally mutated genes')

comparison.subclonal

library(egg)
ggarrange(comparison, comparison.subclonal, nrow = 1)





# make a plot to compare the CCF density of certain genes;
# now I compare CCF obtained from ABSOLUTE/TP purity with Facets data
# look at several genes

plot_list = list()
for(i in c('TP53', 'FOXA1', 'SPOP', 'KMT2D', 'KMT2C', 'ATM')){
  
  data_facets = ccf_facets[which(ccf_facets$Hugo_Symbol == i), ]
  data_facets$Sample = NULL
  data_facets$sample = substr(data_facets$Tumor_Sample_Barcode, start = 1, stop = 15)
  data_facets = data_facets[!duplicated(data_facets$Tumor_Sample_Barcode), ]
  
  # absolute data
  data_absolute = ccf_absolute_out[which(ccf_absolute_out$Hugo_Symbol == i), ]
  data_absolute = data_absolute[!duplicated(data_absolute), ]
  
  data_merged = merge(data_absolute, data_facets[, c('sample', 'purity', 'ccf_expected_copies')],
                      by = 'sample', all.x = T)
  
  data_merged$class = ifelse(data_merged$CCF > data_merged$ccf_expected_copies, 'red', 'green')
  
  # make plot:
  # left_label = paste(FOXA_merged$sample)
  
  ccf_comparison_plot = ggplot(data_merged) +
    geom_segment(aes(x = 1, xend = 2, y = CCF, yend = ccf_expected_copies, col = class), 
                 size = 0.75, show.legend = F) +
    
    geom_vline(xintercept = 1, linetype = 'dashed', size = 0.1) +
    geom_vline(xintercept = 2, linetype = 'dashed', size = 0.1) +
    
    scale_color_manual(labels = c('Up', 'Down'),
                       values = c("green" = '#392C8D', "red" = '#EAB9B8')) +
    
    labs(x = '', y = 'cancer cell fraction', title = i) +
    
    theme_std(base_size = 14) +
    
    theme(panel.background = element_blank(), 
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_blank(),
          aspect.ratio = 1,
          plot.margin = unit(c(1,1,1,1), "lines")) +
    
    xlim(0.5, 2.5) + 
    
    scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, 0.2), limits = c(0, 1.1)) +
    
    geom_text(label="TCGA", x = 1, y = 1.01, hjust = 1.1, size = 5) +
    geom_text(label="Facets", x = 2.2, y = 1.01, hjust = 0.05, size = 5)
  
  plot_list[[i]] = ccf_comparison_plot

}

# arrange all plots in one picture
egg::ggarrange(plots = plot_list, ncol = 3, nrow = 2)


# png()
# plot(1,
#      type = 'n',
#      xlim = c(0, 1.3),
#      ylim = c(0, 2.8),
#      axes = F,
#      xlab = 'CCF estimation',
#      ylab = '',
#      pch = 16,
#      main = '',
#      xaxt = 'n',
#      xaxs = 'i',
#      yaxs = 'i')
# 
# axis(side = 2, at = c(0, 6), labels = FALSE, tck = 0.0, lwd = 2)
# mtext(side = 2, text = "Density (a.u.)", line = 1)
# 
# axis(side = 1,
#      at = c(0, 0.5, 1),
#      labels = c('0', '0.5', '1'),
#      las = 1,
#      lwd = 2)
# 
# # add Facets density of TP53
# lines(density(x$ccf_expected_copies[which(x$Hugo_Symbol == 'SPOP')], na.rm = T), lwd = 3, col = '#392C8D')
# text(x = 1.15, y = 2.4, labels = 'Facets\n estimation', col = '#392C8D')
# lines(density(ccf_absolute_out$CCF[which(ccf_absolute_out$Hugo_Symbol == 'SPOP')], na.rm = T), lwd = 3, col = '#EAB9B8')
# text(x = 0.3, y = 1.1, labels = 'Absolute\n estimation', col = '#EAB9B8')
# 
# title(main = 'SPOP1 (n=61)', line = -1, adj = 0.1, cex.main = 1.5)
# 
# dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# helper functions to estimate CCF

estimate_ccf = function(purity,
                        total_copies,
                        mutant_copies,
                        t_alt_count,
                        t_depth) {
  
  
  ccfs = seq(0.001, 1, 0.001)
  expected_vaf  = function(ccf, purity, total_copies) {
    purity * ccf * mutant_copies / (2 * (1 - purity) + purity * total_copies)
  }
  
  probs = sapply(ccfs, function(c) {
    stats::dbinom(t_alt_count, t_depth, expected_vaf(c, purity, total_copies))
  })
  probs = probs / sum(probs)
  
  ccf_max = which.max(probs)
  if (identical(ccf_max, integer(0))) ccf_max = NA
  ccf_half_max = which(probs > max(probs) / 2)
  ccf_lower = max(ccf_half_max[1] - 1, 1) # closest ccf value before half-max range (within 0-1 range)
  ccf_upper = min(ccf_half_max[length(ccf_half_max)] + 1, length(ccfs)) # closest ccf value after half-max range (within 0-1 range)
  if (is.na(purity)) ccf.upper = NA 
  ccf_max = ccf_max / length(ccfs)
  ccf_lower = ccf_lower / length(ccfs)
  ccf_upper = ccf_upper / length(ccfs)
  prob95 = sum(probs[950:1000])
  prob90 = sum(probs[900:1000])
  
  list(ccf_max, ccf_lower, ccf_upper, prob95, prob90)
}

# Estimate mutant copy number, given observed VAF, purity, and local ploidy
# Based on PMID 28270531
expected_mutant_copies = function(t_var_freq,
                                  total_copies,
                                  purity) {
  
  if (is.na(total_copies)) {
    NA_real_
  } else {
    if (total_copies == 0) total_copies = 1
    mu = t_var_freq * (1 / purity) * (purity * total_copies + (1 - purity) * 2)
    alt_copies = ifelse(mu < 1, 1, abs(mu)) # mu < 1 ~ 1, mu >= 1 ~ abs(mu)
    round(alt_copies)
  }
}








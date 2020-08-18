# Orthogonal validation of purity estimates
# theoretical model: see power point slide

# load facets data:
data.in = read.csv('~/Desktop/ccf_facets_out.txt', sep = '\t')
head(data.in)
View(data.in)

test = data.in[which(data.in$Tumor_Sample_Barcode == 'TCGA-VP-A87D-01A-11D-A34U-08'), ]
test = test[which(test$t_depth >= 50), ]
test = test[which(test$t_var_freq >= 0.10), ]
test = test[!test$tcn == 0, ]
test = test[complete.cases(test$Hugo_Symbol), ]

library(Ryacas)
library(Ryacas0)

merged = data.frame()
for(i in 1:nrow(test)){
  expected_VAF = (test$expected_alt_copies[i] * test$purity[i]) / (2*(1 - test$purity[i]) + test$tcn[i]*test$purity[i])
  formula.syntax = paste0(test$t_var_freq[i], ' == (', test$expected_alt_copies[i], '*x) / (2*(1 - x) + ', test$tcn[i], '*x)')
  expr = yacas(formula.syntax)
  expected_purity = Solve(expr, 'x')
  expected_purity_solved = Eval(expected_purity)
  purity_vector = round(as.numeric(noquote(substr(expected_purity_solved, start = 8, stop = 12))), 2)
  
  out = data.frame(sample = test$Tumor_Sample_Barcode[i],
                   gene = test$Hugo_Symbol[i],
                   expected_VAF = expected_VAF,
                   observed_VAF = test$t_var_freq[i],
                   expected_purity = purity_vector,
                   observed_purity = test$purity[i],
                   total_copies = test$tcn[i],
                   mut_copies = test$expected_alt_copies[i])
  
  merged = rbind(merged, out)
  rm(expected_VAF)
  rm(formula.syntax)
  rm(expr)
  rm(expected_purity)
  rm(expected_purity_solved)
  rm(purity_vector)
}
merged

p = function(VAF, purity, mut, total){
  return((mut*purity) / (2*(1-purity) + total*purity))
}

p(purity = seq(0, 1, 0.1), mut = 2, total = 2)

# make a plot for heterozygous mutations

plot(1,
     type = 'n',
     xlim = c(0, 1.0),
     ylim = c(0, 1),
     axes = F,
     xlab = '<VAF>',
     ylab = '<purity>',
     pch = 16,
     main = '',
     xaxt = 'n',
     xaxs = 'i',
     yaxs = 'i')

axis(side = 2, at = c(0, 0.5, 1), labels = c(0, 0.5, 1), tck = 0.0, lwd = 2, las = 2)
#mtext(side = 2, text = "<purity>", line = 1)

axis(side = 1,
     at = c(0, 0.5, 1),
     labels = c('0', '0.5', '1'),
     las = 1,
     lwd = 2)

# theory
lines(x = p(purity = seq(0, 1, 0.1), mut = 1, total = 2), 
      y = seq(0, 1, 0.1), lwd = 3, col = '#392C8D')
text(x = 0.65, y = 0.95, labels = expression(paste('Heterozygous: <VAF> ~ ', frac(p, 2))), col = 'black', cex = 1.2)
title(main = 'one arbitrary sample', adj = 0, )

# observed
points(x = merged$observed_VAF[which(merged$total_copies == 2 & merged$mut_copies == 1)], 
       y = merged$observed_purity[which(merged$total_copies == 2 & merged$mut_copies == 1)], pch = 19, col = '#91b0a1')
text(x = 0.15, y = 0.75, labels = 'Facets global purity estimate')

# expected
points(x = merged$observed_VAF[which(merged$total_copies == 2 & merged$mut_copies == 1)], 
       y = merged$expected_purity[which(merged$total_copies == 2 & merged$mut_copies == 1)], 
       pch = 19, col = 'red')

abline(h = mean(merged$expected_purity[which(merged$total_copies == 2 & merged$mut_copies == 1)]), lwd = 0.5, lty = 'dashed', col = 'black')



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## LOH

plot(1,
type = 'n',
xlim = c(0, 1.0),
ylim = c(0, 1),
axes = F,
xlab = '<VAF>',
ylab = '<purity>',
pch = 16,
main = '',
xaxt = 'n',
xaxs = 'i',
yaxs = 'i')

axis(side = 2, at = c(0, 0.5, 1), labels = c(0, 0.5, 1), tck = 0.0, lwd = 2, las = 2)
#mtext(side = 2, text = "<purity>", line = 1)

axis(side = 1,
     at = c(0, 0.5, 1),
     labels = c('0', '0.5', '1'),
     las = 1,
     lwd = 2)

# theory
lines(x = p(purity = seq(0, 1, 0.1), mut = 1, total = 1), 
      y = seq(0, 1, 0.1), lwd = 3, col = '#392C8D')
text(x = 0.65, y = 0.95, labels = expression(paste('LOH: <VAF> ~ ', frac(p, 2-p))), col = 'black', cex = 1.2)
title(main = 'one arbitrary sample', adj = 0, )

# observed
points(x = merged$observed_VAF[which(merged$total_copies == 1 & merged$mut_copies == 1)], 
       y = merged$observed_purity[which(merged$total_copies == 1 & merged$mut_copies == 1)], pch = 19, col = '#91b0a1')
text(x = 0.15, y = 0.75, labels = 'Facets global purity estimate')

# expected
points(x = merged$observed_VAF[which(merged$total_copies == 1 & merged$mut_copies == 1)], 
       y = merged$expected_purity[which(merged$total_copies == 1 & merged$mut_copies == 1)], 
       pch = 19, col = 'red')

abline(h = mean(merged$expected_purity[which(merged$total_copies == 1 & merged$mut_copies == 1)]), lwd = 0.5, lty = 'dashed', col = 'black')





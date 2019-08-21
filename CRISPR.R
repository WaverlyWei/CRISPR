library(here)
library(MAGeCKFlute)

countsummary <- read.delim("CRISPR_Count.countsummary.txt")
gene <- read.table('T0puro_ContD20/T0puro_ContD20_rra.gene_summary.txt',header=T)
sgrna <- read.delim('T0puro_ContD20/T0puro_ContD20_rra.sgrna_summary.txt',header=T)

# generate plots in ???./RRA_Flute_Results/???
FluteRRA(gene, sgrna, prefix="RRA", organism="hsa")

# plot count summary 
head(countsummary)
MapRatesView(countsummary)

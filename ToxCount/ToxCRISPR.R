library(here)
library(MAGeCKFlute)

countsummary <- read.delim(here("Desktop/SuperFund/CRISPR/Count/CRISPR_Count.countsummary.txt"))
gene <- read.table(here('Desktop/SuperFund/CRISPR/Count/T0puro_ContD20/T0puro_ContD20_rra.gene_summary.txt'),header=T)
sgrna <- read.delim(here('Desktop/SuperFund/CRISPR/Count/T0puro_ContD20/T0puro_ContD20_rra.sgrna_summary.txt'),header=T)

# generate plots in ???./RRA_Flute_Results/???
FluteRRA(gene, sgrna, prefix="Tox_100FA_D20_RRA", organism="hsa")

# plot count summary 
head(countsummary)
MapRatesView(countsummary)


# TRY: get gene list 
dd.rra <- ReadRRA(gene)
head(dd.rra)
dd.sgrna <- ReadsgRRA(sgrna)
head(dd.sgrna)

geneList <- dd.rra$LFC
names(geneList) <- dd.rra$Official


write.csv(dd.rra, file = here('Desktop/SuperFund/CRISPR/Count/Genomewide_GeneList/T0puro_ContD20_geneList.csv') )

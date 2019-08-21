Sweave("CRISPR_Count_countsummary.Rnw");
library(tools);

texi2dvi("CRISPR_Count_countsummary.tex",pdf=TRUE);


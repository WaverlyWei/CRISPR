library(here)
library(limma)
library(tidyverse)
library(edgeR)

count <- read.delim(here("SuperFund/CRISPR/Count","CRISPR_Count.count.txt"))
D8_40 <- count %>% select(X40.D8.1_S4_L001,
                             X40.D8.2_S5_L001,
                             X40.D8.3_S6_L001,
                            Cont.D8.1_S1_L001,
                            Cont.D8.2_S2_L001,
                            Cont.D8.4_S3_L001)
                         
#design <- model.matrix(~X100.D20.1_S18_L002,data = D20_100)
design <- cbind(Intercept=1,Group=c(1,1,0,0))
fit <- lmFit(D20_100FA,design = design)
res <- decideTests(fit)

genes <- toptable(fit, n = nrow(count))
idx <- as.numeric(rownames(genes))
# incorporate gene names 
genes$gene <- as.character(count$Gene[idx])
# filter significant genes
geneList <- genes %>% filter(adj.P.Val < 0.05)
dim(geneList)

write.csv(geneList,
          file = here("Desktop/SuperFund/CRISPR/ToxCount/Tox_Limma/D20_100FA.csv"
                      ))



# ========= edge R ========= #
group <- factor(c(1,1,0,0))
design <- model.matrix(~group)
y <- DGEList(D20_100FA, group=group)
y_d <- estimateDisp(y, design)
#y_lrt <- glmFit(y)
# y_exact <- exactTest(y_d)
# res <- decideTestsDGE(y_exact, adjust.method="BH", p.value=0.05, lfc=0)
# 
# predlfc<-predFC(y,design,prior.count=1)
# logfc <- predFC(y,design,prior.count=0)
# 
# p.value <- exactTestDoubleTail(D8_0.5AA[,1:2], D20_0.5AA[,3:4])

fit <- glmQLFit(y_d, design)
qlf <- glmQLFTest(fit)
res <- topTags(qlf, n = dim(count)[1])
output <- res$table


# output <- as.data.frame(cbind(as.character(count$Gene),
#                               logfc[,2],
#                               p.value))
                        
# names(output) <- c("gene","LFC","p.value")
# output$gene <- as.character(output$gene)
# output$LFC <- as.numeric(as.character(output$LFC))
# output$p.value <- as.numeric(as.character(output$p.value))

head(output)
idx <- as.numeric(rownames(output))
# incorporate gene names 
output$gene <- as.character(count$Gene[idx])

write.csv(output,
          file = here("Desktop/SuperFund/CRISPR/ToxCount/Tox_edgeR/D20_100FA.csv"
          ))


# ========== edge R normalization ========== #
group <- c(1,1,1,2,2,2)
y <- DGEList(counts=D8_40, group=group)
y$samples
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# model 
group <- factor(c(1,1,1,0,0,0))
design <- model.matrix(~group)
y <- DGEList(y, group=group)
y_d <- estimateDisp(y, design)


fit <- glmQLFit(y_d, design)
qlf <- glmQLFTest(fit)
res <- topTags(qlf, n = dim(count)[1])
output <- res$table


head(output)
idx <- as.numeric(rownames(output))
# incorporate gene names 
output$gene <- as.character(count$Gene[idx])

write.csv(output,
          file = here("SuperFund/CRISPR/ToxCount/Tox_edgeR/D8_40_postpreprocessing.csv"
          ))

# R version 4.2.2 (2022-10-31) -- "Innocent and Trusting"
# R Core Team (2022). R: A language and environment for statistical computing. R Foundation for Statistical
# Computing, Vienna, Austria. URL https://www.R-project.org/.
# Copyright (C) 2022 The R Foundation for Statistical Computing
# Platform: x86_64-apple-darwin17.0 (64-bit)

# Transvar: HG19, UCSC
# Input genome	hg19
# open-cravat	2.3.0
# Gene mapper	UCSC hg38 Gene Mapper (1.10.4)

library(dplyr)
library(ggplot2)
library(ggpubr)

#### results from cravat ####
cravat_output <-
  read.table(
    'cravat/cravat_output_hg19_ucsc.tsv',
    header = T,
    skip = 6,
    sep = '\t'
  )
cravat.headers = colnames(cravat_output)
# expand chasmplus and vest4
cravat.headers.select = cravat.headers[c(1:12, 15, 20, 21,
                                         26, 27, 45:49, 54:64, 66, 67)]

new.headers = c(
  cravat.headers[1:12],
  'CGC_class',
  'CHASMplus_pval',
  'CHASMplus_score',
  'CScape_score',
  'DANN_score',
  'MA_rank_score',
  'MA_functional_impact',
  'MT_rank_score',
  'MT_prediction',
  'MP_rank_score',
  'PhD_SNPg_prediction',
  'PhD_SNPg_score',
  'PP2_hdiv_prediction',
  'PP2_hvar_prediction',
  'PP2_hdiv_rank_score',
  'PP2_hvar_rank_score',
  'PROVEAN_rank_score',
  'PROVEAN_prediction',
  'SIFT_prediction',
  'SIFT_rank_score',
  'Siphy_rank_score',
  'VEST4_score',
  'VEST4_pval'
)
cravat.output.select <- cravat_output[, cravat.headers.select]
colnames(cravat.output.select) <- new.headers
CGI.aa <-
  full_join(CGI, transvar_aa_map, by = c('Alteration' = 'input'))
CGI.cravat = full_join(CGI.aa, cravat.output.select, by = c("id" = "Samples"))

# driver/passenger prediction
# CHASMplus_score: probability of being driver for missense mutations
# VEST4_pval: p-value of being pathogenic for all nonsilent consequence types
CGI.cravat.chasmvest = CGI.cravat %>% # filter(Oncogenic.Summary=="") %>%
  select(Oncogenic.Summary,
         CHASMplus_score,
         CScape_score,
         DANN_score,
         VEST4_pval)

p1 = ggplot(CGI.cravat.chasmvest,
            aes(x = Oncogenic.Summary, y = CHASMplus_score)) + geom_boxplot()
p2 = ggplot(CGI.cravat.chasmvest,
            aes(x = Oncogenic.Summary, y = CScape_score)) + geom_boxplot()
p3 = ggplot(CGI.cravat.chasmvest, aes(x = Oncogenic.Summary, y = DANN_score)) + geom_boxplot()
p4 = ggplot(CGI.cravat.chasmvest, aes(x = Oncogenic.Summary, y = VEST4_pval)) + geom_boxplot()

ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
# chasm, cscape, and vest 4 are better than dann; vest4 is the best

# prediction
CGI.cravat$VEST4_pred <-
  ifelse(CGI.cravat$VEST4_pval < 0.05, 'driver (predicted)', 'passenger')
CGI.cravat$CHASMplus_pred <-
  ifelse(CGI.cravat$CHASMplus_pval < 0.05,
         'driver (predicted)',
         'passenger')
CGI.cravat$CScape_pred <-
  ifelse(CGI.cravat$CScape_score > 0.5,
         'driver (predicted)',
         'passenger')

CGI.cravat.result = CGI.cravat %>% select(
  id,
  Alteration,
  p.aa1infer,
  p.check,
  p.aa3infer,
  Oncogenic.Summary,
  Chrom,
  Position,
  Ref_Base,
  Alt_Base,
  Gene,
  Transcript,
  Sequence_Ontology,
  c.infer,
  cDNA_change,
  Protein_Change,
  VEST4_pred,
  CHASMplus_pred,
  CScape_pred
)

# QC on protein change position
ppos1 <- strsplit(CGI.cravat.result$Alteration, ":")
ppos1 <-
  sapply(strsplit(CGI.cravat.result$Alteration, ":"), tail, 1)
ppos1 <- as.numeric(gsub(".*?([0-9]+).*", "\\1", ppos1))

ppos2 <-
  as.numeric(gsub(".*?([0-9]+).*", "\\1", CGI.cravat.result$Protein_Change))
CGI.cravat.result$ppos.check = ppos1 == ppos2

# remove unmatched multiple mapping with at least one matched
map.count <-
  CGI.cravat.result %>% group_by(id) %>% count() %>% filter(n > 1)
map.true.count <-
  CGI.cravat.result %>% filter(ppos.check) %>% group_by(id) %>% count()
map.count.tofilter <- intersect(map.count$id, map.true.count$id)

CGI.cravat.result$ppos.filter <-
  CGI.cravat$id %in% map.count.tofilter
CGI.cravat.result$ppos.filter <-
  CGI.cravat.result$ppos.filter & CGI.cravat.result$ppos.check == F
CGI.cravat.result2 <- CGI.cravat.result %>% filter(ppos.filter == F)

CGI.cravat.result2$pcheck_transvar_cravat <-
  CGI.cravat.result2$p.aa3infer == CGI.cravat.result2$Protein_Change
CGI.cravat.result2 %>% group_by(pcheck_transvar_cravat) %>% count()
CGI.cravat.result2$ccheck_transvar_cravat <-
  CGI.cravat.result2$c.infer == CGI.cravat.result2$cDNA_change
CGI.cravat.result2 %>% group_by(ccheck_transvar_cravat) %>% count()

CGI.cravat.result2 %>% group_by(ccheck_transvar_cravat, pcheck_transvar_cravat) %>% count()

write.csv(CGI.cravat.result2,
          'cravat/cravat_hg19_ucsc_analysis_result.csv',
          row.names = F)

CGI.cravat.result2 %>% filter(!ccheck_transvar_cravat |
                                !pcheck_transvar_cravat) %>% write.csv('cravat/ucsc_hg19_match_check_descrepancies.csv')

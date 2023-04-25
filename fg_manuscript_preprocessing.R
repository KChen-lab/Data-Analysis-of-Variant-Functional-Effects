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

# input data
hgvsinput438 <- read.table('438hgvs.txt', header = F)
hgvsinput74 <- read.table('74hgvs.txt', header = F)

CGI <- read.csv('438informativeCGI.csv', header = T)
CGI$Oncogenic.Summary = replace(CGI$Oncogenic.Summary, CGI$Oncogenic.Summary ==
                                  "#N/A", "")
tvinput <-
  read.table('transvar/transvar_output_hg19_ucsc.txt',
             sep = '\t',
             header = T) # offline script

#### QC for protein inference ####
tvinput <- tvinput %>% rowwise() %>%
  mutate(
    p.input = unlist(strsplit(input, ':'))[2],
    p.infer = unlist(strsplit(coordinates.gDNA.cDNA.protein., '/p.'))[2],
    p.check = p.input == p.infer
  )
tvinput %>% group_by(p.check) %>% count()
tvinput %>% filter(p.check == F) %>% select(input, p.input, p.infer) %>% distinct -> protein_check

tvinput <- tvinput %>% filter(!p.infer %in% c('A319delA',
                                              'M319delM',
                                              'E557_N562delEEINGN'))

transvar_aa1 <- tvinput %>% rowwise() %>%
  mutate(
    p.input = unlist(strsplit(input, ':'))[2],
    p.infer = unlist(strsplit(coordinates.gDNA.cDNA.protein., '/p.'))[2],
    p.check = p.input == p.infer,
    g.infer = unlist(strsplit(coordinates.gDNA.cDNA.protein., '/'))[1],
    c.infer = unlist(strsplit(coordinates.gDNA.cDNA.protein., '/'))[2],
    p.aa1infer = unlist(strsplit(coordinates.gDNA.cDNA.protein., '/'))[3]
  ) %>%
  select(input, g.infer, c.infer, p.aa1infer, p.check) %>% distinct()
transvar_aa1_count <- transvar_aa1 %>% group_by(input) %>% count()
transvar_aa1 <-
  full_join(transvar_aa1, transvar_aa1_count, by = 'input')
transvar_aa1 %>% filter(n > 1)

# filter match discrepancies
transvar_aa1 <-
  transvar_aa1 %>% filter(!p.aa1infer %in% c('p.A319delA',
                                             'p.M319delM',
                                             'p.E557_N562delEEINGN')) %>% distinct()
# check
transvar_aa1_count <- transvar_aa1 %>% group_by(input) %>% count()
transvar_aa1 %>% filter(n > 1)

transvar_aa3 <-
  read.table('transvar/transvar_output_hg19_ucsc_aa3.txt',
             sep = '\t',
             header = T)
transvar_aa3 <- transvar_aa3 %>% rowwise() %>%
  mutate(
    p.input = unlist(strsplit(input, ':'))[2],
    p.aa3infer = unlist(strsplit(coordinates.gDNA.cDNA.protein., '/'))[3],
    p.check2 = p.input == p.aa3infer
  )
transvar_aa3 %>% group_by(p.check2) %>% count()
transvar_aa3 <-
  transvar_aa3 %>% select(input, p.aa3infer) %>% distinct()

transvar_aa3_count <- transvar_aa3 %>% group_by(input) %>% count()
transvar_aa3 <-
  full_join(transvar_aa3, transvar_aa3_count, by = 'input')
transvar_aa3 %>% filter(n > 1)

# filter match discrepancies
transvar_aa3 <-
  transvar_aa3 %>% filter(
    !p.aa3infer %in% c(
      'p.Ala319delAla',
      'p.Met319delMet',
      'p.Glu557_Asn562delGluGluIleAsnGlyAsn'
    )
  )


transvar_aa_map <-
  full_join(transvar_aa1, transvar_aa3, by = 'input') %>%
  select(input, c.infer, p.aa1infer, p.check, p.aa3infer) %>% distinct()

tvCGI = full_join(tvinput, CGI, by = c("input" = "Alteration"))
vdata <-
  tvCGI %>% select(input, gene, CHROM, strand, POS, REF, ALT, id, Oncogenic.Summary)

#### cravat tsv format ####
cravat_fmt <-
  vdata %>% select(CHROM, POS, strand, REF, ALT, id, input) %>% distinct()
head(cravat_fmt)

# QC
# check indeterminable positions
cravat_fmt_check1 <- cravat_fmt %>% filter(POS == '.')
vdata_check1 <- vdata %>% filter(input %in% cravat_fmt_check1$input)
tvinput_check1 <-
  tvinput %>% filter(input %in% cravat_fmt_check1$input)
cravat_fmt <- cravat_fmt %>% filter(POS != '.')

# check one variant with multiple entries
cravat_fmt_check2 <- cravat_fmt %>% group_by(input) %>% count()
input_multiple <- cravat_fmt_check2 %>% filter(n > 1)
tvinput_check2 <-
  tvinput %>% filter(input %in% input_multiple$input)

# check variants not identified
set1 = setdiff(hgvsinput438$V1, cravat_fmt$input)
set2 = setdiff(hgvsinput74$V1, cravat_fmt$input)
# need to work on the three missing variants

# output
cravat_fmt$POS <- as.integer(cravat_fmt$POS)
summary(cravat_fmt)
write.table(
  cravat_fmt,
  file = 'cravat/cravat_input_hg19_ucsc.tsv',
  col.names = F,
  sep = '\t',
  row.names = F,
  quote = F
)

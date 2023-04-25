# R version 4.2.2 (2022-10-31) -- "Innocent and Trusting"
# R Core Team (2022). R: A language and environment for statistical computing. R Foundation for Statistical
# Computing, Vienna, Austria. URL https://www.R-project.org/.
# Copyright (C) 2022 The R Foundation for Statistical Computing
# Platform: x86_64-apple-darwin17.0 (64-bit)

library(dplyr)

#### Test 1 ####
dat <- data.frame(
  "Potentially" = c(290, 369),
  "Unknown" = c(9, 109),
  row.names = c("FG_Yes", "FG_No"),
  stringsAsFactors = FALSE
)
test.eq <- fisher.test(dat)
dat
test.eq


#### CGI ####
dat <- data.frame(
  "driver" = c(65, 155),
  "passenger" = c(13, 127),
  row.names = c("FG:Yes", "FG:No"),
  stringsAsFactors = FALSE
)
test.eq <- fisher.test(dat)
test.eq

# Oncogenic.Summary vs Annotator
library(readxl)
FG_438 <-
  read_excel("manuscript/cravat_hg19_ucsc_analysis_result_438_20220327_FG.xlsx")
cravat_output_29variants_mapper <-
  cravat_output_29variants %>% select(Tags, chasm.result.final, vest4.result.final) %>% distinct() %>%
  mutate(variant_check = TRUE)

FG_438_corrected <-
  left_join(FG_438,
            cravat_output_29variants_mapper,
            by = c("Alteration" = "Tags"))

# Final update of prediction results
FG_438_corrected <- FG_438_corrected %>% rowwise() %>%
  mutate(
    VEST4_pred.final = ifelse(
      !is.na(variant_check),
      ifelse(
        vest4.result.final == 1,
        'driver (predicted)',
        ifelse(vest4.result.final == 0, 'passenger',
               "NA")
      ),
      VEST4_pred
    ),
    CHASMplus_pred.final = ifelse(
      !is.na(variant_check),
      ifelse(
        chasm.result.final == 1,
        'driver (predicted)',
        ifelse(chasm.result.final ==
                 0, 'passenger',
               "NA")
      ),
      CHASMplus_pred
    ),
    CScape_pred.final = ifelse(!is.na(variant_check), 'NA', CScape_pred)
  )

sum(FG_438_corrected$VEST4_pred == FG_438_corrected$VEST4_pred.final)
sum(FG_438_corrected$CHASMplus_pred == FG_438_corrected$CHASMplus_pred.final)

cravat_result <- FG_438_corrected

#### vs CHASMplus ####
analysis_dat <-
  cravat_result %>% filter(`Functional Genomics Actionability` != 'No Call',
                           CHASMplus_pred.final != 'NA') %>%
  group_by(`Functional Genomics Actionability`, CHASMplus_pred.final) %>% count()
analysis_dat <-
  analysis_dat %>% pivot_wider(id_cols = `Functional Genomics Actionability`,
                               names_from = CHASMplus_pred.final,
                               values_from = n)
analysis_dat <-
  analysis_dat %>% arrange(desc(`Functional Genomics Actionability`))
analysis_dat <- data.frame(analysis_dat[, 2:3])
analysis_dat
rownames(analysis_dat) <- c('Yes:FG', 'No:FG')
test.eq <- fisher.test(analysis_dat)
test.eq

#### vs VEST4 ####
analysis_dat <-
  cravat_result %>% filter(`Functional Genomics Actionability` != 'No Call',
                           VEST4_pred.final != 'NA') %>%
  group_by(`Functional Genomics Actionability`, VEST4_pred.final) %>% count()
analysis_dat <-
  analysis_dat %>% pivot_wider(id_cols = `Functional Genomics Actionability`,
                               names_from = VEST4_pred.final,
                               values_from = n)
analysis_dat <-
  analysis_dat %>% arrange(desc(`Functional Genomics Actionability`))
analysis_dat <- data.frame(analysis_dat[, 2:3])
analysis_dat
rownames(analysis_dat) <- c('Yes:FG', 'No:FG')
test.eq <- fisher.test(analysis_dat)
test.eq

#### vs CSscape ####
analysis_dat <-
  cravat_result %>% filter(`Functional Genomics Actionability` != 'No Call',
                           CScape_pred.final != 'NA') %>%
  group_by(`Functional Genomics Actionability`, CScape_pred.final) %>% count()
analysis_dat <-
  analysis_dat %>% pivot_wider(id_cols = `Functional Genomics Actionability`,
                               names_from = CScape_pred.final,
                               values_from = n)
analysis_dat <-
  analysis_dat %>% arrange(desc(`Functional Genomics Actionability`))
analysis_dat <- data.frame(analysis_dat[, 2:3])
analysis_dat
rownames(analysis_dat) <- c('Yes:FG', 'No:FG')
test.eq <- fisher.test(analysis_dat)
test.eq

write.csv(
  cravat_result,
  'manuscript/cravat_hg19_ucsc_analysis_result_438_final.csv',
  row.names = F
)

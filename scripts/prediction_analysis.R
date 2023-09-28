###librraies 
library(tidyverse)
library(readr)
library(dplyr)


### reading in the data 
#prediction from TB-Profiler 
latest_results <- read_csv("joined_results.csv")
latest_results <- latest_results[!duplicated(latest_results$sample_accession),]

#final phenotypic data 
full_phenotypic_data <- read_csv("full_phenotypic_dataset.csv")
full_phenotypic_data <- full_phenotypic_data[!duplicated(full_phenotypic_data$sample_accession),]

# cleaning values

#Bedaquiline 
full_phenotypic_data$bedaquiline[full_phenotypic_data$bedaquiline == 'N.T'] <- NA
full_phenotypic_data$bedaquiline[full_phenotypic_data$bedaquiline == 'not tested'] <- NA
full_phenotypic_data$bedaquiline[full_phenotypic_data$bedaquiline == 'R'] <- '2'
full_phenotypic_data$bedaquiline[full_phenotypic_data$bedaquiline == 'S'] <- '0.5'


#clofazimine
full_phenotypic_data$clofazimine[full_phenotypic_data$clofazimine == 'R'] <- '2'
full_phenotypic_data$clofazimine[full_phenotypic_data$clofazimine == 'S'] <- '0.5'


# #linezolid
full_phenotypic_data$linezolid[full_phenotypic_data$linezolid == 'R'] <- '2'
full_phenotypic_data$linezolid[full_phenotypic_data$linezolid == 'S'] <- '0.5'


# #delamanid 
full_phenotypic_data$delamanid[full_phenotypic_data$delamanid == 'R'] <- '0.12'
full_phenotypic_data$delamanid[full_phenotypic_data$delamanid == 'S'] <- '0.03'
full_phenotypic_data$delamanid[full_phenotypic_data$delamanid == 'not tested'] <- NA



#checking data has been cleaned

#full_phenotypic_data$bedaquiline <- as.factor(full_phenotypic_data$bedaquiline)
#summary(full_phenotypic_data$bedaquiline)

#full_phenotypic_data$clofazimine <- as.factor(full_phenotypic_data$clofazimine)
#summary(full_phenotypic_data$clofazimine)

#full_phenotypic_data$linezolid <- as.factor(full_phenotypic_data$linezolid)
#summary(full_phenotypic_data$linezolid)

#full_phenotypic_data$delamanid <- as.factor(full_phenotypic_data$delamanid)
#summary(full_phenotypic_data$delamanid)

#full_phenotypic_data$pretomanid <- as.factor(full_phenotypic_data$pretomanid)
#summary(full_phenotypic_data$pretomanid)


#creating long form of dataset 

#pivot phenotypic_data from wide to long
phenotypic_long <- pivot_longer(data = full_phenotypic_data, cols = c(bedaquiline, clofazimine, linezolid, delamanid, pretomanid), 
                                names_to = "Drug",
                                values_to = "MIC")


#making the column numeric 
phenotypic_long$MIC <- as.numeric(phenotypic_long$MIC)


#don convert to 1 just yet, keep as R and S 
## applying MIC cut off's for the drugs in phenotypic dataset 
phenotypic_long <- phenotypic_long %>% 
  mutate(resistance = case_when(
    Drug == 'bedaquiline' & MIC >= 1 ~ "R",
    Drug == 'bedaquiline' & MIC < 1 ~ "S",
    Drug == 'clofazimine' & MIC >= 1 ~ "R",
    Drug == 'clofazimine' & MIC < 1 ~ "S",
    Drug == 'linezolid' & MIC >= 1 ~ "R",
    Drug == 'linezolid' & MIC < 1 ~ "S",
    Drug == 'delamanid' & MIC >= 0.06 ~ "R",
    Drug == 'delamanid' & MIC < 0.06 ~ "S"))
    #is.character(MIC) ~ MIC))


summary(as.factor(phenotypic_long$resistance))


phenotypic_long$resistance[phenotypic_long$resistance == 'R'] <- '1'
phenotypic_long$resistance[phenotypic_long$resistance == 'S'] <- '0'


###prediction dataset 

##selecting relevant columns 
prediction <- latest_results %>%  
  select(sample_id, study_name, PMID, sample_accession, study_accession, bedaquiline_y, clofazimine_y, linezolid_y, delamanid_y, pretomanid_y)


## pivot to long form 

prediction_long <- pivot_longer(data = prediction, cols = c(bedaquiline_y, clofazimine_y, linezolid_y, delamanid_y, pretomanid_y), 
                                names_to = "Drug",
                                values_to = "Mutation")

#both have same dimension, perfect 
dim(phenotypic_long)
dim(prediction_long)

#changing name of drugs
prediction_long$Drug[prediction_long$Drug == 'bedaquiline_y'] <- 'bedaquiline'
prediction_long$Drug[prediction_long$Drug == 'clofazimine_y'] <- 'clofazimine'
prediction_long$Drug[prediction_long$Drug == 'linezolid_y'] <- 'linezolid'
prediction_long$Drug[prediction_long$Drug == 'delamanid_y'] <- 'delamanid'
prediction_long$Drug[prediction_long$Drug == 'pretomanid_y'] <- 'pretomanid'


prediction_long_need <- prediction_long %>% 
  select(sample_accession, resistance)

#if a mutation is present, it is resitant. otherwise it is susceptible 

prediction_long$resistance <- ifelse( is.na(prediction_long$Mutation), 0, 1) 

#sanity check 
identical(prediction_long$sample_accession, phenotypic_long$sample_accession)


#### bedaquiline ####

#create dataset with id, resistance for bedaquline phenotypic and prediction 

bdq_pheno <- phenotypic_long %>%  
  filter(Drug == 'bedaquiline' & !is.na(resistance)) %>% 
  select(sample_accession, resistance) %>% 
  dplyr::semi_join(prediction_long, by = 'sample_accession')


bdq_pred <- prediction_long %>%  
  filter(Drug == 'bedaquiline') %>% 
  select(sample_accession, resistance) %>% 
  dplyr::semi_join(bdq_pheno, by = 'sample_accession')
  

dim(bdq_pheno)
nrow(unique(bdq_pheno))
dim(bdq_pred)
nrow(unique(bdq_pred))

clo_pheno <- phenotypic_long %>%  
  filter(Drug == 'clofazimine' & !is.na(resistance)) %>% 
  select(sample_accession, resistance) %>% 
  dplyr::semi_join(prediction_long, by = 'sample_accession')

clo_pred <- prediction_long %>%  
  filter(Drug == 'clofazimine') %>% 
  select(sample_accession, resistance) %>% 
  dplyr::semi_join(clo_pheno, by = 'sample_accession')



lzd_pheno <- phenotypic_long %>%  
  filter(Drug == 'linezolid' & !is.na(resistance)) %>% 
  select(sample_accession, resistance) %>% 
  dplyr::semi_join(prediction_long, by = 'sample_accession')

lzd_pred <- prediction_long %>%  
  filter(Drug == 'linezolid') %>% 
  select(sample_accession, resistance) %>% 
  dplyr::semi_join(lzd_pheno, by = 'sample_accession')




dlm_pheno <- phenotypic_long %>%  
  filter(Drug == 'delamanid' & !is.na(resistance)) %>% 
  select(sample_accession, resistance) %>% 
  dplyr::semi_join(prediction_long, by = 'sample_accession')

dlm_pred <- prediction_long %>%  
  filter(Drug == 'delamanid') %>% 
  select(sample_accession, resistance) %>% 
  dplyr::semi_join(dlm_pheno, by = 'sample_accession')





pa_pheno <- phenotypic_long %>%  
  filter(Drug == 'pretomanid'& !is.na(resistance)) %>% 
  select(sample_accession, resistance) %>% 
  dplyr::semi_join(prediction_long, by = 'sample_accession')



#how many accesions in common -> all in common!!

ac <- full_phenotypic_data %>% 
  inner_join(latest_results, by= 'sample_accession')
View(ac)


  
  
#### confusion matrix ####
library(epiR)

#bedaquiline

bdq_table <- bdq_pheno %>% 
  inner_join(bdq_pred, by = 'sample_accession') %>% 
  rename( dis = 'resistance.x',
          tes = 'resistance.y' ) 

tmp.df01 <- bdq_table %>%
  mutate(dis = factor(dis, levels = c(1,0), labels = c("Dis+","Dis-"))) %>%
  mutate(tes = factor(tes, levels = c(1,0), labels = c("Test+","Test-"))) %>%
  group_by(tes, dis) %>%
  summarise(n = n())
tmp.df01


table(phenotypic = bdq_table$dis, genotypic= bdq_table$tes)


## View the data in conventional 2 by 2 table format:
pivot_wider(tmp.df01, id_cols = c(tes), names_from = dis, values_from = n)

rval.tes01 <- epi.tests(tmp.df01, method = "wilson", digits = 2, 
                        conf.level = 0.95)
summary(rval.tes01)


#clofazimine
clo_table <- clo_pheno %>% 
  inner_join(clo_pred, by = 'sample_accession') %>% 
  rename( dis = 'resistance.x',
          tes = 'resistance.y' ) 

tmp.df02 <- clo_table %>%
  mutate(dis = factor(dis, levels = c(1,0), labels = c("Dis+","Dis-"))) %>%
  mutate(tes = factor(tes, levels = c(1,0), labels = c("Test+","Test-"))) %>%
  group_by(tes, dis) %>%
  summarise(n = n())
tmp.df02


table(clo_table$dis, clo_table$tes)


## View the data in conventional 2 by 2 table format:
pivot_wider(tmp.df02, id_cols = c(tes), names_from = dis, values_from = n)

rval.tes02 <- epi.tests(tmp.df02, method = "wilson", digits = 2, 
                        conf.level = 0.95)
summary(rval.tes02)


#linezolid
lzd_table <- lzd_pheno %>% 
  inner_join(lzd_pred, by = 'sample_accession') %>% 
  rename( dis = 'resistance.x',
          tes = 'resistance.y' )


tmp.df03 <- lzd_table %>%
  mutate(dis = factor(dis, levels = c(1,0), labels = c("Dis+","Dis-"))) %>%
  mutate(tes = factor(tes, levels = c(1,0), labels = c("Test+","Test-"))) %>%
  group_by(tes, dis) %>%
  summarise(n = n())
tmp.df03


table(lzd_table$dis, lzd_table$tes)


## View the data in conventional 2 by 2 table format:
pivot_wider(tmp.df03, id_cols = c(tes), names_from = dis, values_from = n)

rval.tes03 <- epi.tests(tmp.df03, method = "wilson", digits = 2, 
                        conf.level = 0.95)
summary(rval.tes03)


#delamanid 
dlm_table <- dlm_pheno %>% 
  inner_join(dlm_pred, by = 'sample_accession') %>% 
  rename( dis = 'resistance.x',
          tes = 'resistance.y' )


tmp.df04 <- dlm_table %>%
  mutate(dis = factor(dis, levels = c(1,0), labels = c("Dis+","Dis-"))) %>%
  mutate(tes = factor(tes, levels = c(1,0), labels = c("Test+","Test-"))) %>%
  group_by(tes, dis) %>%
  summarise(n = n())
tmp.df04


table(dlm_table$dis, dlm_table$tes)


## View the data in conventional 2 by 2 table format:
pivot_wider(tmp.df04, id_cols = c(tes), names_from = dis, values_from = n)

rval.tes04 <- epi.tests(tmp.df04, method = "wilson", digits = 2, 
                        conf.level = 0.95)
summary(rval.tes04)


#pretomanid

pa_table <- pa_pheno %>% 
  rename( dis = 'resistance.x',
          tes = 'resistance.y' ) 


tmp.df05 <- pa_table %>%
  mutate(dis = factor(dis, levels = c(1,0), labels = c("Dis+","Dis-"))) %>%
  mutate(tes = factor(tes, levels = c(1,0), labels = c("Test+","Test-"))) %>%
  group_by(tes, dis) %>%
  summarise(n = n())
tmp.df05


table(pa_table$dis, pa_table$tes)

## View the data in conventional 2 by 2 table format:
pivot_wider(tmp.df05, id_cols = c(tes), names_from = dis, values_from = n)

rval.tes05 <- epi.tests(tmp.df05, method = "wilson", digits = 2, 
                        conf.level = 0.95)
summary(rval.tes05)







###exploratory analysis 
View(full_phenotypic_data)
nrow(full_phenotypic_data)

#total not NA for bedaquiline is 8,800 (post-cleaning of NT and not tested etc)
summary(as.factor(full_phenotypic_data$bedaquiline))

#total not NA for clofazimine is 11010
summary(as.factor(full_phenotypic_data$clofazimine))

#total not NA for linezolid = 10206
summary(as.factor(full_phenotypic_data$linezolid))

#total not NA for delmanid = 8150
summary(as.factor(full_phenotypic_data$delamanid))

#all are NA for pretomanid 
summary(as.factor(full_phenotypic_data$pretomanid))


#creating a table 

#creating a df with relavant data 
  
df <- data.frame(Drug  = c("bedaquiline", "clofazimine", "linezolid", "delamanid"),
                R = c("136", "216", "143", "146"), 
                S = c("8664", "10794", "10063", "8003"), 
                Prevalence = c("1.55", "1.96", "1.4", "1.79"), 
                TP = c("136", "212", "62",  "105"),
                FN = c("0", "4", "81", "41"),
                TN = c("87", "138", "10034", "2861"),
                FP = c("8577", "10656", "29", "5142"), 
                Sensitivity = c("1 (0.973 - 1)", "0.98 (0.95 - 0.99)", "0.43 (0.36 - 0.52)", "0.719 (0.641 - 0.786)"),
                Specificty = c("0.01 (0.0081 - 0.012)", "0.0127 (0.011 - 0.015)", "0.997 (0.996 - 0.998)", "0.357 (0.347 - 0.368)"), 
                PPV = c(" 0.0156", "0.0195", "0.681", "0.02"), 
                NPV = c("1", "0.97", "0.992", "0.98")) 

df <- df %>% 
  rename( "Prevalence (%)" = Prevalence)
  
  
  
library(tidyverse)
library(gt)
library(espnscrapeR)
library(ggplot2)

library(gt)
library(webshot2)

tablee <- df %>% 
  gt(auto_align = FALSE) %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()) %>%
  tab_style(style = cell_fill(color = "grey"),
            locations = cells_body(rows = seq(1, 4, 2)))
  
df %>% 
  gt(auto_align = FALSE) %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()) %>% 
  cols_label(
    Prevalence = "Prevalence (%)") %>% 
  tab_options(
    column_labels.border.top.color = "black",
    column_labels.border.bottom.color = "black",
    table.border.bottom.width = px(3),
    table.border.bottom.color = "black",
    table_body.hlines.style = "white") %>% 
  tab_footnote(
    footnote = md("Numbers in brackets refer to the *95% Confidence Interval*."),
    locations = cells_column_labels(columns = c(Sensitivity, Specificty))) %>% 
  tab_footnote(
    footnote = md("R stands for *Number of Resistant Strains*."),
    locations = cells_column_labels(columns = R)) %>% 
  tab_footnote(
    footnote = md("S stands for *Number of Susceptible Strains*."),
    locations = cells_column_labels(columns = S)) %>% 
  tab_footnote(
    footnote = md("Prevalance(%) stands for *Prevalance of Resistant strains in the sample of strains tested for that drug*."),
    locations = cells_column_labels(columns =  Prevalence)) %>% 
  tab_footnote(
    footnote = md("TP stands for *Number of True Positive*."),
    locations = cells_column_labels(columns = TP)) %>% 
  tab_footnote(
    footnote = md("FN stands for *Number of False Negative*."),
    locations = cells_column_labels(columns = FN)) %>% 
  tab_footnote(
    footnote = md("TN stands for *Number of True Negative*."),
    locations = cells_column_labels(columns = TN))  %>% 
  tab_footnote(
    footnote = md("FP stands for *Number of False Positive*."),
    locations = cells_column_labels(columns = FP)) %>% 
  tab_footnote(
    footnote = md("PPV stands for *Postive Predictive Value*."),
    locations = cells_column_labels(columns = PPV)) %>% 
  tab_footnote(
    footnote = md("NPV stands for *Negative Predictive Value*."),
    locations = cells_column_labels(columns = NPV)) %>% 
    tab_options(footnotes.multiline = FALSE) %>% 
  tab_options(data_row.padding = px(10)) %>% 
  cols_width(Sensitivity ~ px(150),
             Specificty ~ px(150),
             Drug ~ px(100),
             Prevalence  ~ px(110),
             everything() ~ px(50)) %>% 
  tab_options(table.font.size = px(12L))
  














gtsave(tablee, "tableee.png")
tablee

gtsave(tablee, "tab_1.docx", expand = 10)

?basic_theme()

?gt()





#### false predictors mutations ####
View(latest_results)

##claening MIC 
#Bedaquiline 
latest_results$bedaquiline_x[latest_results$bedaquiline_x == 'N.T'] <- NA
latest_results$bedaquiline_x[latest_results$bedaquiline_x == 'not tested'] <- NA
latest_results$bedaquiline_x[latest_results$bedaquiline_x == 'R'] <- '2'
latest_results$bedaquiline_x[latest_results$bedaquiline_x == 'S'] <- '0.5'


#clofazimine
latest_results$clofazimine_x[latest_results$clofazimine_x == 'R'] <- '2'
latest_results$clofazimine_x[latest_results$clofazimine_x == 'S'] <- '0.5'


# #linezolid
latest_results$linezolid_x[latest_results$linezolid_x == 'R'] <- '2'
latest_results$linezolid_x[latest_results$linezolid_x == 'S'] <- '0.5'


# #delamanid 
latest_results$delamanid_x[latest_results$delamanid_x == 'R'] <- '0.12'
latest_results$delamanid_x[latest_results$delamanid_x == 'S'] <- '0.03'
latest_results$delamanid_x[latest_results$delamanid_x == 'not tested'] <- NA


latest_results$bedaquiline_x <- as.numeric(latest_results$bedaquiline_x )
latest_results$clofazimine_x <- as.numeric(latest_results$clofazimine_x )
latest_results$linezolid_x <- as.numeric(latest_results$linezolid_x )
latest_results$delamanid_x <- as.numeric(latest_results$delamanid_x )

summary(as.factor(latest_results$bedaquiline_x))

#bdq

#selecting rows where S and mutation 
mutations_bdq <- latest_results %>% 
  #selecting susceptible and mutation present 
  filter(bedaquiline_x < 1 & !is.na(bedaquiline_y)) %>% 
  select(bedaquiline_y)
 

summary(as.factor(mutations_bdq$bedaquiline_y))

library(tidyr)
mutatations_bdq_clean <-  mutations_bdq %>% 
  separate_rows(bedaquiline_y, sep=",")

#stripping whitespace 
library(stringr)
df_new <-as.data.frame(apply(mutatations_bdq_clean,2, str_remove_all, " "))

summary(as.factor(df_new$bedaquiline_y))

df_new <- df_new %>% 
  group_by(bedaquiline_y) %>% 
  dplyr::mutate(count_susceptible = n()) 

df_new <- unique(df_new)

df_new %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    bedaquiline_y = md("**Mutation**"),
    count_susceptible = md("**Number of Susceptible Samples with mutation**")) %>% 
  tab_options(., container.width = 1000, container.height = 1000)



##resistant samples for bdq 
mutations_bdq2 <- latest_results %>% 
  #selecting susceptible and mutation present 
  filter(bedaquiline_x >= 1 & !is.na(bedaquiline_y)) %>% 
  select(bedaquiline_y)

summary(as.factor(mutations_bdq2$bedaquiline_y))

mutatations_bdq_clean2 <-  mutations_bdq2 %>% 
  separate_rows(bedaquiline_y, sep=",")

#stripping whitespace 
library(stringr)
df_new_2 <-as.data.frame(apply(mutatations_bdq_clean2,2, str_remove_all, " "))

summary(as.factor(df_new_2$bedaquiline_y))

df_new_2 <- df_new_2 %>% 
  group_by(bedaquiline_y) %>% 
  dplyr::mutate(count_resistant = n()) 

df_new_2 <- unique(df_new_2)

df_new_2 %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    bedaquiline_y = md("**Mutation**"),
    count_resistant = md("**Number of Resistant Samples with mutation**")) %>% 
  tab_options(., container.width = 1000, container.height = 1000)



final_df <- merge(df_new, df_new_2, by= "bedaquiline_y", all.x = TRUE, all.y=TRUE)


#converting NA to O 
final_df[is.na(final_df)] <- 0

#creating new columns for R and S without mutations
final_df$notin_resis <- 136 - final_df$count_resistant
final_df$notin_sus <- 8664 - final_df$count_susceptible

#column for odds ratio 
final_df$odds_ratio <- (final_df$count_resistant * final_df$notin_sus) / (final_df$count_susceptible * final_df$notin_resis)
final_df$lower_ci <- exp(log(final_df$odds_ratio) - 1.96 * (sqrt(1/final_df$count_resistant + 1/final_df$notin_sus  + 1/final_df$count_susceptible   + 1/final_df$notin_resis)   )) 
final_df$upper_ci <- exp(log(final_df$odds_ratio) + 1.96 * (sqrt(1/final_df$count_resistant + 1/final_df$notin_sus  + 1/final_df$count_susceptible   + 1/final_df$notin_resis)   )) 

f <- chromote::default_chromote_object() #get the f object
f$close()


final_df %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    bedaquiline_y = md("**Mutation**"),
    count_susceptible = md("**Number of Susceptible Samples with mutation**"),
    count_resistant = md("**Number of Resistant Samples with mutation**"),
    notin_sus = md("**Number of Susceptible Samples without mutation**"),
    notin_resis = md("**Number of Resistant Samples without mutation**"),
    odds_ratio = md("**Odds Ratio**"),) %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()) %>% 
  tab_options(
    column_labels.border.top.color = "black",
    column_labels.border.bottom.color = "black",
    table.border.bottom.width = px(3),
    table.border.bottom.color = "black",
    table_body.hlines.style = "white") %>% 
  tab_options(data_row.padding = px(3))



#  tab_options(., container.width = 1000, container.height = 1000) %>% 
tab_options(table.font.size = px(10L)) %>% 
  cols_width(everything() ~ px(150))   gt(auto_align = FALSE) 

#%>%  gtsave(filename = "bdq_problematic2.png")

final_df %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    bedaquiline_y = md("**Mutation**"),
    count_susceptible = md("**Susceptible Samples with mutation**"),
    count_resistant = md("**Resistant Samples with mutation**"),
    notin_sus = md("**Susceptible Samples without mutation**"),
    notin_resis = md("**Resistant Samples without mutation**"),
    odds_ratio = md("**Odds Ratio**"),
    lower_ci = md("**Lower 95% Confidence Interval**"),
    upper_ci = md("**Upper 95% Confidence Interval**"),) %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()) %>% 
  tab_options(
    column_labels.border.top.color = "black",
    column_labels.border.bottom.color = "black",
    table.border.bottom.width = px(3),
    table.border.bottom.color = "black",
    table_body.hlines.style = "white") %>% 
    tab_options(data_row.padding = px(0.3)) %>% 
  cols_width(bedaquiline_y ~ px(170),
             everything() ~ px(100)) %>% 
  tab_options(table.font.size = px(12L))









final_df %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    bedaquiline_y = md("**Mutation**"),
    count_susceptible = md("**Susceptible Samples with mutation**"),
    count_resistant = md("**Resistant Samples with mutation**"),
    notin_sus = md("**Susceptible Samples without mutation**"),
    notin_resis = md("**Resistant Samples without mutation**"),
    odds_ratio = md("**Odds Ratio**"),
    lower_ci = md("**Lower 95% CI**"),
    upper_ci = md("**Upper 95% CI**"),) %>% 
  tab_options(., container.width = 1000, container.height = 1000) %>% 
  tab_options(
    table.font.size = px(10L)) %>% 
  cols_width(everything() ~ px(100)) %>% 
  gtsave(filename = "bdq_problematic2.png")



#clo 
#selecting rows where S and mutation 
mutations_clo <- latest_results %>% 
  #selecting susceptible and mutation present 
  filter(clofazimine_x < 1 & !is.na(clofazimine_y)) %>% 
  select(clofazimine_y)


library(tidyr)
mutatations_clo_clean <-  mutations_clo %>% 
  separate_rows(clofazimine_y, sep=",")

#stripping whitespace 
library(stringr)
df_new <-as.data.frame(apply(mutatations_clo_clean,2, str_remove_all, " "))

summary(as.factor(df_new$clofazimine_y))

df_new <- df_new %>% 
  group_by(clofazimine_y) %>% 
  dplyr::mutate(count_susceptible = n()) 

df_new <- unique(df_new)

df_new %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    clofazimine_y = md("**Mutation in Clofazimine Susceptible Samples**"),
    count = md("**Number of Samples**")) %>% 
  tab_options(
    table.font.size = px(10L)
  )

##resistant samples for clo
mutations_clo2 <- latest_results %>% 
  #selecting susceptible and mutation present 
  filter(clofazimine_x >= 1 & !is.na(clofazimine_y)) %>% 
  select(clofazimine_y)

summary(as.factor(mutations_clo2$clofazimine_y))

mutatations_clo_clean2 <-  mutations_clo2 %>% 
  separate_rows(clofazimine_y, sep=",")

#stripping whitespace 
library(stringr)
df_new_2 <-as.data.frame(apply(mutatations_clo_clean2,2, str_remove_all, " "))


df_new_2 <- df_new_2 %>% 
  group_by(clofazimine_y) %>% 
  dplyr::mutate(count_resistant = n()) 

df_new_2 <- unique(df_new_2)

df_new_2 %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    clofazimine_y = md("**Mutation**"),
    count_resistant = md("**Number of Resistant Samples with mutation**")) %>% 
  tab_options(., container.width = 1000, container.height = 1000)



final_df <- merge(df_new, df_new_2, by= "clofazimine_y", all.x = TRUE, all.y=TRUE)


#converting NA to O 
final_df[is.na(final_df)] <- 0

#creating new columns for R and S without mutations
final_df$notin_resis <- 212 - final_df$count_resistant
final_df$notin_sus <- 10656  - final_df$count_susceptible

#column for odds ratio 
final_df$odds_ratio <- (final_df$count_resistant * final_df$notin_sus) / (final_df$count_susceptible * final_df$notin_resis)
final_df$lower_ci <- exp(log(final_df$odds_ratio) - 1.96 * (sqrt(1/final_df$count_resistant + 1/final_df$notin_sus  + 1/final_df$count_susceptible   + 1/final_df$notin_resis)   )) 
final_df$upper_ci <- exp(log(final_df$odds_ratio) + 1.96 * (sqrt(1/final_df$count_resistant + 1/final_df$notin_sus  + 1/final_df$count_susceptible   + 1/final_df$notin_resis)   )) 




final_df %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    clofazimine_y = md("**Mutation**"),
    count_susceptible = md("**Number of Susceptible Samples with mutation**"),
    count_resistant = md("**Number of Resistant Samples with mutation**"),
    notin_sus = md("**Number of Susceptible Samples without mutation**"),
    notin_resis = md("**Number of Resistant Samples without mutation**"),
    odds_ratio = md("**Odds Ratio**"),) %>% 
  tab_options(., container.width = 1000, container.height = 1000) %>% 
  tab_options(
    table.font.size = px(10L)) %>% 
  cols_width(everything() ~ px(150)) %>% 
  gtsave(filename = "clo_problematic.png")

final_df %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
   clofazimine_y = md("**Mutation**"),
    count_susceptible = md("**Susceptible Samples with mutation**"),
    count_resistant = md("**Resistant Samples with mutation**"),
    notin_sus = md("**Susceptible Samples without mutation**"),
    notin_resis = md("**Resistant Samples without mutation**"),
    odds_ratio = md("**Odds Ratio**"),
    lower_ci = md("**Lower 95% CI**"),
    upper_ci = md("**Upper 95% CI**"),) %>% 
  tab_options(., container.width = 1000, container.height = 1000) %>% 
  tab_options(
    table.font.size = px(10L)) %>% 
  cols_width(everything() ~ px(100)) %>% 
  gtsave(filename = "clo_problematic2.png")


final_df %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    clofazimine_y = md("**Mutation**"),
    count_susceptible = md("**Susceptible Samples with mutation**"),
    count_resistant = md("**Resistant Samples with mutation**"),
    notin_sus = md("**Susceptible Samples without mutation**"),
    notin_resis = md("**Resistant Samples without mutation**"),
    odds_ratio = md("**Odds Ratio**"),
    lower_ci = md("**Lower 95% Confidence Interval**"),
    upper_ci = md("**Upper 95% Confidence Interval**"),) %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()) %>% 
  tab_options(
    column_labels.border.top.color = "black",
    column_labels.border.bottom.color = "black",
    table.border.bottom.width = px(3),
    table.border.bottom.color = "black",
    table_body.hlines.style = "white") %>% 
  tab_options(data_row.padding = px(0.5)) %>% 
  cols_width(clofazimine_y ~ px(170),
             everything() ~ px(100)) %>% 
  tab_options(table.font.size = px(12L))

#lzd
#selecting rows where S and mutation 
mutations_lzd <- latest_results %>% 
  #selecting susceptible and mutation present 
  filter(linezolid_x < 1 & !is.na(linezolid_y)) %>% 
  select(linezolid_y)


library(tidyr)
mutatations_lzd_clean <-  mutations_lzd %>% 
  separate_rows(linezolid_y, sep=",")

#stripping whitespace 
library(stringr)
df_new <-as.data.frame(apply(mutatations_lzd_clean,2, str_remove_all, " "))

summary(as.factor(df_new$linezolid_y))

df_new <- df_new %>% 
  group_by(linezolid_y) %>% 
  dplyr::mutate(count_susceptible = n()) 

df_new <- unique(df_new)

df_new %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    linezolid_y = md("**Mutation in linezolid Susceptible Samples**"),
    count = md("**Number of Samples**")) %>% 
  tab_options(
    table.font.size = px(10L)
  )

##resistant samples for lzd
mutations_lzd2 <- latest_results %>% 
  #selecting susceptible and mutation present 
  filter(linezolid_y >= 1 & !is.na(linezolid_y)) %>% 
  select(linezolid_y)

summary(as.factor(mutations_lzd2$linezolid_y))

mutatations_lzd_clean2 <-  mutations_lzd2 %>% 
  separate_rows(linezolid_y, sep=",")

#stripping whitespace 
library(stringr)
df_new_2 <-as.data.frame(apply(mutatations_lzd_clean2,2, str_remove_all, " "))


df_new_2 <- df_new_2 %>% 
  group_by(linezolid_y) %>% 
  dplyr::mutate(count_resistant = n()) 

df_new_2 <- unique(df_new_2)

df_new_2 %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    linezolid_y = md("**Mutation**"),
    count_resistant = md("**Number of Resistant Samples with mutation**")) %>% 
  tab_options(., container.width = 1000, container.height = 1000)



final_df <- merge(df_new, df_new_2, by= "linezolid_y", all.x = TRUE, all.y=TRUE)


#converting NA to O 
final_df[is.na(final_df)] <- 0

#creating new columns for R and S without mutations
final_df$notin_resis <- 143 - final_df$count_resistant
final_df$notin_sus <- 10063  - final_df$count_susceptible

#column for odds ratio 
final_df$odds_ratio <- (final_df$count_resistant * final_df$notin_sus) / (final_df$count_susceptible * final_df$notin_resis)
final_df$lower_ci <- exp(log(final_df$odds_ratio) - 1.96 * (sqrt(1/final_df$count_resistant + 1/final_df$notin_sus  + 1/final_df$count_susceptible   + 1/final_df$notin_resis)   )) 
final_df$upper_ci <- exp(log(final_df$odds_ratio) + 1.96 * (sqrt(1/final_df$count_resistant + 1/final_df$notin_sus  + 1/final_df$count_susceptible   + 1/final_df$notin_resis)   )) 




final_df %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    linezolid_y = md("**Mutation**"),
    count_susceptible = md("**Number of Susceptible Samples with mutation**"),
    count_resistant = md("**Number of Resistant Samples with mutation**"),
    notin_sus = md("**Number of Susceptible Samples without mutation**"),
    notin_resis = md("**Number of Resistant Samples without mutation**"),
    odds_ratio = md("**Odds Ratio**"),) %>% 
  tab_options(., container.width = 1000, container.height = 1000) %>% 
  tab_options(
    table.font.size = px(10L)) %>% 
  cols_width(everything() ~ px(150)) %>% 
  gtsave(filename = "lzd_problematic.png")

final_df %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    linezolid_y = md("**Mutation**"),
    count_susceptible = md("**Susceptible Samples with mutation**"),
    count_resistant = md("**Resistant Samples with mutation**"),
    notin_sus = md("**Susceptible Samples without mutation**"),
    notin_resis = md("**Resistant Samples without mutation**"),
    odds_ratio = md("**Odds Ratio**"),
    lower_ci = md("**Lower 95% CI**"),
    upper_ci = md("**Upper 95% CI**"),) %>% 
  tab_options(., container.width = 1000, container.height = 1000) %>% 
  tab_options(
    table.font.size = px(10L)) %>% 
  cols_width(everything() ~ px(100)) %>% 
  gtsave(filename = "lzd_problematic2.png")


final_df %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    linezolid_y = md("**Mutation**"),
    count_susceptible = md("**Susceptible Samples with mutation**"),
    count_resistant = md("**Resistant Samples with mutation**"),
    notin_sus = md("**Susceptible Samples without mutation**"),
    notin_resis = md("**Resistant Samples without mutation**"),
    odds_ratio = md("**Odds Ratio**"),
    lower_ci = md("**Lower 95% Confidence Interval**"),
    upper_ci = md("**Upper 95% Confidence Interval**"),) %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()) %>% 
  tab_options(
    column_labels.border.top.color = "black",
    column_labels.border.bottom.color = "black",
    table.border.bottom.width = px(3),
    table.border.bottom.color = "black",
    table_body.hlines.style = "white") %>% 
  tab_options(data_row.padding = px(0.6)) %>% 
  cols_width(linezolid_y ~ px(180),
             everything() ~ px(100)) %>% 
  tab_options(table.font.size = px(12L))


#dlm
#selecting rows where S and mutation 
mutations_dlm <- latest_results %>% 
  #selecting susceptible and mutation present 
  filter(delamanid_x < 1 & !is.na(delamanid_y)) %>% 
  select(delamanid_y)


library(tidyr)
mutatations_dlm_clean <-  mutations_dlm %>% 
  separate_rows(delamanid_y, sep=",")

#stripping whitespace 
library(stringr)
df_new <-as.data.frame(apply(mutatations_dlm_clean,2, str_remove_all, " "))

summary(as.factor(df_new$delamanid_y))

df_new <- df_new %>% 
  group_by(delamanid_y) %>% 
  dplyr::mutate(count_susceptible = n()) 

df_new <- unique(df_new)

df_new %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    delamanid_y = md("**Mutation in delamanid Susceptible Samples**"),
    count = md("**Number of Samples**")) %>% 
  tab_options(
    table.font.size = px(10L)
  )
# 


##resistant samples for dlm
mutations_dlm2 <- latest_results %>% 
  #selecting susceptible and mutation present 
  filter(delamanid_x >= 1 & !is.na(delamanid_y)) %>% 
  select(delamanid_y)

mutatations_dlm_clean2 <-  mutations_dlm2 %>% 
  separate_rows(delamanid_y, sep=",")

#stripping whitespace 
library(stringr)
df_new_2 <-as.data.frame(apply(mutatations_dlm_clean2,2, str_remove_all, " "))


df_new_2 <- df_new_2 %>% 
  group_by(delamanid_y) %>% 
  dplyr::mutate(count_resistant = n()) 

df_new_2 <- unique(df_new_2)

df_new_2 %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    delamanid_y = md("**Mutation**"),
    count_resistant = md("**Number of Resistant Samples with mutation**")) %>% 
  tab_options(., container.width = 1000, container.height = 1000)



final_df <- merge(df_new, df_new_2, by= "delamanid_y", all.x = TRUE, all.y=TRUE)


#converting NA to O 
final_df[is.na(final_df)] <- 0

#creating new columns for R and S without mutations
final_df$notin_resis <- 146 - final_df$count_resistant
final_df$notin_sus <- 8003  - final_df$count_susceptible

#column for odds ratio 
final_df$odds_ratio <- (final_df$count_resistant * final_df$notin_sus) / (final_df$count_susceptible * final_df$notin_resis)
final_df$lower_ci <- exp(log(final_df$odds_ratio) - 1.96 * (sqrt(1/final_df$count_resistant + 1/final_df$notin_sus  + 1/final_df$count_susceptible   + 1/final_df$notin_resis)   )) 
final_df$upper_ci <- exp(log(final_df$odds_ratio) + 1.96 * (sqrt(1/final_df$count_resistant + 1/final_df$notin_sus  + 1/final_df$count_susceptible   + 1/final_df$notin_resis)   )) 




final_df %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    delamanid_y = md("**Mutation**"),
    count_susceptible = md("**Number of Susceptible Samples with mutation**"),
    count_resistant = md("**Number of Resistant Samples with mutation**"),
    notin_sus = md("**Number of Susceptible Samples without mutation**"),
    notin_resis = md("**Number of Resistant Samples without mutation**"),
    odds_ratio = md("**Odds Ratio**"),) %>% 
  tab_options(., container.width = 1000, container.height = 1000) %>% 
  tab_options(
    table.font.size = px(10L)) %>% 
  cols_width(everything() ~ px(150)) %>% 
  gtsave(filename = "dlm_problematic.png")

final_df %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    delamanid_y = md("**Mutation**"),
    count_susceptible = md("**Susceptible Samples with mutation**"),
    count_resistant = md("**Resistant Samples with mutation**"),
    notin_sus = md("**Susceptible Samples without mutation**"),
    notin_resis = md("**Resistant Samples without mutation**"),
    odds_ratio = md("**Odds Ratio**"),
    lower_ci = md("**Lower 95% CI**"),
    upper_ci = md("**Upper 95% CI**"),) %>% 
  tab_options(., container.width = 1000, container.height = 1000) %>% 
  tab_options(
    table.font.size = px(10L)) %>% 
  cols_width(everything() ~ px(100)) %>% 
  gtsave(filename = "dlm_problematic2.png")

final_df %>% 
  ungroup %>% 
  gt(auto_align = FALSE) %>% 
  cols_label(
    delamanid_y = md("**Mutation**"),
    count_susceptible = md("**Susceptible Samples with mutation**"),
    count_resistant = md("**Resistant Samples with mutation**"),
    notin_sus = md("**Susceptible Samples without mutation**"),
    notin_resis = md("**Resistant Samples without mutation**"),
    odds_ratio = md("**Odds Ratio**"),
    lower_ci = md("**Lower 95% Confidence Interval**"),
    upper_ci = md("**Upper 95% Confidence Interval**"),) %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()) %>% 
  tab_options(
    column_labels.border.top.color = "black",
    column_labels.border.bottom.color = "black",
    table.border.bottom.width = px(3),
    table.border.bottom.color = "black",
    table_body.hlines.style = "white") %>% 
  tab_options(data_row.padding = px(0.6)) %>% 
  cols_width(delamanid_y ~ px(170),
             everything() ~ px(100)) %>% 
  tab_options(table.font.size = px(12L))





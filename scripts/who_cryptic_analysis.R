#libraries
library(tidyverse)

#### reading in the data ####
cryptic_data <- read_csv("Cryptic_DST_MEASUREMENTS.csv")
who_phenotypic_data <- read_csv("who_phenotypic_data.csv")


## joining the datasets but only subsettig UNIQUEID ?

##### creating a who phenotypic dataset where at least one of the (4/5 last line drugs available) have phenotypic data ####


#who_phenotypic data, subset where phenotype data is not blank ie actually have phenotypic data present 
summary(who_phenotypic_data)

#columns neeeded from dataframe
require_clumns <- c( 'isolate name', 'ena_project', 'ena_sample', 'phenotype Bedaquiline', 'phenotype Clofazimine', 'phenotype Delamanid', 'phenotype Linezolid')

#columns needed to not be blank 
list_of_cols <- c('phenotype Bedaquiline', 'phenotype Clofazimine', 'phenotype Delamanid', 'phenotype Linezolid'  )

#python
who_phenotypic_data[list_of_cols].dropna(thresh=1).head()


#removing only rows where all values are missing 

who_complete_phenotypic <- who_phenotypic_data %>% 
  select( c('isolate name', 'ena_project', 'ena_sample', 'phenotype Bedaquiline', 
            'phenotype Clofazimine', 'phenotype Delamanid', 'phenotype Linezolid')) %>% 
  # keeping rows where at least one of the required columns is not balnk, with filter(df, rowSums(is.na(df)) != ncol(df))
  filter(rowSums(is.na(who_phenotypic_data[c('phenotype Bedaquiline', 
                                             'phenotype Clofazimine', 'phenotype Delamanid', 'phenotype Linezolid')])) != ncol(who_phenotypic_data[c('phenotype Bedaquiline', 
                                                                                                                                                     'phenotype Clofazimine', 'phenotype Delamanid', 'phenotype Linezolid')]))

#checking missing values 
#all have isolate names but 1597 missing ena_project and 859 missing ena_sample -> should these be removed? 
sum(is.na(who_complete_phenotypic$`isolate name`))
sum(is.na(who_complete_phenotypic$ena_sample))
sum(is.na(who_complete_phenotypic$ena_project))

#isolate name is unique in this dataset 
nrow(unique(who_complete_phenotypic['isolate name']))

# question: should I remove the rows that have no ena_sample or ena_project ?? 

# number of rows when all ena_sample and ena_project are present 
nrow(who_complete_phenotypic %>%  
       drop_na(c('ena_sample', 'ena_project')))

who_prac <- who_phenotypic_data %>% 
  select( c('isolate name', 'ena_project', 'ena_sample', 'phenotype Bedaquiline', 
            'phenotype Clofazimine', 'phenotype Delamanid', 'phenotype Linezolid')) %>% 
  # keeping rows where at least one of the required columns is not balnk, with filter(df, rowSums(is.na(df)) != ncol(df))
  filter(rowSums(is.na(who_phenotypic_data[c('phenotype Bedaquiline', 
                                             'phenotype Clofazimine', 'phenotype Delamanid', 'phenotype Linezolid')])) != ncol(who_phenotypic_data[c('phenotype Bedaquiline', 
                                                                                                                                                    'phenotype Clofazimine', 'phenotype Delamanid', 'phenotype Linezolid')])) %>%
  
  drop_na(c('ena_sample','ena_project'))


View(who_prac)

#### joining with cryptic dataset #### 

joint_who_cryptic <- who_complete_phenotypic %>% 
  rename(UNIQUEID = `isolate name`) %>% 
  inner_join(cryptic_data['UNIQUEID'], by = 'UNIQUEID') 

#it seems that all the who data with at least one phenotypic data present are all covered by the cryptic dataset 
nrow(who_complete_phenotypic)
nrow(unique(joint_who_cryptic))


## understanding variety of drugs in cryptic dataset 
cryptic_data$DRUG <- as.factor(cryptic_data$DRUG)
View(cryptic_data)
summary(cryptic_data$DRUG)



df1[!(df1$name %in% df2$name),]

#there are no isolate names that are not also present in the cryptic dataset 
who_complete_phenotypic[!(who_complete_phenotypic$`isolate name` %in% cryptic_data$UNIQUEID),]

#There are 33,458 rows in cryptic dataset, not included in WHO dataset with information on drugs of interest 
cryptic_data[!(cryptic_data$UNIQUEID %in% who_complete_phenotypic$`isolate name` ),] %>% 
  filter(DRUG %in% c("BDQ", "CFZ", "DLM", "LZD")) %>% 
  select(unique('UNIQUEID'))

View(cryptic_data[!(cryptic_data$UNIQUEID %in% who_complete_phenotypic$`isolate name` ),] %>% 
       filter(DRUG %in% c("BDQ", "CFZ", "DLM", "LZD")))




View(cryptic_data %>% 
       filter(DRUG %in% c("BDQ", "CFZ", "DLM", "LZD")))

#isolate name is unique in this dataset 
nrow(unique(cryptic_data['UNIQUEID']))


##### saving the dataframes as csv files ####


write_csv(who_complete_phenotypic, "who_complete_phenotypic.csv")
write_csv(who_prac, "who_phenotypic_complete_cases.csv")


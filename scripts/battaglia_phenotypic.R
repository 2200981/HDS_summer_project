
### trying to read in all sheets from excel files ####
#libraries 
install.packages("openxlsx")
library(openxlsx)


#reading in excel file with multiple sheets 

# specifying the path name
path <- "Desktop/Summer_Project /phenotypic_battaglia2020.xlsx"

# importing the required library
library(openxlsx)

# getting data from sheets
sheets <- openxlsx::getSheetNames(path)
data_frame <- lapply(sheets, openxlsx::read.xlsx, xlsxFile=path)

# assigning names to data frame
names(data_frame) <- sheets

# printing the data
print(data_frame)
class(data_frame)


head(data_frame$`WHO database`$X6)







#### Prep: reading in the files: approach 2 ####

#libraries
library(tidyverse)


# Make a list of all file paths and file names
files <- list.files("battaglia_2020", full.names = TRUE, pattern = ".csv")
filenames <- str_to_lower(str_remove(str_remove(tools::file_path_sans_ext(basename(files)), "_2020"), "battaglia_"))

# Apply the read_csv function to all files and save the output as one big list
data <- lapply(files, readr::read_csv)
names(data) <- filenames

sapply(data, function(x) sapply(x, class))

View(data$complete_papers)



#### data analysis ####

### bedaquiline 
#subset of isolate 1 

bdq_isolate_1 <- data$bdq %>% 
  select('ISOLATE 1', 'DRS Sample selected (Isolate 1)') %>% 
  rename( MIC = 'ISOLATE 1', 
          SampleID  ='DRS Sample selected (Isolate 1)') %>% 
  left_join(data$who, by="SampleID") %>% 
  mutate(drug = 'Bedaquiline')


View(bdq_isolate_1)

#total rows with MIC ie MIC not empty 
sum(!is.na(bdq_isolate_1$MIC))




#subset of isolate 2

bdq_isolate_2 <- data$bdq %>% 
  select('ISOLATE 2', 'DRS Sample selected (Isolate 2)') %>% 
  rename( MIC = 'ISOLATE 2', 
          SampleID  ='DRS Sample selected (Isolate 2)') %>% 
  left_join(data$who, by="SampleID") %>% 
  mutate(drug = 'Bedaquiline')

View(bdq_isolate_2)

#total rows with MIC ie MIC not empty 
sum(!is.na(bdq_isolate_2$MIC))


## delamanid 


dlm_isolate_1 <- data$dlm %>% 
  select('ISOLATE 1', 'DRS Sample selected (Isolate 1)') %>% 
  rename( MIC = 'ISOLATE 1', 
          SampleID  ='DRS Sample selected (Isolate 1)') %>% 
  left_join(data$who, by="SampleID") %>% 
  mutate(drug = 'Delamanid')


View(dlm_isolate_1)

#total rows with MIC ie MIC not empty 
sum(!is.na(dlm_isolate_1$MIC))


#subset of isolate 2

dlm_isolate_2 <- data$dlm %>% 
  select('ISOLATE 2', 'DRS Sample selected (Isolate 2)') %>% 
  rename( MIC = 'ISOLATE 2', 
          SampleID  ='DRS Sample selected (Isolate 2)') %>% 
  left_join(data$who, by="SampleID") %>% 
  mutate(drug = 'Delamanid')

View(dlm_isolate_2)

#total rows with MIC ie MIC not empty 
sum(!is.na(dlm_isolate_1$MIC))



###joining all the diffrent datasets 
final_dataset <- rbind(bdq_isolate_1, bdq_isolate_2, dlm_isolate_1, dlm_isolate_2 )

#now removing any empty cells 
final_dataset <- final_dataset %>% 
  filter(!is.na(MIC))

#total rows with MIC ie MIC not empty 
sum(!is.na(final_dataset$MIC))

# check if total is correct, should be 0 

sum(!is.na(final_dataset$MIC)) - ( sum(!is.na(bdq_isolate_1$MIC)) + sum(!is.na(bdq_isolate_2$MIC)) +  + 
  sum(!is.na(dlm_isolate_1$MIC)) + sum(!is.na(dlm_isolate_2$MIC)) )




#save df as csv file 
write_csv(final_dataset, "battaglia_phenotypic_data.csv")





n_distinct(final_dataset$SampleID)
#looking for dupliacte sample ID
final_dataset %>% count(SampleID) %>% filter(n > 1)




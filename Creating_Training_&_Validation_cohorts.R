#Reading in datasets
baseline <- read.csv("Covariates.csv", header = TRUE)
Scores <- read.table("Scores.txt", header = TRUE)

#Merging protein levels and MS-DSS and eliminating unwanted patients
baseline <- merge(baseline, Scores, by = "patient")

#Deleting unwanted columns
baseline <- baseline %>%
  select(-c("LP_date", "MS_status", "age_squared", "IgG.index", "serum_CRP", "serum_VitD25OH_D2D3", "platform")) %>%
  select(-c("edss", "ageatedss", "dd", "armss", "armss_std", "msdss_std", "msss", "msss_std")) %>%
  select(patient, msdss_last, everything())

#Fixing misspellings in data
baseline$family_history_MS <- factor(gsub("Yes", "yes", baseline$family_history_MS))
baseline$race <- factor(gsub("White", "white", baseline$race))

#Duplicating baseline dataset for msdss dataset
baselinemsdss <- baseline

#Functions
quant_bin<-function(x,quants=c(0,.5,1),...){
  quantiles<-quantile(x,quants)
  quantiles[1]<-quantiles[1]-0.001
  cut(x,breaks=quantiles,...)
}

fac_strip <- function(x){
  if(!(class(x) %in% c("character","factor"))){return(x)}
  if(class(x) %in% c("character","factor")){
    return(as.factor(as.character(x)))
  }
}

# Training/Validation split - setting up the sampling groups (not for MS-DSS)
baseline <- baseline %>%
  mutate(cat_group = interaction(race,sex,family_history_MS)) %>%
  group_by(cat_group) %>%
  mutate(age_group = quant_bin(age_at_LP)) %>%
  ungroup %>%
  mutate(sample_group = as.character(fac_strip(interaction(cat_group, age_group))))

# Training/Validation split - setting up the sampling groups (for MS-DSS)
baselinemsdss <- baselinemsdss %>%
  mutate(cat_group = interaction(race, sex, family_history_MS)) %>%
  group_by(cat_group) %>%
  mutate(msdss_group = quant_bin(msdss_last),
         age_group = quant_bin(age_at_LP)) %>%
  ungroup %>%
  mutate(msdss_group = fac_strip(interaction(cat_group, msdss_group)),
         age_group = fac_strip(interaction(cat_group, age_group)),
         sample_group = as.character(interaction(msdss_group,age_group)))

#Identifying training patients
set.seed(1085)
train_pats <- baseline %>%
  group_by(sample_group) %>%
  sample_frac(0.66) %>%
  ungroup %>%
  .$patient %>%
  as.character

#Identifying training patients for MSDSS
set.seed(98)
train_pats_msdss <- baselinemsdss %>%
  group_by(sample_group) %>%
  sample_frac(0.66) %>%
  ungroup %>%
  .$patient %>%
  as.character

#Creating column for identification of training and validation
baseline <- baseline %>%
  mutate(cohort = ifelse(patient %in% train_pats,"training","validation")) %>%
  select(-cat_group,-age_group,-sample_group) %>%
  select(patient, msdss_last, cohort, everything())

baselinemsdss <- baselinemsdss %>%
  mutate(cohort = ifelse(patient %in% train_pats_msdss,"training","validation")) %>%
  select(-cat_group,-age_group,-msdss_group, -sample_group) %>%
  select(patient, msdss_last, cohort, everything())

write.csv(baseline, "unadjbaseline.csv")         
write.csv(baselinemsdss, "unadjbaselinemsdss.csv")  

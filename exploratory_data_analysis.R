library(haven)
library(tidyverse)
library(survival)
library(survminer)
library(pheatmap)
library(NNMIS)

synth_cancer <- read_sas("/Users/jeffr/Documents/CCRC_19/CPAC Hack4Cancer/synth_cancer_linked_file.sas7bdat")

subset_synth <- synth_cancer %>% 
  gather(cancer_type, diagnosed, -c(TREGION_SYNTH:SEX_SYNTH)) 

cancer_types <- subset_synth %>% 
  #filter(cancer_type == c("COLOREC_SYNTH","BREAST_SYNTH")) %>%
  filter(diagnosed == 1) %>%
  filter(STAGE_SYNTH !=6) %>%
  filter(SURVINT_SYNTH != 99998) %>%
  filter(SURVINT_SYNTH != 0) %>%
  filter(SURVINT_SYNTH != 99999) %>%
  filter(S_DEAD_SYNTH != 9) 

ethnic_subset_binaryminority <- cancer_types %>%
  mutate(DVISMIN_SYNTH = case_when(DVISMIN_SYNTH==14~0, # aboriginals
                                   DVISMIN_SYNTH==13~1, # not visible minority,
                                   DVISMIN_SYNTH<=12~2)) # visible minority)
# Look at aboriginal status
#ethnic_subset_aboriginal <- cancer_types %>%
#  mutate(aboriginal=case_when(DVISMIN_SYNTH==14~0,
#                              DVISMIN_SYNTH<=13~1))

breast <- ethnic_subset %>%
  filter(cancer_type=="BREAST_SYNTH")

colon <- ethnic_subset %>%
  filter(cancer_type=="COLOREC_SYNTH")

# Run a survival analysis on survival interval and death covariates against all other variables , in particular contorlling for ethnic groups (encoded 0 = indigenous, 1 = non-visible, 2 = visible minority)
fit <- survfit(Surv(time=ethnic_subset_binaryminority$SURVINT_SYNTH, event=ethnic_subset_binaryminority$S_DEAD_SYNTH)~ethnic_subset_binaryminority$DVISMIN_SYNTH, data = ethnic_subset_binaryminority)

# Fit the plot (Nancy's is a fancier version with tidier coloring scheme)
ggsurvplot(fit,
           size =1,
           conf.int = FALSE,
           pval = TRUE,
           risk.table = TRUE,
           risk.table.col="strata")


# Clean data and reduce levels for some variables

all_cancer <- 
  subset_synth %>%
  filter(diagnosed == 1) %>%
  filter(STAGE_SYNTH !=9) %>%
  filter(SURVINT_SYNTH != 99998) %>%
  filter(SURVINT_SYNTH != 0) %>%
  filter(SURVINT_SYNTH != 99999) %>%
  filter(S_DEAD_SYNTH != 9) %>%
  mutate(income = case_when (D_LICORATIO_DA_BEF_SYNTH==1 ~ 1,
                             D_LICORATIO_DA_BEF_SYNTH==2 ~ 1,
                             D_LICORATIO_DA_BEF_SYNTH==3 ~ 2,
                             D_LICORATIO_DA_BEF_SYNTH==4 ~ 2,
                             D_LICORATIO_DA_BEF_SYNTH==5 ~ 3,
                             D_LICORATIO_DA_BEF_SYNTH==6 ~ 3,
                             D_LICORATIO_DA_BEF_SYNTH==7 ~ 4,
                             D_LICORATIO_DA_BEF_SYNTH==8 ~ 4,
                             D_LICORATIO_DA_BEF_SYNTH==9 ~ 5,
                             D_LICORATIO_DA_BEF_SYNTH==9 ~ 5)) %>%
  mutate(education= case_when (HCDD_SYNTH==1~0,
                               HCDD_SYNTH %in% c(2,3,4,5,6,7) ~ 1,
                               HCDD_SYNTH %in% c(8,9,10,11,12,13)~2)) %>%
  mutate(marriage=case_when(MARST_SYNTH==1~0,
                            MARST_SYNTH==2~1,
                            TRUE ~ 0)) %>%
  mutate(difficult_activity=case_when(DISABFL_SYNTH==1 ~ 0,
                                      DISABFL_SYNTH %in% c(2,3) ~ 1)) %>%
  drop_na()

# maybe restrict to lung, prostate, breast, colorec
top_cancer_types <- c("BREAST_SYNTH", "LUNG_SYNTH", "COLOREC_SYNTH", "PROSTATE_SYNTH")

TOP_CANCERS <- all_cancer %>%
  filter(cancer_type %in% top_cancer_types) 

# Sanity check, look for encoded NAs
# quick sanity check to look at the factor levels to determine which to drop, but can also use the data dictionary
TOP_CANCERS %>%
  mutate_all(as.factor) %>%
  str()


TOP_CANCERS_IMPUTE <- TOP_CANCERS %>% 
  select(-SURG_1_SYNTH, -SURG_2_SYNTH) %>%
  mutate_at("S_DEAD_SYNTH", list(~replace(., . == 2, 0))) %>%
  mutate_at("STAGE_SYNTH", list(~replace(., . == 6, NA)))


remove <- c("PREGION_SYNTH","AGE_SYNTH", "CENSOR_SYNTH", 
            "AGE_IMM_REVISED_group_SYNTH", "IMMDER_SYNTH", 
            "CITBIR_SYNTH", "YRIM_GROUP_SYNTH","DISABFL_SYNTH", "KID_GROUP_SYNTH", "COWD_SYNTH", "HCDD_SYNTH","D_LICORATIO_DA_BEF_SYNTH",
            "SURVINT_SYNTH", "S_DEAD_SYNTH")

BREAST <-TOP_CANCERS_IMPUTE %>%
  filter(cancer_type=="BREAST_SYNTH") %>%
  select(-cancer_type)


NO_STAGE_BREAST <- TOP_CANCERS_IMPUTE %>%
  select(-STAGE_SYNTH, -diagnosed, -SEX_SYNTH) %>%
  filter(cancer_type == "BREAST_SYNTH") %>%
  select(-cancer_type)

NO_STAGE_BREAST_RMVAR <- NO_STAGE_BREAST %>%
  select(-remove)

#smu(var1 == n & STAGE_SYNTH == 6)/sum(STAGE_SYNTH == 6)
imp.dat_breast <- NNMIS(BREAST$STAGE_SYNTH, 
                        xa=NO_STAGE_BREAST_RMVAR, 
                        xb=NO_STAGE_BREAST_RMVAR, 
                        time=NO_STAGE_BREAST$SURVINT_SYNTH, 
                        event=NO_STAGE_BREAST$S_DEAD_SYNTH, 
                        imputeCT=TRUE, 
                        Seed = 2016,  # Set for RNG
                        mc.cores = 3) # Set number of cores depending on your computating power

data.table::fwrite(imp.dat_breast, "./CPAC_HACK4Cancer/results/", sep = "\t")
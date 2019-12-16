library(sas7bdat)
library(tidyr)
library(dplyr)
library(magrittr)
library(survival)
library(survminer)
library(tableone)
cancer_link<-read.sas7bdat('./Hackthon/synth_cancer_linked_file.sas7bdat')

#no patient id available, so treat each row as a different individual?

#clean cancer type columns, consolidate to single column, cancer-type:
#filter for cancer patients:
cancer_link<-cancer_link%>%
             gather('cancer_type','value',-c(TREGION_SYNTH:SEX_SYNTH))%>%
             filter(value==1)%>%
             mutate(cancer_type=gsub('_SYNTH$','',cancer_type))

#limit to staging data:

#cancer_link_subset<-cancer_link%>%filter(STAGE_SYNTH!=6)

#end point definition: death from cancer/other causes, or end of data collection (not dead)
#remove outcome not known
outcome_unknown<-cancer_link_subset%>%filter(S_DEAD_SYNTH==9)


#remove 51 individuals who have missing values in end point:
nrow(cancer_link_subset[cancer_link_subset$SURVINT_SYNTH==99999,])

#due to data quality, we remove zero survival time and autopsy
autopsy<-cancer_link_subset%>%filter(SURVINT_SYNTH==99998)
zero_survival<-cancer_link_subset%>%filter(SURVINT_SYNTH==0)

all_cancer<-cancer_link%>%
            filter(SURVINT_SYNTH!=99998)%>%
            filter(STAGE_SYNTH!=6)%>%
            filter(SURVINT_SYNTH!=0)%>%
            filter(SURVINT_SYNTH!=99999)%>%
            filter(S_DEAD_SYNTH!=9)

summary(as.factor(all_cancer$cancer_type))

#summarise cancer type proportion:
summary(as.factor(cancer_link$cancer_type))

#define exposure:
all_cancer<-all_cancer%>%
                      mutate(enthnic=case_when(DVISMIN_SYNTH==14~0, #aboriginal
                                               DVISMIN_SYNTH==13~1, #non-visible minority
                                               TRUE ~ 2))           # visible minority
all_cancer$enthnic<-as.factor(all_cancer$enthnic)

all_cancer$aboriginal<-ifelse(all_cancer$DVISMIN_SYNTH==14,1,0)
summary(as.factor(all_cancer$enthnic))


survfit_graph<-function(fit,title){
ggsurvplot(fit,
           size = 1,                 # change line size
          # palette = 
          #   c("#E7B800", "#2E9FDF","#00FF00"),# custom color palettes
           conf.int = FALSE,          # Add confidence interval
           pval = TRUE,              # Add p-value
           risk.table = TRUE,        # Add risk table
           risk.table.col = "strata",# Risk table color by groups
           # Change legend labels
           risk.table.height = 0.25, # Useful to change when you have multiple groups
           ggtheme = theme_bw(),
           title=title)
}


#recode variables:
all_cancer<-all_cancer%>%
  mutate(income = case_when (D_LICORATIO_DA_BEF_SYNTH==1 ~ 1,
                             D_LICORATIO_DA_BEF_SYNTH==2 ~ 1,
                             D_LICORATIO_DA_BEF_SYNTH==3 ~ 2,
                             D_LICORATIO_DA_BEF_SYNTH==4 ~ 2,
                             D_LICORATIO_DA_BEF_SYNTH==5 ~ 3,
                             D_LICORATIO_DA_BEF_SYNTH==6 ~ 3,
                             D_LICORATIO_DA_BEF_SYNTH==7 ~ 4,
                             D_LICORATIO_DA_BEF_SYNTH==8 ~ 4,
                             D_LICORATIO_DA_BEF_SYNTH==9 ~ 5,
                             D_LICORATIO_DA_BEF_SYNTH==9 ~ 5))

all_cancer<-all_cancer%>%
  mutate(education= case_when (HCDD_SYNTH==1~0,
                               HCDD_SYNTH %in% c(2,3,4,5,6,7) ~ 1,
                               HCDD_SYNTH %in% c(8,9,10,11,12,13)~2))

all_cancer<-all_cancer%>%
  mutate(marriage=case_when(MARST_SYNTH==1~0,
                            MARST_SYNTH==2~1,
                            TRUE ~ 0))



all_cancer<-all_cancer%>%
  mutate(difficult_activity=case_when(DISABFL_SYNTH==1 ~ 0,
                                      DISABFL_SYNTH %in% c(2,3) ~ 1))



#create subset dataset (breast and colectoral cancer)
breast<-all_cancer%>%filter(cancer_type=='BREAST')
colec<-all_cancer%>%filter(cancer_type=='COLOREC')

#produce table one:
factorVars<-c('S_DEAD_SYNTH','income','RUINDFG_SYNTH','education','marriage','difficult_activity','cancer_type',
              'SEX_SYNTH','TREGION_SYNTH','STAGE_SYNTH')

vars<-c('SURVINT_SYNTH','AGEDIAG_SYNTH')


CreateTableOne(vars = vars, strata = "enthnic2", data = breast, factorVars = catVars)
CreateCatTable(vars=factorVars,data=breast,strata='enthnic2')


CreateTableOne(vars = vars, strata = "enthnic2", data = colec, factorVars = catVars)
CreateCatTable(vars=factorVars,data=colec,strata='enthnic2')



#assess overall survival rate from all types of cancer, stratified by  aboriginal:
fit_all<-survfit(Surv(time=SURVINT_SYNTH,event=S_DEAD_SYNTH)~ enthnic,data=all_cancer)
survfit_graph(fit_all,'Survival rate for all types of cancer')

#rs<-glm(as.factor(S_DEAD_SYNTH)~ AGEDIAG_SYNTH + STAGE_SYNTH + aboriginal, data =  all_cancer,family='binomial')
#summary(rs)


#breast cancer:
fit_breast<-survfit(Surv(time=SURVINT_SYNTH,event=S_DEAD_SYNTH)~enthnic,data=all_cancer,
             subset=(cancer_type=='BREAST'))

fit_prostate<-survfit(Surv(time=all_cancer$SURVINT_SYNTH,event=all_cancer$S_DEAD_SYNTH)~enthnic,data=all_cancer,
                    subset=(cancer_type=='PROSTATE'))
fit_colorec<-survfit(Surv(time=all_cancer$SURVINT_SYNTH,event=all_cancer$S_DEAD_SYNTH)~income,data=all_cancer,
                    subset=(cancer_type=='COLOREC'))
fit_lung<-survfit(Surv(time=all_cancer$SURVINT_SYNTH,event=all_cancer$S_DEAD_SYNTH)~enthnic,data=all_cancer,
                    subset=(cancer_type=='LUNG'))

p1<-survfit_graph(fit_breast,title='Survival rate for breast cancer')
p2<-survfit_graph(fit_prostate,title='Survival rate for prostate cancer')
p3<-survfit_graph(fit_colorec,title='Survival rate for colorectal cancer')
p4<-survfit_graph(fit_lung,title='survival rate for lung cancer')


arrange_ggsurvplots(list(p1,p2,p3,p4),nrow=2)

#Cox model: (look at risk factors, especially social economic determinants)

#CoxpH model:

#remove stage==9

#adjust for exposure base level:
all_cancer<-all_cancer%>%
                mutate(enthnic2=case_when(DVISMIN_SYNTH==13~0,
                                          DVISMIN_SYNTH==14~1,
                                          TRUE ~2))

#colec<-colec%>%filter(STAGE_SYNTH!=9)
colec$enthnic_group<-as.factor(colec$enthnic2)
cole_cox<- coxph(Surv(SURVINT_SYNTH,outcome) ~ enthnic_group + AGEDIAG_SYNTH + marriage + STAGE_SYNTH+ TREGION_SYNTH,data = colec)
summary(cole_cox)

cox.zph(cole_cox) #model diagnostic

breast$enthnic_group<-as.factor(breast$enthnic_group)
breast_cox<- coxph(Surv(SURVINT_SYNTH,outcome) ~ enthnic_group + AGEDIAG_SYNTH + marriage + TREGION_SYNTH ,data = breast)
summary(breast_cox)

#plot forest plot for covariates HR
ggforest(cole_cox, data = colec)
ggforest(breast_cox,data=breast)



#assess surgery rate among enthnic group in breast and colorectal cancer patients:
#recode SURG_1 and SURG_2 variable
#create new surgery variable defined as having at least surgery once:

all_cancer$SURG_1_SYNTH[all_cancer$SURG_1_SYNTH %in% c(6,9)]<-NA
all_cancer$SURG_2_SYNTH[all_cancer$SURG_2_SYNTH %in% c(6,9)]<-NA

all_cancer$SURG_1_SYNTH[all_cancer$SURG_1_SYNTH ==2]<-0
all_cancer$SURG_2_SYNTH[all_cancer$SURG_2_SYNTH ==2]<-0


all_cancer_surgery<-all_cancer_surgery%>%mutate(surgery=ifelse(SURG_1_SYNTH==1 | SURG_2_SYNTH==1,1,0))

breast_surgery<-all_cancer_surgery%>%filter(cancer_type=='BREAST')%>%filter(STAGE_SYNTH %in% c(1,2,3))
colec_surgery<-all_cancer_surgery%>%filter(cancer_type=='COLOREC')%>%filter(STAGE_SYNTH %in% c(1,2,3))

CreateCatTable(vars=c('surgery'),data=breast_surgery,strata='enthnic2')
CreateCatTable(vars=c('surgery'),data=colec_surgery,strata='enthnic2')


tb<-table(breast_surgery$surgery,breast_surgery$enthnic2)
chisq.test(tb)   
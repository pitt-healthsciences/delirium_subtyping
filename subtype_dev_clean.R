# Kelly Potter            #
# kelly.potter@pitt.edu   #
# Reviewed 2/21/2023      #

setwd("~/Dropbox/4. T32/Projects/Delirium Phenotypes/Analysis/Analysis_v2")

# Load libraries
library(haven)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(dplyr)
library(base)
library(ggplot2)
library(missRanger)
library(corrplot)
library(gmodels)
library(reshape2)
library(table1)
library(survival)
library(survminer)
library(fastDummies)
library(Hmisc)
library(ggrepel)
library(sparklyr)
library(gridExtra)
library(cobalt)
library(magick)
library(scatterplot3d)
library(ggalluvial)
library(Rtsne)
library(rgl)
library(flextable)
library(openxlsx)


my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x, ), digits = 2),
       c("",
         "median (Q1-Q3)" =
           sprintf(paste("%s (",Q1,"- %s)"), MEDIAN,Q3)))
  }

pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- kruskal.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

# Load data
load("/Users/kelly/Dropbox/4. T32/Projects/Delirium Phenotypes/brainmind_deid.Rdata")

assess<-brainmind.assess # mental status assessments
daily<-brainmind.daily # daily data collection
oneobs<-brainmind.oneobs # baseline and hospital outcomes

# Create linkage variable from first study day/assessment delirium was identified
### Assessments
id<-as.character(assess$id)
study.day<-as.character(assess$study.day)
asmt.num<-as.character(assess$asmt.num)

assess$id.studyday <-paste(id,"_",study.day,".",asmt.num)
assess <- assess[order(assess$id.studyday),]

assess.positive<-filter(assess, mental.stat=='Delirious', preserve=TRUE)
assess.positive <- assess.positive[order(assess.positive$id.studyday),]

del_firstpos <- assess.positive[!duplicated(assess.positive$id), ]  # N = 731
del_firstpos$id.studyday<-paste(del_firstpos$id,"_",del_firstpos$study.day)

### Daily data
id<-as.character(daily$id)
study.day<-as.character(daily$study.day)

daily$id.studyday <-paste(daily$id,"_",daily$study.day)

##### One observation data does not have study day - must link by id

# Merge files on id.studyday (daily) or id (oneobs)
brain_cohort<-merge(del_firstpos, daily, by='id.studyday')

oneobs <- rename(oneobs, id.x='id')
brain_cohort<-merge(brain_cohort, oneobs, by='id.x') # Full dataset = brain_cohort

# Generate CV SOFA Score
brain_cohort$cvsofa.icu<-NA 

brain_cohort$cvsofa.icu[brain_cohort['cv.icu'] == 'No hypotension'] <- 0
brain_cohort$cvsofa.icu[brain_cohort['cv.icu'] == 'MAP <70'] <- 1
brain_cohort$cvsofa.icu[brain_cohort['cv.icu'] == 'Dop<5 or Dob (any dose), milrinone, vasopressin (alone)'] <- 2
brain_cohort$cvsofa.icu[brain_cohort['cv.icu'] == 'Dop >5, epi/norepi <0.1'] <- 3
brain_cohort$cvsofa.icu[brain_cohort['cv.icu'] == 'Dop >15, epi/norepi >0.1'] <- 4
brain_cohort$cvsofa.icu[brain_cohort['cv.icu'] == 'Dop >15, epi/norepi <0.1 or balloon pump'] <- 4

# Generate average daily dose of sedatives
brain_cohort$cum.adj.benz <- brain_cohort$cum.benz / brain_cohort$study.day.x
brain_cohort$cum.adj.dex <- brain_cohort$cum.dex / brain_cohort$study.day.x
brain_cohort$cum.adj.op <- brain_cohort$cum.op / brain_cohort$study.day.x
brain_cohort$cum.adj.prop <- brain_cohort$cum.prop / brain_cohort$study.day.x

# Generate cumulative days of sepsis
brain_cohort$cum.any.sepsis <- brain_cohort$cum.sepsis + brain_cohort$cum.sevsepsis

# Build active dataset
brain_cohort <- as_tibble(brain_cohort)

brain_mod_feats <- subset(brain_cohort, select = c('id.x', 'age.enroll', 'bmi', 'charlson.score', 'edu', 'ses.score',
                                                   'cum.sat.int', 'highhr.icu', 'map.int.icu', 'map.low.icu' , 'resp.icu', 
                                                   'sat.int.icu', 'sat.low.icu','sbp.int.icu', 'sbp.low.icu', 'temp.high.icu', 'temp.low.icu', 
                                                   'lowsf.icu', 'cvsofa.icu', 'alt.icu', 'ast.icu', 'bili.icu', 'bun.icu', 'cr.icu', 'glu.low.icu',
                                                   'inr.icu', 'lactate.icu', 'pao2.icu', 'plt.icu', 'sodium.icu', 'tro.icu', 'wbc.icu',
                                                   'fio2.icu', 'cum.vent', 'gcs.icu', 'rass.high.icu', 'rass.low.icu', 
                                                   'cum.adj.benz', 'cum.adj.dex', 'cum.adj.op', 'cum.adj.prop','sirs.num', 'cum.any.sepsis'  ))

sum(is.na(brain_mod_feats))

# Impute missing data
data<-brain_mod_feats

missing_perc<- colMeans(is.na(data))*100

set.seed = 04241991
n.features = 42 # Change to number of features in your data
data.imputed <- missRanger(data[,c(2:(n.features+1))],pmm.k=5,num.trees=1000,seed=set.seed)
data.imputed$id.x <- data$id.x
data.imputed <- data.imputed[,c((n.features+1),(1:n.features))]

head(data.imputed)
sum(is.na(data.imputed))

missing_perc_imputed<-colMeans(is.na(data.imputed))*100

raw_means<-colMeans(data, na.rm = TRUE)
raw_sd<-apply(data, 2, sd, na.rm=TRUE)

imputed_means<-colMeans(data.imputed, na.rm = TRUE)
imputed_sd<-apply(data.imputed, 2, sd, na.rm=TRUE)

compare_impute<-tibble(variable = c('id.x', 'age.enroll', 'bmi', 'charlson.score', 'edu', 'ses.score',
                                    'cum.sat.int', 'highhr.icu', 'map.int.icu', 'map.low.icu' , 'resp.icu', 
                                    'sat.int.icu', 'sat.low.icu','sbp.int.icu', 'sbp.low.icu', 'temp.high.icu', 'temp.low.icu', 
                                    'lowsf.icu', 'cvsofa.icu', 'alt.icu', 'ast.icu', 'bili.icu', 'bun.icu', 'cr.icu', 'glu.low.icu',
                                    'inr.icu', 'lactate.icu', 'pao2.icu', 'plt.icu', 'sodium.icu', 'tro.icu', 'wbc.icu',
                                    'fio2.icu', 'cum.vent', 'gcs.icu', 'rass.high.icu', 'rass.low.icu', 
                                    'cum.adj.benz', 'cum.adj.dex', 'cum.adj.op', 'cum.adj.prop','sirs.num', 'cum.any.sepsis'))

df1<-as_tibble_col(raw_means)
df1<-rename(df1, "raw_mean"="value")

df2<-as_tibble_col(imputed_means)
df2<-rename(df2, "imputed_mean"="value")

df3<-as_tibble_col(raw_sd)
df3<-rename(df3, "raw_sd"="value")

df4<-as_tibble_col(imputed_sd)
df4<-rename(df4, "imputed_sd"="value")

compare_impute$missing_perc<-missing_perc
compare_impute$raw_means<-df1$raw_mean
compare_impute$raw_sd<-df3$raw_sd
compare_impute$imputed_mean<-df2$imputed_mean
compare_impute$imputed_sd<-df4$imputed_sd

compare_impute = compare_impute[-1,]

compare_impute<-as.data.frame(compare_impute)
write.xlsx(compare_impute, 'comp_imputed.xlsx')

# Correlations
corrs_all_data <- data.imputed
corrs_all_data$id.x <- NULL
title_all_data<-'Correlation Plot of All Variables'
corrplot(cor(corrs_all_data), tl.cex=0.5, tl.col='black',
         title=title_all_data, mar=c(3,0,3,0))

corrs_sub1 <- subset(data.imputed, select = c('id.x', 'age.enroll', 'bmi', 'charlson.score', 'edu', 'ses.score',
                                                    'highhr.icu', 'map.int.icu',  'resp.icu', 
                                                     'sat.int.icu',  'temp.high.icu',  
                                                     'lowsf.icu', 'cvsofa.icu', 'alt.icu', 'ast.icu', 'bili.icu', 'bun.icu', 'cr.icu', 'glu.low.icu',
                                                     'inr.icu', 'lactate.icu',  'plt.icu', 'sodium.icu', 'tro.icu', 'wbc.icu',
                                                      'cum.vent',  'rass.low.icu', 
                                                     'cum.adj.benz', 'cum.adj.dex', 'cum.adj.op', 'cum.adj.prop', 'cum.any.sepsis'  ))
corrs_sub1$id.x <- NULL
title_sub1_data<-'Correlation Plot of All Candidate Variables'
corrplot(cor(corrs_sub1), tl.cex=0.5, tl.col='black',
         title=title_sub1_data, mar=c(3,0,3,0))

corrs_sub2<-subset(data.imputed, select = c('id.x', 'age.enroll', 'bmi', 'charlson.score', 'edu', 'ses.score',
                                       'highhr.icu', 'map.int.icu',  'resp.icu', 
                                       'sat.int.icu',  'temp.high.icu',  
                                       'lowsf.icu', 'cvsofa.icu',  'bili.icu', 'cr.icu', 'glu.low.icu',
                                       'inr.icu', 'lactate.icu',  'plt.icu', 'sodium.icu', 'tro.icu', 'wbc.icu',
                                       'cum.vent',  'rass.low.icu', 
                                       'cum.adj.benz', 'cum.adj.dex', 'cum.adj.op', 'cum.adj.prop' ))

corrs_sub2$id.x <- NULL
title_sub2_data<-'Correlation Plot of Final Model Variables'
corrplot(cor(corrs_sub2), tl.cex=0.5, tl.col='black',
         title=title_sub2_data, mar=c(3,0,3,0))

# Build latent class model
mod_indic<-subset(data.imputed, select = c('id.x', 'age.enroll', 'bmi', 'charlson.score', 'edu', 'ses.score',
                                           'highhr.icu', 'map.int.icu',  'resp.icu', 
                                           'sat.int.icu',  'temp.high.icu',  
                                           'lowsf.icu', 'cvsofa.icu',  'bili.icu', 'cr.icu', 'glu.low.icu',
                                           'inr.icu', 'lactate.icu',  'plt.icu', 'sodium.icu', 'tro.icu', 'wbc.icu',
                                           'cum.vent',  'rass.low.icu', 
                                           'cum.adj.benz', 'cum.adj.dex', 'cum.adj.op', 'cum.adj.prop' ))

# Latent class analysis features
features<-c('age.enroll', 'bmi', 'charlson.score', 'edu', 'ses.score',
            'highhr.icu', 'map.int.icu',  'resp.icu', 
            'sat.int.icu',  'temp.high.icu',  
            'lowsf.icu', 'cvsofa.icu',  'bili.icu', 'cr.icu', 'glu.low.icu',
            'inr.icu', 'lactate.icu',  'plt.icu', 'sodium.icu', 'tro.icu', 'wbc.icu',
            'cum.vent',  'rass.low.icu', 
            'cum.adj.benz', 'cum.adj.dex', 'cum.adj.op', 'cum.adj.prop')

mod_indic<-as_tibble(mod_indic)
write.csv(mod_indic,"./mod_indic_imputed.csv")
write_sav(mod_indic, "./mod_indic_imputed.sav")

### Latent class analysis with Latent Gold 6.0 was used to identify latent subtypes of patients with delirium. 
### 1 through 10 hypothesized classes were fit and global fit was examined using BIC, AIC, and entropy R^2. 

### Test models with variables >25% missingness (troponin, lactate, INR, bilirubin) removed and compare to model including all indicators.
mod_fit_0801 <- read_csv('mod_fit_08012022 - mod_fit_c.csv')

mod_fit_0801$Classes<-as.numeric(mod_fit_0801$Classes)

ggplot(data=mod_fit_0801, aes(x=Classes, y=BIC_LL, group=Model_type)) +
  geom_line(aes(linetype=Model_type, color=Model_type))+
  geom_point(aes(color=Model_type))+
  theme(legend.position="top")+
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

mod_fit<-as.data.frame(mod_fit_0801)
mod_fit$LL<-NULL
mod_fit$AIC_LL<-NULL
mod_fit$AIC3_LL<-NULL
mod_fit$Npar<-NULL
mod_fit$Max_BVR<-NULL

mod_fit

## Compare 3 vs 3
all_indic_3_0801<-read_sav('/Users/kelly/Dropbox/4. T32/Projects/Delirium Phenotypes/Analysis/Analysis_v2/all_indic_3_08012022_b.sav')
high_miss_rem_3_0801<-read_sav('/Users/kelly/Dropbox/4. T32/Projects/Delirium Phenotypes/Analysis/Analysis_v2/high_miss_rem_3_08012022_b.sav')

mod_indic$all_indic_clu_3<-all_indic_3_0801$cluster
mod_indic$high_miss_rem_clu_3<-high_miss_rem_3_0801$cluster

mod_indic %>%
  count(high_miss_rem_clu_3)

CrossTable(mod_indic$all_indic_clu_3, mod_indic$high_miss_rem_clu_3, format=c("SAS"))

cluster_crosstabs_3_class<-read_csv('ClusterCrossTabs3Class.csv')
cluster_crosstabs_3_class$three_class_all<-as.factor(cluster_crosstabs_3_class$three_class_all)
cluster_crosstabs_3_class$three_class_miss_rem<-as.factor(cluster_crosstabs_3_class$three_class_miss_rem)

ggplot(data = cluster_crosstabs_3_class,
       aes(axis1 = three_class_all, axis2 = three_class_miss_rem,
           y = frequency)) +
  scale_x_discrete(limits = c("All Indicator Model", "High Missing-Removed Model"), expand = c(.2, .05)) +
  xlab("Class Assignment") +
  geom_alluvium(aes(fill = three_class_all)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  xlab(label = "Delirium Subtypes")+
  ylab(label='Frequency')+
  theme_minimal() +
  ggtitle("Crosstabulations between models")

## Features
mod_indic_z <-subset(mod_indic, select = c('age.enroll', 'bmi', 'charlson.score', 'edu', 'ses.score','highhr.icu', 'map.int.icu',  'resp.icu',  'sat.int.icu',  'temp.high.icu',  
                                           'lowsf.icu', 'cvsofa.icu',  'bili.icu', 'cr.icu', 'glu.low.icu',
                                           'inr.icu', 'lactate.icu',  'plt.icu', 'sodium.icu', 'tro.icu', 'wbc.icu',
                                           'cum.vent',  'rass.low.icu',  'cum.adj.benz', 'cum.adj.dex', 'cum.adj.op', 'cum.adj.prop' ))

mod_indic_z<-scale(mod_indic_z, center=TRUE,scale=TRUE)
mod_indic_z<-as_tibble(mod_indic_z)

colnames(mod_indic_z) <- c('age', 'bmi', 'cci', 'edu', 'ses', 'hr', 'map', 'rr', 'sat',
                           'temp', 'sf', 'cvsof', 'bili', 'cr', 'glu', 'inr', 'lac', 'plt',
                           'na','tro', 'wbc', 'vent', 'rass', 'benz', 'dex', 'op', 'prop')

mod_indic_z$id.x<-mod_indic$id.x

mod_indic_z_all_3<-mod_indic_z
mod_indic_z_all_3$all_indic_clu_3<-mod_indic$all_indic_clu_3

mod_indic_z_high_miss_rem_3<-mod_indic_z
mod_indic_z_high_miss_rem_3$high_miss_rem_clu_3<-mod_indic$high_miss_rem_clu_3

## Compare 4 vs 4
all_indic_4_0801<-read_sav('/Users/kelly/Dropbox/4. T32/Projects/Delirium Phenotypes/Analysis/Analysis_v2/all_indic_4_08012022_b.sav')
high_miss_rem_4_0801<-read_sav('/Users/kelly/Dropbox/4. T32/Projects/Delirium Phenotypes/Analysis/Analysis_v2/high_miss_rem_4_08012022_b.sav')

mod_indic$all_indic_clu_4<-all_indic_4_0801$cluster
mod_indic$high_miss_rem_clu_4<-high_miss_rem_4_0801$cluster

all_indic_4_probs<-subset(all_indic_4_0801, select=c('id.x', 'clu1', 'clu2', 'clu3', 'clu4', 'cluster'))
high_miss_rem_4_probs<-subset(high_miss_rem_4_0801, select=c('id.x', 'clu1', 'clu2', 'clu3', 'clu4', 'cluster'))
compare_reclass<-merge(all_indic_4_probs, high_miss_rem_4_probs, by='id.x')
compare_reclass$change<-all_indic_4_probs$cluster - high_miss_rem_4_probs$cluster
write.csv(compare_reclass, './compare_reclass.csv')

CrossTable(mod_indic$all_indic_clu_4, mod_indic$high_miss_rem_clu_4, format=c("SAS"))

cluster_crosstabs_4_class<-read_csv('ClusterCrossTabs4Class.csv')
cluster_crosstabs_4_class$four_class_all<-as.factor(cluster_crosstabs_4_class$four_class_all)
cluster_crosstabs_4_class$four_class_miss_rem<-as.factor(cluster_crosstabs_4_class$four_class_miss_rem)

ggplot(data = cluster_crosstabs_4_class,
       aes(axis1 = four_class_all, axis2 = four_class_miss_rem,
           y = frequency)) +
  scale_x_discrete(limits = c("All Variable Model", "> 25% Missing-Removed Model"), expand = c(.2, .05)) +
  xlab("Class Assignment") +
  geom_alluvium(aes(fill = four_class_all)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_y_continuous(breaks = c(0, 100, 200, 300, 400, 500, 600, 700, 800))+
  xlab(label = "Class Assignment")+
  ylab(label='Frequency')+
  theme_minimal() +
  ggtitle("Reclassification After Variable Removal")+
  scale_fill_discrete(name = "Four Class Model")

# Compare 3 vs 4 in selected model
CrossTable(mod_indic$all_indic_clu_3, mod_indic$all_indic_clu_4, format=c("SAS"))

crosstabs = read.csv("ClusterCrossTabs.csv")  # read csv file

crosstabs$three_class<-as.factor(crosstabs$three_class)
crosstabs$four_class<-as.factor(crosstabs$four_class)

ggplot(data = crosstabs,
       aes(axis1 = four_class, axis2 = three_class,
           y = frequency)) +
  scale_x_discrete(limits = c("Four Class", "Three Class"), expand = c(.2, .05)) +
  scale_y_continuous(breaks = c(0, 100, 200, 300, 400, 500, 600, 700, 800))+
  xlab("Class Assignment") +
  ylab(label='Frequency')+
  geom_alluvium(aes(fill = four_class)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  ggtitle("Reclassification from Three- to Four-Class Model")+
  scale_fill_discrete(name = "Four Class Model")

# Highly missing removed model
CrossTable(mod_indic$high_miss_rem_clu_3, mod_indic$high_miss_rem_clu_4, format=c("SAS"))

### Cluster link
cluster.link<-as_tibble(mod_indic$id.x)
cluster.link<-rename(cluster.link, "id.x"="value")
cluster.link$cluster_3<-all_indic_3_0801$cluster
cluster.link$cluster_4<-all_indic_4_0801$cluster

write.csv(cluster.link, './cluster.link.csv')

# Examine hospital outcomes
brain_cohort$cluster_4<-mod_indic$all_indic_clu_4

brain_cohort$mort.30day<-brain_cohort$days.deathlast.30
brain_cohort$mort.30day<-as.numeric(brain_cohort$mort.30day)
brain_cohort["mort.30day"][brain_cohort["days.deathlast.30"] < 30] <- 1
brain_cohort["mort.30day"][brain_cohort["days.deathlast.30"] > 29] <- 0

brain_cohort <-brain_cohort %>%
  reorder_levels(cluster_4, order = c('1', '2', '3', '4'))

hosp_outs<-c('dfd.s', 'cfd.s', 'dcfd.s', 'mort.30day')

kw_tests_hosp_4 <- lapply(hosp_outs, function(x) kruskal_test(reformulate("cluster_4", x), data = brain_cohort))
kw_tests_hosp_4

eff_size_hosp_4 <- lapply(hosp_outs, function(x) kruskal_effsize(reformulate("cluster_4", x), data=brain_cohort))
eff_size_hosp_4

### 4 x Delirium-Free Days
res.kruskal <- brain_cohort %>% kruskal_test(dfd.s ~ cluster_4) ### CHANGE DEP VARIABLE ###

pwc <- brain_cohort %>% 
  dunn_test(dfd.s ~ cluster_4, p.adjust.method = 'none') ### CHANGE DEP VARIABLE ###
pwc

pwc <- pwc %>% add_xy_position(x = "cluster_4")
del <- ggboxplot(brain_cohort, x = "cluster_4", y = "dfd.s") + ### CHANGE DEP VARIABLE ###
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc))+ 
  ggtitle("Delirium-Free Days Among Subtypes (Non-Survivors Included)") +
  xlab("Classes") + ylab("Delirium-Free Days")

### 4 x Coma-Free Days
res.kruskal <- brain_cohort %>% kruskal_test(cfd.s ~ cluster_4) ### CHANGE DEP VARIABLE ###

pwc <- brain_cohort %>% 
  dunn_test(cfd.s ~ cluster_4, p.adjust.method = 'none') ### CHANGE DEP VARIABLE ###
pwc

pwc <- pwc %>% add_xy_position(x = "cluster_4")
coma <- ggboxplot(brain_cohort, x = "cluster_4", y = "cfd.s") + ### CHANGE DEP VARIABLE ###
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)) +
  ggtitle("Coma-Free Days Among Subtypes (Non-Survivors Included)") +
  xlab("Classes") + ylab("Coma-Free Days")

### 4 x Delirium/Coma-Free Days
res.kruskal <- brain_cohort %>% kruskal_test(dcfd.s ~ cluster_4) ### CHANGE DEP VARIABLE ###

pwc <- brain_cohort %>% 
  dunn_test(dcfd.s ~ cluster_4, p.adjust.method = 'none') ### CHANGE DEP VARIABLE ###
pwc

pwc <- pwc %>% add_xy_position(x = "cluster_4")
delcoma <- ggboxplot(brain_cohort, x = "cluster_4", y = "dcfd.s") + ### CHANGE DEP VARIABLE ###
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)) + 
  ggtitle("Delirium- or Coma-Free Days Among Subtypes (Non-Survivors Included)") +
  xlab("Classes") + ylab("Delirium/Coma-Free Days")

### 4 x Duration of Delirium
del.res.kruskal <- brain_cohort %>% kruskal_test(del.s ~ cluster_4) ### CHANGE DEP VARIABLE ###

del.pwc <- brain_cohort %>% 
  dunn_test(del.s ~ cluster_4, p.adjust.method = 'none') ### CHANGE DEP VARIABLE ###
del.pwc

del.pwc <- del.pwc %>% add_xy_position(x = "cluster_4")
del.class4 <- ggboxplot(brain_cohort, x = "cluster_4", y = "del.s") + ### CHANGE DEP VARIABLE ###
  stat_pvalue_manual(del.pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(del.res.kruskal, detailed = TRUE),
    caption = get_pwc_label(del.pwc)) + 
  ggtitle("Delirium Among Subtypes (Non-Survivors Included)") +
  xlab("Classes") + ylab("Days of Delirium")

# 30-Day Mortality
surv_object_4 = survival::Surv(time=brain_cohort$days.deathlast.30, event=brain_cohort$mort.30day)
surv_object_4
fit_4=survminer::surv_fit(surv_object_4 ~ cluster_4, data=brain_cohort)
summary(fit_4)
ggsurvplot(fit_4, pval=TRUE,   risk.table = TRUE)

# Hospital outcomes restricted to survivors
brain_cohort_surv<-brain_cohort
brain_cohort_surv<-filter(brain_cohort, died.inhosp.s=='Alive through or withdrew during 30-day study period', preserve=TRUE)
head(brain_cohort_surv)

brain_cohort_surv <-brain_cohort_surv %>%
  reorder_levels(cluster_4, order = c('1', '2', '3', '4'))

# 4 x Days of Delirium (Survivors Only)
res.kruskal <- brain_cohort_surv %>% kruskal_test(del.s ~ cluster_4) ### CHANGE DEP VARIABLE ###
res.kruskal

res.eff<- brain_cohort_surv %>% kruskal_effsize(del.s ~ cluster_4) ### CHANGE DEP VARIABLE ###
res.eff

pwc <- brain_cohort_surv %>% 
  dunn_test(del.s ~ cluster_4, p.adjust.method = 'none') ### CHANGE DEP VARIABLE ###
pwc

pwc <- pwc %>% add_xy_position(x = "cluster_4")
del_surv <- ggboxplot(brain_cohort_surv, x = "cluster_4", y = "del.s") + ### CHANGE DEP VARIABLE ###
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc))+ 
  ggtitle("Days of Delirium Among Subtypes (Survivors Only)") +
  xlab("Subtype") + ylab("Days of Delirium")

# 4 x Days of Coma (Survivors Only)
res.kruskal <- brain_cohort_surv %>% kruskal_test(coma.s ~ cluster_4) ### CHANGE DEP VARIABLE ###
res.kruskal

res.eff<- brain_cohort_surv %>% kruskal_effsize(coma.s ~ cluster_4) ### CHANGE DEP VARIABLE ###
res.eff

pwc <- brain_cohort_surv %>% 
  dunn_test(coma.s ~ cluster_4, p.adjust.method = 'none') ### CHANGE DEP VARIABLE ###
pwc

pwc <- pwc %>% add_xy_position(x = "cluster_4")
coma_surv <- ggboxplot(brain_cohort_surv, x = "cluster_4", y = "coma.s") + ### CHANGE DEP VARIABLE ###
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)) + 
  ggtitle("Days of Coma Among Subtypes (Survivors Only)") +
  xlab("Subtype") + ylab("Days of Coma")

# Impute missing follow-up data
### Identify patients with partial follow-up - 3mo
fu_3mo<-fu_3mo<-filter(brainmind.fu, fu.period=='3 Month')
fu_3mo_ptl<-filter(fu_3mo, status=='Living-Active in the study'| status== 'Living')
fu_3mo_ptl$id.x<-fu_3mo_ptl$id

fu_3mo_ptl$frailty.fu.num<-fu_3mo_ptl$frailty.fu
fu_3mo_ptl$frailty.fu.num<-NA
fu_3mo_ptl$frailty.fu.num<-as.numeric(fu_3mo_ptl$frailty.fu.num)

fu_3mo_ptl["frailty.fu.num"][fu_3mo_ptl["frailty.fu"] == '1. Very Fit'] <- 1
fu_3mo_ptl["frailty.fu.num"][fu_3mo_ptl["frailty.fu"] == '2. Well'] <- 2
fu_3mo_ptl["frailty.fu.num"][fu_3mo_ptl["frailty.fu"] == '3. Well, with treated co-morbid ds'] <- 3
fu_3mo_ptl["frailty.fu.num"][fu_3mo_ptl["frailty.fu"] == '4. Apparently Vulnerable'] <- 4
fu_3mo_ptl["frailty.fu.num"][fu_3mo_ptl["frailty.fu"] == '5. Mildly Frail'] <- 5
fu_3mo_ptl["frailty.fu.num"][fu_3mo_ptl["frailty.fu"] == '6. Moderately Frail'] <- 6
fu_3mo_ptl["frailty.fu.num"][fu_3mo_ptl["frailty.fu"] == '7. Severely Frail'] <- 7

fu_3mo_ptl<-subset(fu_3mo_ptl, select=c('id.x', 'bdi.totscore', 'mmse.rawscore', 'mmse.tscore', 
                                        'trail.a.tscore', 'trail.b.tscore', 'list.learn.tot', 'story.mem.tot', 'rbans.immmemory.tscore',
                                        'fig.copy.tot', 'line.orient.tot', 'rbans.visuo.tscore', 'rbans.sem.score', 
                                        'pic.name.tot', 'rbans.language.tscore', 'rbans.code.score', 'digit.span.tot',
                                        'rbans.attention.tscore', 'list.recall.tot', 'story.recall.tot', 'fig.recall.tot',
                                        'list.recog.tot', 'rbans.delayedmem.tscore', 'rbans.domains.tot', 'rbans.global.score',
                                        'rbans.global.cat', 'rbans.immmemory.cat', 'rbans.visuo.cat', 'rbans.language.cat',
                                        'rbans.attention.cat', 'rbans.delayedmem.cat', 'cogimp', 'cilevel',
                                        'cort.subcort.score', 'faq.totscore', 'adl.totscore', 'sf36.pcs', 'sf36.mcs', 'frailty.fu.num' ))

fu_3mo_ptl<-merge(fu_3mo_ptl, data.imputed, by='id.x')

#### Impute missing follow-up data
colMeans(is.na(fu_3mo_ptl))
data.3mo<-fu_3mo_ptl
data.3mo$id<-NULL
data.3mo$id.x<-NULL

set.seed = 04241991
data.imputed.3mofu <- missRanger(data.3mo,pmm.k=5,num.trees=1000,seed=set.seed)
data.imputed.3mofu$id.x <- data.3mo$id.x

sum(is.na(data.imputed.3mofu))

### Identify patients with partial follow-up - 12mo
fu_12mo<-filter(brainmind.fu, fu.period=='12 Month')
fu_12mo_ptl<-filter(fu_12mo, status=='Living-Active in the study' | status== 'Living')
fu_12mo_ptl$id.x<-fu_12mo_ptl$id

fu_12mo_ptl$frailty.fu.num<-fu_12mo_ptl$frailty.fu
fu_12mo_ptl$frailty.fu.num<-NA
fu_12mo_ptl$frailty.fu.num<-as.numeric(fu_12mo_ptl$frailty.fu.num)

fu_12mo_ptl["frailty.fu.num"][fu_12mo_ptl["frailty.fu"] == '1. Very Fit'] <- 1
fu_12mo_ptl["frailty.fu.num"][fu_12mo_ptl["frailty.fu"] == '2. Well'] <- 2
fu_12mo_ptl["frailty.fu.num"][fu_12mo_ptl["frailty.fu"] == '3. Well, with treated co-morbid ds'] <- 3
fu_12mo_ptl["frailty.fu.num"][fu_12mo_ptl["frailty.fu"] == '4. Apparently Vulnerable'] <- 4
fu_12mo_ptl["frailty.fu.num"][fu_12mo_ptl["frailty.fu"] == '5. Mildly Frail'] <- 5
fu_12mo_ptl["frailty.fu.num"][fu_12mo_ptl["frailty.fu"] == '6. Moderately Frail'] <- 6
fu_12mo_ptl["frailty.fu.num"][fu_12mo_ptl["frailty.fu"] == '7. Severely Frail'] <- 7

fu_12mo_ptl<-subset(fu_12mo_ptl, select=c('id.x', 'bdi.totscore', 'mmse.rawscore', 'mmse.tscore', 
                                          'trail.a.tscore', 'trail.b.tscore', 'list.learn.tot', 'story.mem.tot', 'rbans.immmemory.tscore',
                                          'fig.copy.tot', 'line.orient.tot', 'rbans.visuo.tscore', 'rbans.sem.score', 
                                          'pic.name.tot', 'rbans.language.tscore', 'rbans.code.score', 'digit.span.tot',
                                          'rbans.attention.tscore', 'list.recall.tot', 'story.recall.tot', 'fig.recall.tot',
                                          'list.recog.tot', 'rbans.delayedmem.tscore', 'rbans.domains.tot', 'rbans.global.score',
                                          'rbans.global.cat', 'rbans.immmemory.cat', 'rbans.visuo.cat', 'rbans.language.cat',
                                          'rbans.attention.cat', 'rbans.delayedmem.cat', 'cogimp', 'cilevel',
                                          'cort.subcort.score', 'faq.totscore', 'adl.totscore', 'sf36.pcs', 'sf36.mcs', 'frailty.fu.num' ))

fu_12mo_ptl<-merge(fu_12mo_ptl, data.imputed, by='id.x')

#### Impute missing follow-up data
data.12mo<-fu_12mo_ptl
data.12mo$id<-NULL
data.12mo$id.x<-NULL

set.seed = 04241991
data.imputed.12mofu <- missRanger(data.12mo,pmm.k=5,num.trees=1000,seed=set.seed)
data.imputed.12mofu$id.x <- data.12mo$id.x

head(data.imputed.12mofu)
sum(is.na(data.imputed.12mofu))

# Examine long-term outcomes
### 3 Month
data.imputed.3mofu$id.x<-fu_3mo_ptl$id.x
data.imputed.12mofu$id.x<-fu_12mo_ptl$id.x

ltfu_3mo_cluster<-merge(data.imputed.3mofu, cluster.link, by='id.x')
ltfu_12mo_cluster<-merge(data.imputed.12mofu, cluster.link, by='id.x')

res.kruskal <- ltfu_3mo_cluster %>% kruskal_test(rbans.global.score ~ cluster_4)

pwc <- pwc %>% add_xy_position(x = "cluster_4")
rbans_3mo<- ggboxplot(ltfu_3mo_cluster, x = "cluster_4", y = "rbans.global.score") + ### CHANGE DEP VARIABLE ###
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc), y.position=120)+ 
  ggtitle("RBANS Global Score (Cognition) - 3 Months") +
  xlab("Classes") + ylab("RBANS Global Score")

res.kruskal <- ltfu_3mo_cluster %>% kruskal_test(bdi.totscore ~ cluster_4)

pwc <- pwc %>% add_xy_position(x = "cluster_4")
bdi_3mo <- ggboxplot(ltfu_3mo_cluster, x = "cluster_4", y = "bdi.totscore") + ### CHANGE DEP VARIABLE ###
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)) +
  ggtitle("BDI Total Score (Depression) - 3 Months") +
  xlab("Classes") + ylab("BDI Total Score")

res.kruskal <- ltfu_3mo_cluster %>% kruskal_test(faq.totscore ~ cluster_4)

pwc <- pwc %>% add_xy_position(x = "cluster_4")
faq_3mo <- ggboxplot(ltfu_3mo_cluster, x = "cluster_4", y = "faq.totscore") + ### CHANGE DEP VARIABLE ###
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)) +
  ggtitle("FAQ Score (Physical Ability) - 3 Months") +
  xlab("Classes") + ylab("FAQ Score")

res.kruskal <- ltfu_3mo_cluster %>% kruskal_test(adl.totscore ~ cluster_4)

pwc <- pwc %>% add_xy_position(x = "cluster_4")
katz_3mo <- ggboxplot(ltfu_3mo_cluster, x = "cluster_4", y = "adl.totscore") + ### CHANGE DEP VARIABLE ###
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)) +
  ggtitle("Katz ADL Score (Physical Ability) - 3 Months") +
  xlab("Classes") + ylab("Katz ADL Score")

### 12 Month
res.kruskal <- ltfu_12mo_cluster %>% kruskal_test(rbans.global.score ~ cluster_4)

pwc <- pwc %>% add_xy_position(x = "cluster_4")
rbans_12mo<- ggboxplot(ltfu_12mo_cluster, x = "cluster_4", y = "rbans.global.score") + ### CHANGE DEP VARIABLE ###
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)) + 
  ggtitle("RBANS Global Score (Cognition) - 12 Months") +
  xlab("Classes") + ylab("RBANS Global Score")

res.kruskal <- ltfu_12mo_cluster %>% kruskal_test(bdi.totscore ~ cluster_4)

pwc <- pwc %>% add_xy_position(x = "cluster_4")
bdi_12mo <- ggboxplot(ltfu_12mo_cluster, x = "cluster_4", y = "bdi.totscore") + ### CHANGE DEP VARIABLE ###
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)) +
  ggtitle("BDI Total Score (Depression) - 12 Months") +
  xlab("Classes") + ylab("BDI Total Score")

res.kruskal <- ltfu_12mo_cluster %>% kruskal_test(faq.totscore ~ cluster_4)

pwc <- pwc %>% add_xy_position(x = "cluster_4")
faq_12mo <- ggboxplot(ltfu_12mo_cluster, x = "cluster_4", y = "faq.totscore") + ### CHANGE DEP VARIABLE ###
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)) +
  ggtitle("FAQ Score (Physical Ability) - 12 Months") +
  xlab("Classes") + ylab("FAQ Score")

res.kruskal <- ltfu_12mo_cluster %>% kruskal_test(adl.totscore ~ cluster_4)

pwc <- pwc %>% add_xy_position(x = "cluster_4")
katz_12mo <- ggboxplot(ltfu_12mo_cluster, x = "cluster_4", y = "adl.totscore") + ### CHANGE DEP VARIABLE ###
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)) +
  ggtitle("Katz ADL Score (Physical Ability) - 12 Months") +
  xlab("Classes") + ylab("Katz ADL Score")

# Compare acuity score quartiles to class quartiles for 4-class model
acc.data<-data.imputed

acc.data$cum.avg.sofa.mod<-brain_cohort$cum.avg.sofa.mod
acc.data$sofa.cns<-brain_cohort$sofa.cns
acc.data$sofa.coag<-brain_cohort$sofa.coag
acc.data$sofa.cv<-brain_cohort$sofa.cv
acc.data$sofa.liver<-brain_cohort$sofa.liver
acc.data$sofa.renal<-brain_cohort$sofa.renal
acc.data$sofa.resp<-brain_cohort$sofa.resp

acc.data<-merge(acc.data, brain_cohort[c('id.x', 'apache.aps', 'mod.apache', 'sex.pp', 'race.pp', 'race.wb', 'hisp', 'ins.pp',
                                        'frailty', 'hx.depression', 'hx.alcohol', 'hx.bipolar',
                                        'hx.personality', 'hx.schizo', 'hx.schizo.affect', 'hx.ptsd',
                                        'hx.other', 'hx.psychill', 'hx.concussion', 'hx.stroke', 'hx.isch',
                                        'hx.chf', 'hx.cancer', 'hx.mildliver', 'hx.sevliver', 'hx.liver',
                                        'hx.renal', 'hx.cvd', 'hx.copd', 'hx.diabetes', 'hx.pvd', 'hx.hiv',
                                        'surgery.emergent', 'surgery.elective', 'icu.type', 'admit.dx')], by='id.x')
acc.data$id.x<-NULL

set.seed = 04241991
n.features = 94 # Change to number of features in your data
acc.data.imputed <- missRanger(acc.data,pmm.k=5,num.trees=1000,seed=set.seed)
acc.data.imputed$id.x <- brain_cohort$id.x

head(acc.data.imputed)
sum(is.na(acc.data.imputed))

acc.data.imputed$daily.sofa <- acc.data.imputed$sofa.cns + acc.data.imputed$sofa.coag + acc.data.imputed$sofa.cv + acc.data.imputed$sofa.liver + acc.data.imputed$sofa.renal + acc.data.imputed$sofa.resp

### Assign SOFA quartiles
sofa_quart<-merge(acc.data.imputed, cluster.link[c('id.x', 'cluster_3', 'cluster_4' )], by='id.x')
min(sofa_quart$daily.sofa)
max(sofa_quart$daily.sofa)

sofa_quart$sofa_quart<-1
sofa_quart$sofa_quart[sofa_quart['daily.sofa'] < 7 ] <- 1
sofa_quart$sofa_quart[sofa_quart['daily.sofa'] == 7 ] <- 2
sofa_quart$sofa_quart[sofa_quart['daily.sofa'] == 8 ] <- 2
sofa_quart$sofa_quart[sofa_quart['daily.sofa'] == 9 ] <- 3
sofa_quart$sofa_quart[sofa_quart['daily.sofa'] == 10 ] <- 3
sofa_quart$sofa_quart[sofa_quart['daily.sofa'] == 11 ] <- 3
sofa_quart$sofa_quart[sofa_quart['daily.sofa'] > 11 ] <- 4

sofa_quart$sofa_quart<-as.numeric(sofa_quart$sofa_quart)
sofa_quart$cluster_4<-as.numeric(sofa_quart$cluster_4)

sofa_quart$ins.conds<-sofa_quart$ins.pp
sofa_quart$ins.conds<-as.character(sofa_quart$ins.conds)
sofa_quart$ins.pp<-as.character(sofa_quart$ins.pp)

sofa_quart$ins.conds[sofa_quart$ins.pp == 'Medicaid' ] <- 'Medicaid or dual enrolled'
sofa_quart$ins.conds[sofa_quart$ins.pp == 'Medicare + Medicaid' ] <- 'Medicaid or dual enrolled'
sofa_quart$ins.conds[sofa_quart$ins.pp == 'Medicare' ] <- 'Medicare or medicare + private'
sofa_quart$ins.conds[sofa_quart$ins.pp == 'Medicare + Private' ] <- 'Medicare or medicare + private'
sofa_quart$ins.conds[sofa_quart$ins.pp == 'None' ] <- 'No insurance'
sofa_quart$ins.conds[sofa_quart$ins.pp == 'Private' ] <- 'Private'

write.csv(sofa_quart,"./sofa_quart.csv")
options(stringsAsFactors=F) 
sofa_quart<-read_csv("./sofa_quart.csv")

sofa_quart = select(sofa_quart, -1)

colnames(sofa_quart) <- c('id.x', 'age', 'bmi', 'cci', 'edu', 'ses', 'cum.sat.int', 'high.hr', 'map.int.65',
                          'map.low.icu', 'rr', 'sat.int.90', 'sat.low.icu', 'sbp.int.icu', 'sbp.low.icu', 
                           'high.temp', 'temp.low.icu', 'low.sf', 'cvsofa', 'alt', 'ast', 'bili', 'bun', 'cr', 'glu', 'inr', 'lac', 
                          'pao2.icu', 'plt', 'na', 'tro', 'wbc', 'fio2.icu', 'cumul.vent', 
                          'gcs.icu', 'rass.high.icu', 'low.rass', 'avg.daily.benz', 'avg.daily.dex', 'avg.daily.op', 'avg.daily.prop',
                          'sirs.num', 'cum.any.sepsis', 'cum.avg.sofa.mod', 'sofa.cns', 'sofa.coag', 'sofa.cv', 'sofa.liver','sofa.renal',
                          'sofa.resp', 'apache.aps', 'mod.apache', 'sex.pp', 'race.pp', 'race.wb', 'hisp', 'ins.pp', 
                          'frailty', 'hx.depression', 'hx.alcohol', 'hx.bipolar', 'hx.personality', 'hx.schizo', 
                          'hx.schizo.affect', 'hx.ptsd', 'hx.other', 'hx.psychill', 'hx.concussion', 'hx.stroke', 'hx.isch', 
                          'hx.chf', 'hx.cancer', 'hx.mildliver', 'hx.sevliver', 'hx.liver', 'hx.renal', 'hx.cvd',
                          'hx.copd', 'hx.diabetes', 'hx.pvd', 'hx.hiv', 'surgery.emergent', 'surgery.elective', 'icu.type', 'admit.dx', 
                          'daily.sofa', 'cluster_3', 'cluster_4', 'sofa_quartile', 'ins.conds')

view(sofa_quart)

# Baseline & Vital signs
sofa_quart_age<-ggplot(sofa_quart, aes(y=age)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") +  # second layer
  xlab(label = "")+
  ylab(label = 'Age')

sofa_quart_cci<-ggplot(sofa_quart, aes(y=cci)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") +  # second layer
  xlab(label = "")+
  ylab(label = 'CCI')


sofa_quart_hr<-ggplot(sofa_quart, aes(y=high.hr)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") +  # second layer
  xlab(label = "")+
  ylab(label = 'High HR')


sofa_quart_sf<-ggplot(sofa_quart, aes(y=low.sf)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") +  # second layer
  xlab(label = "")+
  ylab(label = 'Low SF')


sofa_quart_map<-ggplot(sofa_quart, aes(y=map.int.65)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") +  # second layer
  xlab(label = "")+
  ylab(label = 'Map Int < 65')

sofa_quart_spo2<-ggplot(sofa_quart, aes(y=sat.int.90)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") +  # second layer
  xlab(label = "")+
  ylab(label = 'Sat Int < 90')

sofa_quart_cvsofa<-ggplot(sofa_quart, aes(y=cvsofa)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") +  # second layer
  xlab(label = "")+
  ylab(label = 'CV SOFA')

### Laboratory Results
sofa_quart_bili<-ggplot(sofa_quart, aes(y=bili)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") +  # second layer
  xlab(label = "")+
  ylab(label = 'Bilirubin')

sofa_quart_cr<-ggplot(sofa_quart, aes(y=cr)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") +  # second layer
  xlab(label = "")+
  ylab(label = 'Creatinine')

sofa_quart_lactate<-ggplot(sofa_quart, aes(y=lac)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") + # second layer
  xlab(label = "")+
  ylab(label = 'Lactate')


sofa_quart_plt<-ggplot(sofa_quart, aes(y=plt)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") +  # second layer
  xlab(label = "")+
  ylab(label = 'Platelet count')

sofa_quart_tro<-ggplot(sofa_quart, aes(y=tro)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") +  # second layer
  xlab(label = "")+
  ylab(label = 'Troponin')

sofa_quart_inr<-ggplot(sofa_quart, aes(y=inr)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red")  + # second layer
  xlab(label = "")+
  ylab(label = 'INR')

### Ventilation and Sedation
sofa_quart_vent<-ggplot(sofa_quart, aes(y=cumul.vent)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") +  # second layer
  xlab(label = "")+
  ylab(label = 'Cumul Days Vent')

sofa_quart_rass<-ggplot(sofa_quart, aes(y=low.rass)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") +  # second layer
  xlab(label = "")+
  ylab(label = 'Low RASS')

sofa_quart_sepsis<-ggplot(sofa_quart, aes(y=cum.any.sepsis)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") +  # second layer
  xlab(label = "")+
  ylab(label = 'Cumul Days Sepsis')

### Sedatives
sofa_quart_benz<-ggplot(sofa_quart, aes(y=avg.daily.benz)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") +  # second layer
  xlab(label = "")+
  ylab(label = 'Avg Daily Benz')

sofa_quart_dex<-ggplot(sofa_quart, aes(y=avg.daily.dex)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") +  # second layer
  xlab(label = "")+
  ylab(label = 'Avg Daily Dex')

sofa_quart_prop<-ggplot(sofa_quart, aes(y=avg.daily.prop)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red")  + # second layer
  xlab(label = "")+
  ylab(label = 'Avg Daily Prop')

sofa_quart_op<-ggplot(sofa_quart, aes(y=avg.daily.op)) + # basic graphical object
  geom_smooth(aes(x=sofa_quartile), colour="blue") +  # first layer
  geom_smooth(aes(x=cluster_4), colour="red") + # second layer
  xlab(label = "")+
  ylab(label = 'Avg Daily Op')

sofa_plots_base_quart <- grid.arrange(sofa_quart_age, 
                                sofa_quart_cci,
                                sofa_quart_vent,
                                sofa_quart_rass, 
                                sofa_quart_sepsis, nrow = 3,
                                top = "Delirium Subtypes versus SOFA Quartiles Among Variables",
                                bottom = "Red = Delirium Subtype, Blue = SOFA Quartile")

sofa_plots_vs_quart <- grid.arrange(
  sofa_quart_hr,
  sofa_quart_map,
  sofa_quart_cvsofa,
  sofa_quart_spo2,
  sofa_quart_sf,
  nrow = 3,
  top = "Delirium Subtypes versus SOFA Quartiles Among Variables",
  bottom = "Red = Delirium Subtype, Blue = SOFA Quartile")

sofa_plots_labs_quart <- grid.arrange(sofa_quart_bili,
                                sofa_quart_cr,
                                sofa_quart_lactate,
                                sofa_quart_plt,
                                sofa_quart_tro,
                                sofa_quart_inr, nrow = 3,
                                top = "Delirium Subtypes versus SOFA Quartiles Among Variables",
                                bottom = "Red = Delirium Subtype, Blue = SOFA Quartile")

sofa_plots_meds_quart <- grid.arrange(sofa_quart_benz,
                                sofa_quart_dex,
                                sofa_quart_prop,
                                sofa_quart_op, nrow = 2,
                                top = "Delirium Subtypes versus SOFA Quartiles Among Variables",
                                bottom = "Red = Delirium Subtype, Blue = SOFA Quartile")

sofa_plots_quart <- grid.arrange(sofa_quart_age,
             sofa_quart_cci,
             sofa_quart_hr,
             sofa_quart_sf,
             sofa_quart_map,
             sofa_quart_spo2,
             sofa_quart_cvsofa,
             sofa_quart_bili,
             sofa_quart_cr,
             sofa_quart_lactate,
             sofa_quart_plt,
             sofa_quart_tro,
             sofa_quart_inr,
             sofa_quart_vent,
             sofa_quart_rass,
             sofa_quart_sepsis,
             sofa_quart_benz,
             sofa_quart_dex,
             sofa_quart_prop,
             sofa_quart_op, ncol = 3,
             top = "Delirium Subtypes versus SOFA Quartiles Among Variables",
             bottom = "Red = Delirium Subtype, Blue = SOFA Quartile")

# Table 1
library(devtools)
library(table1)

descript <- table1(~ age + sex.pp + race.wb + ins.conds + bmi + cci + edu + ses + icu.type + surgery.elective + surgery.emergent + high.hr + rr +  map.int.65 + sat.int.90 + high.temp +
         low.sf + cvsofa + alt + ast +  bili + na + bun + cr + glu + inr + plt + wbc +
         tro + lac + cumul.vent + low.rass + avg.daily.benz + avg.daily.dex + avg.daily.op + avg.daily.prop + 
         daily.sofa | cluster_4, data=sofa_quart, overall=T, extra.col=list(`P-value`=pvalue), render.continuous=my.render.cont)

st_outs<-table1( ~ dcfd.s + del.s + coma.s + died.inhosp + died.inicu + died.study.30 | cluster_4, data=brain_cohort, overall=T, extra.col=list(`P-value`=pvalue), render.continuous=my.render.cont)

# Compare heat maps of high missing removed model and final model
t1flex(st_outs) %>% 
  save_as_docx(path="del_subtypes_st_outs.docx")

colnames(mod_indic_rob_scale_highmiss_rem)

rob_high_miss<- table1(~age + bmi + cci + edu + ses + hr + map +
                         rr + sat + temp + sf + cvsof + bili +
                         cr + glu + inr + lac + plt + na + 
                         tro + wbc + vent + rass + benz + dex +
                         op + prop + alt + ast + bun + lowsat| high_miss_rem_clu_4, data=mod_indic_rob_scale_highmiss_rem,
                       overall=F, extra.col=list(`P-value`=pvalue), render.continuous=my.render.cont)

t1flex(rob_high_miss) %>% 
  save_as_docx(path="mod_indic_high_miss_rem_descrip.docx")

# Compare knowledge-based subtypes with data-driven subtypes
clin_phen<-subset(brain_cohort, select=c('id.x', 'septic.del', 'hypoxic.del', 'metab.del', 'sed.del', 'other.del', 'days.septic.del', 'days.hypoxic.del', 'days.metab.del', 'days.sed.del', 'days.other.del'))
clin_phen<-merge(clin_phen, cluster.link, by='id.x')

table1(~days.septic.del + days.hypoxic.del + days.metab.del + days.sed.del + days.other.del | cluster_4, data=clin_phen,
       overall=F, extra.col=list(`P-value`=pvalue), render.continuous=my.render.cont)

table1(~septic.del + hypoxic.del + metab.del + sed.del + other.del | cluster_4, data=clin_phen, overall=F,
       extra.col=list(`P-value`=pvalue), render.continuous=my.render.cont)


clin_phen_crosstabs<-read_csv('clin_phen_crosstabs.csv')
clin_phen_crosstabs$clinical_phenotype<-as.factor(clin_phen_crosstabs$clinical_phenotype)
clin_phen_crosstabs$four_class<-as.factor(clin_phen_crosstabs$four_class)

clin_alluvial<- ggplot(data = clin_phen_crosstabs,
       aes(axis1 = clinical_phenotype, axis2 = four_class,
           y = frequency)) +
  scale_x_discrete(limits = c("Clinical Phenotype", "Data-Driven Subtype"), expand = c(.2, .05)) +
  xlab("Class Assignment") +
  geom_alluvium(aes(fill = clinical_phenotype)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  xlab(label = "Clinical vs. Data-Driven Subtypes")+
  ylab(label='Frequency')+
  theme_minimal() +
  ggtitle("Crosstabulations between clinical and data-driven phenotypes")

clin_alluvial<-clin_alluvial + scale_fill_discrete(name = "Clinical Phenotype")
clin_alluvial

clin_phen$septic.del.y<-0
clin_phen$hypoxic.del.y<-0
clin_phen$metab.del.y<-0
clin_phen$sed.del.y<-0
clin_phen$other.del.y<-0

clin_phen["septic.del.y"][clin_phen["septic.del"] == "Septic delirium"] <- 1
clin_phen["hypoxic.del.y"][clin_phen["hypoxic.del"] == "Hypoxic delirium"] <- 1
clin_phen["metab.del.y"][clin_phen["metab.del"] == "Metabolic delirium"] <- 1
clin_phen["sed.del.y"][clin_phen["sed.del"] == "Sedative delirium"] <- 1
clin_phen["other.del.y"][clin_phen["other.del"] == "Other delirium"] <- 1

clin_phen$add.del<- clin_phen$septic.del.y + clin_phen$hypoxic.del.y + clin_phen$metab.del.y + clin_phen$sed.del.y + clin_phen$other.del.y

adj_clin_phen_crosstabs<-read_csv('adj_clin_phen_crosstabs_freq.csv')
adj_clin_phen_crosstabs$clinical_phenotype<-as.factor(adj_clin_phen_crosstabs$clinical_phenotype)
adj_clin_phen_crosstabs$four_class<-as.factor(adj_clin_phen_crosstabs$four_class)

clin_alluvial_adj<- ggplot(data = adj_clin_phen_crosstabs,
                       aes(axis1 = clinical_phenotype, axis2 = four_class,
                           y = frequency)) +
  scale_x_discrete(limits = c("Clinical Phenotype", "Data-Driven Subtype"), expand = c(.2, .05)) +
  xlab("Class Assignment") +
  geom_alluvium(aes(fill = clinical_phenotype)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_y_continuous(breaks = c(0, 100, 200, 300, 400, 500, 600, 700, 800))+
  xlab(label = "Clinical Phenotypes versus Data-Driven Subtypes")+
  ylab(label='Frequency')+
  theme_minimal() +
  theme(legend.position="top")+
  scale_fill_manual(values=c("#0F2080", "#85C0F9", "#A95AA0", "#F5793A", "#44AA99" ), name='Clinical Phenotype') + 
  ggtitle("Crosstabulations Between Clinical Phenotypes and Data-Driven Subtypes")

#clin_alluvial_adj<-clin_alluvial_adj + scale_fill_discrete(name = "Data-Driven Subtype")+theme(legend.position="top")
clin_alluvial_adj

clin_phen$sofa_quart<-sofa_quart$sofa_quartile
clin_phen$rass.high.icu<-data.imputed$rass.high.icu
clin_phen$rass.low.icu<-data.imputed$rass.low.icu

clin_phen$sofa_quartile_str<-''
clin_phen["sofa_quartile_str"][clin_phen["sofa_quart"] == 1] <- 'Q1'
clin_phen["sofa_quartile_str"][clin_phen["sofa_quart"] == 2] <- 'Q2'
clin_phen["sofa_quartile_str"][clin_phen["sofa_quart"] == 3] <- 'Q3'
clin_phen["sofa_quartile_str"][clin_phen["sofa_quart"] == 4] <- 'Q4'

clin_phen$hyperactive.del<-0
clin_phen$hypoactive.del<-0
clin_phen$mixed.del<-0

clin_phen$psychomotor<-''

clin_phen["hyperactive.del"][clin_phen["rass.high.icu"] > 0] <- 1
clin_phen["hypoactive.del"][clin_phen["rass.low.icu"] < 1] <- 1
clin_phen$mixed.del<-clin_phen$hyperactive.del+clin_phen$hypoactive.del
clin_phen["psychomotor"][clin_phen["hyperactive.del"] == 1] <- 'Hyperactive'
clin_phen["psychomotor"][clin_phen["hypoactive.del"] == 1] <- 'Hypoactive'
clin_phen["psychomotor"][clin_phen["mixed.del"] == 2] <- 'Mixed'
clin_phen["psychomotor"][clin_phen["mixed.del"] == 2] <- 'Mixed'
clin_phen["psychomotor"][clin_phen["psychomotor"] == ''] <- 'Mixed'

write.csv(clin_phen, './clin_phen.csv')

adj_clin_phen_crosstabs_all<-read_csv('clin_phen.csv')
adj_clin_phen_crosstabs$clinical_phenotype<-as.factor(adj_clin_phen_crosstabs$clinical_phenotype)
adj_clin_phen_crosstabs$four_class<-as.factor(adj_clin_phen_crosstabs$four_class)

combinations<-crossing(var1 = c('Class 1', 'Class 2', 'Class 3', 'Class'), var2 = c('SOFA 1', 'SOFA 2', 'SOFA 3', 'SOFA 4'), var3 = c('Septic Delirium', 'Hypoxic Delirium', 'Sedative-Associated Delirium', 'Metabolic Delirium', 'Unclassified Delirium'), var4 = c('Hyperactive Delirium', 'Hypoactive Delirium', 'Mixed Delirium'))

write.csv(combinations, './combinations.csv')

psychomotor_crosstabs<-read_csv('psychomotor_crosstabs_freq.csv')
psychomotor_crosstabs$psychomotor<-as.factor(psychomotor_crosstabs$psychomotor)
psychomotor_crosstabs$four_class<-as.factor(psychomotor_crosstabs$four_class)

psychomotor_alluvial<- ggplot(data = psychomotor_crosstabs,
                           aes(axis1 = psychomotor, axis2 = four_class,
                               y = frequency)) +
  scale_x_discrete(limits = c("Psychomotor Subclass", "Data-Driven Subtype"), expand = c(.2, .05)) +
  xlab("Class Assignment") +
  geom_alluvium(aes(fill = psychomotor)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_y_continuous(breaks = c(0, 100, 200, 300, 400, 500, 600, 700, 800))+
  xlab(label = "Psychomotor Subclass versus Data-Driven Subtypes")+
  ylab(label='Frequency')+
  scale_fill_manual(values=c("firebrick3", "royalblue", "purple4" ), name='Psychomotor Subclass') + 
  theme_minimal() +
  ggtitle("Crosstabulations Between Psychomotor Subclasses and Data-Driven Subtypes")+
  theme(legend.position="top")

psychomotor_alluvial

table1(~ psychomotor | cluster_4, data=clin_phen, overall=F,
       extra.col=list(`P-value`=pvalue), render.continuous=my.render.cont)

sofa_quart_crosstabs<-read_csv('sofa_quart_crosstabs_freq.csv')
sofa_quart_crosstabs$sofa_quartile<-as.factor(sofa_quart_crosstabs$sofa_quartile)
sofa_quart_crosstabs$four_class<-as.factor(sofa_quart_crosstabs$four_class)

sofa_quart_alluvial<- ggplot(data = sofa_quart_crosstabs,
                              aes(axis1 = sofa_quartile, axis2 = four_class,
                                  y = frequency)) +
  scale_x_discrete(limits = c("SOFA Quartile", "Data-Driven Subtype"), expand = c(.2, .05)) +
  xlab("Class Assignment") +
  geom_alluvium(aes(fill = sofa_quartile)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_y_continuous(breaks = c(0, 100, 200, 300, 400, 500, 600, 700, 800))+
  xlab(label = "SOFA Quartile versus Data-Driven Subtypes")+
  ylab(label='Frequency')+
  scale_fill_manual(values=c("gray47", "pink", "indianred", 'firebrick4' ), name='SOFA Quartile ') + 
  theme_minimal() +
  ggtitle("Crosstabulations Between SOFA Quartiles and Data-Driven Subtypes")+
  theme(legend.position="top")

sofa_quart_alluvial

table1(~  sofa_quartile_str | cluster_4, data=clin_phen, overall=F,
       extra.col=list(`P-value`=pvalue), render.continuous=my.render.cont)

# Class Membership Probabilities
clu4_pp_1<-filter(all_indic_4_0801, cluster==1, preserve=TRUE)
clu4_pp_2<-filter(all_indic_4_0801, cluster==2, preserve=TRUE)
clu4_pp_3<-filter(all_indic_4_0801, cluster==3, preserve=TRUE)
clu4_pp_4<-filter(all_indic_4_0801, cluster==4, preserve=TRUE)

pp_4_clu_1<- ggplot(clu4_pp_1, aes(x=clu1)) + 
  geom_histogram(binwidth=0.1, fill='white', colour='black')+
  xlab(label = "Probability")+
  scale_x_continuous(breaks = c(0, 0.2, .04, 0.6, 0.8, 1.0))+
  ylab(label='Frequency')+
  ggtitle(label = "Subtype 1")

pp_4_clu_2<- ggplot(clu4_pp_2, aes(x=clu2)) + 
  geom_histogram(binwidth=0.1, fill='white', colour='black')+
  scale_x_continuous(breaks = c(0, 0.2, .04, 0.6, 0.8, 1.0))+
  xlab(label = "Probability")+
  ylab(label='Frequency')+
  ggtitle(label = "Subtype 2")

pp_4_clu_3<- ggplot(clu4_pp_3, aes(x=clu3)) + 
  geom_histogram(binwidth=0.1, fill='white', colour='black')+
  scale_x_continuous(breaks = c(0, 0.2, .04, 0.6, 0.8, 1.0))+
  xlab(label = "Probability")+
  ylab(label='Frequency')+
  ggtitle(label = "Subtype 3")

pp_4_clu_4<- ggplot(clu4_pp_4, aes(x=clu4)) + 
  geom_histogram(binwidth=0.1, fill='white', colour='black')+
  scale_x_continuous(breaks = c(0, 0.2, .04, 0.6, 0.8, 1.0))+
  xlab(label = "Probability")+
  ylab(label='Frequency')+
  ggtitle(label = "Subtype 4")

grid.arrange(pp_4_clu_1, pp_4_clu_2, pp_4_clu_3,pp_4_clu_4, nrow = 2)

# tSNE
tsne.data<-mod_indic
tsne.data$id.x<-NULL
tsne.data$all_indic_clu_3<-NULL
tsne.data$all_indic_clu_4<-NULL
tsne.data$high_miss_rem_clu_3<-NULL
tsne.data$high_miss_rem_clu_4<-NULL

tsne.data
tsne.data.z=scale(tsne.data,T,T)

#Run PCA
tsne <- Rtsne(tsne.data.z, perplexity=20, check_duplicates = FALSE, dims=3)

load1<- as.matrix(tsne$Y[,1])
load2<- as.matrix(tsne$Y[,2])
load3<- as.matrix(tsne$Y[,3])

# combine load1 and load2 and clustering assignment to be a new data for plotting
plot1.d<- data.frame(load1=load1,load2=load2,load3=load3,id=1:nrow(tsne.data))
plot1.d$id.x<-mod_indic$id.x
plot1.d$cluster_4<-cluster.link$cluster_4

#Plot of All 4 Clusters, Not Separating for Core

ggplot(plot1.d) + 
  geom_point(aes(x = load1, y = load2, color = as.factor(cluster_4))) + 
  scale_color_manual(values=c("#0F2080", "#85C0F9", "#A95AA0", "#F5793A")) + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
 ) + coord_fixed() +
  labs(color="Delirium Subtype")+
  xlab(label = "t-SNE Loading 1")+
  ylab(label='t-SNE Loading 2')
  

ggplot(plot1.d) + 
  geom_point(aes(x = load1, y = load3, color = as.factor(cluster_4))) + 
  scale_color_manual(values=c("#0F2080", "#85C0F9", "#A95AA0", "#F5793A")) + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) + coord_fixed() +
  labs(color="Delirium Subtype")+
  xlab(label = "t-SNE Loading 1")+
  ylab(label='t-SNE Loading 3')


ggplot(plot1.d) + 
  geom_point(aes(x = load2, y = load3, color = as.factor(cluster_4))) + 
  scale_color_manual(values=c("#0F2080", "#85C0F9", "#A95AA0", "#F5793A")) + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) + coord_fixed() +
  labs(color="Delirium Subtype")+
  xlab(label = "t-SNE Loading 2")+
  ylab(label='t-SNE Loading 3')

mycolors <- c('royalblue1', 'coral1', 'goldenrod1', 'lightblue1')
mycolors<- c('#0F2080', '#85C0F9',  '#A95AA1', '#F5793A')
plot1.d$color <- mycolors[ as.numeric(plot1.d$cluster_4) ]

#plot3d( 
  #x=plot1.d$load1, y=plot1.d$load2, z=plot1.d$load3, 
  #col = plot1.d$color, 
  #type = 's', 
  #radius = .5,
  #xlab="Load 1", ylab="Load 2", zlab="Load 3")
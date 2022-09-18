#######4.1 GMB was also applicable to multiple cancers with proper gene mutation sets
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/08_GMB_riskscore"
# 设置工作目录
setwd(work_dir)

##########
## Risk model
install.packages("rms")
library(rms)
library(survival)
clinical<-read.csv("01_clinical.csv",row.names = 1,check.names=FALSE)

library(survminer)
#library(ggplot2)
fit<- survfit(Surv(Overall_survival,Survival_status) ~ Risk_group, data = clinical)
surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           break.x.by=24,
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("High_risk", "Low_risk"),
           palette="jto", ylab = "Overall survival" )


###GMB
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/08_GMB_riskscore"
# 设置工作目录
setwd(work_dir)

clinical <- read.csv("01_clinical.csv",row.names = 1,check.names=FALSE)
Gene <- read.csv("01_data_mutations.csv")
gene <- read.csv("01_gene_select.csv")
attach(Gene)

###挑选特定基因所有突变样本
D <- data.frame()
for (i in (1:11)) {
  g <-gene[i,1]
  r<-which(Hugo_Symbol==g)
  Gene_r<-Gene[r,]
  D <- rbind(D,Gene_r)
}

###每个样本含有特定基因突变数目
Da <- data.frame()
for (i in (1:1572)) {
  samp <-rownames(clinical)[i]
  r<-which(D$Tumor_Sample_Barcode==samp)
  Gene_r<-D[r,]
  n=length(r)
  
  Result <-data.frame(n)
  Da <- rbind(Da, Result)
}
write.csv(Da,file = "09_mutation_number.csv")

#####surv-cutpoint OS决定最佳分界值
clinical <- read.csv("09_clinical.csv",row.names = 1,check.names=FALSE)
library(survminer)
res.cut <- surv_cutpoint(clinical,time="Overall_survival",event = "Survival_status",
                         variables = "n")
summary(res.cut)   
plot(res.cut,"n",palette="npg")
res.cat <- surv_categorize(res.cut)

library("survival")
fit<- survfit(Surv(Overall_survival, Survival_status) ~ n, data = res.cat)
surv_pvalue(fit)$pval.txt
clinical$cut<-  res.cat$n
write.csv(clinical,file = "09_clinical.csv")

ggsurvplot(fit, data = res.cat,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           break.x.by=24,
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("Gene_mutation_burden ≥ 1", "Gene_mutation_burden < 1"),
           palette="jto", ylab = "Overall survival" )


##GMB risk model contrast
library(survival)
#Kaplan-Meier 分析，详情 ?survfit
clinical <- read.csv("09_clinical.csv",row.names = 1,check.names=FALSE)
attach(clinical)
KM <- survfit(Surv(Overall_survival,Survival_status) ~ Risk_Mut, data =clinical, type = 'kaplan-meier', conf.type = 'log')
KM
#####
#绘制生存曲线，反映了尚在世的患者数量比例和时间的关系
library(survminer)
#生存曲线，详情 ?ggsurvplot
ggsurvplot(KM, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           palette="jto", ylab = "Overall survival" )
#对数秩检验，详情 ?survdiff
survdiff(Surv(Overall_survival,Survival_status) ~ Risk_Mut, data = clinical)
restest <- pairwise_survdiff(Surv(Overall_survival,Survival_status) ~ Risk_Mut,
                             data = clinical)
restest



#######4.2 GMB was better than risk model in NSCLC
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/09_GMB_riskscore_NSCLC"
# 设置工作目录
setwd(work_dir)

##risk model OS
library(rms)
library(survival)
clinical<-read.csv("01_clinical.csv",row.names = 1,check.names=FALSE)

table(clinical$Cancer_type)
r<-which(clinical$Cancer_type=="Non-small cell lung cancer")
clinical_NSCLC<-clinical[r,]

library(survminer)
#library(ggplot2)
fit<- survfit(Surv(Overall_survival,Survival_status) ~ Risk_group, data = clinical_NSCLC)
surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = clinical_NSCLC,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           break.x.by=24,
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("High_risk", "Low_risk"),
           palette="jto", ylab = "Overall survival" )


###GMB OS
clinical <- read.csv("01_clinical.csv",row.names = 1,check.names=FALSE)
Gene <- read.csv("01_data_mutations.csv")
gene <- read.csv("01_gene_select.csv")
attach(Gene)

##肺癌中突变基因
D <- data.frame()
for (i in (1:350)) {##
  samp <- rownames(clinical_NSCLC)[i]##
  r<-which(Gene$Tumor_Sample_Barcode==samp)
  Gene_r<-Gene[r,]
  D <- rbind(D,Gene_r)
}

###挑选特定基因所有突变样本
Da <- data.frame()
for (i in (1:11)) {
  g <-gene[i,1]
  r<-which(D$Hugo_Symbol==g)
  Gene_r<-D[r,]
  Da <- rbind(Da,Gene_r)
}

###每个样本含有特定基因突变数目
Dao <- data.frame()
for (i in (1:350)) {
  samp <-rownames(clinical_NSCLC)[i]
  r<-which(Da$Tumor_Sample_Barcode==samp)
  Gene_r<-Da[r,]
  n=length(r)
  
  Result <-data.frame(n)
  Dao <- rbind(Dao, Result)
}
write.csv(Dao,file = "09_mutation_number.csv")
write.csv(clinical_NSCLC,file = "09_clinical_NSCLC.csv")

#####surv-cutpoint OS决定最佳分界值
clinical <- read.csv("09_clinical_NSCLC.csv",row.names = 1,check.names=FALSE)
library(survminer)
res.cut <- surv_cutpoint(clinical,time="Overall_survival",event = "Survival_status",
                         variables = "n")
summary(res.cut)   
plot(res.cut,"n",palette="npg")

res.cat <- surv_categorize(res.cut)
#head(res.cat)

library("survival")
fit<- survfit(Surv(Overall_survival, Survival_status) ~ n, data = res.cat)
surv_pvalue(fit)$pval.txt
clinical$cut<-  res.cat$n
write.csv(clinical,file = "09_clinical_NSCLC.csv")

ggsurvplot(fit, data = res.cat,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           break.x.by=24,
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("GMB_H", "GMB_L"),
           palette="jto", ylab = "Overall survival" )


####risk model VS. GMB 
library(survival)
#Kaplan-Meier 分析，详情 ?survfit
clinical <- read.csv("09_clinical_NSCLC_time0.csv",row.names = 1,check.names=FALSE)
attach(clinical)
KM <- survfit(Surv(Overall_survival,Survival_status) ~ Risk_Mut, data =clinical, type = 'kaplan-meier', conf.type = 'log')
KM
#####
#绘制生存曲线，反映了尚在世的患者数量比例和时间的关系
library(survminer)
#生存曲线，详情 ?ggsurvplot
ggsurvplot(KM, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           #legend.labs = c("TMB < 10", "TMB > 10"),
           palette="jto", ylab = "Overall survival" )

#对数秩检验，详情 ?survdiff
survdiff(Surv(Overall_survival,Survival_status) ~ Risk_Mut, data = clinical)
restest <- pairwise_survdiff(Surv(Overall_survival,Survival_status) ~ Risk_Mut,
                             data = clinical)
restest

##C-index
library(survcomp)
f0 <- coxph(Surv(Overall_survival,Survival_status)~Risk_group,data = clinical)
f1 <- coxph(Surv(Overall_survival,Survival_status)~cut,data = clinical)

ci4 <- concordance.index(
  predict(f1), surv.time=Overall_survival, surv.event=Survival_status,
  method="noether"
)
ci4 <- concordance.index(
  predict(f0), surv.time=Overall_survival, surv.event=Survival_status,
  method="noether"
)
ci4$c.index
ci4$lower
ci4$upper









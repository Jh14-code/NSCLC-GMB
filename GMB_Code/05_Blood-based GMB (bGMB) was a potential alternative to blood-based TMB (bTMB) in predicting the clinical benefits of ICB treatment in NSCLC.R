#######5.1 Blood-based GMB (bGMB) predicting the clinical benefits of ICB treatment in NSCLC
##POPLAR GMB
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/15_bGMB"
# 设置工作目录
setwd(work_dir)


clinical <- read.csv("15_POPLAR_clinical.csv",row.names = 1,check.names=FALSE)
Gene <- read.csv("15_POPLAR_Mut.csv")
gene <- read.csv("03_gene_select_NSCLC.csv")
attach(Gene)

###挑选特定基因所有突变样本
D <- data.frame()
for (i in (1:41)) {
  g <-gene[i,2]
  r<-which(Gene$gene_name==g)
  Gene_r<-Gene[r,]
  D <- rbind(D,Gene_r)
}

###每个样本含有特定基因突变数目
Da <- data.frame()
for (i in (1:105)) {
  samp <-rownames(clinical)[i]
  r<-which(D$PtID==samp)
  Gene_r<-D[r,]
  n=length(r)
  
  Result <-data.frame(n)
  Da <- rbind(Da, Result)
}
write.csv(Da,file = "15_POPLAR_clinical_number.csv")

#####surv-cutpoint OS决定最佳分界值
rm(list=ls())
clinical <- read.csv("16_POPLAR_clinical.csv",row.names = 1,check.names=FALSE)
library(survminer)
res.cut <- surv_cutpoint(clinical,time="OS",event = "OS_status",
                         variables = "n")
summary(res.cut)   
plot(res.cut,"n",palette="npg")

res.cat <- surv_categorize(res.cut)
#head(res.cat)

library("survival")
fit<- survfit(Surv(OS, OS_status) ~ n, data = res.cat)
surv_pvalue(fit)$pval.txt
clinical$cut<-  res.cat$n
write.csv(clinical,file = "16_POPLAR_clinical.csv")

ggsurvplot(fit, data = res.cat,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("GMB-H", "GMB-L"),
           palette="jto", ylab = "Overall survival" )

###添加HR ,CI ,P
p3 <- ggsurvplot(fit, data = res.cat,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("GMB-H", "GMB-L"),
                 palette="jto", ylab = "Overall survival" )

res_cox<-coxph(Surv(OS, OS_status) ~ Cut, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3

##OAK GMB
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/15_bGMB"
# 设置工作目录
setwd(work_dir)

clinical <- read.csv("15_OAK_clinical.csv",row.names = 1,check.names=FALSE)
Gene <- read.csv("15_OAK_Mut.csv")
gene <- read.csv("03_gene_select_NSCLC.csv")
attach(Gene)

###挑选特定基因所有突变样本
D <- data.frame()
for (i in (1:41)) {
  g <-gene[i,2]
  r<-which(Gene$gene_name==g)
  Gene_r<-Gene[r,]
  D <- rbind(D,Gene_r)
}

###每个样本含有特定基因突变数目
Da <- data.frame()
for (i in (1:324)) {
  samp <-rownames(clinical)[i]
  r<-which(D$PtID==samp)
  Gene_r<-D[r,]
  n=length(r)
  
  Result <-data.frame(n)
  Da <- rbind(Da, Result)
}
write.csv(Da,file = "15_OAK_clinical_number.csv")

#####
clinical <- read.csv("16_OAK_clinical.csv",row.names = 1,check.names=FALSE)
library(survminer)
library("survival")
fit<- survfit(Surv(OS, OS_status) ~ cut, data = clinical)
surv_pvalue(fit)$pval.txt

ggsurvplot(fit, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("GMB-H", "GMB-L"),
           palette="jto", ylab = "Overall survival" )

###添加HR ,CI ,P
p3 <- ggsurvplot(fit, data = clinical,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("GMB-H", "GMB-L"),
                 palette="jto", ylab = "Overall survival" )

res_cox<-coxph(Surv(OS, OS_status) ~ Cut, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



#######5.2 Blood-based TMB (bTMB) predicting the clinical benefits of ICB treatment in NSCLC
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/19_SPSS"
# 设置工作目录
setwd(work_dir)

#####POPLAR surv-cutpoint OS决定最佳分界值
clinical <- read.csv("16_POPLAR_clinical.csv",row.names = 1,check.names=FALSE)
library(survminer)
res.cut <- surv_cutpoint(clinical,time="OS",event = "OS_status",
                         variables = "btmb")
summary(res.cut)   
plot(res.cut,"btmb",palette="npg")
res.cat <- surv_categorize(res.cut)

library("survival")
fit<- survfit(Surv(OS, OS_status) ~ btmb, data = res.cat)
surv_pvalue(fit)$pval.txt
clinical$cut_bTMB<-  res.cat$btmb
write.csv(clinical,file = "19_POPLAR_clinical.csv")

ggsurvplot(fit, data = res.cat,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("TMB-H", "TMB-L"),
           palette="jto", ylab = "Overall survival" )
###添加HR ,CI ,P
p3 <- ggsurvplot(fit, data = res.cat,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("TMB-H", "TMB-L"),
                 palette="jto", ylab = "Overall survival" )
clinical <- read.csv("19_POPLAR_clinical.csv",row.names = 1,check.names=FALSE)
res_cox<-coxph(Surv(OS, OS_status) ~ Cut_bTMB, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3


#####OAK validation
rm(list=ls())
clinical <- read.csv("19_OAK_clinical.csv",row.names = 1,check.names=FALSE)
library(survminer)
library("survival")
fit<- survfit(Surv(OS, OS_status) ~ cut_btmb, data = clinical)
surv_pvalue(fit)$pval.txt

ggsurvplot(fit, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("TMB-H", "TMB-L"),
           palette="jto", ylab = "Overall survival" )

###添加HR ,CI ,P
p3 <- ggsurvplot(fit, data = clinical,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("TMB-H", "TMB-L"),
                 palette="jto", ylab = "Overall survival" )

res_cox<-coxph(Surv(OS, OS_status) ~ Cut_btmb, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



#######5.3 clinical factors associate with OS
###RACE OS POPLAR
rm(list=ls())
clinical <- read.csv("16_POPLAR_clinical.csv")
library(survminer)
library("survival")
fit<- survfit(Surv(OS, OS_status) ~ race2, data = clinical)
surv_pvalue(fit)$pval.txt

ggsurvplot(fit, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           #legend.labs = c("OTHER", "WHITE"),
           palette="jto", ylab = "Overall survival" )

###Histology EGFR SMOKER SEX ECOG OS OAK
rm(list=ls())
clinical <- read.csv("16_OAK_clinical.csv")
library(survminer)
library("survival")
fit<- survfit(Surv(OS, OS_status) ~ SEX, data = clinical)
surv_pvalue(fit)$pval.txt

ggsurvplot(fit, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           palette="jto", ylab = "Overall survival" )

##Multi Cox
clinical<-read.csv("16_OAK_clinical.csv",row.names = 1,check.names=FALSE)
res.cox <- coxph(Surv(OS, OS_status) ~ HIST +ECOGGR + cut, data = Type)
x <- summary(res.cox)
pvalue=round(as.matrix(x$coefficients)[,5],4)
HR=round(as.matrix(x$coefficients)[,2],4)
low=round(x$conf.int[,3],4)
high=round(x$conf.int[,4],4)
multi_res=data.frame(p.value=pvalue,
                     HR,
                     low,
                     high,
                     stringsAsFactors = F
)
multi_res
multi_res$HR.CI95<-paste0(multi_res$HR," (",multi_res$low,"-",multi_res$high,")")
write.csv(multi_res,file = "19_multi_.csv")



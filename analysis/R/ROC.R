library(ROCR)
options(stringsAsFactors=FALSE)
x <- read.table("~/Desktop/ROC.data", sep=",", row.names=1, header=FALSE)
df <- as.data.frame(t(x))
names(df) <- row.names(x)
rownames(df)<-NULL
df1 <- data.matrix(df,rownames.force = NA)

base.predict <- prediction(df1[,"base.predictions"],df1[,"base.labels"])
base.perform <- performance(base.predict, 'tpr', 'fpr')

partial.predict <- prediction(df1[,"partial.predictions"],df1[,"partial.labels"])
partial.perform <- performance(partial.predict, 'tpr', 'fpr')

weighted.predict <- prediction(df1[,"weighted.predictions"],df1[,"weighted.labels"])
weighted.perform <- performance(weighted.predict, 'tpr', 'fpr')
 
super.predict <- prediction(df1[,"super.predictions"],df1[,"super.labels"])
super.perform <- performance(super.predict, 'tpr', 'fpr')


# AUC Values
base.auc <- performance(base.predict,"auc")
base.auc <- unlist(slot(base.auc, "y.values"))
base.auc<-min(round(base.auc, digits = 3))
base.auc.t <- paste(c("base MSA \t\t\t\t\t(AUC)  = "),base.auc,sep="")

partial.auc <- performance(partial.predict,"auc")
partial.auc <- unlist(slot(partial.auc, "y.values"))
partial.auc<-min(round(partial.auc, digits = 3))
partial.auc.t <- paste(c("Partial SuperMSA \t\t\t(AUC)  = "),partial.auc,sep="")

weighted.auc <- performance(weighted.predict,"auc")
weighted.auc <- unlist(slot(weighted.auc, "y.values"))
weighted.auc<-min(round(weighted.auc, digits = 3))
weighted.auc.t <- paste(c("Partial Weighted SuperMSA \t(AUC)  = "),weighted.auc,sep="")

super.auc <- performance(super.predict,"auc")
super.auc <- unlist(slot(super.auc, "y.values"))
super.auc<-min(round(super.auc, digits = 3))
super.auc.t <- paste(c("SuperMSA \t\t\t\t\t(AUC)  = "),super.auc,sep="")



plot (base.perform, col="#e41a1c", lwd=3, main='ROC curve accuracy of bootstrap support for 10 Simulated Datasets (16, 32, 64 Tips)')
plot(partial.perform, add= TRUE, col="#965628",lwd=3)
plot(weighted.perform, add= TRUE, col="#ffff33",lwd=3)
plot(super.perform, add= TRUE, col="#ff7f00",lwd=3)
legend(0.4,0.4,c(base.auc.t,partial.auc.t,weighted.auc.t,super.auc.t), col=c('#e41a1c','#965628', '#ffff33', '#ff7f00'),lwd=3)


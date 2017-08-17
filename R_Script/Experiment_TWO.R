#### Package for foldchange ratio
install.packages("gtools")
install.packages("dplyr")
library(gtools)
library(survival)
library(survminer)
library(dplyr)


#### Reading TP and NT data
Merged_Dataset <- read.csv(file="/Users/jiteshpatel/Documents/Ryerson/NEW_MRP/GENE_EXPRESSIOn/MERGED_TP_Data.csv", header=TRUE, sep=",")
Merged_NT_Dataset <- read.csv(file="/Users/jiteshpatel/Documents/Ryerson/NEW_MRP/GENE_EXPRESSIOn/MERGED_NT_Data.csv", header=TRUE, sep=",")


Merged_Dataset <- read.csv(file="C:/Users/jpatel/Documents/SUBMISSION/DATA/MERGED_TP_Data.csv", header=TRUE, sep=",")
Merged_NT_Dataset <- read.csv(file="C:/Users/jpatel/Documents/SUBMISSION/DATA/MERGED_NT_Data.csv", header=TRUE, sep=",")


test_TP_data <-Merged_Dataset[,c("tcga_participant_barcode","AKT1","AKT2","AKT3","MYC","MYCL1","MYCN","PRKAA1","PRKAA2","PRKAB1","PRKAB2","RB1","STK11","CAMKK2","PTEN","PIK3CA","MTOR", "ATM", "TP53", "CHEK2")]
head(test_TP_data)
test_NT_data <- Merged_NT_Dataset[,c("tcga_participant_barcode","AKT1","AKT2","AKT3","MYC","MYCL1","MYCN","PRKAA1","PRKAA2","PRKAB1","PRKAB2","RB1","STK11","CAMKK2","PTEN","PIK3CA","MTOR", "ATM", "TP53", "CHEK2")]

head(test_NT_data)
summary(Merged_Dataset$STK11)
summary(test_NT_data$AKT1)

#### Removing four respondents without mRNA gene expressions

#test_TP_data <- test_TP_data[-c(382,579,645,1063),]
#test_NT_data <- test_NT_data[-c(382,579,645,1063),]

#### Replacing NA with median in NT dataset

AKT1.median <- median(test_NT_data$AKT1, na.rm=TRUE)
AKT2.median <- median(test_NT_data$AKT2, na.rm=TRUE)
AKT3.median <- median(test_NT_data$AKT3, na.rm=TRUE)
MYC.median <- median(test_NT_data$MYC, na.rm=TRUE)
MYCL1.edian <- median(test_NT_data$MYCL1, na.rm=TRUE)
MYCN.median <- median(test_NT_data$MYCN, na.rm=TRUE)
PRKAA1.median <- median(test_NT_data$PRKAA1, na.rm=TRUE)
PRKAA2.median <- median(test_NT_data$PRKAA2, na.rm=TRUE)
PRKAB1.median <- median(test_NT_data$PRKAB1, na.rm=TRUE)
PRKAB2.median <- median(test_NT_data$PRKAB2, na.rm=TRUE)
RB1.median <- median(test_NT_data$RB1, na.rm=TRUE)
STK11.median <- median(test_NT_data$STK11, na.rm=TRUE)
CAMKK2.median <- median(test_NT_data$CAMKK2, na.rm=TRUE)
PTEN.median <- median(test_NT_data$PTEN, na.rm=TRUE)
PIK3CA.median <- median(test_NT_data$PIK3CA, na.rm=TRUE)
MTOR.median <- median(test_NT_data$MTOR, na.rm=TRUE)

ATM.median <- median(test_NT_data$ATM, na.rm=TRUE)
TP53.median <- median(test_NT_data$TP53, na.rm=TRUE)
CHEK2.median <- median(test_NT_data$CHEK2, na.rm=TRUE)


New_NT_data <- test_NT_data

New_NT_data$AKT1[is.na(New_NT_data$AKT1)] = AKT1.median
New_NT_data$AKT2[is.na(New_NT_data$AKT2)] = AKT2.median
New_NT_data$AKT3[is.na(New_NT_data$AKT3)] = AKT3.median
New_NT_data$MYC[is.na(New_NT_data$MYC)] = MYC.median
New_NT_data$MYCL1[is.na(New_NT_data$MYCL1)] = MYCL1.edian
New_NT_data$MYCN[is.na(New_NT_data$MYCN)] = MYCN.median
New_NT_data$PRKAA1[is.na(New_NT_data$PRKAA1)] = PRKAA1.median
New_NT_data$PRKAA2[is.na(New_NT_data$PRKAA2)] = PRKAA2.median
New_NT_data$PRKAB1[is.na(New_NT_data$PRKAB1)] = PRKAB1.median
New_NT_data$PRKAB2[is.na(New_NT_data$PRKAB2)] = PRKAB2.median
New_NT_data$RB1[is.na(New_NT_data$RB1)] = RB1.median
New_NT_data$STK11[is.na(New_NT_data$STK11)] = STK11.median
New_NT_data$CAMKK2[is.na(New_NT_data$CAMKK2)] = CAMKK2.median
New_NT_data$PTEN[is.na(New_NT_data$PTEN)] = PTEN.median
New_NT_data$PIK3CA[is.na(New_NT_data$PIK3CA)] = PIK3CA.median
New_NT_data$MTOR[is.na(New_NT_data$MTOR)] = MTOR.median

New_NT_data$ATM[is.na(New_NT_data$ATM)] = ATM.median
New_NT_data$TP53[is.na(New_NT_data$TP53)] = TP53.median
New_NT_data$CHEK2[is.na(New_NT_data$CHEK2)] = CHEK2.median

head(New_NT_data)
head(test_NT_data)
head(test_TP_data)

write.csv(New_NT_data, file = "C:/Users/jpatel/Documents/SUBMISSION/DATA/NT_with_median.csv")

write.csv(New_NT_data, file = "/Users/jiteshpatel/Documents/Ryerson/NEW_MRP/GENE_EXPRESSIOn/NT_with_median.csv")

#### Fold Change Calculation

FC <- cbind(test_TP_data[1],
            foldchange(test_TP_data$AKT1,New_NT_data$AKT1),
            foldchange(test_TP_data$AKT2,New_NT_data$AKT2),
            foldchange(test_TP_data$AKT3,New_NT_data$AKT3),
            foldchange(test_TP_data$MYC,New_NT_data$MYC),
            foldchange(test_TP_data$MYCL1,New_NT_data$MYCL1),
            foldchange(test_TP_data$MYCN,New_NT_data$MYCN),
            foldchange(test_TP_data$PRKAA1,New_NT_data$PRKAA1),
            foldchange(test_TP_data$PRKAA2,New_NT_data$PRKAA2),
            foldchange(test_TP_data$PRKAB1,New_NT_data$PRKAB1),
            foldchange(test_TP_data$PRKAB2,New_NT_data$PRKAB2),
            foldchange(test_TP_data$RB1,New_NT_data$RB1),
            foldchange(test_TP_data$STK11,New_NT_data$STK11),
            foldchange(test_TP_data$CAMKK2,New_NT_data$CAMKK2),
            foldchange(test_TP_data$PTEN,New_NT_data$PTEN),
            foldchange(test_TP_data$PIK3CA,New_NT_data$PIK3CA),
            foldchange(test_TP_data$MTOR,New_NT_data$MTOR),
            foldchange(test_TP_data$ATM,New_NT_data$ATM),
            foldchange(test_TP_data$TP53,New_NT_data$TP53),
            foldchange(test_TP_data$CHEK2,New_NT_data$CHEK2))


head(FC)

colnames(FC) <- c("tcga_participant_barcode","FC.AKT1","FC.AKT2","FC.AKT3","FC.MYC","FC.MYCL1","FC.MYCN","FC.PRKAA1","FC.PRKAA2","FC.PRKAB1","FC.PRKAB2","FC.RB1","FC.STK11","FC.CAMKK2","FC.PTEN","FC.PIK3CA","FC.MTOR","FC.ATM","FC.TP53","FC.CHEK2")
write.csv(FC, file = "C:/Users/jpatel/Documents/SUBMISSION/DATA/FOLD_CHANGE.csv")
FC <- read.csv(file="C:/Users/jpatel/Documents/SUBMISSION/DATA/FOLD_CHANGE.csv")

summary(FC$FC.AKT1)
#### Fold Change Log Ratio Calculation
FC_log <-0
data.frame(FC_log)

FC_log <- data.frame(FC$tcga_participant_barcode,foldchange2logratio(FC$FC.AKT1),foldchange2logratio(FC$FC.AKT2),foldchange2logratio(FC$FC.AKT3),foldchange2logratio(FC$FC.MYC),foldchange2logratio(FC$FC.MYCL1),foldchange2logratio(FC$FC.MYCN),foldchange2logratio(FC$FC.PRKAA1),foldchange2logratio(FC$FC.PRKAA2),foldchange2logratio(FC$FC.PRKAB1),foldchange2logratio(FC$FC.PRKAB2),foldchange2logratio(FC$FC.RB1),foldchange2logratio(FC$FC.STK11),foldchange2logratio(FC$FC.CAMKK2),foldchange2logratio(FC$FC.PTEN),foldchange2logratio(FC$FC.PIK3CA),foldchange2logratio(FC$FC.MTOR),foldchange2logratio(FC$FC.ATM),foldchange2logratio(FC$FC.TP53),foldchange2logratio(FC$FC.CHEK2))
colnames(FC_log) <- c("tcga_participant_barcode","L2.AKT1","L2.AKT2","L2.AKT3","L2.MYC","L2.MYCL1","L2.MYCN","L2.PRKAA1","L2.PRKAA2","L2.PRKAB1","L2.PRKAB2","L2.RB1","L2.STK11","L2.CAMKK2","L2.PTEN","L2.PIK3CA","L2.MTOR","L2.ATM","L2.TP53","L2.CHEK2")

#### Saving dataset with Fold Change Log 2 Ratio values
write.csv(FC_log, file = "C:/Users/jpatel/Documents/SUBMISSION/DATA/FOLDCHANGE_LOGRATIO.csv")

summary(FC_log$L2.AKT1)
summary(FC_log$L2.AKT2)
summary(FC_log$L2.AKT3)
summary(FC_log$L2.PRKAA1)
summary(FC_log$L2.STK11)


#### Read clinical data
clinical <- read.csv(file="C:/Users/jpatel/Documents/SUBMISSION/DATA/clinical.csv", header=TRUE, sep=",")

#### Removing four records without mRNA expression data from the dataset and merging FoldChange log ratio

clinical <- clinical[-c(382,439,579,645,1063),]
Clin_log <- left_join(clinical, FC_log, by = "tcga_participant_barcode")
Clin_log <- left_join(Clin_log, FC, by = "tcga_participant_barcode")

head(Clin_log)
nrow(Clin_log)
#### SAving dataset with clinical + Log2FoldChange

write.csv(Clin_log, file = "C:/Users/jpatel/Documents/SUBMISSION/DATA/Clinical_logratio2FC.csv")

write.csv(Clin_log, file = "/Users/jiteshpatel/Documents/Ryerson/NEW_MRP/GENE_EXPRESSIOn/Clinical_logratio2FC.csv")

################ Multivariate COX Regression Analysis ##########################


test1 <- Clin_log[-c(254),]

test1$ER <- 'none'
test1$ER[test1$er_level_cell_percentage_category == "<10%"] <- 'negative'
test1$ER[test1$er_level_cell_percentage_category == "10-19%"] <- 'positive'
test1$ER[test1$er_level_cell_percentage_category == "20-29%"] <- 'positive'
test1$ER[test1$er_level_cell_percentage_category == "30-39%"] <- 'positive'
test1$ER[test1$er_level_cell_percentage_category == "40-49%"] <- 'positive'
test1$ER[test1$er_level_cell_percentage_category == "50-59%"] <- 'positive'
test1$ER[test1$er_level_cell_percentage_category == "60-69%"] <- 'positive'
test1$ER[test1$er_level_cell_percentage_category == "70-79%"] <- 'positive'
test1$ER[test1$er_level_cell_percentage_category == "80-89%"] <- 'positive'
test1$ER[test1$er_level_cell_percentage_category == "90-99%"] <- 'positive'

test1$PG <- 'none'
test1$PG[test1$progesterone_receptor_level_cell_percent_category == "<10%"] <- 'negative'
test1$PG[test1$progesterone_receptor_level_cell_percent_category == "<10%"] <- 'negative'
test1$PG[test1$progesterone_receptor_level_cell_percent_category == "10-19%"] <- 'positive'
test1$PG[test1$progesterone_receptor_level_cell_percent_category == "20-29%"] <- 'positive'
test1$PG[test1$progesterone_receptor_level_cell_percent_category == "30-39%"] <- 'positive'
test1$PG[test1$progesterone_receptor_level_cell_percent_category == "40-49%"] <- 'positive'
test1$PG[test1$progesterone_receptor_level_cell_percent_category == "50-59%"] <- 'positive'
test1$PG[test1$progesterone_receptor_level_cell_percent_category == "60-69%"] <- 'positive'
test1$PG[test1$progesterone_receptor_level_cell_percent_category == "70-79%"] <- 'positive'
test1$PG[test1$progesterone_receptor_level_cell_percent_category == "80-89%"] <- 'positive'
test1$PG[test1$progesterone_receptor_level_cell_percent_category == "90-99%"] <- 'positive'

test1$HER <- 'none'
test1$HER[test1$her2_erbb_pos_finding_cell_percent_category == "<10%"] <- 'negative'
test1$HER[test1$her2_erbb_pos_finding_cell_percent_category == "<10%"] <- 'negative'
test1$HER[test1$her2_erbb_pos_finding_cell_percent_category == "10-19%"] <- 'positive'
test1$HER[test1$her2_erbb_pos_finding_cell_percent_category == "20-29%"] <- 'positive'
test1$HER[test1$her2_erbb_pos_finding_cell_percent_category == "30-39%"] <- 'positive'
test1$HER[test1$her2_erbb_pos_finding_cell_percent_category == "40-49%"] <- 'positive'
test1$HER[test1$her2_erbb_pos_finding_cell_percent_category == "50-59%"] <- 'positive'
test1$HER[test1$her2_erbb_pos_finding_cell_percent_category == "60-69%"] <- 'positive'
test1$HER[test1$her2_erbb_pos_finding_cell_percent_category == "70-79%"] <- 'positive'
test1$HER[test1$her2_erbb_pos_finding_cell_percent_category == "80-89%"] <- 'positive'
test1$HER[test1$her2_erbb_pos_finding_cell_percent_category == "90-99%"] <- 'positive'


write.csv(test1, file = "C:/Users/jpatel/Documents/SUBMISSION/DATA/Clin_logratio_class.csv")


#### Read clinical+log2change ratio data
test1 <- read.csv(file="C:/Users/jpatel/Documents/SUBMISSION/DATA/Clin_logratio_class.csv", header=TRUE, sep=",")


### Triple Negative + HER2
test1 <- test1[which((test1$ER == 'negative' & test1$PG == 'negative' & test1$HER == 'negative') |
                       (test1$ER == 'negative' & test1$PG == 'negative' & test1$HER == 'positive')),]

### LB
test1 <- test1[which((test1$ER == 'positive' & test1$PG == 'positive' & test1$HER == 'positive') |
                       (test1$ER == 'negative' & test1$PG == 'positive' & test1$HER == 'positive') |
                       (test1$ER == 'positive' & test1$PG == 'negative' & test1$HER == 'positive')),]

### LA
test1 <- test1[which((test1$ER == 'positive' & test1$PG == 'positive' & test1$HER == 'negative') |
                       (test1$ER == 'positive' & test1$PG == 'negative' & test1$HER == 'negative') |
                       (test1$ER == 'negative' & test1$PG == 'positive' & test1$HER == 'negative')),]

nrow(test1)
tail(test1)


#### Sorting 

#test1 <- test1[with(test1,order(L2.STK11)),]
#test1 <- test1[with(test1,order(L2.MTOR)),]
#test1 <- test1[with(test1,order(L2.RB1)),]
#test1 <- test1[with(test1,order(L2.PRKAA1)),]
#test1 <- test1[with(test1,order(L2.PRKAA2)),]
#test1 <- test1[with(test1,order(L2.PRKAB1)),]
#test1 <- test1[with(test1,order(L2.PRKAB2)),]
#test1 <- test1[with(test1,order(L2.AKT1)),]
#test1 <- test1[with(test1,order(L2.AKT2)),]
#test1 <- test1[with(test1,order(L2.AKT3)),]
#test1 <- test1[with(test1,order(L2.MYC)),]
#test1 <- test1[with(test1,order(L2.MYCL1)),]
#test1 <- test1[with(test1,order(L2.MYCN)),]
#test1 <- test1[with(test1,order(L2.CAMKK2)),]
#test1 <- test1[with(test1,order(L2.PTEN)),]
#test1 <- test1[with(test1,order(L2.PIK3CA)),]
#test1 <- test1[with(test1,order(L2.ATM)),]
#test1 <- test1[with(test1,order(L2.TP53)),]
#test1 <- test1[with(test1,order(L2.CHEK2)),]


#### COXPH regression model

res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.STK11, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.MTOR, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.RB1, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.PRKAA1, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.PRKAA2, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.PRKAB1, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.PRKAB2, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.AKT1, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.AKT2, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.AKT3, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.MYC, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.MYCL1, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.MYCN, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.CAMKK2, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.PTEN, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.PIK3CA, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.ATM, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.TP53, data = test1)
res.cox1 <- coxph(Surv(new_death,death_event) ~ L2.CHEK2, data = test1)

summary(res.cox1)

results <- test1[c('L2.STK11','L2.MTOR','L2.RB1','L2.PRKAA1','L2.PRKAA2','L2.PRKAB1','L2.PRKAB2','L2.AKT1','L2.AKT2','L2.AKT3','L2.MYC','L2.MYCL1','L2.MYCN','L2.CAMKK2','L2.PTEN','L2.PIK3CA','L2.ATM','L2.TP53','L2.CHEK2')]


results <- cbind(results, lp=predict(res.cox1, type="lp"))
results <- cbind(results, risk=predict(res.cox1, type="risk"))
results <- cbind(results, expected=predict(res.cox1, type="expected"))
results$exp.expected <- exp(-(results$expected))
l.coxph.age.results <- data.frame(results)
l.coxph.age.results$ID <- seq.int(nrow(l.coxph.age.results))
 

#################################################


attach(l.coxph.age.results)

plot(L2.STK11, risk, main="STK11 vs RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.MTOR, risk, main="MTOR vs RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.RB1, risk, main="RB1 vs RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.PRKAA1, risk, main="PRKAA1 vs RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.PRKAA2, risk, main="PRKAA2 VS RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.PRKAB1, risk, main="PRKAB1 VS RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.PRKAB2, risk, main="PRKAB2 VS RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.AKT1, risk, main="AKT1 VS RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.AKT2, risk, main="AKT2 VS RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.AKT3, risk, main="AKT3 VS RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.MYC, risk, main="MYC VS RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.MYCL1, risk, main="MYCL1 VS RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.MYCN, risk, main="MYCN VS RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.CAMKK2, risk, main="CAMKK2 VS RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.PTEN, risk, main="PTEN VS RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.PIK3CA, risk, main="PIK3CA VS RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.ATM, risk, main="ATM VS RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.TP53, risk, main="TP53 VS RISK", 
     xlab="FC", ylab="RISK", pch=19)

plot(L2.CHEK2, risk, main="CHEK2 VS RISK", 
     xlab="FC", ylab="RISK", pch=19)






### Kaplan Miere
pos <- 0
mygroup_STK11 <- NULL
for(i in length(test1$fc.STK11)){
  qq <- test1$fc.STK11[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.33, 0.66), na.rm=TRUE)
  qq[test1$fc.STK11[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[test1$fc.STK11[(pos+1):(pos+i)] >= myq[1] & test1$fc.STK11[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[test1$fc.STK11[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup_STK11 <- c(mygroup_STK11,qq)
  pos <- pos + i
}
test1$group_STK11 <-mygroup_STK11
tail(test1)

sfit <- survfit(Surv(new_death,death_event) ~ group_STK11, data = test1) #### fit on STK11
ggsurvplot ( sfit,
             pval = TRUE, conf.int = TRUE,
             risk.table = TRUE,               # Add risk table
             risk.table.col = "strata",       # Change risk table color by groups
             linetype = "strata",             # Change line type by groups
             surv.median.line = "hv",         # Specify median survival
             ggtheme = theme_bw(),            # Change ggplot2 theme
             legend.labs = c("LOW  ", "INTERMEDIATE   ", "HIGH   "),     # change legend labels.
             palette = c("#E7B800", "#2E9FDF","#0000FF"))
####


##################################################
install.packages("ggplot2")
require(ggplot2)
##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
l.coxph.age.results$threshold = as.factor(abs(l.coxph.age.results$fc.STK11) > 2)

##Construct the plot object
g = ggplot(data=l.coxph.age.results, aes(x=fc.STK11, y=expected, colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("Survival Probability")
g

## Check for violation of proportional hazard (constant HR over time)
(res.zph1 <- cox.zph(res.cox1))
plot(res.zph1)

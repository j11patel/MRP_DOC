###### INSTALLATION OF REQUIRED PACKAGE
#install.packages("survival")
#install.packages("survminer")
#install.packages("caret")
#library(survival)
#library(survminer)
#library(caret)

###### Read clinical data
clinical <- read.csv(file="C:/Users/jpatel/Documents/Ryerson/1_MRP/clinical.csv", header=TRUE, sep=",")
#write.csv(clinical, file = "C:/Users/jpatel/Documents/Ryerson/1_MRP/test.csv",row.names=TRUE,col.names = TRUE, na="")

###### Removing four records without mRNA expression data from the dataset
clinical <- clinical[-c(382,579,645,1063),]

###### CREATING HIGH, INTERMEDIATE and LOW categories using quantile  function

###### STK11 groups
pos <- 0
mygroup <- NULL
for(i in length(clinical$STK11)){
  qq <- clinical$STK11[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.30, 0.60), na.rm=TRUE)
  qq[clinical$STK11[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[clinical$STK11[(pos+1):(pos+i)] >= myq[1] & clinical$STK11[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[clinical$STK11[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup <- c(mygroup,qq)
  pos <- pos + i
}

###### RB1 groups
pos <- 0
mygroup_RB1 <- NULL
for(i in length(clinical$RB1)){
  qq <- clinical$RB1[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.30, 0.60), na.rm=TRUE)
  qq[clinical$RB1[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[clinical$RB1[(pos+1):(pos+i)] >= myq[1] & clinical$RB1[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[clinical$RB1[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup_RB1 <- c(mygroup_RB1,qq)
  pos <- pos + i
}

###### PRKAA2 groups
pos <- 0
mygroup_PRKAA2 <- NULL
for(i in length(clinical$PRKAA2)){
  qq <- clinical$PRKAA2[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.30, 0.60), na.rm=TRUE)
  qq[clinical$PRKAA2[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[clinical$PRKAA2[(pos+1):(pos+i)] >= myq[1] & clinical$PRKAA2[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[clinical$PRKAA2[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup_PRKAA2 <- c(mygroup_PRKAA2,qq)
  pos <- pos + i
}

###### Merginf STK11, PRKAA2, RB1 Categories to CLINICAL dataset

clinical$group_STK11 <- mygroup
clinical$group_PRKAA2 <- mygroup_PRKAA2
clinical$group_RB1 <- mygroup_RB1


###### MTOR groups
pos <- 0
mygroup_MTOR <- NULL
for(i in length(clinical$MTOR)){
  qq <- clinical$MTOR[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.30, 0.60), na.rm=TRUE)
  qq[clinical$MTOR[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[clinical$MTOR[(pos+1):(pos+i)] >= myq[1] & clinical$MTOR[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[clinical$MTOR[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup_MTOR <- c(mygroup_MTOR,qq)
  pos <- pos + i
}

clinical$group_MTOR <- mygroup_MTOR


################## HIgh EXpression of two genes #####################

#################### STK11_PRKAA2

clinical$STK11_PRKAA2[clinical$group_STK11 == 3 & clinical$group_PRKAA2 == 3] <- 3 
clinical$STK11_PRKAA2[clinical$group_STK11 == 1 & clinical$group_PRKAA2 == 1] <- 1 
clinical$STK11_PRKAA2[clinical$group_STK11 == 2 & clinical$group_PRKAA2 == 3] <- 3 
clinical$STK11_PRKAA2[clinical$group_STK11 == 3 & clinical$group_PRKAA2 == 2] <- 3 
clinical$STK11_PRKAA2[is.na(clinical$STK11_PRKAA2)] <- 2

#################### MTOR_STK11

clinical$MTOR_STK11[clinical$group_STK11 == 3 & clinical$group_MTOR == 3] <- 3 
clinical$MTOR_STK11[clinical$group_STK11 == 1 & clinical$group_MTOR == 1] <- 1 
clinical$MTOR_STK11[clinical$group_STK11 == 2 & clinical$group_MTOR == 3] <- 3 
clinical$MTOR_STK11[clinical$group_STK11 == 3 & clinical$group_MTOR == 2] <- 3 
clinical$MTOR_STK11[is.na(clinical$MTOR_STK11)] <- 2


#################### MTOR_PRKAA2

clinical$MTOR_PRKAA2[clinical$group_PRKAA2 == 3 & clinical$group_MTOR == 3] <- 3 
clinical$MTOR_PRKAA2[clinical$group_PRKAA2 == 1 & clinical$group_MTOR == 1] <- 1 
clinical$MTOR_PRKAA2[clinical$group_PRKAA2 == 2 & clinical$group_MTOR == 3] <- 3 
clinical$MTOR_PRKAA2[clinical$group_PRKAA2 == 3 & clinical$group_MTOR == 2] <- 3 
clinical$MTOR_PRKAA2[is.na(clinical$MTOR_PRKAA2)] <- 2

################## SURVIVAL RATE ###################################

sfit <- survfit(Surv(new_death,death_event) ~ group_STK11, data = clinical) #### fit on STK11
sfit <- survfit(Surv(new_death,death_event) ~ group_PRKAA2, data = clinical) #### fit on PRKAA2
sfit <- survfit(Surv(new_death,death_event) ~ group_RB1, data = clinical) #### fit on RB1
sfit <- survfit(Surv(new_death,death_event) ~ group_test_STK11, data = clinical) #### fit on convoluted STK11 & PRKAA2
sfit <- survfit(Surv(new_death,death_event) ~ group_test_MTOR_STK11, data = clinical) #### fit on convoluted MTOR & STK11
sfit <- survfit(Surv(new_death,death_event) ~ group_test_MTOR_PRKAA2, data = clinical) #### fit on convoluted MTOR & PRKAA2

sfit <- survfit(Surv(new_death,death_event) ~ STK11_PRKAA2, data = clinical) #### fit on STK11_PRKAA2
sfit <- survfit(Surv(new_death,death_event) ~ MTOR_STK11, data = clinical) #### fit on MTOR_STK11
sfit <- survfit(Surv(new_death,death_event) ~ MTOR_PRKAA2, data = clinical) #### fit on MTOR_PRKAA2


ggsurvplot ( sfit,
             pval = TRUE, conf.int = TRUE,
             risk.table = TRUE,               # Add risk table
             risk.table.col = "strata",       # Change risk table color by groups
             linetype = "strata",             # Change line type by groups
             surv.median.line = "hv",         # Specify median survival
             ggtheme = theme_bw(),            # Change ggplot2 theme
             legend.labs = c("LOW  ", "INTERMEDIATE   ", "HIGH   "),     # change legend labels.
             palette = c("#E7B800", "#2E9FDF","#0000FF"))

#### GRAPH with more details
ggsurvplot(
  sfit,                      # survfit object with calculated statistics.
  pval = TRUE,               # show p-value of log-rank test.
  conf.int = TRUE,           # show confidence intervals for point estimaes of survival curves.
  conf.int.style = "step",   # customize style of confidence intervals
  xlab = "Time in days",     # customize X axis label.
  break.time.by = 2000,      # break X axis in time intervals by 200.
  ggtheme = theme_bw(),   # customize plot and risk table with a theme.
  risk.table.col = "strata",
  risk.table = "abs_pct",    # absolute number and percentage at risk.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE, # show bars instead of names in text annotations in legend of risk table.
  linetype = "strata",             # Change line type by groups
  ncensor.plot = FALSE,       # plot the number of censored subjects at time t
  surv.median.line = "hv",   # add the median survival pointer.
  legend.labs = c("LOW  ", "INTERMEDIATE   ", "HIGH   "),     # change legend labels.
  palette = c("#E7B800", "#2E9FDF","#0000FF")      # custom color palettes.
)



################## COMBINATION using CONVOLUTION ####################

x= clinical$STK11
y = clinical$PRKAA2
z = clinical$RB1
m = clinical$MTOR

clinical$test_STK11 <- convolve(x, y, conj = TRUE, type = "circular")
clinical$test_PRKAA2 <- convolve(y, z, conj = TRUE, type = "circular")

clinical$test_MTOR_STK11 <- convolve(x, m, conj = TRUE, type = "circular")
clinical$test_MTOR_PRKAA2 <- convolve(y, m, conj = TRUE, type = "circular")
  
################## CONVOLUTION GROUP ################################

###### STK11 and PRKAA
pos <- 0
mygroup_test_STK11 <- NULL
for(i in length(clinical$test_STK11)){
  qq <- clinical$test_STK11[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.30, 0.60), na.rm=TRUE)
  qq[clinical$test_STK11[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[clinical$test_STK11[(pos+1):(pos+i)] >= myq[1] & clinical$test_STK11[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[clinical$test_STK11[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup_test_STK11 <- c(mygroup_test_STK11,qq)
  pos <- pos + i
}
clinical$group_test_STK11 <- mygroup_test_STK11

###### MTOR and STK11
pos <- 0
mygroup_test_MTOR_STK11 <- NULL
for(i in length(clinical$test_MTOR_STK11)){
  qq <- clinical$test_MTOR_STK11[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.30, 0.60), na.rm=TRUE)
  qq[clinical$test_MTOR_STK11[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[clinical$test_MTOR_STK11[(pos+1):(pos+i)] >= myq[1] & clinical$test_MTOR_STK11[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[clinical$test_MTOR_STK11[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup_test_MTOR_STK11 <- c(mygroup_test_MTOR_STK11,qq)
  pos <- pos + i
}
clinical$group_test_MTOR_STK11 <- mygroup_test_MTOR_STK11

###### MTOR and PRKAA
pos <- 0
mygroup_test_MTOR_PRKAA2 <- NULL
for(i in length(clinical$test_MTOR_PRKAA2)){
  qq <- clinical$test_MTOR_PRKAA2[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.30, 0.60), na.rm=TRUE)
  qq[clinical$test_MTOR_PRKAA2[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[clinical$test_MTOR_PRKAA2[(pos+1):(pos+i)] >= myq[1] & clinical$test_MTOR_PRKAA2[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[clinical$test_MTOR_PRKAA2[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup_test_MTOR_PRKAA2 <- c(mygroup_test_MTOR_PRKAA2,qq)
  pos <- pos + i
}
clinical$group_test_MTOR_PRKAA2 <- mygroup_test_MTOR_PRKAA2

########### COXPH for (MTOR and STK11) AND (MTOR AND PRKAA2)###########

res.cox1 <- coxph(Surv(new_death,death_event) ~ STK11+MTOR, data = clinical)
summary(res.cox1)
ggsurvplot(survfit(res.cox1), palette= '#2E9FDF',
           ggtheme = theme_minimal(), title="STK11 + MTOR")

res.cox2 <- coxph(Surv(new_death,death_event) ~ PRKAA2+MTOR, data = clinical)
summary(res.cox2)
ggsurvplot(survfit(res.cox2), palette= '#2E9F00',
           ggtheme = theme_minimal(), title="PRKAA2 + MTOR")

############### Extracting Survival Rate for each records ############

ggplot(test1, aes(x=PRKAA2, y=STK11)) + geom_point()

results <- test1[c('new_death','death_event','test_PRKAA2')]

results <- cbind(results, lp=predict(res.cox1, type="lp"))
results <- cbind(results, risk=predict(res.cox1, type="risk"))
results <- cbind(results, expected=predict(res.cox1, type="expected"))
#exp(-0.0122572336)
results$exp.expected <- exp(-(results$expected))
results <- cbind(results, terms=predict(res.cox1, type="terms"))
l.coxph.age.results <- results
head(l.coxph.age.results)
tail(l.coxph.age.results)
nrow(l.coxph.age.results)

l.coxph.age.results$ID <- seq.int(nrow(l.coxph.age.results))

ggplot(l.coxph.age.results, aes(x=exp.expected, y=test_PRKAA2)) + geom_point() + geom_smooth(method="lm")


##################### Kaplan Miere plot for COXPH #####################

source("http://bioconductor.org/biocLite.R") 
biocLite("survcomp")
library(survcomp)
dd <- data.frame("time"=clinical$new_death, "event"=clinical$death_event, "gg"=clinical$group_STK11)
km.coxph.plot(formula.s=formula(Surv(time,event) ~ mygroup), data.s = dd, x.label = "days", y.label = "SP", main.title = "STK11")

#######################################################################

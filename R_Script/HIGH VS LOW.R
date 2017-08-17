test1 <- read.csv(file="C:/Users/jpatel/Documents/FINAL_MRP/DATASETS/Clin_logratio_class_two.csv", header=TRUE, sep=",")

test1 <- read.csv(file="C:/Users/jpatel/Documents/SUBMISSION/DATA/Clin_logratio_class.csv", header=TRUE, sep=",")


#### MEDIAN FOLD CHANGE EXPRESSION OF AMPK (WITHOUT PRKAA2 BECAUSE ITS OPPOSITE TO REST)

for(i in 1:nrow(test1)){
  test1$AMPK_MEDIAN[i] <- median(c(test1$FC.PRKAA1[i],test1$FC.PRKAB1[i],test1$FC.PRKAB2[i],test1$FC.MTOR[i], test1$FC.STK11[i], test1$FC.AKT2[i], test1$FC.AKT3[i],test1$FC.MYC[i],test1$FC.MYCL1[i],test1$FC.MYCN[i],test1$FC.CAMKK2[i],test1$FC.PTEN[i],test1$FC.TP53[i],test1$FC.ATM[i]))  
}

#test123 <- median(c(test1$FC.PRKAA1[1],test1$FC.PRKAB1[1],test1$FC.PRKAB2[1],test1$FC.MTOR[1], test1$FC.STK11[1], test1$FC.AKT2[1], test1$FC.AKT3[1],test1$FC.MYC[1],test1$FC.MYCL1[1],test1$FC.MYCN[1],test1$FC.CAMKK2[1],test1$FC.PTEN[1]))  

for(i in 1:nrow(test1)){
  test1$LOW_MEDIAN[i] <- median(c(test1$FC.RB1[i],test1$FC.AKT1[i],test1$FC.PIK3CA[i],test1$FC.PRKAA2[i],test1$FC.CHEK2[i]))  
}

###############################################USING log2FOLDCHANGE
for(i in 1:nrow(test1)){
  test1$AMPK_MEDIAN[i] <- median(c(test1$L2.PRKAA1[i],test1$L2.PRKAB1[i],test1$L2.PRKAB2[i],test1$L2.MTOR[i], test1$L2.STK11[i], test1$L2.AKT2[i], test1$L2.AKT3[i],test1$L2.MYC[i],test1$L2.MYCL1[i],test1$L2.MYCN[i],test1$L2.CAMKK2[i],test1$L2.PTEN[i],test1$L2.TP53[i],test1$L2.ATM[i]))  
}


for(i in 1:nrow(test1)){
  test1$LOW_MEDIAN[i] <- median(c(test1$L2.RB1[i],test1$L2.AKT1[i],test1$L2.PIK3CA[i],test1$L2.PRKAA2[i],test1$L2.CHEK2))  
}
###################################################################

for (i in 1:nrow(test1)){
  test1$fc.diff[i] <- (test1$AMPK_MEDIAN[i]/test1$LOW_MEDIAN[i])
}



###### CREATING HIGH, INTERMEDIATE and LOW categories based on upper lower quartile
###### fc.diff 

pos <- 0
mygroup <- NULL
for(i in length(test1$fc.diff)){
  qq <- test1$fc.diff[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.50), na.rm=TRUE)
  qq[test1$fc.diff[(pos+1):(pos+i)] <= myq[1]] <- 1
  qq[test1$fc.diff[(pos+1):(pos+i)] > myq[1]] <- 2
  qq <- factor(x=qq, levels=1:2)
  mygroup <- c(mygroup,qq)
  pos <- pos + i
}

test1$fc.diff_group <- mygroup

sfit <- survfit(Surv(new_death,death_event) ~ fc.diff_group, data = test1) #### fit on STK11

ggsurvplot ( sfit,
             pval = TRUE, conf.int = TRUE,
             xscale = "d_y",
             xlab = "Time in Years",
             risk.table = TRUE,               # Add risk table
             risk.table.col = "strata",       # Change risk table color by groups
             linetype = "strata",             # Change line type by groups
             surv.median.line = "hv",         # Specify median survival
             ggtheme = theme_bw(),            # Change ggplot2 theme
             legend.labs = c("High diiference",   "Low difference   "),     # change legend labels.
             palette = c("#E7B800", "#2E9FDF"),
             title = "(OVERALL) (MEDIAN OF HIGH EXPRESSED GENES vs MEDIAN OF LOW EXPRESSED GENES")



#####################################TRIPLE NEGATIVE
test2 <- test1[which((test1$ER == 'negative' & test1$PG == 'negative' & test1$HER == 'negative')),]

sfit <- survfit(Surv(new_death,death_event) ~ fc.diff_group, data = test2) #### fit on STK11

ggsurvplot ( sfit,
             pval = TRUE, conf.int = TRUE,
             xscale = "d_y",
             xlab = "Time in Years",
             risk.table = TRUE,               # Add risk table
             risk.table.col = "strata",       # Change risk table color by groups
             linetype = "strata",             # Change line type by groups
             surv.median.line = "hv",         # Specify median survival
             ggtheme = theme_bw(),            # Change ggplot2 theme
             legend.labs = c("High difference   ", "Low difference   "),     # change legend labels.
             palette = c("#E7B800", "#2E9FDF"),
             title = "(TRIPLE NEGATIVE) - (MEDIAN OF HIGH EXPRESSED GENES vs MEDIAN OF LOW EXPRESSED GENES")              # Main Title
#title = "(LUMINAL) - (PRKAA1,PRKAA2,PRKAB1,PRKAB2,MTOR,STK11) vs RB1 FOLD CHANGE")              # Main Title



#####################################LUMINAL

test2 <- test1[which((test1$ER == 'positive' & test1$PG == 'positive' & test1$HER == 'positive') |
                       (test1$ER == 'negative' & test1$PG == 'positive' & test1$HER == 'positive') |
                       (test1$ER == 'positive' & test1$PG == 'negative' & test1$HER == 'positive') |
                       (test1$ER == 'positive' & test1$PG == 'positive' & test1$HER == 'negative') |
                       (test1$ER == 'positive' & test1$PG == 'negative' & test1$HER == 'negative') |
                       (test1$ER == 'negative' & test1$PG == 'positive' & test1$HER == 'negative')),]


sfit <- survfit(Surv(new_death,death_event) ~ fc.diff_group, data = test2) #### fit on STK11

ggsurvplot ( sfit,
             pval = TRUE, conf.int = TRUE,
             xscale = "d_y",
             xlab = "Time in Years",
             risk.table = TRUE,               # Add risk table
             risk.table.col = "strata",       # Change risk table color by groups
             linetype = "strata",             # Change line type by groups
             surv.median.line = "hv",         # Specify median survival
             ggtheme = theme_bw(),            # Change ggplot2 theme
             legend.labs = c("High diiference   ", "Low difference   "),     # change legend labels.
             palette = c("#E7B800", "#2E9FDF"),
             title = "(LUMINAL) - (MEDIAN OF HIGH EXPRESSED GENES vs MEDIAN OF LOW EXPRESSED GENES")              # Main Title




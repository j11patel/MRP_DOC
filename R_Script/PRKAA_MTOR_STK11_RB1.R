test1 <- read.csv(file="C:/Users/jpatel/Documents/SUBMISSION/DATA/Clin_logratio_class.csv", header=TRUE, sep=",")


#### MEDIAN FOLD CHANGE EXPRESSION OF AMPK (EITHOUT PRKAA2 BECAUSE ITS OPPOSITE TO REST)

for(i in 1:nrow(test1)){
  test1$AMPK_MEDIAN[i] <- median(c(test1$FC.PRKAA1[i],test1$FC.PRKAB1[i],test1$FC.PRKAB2[i],test1$FC.MTOR[i], test1$FC.STK11[i]))  
}

for (i in 1:nrow(test1)){
  test1$fc.diff[i] <- (test1$AMPK_MEDIAN[i]/test1$FC.RB1[i])
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
             title = "(OVERALL) (PRKAA1, PRKAB1, PRKAB2, MTOR, STK11 MEDIAN) vs RB1 FOLD CHANGE")



##################### TRIPLE NEGATIVE

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
             legend.labs = c("High diiference   ", "Low Difference   "),     # change legend labels.
             palette = c("#E7B800", "#2E9FDF"),
             title = "(TRIPLE NEGATIVE) - (PRKAA1,PRKAA2,PRKAB1,PRKAB2,MTOR, STK11 MEDIAN) vs RB1 FOLD CHANGE")              # Main Title
             #title = "(LUMINAL) - (PRKAA1,PRKAA2,PRKAB1,PRKAB2,MTOR,STK11) vs RB1 FOLD CHANGE")              # Main Title




### LUMINAL
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
             legend.labs = c("High diiference   ", "Low Difference   "),     # change legend labels.
             palette = c("#E7B800", "#2E9FDF"),
             title = "(LUMINAL) - (PRKAA1,PRKAA2,PRKAB1,PRKAB2,MTOR, STK11 MEDIAN) vs RB1 FOLD CHANGE")              # Main Title




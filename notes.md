

################ To Figureout correlation between mRNA expression of genes and survival rate ###############
I have categorized mRNA expression of STK11, RB1, PRKAA2 and MTOR into three categories (1)Low, (2)Intermediate and (3)High using quantile function to check effect of High mRNA expression of respective gene to Survival Rate.
Survival Plots for STK11, PRKAA2 and RB1 are available in MRP_DOC/PLOTS
I have created new category variable to check effect of HIGH expression of two genes (STK11 and PRKAA2), (MTOR and STK11) and (MTOR and PRKAA2). 
For example, i have prepared a new variable STK11_PRKAA2 (with 3 as high where (STK11 = 3 and PRKAA2 = 3, STK11 = 3 and PRKAA = 2, STK11 =2 and PRKAA = 3)  and prepared survival plot check high mRNA expression of STk11 and PRKAA2 on survival rate.
STk11_PRKAA2 => Survival rate based on STK11 and PRKAA mRNA expression, MTOR_PRKAA2 => Survival rate based on MTOR and PRKAA mRNA expression, MTOR_STK11 => Survival rate based on MTOR and PRKAA mRNA expression are available in MRP_DOC/PLOTS

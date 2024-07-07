## Task 3

library(tidyverse)
library(survival)
library(survminer)

df <- read.delim("Practical2Data.tsv") %>%
  filter(RT %in% c("YES", "NO"))
head(df)

# Following the methodology from Hosmer et al, we start by looking at what
# covariates are significant at a 20% level in a univariate analysis
## STEP 1 ##
# We will use the log-rank test for categorical variables, and Wald tests for
# continuous variables

head(df)

# Categorical: RT and Sex
dfCat <- df %>% select(Surv, Time, RT, Sex)
pvalsCat <- data.frame("Variable" = c("RT", "Sex"),
                       "Test" = c("Log-rank", "Log-rank"),
                       "P.value" = rep(NA, 2))
for (i in 1:nrow(pvalsCat)) {
  diffobj <- survdiff(Surv(Time, Surv) ~ dfCat[,i+2], data = dfCat)
  pval <- 1 - pchisq(diffobj$chisq, df = 1)
  pvalsCat$P.value[i] <- pval
}

# Continuous: Age and gene expressions
dfCont <- df %>% select(-RT, -Sex)
pvalsCont <- data.frame("Variable" = colnames(dfCont)[3:ncol(dfCont)],
                        "Test" = "Wald",
                        "P.value" = NA)
for (i in 1:nrow(pvalsCont)) {
  fitPH <- coxph(Surv(Time, Surv) ~ dfCont[,i+2], data = dfCont)
  sum <- summary(fitPH)
  pval <- sum$waldtest[3]
  pvalsCont$P.value[i] <- pval
}

# Combine these two and order on p-value
pvals <- rbind(pvalsCat, pvalsCont)
pvals$P.value <- round(pvals$P.value, 4)
pvals[order(pvals$P.value),]

# We will be keeping those with p-values lower than 0.2. We will also keep
# KDM5C since we have reason to believe it has an interaction with RT.

## STEPS 2 & 3 ##
# Fit an initial multivariate model
mdlInit <- coxph(Surv(Time, Surv) ~ RT + Age + BIRC3 + CD79B + MALT1 + JAK3 +
                   HOXA13 + TBL1XR1 + TSHR + GATA3 + PIM1 + CSF3R + AR + KDM5C,
                 data = df, ties = "exact")
summary(mdlInit)
# We remove KDM5C (we can look at its interaction later)
mdlRemoveKDM5C <- coxph(Surv(Time, Surv) ~ RT + Age + BIRC3 + CD79B + MALT1 + 
                          JAK3 + HOXA13 + TBL1XR1 + TSHR + GATA3 + PIM1 + 
                          CSF3R + AR,
                        data = df, ties = "exact")
# We need to take care that this exclusion is justified by the partial LRT
# and that it doesn't have a big effect on other covariates
# The partial likelihood ratio test can be carried out using anova
anova(mdlInit, mdlRemoveKDM5C)
# To look for big changes in coefficients, I will write a function to make this
# easier
compareCoeffs <- function(model, submodel) {
  coeffs <- labels(submodel$coefficients)
  df <- data.frame("Coefficient" = coeffs,
                   "Old.value" = NA,
                   "New.value"= NA,
                   "Percentage.change" = NA)
  for (i in 1:nrow(df)) {
    df$Old.value[i] <- model$coefficients[df$Coefficient[i]]
    df$New.value[i] <- submodel$coefficients[df$Coefficient[i]]
    df$Percentage.change[i] <- round(((df$Old.value[i] -  df$New.value[i]) / 
                                        df$Old.value[i]) * 100, 4)
  }
  return(df)
}
compareCoeffs(mdlInit, mdlRemoveKDM5C)
# All looks good
summary(mdlRemoveKDM5C)
# Now remove PIM1
mdlRemovePIM1 <- coxph(Surv(Time, Surv) ~ RT + Age + BIRC3 + CD79B + MALT1 + 
                         JAK3 + HOXA13 + TBL1XR1 + TSHR + GATA3 + CSF3R + AR,
                                         data = df, ties = "exact")
# Again check PLRT and change in covariates
anova(mdlRemoveKDM5C, mdlRemovePIM1)
compareCoeffs(mdlRemoveKDM5C, mdlRemovePIM1)
# All ok.
summary(mdlRemovePIM1)
# Remove GATA3
mdlRemoveGATA3 <- coxph(Surv(Time, Surv) ~ RT + Age + BIRC3 + CD79B + MALT1 + 
                          JAK3 + HOXA13 + TBL1XR1 + TSHR + CSF3R + AR,
                        data = df, ties = "exact")
anova(mdlRemovePIM1, mdlRemoveGATA3)
compareCoeffs(mdlRemovePIM1, mdlRemoveGATA3)
# Some bigger percentage changes but still none bigger than 20%
summary(mdlRemoveGATA3)
# Remove JAK3
mdlRemoveJAK3 <- coxph(Surv(Time, Surv) ~ RT + Age + BIRC3 + CD79B + MALT1 +
                         HOXA13 + TBL1XR1 + TSHR + CSF3R + AR,
                       data = df, ties = "exact")
anova(mdlRemoveGATA3, mdlRemoveJAK3)
compareCoeffs(mdlRemoveGATA3, mdlRemoveJAK3)
# There are some large percentage changes (on BIRC3 and CSF3R) so JAK3 is an 
# important confounder and so we should keep it in the model
# Instead remove TSHR
mdlRemoveTSHR <- coxph(Surv(Time, Surv) ~ RT + Age + BIRC3 + CD79B + MALT1 + 
                          JAK3 + HOXA13 + TBL1XR1 + CSF3R + AR,
                        data = df, ties = "exact")
anova(mdlRemoveGATA3, mdlRemoveTSHR)
compareCoeffs(mdlRemoveGATA3, mdlRemoveTSHR)
# This is ok
summary(mdlRemoveTSHR)
# Remove BIRC3
mdlRemoveBIRC3 <- coxph(Surv(Time, Surv) ~ RT + Age + CD79B + MALT1 + JAK3 + 
                         HOXA13 + TBL1XR1 + CSF3R + AR,
                       data = df, ties = "exact")
anova(mdlRemoveTSHR, mdlRemoveBIRC3)
compareCoeffs(mdlRemoveTSHR, mdlRemoveBIRC3)
# This causes a big change in JAK3, but considering this wasn't significant at 
# the 5% level, we could simply remove it along with CSF3R which
# are confounders and not significant
mdlRemoveConfounders <- coxph(Surv(Time, Surv) ~ RT + Age + CD79B + MALT1 + 
                                HOXA13 + TBL1XR1 + AR,
                              data = df, ties = "exact")
anova(mdlRemoveTSHR, mdlRemoveConfounders)
compareCoeffs(mdlRemoveTSHR, mdlRemoveConfounders)
# This is ok
summary(mdlRemoveConfounders)
# Remove CD79B
mdlRemoveCD79B <- coxph(Surv(Time, Surv) ~ RT + Age + MALT1 + HOXA13 + 
                          TBL1XR1 + AR,
                        data = df, ties = "exact")
anova(mdlRemoveConfounders, mdlRemoveCD79B)
compareCoeffs(mdlRemoveConfounders, mdlRemoveCD79B)
# This is ok
summary(mdlRemoveCD79B)
# Remove Age
mdlRemoveAge <- coxph(Surv(Time, Surv) ~ RT + MALT1 + HOXA13 + TBL1XR1 + AR,
                        data = df, ties = "exact")
anova(mdlRemoveCD79B, mdlRemoveAge)
compareCoeffs(mdlRemoveCD79B, mdlRemoveAge)
# This is ok
summary(mdlRemoveAge)
# Remove AR
mdlRemoveAR <- coxph(Surv(Time, Surv) ~ RT + MALT1 + HOXA13 + TBL1XR1,
                      data = df, ties = "exact")
anova(mdlRemoveAge, mdlRemoveAR)
compareCoeffs(mdlRemoveAge, mdlRemoveAR)
# This is ok
summary(mdlRemoveAR)
# Everything is now significant at the 5% level

## STEP 4 ##
# We want to look at all the other variables we excluded initially to see if 
# they are significant in the current model

excluded <- pvals[pvals$P.value >= 0.2,]$Variable
incPVals <- data.frame(excluded,
                       "Included.P.value" = NA)
for (i in 1:length(excluded)) {
  model <- coxph(Surv(Time, Surv) ~ RT + MALT1 + HOXA13 + 
                   TBL1XR1 + df[,excluded[i]],
                 data = df, ties = "exact")
  test <- anova(mdlRemoveAR, model)
  incPVals$Included.P.value[i] <- test$`P(>|Chi|)`[2]
}
# Let's look if any are more significant when added to the model
incPVals
# Sex is now a highly significant variable so we will include that
# We now have our main effects model:
mdl <- coxph(Surv(Time, Surv) ~ RT + MALT1 + HOXA13 + TBL1XR1 + Sex,
             data = df, ties = "exact")
summary(mdl)

## STEP 6 ##
# We now consider interaction terms. First we will try all of the continuous 
# variables interacting with RT, then with Sex.

CatVars <- colnames(df)[c(4, 6:25)]
# With RT:
interactions <- data.frame(Variable = CatVars,
                           "RT" = NA,
                           "Sex" = NA)
for (i in 1:length(CatVars)) {
  modelRT <- coxph(Surv(Time, Surv) ~ RT*df[,CatVars[i]] + MALT1 + HOXA13 + 
                     TBL1XR1 + Sex,
                   data = df, ties = "exact")
  modelSex <- coxph(Surv(Time, Surv) ~ RT + MALT1 + HOXA13 + 
                      TBL1XR1 + Sex*df[,CatVars[i]],
                    data = df, ties = "exact")
  testRT <- anova(mdl, modelRT)
  testSex <- anova(mdl, modelSex)
  interactions$RT[i] <- testRT$`P(>|Chi|)`[2]
  interactions$Sex[i] <- testSex$`P(>|Chi|)`[2]
}

# To avoid adding too many interaction terms and making interpretation too 
# difficult (and also unnecessarily widening HR CIs), we will only add those
# that are significant at a 1% level. This means adding CREB3L1, CREB3L1.1,
# IDH2, KDM5C (as predicted is task 2), STK11, and TBL1XR1, all interacting
# with RT.
# This gives us the preliminary model (we must check model assumptions before
# making this the final model)

mdlInt <- coxph(Surv(Time, Surv) ~ MALT1 + HOXA13 + TBL1XR1 + Sex +
                  RT*CREB3L1 + RT*CREB3L1.1 + RT*IDH2 + RT*KDM5C + 
                  RT*STK11 + RT*TBL1XR1,
                data = df, ties = "exact")
summary(mdlInt)
# It appears CREB3L1.1 is just a copy of CREB3L1 so we need to delete that
mdlInt <- coxph(Surv(Time, Surv) ~ MALT1 + HOXA13 + TBL1XR1 + Sex +
                  RT*CREB3L1 + RT*IDH2 + RT*KDM5C + RT*STK11 + RT*TBL1XR1,
                data = df, ties = "exact")
summary(mdlInt)

# It turns out when added together, some of the p-values are rather large. 
# Let's remove them one at a time until they are sufficiently small

mdlIntRemoveSTK11 <- coxph(Surv(Time, Surv) ~ MALT1 + HOXA13 + TBL1XR1 + Sex +
                             RT*CREB3L1 + RT*IDH2 + RT*KDM5C + RT*TBL1XR1,
                           data = df, ties = "exact")
summary(mdlIntRemoveSTK11)
# Now remove TBL1XR1
mdlIntRemoveTBL1XR1 <- coxph(Surv(Time, Surv) ~ MALT1 + HOXA13 + TBL1XR1 + Sex +
                               RT*CREB3L1 + RT*IDH2 + RT*KDM5C,
                             data = df, ties = "exact")
summary(mdlIntRemoveTBL1XR1)
# Now remove CREB3L1
mdlIntRemoveCREB3L1 <- coxph(Surv(Time, Surv) ~ MALT1 + HOXA13 + TBL1XR1 + Sex +
                               RT*IDH2 + RT*KDM5C,
                             data = df, ties = "exact")
summary(mdlIntRemoveCREB3L1)
# Everything is now significant at the 10% level

# Now we have a preliminary model and just need to check it meets model 
# assumptions

mdlPrelim <- mdlIntRemoveCREB3L1

## STEP 7 ##
# Assessment of fit
cox.zph(mdlPrelim)
# None of the values are significant at the 5% level, therefore we can assume
# proportional hazards
plot(cox.zph(mdlPrelim))

# We should also look for influential individuals. We will do this by plotting
# influence against MArtingale Residuals
Mres <- residuals(mdlPrelim, type = "martingale")
C1 <- residuals(coxph(Surv(Time, Surv) ~ MALT1 + HOXA13 + TBL1XR1 + Sex +
                        RT*IDH2 + RT*KDM5C,
                      data = df, ties = "efron"), 
                type = "dfbeta")
Vars <- labels(coefficients(mdlPrelim))
Plots <- vector(mode = "list", length = 9)
for (i in 1:9) {
  Plots[[i]] <- local({
    i <- i
    p1 <- ggplot() + geom_point(aes(x = Mres, y = C1[, i])) + 
      xlab("Martingale Residuals") + ylab(paste0(Vars[i], " Influence"))
    print(p1)
  })
}
require(gridExtra)
grid.arrange(Plots[[1]], Plots[[2]], Plots[[3]], Plots[[4]], Plots[[5]],
             Plots[[6]], Plots[[7]], Plots[[8]], Plots[[9]],
             ncol = 3, nrow = 3)
# There are two individuals with very low Martingale residuals but they don't 
# to exert any undue influence. There are a few outliers in some of the plots
# but no points are consistently out of place so this looks fine.

### Interpretation
finalmodel <- mdlPrelim
summary(finalmodel)

# Create a table of HRs and CIs and interpret them
# Interpret the condordance

Finally we can check whether we could fit a fully parametric model to these covariates (i.e one that assumes constant hazards). This can be done by producing the following plot:
  ```{r, warning = FALSE}
# Calculate Kaplan-Meier estimates of the survival function
fit <- survfit(Surv(Time, Surv) ~ 1, data = df)
# Create plot
ggplot(data.frame(surv=fit$surv, time=fit$time),
       aes(log(time),log(-log(surv)))) + geom_point() +
  geom_smooth(aes(log(time), log(-log(surv))), color="red", method = "lm")

```
This plot fairly closely follows a straight line, so the constant hazards assumption is met. This means we can fit a Weibull model. We should also check that Weibell over Exponential is justified by performing a Wald test on the scale parameter. In this case the inclusion is justified.
```{r}
mdlParam <- survreg(Surv(Time, Surv) ~ MALT1 + HOXA13 + Sex + RT*IDH2 + RT*KDM5C,
                    data = df, dist = "weibull")
summary(mdlParam)
```


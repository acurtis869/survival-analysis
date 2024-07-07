## Task 2

library(tidyverse)
library(survival)
library(survminer)

df <- read.delim("Practical2Data.tsv") %>%
  filter(RT %in% c("YES", "NO"))
head(df)

# Create new column in df to group levels of KDM5C
range(df$KDM5C)
# Let's create 4 categories: 1-2, 2-3, 3-4, 4-5
df <- df %>%
  mutate(KDM5C_cat = cut(KDM5C, breaks = c(-Inf, median(df$KDM5C), Inf), 
                         labels = c("Low", "High")))
# Now perform a basic log-rank test
survdiff(Surv(Time, Surv) ~ KDM5C_cat, data = dfK)
# What about with 3 categories?
df <- df %>%
  mutate(KDM5C_cat2 = cut(KDM5C, breaks = c(-Inf, 
                                            quantile(df$KDM5C, probs = 1/3),
                                            quantile(df$KDM5C, probs = 2/3), 
                                            Inf), 
                          labels = c("Low", "Medium", "High")))
survdiff(Surv(Time, Surv) ~ KDM5C_cat2, data = df)
# Still not significant. What about if we introduce the trend?
# Use weights -1, 0, 1:
diffobject <- survdiff(Surv(Time, Surv) ~ KDM5C_cat2, data = df)
# Calculate weighted difference between observed and expected numbers
W1 <- c(-1, 0, 1) * (diffobject$obs - diffobject$exp)
# Calculate weighted difference of the mean
W2 <- c(-1 ,0 ,1) * diffobject$exp
# Calculate weighted estimate of the variance
W3 <- c(-1, 0, 1)^2 * diffobject$exp
# Calculate VT from the formula
VT <- sum(W3) - sum(W2)^2/sum(diffobject$exp)
# Calculate test statistic
X2T<-sum(W1)^2/VT
X2T
# Obtain p-value
1 - pchisq(X2T, df = 1)
# There is no evidence for a trend

# We should check whether the data are skewed around Age or Sex in case 
# stratified tests might be required. We can do this by conisdering the cross-
# tabulation
# First sort Age into high medium and low
hist(df$Age)
df <- df %>%
  mutate(Age_cat = cut(Age, breaks = c(-Inf, 
                                       quantile(df$Age, probs = 1/3),
                                       quantile(df$Age, probs = 2/3), 
                                       Inf), 
                       labels = c("Low", "Medium", "High")))
table(df$KDM5C_cat2, df$Age_cat)
# This looks evenly spread
table(df$KDM5C_cat2, df$Sex)
# This suggests it might be worth stratifying over sex
survdiff(Surv(Time, Surv) ~ KDM5C_cat2 + strata(Sex), data = df)
# There is still no strong evidence to suggest KDM5C affects survival
table(df$KDM5C_cat2, df$RT)
# This is fairly balanced within the group

# What about factorial design?

survdiff(Surv(Time, Surv) ~ KDM5C_cat + Sex, data = df)
# Not much evidence
survdiff(Surv(Time, Surv) ~ KDM5C_cat + Age_cat, data = df)
# Less evidence
survdiff(Surv(Time, Surv) ~ KDM5C_cat + RT, data = df)
# Strong evidence for difference. This however is to be expected, since
# we already know there is a strong association with RT. It does, however look
# interesting that RT seems to be somewhat more effective for patients with low
# levels of KDM5C rather than high (compare difference between exp. and obs.
# deaths when RT = YES). This relationship is flipped when there was no RT 
# given. This suggests there may be some kind of interaction that could be 
# investigated further using parametric methods.




































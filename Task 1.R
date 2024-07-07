## Task 1

library(tidyverse)
library(survival)
library(survminer)

df <- read.delim("Practical2Data.tsv") %>%
  filter(RT %in% c("YES", "NO"))
head(df)

# Investigate the association between radiation therapy and survival.
# We can do this with a Log-rank test. This is non-parametric, and just 
# assumes the data come from two iid groups.

survdiff(Surv(Time, Surv) ~ RT, data = df)

# The non-parametric test uses the standard chi-squared statistic
# The p-value is calculated by
X_2 <- 59.5 + 2
1 - pchisq(X_2, df = 1)
# This suggests there is a clear difference
# The parametric test (which assumes a hyperheometric distribution) has
# the following result
1 - pchisq(63.2, df = 1)
# This is again a very small p-value so suggests there is a significant 
# difference between the two groups

# Illustrate this using a plot

# First create a survival fit object (i.e a high resolution life table)
# This uses the Kaplan-Meier formula
fit <- survfit(Surv(Time, Surv) ~ RT, data = df, conf.type = "log")
# Now create a plot using ggsurvplot
ggsurvplot(fit, data = df, surv.median.line = "hv")
# It is pretty clear from this that they are distinct. It is also clear that
# RT increases survival time.

# Report confidence intervals for median survival in each radiation therapy 
# category

# Confidence intervals are calculated using log confidence intervals detailed by
# Collett (1994), and convention by Hosmer and Lemeshow
fit
# Plot these?
# Create dataframe
CIs <- data.frame(RT = c("NO", "YES"),
                  Median = c(83, 448),
                  UpperCL = c(124, 489),
                  LowerCL = c(33, 399))
qplot(x = RT, y = Median, data = CIs) +
  geom_errorbar(aes(ymin = LowerCL, ymax = UpperCL, width = 0.15))
# These confidence clearly to not overlap so once again we have strong evidence
# that there is a difference between the survivals of these two groups.

# We can also look at this relationship by calculating the hazard ratio and
# its confidence interval

HR <- (20 / 4.2) / (109 / 124.8)

Var_log_HR <- 1 / 4.2 + 1 / 124.8

round(exp(log(HR + c(-1, 0, 1) * pnorm(0.975) * sqrt(Var_log_HR))), 2)
# This gives strong evidence that the hazard rates of whether RT is recieved 
# are significantly different (as the hazard ratio doesn't include 1). 
# Furthermore, it suggests that not having RT increases the hazard rate 
# by a factor of somewhere between 5 and 6, which is quite considerable.


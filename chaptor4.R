#install.packages("survival")
#install.packages("aod")
library(aod)
library(survival)
library(Hmisc)

data(cancer, package="survival")

lung2 <- lung[complete.cases(lung),]
# Fit a Cox PH model with age and sex as covariates - different ties
cox_model_efron <- coxph(Surv(time, status) ~ age + sex, data = lung2,
                   ties=c("efron"))
summary(cox_model_efron)

cox_model_breslow <- coxph(Surv(time, status) ~ age + sex, data = lung2,
                   ties=c("breslow"))
summary(cox_model_breslow)

cox_model_exact <- coxph(Surv(time, status) ~ age + sex, data = lung2,
                   ties=c("exact"))
summary(cox_model_exact)

# continuous ph.ecog
# Fit the full model
cox_full <- coxph(Surv(time, status) ~ age + sex + ph.ecog, data = lung2)

# Fit the reduced model (nested)
cox_reduced <- coxph(Surv(time, status) ~ age + sex, data = lung2)

# Compare the models using a likelihood ratio test
anova(cox_reduced, cox_full)

# categorical ph.ecog
cox_full2 <- coxph(Surv(time, status) ~ age + sex + I(ph.ecog==1) + I(ph.ecog==2) + I(ph.ecog==3), data = lung2)

# wald
wald.test(b = coef(cox_full2), Sigma = vcov(cox_full2), Terms = 3:5)

# LRT
anova(cox_reduced, cox_full2)


surv_fit <- survfit(cox_full)
summary(surv_fit)

# Plot the survival curve
plot(surv_fit, xlab = "Time", ylab = "Survival Probability", main = "Survival Function Based on Cox Model")

# Predict survival curves for males (sex = 1) and females (sex = 2)
# For males (sex = 1)
newdata_male <- data.frame(age = median(lung2$age), sex = 1, ph.ecog = median(lung2$ph.ecog, na.rm = TRUE))
surv_fit_male <- survfit(cox_full, newdata = newdata_male)

# For females (sex = 2)
newdata_female <- data.frame(age = median(lung2$age), sex = 2, ph.ecog = median(lung2$ph.ecog, na.rm = TRUE))
surv_fit_female <- survfit(cox_full, newdata = newdata_female)

# Plot survival curves for both subgroups
plot(surv_fit_male, col = "blue", lty = 1, xlab = "Time", ylab = "Survival Probability",
     main = "Survival Curves by Sex", ylim = c(0, 1))
lines(surv_fit_female, col = "red", lty = 2)

# Add legend
legend("bottomleft", legend = c("Male", "Female"), col = c("blue", "red"), lty = c(1, 2))


# Fit the Kaplan-Meier estimator for subgroups based on sex
km_fit <- survfit(Surv(time, status) ~ sex, data = lung2)

# Plot the Kaplan-Meier survival curves for males and females
plot(km_fit, col = c("blue", "red"), lty = 1:2, xlab = "Time", ylab = "Survival Probability",
     main = "Kaplan-Meier Survival Curves by Sex", ylim = c(0, 1))

# Add legend
legend("bottomleft", legend = c("Male", "Female"), col = c("blue", "red"), lty = 1:2)

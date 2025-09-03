
################Chapter 1################

rasbc
names(rasbc)


library(tidyverse)
library(survival)
library(survminer)
library(rms)



ref=rasbc$V001
age=rasbc$V015
menop=rasbc$V016
presite=rasbc$V017
side=rasbc$V018
maxgr=rasbc$V020
nodes=rasbc$V023
t=rasbc$V026
n=rasbc$V027
m=rasbc$V028
stage=rasbc$V029
dateprisur=rasbc$V034
typeprisur=rasbc$V035
radio=rasbc$V040
histo=rasbc$V044
involve=rasbc$V045
remove=rasbc$V046
pagr=rasbc$V052
grade=rasbc$V053
estr=rasbc$V054
progest=rasbc$V057
datelast=rasbc$V060
curr=rasbc$V061
datere=rasbc$V075
s=rasbc$survival_time_days
w=rasbc$censored

nrat=involve/remove

boxplot(age~w,main="age")
table(menop,w)
round(prop.table(table(menop,w), margin = 1), 2) 

table(maxgr,w)
round(prop.table(table(maxgr,w), margin = 1), 2) 

table(grade,w)
round(prop.table(table(grade,w), margin = 1), 2) 





################Chapter 2################
surv_obj <- Surv(time = s, event = w)
fit <- survfit(surv_obj ~ 1, data = rasbc)
ggsurvplot(fit, 
           ggtheme = theme_minimal(),
           title = "Overall Kaplan-Meier Survival Curve",
           xlab = "Time", 
           ylab = "Survival Probability")


# categorising age into 3 categories
agegr=rep(0,572)
for(i in 1:572)
{if(age[i]<=44)agegr[i]=1
if(age[i]>=45&age[i]<=64)agegr[i]=2
if(age[i]>=65)agegr[i]=3
}


fit_age <- survfit(surv_obj ~ agegr, data = rasbc)
ggsurvplot(fit_age, 
           legend.title = "Age Group",
           ggtheme = theme_minimal(),
           title = "Kaplan-Meier Curve by Age Group",
           xlab = "Time", 
           ylab = "Survival Probability")


fit_menop <- survfit(surv_obj ~ factor(menop), data = rasbc)
ggsurvplot(fit_menop,
           legend.title = "Menopausal Status",
           ggtheme = theme_minimal(),
           title = "Kaplan-Meier Curve by Menopausal Status",
           xlab = "Time", 
           ylab = "Survival Probability")


fit_maxgr <- survfit(surv_obj ~ factor(maxgr), data = rasbc)
ggsurvplot(fit_maxgr,
           legend.title = "Maximum Tumour Diameter Group",
           ggtheme = theme_minimal(),
           title = "Kaplan-Meier Curve by Maximum Tumour Diameter",
           xlab = "Time", 
           ylab = "Survival Probability")

fit_grade <- survfit(surv_obj ~ factor(grade), data = rasbc)
ggsurvplot(fit_grade,
           legend.title = "Tumour Grade",
           ggtheme = theme_minimal(),
           title = "Kaplan-Meier Curve by Tumour Grade",
           xlab = "Time", 
           ylab = "Survival Probability")



##log-rank
logrank_age <- survdiff(surv_obj ~ agegr, data = rasbc)
logrank_age

logrank_menop <- survdiff(surv_obj ~ factor(menop), data = rasbc)
logrank_menop

logrank_maxgr <- survdiff(surv_obj ~ factor(maxgr), data = rasbc)
logrank_maxgr

logrank_grade <- survdiff(surv_obj ~ factor(grade), data = rasbc)
logrank_grade






################Chapter 3################

##group the variables
# a categorisation of nratio into 7 groups by censored and uncensored
nr=rep(0,572)
for(i in 1:572)
{if(nrat[i]==0)nr[i]=1
if(nrat[i]>0&nrat[i]<=0.2)nr[i]=2
if(nrat[i]>0.2&nrat[i]<=0.4)nr[i]=3
if(nrat[i]>0.4&nrat[i]<=0.6)nr[i]=4
if(nrat[i]>0.6&nrat[i]<=0.8)nr[i]=5
if(nrat[i]>0.8&nrat[i]<=1)nr[i]=6
if(nrat[i]==Inf)nr[i]=7
}
# nr
table(nr,w)

# categorising nratio into 3 categories

nrr=rep(0,572)
for(i in 1:572)
{if(nrat[i]==0 || nrat[i]==Inf )nrr[i]=1
if(nrat[i]>0&nrat[i]<=0.4)nrr[i]=2
if(nrat[i]>0.4&nrat[i]<=1)nrr[i]=3
}

table(nrr,w)

# categorising age into 3 categories
agegr=rep(0,572)
for(i in 1:572)
{if(age[i]<=44)agegr[i]=1
if(age[i]>=45&age[i]<=64)agegr[i]=2
if(age[i]>=65)agegr[i]=3
}

table(agegr,w)


#################################

# combining categories in some of the variables

# 'MENOPAUSAL STATUS'
table(menop,w)

# categories 1 and 2 comnined into a new category 1 corresponding to last menstruation =< 3 years
# the old category 3 becomes the new category 2 corresponding to last menstruation > 3 years
# category 9, unknown will remain unchanged

menop.new=menop
for(i in 1:572){
  if(menop[i]==2) menop.new[i]=1
  if(menop[i]==3) menop.new[i]=2
}

table(menop.new, w)


# 'PREDOMINANT SITE in the Breast'
table(presite,w)

# we will combine the first 4 categories into a new bew category 1 named "non-central sites"
# category 5 becomes the new category 2 called "subareola (middle)"
# category 9, site unknown, will remain the same.

presite.new=presite
for(i in 1:572){
  if(presite[i]==1 || presite[i]==2 || presite[i]==3 || presite[i]==4) presite.new[i]=1
  if(presite[i]==5) presite.new[i]=2
}

table(presite.new, w)



# 'TYPE OF PRIMARY SURGERY'
table(typeprisur,w)

typeprisurfit=survfit(Surv(s, w)~typeprisur,type="kaplan-meier",error="greenwood",conf.type="plain")
plot(typeprisurfit,col=1:6,lwd=2,lty=1:6,xlab="Survival days",ylab="Survival Probability",
     main="Type of Primary Surgery")
legend("bottomleft",legend=c("3(excision biopsy)","4(simple mastectomy)",
                             "6(wide local excision with axillary nodal clearance)","7(mastectomy with axillary nodal clearance)",
                             "8(quadratectomy with axillary nodal clearance)","99(unknown)"),col=1:6,lwd=2,lty=1:6)

# the unknown (99) and excision biopsy (3) categories have very similar survival curves 
# and the two groups are to be combined into a new category 3.

typeps.new = typeprisur
for(i in 1:572){
  if(typeprisur[i]==99) typeps.new[i]=3 
}

table(typeps.new, w)


# 'HISTOLOGY'

table(histo,w)

# categories 3, 4, 5 are to be combined into a new category 3 named "Other"
# the existing categories 1 and 2 are to remain unchanged.

histo.new=histo
for(i in 1:572){
  if(histo[i]== 4 || histo[i]==5) histo.new[i]=3
}

table(histo.new, w)


# 'TUMOUR GRADE'


table(grade,w)

# combine categories 1 and 2 into a new category 2 which will now comprise grades 1 and 2
# categories 3 and 9 are to be unchanged

tgrade.new = grade
for(i in 1:572){
  if(grade[i]==1) tgrade.new[i]= 2
}

table(tgrade.new, w)


# RECEPTOR LEVEL AS A NEW VARIABLE - A COMBINATION OF
# ESTROGEN AND PROGESTERONE LEVELS


bcs=Surv(time=rasbc[,25], event=rasbc[,26])

# we can estimate all the paarmeters in Cox PH models just containg a simgle covariate of either estrogen
# or progesterone levels. ie.

coxph(bcs~factor(estr))  #  estrogen category 0 is the reference category

coxph(bcs~factor(progest))  # progest category 0 is the reference category

# but when we include both factors in a single additive linear predictor
# we cannot estimate the effect for progest level 9.  ie.

coxph(bcs~factor(estr)+factor(progest))


# To see where the problem is tabulate these two factors

addmargins(table(progest, estr, w))
addmargins(table( progest, estr))
addmargins(table( progest, w))
addmargins(table( estr, w))

# the same 289 cases have both estr=9 and progest=9.  
# ie. both are missing values for these cases

# hence in the PH model containg both factors only the effect for estrogen level 9 is estimated

# we therefore combine estrogen and progesterone into a new varaible "recept" corresponding 
# to receptor level with five categories.  ie.

recept=0L
for(i in 1:572){
  if(estr[i]==0 & progest[i]==0) recept[i]=1
  if(estr[i]==0 & progest[i]==1) recept[i]=2
  if(estr[i]==1 & progest[i]==0) recept[i]=3
  if(estr[i]==1 & progest[i]==1) recept[i]=4
  if(estr[i]==9 & progest[i]==9) recept[i]=9
}

addmargins(table(recept, w))

coxph(bcs~factor(recept)) # recept level 1 is the reference level

# parameter estimation works fine now

# MANCHESTER STAGE - only include this stage variable in Cox PH models.  ie. stage=rasbc$V029

# Ignore the N, M, and T stage factors



##cox

cox.full <- coxph(bcs ~ 
                    factor(agegr) +
                    factor(menop.new) +
                    factor(presite.new) +
                    factor(side) +
                    factor(maxgr) +
                    factor(nodes) +
                    factor(stage) +
                    factor(radio) +
                    factor(histo.new) +
                    factor(nrr) +
                    factor(pagr) +
                    factor(tgrade.new) +
                    factor(recept) +
                    factor(typeps.new),
                  data =rasbc) 
summary(cox.full)



##diagnostics

# Cox-Snell

cox.snell <- resid(cox.full, type = "martingale")
cox.snell_res <- w - cox.snell

surv.coxsnell <- Surv(cox.snell_res, w)
fit.coxsnell <- survfit(surv.coxsnell ~ 1)

surv.coxsnell
cumhaz.coxsnell <- -log(fit.coxsnell$surv)

plot(fit.coxsnell$time, cumhaz.coxsnell, 
     type = "s", 
     xlab = "Cox-Snell Residuals", 
     ylab = "Cumulative Hazard",
     main = "Cox-Snell Residuals Diagnostic Plot")

# y=x line
abline(0, 1, col = "red", lwd = 2)

legend("bottomright", legend = c("Observed", "Expected"), 
       col = c("black", "red"), lty = c(1, 1), lwd = c(1, 2))



# Martingale
martingale_res <- resid(cox.full, type = "martingale")


time <- rasbc$survival_time_days
status <- rasbc$censored  # 1 = died / event, 0 = censored


color_vec <- ifelse(status == 1, "red", "blue")
label_vec <- ifelse(status == 1, "Died / Liver failure", "Censored")


plot(time, martingale_res,
     type = "n",  
     xlab = "Time (Days)",
     ylab = "Martingale Residuals",
     main = "Martingale Residuals vs Time (by Event Status)")


points(time[status == 0], martingale_res[status == 0],
       col = "blue", pch = 1)  # Censored
points(time[status == 1], martingale_res[status == 1],
       col = "red", pch = 16)  # Died
abline(h = 0, col = "gray40", lty = 2)

legend("bottomleft",
       legend = c("Censored", "Uncensored", "Zero Line"),
       col = c("blue", "red", "gray40"),
       pch = c(1, 16, NA),
       lty = c(NA, NA, 2),
       lwd = c(NA, NA, 1),
       bty = "n")


df_res <- data.frame(
  index = 1:length(martingale_res),
  time = time,
  status = status,
  residual = martingale_res
)


min_blue <- df_res[df_res$status == 0, ]
min_blue <- min_blue[order(min_blue$residual), ][1:5, ]


min_red <- df_res[df_res$status == 1, ]
min_red <- min_red[order(min_red$residual), ][1:5, ]


cat("---- the most negative blue ones ----\n")
print(min_blue)

cat("\n---- the most negative red ones ----\n")
print(min_red)




# Deviance residuals
deviance_res <- resid(cox.full, type = "deviance")

time <- rasbc$survival_time_days

status <- rasbc$censored

plot(time, deviance_res,
     type = "n",
     xlab = "Time (Days)",
     ylab = "Deviance Residuals",
     main = "Deviance Residuals vs Time")

points(time[status == 0], deviance_res[status == 0],
       col = "blue", pch = 1)
points(time[status == 1], deviance_res[status == 1],
       col = "red", pch = 16)

lines(lowess(time, deviance_res), col = "darkgreen", lwd = 2)

abline(h = 0, lty = 2, col = "gray40")

legend("bottomleft",
       legend = c("Censored", "Uncensored", "LOWESS Smooth", "Zero Line"),
       col = c("blue", "red", "darkgreen", "gray40"),
       pch = c(1, 16, NA, NA),
       lty = c(NA, NA, 1, 2),
       lwd = c(NA, NA, 2, 1),
       bty = "n")
qqnorm(deviance_res, main = "QQ Plot of Deviance Residuals")
qqline(deviance_res, col = "red", lwd = 2)




# Schoenfeld residuals
schoenfeld_res <- resid(cox.full, type = "schoenfeld")


ph_test <- cox.zph(cox.full)
print(ph_test)

par(mfrow = c(3, 5))  

for(i in 1:ncol(schoenfeld_res)){
  plot(ph_test, var = i, main = colnames(schoenfeld_res)[i])
}




##Model Reduction
# Model Reduction using Stepwise Selection
library(MASS)

# Stepwise selection methods
step_backward <- step(cox.full, direction = "backward")
step_both <- step(cox.full, direction = "both")

# Compare models
cat("=== Model Comparison ===\n")
cat("Full model AIC:", AIC(cox.full), "\n")
cat("Backward selection AIC:", AIC(step_backward), "\n")
cat("Both directions AIC:", AIC(step_both), "\n")

# Use the best model (lowest AIC)
final_model <- step_both

# Show selected variables
cat("\n=== Selected Variables ===\n")
print(names(final_model$coefficients))

# Model diagnostics
cat("\n=== Proportional Hazards Test ===\n")
ph_test_final <- cox.zph(final_model)
print(ph_test_final)


# Likelihood ratio test
if (length(final_model$coefficients) < length(cox.full$coefficients)) {
  lrt <- anova(final_model, cox.full, test = "Chisq")
  cat("\n=== Likelihood Ratio Test ===\n")
  print(lrt)
}


summary(final_model)



cox.reduced <- coxph(bcs ~ 
                    factor(maxgr) +
                    factor(nodes) +
                    factor(radio) +
                    factor(histo.new) +
                    factor(nrr) +
                    factor(pagr) +
                    factor(tgrade.new) +
                    factor(recept),
                  data =rasbc) 
summary(cox.reduced)


# Cox-Snell residuals
cox.snell_reduced <- resid(cox.reduced, type = "martingale")
cox.snell_res_reduced <- rasbc$censored - cox.snell_reduced

surv.coxsnell_reduced <- Surv(cox.snell_res_reduced, rasbc$censored)
fit.coxsnell_reduced <- survfit(surv.coxsnell_reduced ~ 1)

cumhaz.coxsnell_reduced <- -log(fit.coxsnell_reduced$surv)

plot(fit.coxsnell_reduced$time, cumhaz.coxsnell_reduced,
     type = "s",
     xlab = "Cox-Snell Residuals",
     ylab = "Cumulative Hazard",
     main = "Cox-Snell Residuals Diagnostic Plot (Reduced Model)")

abline(0, 1, col = "red", lwd = 2)

legend("bottomright", legend = c("Observed", "Expected"),
       col = c("black", "red"), lty = c(1, 1), lwd = c(1, 2))


# Martingale
martingale_res_reduced <- resid(cox.reduced, type = "martingale")

time <- rasbc$survival_time_days
status <- rasbc$censored

plot(time, martingale_res_reduced,
     type = "n",
     xlab = "Time (Days)",
     ylab = "Martingale Residuals",
     main = "Martingale Residuals vs Time (Reduced Model)")

points(time[status == 0], martingale_res_reduced[status == 0],
       col = "blue", pch = 1)
points(time[status == 1], martingale_res_reduced[status == 1],
       col = "red", pch = 16)
abline(h = 0, col = "gray40", lty = 2)

legend("bottomleft",
       legend = c("Censored", "Uncensored", "Zero Line"),
       col = c("blue", "red", "gray40"),
       pch = c(1, 16, NA),
       lty = c(NA, NA, 2),
       lwd = c(NA, NA, 1),
       bty = "n")

df_res_reduced <- data.frame(
  index = seq_along(martingale_res_reduced),
  time = time,
  status = status,
  residual = martingale_res_reduced
)

# Entire sample: 5 most negative residuals 
worst5_overall <- df_res_reduced[order(df_res_reduced$residual), ][1:5, ]
print(worst5_overall)



# Deviance Residuals
deviance_res_reduced <- resid(cox.reduced, type = "deviance")

plot(time, deviance_res_reduced,
     type = "n",
     xlab = "Time (Days)",
     ylab = "Deviance Residuals",
     main = "Deviance Residuals vs Time (Reduced Model)")

points(time[status == 0], deviance_res_reduced[status == 0],
       col = "blue", pch = 1)
points(time[status == 1], deviance_res_reduced[status == 1],
       col = "red", pch = 16)
lines(lowess(time, deviance_res_reduced), col = "darkgreen", lwd = 2)
abline(h = 0, lty = 2, col = "gray40")

legend("bottomleft",
       legend = c("Censored", "Uncensored", "LOWESS Smooth", "Zero Line"),
       col = c("blue", "red", "darkgreen", "gray40"),
       pch = c(1, 16, NA, NA),
       lty = c(NA, NA, 1, 2),
       lwd = c(NA, NA, 2, 1),
       bty = "n")

# QQ plot
qqnorm(deviance_res_reduced, main = "QQ Plot of Deviance Residuals (Reduced Model)")
qqline(deviance_res_reduced, col = "red", lwd = 2)


# Schoenfeld residuals 
schoenfeld_res_reduced <- resid(cox.reduced, type = "schoenfeld")

ph_test1 <- cox.zph(cox.reduced)
print(ph_test1)
par(mfrow = c(3,3))  

for(i in 1:ncol(schoenfeld_res_reduced)){
  plot(ph_test1, var = i, main = colnames(schoenfeld_res_reduced)[i])
}







################Chapter 4################

library(survival)
library(rpart.plot) 
library(randomForest)
library(party)
library(ggplot2)
library(survminer)
library(tidyverse)
library(caret)


# ============================================================================
# SECTION 1: LIBRARY LOADING AND DATA PREPARATION
# ============================================================================

# Response: survival time and event indicator (1=event, 0=censored)
# Predictors: all full model variables (factors)

tree_data <- data.frame(
  surv_time   = s,
  censored    = w,             # 1 = event, 0 = censored
  agegr       = factor(agegr),
  menop_new   = factor(menop.new),
  presite_new = factor(presite.new),
  side        = factor(side),
  maxgr       = factor(maxgr),
  nodes       = factor(nodes),
  stage       = factor(stage),
  radio       = factor(radio),
  histo_new   = factor(histo.new),
  nrr         = factor(nrr),
  pagr        = factor(pagr),
  tgrade_new  = factor(tgrade.new),
  recept      = factor(recept),
  typeps_new  = factor(typeps.new)
)

# Clean 
tree_data <- na.omit(tree_data)
cat("Dataset:", nrow(tree_data), "rows,", ncol(tree_data), "columns\n\n")
print(summary(tree_data))

# Train/Test split (70/30) 
set.seed(135)
n <- nrow(tree_data)
test_id <- sample(seq_len(n), size = round(0.30 * n), replace = FALSE)
train <- tree_data[-test_id, ]
test <- tree_data[ test_id, ]

cat("\nTrain:", nrow(train), " | Test:", nrow(test), "\n")



# ============================================================================
# SECTION 2: FITTING PROCEDURE (rpart, method="exp")
# ============================================================================

surv_train <- Surv(train$surv_time, train$censored)


# Fit initial generous tree
fit_exp <- rpart(
  formula = surv_train ~ . - surv_time - censored,
  data    = train,
  method  = "exp",    # survival deviance-based splitting
  control = rpart.control(
    minsplit  = 20,
    minbucket = 10,
    cp        = 0.001,  # liberal initial cp
    maxdepth  = 30
  )
)


cat("\n=== FITTING PROCEDURE ===\n")
cat("Initial controls: minsplit =", fit_exp$control$minsplit,
    "| minbucket =", fit_exp$control$minbucket,
    "| cp =", fit_exp$control$cp, "\n\n")


# Get CP table
cat("=== Cross-Validation Results ===\n")
cp_tab <- fit_exp$cptable
print(cp_tab)

# Plot CV curve
plotcp(fit_exp, main = "Cross-Validation Error vs Complexity Parameter")

# Plot unpruned survival tree
rpart.plot(fit_exp,main = "Unpruned Survival Tree (method = 'exp')",
type = 4, extra = 0, fallen.leaves = TRUE, cex = 0.7)





# ============================================================================
# SECTION 3: MODEL SELECTION - Minimum-CV Tree vs 1-SE Tree
# ============================================================================


# Find minimum CV error
min_row <- which.min(cp_tab[, "xerror"])
min_xerror <- cp_tab[min_row, "xerror"]
min_xstd <- cp_tab[min_row, "xstd"]
min_cp <- cp_tab[min_row, "CP"]
min_nsplit <- cp_tab[min_row, "nsplit"]

cat("Minimum CV error at nsplit =", min_nsplit, ", CP =", round(min_cp, 6), 
    ", xerror =", round(min_xerror, 4), "\n")

# Apply 1-SE rule
se_threshold <- min_xerror + min_xstd
candidates <- which(cp_tab[, "xerror"] <= se_threshold)
se_row <- min(candidates)  
se_cp <- cp_tab[se_row, "CP"]
se_nsplit <- cp_tab[se_row, "nsplit"]

cat("1-SE rule selects nsplit =", se_nsplit, ", CP =", round(se_cp, 6), 
    ", xerror =", round(cp_tab[se_row, "xerror"], 4), "\n\n")

# Create both subtrees
tree_min_cv <- prune(fit_exp, cp = min_cp)
tree_1se <- prune(fit_exp, cp = se_cp)





# ============================================================================
# SECTION 4: SELECTED SUBTREES AND THEIR STRUCTURES
# ============================================================================

# Plot minimum CV error tree
rpart.plot(
  tree_min_cv,
  main = paste("Minimum CV Error Subtree (nsplit =", min_nsplit, ")"),
  type = 4, extra = 0, fallen.leaves = TRUE,
  cex = 0.8
)

# Plot 1-SE tree
rpart.plot(
  tree_1se,
  main = paste("1-SE Rule Subtree (nsplit =", se_nsplit, ")"),
  type = 4, extra = 0, fallen.leaves = TRUE,
  cex = 0.8
)



# ============================================================================
# SECTION 5: LEAF-WISE SURVIVAL CURVES (Kaplan-Meier)
# ============================================================================

library(survival)
library(survminer)

# Min-CV Tree
# Get predicted risk values (represent different leaf nodes)
risk_min_cv <- predict(tree_min_cv, newdata = test)

# Create risk groups based on unique predicted values
risk_groups_min_cv <- as.factor(risk_min_cv)

# Create survival object
surv_obj <- Surv(test$surv_time, test$censored)

# Create data frame for analysis
data_min_cv <- data.frame(
  time = test$surv_time,
  status = test$censored,
  group = risk_groups_min_cv
)

# Fit Kaplan-Meier curves
km_min_cv <- survfit(Surv(time, status) ~ group, data = data_min_cv)

# Plot
ggsurvplot(
  km_min_cv,
  data = data_min_cv,
  title = "Min-CV Tree: Survival by Risk Groups",
  xlab = "Time (days)",
  ylab = "Survival Probability",
  legend.title = "Risk Group",
  risk.table = TRUE,
  conf.int = TRUE
)

# Print group information
cat("Min-CV Tree:\n")
cat("Number of risk groups:", length(unique(risk_min_cv)), "\n")
print(table(risk_groups_min_cv))




# 1-SE Tree
# Get predicted risk values
risk_1se <- predict(tree_1se, newdata = test)

# Create risk groups
risk_groups_1se <- as.factor(risk_1se)

# Create data frame for analysis
data_1se <- data.frame(
  time = test$surv_time,
  status = test$censored,
  group = risk_groups_1se
)

# Fit Kaplan-Meier curves
km_1se <- survfit(Surv(time, status) ~ group, data = data_1se)

# Plot
ggsurvplot(
  km_1se,
  data = data_1se,
  title = "1-SE Tree: Survival by Risk Groups",
  xlab = "Time (months)",
  ylab = "Survival Probability",
  legend.title = "Risk Group",
  risk.table = TRUE,
  conf.int = TRUE
)

# Print group information
cat("\n1-SE Tree:\n")
cat("Number of risk groups:", length(unique(risk_1se)), "\n")
print(table(risk_groups_1se))




# ============================================================================
# SECTION 6: VARIABLE IMPORTANCE
# ============================================================================

# Extract variable importance from the minimum CV tree 
var_imp_cv <- tree_min_cv$variable.importance
var_imp_cv <- sort(var_imp_cv, decreasing = TRUE)

# Plot variable importance
barplot(var_imp_cv, 
        las = 2, 
        main = "Variable Importance (Min-CV Tree)",
        ylab = "Importance Score",
        cex.names = 0.8)

cat("Top 5 most important variables:\n")
print(head(var_imp_cv, 5))



# Extract variable importance from the 1-SE tree 
var_imp_1se <- tree_1se$variable.importance
var_imp_1se <- sort(var_imp_1se, decreasing = TRUE)

# Plot variable importance
barplot(var_imp_1se, 
        las = 2, 
        main = "Variable Importance (1-SE Tree)",
        ylab = "Importance Score",
        cex.names = 0.8)

cat("Top 5 most important variables:\n")
print(head(var_imp_1se, 5))




# ============================================================================
# SECTION 7: RANDOM FOREST
# ============================================================================


# Random Forest (for survival time prediction)

library(randomForestSRC)
rf_model <- rfsrc(Surv(surv_time, censored) ~ ., 
                  data = train, 
                  ntree = 500,
                  importance = TRUE)

print(rf_model)

# RF Variable importance
importance_scores <- rf_model$importance
print(importance_scores)

rf_vimp <- vimp(rf_model)
print(rf_vimp)

# Covert to data frame
vimp_df <- data.frame(
  variable = names(rf_vimp$importance),
  importance = rf_vimp$importance
)

# Sort by importance
vimp_df <- vimp_df[order(vimp_df$importance, decreasing = TRUE), ]

ggplot(vimp_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Variable Importance (Random Survival Forest)",
    x = "Variable",
    y = "Importance"
  ) +
  theme_minimal()




# ============================================================================
# SECTION 8: COX MODEL AND THREE-MODEL COMPARISON
# ============================================================================

# Cox model using training set
cox_model <- coxph(Surv(surv_time, censored) ~ 
                     agegr + menop_new + presite_new + side +
                     maxgr + nodes + stage + radio + histo_new + 
                     nrr + pagr + tgrade_new + recept + typeps_new,
                   data = train)

print(summary(cox_model))


# ============================================================================
# SECTION 9: MODEL EVALUATION AND COMPARISON
# ============================================================================

# Generate predictions
cox_pred <- predict(cox_model, newdata = test, type = "risk")
tree_pred_cv <- predict(tree_min_cv, newdata = test)
tree_pred_1se <- predict(tree_1se, newdata = test) 
rf_pred <- predict(rf_model, newdata = test)$predicted


# Test set survival object
test_surv <- Surv(test$surv_time, test$censored)


# C-index calculation
cox_c <- concordance(test_surv ~ cox_pred)$concordance
tree_cv_c <- concordance(test_surv ~ tree_pred_cv)$concordance  
tree_1se_c <- concordance(test_surv ~ tree_pred_1se)$concordance
rf_c <- concordance(test_surv ~ rf_pred)$concordance

cat("C-index Results:\n")
cat("Cox Model:", round(cox_c, 4), "\n")
cat("Tree (Min-CV):", round(tree_cv_c, 4), "\n") 
cat("Tree (1-SE):", round(tree_1se_c, 4), "\n")
cat("Random Forest:", round(rf_c, 4), "\n\n")



library(timeROC)

# Time-dependent AUC at median survival time
median_time <- quantile(test$surv_time[test$censored == 1], 0.5, na.rm = TRUE)

# AUC calculations
cox_auc <- timeROC(T = test$surv_time, delta = test$censored, 
                   marker = cox_pred, times = median_time, cause = 1)$AUC[2]
tree_cv_auc <- timeROC(T = test$surv_time, delta = test$censored,
                       marker = tree_pred_cv, times = median_time, cause = 1)$AUC[2]
tree_1se_auc <- timeROC(T = test$surv_time, delta = test$censored,
                        marker = tree_pred_1se, times = median_time, cause = 1)$AUC[2]
rf_auc <- timeROC(T = test$surv_time, delta = test$censored,
                  marker = rf_pred, times = median_time, cause = 1)$AUC[2]

cat("AUC Results (at", round(median_time), "days):\n")
cat("Cox Model:", round(cox_auc, 4), "\n")
cat("Tree (Min-CV):", round(tree_cv_auc, 4), "\n")
cat("Tree (1-SE):", round(tree_1se_auc, 4), "\n") 
cat("Random Forest:", round(rf_auc, 4), "\n\n")




# C-Index comparison for all four models
c_values <- c(cox_c, tree_cv_c, tree_1se_c, rf_c)
barplot(c_values, 
        names.arg = c("Cox", "Tree-CV", "Tree-1SE", "RF"),
        main = "C-Index Comparison", ylab = "C-Index",
        col = c("lightblue", "lightgreen", "orange", "lightcoral"),
        ylim = c(0, max(c_values, na.rm = TRUE) * 1.1))
text(1:length(c_values),
     c_values + 0.02,
     labels = format(c_values, digits = 3),
     cex = 0.85)

# AUC comparison for all four models
auc_values <- c(cox_auc, tree_cv_auc, tree_1se_auc, rf_auc)
barplot(auc_values,
        names.arg = c("Cox", "Tree-CV", "Tree-1SE", "RF"), 
        main = "AUC Comparison", ylab = "AUC",
        col = c("lightblue", "lightgreen", "orange", "lightcoral"),
        ylim = c(0, max(auc_values, na.rm = TRUE) * 1.1))
text(1:length(auc_values),
     auc_values + 0.02,
     labels = format(auc_values, digits = 3),
     cex = 0.85)



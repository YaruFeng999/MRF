library(survival)
library(survminer)
library(eoffice)
library(ggsci)
library(Hmisc)
library(ggplot2)
library(ggcorrplot)
library(tidyr, help, pos = 2, lib.loc = NULL)
library(data.table)
library(rstatix)
library(ggcorrplot, help, pos = 2, lib.loc = NULL)
library(dplyr)
library(openxlsx)
library(caret)
library(Metrics)
library(timeROC)
library(openxlsx)

# KM
mrf$EMVI_invasion=factor(mrf$EMVI_invasion,levels=c("1","0"))
fit_tumor <- survfit(Surv(survival, status) ~ EMVI_invasion, data = mrf)

ggsurvplot(fit_tumor, 
           pval=TRUE, 
           pval.size =5,
           risk.table=TRUE,
           palette =c("#DC0000FF","#3C5488FF"),
) 




# multi-cox
mySurv <- Surv(mrf$dfs, mrf$dfs_status)
n <- coxph(
  mySurv ~ Tumor_invasion + LNM_invasion + TDs_invasion + EMVI_invasion + T +
    N + distance + length+gender+age+CEA+CA199+histology,
  data = mrf
)



results_table <- tidy(n, exponentiate = TRUE, conf.int = TRUE) %>% 
  dplyr::select(    # 明确指定使用dplyr的select
    term, estimate, conf.low, conf.high, p.value
  ) %>% 
  mutate(
    CI = sprintf("%.2f - %.2f", conf.low, conf.high),
    p.value = ifelse(p.value < 0.001, "<0.001", format.pval(p.value, digits = 2)),
    estimate = round(estimate, 2)
  ) %>% 
  rename(
    Variable = term,
    HR = estimate,
    `95% CI` = CI,
    `P-value` = p.value
  )



# Lasso

set.seed(80694571)
sam <- createDataPartition(mrf$survival, p = .7, list = FALSE)
train <- mrf[sam, ]
valid <- mrf[-sam, ]
valid_answer <- mrf[-sam, "three_dfs"]
train_answer <- mrf[sam, "three_dfs"]
x <- as.matrix(train[, c(
  "Tumor_invasion",
  "LNM_invasion",
  "TDs_invasion",
  "EMVI_invasion",
  "gender",
  "age",
  "CEA","CA199",
  "T",
  "N",
  "distance","length"
  ,"histology"
  ,"grade"
)])
y <- unlist(train$three_dfs)


set.seed(80694571)
model_lasso <- glmnet(x, y, alpha = 1)
cv_model <- cv.glmnet(x, y, alpha = 1)
model_lasso_min <- glmnet(x, y, alpha = 1, lambda = cv_model$lambda.min)
choose_signature_min <- row.names(model_lasso_min$beta)[as.numeric(model_lasso_min$beta) != 0]
index <- model_lasso_min$beta[as.numeric(model_lasso_min$beta) != 0]
score_list <- data.frame(choose_signature_min, index)

valid <- valid[, choose_signature_min]
valid <- as.data.frame(valid)
riskscore <- c()

for (i in 1:nrow(valid)) {
  score <- sum(as.numeric(valid[i, ]) * as.numeric(score_list$index))
  riskscore[i] <- score
}
valid$riskscore <- riskscore


riskscore <- c()
for (i in 1:nrow(train)) {
  score <- sum(as.numeric(train[i, ]) * as.numeric(score_list$index))
  riskscore[i] <- score
}
train$riskscore <- riskscore


gfit <- roc(train$three_dfs, train$riskscore)
ci <- ci(gfit)

pred <- prediction(train$riskscore, train$three_dfs)
auc <- performance(pred, measure = "auc")@y.values[[1]]
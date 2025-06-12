


remove(list = ls())
#### moedel
######## RF ##########


library(dplyr)
library(data.table)
library(randomForest) 
library(caret) 
library(pROC)
library(ggplot2) 
library(ggpubr)
library(ggprism)
library(caTools)

data1 <- readRDS("./result/ML/ML_usedata.rds")
data <- data1[,-1]
colnames(data)[1] <- "type"
set.seed(123)  
split <- sample.split(data$type, SplitRatio = 0.7)  
train_data <- subset(data, split == TRUE)  
test_data <- subset(data, split == FALSE)  

X_train <- train_data[, -1]
y_train <- as.factor(train_data[, 1])

model <- randomForest(x = X_train, y = y_train, ntree = 100)
print(model)

set.seed(1234)
ctrl <- trainControl(method = "cv", number = 10) 
grid <- expand.grid(mtry = c(2, 4, 6, 8 ,10))

rf_model <- train(x = X_train, y = y_train,
                  method = "rf",
                  trControl = ctrl,
                  tuneGrid = grid)


print(rf_model)

grid <- expand.grid(mtry = c(2))  
modellist <- list()

for (ntree in c(100, 200, 300, 400, 500)) {
  set.seed(12346)
  fit <- train(x = X_train, y = y_train, method="rf", 
               metric="Accuracy", tuneGrid=grid, 
               trControl=ctrl, ntree=ntree)
  key <- toString(ntree)
  modellist[[key]] <- fit
}

results <- resamples(modellist)
summary(results)


final_model <- randomForest(x = X_train, y = y_train,mtry = 2,ntree = 100)
# final_model <- readRDS("./result/ML/RF_model.rds")

print(final_model)
saveRDS(final_model,"./result/ML/RF_model.rds")

X_test <- test_data[, -1]
y_test <- as.factor(test_data[, 1])
test_predictions <- predict(final_model, newdata = test_data)

confusion_matrix <- confusionMatrix(test_predictions, y_test)
accuracy <- confusion_matrix$overall["Accuracy"]
precision <- confusion_matrix$byClass["Pos Pred Value"]
recall <- confusion_matrix$byClass["Sensitivity"]
f1_score <- confusion_matrix$byClass["F1"]

print(confusion_matrix)
print(paste("Accuracy:", accuracy))
print(paste("Precision:", precision))
print(paste("Recall:", recall))
print(paste("F1 Score:", f1_score))


confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("NonPE","PE")
rownames(confusion_matrix_df) <- c("NonPE","PE")
draw_data <- round(confusion_matrix_df / rowSums(confusion_matrix_df),2)
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)

ggplot(draw_data, aes(real,variable, fill = value)) +
  geom_tile() +
  geom_step(size = 1.2) + 
  geom_text(aes(label = scales::percent(value))) +
  scale_fill_gradient(low = "#F0F0F0", high = "#3575b5") +
  labs(x = "True", y = "Guess", title = "Confusion matrix") +
  theme_prism(border = T)+
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position="none")
test_predictions <- predict(final_model, newdata = test_data,type = "prob")

roc_obj <- roc(response = y_test, predictor = test_predictions[, 2])
roc_auc <- auc(roc_obj)

plot(roc_obj,
     add = FALSE,   
     col = 'red',   
     legacy.axes = TRUE,  
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "RF AUC: %.2f",  
     print.auc.x = 0.5,   
     print.auc.y = 0.5,   
     print.thres = F)

saveRDS(roc_obj,"./result/ML/ROC_result/RF_ROC.rds")



######## XGBoost ##########
## 
library(xgboost) 
library(Matrix) 

data1 <- readRDS("./result/ML/ML_usedata.rds")
data <- data1[,-1]
colnames(data)[1] <- "type"

set.seed(123)  
split <- sample.split(data$type, SplitRatio = 0.7)  
train_data <- subset(data, split == TRUE)  
test_data <- subset(data, split == FALSE)  

X_train <- train_data[, -1]
y_train <- train_data[, 1]

dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
params <- list(objective = "binary:logistic", eval_metric = "logloss", eta = 0.1, max_depth = 3)
nrounds <- 100
xgb_model <- xgboost(params = params, data = dtrain, nrounds = nrounds)

train_predictions <- predict(xgb_model, newdata = dtrain)
train_predictions <- ifelse(train_predictions > 0.5,1,0)

accuracy <- mean(train_predictions == y_train)
print(paste("accuracy:", accuracy))

X_test <- test_data[, -1]
y_test <- as.factor(test_data[, 1])

dtest <- xgb.DMatrix(data = as.matrix(X_test))
test_predictions <- predict(xgb_model, newdata = dtest)
test_predictions <- ifelse(test_predictions > 0.5,1,0)

ctrl <- trainControl(
  method = "cv",   
  number = 10,    
  verboseIter = FALSE)

param_grid <- expand.grid(
  nrounds = c(100, 200), 
  max_depth = c(3, 6), 
  eta = c(0.1), 
  gamma = c(0, 0.1), 
  colsample_bytree = c(0.8), 
  min_child_weight = c(1, 3), 
  subsample = c(0.8)) 

xgb_model <- train(
  x = X_train,
  y = y_train,
  method = "xgbTree",
  trControl = ctrl,
  tuneGrid = param_grid)

print(xgb_model$bestTune)

params <- list(objective = "binary:logistic", eval_metric = "logloss",
               nrounds = 100, eta = 0.1, max_depth = 6, gamma = 0.1,
               colsample_bytree = 0.8,
               min_child_weight = 3,
               subsample = 0.8)

xgb_model_final <- xgb.train(params = params, data = dtrain, nrounds = 100)
saveRDS(xgb_model_final,"./result/ML/XGB_model.rds")

train_predictions <- predict(xgb_model_final, newdata = dtrain)
train_predictions <- ifelse(train_predictions > 0.5,1,0)
accuracy <- mean(train_predictions == y_train)
print(paste("accuracy:", accuracy))

X_test <- test_data[, -1]
y_test <- as.factor(test_data[, 1])
dtest <- xgb.DMatrix(data = as.matrix(X_test))
test_predictions <- predict(xgb_model_final, newdata = dtest)
test_predictions <- ifelse(test_predictions > 0.5,1,0)

accuracy <- mean(test_predictions == y_test)
print(paste("accuracy:", accuracy))

X_test <- test_data[, -1]
y_test <- as.factor(test_data[, 1])
dtest <- xgb.DMatrix(data = as.matrix(X_test))
test_predictions <- predict(xgb_model_final, newdata = dtest)


roc_obj <- roc(response = y_test, predictor = test_predictions)
roc_auc <- auc(roc_obj)

plot(roc_obj,
     add = FALSE,  
     col = 'red',   
     legacy.axes = TRUE,   
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "RF AUC: %.2f",  
     print.auc.x = 0.5,   
     print.auc.y = 0.5,   
     print.thres = F)

saveRDS(roc_obj,"./result/ML/ROC_result/XGBoost_ROC.rds")



######## SVM ##########
library(e1071)

data1 <- readRDS("./result/ML/ML_usedata.rds")
data <- data1[,-1]
colnames(data)[1] <- "type"

set.seed(123)  
split <- sample.split(data$type, SplitRatio = 0.7)  
train_data <- subset(data, split == TRUE) 
test_data <- subset(data, split == FALSE)  

X_train <- train_data[, -1]
y_train <- as.factor(train_data[, 1])

svm_model <- svm(x = X_train, y = y_train)

train_predictions <- predict(svm_model, newdata = X_train)

table(y_train, train_predictions)

param_grid <- expand.grid(
  sigma = c(0.1, 1, 10),
  C = c(0.1, 1, 10))

ctrl <- trainControl(method = "cv", number = 10, verboseIter = FALSE)

tuned_model <- train(
  x = X_train,
  y = y_train,
  method = "svmRadial",
  tuneGrid = param_grid,
  trControl = ctrl)
print(tuned_model)


svm_model <- svm(x = X_train, y = y_train, probability = TRUE,sigma = 0.1,coef0 = 0.1)
saveRDS(svm_model,"./result/ML/SVM_model.rds")
test_predictions <- predict(svm_model, newdata = X_test,probability = T)

prob_estimates <- attr(test_predictions, "probabilities")

print(prob_estimates)

roc_obj <- roc(response = y_test, predictor = prob_estimates[,1])
roc_auc <- auc(roc_obj)

plot(roc_obj,
     add = FALSE,   
     col = 'red',  
     legacy.axes = TRUE,  
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "RF AUC: %.2f", 
     print.auc.x = 0.5,   
     print.auc.y = 0.5,   
     print.thres = F)

saveRDS(roc_obj,"./result/ML/ROC_result/SVM_ROC.rds")


######## LR ##########

data1 <- readRDS("./result/ML/ML_usedata.rds")
data <- data1[,-1]
colnames(data)[1] <- "type"

library(caTools) 

set.seed(123)  
split <- sample.split(data$type, SplitRatio = 0.7)  
train_data <- subset(data, split == TRUE)  
test_data <- subset(data, split == FALSE) 

model <- glm(type ~ ., data = train_data, family = binomial)
saveRDS(model,"./result/ML/LR_model.rds")
summary(model)

predictions <- predict(model, newdata = test_data, cluster = "response")

threshold <- 0.5  
predicted_classes <- ifelse(predictions >= threshold, 1, 0)

roc_obj <- roc(test_data$type, predictions)

plot(roc_obj,
     add = FALSE,   
     col = 'red',   
     legacy.axes = TRUE,   
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "LR AUC: %.3f",  
     print.auc.x = 0.5,   
     print.auc.y = 0.5,   
     print.thres = F)

saveRDS(roc_obj,"./result/ML/ROC_result/LR_ROC.rds")


remove(list = ls())

XGBoost_ROC <- readRDS("./result/ML/ROC_result/XGBoost_ROC.rds")
LR_ROC <- readRDS("./result/ML/ROC_result/LR_ROC.rds")
SVM_ROC <- readRDS("./result/ML/ROC_result/SVM_ROC.rds")
RF_ROC <- readRDS("./result/ML/ROC_result/RF_ROC.rds")


plot(XGBoost_ROC,
     add = FALSE,   
     col = 'red',   
     legacy.axes = TRUE,   
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "XGBoost: AUC = %.2f (GSE75010_test)", 
     print.auc.x = 0.65,   
     print.auc.y = 0.2,  
     print.thres = F)

plot(LR_ROC,
     add = FALSE,   
     col = 'red',  
     legacy.axes = TRUE,   
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "LR: AUC = %.2f (GSE75010_test)",  
     print.auc.x = 0.65,   
     print.auc.y = 0.2,  
     print.thres = F)

plot(RF_ROC,
     add = FALSE,   
     col = 'red',   
     legacy.axes = TRUE,   
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "RF: AUC = %.2f (GSE75010_test)",  
     print.auc.x = 0.65,   
     print.auc.y = 0.2,   
     print.thres = F)


plot(SVM_ROC,
     add = FALSE,   
     col = 'red',   
     legacy.axes = TRUE,   
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "SVM: AUC = %.2f (GSE75010_test)",  
     print.auc.x = 0.65,   
     print.auc.y = 0.2,   
     print.thres = F)



library(ggplot2) 
library(RColorBrewer) 
library(tidyverse) 
library(randomForest) 
library(rfPermute) 
library(ggthemes)
library(viridis)
rf_data <- readRDS("./result/data/")
rf_data <- rf_data[, -1]
rf_data$group <- as.factor(rf_data$group)
set.seed(1234)
df.rf <- rfPermute(group ~ ., 
                   data= rf_data, 
                   ntree = 300,
                   nrep = 299, 
                   num.cores = 4)
importance_scores <- importance(df.rf)
feature_importance <- data.frame(Gene = rownames(importance_scores), Importance = importance_scores[, "MeanDecreaseGini"])
ordered_features <- feature_importance[order(-feature_importance$Importance), ]
importance_gini <- importance(df.rf, type = 2)
importance_accuracy <- importance(df.rf, type = 1)
top_genes <- head(ordered_features, 50)
write.csv(top_genes, file = "./use_data/rf_top50_genes.csv", row.names = F)

top_genes <-read.table("./use_data/rf_top50_genes.csv",sep=",",header=T)
color_range <- c("#D8BFD8", "#8B008B")
ggplot(top_genes, aes(x = reorder(Gene, Importance), y = Importance, fill = Importance, label = round(Importance, 2))) +
  geom_bar(stat = "identity") +
  geom_text(size = 2.5, position = position_stack(vjust = 0.5), color = "white") + 
  scale_fill_gradient(low = color_range[1], high = color_range[2]) + 
  theme_few() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 16), 
    plot.caption = element_text(size = 12), 
    axis.text = element_text(size = 8), 
    axis.title = element_text(size = 15)
  ) +
  labs(title = "Top 50 Important Genes in Heart Disease Prediction", x = "Gene", y = "Importance") +
  coord_flip() 


## Validation
remove(list = ls())
exp <- readRDS("./use_data/GSE25906_exp.rds")
use_exp <- as.data.frame(t(exp[colnames(model_data)[3:13],]))
use_exp$sample <- rownames(use_exp)
cli <- readRDS("./use_data/GSE25906_clin.rds")
use_clin <- cli[,1:2]
use_data <- dplyr::inner_join(use_clin,use_exp)
use_data <- use_data %>% dplyr::mutate(status = ifelse(status == "NonPE",0,1))
saveRDS(use_data,"./result/Validation/GSE25906_model_data.rds")

use_data <- readRDS("./result/Validation/GSE25906_model_data.rds")
use_data <- use_data[,-1]

XGB <- readRDS("./result/ML/XGB_model.rds")

X_test <- use_data[, -1]
y_test <- as.factor(use_data[, 1])
dtest <- xgb.DMatrix(data = as.matrix(X_test))
test_predictions <- predict(XGB, newdata = dtest)

roc_obj <- roc(response = y_test, predictor = test_predictions)
roc_auc <- auc(roc_obj)

plot(roc_obj,
     add = FALSE,   
     col = 'red',   
     legacy.axes = TRUE,   
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "XGBoost AUC: %.2f",  
     print.auc.x = 0.5,   
     print.auc.y = 0.5,   
     print.thres = F)

saveRDS(roc_obj,"./result/Validation/GSE25906_XGB_ROC.rds")


LR <- readRDS("./result/ML/LR_model.rds")
predictions <- predict(LR, newdata = use_data, cluster = "response")

threshold <- 0.5 
predicted_classes <- ifelse(predictions >= threshold, 1, 0)

roc_obj <- roc(use_data$status, predictions)

plot(roc_obj,
     add = FALSE,   
     col = 'red',  
     legacy.axes = TRUE,  
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "LR AUC: %.3f",  
     print.auc.x = 0.5,   
     print.auc.y = 0.5,   
     print.thres = F)

saveRDS(roc_obj,"./result/Validation/GSE25906_LR_ROC.rds")

SVM <- readRDS("./result/ML/SVM_model.rds")
test_predictions <- predict(SVM, newdata = X_test,probability = T)
prob_estimates <- attr(test_predictions, "probabilities")
print(prob_estimates)
roc_obj <- roc(response = y_test, predictor = prob_estimates[,1])
roc_auc <- auc(roc_obj)

plot(roc_obj,
     add = FALSE,  
     col = 'red',   
     legacy.axes = TRUE,  
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "RF AUC: %.2f",  
     print.auc.x = 0.5,  
     print.auc.y = 0.5,  
     print.thres = F)

saveRDS(roc_obj,"./result/Validation/GSE25906_RF_ROC.rds")

RF <- readRDS("./result/ML/RF_model.rds")
test_predictions <- predict(RF, newdata = use_data,type = "prob")

roc_obj <- roc(response = y_test, predictor = test_predictions[, 2])
roc_auc <- auc(roc_obj)

plot(roc_obj,
     add = FALSE,   
     col = 'red',   
     legacy.axes = TRUE,   
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "RF AUC: %.2f", 
     print.auc.x = 0.5,   
     print.auc.y = 0.5,  
     print.thres = F)

saveRDS(roc_obj,"./result/Validation/GSE25906_SVM_ROC.rds")


####### Draw ########

remove(list = ls())

XGBoost_ROC <- readRDS("./result/ML/ROC_result/XGBoost_ROC.rds")
LR_ROC <- readRDS("./result/ML/ROC_result/LR_ROC.rds")
SVM_ROC <- readRDS("./result/ML/ROC_result/SVM_ROC.rds")
RF_ROC <- readRDS("./result/ML/ROC_result/RF_ROC.rds")


XGBoost_ROC1 <- readRDS("./result/Validation/GSE25906_XGB_ROC.rds")
LR_ROC1 <- readRDS("./result/Validation/GSE25906_LR_ROC.rds")
SVM_ROC1 <- readRDS("./result/Validation/GSE25906_SVM_ROC.rds")
RF_ROC1 <- readRDS("./result/Validation/GSE25906_RF_ROC.rds")


plot(XGBoost_ROC,
     add = FALSE,  
     col = 'red',  
     legacy.axes = TRUE,  
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "XGBoost: AUC = %.2f (GSE75010_test)", 
     print.auc.x = 0.7,   
     print.auc.y = 0.3,   
     print.thres = F,
     type = "S")


plot(XGBoost_ROC1,
     add = T,   
     col = 'blue',  
     legacy.axes = TRUE,   
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "XGBoost: AUC = %.2f (GSE25906_Validation)",  
     print.auc.x = 0.8,  
     print.auc.y = 0.2,   
     print.thres = F,
     type = "S")



plot(LR_ROC,
     add = FALSE,   
     col = 'red',   
     legacy.axes = TRUE,  
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "LR: AUC = %.2f (GSE75010_test)", 
     print.auc.x = 0.7,  
     print.auc.y = 0.3, 
     print.thres = F,
     type = "S")


plot(LR_ROC1,
     add = T, 
     col = 'blue', 
     legacy.axes = TRUE, 
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "LR: AUC = %.2f (GSE25906_Validation)",
     print.auc.x = 0.8,  
     print.auc.y = 0.2,  
     print.thres = F,
     type = "S")


plot(SVM_ROC,
     add = FALSE,
     col = 'red', 
     legacy.axes = TRUE,  
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "SVM: AUC = %.2f (GSE75010_test)", 
     print.auc.x = 0.7,  
     print.auc.y = 0.3, 
     print.thres = F,
     type = "S")


plot(SVM_ROC1,
     add = T,  
     col = 'blue',  
     legacy.axes = TRUE,  
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "SVM: AUC = %.2f (GSE25906_Validation)", 
     print.auc.x = 0.8,  
     print.auc.y = 0.2,  
     print.thres = F,
     type = "S")


plot(RF_ROC,
     add = FALSE,  
     col = 'red', 
     legacy.axes = TRUE,  
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "XGBoost: AUC = %.2f (GSE75010_test)", 
     print.auc.x = 0.7,  
     print.auc.y = 0.3,  
     print.thres = F,
     type = "S")


plot(RF_ROC1,
     add = T,  
     col = 'blue',  
     legacy.axes = TRUE,  
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.pattern = "RF: AUC = %.2f (GSE25906_Validation)",  
     print.auc.x = 0.8,  
     print.auc.y = 0.2,  
     print.thres = F,
     type = "S")



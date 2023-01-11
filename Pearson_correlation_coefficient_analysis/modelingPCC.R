#PCC response variable prediction
getwd()
setwd("/Users/xung0001/Desktop")
non <- read.table("cor_noncoding06292021.txt", sep="\t", header=T)
x <- non[,-c(13:14,6,2)]
# because I have same names for colum x$id1, I have to rename it
x$id1[59] <- "AT1G26220x"
row.names(x) <- x$id1
x <- x[,-1]
colnames(x)[3] <- "PCC"
# libraries loading
library(magrittr)
library(ggplot2)
library(cluster)
# heatmap function 
par(oma=c(0.8,0,0,0))
x %>% as.matrix() %>% scale() %>% heatmap(Rowv = T,Colv = T,cexCol = 0.8,cexRow = 0.5)#A4 
# variables importance using lasso 
# #find best lamda through regularisation
library(glmnet)
library(dplyr)
library(caret)
# standarize data x
x1 <- scale(x)
x1 <- as.data.frame(x1)
y <- x1%>%select(PCC)%>%as.matrix()
X <- x1%>%select(-PCC)%>%as.matrix()
# chosing lamda based on ridge regression
m_cv <- cv.glmnet(X, y, alpha = 0)
plot(m_cv)
m <- glmnet(X, y, alpha = 0,lambda = m_cv$lambda.min)
coef(m)
m_cv#lamda.min 0.3248448
# chosing lamda based on lasso regression
m_cv <- cv.glmnet(X, y, alpha = 1)
plot(m_cv)
m_cv$lambda.min #0.26
# final model
m <- glmnet(X, y, alpha = 0,lambda = m_cv$lambda.min)
coef(m)
#  lasso modeling
tc <- trainControl(method = "repeatedcv",
                   number = 10, repeats = 100)
m3 <- train(PCC_value~ .,
            data = x1,
            method = "glmnet", 
            tuneGrid = expand.grid(alpha = 1,
                                   lambda = seq(0, 10, 0.1)),
            metric = "RMSE",
            trControl = tc)
m3$bestTune
m3$results[which(rownames(m3$results) == rownames(m3$bestTune)),]
coef(m3$finalModel, m3$bestTune$lambda)
ggplot(m3)
rownames(coef(m3$finalModel, m3$finalModel$lambdaOpt))[
  coef(m3$finalModel, m3$finalModel$lambdaOpt)[,1]!= 0]

# https://www.guttmacher.org/public-use-datasets
# https://www.guttmacher.org/adding-it-up

# Support2
#“Data obtained from http://hbiostat.org/data courtesy of the Vanderbilt University Department of Biostatistics.”


#----------------Setup--------------------
if(!require(tidyverse)) install.packages("tidyverse", repos = "https://cran.rstudio.com/")
if(!require(Hmisc)) install.packages("Hmisc", repos = "https://cran.rstudio.com/")
if(!require(kamila)) install.packages("kamila", repos = "https://cran.rstudio.com/")
if(!require(caret)) install.packages("caret", repos = "https://cran.rstudio.com/")
if(!require(recipes)) install.packages("recipes", repos = "https://cran.rstudio.com/")
if(!require(MLmetrics)) install.packages("MLmetrics", repos = "https://cran.rstudio.com/")
if(!require(GGally)) install.packages("GGally", repos = "https://cran.rstudio.com/")
if(!require(ggbiplot)) install.packages("ggbiplot", repos = "https://cran.rstudio.com/")
if(!require(pls)) install.packages("pls", repos = "https://cran.rstudio.com/")
if(!require(gbm)) install.packages("gbm", repos = "https://cran.rstudio.com/")



#library(readr)
library(tidyverse)
library(Hmisc) # data set is available through this package
library(kamila)
library(caret)
library(recipes)
library(MLmetrics) #LogLoss/Cross-Entropy Loss
library(GGally)
library(ggbiplot)
library(pls)
library(gbm)

#color blind friendly colors
cb_cols <- 
  c("#999999", "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#options(timeout = 120, digits = 2)

#EXPLORE BINARY CROSS ENTROPY instead of RMSE
# https://www.geeksforgeeks.org/binary-cross-entropy-in-r/
# https://medium.com/@jwbtmf/loss-metrics-for-deep-learning-binary-cross-entropy-vs-rmse-8cd5fa80a1e5#:~:text=One%20key%20difference%20between%20binary%20cross-entropy%20and%20RMSE,sensitive%20to%20large%20errors%20in%20the%20predicted%20values.
# https://pareto.ai/blog/cross-entropy-loss-function
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5358992/
# Function to calculate BCE


getHdata(support2)
#remove fields that are findings from previous models and those deceased before discharge 
support <- data.frame(support2) %>% 
  select(-aps, -sps, -surv2m, -surv6m, -prg2m, -prg6m, -dnr, -dnrday) %>%
  filter(hospdead == 0)

# Final hold-out test set will be 10% of patient's data
set.seed(365, sample.kind="Rounding") # if using R 3.6 or later
# set.seed(365) # if using R 3.5 or earlier
test_index <- createDataPartition(y = support$death, times = 1, p = 0.1, list = FALSE)
main <- support[-test_index,]
final <- support[test_index,]

#-----------Exploration-------------------

#function from here - https://stackoverflow.com/questions/2394902/remove-variable-labels-attached-with-foreign-hmisc-spss-import-functions
#Clears the labels from the data frame, since it was causing issues with ggplot
clear.labels <- function(x) {
  if(is.list(x)) {
    for(i in seq_along(x)) {
      class(x[[i]]) <- setdiff(class(x[[i]]), 'labelled') 
      attr(x[[i]],"label") <- NULL
    } 
  } else {
    class(x) <- setdiff(class(x), "labelled")
    attr(x, "label") <- NULL
  }
  return(x)
}

main <- clear.labels(main)

summary(main)

# checking which columns have NA and how many
x <- colSums(is.na(main))
# less than 25% missing
x[which(x < (nrow(main)/4) & x > 0)]
# more than 25% missing values
x[which(x >= (nrow(main)/4))]
# no missing values
x[which(x == 0)]

#effect of age, higher spread 
main %>%  ggplot(aes(age, fill = (death == 1))) + geom_histogram(binwidth = 10)


#effect of sex
main %>% group_by(sex) %>% summarise(n())
main %>% group_by(sex, death) %>% ggplot(aes(death == 1,age)) + geom_boxplot(aes(color = sex))

#exploring race as a factor
main %>% group_by(race, death) %>% summarise(age = mean(age)) %>% arrange(desc(death))
main %>% group_by(race, death) %>% summarise(age = mean(age)) %>% 
  ggplot(aes(race,age, fill = death == 1)) + 
  geom_bar(position="dodge", stat="identity") 
main %>% group_by(race, death)  %>% 
  ggplot(aes(race, fill = death == 1)) + 
  geom_histogram(stat="count", position="dodge") 

# Is there is a large difference in survival rates when insurance covers more of the cost
main %>% 
  ggplot(aes(death == 1,totcst)) + 
  geom_boxplot() + scale_y_log10()

#dzclass and dzgroup(subgroup of dzclass) are highly correlated, keep only one

# The number of simultaneous diseases (or comorbidities)
main %>% ggplot(aes(num.co,  fill = death == 1)) + 
  geom_histogram(stat="count", position="dodge") 


# Possible causes columns based on too many missing values or is a dependent var
modeling_index = sort(c(c(1,5,9,10,12,18,seq(22,31),34,35), c(3,7,17,19,20,21)))
# adding column 2 as predictor, adjust the indicies to view different comparisons
# visual correlation matrix
# Can take several minutes to run, my system would fail on over 12 at a time. 
# 1-10 seemed to have the best results
reference <- main[,c(2,modeling_index[1:10])]
ggpairs(reference, aes(color = as.factor(death), alpha = 0.4))


#------------ Cleaning up NAs



# Baseline Normal Fill-in Values (recommened by source)
# Serum albumin	3.5
# PaO2/FiO2 ratio (pafi)	333.3
# Bilirubin	1.01
# Creatinine	1.01
# BUN	6.51
# White blood count	9 (thousands)
# Urine output	2502
#
# Unknown race set at other
# education assume finished high school = 12 (may create a bias)
main <- main %>% replace_na(list(alb = 3.5, pafi = 333.3,
                                 bili = 1.01, crea = 1.01,
                                 bun = 6.51, wblc = 9,
                                 urine = 2502,
                                 race = "other",
                                 edu = 12))





#use mean for NA, numeric
main <- main %>% mutate_at(vars("charges", "totcst", "avtisst"),
                           ~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))
#mode for catagorical
main <- main %>% mutate_at(vars("sfdm2"),
                           ~ifelse(is.na(.x), mode(.x), .x))


# check correlation in numeric data
sort(abs(cor(main %>% select(where(is.numeric)), use = "complete.obs"))[,2], decreasing=TRUE)



#----------- METHODS ------------
# create a testing and training set
set.seed(713)
test_index <- createDataPartition(y = main$death, times = 1, p = 0.2, list = FALSE)
train <- main[-test_index,]
test <- main[test_index,]




# just using mean
mu <- train %>% summarise(x = mean(death)) %>% pull(x)
RMSE(mu,test$death)
# Calculate BCE
LogLoss(mu, test$death)

#Linear Regression
relevant_data <- train[,c(2,modeling_index)]
fit <- lm(death ~ ., data = relevant_data)
RMSE(predict(fit,test),test$death)
LogLoss(predict(fit,test), test$death)

# graident boosted method
# significantly worse


lambdas <- seq(2, 10)
trees <- sapply(lambdas, function(lambda){
  fit <- gbm(death ~ ., data = relevant_data, cv.folds = 2, 
             n.trees = lambda,
             distribution = "bernoulli")
  LogLoss(predict(fit,test),test$death)
})
qplot(lambdas, trees, geom = "line")
lambda_u <- lambdas[which.min(trees)]
print(lambda_u)


fit <- gbm(death ~ ., data = relevant_data, cv.folds = 5, n.trees = lambda_u)
RMSE(predict(fit,test),test$death)
LogLoss(predict(fit,test), test$death)



#Exploring PCA----------------------------------------------------------------

#explaining PCA https://www.turing.com/kb/guide-to-principal-component-analysis

#calculate principal components
PCA_data <- relevant_data %>% mutate(death_info = paste(row_number(),death))
row.names(PCA_data) =   PCA_data$death_info
# remove categorical data
results <- prcomp(PCA_data[,-c(3,5,9,13,26)], scale = TRUE)
summary(results)

#scree plot
var_explained <- results$sdev^2/sum(results$sdev^2)
qplot(c(1:21), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 0.15)

# top results in relation to death
sort(abs(results$rotation[1,]), decreasing = TRUE)

plotting_results <- results$x %>% as.data.frame %>% rownames_to_column("death_info") %>%
  separate("death_info",c("row","death")) 

# overall plot with the components that explain overall variance
plotting_results %>%
  ggplot(aes(PC1, PC2)) + geom_point(aes(color = death, alpha = 0.2))
# top 2 components that relate to death instead
plotting_results %>%
  ggplot(aes(PC11, PC9)) + geom_point(aes(color = death, alpha = 0.2))

# ordering features based on importance for overall top and death top
y<-data.frame(results$rotation) %>% mutate(feature = row.names(.))
y %>% ggplot(aes(PC1, reorder(feature, PC1))) + geom_point()
y %>% ggplot(aes(PC11, reorder(feature, PC11))) + geom_point()


#biplot with top 2 
ggbiplot(results, choices = c(11,9), obs.scale = 1, var.scale = 1, 
         groups = as.factor(relevant_data$death)) + 
  scale_colour_manual(values = c("0" = cb_cols[1], "1" = cb_cols[2])) +
  theme(legend.direction = 'horizontal', legend.position = 'top') + 
  ylim(-6, 6)



y %>% ggplot(aes(PC11, PC9, label = feature)) +
  geom_text()


model <- lm(response ~ results$x[, 1] + results$x[, 2])
summary(model)
predict(results,test)

#make this example reproducible
set.seed(1)

#fit PCR model
model <- pcr(death~., data=relevant_data, scale=TRUE, validation="CV")
summary(model)

#visualize cross-validation plots
validationplot(model)

pcr_pred <- predict(model, test, ncomp=30)
RMSE(pcr_pred,test$death)
LogLoss(pcr_pred, test$death)
#---------------------------------



# Tree method
# Improvement on binary but loss on RMSE
# Min 5 observations per node
# 10 cross validations
library(rpart)
fit <- rpart(death ~ ., data = select(train, -c(d.time, sfdm2, hospdead)), 
             control = rpart.control( minsplit = 5, xval = 10))
plot(fit, margin = 0.1)
text(fit, cex = 0.75)
RMSE(predict(fit,test),test$death)
LogLoss(predict(fit,test), test$death)


#----RECIPES

# -1 due death being removed from the indexes when creating recipes
id_index = sapply(c(4,6,8,11,seq(13,16),32,33,seq(36,39)),function(x){x-1}) 

my_recipe <- recipe(death ~ ., data = train) %>%
  update_role(all_of(id_index), new_role = "ID")  #sets the title feature to be kept but not used for modeling
  
  #  step_dummy(genres, one_hot = TRUE)  #One Hot encoding the genres

my_recipe
print(summary(my_recipe),n = 40)
#library(kernlab)


#loess method
# No improvement
b <- 5
control <- trainControl(method = "cv", number = b, p = .9)
train_Loess <- train(my_recipe, train,
                     method = "gamLoess",
                     #tuneGrid = data.frame(k = seq(1,60,2)),
                     trControl = control)

RMSE(predict(train_Loess, test), test$death)
LogLoss(predict(train_Loess, test), test$death)

# #knn method
# b <- 5
# control <- trainControl(method = "cv", number = b, p = .9)
# train_knn <- train(my_recipe, train,
#                      method = "knn",
#                      tuneGrid = data.frame(k = seq(1,60,2)),
#                      trControl = control)
# 
# RMSE(predict(train_Loess, testing), testing$death)
# LogLoss(predict(train_Loess, test), test$death)



#-----------------


 
 
 
 #------------------
 # # Kamila method, ran for 15-30 min on my device and pushed memory limits but is possible
 # # supposed to be good for lots of categorical data (genres) with better performance. 
 # # But it did not beat regularized data above
 set.seed(52)
 kamind <- createDataPartition(main$death, p = 0.7, list = FALSE)
 kamtrain <- main[kamind,]
 kamtest <- main[-kamind,]
 
 
 #continous variables (numbers) 
 conInd <- c(1,5,9,10,12,18,seq(22,31),34,35)
 conVars <- kamtrain[,conInd]
 conTest <- kamtest[,conInd]
 
 #categorical variables
 # (dzgroup and dzclass highly correlated, just using dzgroup)
 catInd <- c(3,7,17,19,20,21)
 #catInd2 <- c(1,2,6) #genre groups in tact
 catVarsFac <- kamtrain[,catInd]
 catVarsFac[] <- lapply(catVarsFac, factor)

 #make sure cat variables are factors
 catTestFac <- kamtest[,catInd]
 catTestFac[] <- lapply(catTestFac, factor)
 
 
 
 
 
 # KAMILA
 kamRes <- kamila(conVars, catVarsFac, numClust=8, numInit=5)
 kamPred <- classifyKamila(kamRes, list(conTest, catTestFac))
 table(kamtest$death, kamPred)
 
 RMSE(kamPred,kamtest$death)
 LogLoss(kamtest$death,kamPred)
 
 # FIGURE OUT HOW TO DO THIS 
 # Plot KAMILA results
 plottingData <- cbind(
   conVars,
   catVarsFac,
   KamilaCluster = factor(kamRes$finalMemb),
   death = factor(kamtrain$death))
 plottingData$result <- ifelse(
   plottingData$KamilaCluster == '1', yes='Yes',no='No')
 
 plottingData %>% ggplot(aes(KamilaCluster, death)) + geom_histogram()
 
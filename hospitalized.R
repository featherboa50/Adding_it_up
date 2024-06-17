
# https://www.guttmacher.org/public-use-datasets
# https://www.guttmacher.org/adding-it-up

# Support2
#“Data obtained from http://hbiostat.org/data courtesy of the Vanderbilt University Department of Biostatistics.”

if(!require(tidyverse)) install.packages("tidyverse", repos = "https://cran.rstudio.com/")
if(!require(Hmisc)) install.packages("Hmisc", repos = "https://cran.rstudio.com/")
if(!require(kamila)) install.packages("kamila", repos = "https://cran.rstudio.com/")
if(!require(caret)) install.packages("caret", repos = "https://cran.rstudio.com/")


#library(readr)
library(tidyverse)
library(Hmisc) # data set is available through this package
library(kamila)
library(caret)

#color blind friendly colors
cb_cols <- 
  c("#999999", "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#options(timeout = 120, digits = 2)




getHdata(support2)
#remove fields that are findings from previous models 
support <- data.frame(support2) %>% select(-aps, -sps, -surv2m, -surv6m, -prg2m, -prg6m, -dnr, -dnrday)

# Final hold-out test set will be 10% of patient's data
set.seed(365, sample.kind="Rounding") # if using R 3.6 or later
# set.seed(365) # if using R 3.5 or earlier
test_index <- createDataPartition(y = support$death, times = 1, p = 0.1, list = FALSE)
main <- support[-test_index,]
final <- support[test_index,]



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

main %>% group_by(sex) %>% summarise(n())

main %>% group_by(sex, death) %>% ggplot(aes(death == 1,age)) + geom_boxplot(aes(color = sex))

#exploring race as
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


# checking which columns have NA and how many
x <- colSums(is.na(main))
# less than 25% missing
x[which(x < 2048 & x > 0)]

# more than 25% missing values
x[which(x >= 2048)]


# Baseline Normal Fill-in Values
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
main <- main %>% mutate_at(vars("meanbp", "hrt", "resp", "scoma", "sod", "temp"),
                                   ~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5358992/



#mode for NA, categorical data
support2 <- support2 %>% mutate_at(vars("ca"),
                                   ~ifelse(is.na(.x), mode(.x), .x))
# check correlation in numeric data
sort(abs(cor(support2 %>% select(where(is.numeric)), use = "complete.obs"))[,2], decreasing=TRUE)

support2 %>% filter(death == 1) %>% pairs(~sex + bun  + sod + race, data = .)
 charges totcst avtisst race  temp sfdm2 
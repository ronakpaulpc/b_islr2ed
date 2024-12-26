# Contains code of "Introduction to Statistical Learning with Applications in R".
# The chapter-wise code sections contains Lab codes from each chapter end.
# Let's get going!!!


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C0 - Libraries and packages ---------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Working directory check
getwd()
# "C:/Users/Ronak/Documents/ALL Research/Rsoftware/b_islr/islr2ed"

# Installing required pkgs
# install.packages("ISLR2")

# Loading the pkgs
library(easypackages)
libraries(
    
    # Book package
    "ISLR2",
    
    # Data io
    "rio",
    "here",
    
    # Data manipulate
    "janitor",
    "tidyverse",
    "skimr",
    
    # Analysis
    "survival"
)
# ctrl+shift+2 to see console in full view.



#_====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C1 - Introduction -------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NO CODE.



#_====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C2 - Statistical Learning -----------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Code for Lab on Introduction to R
# In this lab, we will introduce some simple R commands. 
# The best way to learn a new language is to try out the commands. 

# Creating a vector x
x <- c(1, 3, 2, 5)
x
# We can also save things using = rather than <-
x = c(1, 6, 2)
x
y = c(1, 4, 3)
y

# 2.3.1 Basic Commands ----------------------------------------------------
# object length and arithmetic operation
length(x)
length(y)
x + y

# list of all of the objects
ls()
# delete objects
rm(x, y)
ls()

# remove all objects at once
x <- c(1, 6, 2)
y <- c(1, 4, 3)
rm(list = ls())


# Using the matrix() fn to create matrices
?matrix
# we create a simple matrix
x <- matrix(
    data = c(1, 2, 3, 4),
    nrow = 2,
    ncol = 2
)
print(x)

# NOTE: We could re-run the matrix cmd with no args and get the same result.
x <- matrix(c(1, 2, 3, 4), 2, 2)
print(x)
# However, it can sometimes be useful to specify the arg names.
# Otherwise R will assume that the fn args are passed into the function in 
# - the same order that is given in the function’s help fle.

# By default R creates matrices by successively filling in columns.
# We can fill matrices by order of rows using byrow = T arg.
matrix(
    data = c(1, 2, 3, 4), 
    nrow = 2, 
    ncol = 2, 
    byrow = T
)

# Square root of matrix elements
sqrt(x)
# Square of matrix elements
x ^ 2


# Generating vector of random normal var
?rnorm
# NOTE: By default, rnorm() creates std normal rdvars with mean 0 and sd 1.
# However, the mean and sd can be altered using the mean and sd arguments.
x <- rnorm(50)
y <- x + rnorm(50, mean = 50, sd = 0.1)
cor(x, y)

# We can reproduce the exact same set of random numbers using set.seed() fn.
# Before set.seed()
rnorm(5)
rnorm(5)
# After set.seed()
set.seed(5)
rnorm(5)
set.seed(5)
rnorm(5)

# Compute the mean and var of rdvar
set.seed(3)
y <- rnorm(100)
mean(y)
var(y)
sqrt(var(y))
sd(y)


# 2.3.2 Graphics ----------------------------------------------------------
# The plot() fn is the primary way to plot data in R.
?plot()
# plot(x, y) produces a scatter plot of the numbers in x vs y.
x <- rnorm(100)
y <- rnorm(100)
plot(x, y)
# Scatter plot with labels
plot(
    x, y,
    xlab = "This is the x-axis",
    ylab = "This is the y-axis",
    main = "Plot of X vs Y"
)

# We will often want to save the output of an R plot. 
# The cmd that we will use depend on the file type that we would create.
# to create a pdf, we use the pdf() function.
pdf("figure_ch2-3.pdf")
plot(x, y, col = "green")
dev.off()
# to create a jpeg, we use the jpeg() function.
jpeg("figure_ch2-3.jpeg")
plot(x, y, col = "orange")
dev.off()
# The function dev.off() indicates to R that we are done creating the plot.
# Altly, we can copy the plot window and paste it to a file, like a Word doc.

# The seq() fn can be used to create a sequence of numbers.
x <- seq(1, 10)
print(x)
x <- 1:10
print(x)
x <- seq(-pi, pi, length = 50)
print(x)

# Now we will create contour plot using countour() fn. 
# It represents 3-D data and is like a topographical map.
y <- x
f <- outer(x, y, function(x, y) cos(y)/(1 + x^2) )
contour(x, y, f)
contour(x, y, f, nlevels = 45, add = T)
fa <- (f - t(f))/2
contour(x, y, fa, nlevels = 15)

# Creating a heat map using image(). It works the same way as contour().
image(x, y, fa)

# Creating a 3-dimensional plot using persp().
persp(x, y, fa)
# The args theta and phi control the angles at which the plot is viewed.
persp(x, y, fa, theta = 30)
persp(x, y, fa, theta = 30, phi = 20)
persp(x, y, fa, theta = 30, phi = 70)
persp(x, y, fa, theta = 30, phi = 40)


# 2.3.3 Indexing Data -----------------------------------------------------
# We often wish to examine part of a set of data.
A <- matrix(1:16, 4, 4)
print(A)
A[2, 3]

# Selecting multiple rows and cols at a time
A[c(1, 3), c(2, 4)]
A[1:3, 2:4]
A[1:2, ]
A[ , 1:2]
# R treats a single row or column of a matrix as a vector.
A[1, ]
# Using "-" tells R to keep all rows & cols except those in the index.
A[-c(1, 3), ]
A[-c(1, 3), -c(1, 3, 4)]
# The dim() fn gives the no of rows and cols of a matrix.
dim(A)


# 2.3.4 Loading Data ------------------------------------------------------
# For most analyses, the first step involves importing a data set into R.
# The read.table() fn is one of the primary ways to do this.
# We can use the write.table() fn to export data.
?read.table
# Reading a file
auto <- read.table(
    file = here("data", "Auto.data")
)
head(auto)
View(auto)
# NOTE: The data is loaded incorrectly, because R has assumed that varnames
# - are part of the data and so has included them in the first row.
auto <- read.table(
    file = here("data", "Auto.data"),
    header = T,
    stringsAsFactors = T,
    na.strings = "?"
)
View(auto)
glimpse(auto)

# An easy way to load data from Excel into R is to save it as a csv file
# - and then use the read.csv() fn for importing in R.
auto <- read.csv(
    file = here("data", "Auto.csv"),
    header = T,
    stringsAsFactors = T,
    na.strings = "?"
)
glimpse(auto)
View(auto)
dim(auto)       # 397 X 9
# Checking part of dataset
auto[1:4, ]

# Finding the total missing values in dataset
sum(is.na(auto))        # 5 missing values
# Missing values in dataset variables
colSums(is.na(auto))
# Missing values from data summary
summary(auto)           # All missval in "horsepower" var
skim_without_charts(auto)
# Since there only five of the rows with missing obs we omit them.
auto <- na.omit(auto)
dim(auto)       # 392 X 9
# Checking the var names
names(auto)


# 2.3.5 Additional Graphical and Numerical Summaries ----------------------
# In plot() simply typing varnames will give an error as dataset is not 
# specified.
plot(cylinders, mpg)
# We need to specify dataset followed by $ and varname.
plot(auto$cylinders, auto$mpg)
# Altly, we can use attach() fn and make vars in dataframe available by name.
attach(auto)
plot(cylinders, mpg)

# The as.factor() fn converts quantitative vars into qualitative vars.
cylinders <- as.factor(cylinders)
class(cylinders)

# If x-axis var is qualitative then plot() will automatically give box plots.
plot(cylinders, mpg)
# Customizing the plot
plot(cylinders, mpg, col = "red")
plot(cylinders, mpg, col = "red", varwidth = T)
plot(cylinders, mpg, col = "red", varwidth = T, horizontal = T)
plot(
    cylinders, mpg, 
    col = "red", 
    varwidth = T, 
    xlab = "cylinders", 
    ylab = "MPG"
)

# Plotting histograms
hist(mpg)
hist(mpg, col = "red")
hist(mpg, col = "red", breaks = 15)

# pairs() creates a scatterplot matrix, i.e. a scatterplot for all var pairs.
pairs(auto)
# We can also produce scatterplots for just a subset of the vars.
pairs( ~ mpg + displacement + horsepower + weight + acceleration, data = auto)

# identify() allows for identifying the value of a var for points on a plot.
plot(horsepower, mpg)
identify(horsepower, mpg, name)

# summary() fn produces a numerical summary of each var in a dataset.
summary(auto)
# NOTE: For qualitative vars like name, R will list the num of obs in each 
# category. We can also produce a summary of just a single var in dataframe.
summary(auto$mpg)



#_====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C3 - Linear Regression --------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Lab on Linear Regression
# In this lab, we will learn how to perform linear regression in R. 
# Let's roll!


# 3.6.1 Libraries ---------------------------------------------------------
library(car)
library(MASS)
library(ISLR2)


# 3.6.2 Simple Linear Regression ------------------------------------------
head(Boston)
glimpse(Boston)
?Boston

# Incorrect cmd.
# The dataset to which the vars come from is not given.
lm.fit <- lm(medv ~ lstat)
# Gives ERROR

# Correct cmd.
lm.fit <- lm(medv ~ lstat, data = Boston)
# We need not specify data if it is attached.
# attach(Boston)
# lm.fit <- lm(medv ~ lstat)

# See basic information about the model.
print(lm.fit)
# See detailed information about the model.
summary(lm.fit)
# Using names() to see what other pieces of info are stored in the model.
names(lm.fit)

# Checking the model coefficients.
lm.fit$coefficients
# Altly, for model coefficients.
coef(lm.fit)
# Checking the confidence interval of model coefficients.
confint(lm.fit)

# Using predict() to predict medv values and CI for custom lstat values.
predict(
    lm.fit, 
    data.frame(lstat = c(5, 10, 15)),
    interval = "confidence"
)
# Using predict() to predict medv and prediction int for custom lstat values.
predict(
    lm.fit,
    data.frame(lstat = c(5, 10, 15)),
    interval = "prediction"
)

# Plotting medv and lstat along with the least squares regression line.
plot(Boston$lstat, Boston$medv)
abline(lm.fit)
# Increasing regression line width.
abline(lm.fit, lwd = 3)
# Changing regression line colour.
abline(lm.fit, lwd = 3, col = "red")
# Changing data point colour.
plot(Boston$lstat, Boston$medv, col = "red")
# Modifying the data points.
plot(Boston$lstat, Boston$medv, pch = 20)
plot(Boston$lstat, Boston$medv, pch = "+")
plot(1:20, 1:20, pch = 1:20)

# We examine some diagnostic plots of the regression model.
# NOTE: 4 diagnostic plots are produced by plot() fn directly from lm(). 
# We split the plot panel in 4 parts using below code.
par(mfrow = c(2, 2))
# Seeing all 4 diagnostic plots together.
plot(lm.fit)

# Altly, we can graph the diagnostic plots manually.
predval <- predict(lm.fit)
resdval <- residuals(lm.fit)
rstudval <- rstudent(lm.fit)
# Undivding the plot area to view single plot.
par(mfrow = c(1, 1))
plot(predval, resdval)
# plot(predict(lm.fit), residuals(lm.fit))
plot(predval, rstudval)
# plot(predict(lm.fit), rstudent(lm.fit))

# Compute and plot Leverage statistics for predictors.
plot(hatvalues(lm.fit))
# Identifying the index of the largest element of a vector.
which.max(hatvalues(lm.fit))


# 3.6.3 Multiple Linear Regression ----------------------------------------
# We use the same lm() to fit multiple linear regerssion models.
mulmodel <- lm(formula = medv ~ lstat + age, data = Boston)
summary(mulmodel)

# The Boston data set contains 12 variables. 
# It would be cumbersome to type all the vars to perform a regression.
# Instead, we can use the following short-hand.
mulmodel <- lm(formula = medv ~ ., data = Boston)
summary(mulmodel)

# We can access the individual components of a summary object by name.
?summary.lm
# See the R-square
summary(mulmodel)$r.squared
# See the RSE
summary(mulmodel)$sigma

# Computing the variance infation factors. 
vif(mulmodel)

# Run regression model using all of the variables but "age" var.
mulmodel1 <- lm(medv ~ . - age, data = Boston)
summary(mulmodel1)
# Altly, we use update() to modify existing model.
mulmodel2 <- update(mulmodel, ~ . - age)
summary(mulmodel2)


# 3.6.4 Interaction Terms -------------------------------------------------


# TBC ====


#_====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C11 - Survival Analysis and Censored Data -------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Lab on Survival Analysis.
# In this lab, we perform survival analyses on three separate datasets.
# Let's go!

# We load the required pkg which also contains the example datasets.
library(ISLR2)
library(survival)


# 11.8.1 Brain cancer data ------------------------------------------------
# We begin with the BrainCancer dataset.
names(BrainCancer)
glimpse(BrainCancer)
help("BrainCancer")
# Description: A dataset consisting of survival times for patients 
# diagnosed with brain cancer.




# TBC ====






























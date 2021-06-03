library(ggplot2)
library(rcompanion)
library(car) 
library(lmtest)
library(ggiraphExtra)
library(broom)
library(ggfortify)
library(ggnewscale)
library(RColorBrewer)
library(nplr) #this package fits logistic functions and allows you to fit parameters, like the asymptotes.
library(drc) #this package also fits 3 parameter logistic functions; dose response curves.
library(medrc) #this package extends the 'drc' package for mixed effects models, and does this via 'nlme.' 

survival_data <- read.delim(here::here("ox_stress/survival/Raw_survival_data.txt"))
survival_data <- survival_data[,1:6] #removes extra columns that, for some reason, were imported in from the .txt file.
str(survival_data)
# 'data.frame':	933 obs. of  6 variables:
# $ temperature: int  20 20 20 20 20 25 25 25 25 30 ...
# $ genotype   : chr  "CH" "CH" "CH" "CH" ...
# $ region     : chr  "tropical" "tropical" "tropical" "tropical" ...
# $ eggs       : int  20 20 20 20 20 20 20 20 20 20 ...
# $ hatched    : num  20 17 20 19 19 17 15 12 20 12 ...
# $ survival   : num  1 0.85 1 0.95 0.95 0.85 0.75 0.6 1 0.6 ...

survival_data$genotype <- as.factor(survival_data$genotype)
survival_data$region <- as.factor(survival_data$region)
survival_data$hatched <- as.integer(survival_data$hatched)
str(survival_data)
# 'data.frame':	933 obs. of  6 variables:
# $ temperature: int  20 20 20 20 20 25 25 25 25 30 ...
# $ genotype   : Factor w/ 26 levels "BEA_16","BEA_21",..: 6 6 6 6 6 6 6 6 6 6 ...
# $ region     : Factor w/ 4 levels "temperate","temperate.X.tropical",..: 3 3 3 3 3 3 3 3 3 3 ...
# $ eggs       : int  20 20 20 20 20 20 20 20 20 20 ...
# $ hatched    : int  20 17 20 19 19 17 15 12 20 12 ...
# $ survival   : num  1 0.85 1 0.95 0.95 0.85 0.75 0.6 1 0.6 ...


# TSO Quick data for exp
require(tidyverse)

df <- survival_data %>% 
  filter(genotype %in% c("MU", "CH", "VTECK_12", "VTECK_10")) %>%
  group_by(temperature, genotype) %>%
  summarise(survival_avg = sum(hatched)/sum(eggs)) %>%
  arrange(genotype, temperature) %>%
  filter(temperature > 28 & temperature < 36)

df %>%
  ggplot() +
  geom_line(aes(x = temperature, y = survival_avg, color = genotype)) +
  theme_classic()
# ----------

#Mixed effects logistic regression via the 'medrc' package...I can't get this to work! I think this is because,
#even though it makes logical sense to model genotype as a random effect, this dataset is unbalanced in terms
#of replication across temperatures and replication of number of lines per region (only 5 lines for tropics). It also
#could be that this package is not particularly user friendly; it is straight from the developer's Github. It's also
#tricky to run nonlinear mixed effects models to begin with.
#Thus, I think the best approach is to implement dose-response curve analysis in the 'drc' package, treating factors 
#as fixed effects. See below far...

drc.me.1 <- medrm(survival~temperature,  
      data = subset(survival_data, region == "tropical" | region == "temperate"), 
      fct=LL.3(), random = b + d + e ~ 1|genotype, start =  c(47,.9,35))
summary(drc.me.1)
# Nonlinear mixed-effects model fit by maximum likelihood
# Model: survival ~ meLL.3(temperature, b, d, e) 
# Data: data 
#       AIC       BIC   logLik
# -935.1433 -888.0747 477.5717
# 
# Random effects:
#   Formula: list(b ~ 1, d ~ 1, e ~ 1)
# Level: genotype
# Structure: General positive-definite, Log-Cholesky parametrization
#          StdDev      Corr         
# b        18.48611936 b      d     
# d         0.06759829 -0.032       
# e         0.48975223  0.963 -0.295
# Residual  0.12765249              
# 
# Fixed effects: b + d + e ~ 1 
#      Value Std.Error  DF  t-value p-value
# b 74.16965  5.143291 792  14.4207       0
# d  0.85978  0.016202 792  53.0663       0
# e 35.04129  0.109523 792 319.9446       0
# Correlation: 
#   b      d     
# d -0.142       
# e  0.702 -0.333
# 
# Standardized Within-Group Residuals:
#         Min          Q1         Med          Q3         Max 
# -3.27154051 -0.32161916  0.04305892  0.63904482  7.81854920 
# 
# Number of Observations: 818
# Number of Groups: 24 
plot(drc.me.1, ndose=250, ranef=TRUE) + theme_classic()

drc.me.2 <- medrm(survival~temperature,  
                  data = subset(survival_data, region == "tropical" | region == "temperate"), 
                  fct=LL.3(), random =  b + d + e ~ 1|region, start =  c(47,.9,35))
summary(drc.me.2)
# Nonlinear mixed-effects model fit by maximum likelihood
# Model: survival ~ meLL.3(temperature, b, d, e) 
# Data: data 
# AIC       BIC   logLik
# -827.3913 -780.3227 423.6957
# 
# Random effects:
#   Formula: list(b ~ 1, d ~ 1, e ~ 1)
# Level: region
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev       Corr   
# b        1.190905e+01 b    d 
# d        5.244123e-04  1     
# e        4.485635e-01 -1   -1
# Residual 1.434490e-01        
# 
# Fixed effects: b + d + e ~ 1 
# Value Std.Error  DF   t-value p-value
# b 53.93608  8.923913 814   6.04399       0
# d  0.85677  0.008617 814  99.42865       0
# e 35.27174  0.321908 814 109.57084       0
# Correlation: 
#   b      d     
# d -0.069       
# e -0.934 -0.110
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -4.25908209 -0.38534532  0.08797917  0.66911165  6.95176937 
# 
# Number of Observations: 818
# Number of Groups: 2
plot(drc.me.2, ndose=250, ranef=TRUE) + theme_classic()

drc.meta.1 <- metadrm(hatched/eggs ~ temperature, data = subset(survival_data, region == "tropical" | region == "temperate"),
                      fct = LL.3(), ind = genotype, cid2 = region, struct = "UN", type = "binomial")
summary(drc.meta.1)
# Two-stage meta-analysis dose-response model
# Model fitted: Log-logistic (ED50 as parameter) with lower limit at 0
# 
# Call:
#   metadrm(formula = hatched/eggs ~ temperature, fct = LL.3(), ind = genotype, 
#           data = subset(survival_data, region == "tropical" | region == 
#                           "temperate"), type = "binomial", cid2 = region, struct = "UN")
# 
# Variance estimates:
#   estim    sqrt
# tau^2.1    0.0000  0.0000
# tau^2.2    0.0000  0.0000
# tau^2.3    0.0000  0.0000
# 
# rho.b:(I  rho.d:(I  rho.e:(I
#                            b:(Intercept)         1    0.3612    0.7103
#                            d:(Intercept)    0.3612         1   -0.3684
#                            e:(Intercept)    0.7103   -0.3684         1
#                            
#                            
#                            Coefficients:
#                                           Estimate   Std.Err  z value  Pr(>|z|)    
#                              b:temperate 34.455156  4.168664   8.2653 < 2.2e-16 ***
#                              b:tropical  28.348034  6.467104   4.3834 1.168e-05 ***
#                              d:temperate  0.925017  0.019150  48.3037 < 2.2e-16 ***
#                              d:tropical   0.869210  0.042102  20.6451 < 2.2e-16 ***
#                              e:temperate 34.944252  0.173175 201.7861 < 2.2e-16 ***
#                              e:tropical  35.890170  0.338679 105.9712 < 2.2e-16 ***
#                              ---
#                              Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
drc.me.3 <- medrm(survival~temperature,  curveid = b + d + e ~ region,
                  data = subset(survival_data, region == "tropical" | region == "temperate"), 
                  fct=LL.3(), random = b + d + e ~ 1|genotype, start = c(34,28,.92,.85,35,36),
                  control = nlmeControl(tolerance = .1, pnlsTol=1, msMaxIter = 5000))
#Error in indmat %*% vcov(fmmixed) : non-conformable arguments

#####BELOW IS THE BEST APPROACH!!!##########
#Fit a three parameter logistic model with the drm() funnction in the "drc" package. Note that fct = LL.3() is for a fixed lower limit of zero.
#Note that the pmodels indicates whether you want to fix any of the three variables, those being the slope between the LT10 and LT90, the upper limit, and the 
#inflection point (LT50). If you want to fix any of these, then you would write ~1 in the appropriate spot, ordered as indicated in the LL.3(). As
#coded below, ~region-1 indicates that we're interested in comparing by region.
drc.2 <- drm(hatched/eggs~temperature, region, weights = eggs, data = subset(survival_data, region == "tropical" | region == "temperate"), 
             fct = LL.3(names = c("slope", "upper limit", "LT50")), type="binomial")
summary(drc.2)
# Model fitted: Log-logistic (ED50 as parameter) with lower limit at 0 (3 parms)
# 
# Parameter estimates:
#   
#                           Estimate Std. Error t-value   p-value    
#   slope:tropical        29.2585239  1.3244476  22.091 < 2.2e-16 ***
#   slope:temperate       36.7105073  0.9737284  37.701 < 2.2e-16 ***
#   upper limit:tropical   0.8624206  0.0096524  89.348 < 2.2e-16 ***
#   upper limit:temperate  0.8836792  0.0056657 155.970 < 2.2e-16 ***
#   LT50:tropical         35.8939982  0.0783013 458.409 < 2.2e-16 ***
#   LT50:temperate        34.8613019  0.0406156 858.323 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
modelFit(drc.2)
# Goodness-of-fit test
# 
#              Df Chisq value p value
# 
# DRC model   812      4149.9       0
coef(drc.2)
# slope:tropical       slope:temperate  upper limit:tropical upper limit:temperate         LT50:tropical        LT50:temperate 
#     29.2585239            36.7105073             0.8624206             0.8836792            35.8939982            34.8613019

ED(drc.2, c(50), interval = "delta") #this gives you estimates of LT50 and unadjusted 95% confidence intervals.
#                 Estimate Std. Error     Lower     Upper
# e:temperate:50 34.861302   0.040616 34.781697 34.940907
# e:tropical:50  35.893998   0.078301 35.740531 36.047466

#This calculates adjusted 95% confidence intervals, in collaboration with the multcomp package.
library(multcomp)
drc.2.LT50res <- ED(drc.2, c(50), interval = "delta",
                     multcomp = TRUE, display = T)
confint(glht(drc.2.LT50res[["EDmultcomp"]]))
# Simultaneous Confidence Intervals
# 
# Fit: NULL
# 
# Quantile = 2.2364
# 95% family-wise confidence level
# 
# 
# Linear Hypotheses:
#                     Estimate lwr     upr    
# e:temperate:50 == 0 34.8613  34.7705 34.9521
# e:tropical:50 == 0  35.8940  35.7189 36.0691

drc.2.2 <- drm(hatched/eggs~temperature, region, weights = eggs, data = subset(survival_data, region == "tropical" | region == "temperate"), 
               fct = LL.3(names = c("slope", "upper limit", "LT50")), pmodels = list(~region-1, ~region-1, ~1), type="binomial")
summary(drc.2.2)
# Model fitted: Log-logistic (ED50 as parameter) with lower limit at 0 (3 parms)
# 
# Parameter estimates:
#   
#                                         Estimate Std. Error t-value   p-value    
#   slope:factor(region)temperate       38.1280297  0.9891567  38.546 < 2.2e-16 ***
#   slope:factor(region)tropical        23.7526388  1.1207829  21.193 < 2.2e-16 ***
#   upper limit:factor(region)temperate  0.8709401  0.0056818 153.287 < 2.2e-16 ***
#   upper limit:factor(region)tropical   0.9018754  0.0077558 116.283 < 2.2e-16 ***
#   LT50:(Intercept)                    35.0570840  0.0374659 935.708 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(drc.2.2,drc.2) #this is comparing models to test for just the difference of LT50 among regions (i.e., by modeling just one LT50 in drc.2.2).
# 1st model
# fct:     LL.3(names = c("slope", "upper limit", "LT50"))
# pmodels: ~region - 1, ~region - 1, ~1
# 2nd model
# fct:     LL.3(names = c("slope", "upper limit", "LT50"))
# pmodels: region (for all parameters)
# 
# ANOVA-like table
# 
#           ModelDf  Loglik Df LR value p value
# 1st model       5 -2340.3                    
# 2nd model       6 -2273.3  1   134.08       0

#Below we are directly comparing the three parameter estimates (LT50, upper limit, and slope of line at LT50) between regions. 
#Note that "-" tells it that you want to compare differences, and the estimate is the difference in that parameter.
#You can also indicate "/", which tells it to compare parameter ratios.
compParm(drc.2, "LT50", "-")
# Comparison of parameter 'LT50' 
# 
#                    Estimate Std. Error t-value   p-value    
# tropical-temperate 1.032696   0.088208  11.707 < 2.2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

compParm(drc.2, "upper limit", "-")
# Comparison of parameter 'upper limit' 
# 
#                     Estimate Std. Error t-value p-value  
# tropical-temperate -0.021259   0.011192 -1.8994 0.05751 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

compParm(drc.2, "slope", "-")
# Comparison of parameter 'slope' 
# 
#                    Estimate Std. Error t-value  p-value    
# tropical-temperate  -7.4520     1.6439 -4.5332 5.81e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Below is another way of comparing differences among regions, via fitting a reduced null model and comparing models via LRT.
#So, basically this compares the entire survival curves.
drc.null <- drm(hatched/eggs~temperature, weights = eggs, data = subset(survival_data, region == "tropical" | region == "temperate"), 
                fct = LL.3(names = c("slope", "upper limit", "LT50")), type="binomial")

anova(drc.null,drc.2)
# 1st model
# fct:     LL.3(names = c("slope", "upper limit", "LT50"))
# pmodels: 1 (for all parameters)
# 2nd model
# fct:     LL.3(names = c("slope", "upper limit", "LT50"))
# pmodels: ~region - 1, ~region - 1, ~region - 1
# 
# ANOVA-like table
# 
# Model          Df  Loglik Df LR value p value
# 1st model       3 -2387.6                    
# 2nd model       6 -2273.3  3   228.72       0

drc.3 <- drm(hatched/eggs~temperature, genotype, weights = eggs, data = subset(survival_data, region == "tropical" | region == "temperate"), 
             fct = LL.3(names = c("slope", "upper limit", "LT50")), pmodels = list(~genotype-1, ~genotype-1, ~genotype-1), type="binomial")
summary(drc.3)
# Model fitted: Log-logistic (ED50 as parameter) with lower limit at 0 (3 parms)
# 
# Parameter estimates:
#   
#                         Estimate Std. Error  t-value   p-value    
#   slope:CH             39.505903   4.018677   9.8306 < 2.2e-16 ***
#   slope:GH             42.199269   5.374300   7.8520 4.045e-15 ***
#   slope:MU             20.209212   2.120222   9.5316 < 2.2e-16 ***
#   slope:SK             32.949505   2.746770  11.9957 < 2.2e-16 ***
#   slope:GU             28.922348   5.957246   4.8550 1.204e-06 ***
#   slope:VTECK_10       44.382398   5.698343   7.7886 6.828e-15 ***
#   slope:VTECK_12       47.531203   7.093669   6.7005 2.077e-11 ***
#   slope:VTECK_14       43.152130   4.141548  10.4193 < 2.2e-16 ***
#   slope:VTECK_2        42.910015   4.483628   9.5704 < 2.2e-16 ***
#   slope:VTECK_4        76.580374   8.906966   8.5978 < 2.2e-16 ***
#   slope:VTECK_5        38.421583   4.893132   7.8521 4.043e-15 ***
#   slope:VTECK_8        46.648058   6.717098   6.9447 3.793e-12 ***
#   slope:VTECK_9        43.159568   4.543908   9.4983 < 2.2e-16 ***
#   slope:BEA_16         31.134337   4.708593   6.6122 3.786e-11 ***
#   slope:BEA_21         40.351895   4.694713   8.5952 < 2.2e-16 ***
#   slope:BEA_32         23.294604   4.416480   5.2745 1.331e-07 ***
#   slope:BEA_36         27.661977   2.286518  12.0979 < 2.2e-16 ***
#   slope:BEA_5          58.780773   7.268522   8.0870 6.387e-16 ***
#   slope:RFM_16         21.806102   2.781721   7.8391 4.490e-15 ***
#   slope:RFM_19         36.546133   3.782244   9.6626 < 2.2e-16 ***
#   slope:RFM_34         35.584281   5.912944   6.0180 1.766e-09 ***
#   slope:RFM_4          44.515784   4.831713   9.2133 < 2.2e-16 ***
#   slope:RFM_48         39.827121   5.132241   7.7602 8.459e-15 ***
#   slope:RFM_6          39.961940   4.465904   8.9482 < 2.2e-16 ***
#   upper limit:CH        0.873273   0.018198  47.9872 < 2.2e-16 ***
#   upper limit:GH        0.755429   0.028144  26.8411 < 2.2e-16 ***
#   upper limit:MU        0.897561   0.021040  42.6605 < 2.2e-16 ***
#   upper limit:SK        0.883568   0.014779  59.7874 < 2.2e-16 ***
#   upper limit:GU        0.858134   0.091743   9.3537 < 2.2e-16 ***
#   upper limit:VTECK_10  0.939936   0.016702  56.2785 < 2.2e-16 ***
#   upper limit:VTECK_12  0.926866   0.043975  21.0773 < 2.2e-16 ***
#   upper limit:VTECK_14  0.955610   0.011463  83.3679 < 2.2e-16 ***
#   upper limit:VTECK_2   0.966335   0.010441  92.5478 < 2.2e-16 ***
#   upper limit:VTECK_4   0.966633   0.016770  57.6423 < 2.2e-16 ***
#   upper limit:VTECK_5   0.858941   0.023550  36.4730 < 2.2e-16 ***
#   upper limit:VTECK_8   0.825752   0.042097  19.6157 < 2.2e-16 ***
#   upper limit:VTECK_9   0.968462   0.010445  92.7230 < 2.2e-16 ***
#   upper limit:BEA_16    0.817745   0.028703  28.4900 < 2.2e-16 ***
#   upper limit:BEA_21    0.923553   0.016940  54.5180 < 2.2e-16 ***
#   upper limit:BEA_32    0.921895   0.083980  10.9776 < 2.2e-16 ***
#   upper limit:BEA_36    0.772503   0.018380  42.0306 < 2.2e-16 ***
#   upper limit:BEA_5     0.915851   0.025470  35.9582 < 2.2e-16 ***
#   upper limit:RFM_16    0.972976   0.080973  12.0161 < 2.2e-16 ***
#   upper limit:RFM_19    0.878875   0.019797  44.3943 < 2.2e-16 ***
#   upper limit:RFM_34    0.927899   0.039055  23.7585 < 2.2e-16 ***
#   upper limit:RFM_4     0.760235   0.032170  23.6318 < 2.2e-16 ***
#   upper limit:RFM_48    0.920222   0.019597  46.9572 < 2.2e-16 ***
#   upper limit:RFM_6     0.840943   0.033251  25.2908 < 2.2e-16 ***
#   LT50:CH              36.317847   0.143588 252.9318 < 2.2e-16 ***
#   LT50:GH              35.553656   0.148804 238.9292 < 2.2e-16 ***
#   LT50:MU              36.463928   0.234069 155.7831 < 2.2e-16 ***
#   LT50:SK              35.628248   0.131035 271.8979 < 2.2e-16 ***
#   LT50:GU              35.771367   0.497642  71.8817 < 2.2e-16 ***
#   LT50:VTECK_10        34.186988   0.168842 202.4790 < 2.2e-16 ***
#   LT50:VTECK_12        34.226944   0.183236 186.7916 < 2.2e-16 ***
#   LT50:VTECK_14        35.186276   0.136338 258.0815 < 2.2e-16 ***
#   LT50:VTECK_2         35.283804   0.141682 249.0349 < 2.2e-16 ***
#   LT50:VTECK_4         35.027196   0.121852 287.4568 < 2.2e-16 ***
#   LT50:VTECK_5         34.758982   0.190079 182.8657 < 2.2e-16 ***
#   LT50:VTECK_8         34.773832   0.194999 178.3281 < 2.2e-16 ***
#   LT50:VTECK_9         35.347546   0.150608 234.6986 < 2.2e-16 ***
#   LT50:BEA_16          34.980136   0.258495 135.3224 < 2.2e-16 ***
#   LT50:BEA_21          34.897073   0.165034 211.4536 < 2.2e-16 ***
#   LT50:BEA_32          34.935776   0.448813  77.8404 < 2.2e-16 ***
#   LT50:BEA_36          35.357653   0.174524 202.5947 < 2.2e-16 ***
#   LT50:BEA_5           35.072319   0.140477 249.6652 < 2.2e-16 ***
#   LT50:RFM_16          33.572025   0.385148  87.1665 < 2.2e-16 ***
#   LT50:RFM_19          34.797191   0.169913 204.7936 < 2.2e-16 ***
#   LT50:RFM_34          35.046900   0.231806 151.1903 < 2.2e-16 ***
#   LT50:RFM_4           34.866853   0.154836 225.1862 < 2.2e-16 ***
#   LT50:RFM_48          35.193429   0.196533 179.0710 < 2.2e-16 ***
#   LT50:RFM_6           34.859729   0.173672 200.7211 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
modelFit(drc.3)
# Goodness-of-fit test
# 
#              Df Chisq value p value
# 
# DRC model   746      3002.9       0

anova(drc.null,drc.3)
# 1st model
# fct:     LL.3(names = c("slope", "upper limit", "LT50"))
# pmodels: 1 (for all parameters)
# 2nd model
# fct:     LL.3(names = c("slope", "upper limit", "LT50"))
# pmodels: ~genotype - 1, ~genotype - 1, ~genotype - 1
# 
# ANOVA-like table
# 
# Model          Df  Loglik Df LR value p value
# 1st model       3 -2387.6                    
# 2nd model      72 -1976.5 69    822.2       0
compParm(drc.3, "LT50", "-")
# Comparison of parameter 'LT50' 
# 
#                       Estimate Std. Error t-value   p-value    
#   CH-GH              0.7641917  0.2067850  3.6956 0.0002194 ***
#   CH-MU             -0.1460810  0.2746006 -0.5320 0.5947424    
#   CH-SK              0.6895994  0.1943905  3.5475 0.0003889 ***
#   CH-GU              0.5464807  0.5179432  1.0551 0.2913806    
#   CH-VTECK_10        2.1308589  0.2216417  9.6140 < 2.2e-16 ***
#   CH-VTECK_12        2.0909032  0.2327935  8.9818 < 2.2e-16 ***
#   CH-VTECK_14        1.1315717  0.1980035  5.7149 1.098e-08 ***
#   CH-VTECK_2         1.0340438  0.2017206  5.1261 2.958e-07 ***
#   CH-VTECK_4         1.2906518  0.1883223  6.8534 7.211e-12 ***
#   CH-VTECK_5         1.5588649  0.2382174  6.5439 5.994e-11 ***
#   CH-VTECK_8         1.5440155  0.2421612  6.3760 1.818e-10 ***
#   CH-VTECK_9         0.9703012  0.2080870  4.6630 3.117e-06 ***
#   CH-BEA_16          1.3377118  0.2956973  4.5239 6.070e-06 ***
#   CH-BEA_21          1.4207746  0.2187548  6.4948 8.313e-11 ***
#   CH-BEA_32          1.3820712  0.4712222  2.9329 0.0033576 ** 
#   CH-BEA_36          0.9601941  0.2260001  4.2486 2.151e-05 ***
#   CH-BEA_5           1.2455285  0.2008763  6.2005 5.629e-10 ***
#   CH-RFM_16          2.7458223  0.4110434  6.6801 2.387e-11 ***
#   CH-RFM_19          1.5206559  0.2224589  6.8357 8.162e-12 ***
#   CH-RFM_34          1.2709476  0.2726749  4.6610 3.146e-06 ***
#   CH-RFM_4           1.4509942  0.2111669  6.8713 6.361e-12 ***
#   CH-RFM_48          1.1244185  0.2433984  4.6197 3.844e-06 ***
#   CH-RFM_6           1.4581187  0.2253431  6.4707 9.757e-11 ***
#   GH-MU             -0.9102727  0.2773639 -3.2819 0.0010312 ** 
#   GH-SK             -0.0745923  0.1982749 -0.3762 0.7067635    
#   GH-GU             -0.2177110  0.5194136 -0.4191 0.6751082    
#   GH-VTECK_10        1.3666671  0.2250563  6.0726 1.259e-09 ***
#   GH-VTECK_12        1.3267114  0.2360468  5.6205 1.904e-08 ***
#   GH-VTECK_14        0.3673800  0.2018184  1.8203 0.0687059 .  
#   GH-VTECK_2         0.2698520  0.2054666  1.3134 0.1890609    
#   GH-VTECK_4         0.5264601  0.1923294  2.7373 0.0061949 ** 
#   GH-VTECK_5         0.7946732  0.2413976  3.2920 0.0009949 ***
#   GH-VTECK_8         0.7798238  0.2452904  3.1792 0.0014769 ** 
#   GH-VTECK_9         0.2061095  0.2117204  0.9735 0.3303055    
#   GH-BEA_16          0.5735201  0.2982653  1.9229 0.0544987 .  
#   GH-BEA_21          0.6565829  0.2222138  2.9547 0.0031294 ** 
#   GH-BEA_32          0.6178795  0.4728379  1.3067 0.1912987    
#   GH-BEA_36          0.1960024  0.2293498  0.8546 0.3927724    
#   GH-BEA_5           0.4813367  0.2046377  2.3521 0.0186657 *  
#   GH-RFM_16          1.9816306  0.4128946  4.7994 1.592e-06 ***
#   GH-RFM_19          0.7564642  0.2258611  3.3492 0.0008103 ***
#   GH-RFM_34          0.5067558  0.2754577  1.8397 0.0658142 .  
#   GH-RFM_4           0.6868025  0.2147481  3.1982 0.0013830 ** 
#   GH-RFM_48          0.3602267  0.2465118  1.4613 0.1439342    
#   GH-RFM_6           0.6939269  0.2287024  3.0342 0.0024118 ** 
#   MU-SK              0.8356804  0.2682506  3.1153 0.0018376 ** 
#   MU-GU              0.6925617  0.5499417  1.2593 0.2079088    
#   MU-VTECK_10        2.2769399  0.2886100  7.8893 3.073e-15 ***
#   MU-VTECK_12        2.2369842  0.2972600  7.5253 5.260e-14 ***
#   MU-VTECK_14        1.2776527  0.2708802  4.7167 2.397e-06 ***
#   MU-VTECK_2         1.1801248  0.2736091  4.3132 1.609e-05 ***
#   MU-VTECK_4         1.4367328  0.2638863  5.4445 5.195e-08 ***
#   MU-VTECK_5         1.7049459  0.3015265  5.6544 1.564e-08 ***
#   MU-VTECK_8         1.6900965  0.3046519  5.5476 2.896e-08 ***
#   MU-VTECK_9         1.1163822  0.2783360  4.0109 6.048e-05 ***
#   MU-BEA_16          1.4837928  0.3487228  4.2549 2.091e-05 ***
#   MU-BEA_21          1.5668556  0.2863990  5.4709 4.478e-08 ***
#   MU-BEA_32          1.5281522  0.5061829  3.0190 0.0025363 ** 
#   MU-BEA_36          1.1062752  0.2919704  3.7890 0.0001513 ***
#   MU-BEA_5           1.3916095  0.2729872  5.0977 3.438e-07 ***
#   MU-RFM_16          2.8919033  0.4506965  6.4165 1.394e-10 ***
#   MU-RFM_19          1.6667369  0.2892381  5.7625 8.287e-09 ***
#   MU-RFM_34          1.4170286  0.3294273  4.3015 1.697e-05 ***
#   MU-RFM_4           1.5970752  0.2806460  5.6907 1.265e-08 ***
#   MU-RFM_48          1.2704995  0.3056362  4.1569 3.226e-05 ***
#   MU-RFM_6           1.6041997  0.2914622  5.5040 3.713e-08 ***
#   SK-GU             -0.1431187  0.5146048 -0.2781 0.7809250    
#   SK-VTECK_10        1.4412594  0.2137240  6.7436 1.546e-11 ***
#   SK-VTECK_12        1.4013037  0.2252680  6.2206 4.952e-10 ***
#   SK-VTECK_14        0.4419723  0.1890986  2.3373 0.0194258 *  
#   SK-VTECK_2         0.3444443  0.1929873  1.7848 0.0742933 .  
#   SK-VTECK_4         0.6010524  0.1789363  3.3590 0.0007822 ***
#   SK-VTECK_5         0.8692655  0.2308688  3.7652 0.0001664 ***
#   SK-VTECK_8         0.8544161  0.2349361  3.6368 0.0002760 ***
#   SK-VTECK_9         0.2807018  0.1996324  1.4061 0.1596965    
#   SK-BEA_16          0.6481124  0.2898100  2.2363 0.0253298 *  
#   SK-BEA_21          0.7311752  0.2107287  3.4697 0.0005209 ***
#   SK-BEA_32          0.6924718  0.4675503  1.4811 0.1385896    
#   SK-BEA_36          0.2705947  0.2182405  1.2399 0.2150153    
#   SK-BEA_5           0.5559290  0.1921046  2.8939 0.0038050 ** 
#   SK-RFM_16          2.0562229  0.4068287  5.0543 4.320e-07 ***
#   SK-RFM_19          0.8310565  0.2145713  3.8731 0.0001075 ***
#   SK-RFM_34          0.5813481  0.2662790  2.1832 0.0290189 *  
#   SK-RFM_4           0.7613948  0.2028407  3.7537 0.0001743 ***
#   SK-RFM_48          0.4348190  0.2362111  1.8408 0.0656498 .  
#   SK-RFM_6           0.7685192  0.2175601  3.5324 0.0004117 ***
#   GU-VTECK_10        1.5843781  0.5255050  3.0150 0.0025701 ** 
#   GU-VTECK_12        1.5444224  0.5303048  2.9123 0.0035874 ** 
#   GU-VTECK_14        0.5850910  0.5159804  1.1339 0.2568196    
#   GU-VTECK_2         0.4875630  0.5174182  0.9423 0.3460393    
#   GU-VTECK_4         0.7441711  0.5123434  1.4525 0.1463668    
#   GU-VTECK_5         1.0123842  0.5327081  1.9004 0.0573743 .  
#   GU-VTECK_8         0.9975348  0.5344834  1.8664 0.0619920 .  
#   GU-VTECK_9         0.4238205  0.5199333  0.8151 0.4149899    
#   GU-BEA_16          0.7912311  0.5607739  1.4110 0.1582556    
#   GU-BEA_21          0.8742939  0.5242939  1.6676 0.0954022 .  
#   GU-BEA_32          0.8355905  0.6701349  1.2469 0.2124346    
#   GU-BEA_36          0.4137134  0.5273580  0.7845 0.4327456    
#   GU-BEA_5           0.6990477  0.5170897  1.3519 0.1764109    
#   GU-RFM_16          2.1993416  0.6292751  3.4950 0.0004740 ***
#   GU-RFM_19          0.9741752  0.5258502  1.8526 0.0639438 .  
#   GU-RFM_34          0.7244668  0.5489827  1.3197 0.1869508    
#   GU-RFM_4           0.9045135  0.5211736  1.7355 0.0826466 .  
#   GU-RFM_48          0.5779377  0.5350450  1.0802 0.2800680    
#   GU-RFM_6           0.9116379  0.5270768  1.7296 0.0836998 .  
#   VTECK_10-VTECK_12 -0.0399557  0.2491648 -0.1604 0.8725987    
#   VTECK_10-VTECK_14 -0.9992872  0.2170154 -4.6047 4.131e-06 ***
#   VTECK_10-VTECK_2  -1.0968151  0.2204121 -4.9762 6.484e-07 ***
#   VTECK_10-VTECK_4  -0.8402071  0.2082200 -4.0352 5.456e-05 ***
#   VTECK_10-VTECK_5  -0.5719940  0.2542397 -2.2498 0.0244603 *  
#   VTECK_10-VTECK_8  -0.5868433  0.2579387 -2.2751 0.0228983 *  
#   VTECK_10-VTECK_9  -1.1605576  0.2262532 -5.1295 2.906e-07 ***
#   VTECK_10-BEA_16   -0.7931471  0.3087510 -2.5689 0.0102025 *  
#   VTECK_10-BEA_21   -0.7100843  0.2361016 -3.0075 0.0026337 ** 
#   VTECK_10-BEA_32   -0.7487876  0.4795213 -1.5615 0.1183984    
#   VTECK_10-BEA_36   -1.1706647  0.2428298 -4.8209 1.429e-06 ***
#   VTECK_10-BEA_5    -0.8853304  0.2196396 -4.0308 5.558e-05 ***
#   VTECK_10-RFM_16    0.6149634  0.4205317  1.4623 0.1436461    
#   VTECK_10-RFM_19   -0.6102030  0.2395376 -2.5474 0.0108522 *  
#   VTECK_10-RFM_34   -0.8599113  0.2867785 -2.9985 0.0027129 ** 
#   VTECK_10-RFM_4    -0.6798646  0.2290890 -2.9677 0.0030005 ** 
#   VTECK_10-RFM_48   -1.0064404  0.2591005 -3.8844 0.0001026 ***
#   VTECK_10-RFM_6    -0.6727402  0.2422185 -2.7774 0.0054794 ** 
#   VTECK_12-VTECK_14 -0.9593315  0.2283931 -4.2004 2.665e-05 ***
#   VTECK_12-VTECK_2  -1.0568594  0.2316231 -4.5628 5.047e-06 ***
#   VTECK_12-VTECK_4  -0.8002514  0.2200530 -3.6366 0.0002762 ***
#   VTECK_12-VTECK_5  -0.5320383  0.2640181 -2.0152 0.0438881 *  
#   VTECK_12-VTECK_8  -0.5468876  0.2675820 -2.0438 0.0409720 *  
#   VTECK_12-VTECK_9  -1.1206019  0.2371882 -4.7245 2.307e-06 ***
#   VTECK_12-BEA_16   -0.7531914  0.3168516 -2.3771 0.0174489 *  
#   VTECK_12-BEA_21   -0.6701286  0.2466003 -2.7175 0.0065783 ** 
#   VTECK_12-BEA_32   -0.7088319  0.4847767 -1.4622 0.1436912    
#   VTECK_12-BEA_36   -1.1307090  0.2530495 -4.4683 7.883e-06 ***
#   VTECK_12-BEA_5    -0.8453747  0.2308881 -3.6614 0.0002508 ***
#   VTECK_12-RFM_16    0.6549191  0.4265146  1.5355 0.1246576    
#   VTECK_12-RFM_19   -0.5702473  0.2498920 -2.2820 0.0224908 *  
#   VTECK_12-RFM_34   -0.8199556  0.2954821 -2.7750 0.0055206 ** 
#   VTECK_12-RFM_4    -0.6399089  0.2398948 -2.6675 0.0076428 ** 
#   VTECK_12-RFM_48   -0.9664847  0.2687021 -3.5969 0.0003221 ***
#   VTECK_12-RFM_6    -0.6327845  0.2524629 -2.5064 0.0121952 *  
#   VTECK_14-VTECK_2  -0.0975279  0.1966261 -0.4960 0.6198895    
#   VTECK_14-VTECK_4   0.1590801  0.1828549  0.8700 0.3843115    
#   VTECK_14-VTECK_5   0.4272932  0.2339191  1.8267 0.0677493 .  
#   VTECK_14-VTECK_8   0.4124438  0.2379342  1.7334 0.0830181 .  
#   VTECK_14-VTECK_9  -0.1612705  0.2031523 -0.7938 0.4272884    
#   VTECK_14-BEA_16    0.2061401  0.2922457  0.7054 0.4805826    
#   VTECK_14-BEA_21    0.2892029  0.2140661  1.3510 0.1766960    
#   VTECK_14-BEA_32    0.2504995  0.4690640  0.5340 0.5933130    
#   VTECK_14-BEA_36   -0.1713776  0.2214648 -0.7738 0.4390275    
#   VTECK_14-BEA_5     0.1139568  0.1957598  0.5821 0.5604822    
#   VTECK_14-RFM_16    1.6142506  0.4085674  3.9510 7.782e-05 ***
#   VTECK_14-RFM_19    0.3890842  0.2178499  1.7860 0.0740961 .  
#   VTECK_14-RFM_34    0.1393759  0.2689280  0.5183 0.6042736    
#   VTECK_14-RFM_4     0.3194225  0.2063058  1.5483 0.1215510    
#   VTECK_14-RFM_48   -0.0071532  0.2391932 -0.0299 0.9761423    
#   VTECK_14-RFM_6     0.3265470  0.2207943  1.4790 0.1391498    
#   VTECK_2-VTECK_4    0.2566080  0.1868736  1.3732 0.1697015    
#   VTECK_2-VTECK_5    0.5248212  0.2370738  2.2137 0.0268462 *  
#   VTECK_2-VTECK_8    0.5099718  0.2410363  2.1157 0.0343664 *  
#   VTECK_2-VTECK_9   -0.0637425  0.2067769 -0.3083 0.7578790    
#   VTECK_2-BEA_16     0.3036681  0.2947768  1.0302 0.3029337    
#   VTECK_2-BEA_21     0.3867308  0.2175089  1.7780 0.0754038 .  
#   VTECK_2-BEA_32     0.3480275  0.4706451  0.7395 0.4596223    
#   VTECK_2-BEA_36    -0.0738496  0.2247943 -0.3285 0.7425179    
#   VTECK_2-BEA_5      0.2114847  0.1995188  1.0600 0.2891564    
#   VTECK_2-RFM_16     1.7117786  0.4103817  4.1712 3.030e-05 ***
#   VTECK_2-RFM_19     0.4866121  0.2212338  2.1995 0.0278397 *  
#   VTECK_2-RFM_34     0.2369038  0.2716764  0.8720 0.3832044    
#   VTECK_2-RFM_4      0.4169505  0.2098759  1.9867 0.0469610 *  
#   VTECK_2-RFM_48     0.0903747  0.2422792  0.3730 0.7091344    
#   VTECK_2-RFM_6      0.4240749  0.2241338  1.8921 0.0584828 .  
#   VTECK_4-VTECK_5    0.2682131  0.2257832  1.1879 0.2348637    
#   VTECK_4-VTECK_8    0.2533637  0.2299405  1.1019 0.2705196    
#   VTECK_4-VTECK_9   -0.3203506  0.1937286 -1.6536 0.0982078 .  
#   VTECK_4-BEA_16     0.0470600  0.2857752  0.1647 0.8691998    
#   VTECK_4-BEA_21     0.1301228  0.2051444  0.6343 0.5258859    
#   VTECK_4-BEA_32     0.0914194  0.4650601  0.1966 0.8441597    
#   VTECK_4-BEA_36    -0.3304576  0.2128534 -1.5525 0.1205395    
#   VTECK_4-BEA_5     -0.0451233  0.1859619 -0.2426 0.8082779    
#   VTECK_4-RFM_16     1.4551705  0.4039644  3.6022 0.0003155 ***
#   VTECK_4-RFM_19     0.2300041  0.2090897  1.1000 0.2713208    
#   VTECK_4-RFM_34    -0.0197042  0.2618820 -0.0752 0.9400231    
#   VTECK_4-RFM_4      0.1603424  0.1970330  0.8138 0.4157683    
#   VTECK_4-RFM_48    -0.1662333  0.2312430 -0.7189 0.4722218    
#   VTECK_4-RFM_6      0.1674669  0.2121557  0.7894 0.4299026    
#   VTECK_5-VTECK_8   -0.0148494  0.2723138 -0.0545 0.9565126    
#   VTECK_5-VTECK_9   -0.5885637  0.2425139 -2.4269 0.0152273 *  
#   VTECK_5-BEA_16    -0.2211531  0.3208577 -0.6893 0.4906622    
#   VTECK_5-BEA_21    -0.1380903  0.2517269 -0.5486 0.5832993    
#   VTECK_5-BEA_32    -0.1767937  0.4874045 -0.3627 0.7168105    
#   VTECK_5-BEA_36    -0.5986708  0.2580480 -2.3200 0.0203410 *  
#   VTECK_5-BEA_5     -0.3133364  0.2363558 -1.3257 0.1849396    
#   VTECK_5-RFM_16     1.1869574  0.4294990  2.7636 0.0057170 ** 
#   VTECK_5-RFM_19    -0.0382090  0.2549524 -0.1499 0.8808693    
#   VTECK_5-RFM_34    -0.2879173  0.2997739 -0.9604 0.3368296    
#   VTECK_5-RFM_4     -0.1078707  0.2451616 -0.4400 0.6599384    
#   VTECK_5-RFM_48    -0.4344464  0.2734146 -1.5890 0.1120681    
#   VTECK_5-RFM_6     -0.1007462  0.2574728 -0.3913 0.6955837    
#   VTECK_8-VTECK_9   -0.5737143  0.2463890 -2.3285 0.0198861 *  
#   VTECK_8-BEA_16    -0.2063037  0.3237966 -0.6371 0.5240337    
# [ reached getOption("max.print") -- omitted 76 rows ]
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
LT50.comp <- as.data.frame(compParm(drc.3, "LT50", "-"))
write.table(LT50.comp, "LT50_pairwise_comparisons.txt", sep = "\t")

parameter.estimates.region <- as.data.frame(drc.2$coefficients)
parameter.estimates.genotype <- as.data.frame(drc.3$coefficients)
write.table(parameter.estimates.region, "parameters_est_region.txt", sep = "\t")
write.table(parameter.estimates.genotype, "parameters_est_genotype.txt", sep = "\t")

#Now compare genotype and region effects...
anova(drc.3,drc.2)
# 1st model
# fct:     LL.3(names = c("slope", "upper limit", "LT50"))
# pmodels: ~genotype - 1, ~genotype - 1, ~genotype - 1
# 2nd model
# fct:     LL.3(names = c("slope", "upper limit", "LT50"))
# pmodels: ~region - 1, ~region - 1, ~region - 1
# 
# ANOVA-like table
# 
#           ModelDf  Loglik Df LR value p value
# 1st model      72 -1976.5                    
# 2nd model       6 -2273.3 66   593.48       0

AIC(drc.2)
#[1] 4558.553
AIC(drc.3)
#[1] 4097.07
AIC(drc.null)
#[1] 4781.269

#Base package plotting...
plot(drc.2, ylim = c(0, 1), xlim = c(20,45), type = "all", xlab = "Temperature °C", ylab = "Proportion hatched", legendPos = c(28,.4))

plot(drc.null, ylim = c(0, 1), xlim = c(20,45), type = "all", xlab = "Temperature °C", ylab = "Proportion hatched")

plot(drc.3, ylim = c(0, 1), xlim = c(20,45), type = "all", xlab = "Temperature °C", ylab = "Proportion hatched", legend = F)

#ggplot-ting
#First plot the logistic regressions by region...I like this plot best, as it reveals the difference between regions,
#and hence gets at the core question of the study.

#new temperature sequence as support for the line
newdata <- expand.grid(temperature=seq(20, 45, length=250), region=c("tropical","temperate"))
# predictions and confidence intervals
pm <- predict(drc.2, newdata=newdata, interval="confidence", type="response")
# new data with predictions
newdata$p <- pm[,1]
newdata$pmin <- pm[,2]
newdata$pmax <- pm[,3]
# plotting the curves of logistic regressions by region
ggplot(subset(survival_data, region == "tropical" | region == "temperate"), aes(x = temperature, y = survival)) +
  geom_ribbon(data=newdata, aes(x=temperature, y=p, ymin=pmin, ymax=pmax, fill=factor(region, levels = c("temperate","tropical"))), alpha=0.5) +
  #geom_line(data=newdata, size = 0.7, aes(x=temperature, y=p, color = factor(region, levels = c("temperate","tropical")), linetype = factor(region, levels = c("temperate","tropical")))) +
  geom_point(aes(color = factor(region, levels = c("temperate","tropical")), shape = factor(region, levels = c("temperate","tropical"))), alpha=0.2, size = 2) +
  geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.5) +
  xlab("Temperature (°C)") + 
  ylab("Proportion hatched") +
  scale_color_manual(values = c("#009E73", "#D55E00")) +
  #scale_shape_manual(values = c(1,4)) +
  scale_fill_manual(values = c("#009E73","#D55E00")) +
  xlim(25,45) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), legend.title = element_blank(), legend.position = c(0.85,0.88)) + 
  guides(colour = guide_legend(override.aes = list(alpha = .6)))

#Now make plots of the logistic regressions by genotype...
#new temperature sequence as support for the line
surv_dat <- na.omit(read.delim("Parental_survival_data.txt")) #parental lines survival data
surv_dat_trop <- na.omit(read.delim("tropical_survival_data.txt")) #tropical parental lines survival data
surv_dat_temp <- na.omit(read.delim("temperate_survival_data.txt")) #temperate parental lines survival data

surv_dat$genotype <- as.factor(surv_dat$genotype)
surv_dat$region <- as.factor(surv_dat$region)
surv_dat$hatched <- as.integer(surv_dat$hatched)
surv_dat_trop$genotype <- as.factor(surv_dat_trop$genotype)
surv_dat_trop$region <- as.factor(surv_dat_trop$region)
surv_dat_trop$hatched <- as.integer(surv_dat_trop$hatched)
surv_dat_temp$genotype <- as.factor(surv_dat_temp$genotype)
surv_dat_temp$region <- as.factor(surv_dat_temp$region)
surv_dat_temp$region <- as.factor(surv_dat_temp$population)
surv_dat_temp$hatched <- as.integer(surv_dat_temp$hatched)

newdata.trop <- expand.grid(temperature=seq(20, 45, length=250), genotype=c(levels(surv_dat_trop$genotype)))
newdata.temp <- expand.grid(temperature=seq(20, 45, length=250), genotype=c(levels(surv_dat_temp$genotype)))

# predictions and confidence intervals
pm.trop <- predict(drc.3, newdata=newdata.trop, interval="confidence", type="response")
pm.temp <- predict(drc.3, newdata=newdata.temp, interval="confidence", type="response")


# new data with predictions
newdata.trop$p <- pm.trop[,1]
newdata.trop$pmin <- pm.trop[,2]
newdata.trop$pmax <- pm.trop[,3]

newdata.temp$p <- pm.temp[,1]
newdata.temp$pmin <- pm.temp[,2]
newdata.temp$pmax <- pm.temp[,3]

# plotting the curve of the logistic regressions by genotype
ggplot(surv_dat, aes(x = temperature, y = survival)) +
  #geom_point(aes(color = genotype), alpha=0.5, size = 2) +
  geom_line(data=newdata.temp, size = .7, alpha=0.5, aes(x=temperature, y=p, color=genotype), linetype="dotted") +
  scale_color_grey() +
  new_scale_color() + #this resets the color scale to have independent color scales for multiple plotting features
  geom_line(data=newdata.trop, size = 0.7, alpha=0.7, aes(x=temperature, y=p, color= genotype)) +
  scale_color_brewer(palette = "Dark2") +  
  geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.5) +
  xlim(25,45) +
  xlab("Temperature (°C)") + 
  ylab("Proportion hatched") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), legend.title = element_blank(), legend.position = "none") 

#Make a combined plot of the logistic regressions by region and by genotype...

ggplot(surv_dat, aes(x = temperature, y = survival)) +
  geom_line(data=newdata.temp, size = .7, alpha=0.3, aes(x=temperature, y=p, color=genotype), linetype="dotted") +
  scale_color_grey() +
  new_scale_color() + #this resets the color scale to have independent color scales for multiple plotting features
  geom_line(data=newdata.trop, size = .7, alpha=0.5, aes(x=temperature, y=p, color= genotype)) +
  scale_color_grey() +
  new_scale_color() +
  geom_ribbon(data=newdata, aes(x=temperature, y=p, ymin=pmin, ymax=pmax, fill=factor(region, levels = c("temperate","tropical"))), alpha=0.5) +
  geom_line(data=newdata, size = 1, aes(x=temperature, y=p, color = factor(region, levels = c("temperate","tropical")), linetype = factor(region, levels = c("tropical","temperate")))) +
  scale_color_manual(values = c("#009E73", "#D55E00")) +
  scale_fill_manual(values = c("#009E73","#D55E00")) +
  geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.5) +
  xlim(25,45) +
  xlab("Temperature (°C)") + 
  ylab("Proportion hatched") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), legend.title = element_blank(), legend.position = "none") 

#Mixed effects models with lme4. Note here that treating genotype as a random effect produces a singular fit,
#which is bad. Treating genotype as a random effect and testing for the effect
#of region is just illogical.

library(lme4)
# baseline model glm
m0.glm <- glm(hatched/eggs ~ temperature, family = binomial, data = subset(survival_data, region == "tropical" | region == "temperate")) 
m1.glm <- glm(hatched/eggs ~ temperature*region, family = binomial, data = subset(survival_data, region == "tropical" | region == "temperate"))
summary(m1.glm)
Anova(m1.glm)
dose.p(m1.glm)
# base-line mixed-model
m0.glmer <- glmer(hatched/eggs ~ temperature + (1|genotype), data = subset(survival_data, region == "tropical" | region == "temperate"), family = binomial, control=glmerControl(optimizer = "bobyqa")) 
AIC(m0.glm)
#[1] 498.0415
AIC(m0.glmer)
#[1] 324.5691
m1.glmer <- glmer(hatched/eggs ~ temperature*region + (1|genotype), data = subset(survival_data, region == "tropical" | region == "temperate"), family = binomial, control=glmerControl(optimizer = "bobyqa")) 
isSingular(m1.glmer)
#[1] TRUE
AIC(m1.glmer)
#[1] 311.5247
m2.glmer <- glmer(hatched/eggs ~ region + (1|genotype), data = subset(survival_data, region == "tropical" | region == "temperate"), family = binomial, control=glmerControl(optimizer = "bobyqa")) 
isSingular(m2.glmer)
#[1] TRUE

#Plot the CH x VTECK 12 hybrid survival data.
hybrid_data <- subset(survival_data, genotype == "CH_F_X_VTECK12_M" | genotype == "VTECK12_F_X_CH_M")

ggplot(hybrid_data, aes(x = temperature, y = survival, color = genotype, shape = genotype)) +
  geom_point(alpha = 0.4, size = 2.5) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE, size = .7) +
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  scale_color_manual(values = c("#D55E00", "#009E73")) 

glm.all.hybrid <- glm(survival ~ temperature, family=binomial, data = hybrid_data)
glm.hybrid.effect <- glm(survival ~ temperature*genotype, family = binomial, data=hybrid_data)
anova(glm.all.hybrid, glm.hybrid.effect)
summary(glm.hybrid.effect)

glm.CH.vs.VT10 <- glm(survival ~ temperature*genotype, family = binomial, data=subset(survival_data, genotype == "CH" | genotype == "VTECK_10"))
glm.CH.vs.VT10.2 <- glm(survival ~ temperature, family = binomial, data=subset(survival_data, genotype == "CH" | genotype == "VTECK_10"))
anova(glm.CH.vs.VT10, glm.CH.vs.VT10.2, test = "Chisq")
# Analysis of Deviance Table
# 
# Model 1: survival ~ temperature * genotype
# Model 2: survival ~ temperature
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
# 1        74     16.563                       
# 2        76     22.143 -2  -5.5805  0.06141 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

compareGLM(glm.CH.vs.VT10, glm.CH.vs.VT10.2)
# $Models
# Formula                            
# 1 "survival ~ temperature * genotype"
# 2 "survival ~ temperature"           
# 
# $Fit.criteria
# Rank Df.res   AIC  AICc   BIC McFadden Cox.and.Snell Nagelkerke   p.value
# 1    4     74 51.51 52.35 63.30   0.6142        0.5715     0.7637 1.889e-10
# 2    2     76 53.02 53.34 60.09   0.5631        0.5402     0.7218 7.007e-11

summary(glm.CH.vs.VT10)  
# Call:
#   glm(formula = survival ~ temperature * genotype, family = binomial, 
#       data = subset(survival_data, genotype == "CH" | genotype == 
#                       "VTECK_10"))
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.0855  -0.3824  -0.1091   0.2480   0.9191  
# 
# Coefficients:
#                              Estimate Std. Error z value Pr(>|z|)   
# (Intercept)                   10.4297     3.7044   2.815  0.00487 **
# temperature                   -0.2957     0.1036  -2.854  0.00432 **
# genotypeVTECK_10              16.2909     9.5046   1.714  0.08653 . 
# temperature:genotypeVTECK_10  -0.4941     0.2781  -1.776  0.07565 . 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 63.350  on 77  degrees of freedom
# Residual deviance: 16.563  on 74  degrees of freedom
# AIC: 49.514
# 
# Number of Fisher Scoring iterations: 6
qchisq(0.95, df.residual(glm.CH.vs.VT10))

#Based on R Companion of The Handbook of Biological Statistics for Multiple Logistic Regression...
#You'll see that this approach is not the best for embryos because it doesn't allow you to estimate the upper limit,
#and instead forces the upper limit at 1. Because flies lay unfertilized eggs sometimes, particularly among the tropical
#genotypes in the lab, this leads to shallow slopes and inaccurate fits of the logistic regression, to the detrimment
#of not being able to discern the differential effects of temperature.
#Thus, the better approach (by far) is to use the 'drc' package, as above.

#Determine the correct model with a stepwise model comparison of nested models:

#Null model:
glm.null <- glm(survival ~ 1, weights = eggs, family = binomial, data = subset(survival_data, region == "tropical" | region == "temperate"))

#Full saturated model:
glm.full <- glm(hatched/eggs ~ temperature*region + genotype, weights = eggs, family = binomial, data = subset(survival_data, region == "tropical" | region == "temperate"))
summary(glm.full)

glm.SK <- glm(hatched/eggs ~ temperature, weights = eggs, family = binomial, data = subset(survival_data, genotype == "SK"))
summary(glm.SK)
Anova(glm.SK)
dose.p(glm.SK, p=0.5)

#stepwise model comparison of nested models:
step(glm.null, scope = list(upper=glm.full), direction="both", test="Chisq", data = subset(survival_data, region == "tropical" | region == "temperate"))
# Start:  AIC=13153.53
# survival ~ 1
# 
# Df Deviance     AIC    LRT Pr(>Chi)    
# + temperature  1   4382.4  6031.3 7124.2  < 2e-16 ***
# + genotype    23  11085.1 12778.1  421.5  < 2e-16 ***
# + region       1  11502.8 13151.8    3.7  0.05283 .  
# <none>            11506.6 13153.5                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Step:  AIC=6031.31
# survival ~ temperature
# 
# Df Deviance     AIC    LRT  Pr(>Chi)    
# + genotype    23   4073.5  5768.4  308.9 < 2.2e-16 ***
# + region       1   4324.3  5975.3   58.0 2.575e-14 ***
# <none>             4382.4  6031.3                     
# - temperature  1  11506.6 13153.5 7124.2 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Step:  AIC=5768.45
# survival ~ temperature + genotype
# 
# Df Deviance     AIC    LRT  Pr(>Chi)    
# + temperature:genotype 23   3218.4  4959.4  855.1 < 2.2e-16 ***
#   <none>                      4073.5  5768.4                     
# - genotype             23   4382.4  6031.3  308.9 < 2.2e-16 ***
#   - temperature           1  11085.1 12778.1 7011.6 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Step:  AIC=4959.36
# survival ~ temperature + genotype + temperature:genotype
# 
# Df Deviance    AIC    LRT  Pr(>Chi)    
# <none>                      3218.4 4959.4                     
# - temperature:genotype 23   4073.5 5768.4 855.09 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Call:  glm(formula = survival ~ temperature + genotype + temperature:genotype, 
#            family = binomial, data = subset(survival_data, region == 
#                                               "tropical" | region == "temperate"), weights = eggs)
# 
# Coefficients:
#   (Intercept)                   temperature                genotypeBEA_21                genotypeBEA_32                genotypeBEA_36                 genotypeBEA_5  
# 14.69890                      -0.43906                       7.88342                       5.51008                      -2.96303                      24.29807  
# genotypeCH                    genotypeGH                    genotypeGU                    genotypeMU                genotypeRFM_16                genotypeRFM_19  
# -4.26918                      -8.24141                       8.30981                      -3.08606                       5.88625                       3.29838  
# genotypeRFM_34                 genotypeRFM_4                genotypeRFM_48                 genotypeRFM_6                    genotypeSK              genotypeVTECK_10  
# 13.75218                      10.35744                       9.82255                      12.49415                      -4.29399                      12.19823  
# genotypeVTECK_12              genotypeVTECK_14               genotypeVTECK_2               genotypeVTECK_4               genotypeVTECK_5               genotypeVTECK_8  
# 24.61347                      11.00214                      13.06056                      44.51360                       3.69960                      14.65218  
# genotypeVTECK_9    temperature:genotypeBEA_21    temperature:genotypeBEA_32    temperature:genotypeBEA_36     temperature:genotypeBEA_5        temperature:genotypeCH  
# 14.82766                      -0.21805                      -0.14477                       0.08737                      -0.68349                       0.14340  
# temperature:genotypeGH        temperature:genotypeGU        temperature:genotypeMU    temperature:genotypeRFM_16    temperature:genotypeRFM_19    temperature:genotypeRFM_34  
# 0.24077                      -0.21561                       0.11349                      -0.17551                      -0.09177                      -0.37903  
# temperature:genotypeRFM_4    temperature:genotypeRFM_48     temperature:genotypeRFM_6        temperature:genotypeSK  temperature:genotypeVTECK_10  temperature:genotypeVTECK_12  
# -0.30131                      -0.26262                      -0.35663                       0.13932                      -0.35869                      -0.71718  
# temperature:genotypeVTECK_14   temperature:genotypeVTECK_2   temperature:genotypeVTECK_4   temperature:genotypeVTECK_5   temperature:genotypeVTECK_8   temperature:genotypeVTECK_9  
# -0.29799                      -0.35170                      -1.25794                      -0.10689                      -0.42207                      -0.40190  
# 
# Degrees of Freedom: 817 Total (i.e. Null);  770 Residual
# Null Deviance:	    11510 
# Residual Deviance: 3218 	AIC: 4959

#Stepwise model comparison indicates that including region did not improve model fits. 
#So the best model is one that includes both genotpye and temperature and their interactions. 
glm.final <- glm(survival ~ temperature*genotype, weights = eggs, family = binomial, data = subset(survival_data, region == "tropical" | region == "temperate"))

nagelkerke(glm.final) #This computes a Pseudo R-squared of lm, gls, lme, lmer, glmer, glm, negbin, zeroinfl, nls, clm, clmm, or vglm objects.
# $Pseudo.R.squared.for.model.vs.null
# Pseudo.R.squared
# McFadden                             0.630206
# Cox and Snell (ML)                   0.999960
# Nagelkerke (Cragg and Uhler)         0.999960
# 
# $Likelihood.ratio.test
# Df.diff LogLik.diff  Chisq p.value
# -47     -4144.1 8288.2       0

summary(glm.final)
# Call:
#   glm(formula = survival ~ temperature * genotype, family = binomial, 
#       data = subset(survival_data, region == "tropical" | region == 
#                       "temperate"), weights = eggs)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -5.5313  -1.6036  -0.4011   1.1054   9.1790  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                  14.69890    1.20235  12.225  < 2e-16 ***
#   temperature                  -0.43906    0.03560 -12.333  < 2e-16 ***
#   genotypeBEA_21                7.88342    2.05481   3.837 0.000125 ***
#   genotypeBEA_32                5.51008    2.43571   2.262 0.023685 *  
#   genotypeBEA_36               -2.96303    1.36777  -2.166 0.030286 *  
#   genotypeBEA_5                24.29807    3.79818   6.397 1.58e-10 ***
#   genotypeCH                   -4.26918    1.46002  -2.924 0.003455 ** 
#   genotypeGH                   -8.24141    1.36194  -6.051 1.44e-09 ***
#   genotypeGU                    8.30981    2.75068   3.021 0.002519 ** 
#   genotypeMU                   -3.08606    1.51592  -2.036 0.041774 *  
#   genotypeRFM_16                5.88625    2.07187   2.841 0.004497 ** 
#   genotypeRFM_19                3.29838    1.75816   1.876 0.060650 .  
#   genotypeRFM_34               13.75218    3.21742   4.274 1.92e-05 ***
#   genotypeRFM_4                10.35744    2.14629   4.826 1.39e-06 ***
#   genotypeRFM_48                9.82255    2.57113   3.820 0.000133 ***
#   genotypeRFM_6                12.49415    2.23425   5.592 2.24e-08 ***
#   genotypeSK                   -4.29399    1.38853  -3.092 0.001985 ** 
#   genotypeVTECK_10             12.19823    2.41306   5.055 4.30e-07 ***
#   genotypeVTECK_12             24.61347    3.97924   6.185 6.19e-10 ***
#   genotypeVTECK_14             11.00214    2.23381   4.925 8.42e-07 ***
#   genotypeVTECK_2              13.06056    2.55334   5.115 3.14e-07 ***
#   genotypeVTECK_4              44.51360    5.87169   7.581 3.43e-14 ***
#   genotypeVTECK_5               3.69960    1.81176   2.042 0.041153 *  
#   genotypeVTECK_8              14.65218    3.05099   4.802 1.57e-06 ***
#   genotypeVTECK_9              14.82766    2.55659   5.800 6.64e-09 ***
#   temperature:genotypeBEA_21   -0.21805    0.06024  -3.620 0.000295 ***
#   temperature:genotypeBEA_32   -0.14477    0.07045  -2.055 0.039870 *  
#   temperature:genotypeBEA_36    0.08737    0.04050   2.157 0.030985 *  
#   temperature:genotypeBEA_5    -0.68349    0.10960  -6.236 4.48e-10 ***
#   temperature:genotypeCH        0.14340    0.04247   3.376 0.000734 ***
#   temperature:genotypeGH        0.24077    0.04000   6.020 1.75e-09 ***
#   temperature:genotypeGU       -0.21561    0.07766  -2.776 0.005498 ** 
#   temperature:genotypeMU        0.11349    0.04396   2.582 0.009831 ** 
#   temperature:genotypeRFM_16   -0.17551    0.06092  -2.881 0.003967 ** 
#   temperature:genotypeRFM_19   -0.09177    0.05187  -1.769 0.076843 .  
#   temperature:genotypeRFM_34   -0.37903    0.09318  -4.068 4.75e-05 ***
#   temperature:genotypeRFM_4    -0.30131    0.06294  -4.787 1.69e-06 ***
#   temperature:genotypeRFM_48   -0.26262    0.07575  -3.467 0.000527 ***
#   temperature:genotypeRFM_6    -0.35663    0.06527  -5.464 4.65e-08 ***
#   temperature:genotypeSK        0.13932    0.04082   3.413 0.000643 ***
#   temperature:genotypeVTECK_10 -0.35869    0.07140  -5.024 5.07e-07 ***
#   temperature:genotypeVTECK_12 -0.71718    0.11673  -6.144 8.06e-10 ***
#   temperature:genotypeVTECK_14 -0.29799    0.06468  -4.607 4.08e-06 ***
#   temperature:genotypeVTECK_2  -0.35170    0.07370  -4.772 1.82e-06 ***
#   temperature:genotypeVTECK_4  -1.25794    0.16859  -7.461 8.56e-14 ***
#   temperature:genotypeVTECK_5  -0.10689    0.05353  -1.997 0.045852 *  
#   temperature:genotypeVTECK_8  -0.42207    0.08911  -4.737 2.17e-06 ***
#   temperature:genotypeVTECK_9  -0.40190    0.07353  -5.466 4.60e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 11506.6  on 817  degrees of freedom
# Residual deviance:  3218.4  on 770  degrees of freedom
# AIC: 4959.4
# 
# Number of Fisher Scoring iterations: 5

Anova(glm.final, type="II", test="Wald") #ANOVA for the individual terms in the model.
# Analysis of Deviance Table (Type II tests)
# 
# Response: survival
#                      Df   Chisq Pr(>Chisq)    
# temperature           1 3018.53  < 2.2e-16 ***
# genotype             23  219.50  < 2.2e-16 ***
# temperature:genotype 23  766.35  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(glm.final, glm.null, test="Chisq") #overall p-value for the model, via Chi-square test.
# Analysis of Deviance Table
# 
# Model 1: survival ~ temperature * genotype
# Model 2: survival ~ 1
#   Resid. Df Resid. Dev  Df Deviance  Pr(>Chi)    
# 1       770     3218.4                           
# 2       817    11506.6 -47  -8288.2 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lrtest(glm.final) #overall p-value for the model, via likelihood ratio test.
# Likelihood ratio test
# 
# Model 1: survival ~ temperature * genotype
# Model 2: survival ~ 1
#   #Df  LogLik  Df  Chisq Pr(>Chisq)    
# 1  48 -2431.7                          
# 2   1 -6575.8 -47 8288.2  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Now let's plot the results:
ggplot(data=subset(survival_data, region == "tropical" | region == "temperate"), aes(x = temperature, y = survival, color = genotype)) +
  geom_point(alpha = 0.4, size = 2.5) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, size = .7, fullrange=TRUE) +
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) 
  
#Below is plotting a Kaplan-Meier curve:
library(survival)
survival.plot <- survfit(Surv(temperature, survival) ~ region, data = subset(survival_data, region == "tropical" | region == "temperate"))
autoplot(survival.plot)

#playing around with count-based methods...

#add column of hatched egg counts to the data frame:
hatched <- survival_data$survival*20
survival_data <- cbind(survival_data, hatched)
survival_data$hatched <- as.integer(survival_data$hatched)  

#add column of not hatched egg counts to the data frame:
dead <- 20 - survival_data$survival*20
survival_data <- cbind(survival_data, dead)
survival_data$dead <- as.integer(survival_data$dead) 

#modeling alive eggs...
glm.CH.vs.VT12.poisson <- glm(hatched ~ temperature*genotype, family = poisson, data=subset(survival_data, genotype == "CH" | genotype == "VTECK_12"))
summary(glm.CH.vs.VT12.poisson)


library(MASS)
nb.glm.1 <- glm.nb(dead ~ temperature*genotype, data=subset(survival_data, genotype == "CH" | genotype == "VTECK_10"))
nb.glm.2 <- glm.nb(dead ~ temperature, data=subset(survival_data, genotype == "CH" | genotype == "VTECK_10"))
anova(nb.glm.1,nb.glm.2)
# Likelihood ratio tests of Negative Binomial Models
# 
# Response: dead
# Model    theta Resid. df    2 x log-lik.   Test    df LR stat.      Pr(Chi)
# 1            temperature 3.678593        76       -452.3001                                   
# 2 temperature * genotype 5.643835        74       -432.2517 1 vs 2     2 20.04849 4.431239e-05

summary(nb.glm.1)
# Call:
#   glm.nb(formula = dead ~ temperature * genotype, data = subset(survival_data, 
#                                                                 genotype == "CH" | genotype == "VTECK_10"), init.theta = 5.643834604, 
#          link = log)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.9110  -0.9678  -0.1319   0.5498   2.2079  
# 
# Coefficients:
#                                Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)                  -2.16221    0.55488  -3.897 9.75e-05 ***
#   temperature                   0.12136    0.01517   7.998 1.27e-15 ***
#   genotypeVTECK_10             -5.10298    1.13112  -4.511 6.44e-06 ***
#   temperature:genotypeVTECK_10  0.14725    0.03135   4.697 2.63e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(5.6438) family taken to be 1)
# 
# Null deviance: 272.464  on 77  degrees of freedom
# Residual deviance:  97.669  on 74  degrees of freedom
# AIC: 442.25
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  5.64 
# Std. Err.:  1.87 
# 
# 2 x log-likelihood:  -432.252
qchisq(0.95, df.residual(nb.glm.1)) #this provides the five-percent critical value for a chi-square with 38 degrees of freedom. 
#[1] 95.08147
deviance(nb.glm.1) #Note that the deviance of the negative binomial model is more than the critical value, so this indicates not a good fit.
# [1] 97.66883
df.residual(nb.glm.1)
# [1] 74
pchisq(deviance(nb.glm.1), df.residual(nb.glm.1), lower.tail = FALSE)
# [1] 0.03412109

ggplot(data=subset(survival_data, genotype == "CH" | genotype == "VTECK_10"), aes(x = temperature, y = hatched, color = genotype, shape = genotype)) +
  geom_point(alpha = 0.4, size = 2.5) +
  geom_smooth(method = "glm.nb", se = TRUE, size = .7) +
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  scale_color_manual(values = c("#D55E00", "#009E73")) 

ggplot(data=subset(survival_data, genotype == "CH" | genotype == "VTECK_10"), aes(x = temperature, y = hatched, color = genotype, shape = genotype)) +
  geom_point(alpha = 0.4, size = 2.5) +
  geom_smooth(method = "glm", method.args = list(family = "poisson"), se = TRUE, size = .7) +
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  scale_color_manual(values = c("#D55E00", "#009E73")) 

ggplot(data=subset(survival_data, genotype == "CH" | genotype == "VTECK_2"), aes(x = temperature, y = hatched, color = genotype, shape = genotype)) +
  geom_point(alpha = 0.4, size = 2.5) +
  geom_smooth(method = "glm", method.args = list(family = "quasipoisson"), se = TRUE, size = .7) +
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  scale_color_manual(values = c("#D55E00", "#009E73")) 

ggplot(data=subset(survival_data, genotype == "CH" | genotype == "VTECK_10"), aes(x = temperature, y = survival, color = genotype, shape = genotype)) +
  geom_point(alpha = 0.4, size = 2.5) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE, size = .7) +
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  scale_color_manual(values = c("#D55E00", "#009E73")) 

#modeling dead eggs...
ggplot(data=subset(survival_data, genotype == "CH" | genotype == "VTECK_10"), aes(x = temperature, y = dead, color = genotype, shape = genotype)) +
  geom_point(alpha = 0.4, size = 2.5) +
  geom_smooth(method = "glm", method.args = list(family = "poisson"), se = TRUE, size = .7) +
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  scale_color_manual(values = c("#D55E00", "#009E73")) 

ggplot(data=subset(survival_data, genotype == "CH" | genotype == "VTECK_2"), aes(x = temperature, y = dead, color = genotype, shape = genotype)) +
  geom_point(alpha = 0.4, size = 2.5) +
  geom_smooth(method = "glm", method.args = list(family = "quasipoisson"), se = TRUE, size = .7) +
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  scale_color_manual(values = c("#D55E00", "#009E73")) 

ggplot(data=subset(survival_data, genotype == "CH" | genotype == "VTECK_10"), aes(x = temperature, y = dead, color = genotype, shape = genotype)) +
  geom_point(alpha = 0.4, size = 2.5) +
  geom_smooth(method = "glm.nb", se = TRUE, size = .7, fullrange=TRUE) +
  #xlim(0,45) +
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  scale_color_manual(values = c("#D55E00", "#009E73")) 



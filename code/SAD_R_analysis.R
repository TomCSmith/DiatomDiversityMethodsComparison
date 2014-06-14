library(car)
## Doing some model fitting on SAD parameters

# Read in Holmes data...I think this is Windows compatible
lake_data <- read.csv(file.path("..", "data", "formatted",
                "Lake_SAD_Analysis.csv"), header=T)

# Correlation matrix
param_data = lake_data[c("sigma", "fishers_alpha", "alpha", "theta")]

# Look at the correlation matrix...are these parameters telling us the same
# things? These parameters are the results of fitting lognormal, logseries,
# and discrete gamma distributions to the SADs
cor_mat = cor(param_data)

# They are pretty tightly correlated
#                    sigma fishers_alpha      alpha      theta
# sigma          1.0000000    -0.9076855 -0.7598571  0.8998177
# fishers_alpha -0.9076855     1.0000000  0.4961870 -0.8244380
# alpha         -0.7598571     0.4961870  1.0000000 -0.8614348
# theta          0.8998177    -0.8244380 -0.8614348  1.0000000

# Looks like fisher's alpha nd the lognormal sigma are very tightly correlated

# Look at the correlation plots...nothing looks shockingly non-linear
pairs(param_data)

# Now if I fit a model with S, N, Shannons Diversity how much variance would be
# left unexplained?


# Fitting a model fisher's alpha as the response and S, N, and Shannons as
# predictors
fit1 = lm(log(sigma) ~ log(S), data=lake_data)

# Look at the residuals and a qqplot...they are looking pretty good
par(mfrow=c(1, 3))
plot(fit1$fitted.values, fit1$resid)
abline(h=0)
qqnorm(fit1$residuals)
qqline(fit1$residuals)

plot(log(lake_data$S), log(lake_data$sigma), xlab="Log S", ylab="Log Sigma")

summary(fit1)

# Call:
# lm(formula = log(sigma) ~ log(S), data = lake_data)

# Residuals:
#       Min        1Q    Median        3Q       Max
# -0.114034 -0.040810  0.003423  0.038664  0.084968

# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)  1.69212    0.11871   14.25 2.32e-14 ***
# log(S)      -0.35164    0.03054  -11.51 3.93e-12 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.05248 on 28 degrees of freedom
# Multiple R-squared:  0.8256,    Adjusted R-squared:  0.8194
# F-statistic: 132.6 on 1 and 28 DF,  p-value: 3.926e-12

# Looks like about 80% of the variance in sigma can be explained by Speices
# richness...let's looks at how other environmental variables are affecting it

fit2 = lm(log(sigma) ~ S + pHLab + Cl + HCO3 + SO4 + Ca + Mg + Na + K + SedDryWt,
        data = lake_data)

# Backward selection to figure out what variables I want to keep
step(fit2, direction="backward")
# Start:  AIC=-182.46
# log(sigma) ~ S + pHLab + Cl + HCO3 + SO4 + Ca + Mg + Na + K +
#     SedDryWt

#            Df Sum of Sq      RSS     AIC
# - SedDryWt  1  0.000005 0.032907 -184.46
# - Mg        1  0.000253 0.033156 -184.23
# - Na        1  0.000949 0.033852 -183.61
# - K         1  0.000992 0.033895 -183.57
# - Ca        1  0.001035 0.033937 -183.53
# - HCO3      1  0.001075 0.033978 -183.50
# - SO4       1  0.001128 0.034031 -183.45
# <none>                  0.032903 -182.46
# - Cl        1  0.003132 0.036035 -181.73
# - pHLab     1  0.006220 0.039123 -179.27
# - S         1  0.135022 0.167924 -135.56

# Step:  AIC=-184.46
# log(sigma) ~ S + pHLab + Cl + HCO3 + SO4 + Ca + Mg + Na + K

#         Df Sum of Sq      RSS     AIC
# - Mg     1  0.000268 0.033175 -186.22
# - Na     1  0.001071 0.033979 -185.50
# - K      1  0.001257 0.034164 -185.33
# - Ca     1  0.001355 0.034263 -185.25
# - HCO3   1  0.001412 0.034319 -185.20
# - SO4    1  0.001625 0.034532 -185.01
# <none>               0.032907 -184.46
# - Cl     1  0.003384 0.036291 -183.52
# - pHLab  1  0.006228 0.039135 -181.26
# - S      1  0.221916 0.254823 -125.05

# Step:  AIC=-186.21
# log(sigma) ~ S + pHLab + Cl + HCO3 + SO4 + Ca + Na + K

#         Df Sum of Sq      RSS     AIC
# - K      1  0.001002 0.034177 -187.32
# - Na     1  0.001204 0.034379 -187.15
# - Ca     1  0.002231 0.035406 -186.26
# - SO4    1  0.002281 0.035456 -186.22
# <none>               0.033175 -186.22
# - HCO3   1  0.002378 0.035554 -186.14
# - Cl     1  0.003688 0.036863 -185.05
# - pHLab  1  0.006080 0.039256 -183.17
# - S      1  0.245101 0.278277 -124.41

# Step:  AIC=-187.32
# log(sigma) ~ S + pHLab + Cl + HCO3 + SO4 + Ca + Na

#         Df Sum of Sq     RSS     AIC
# - Na     1  0.001466 0.03564 -188.06
# - SO4    1  0.001490 0.03567 -188.04
# - Ca     1  0.001579 0.03576 -187.97
# - HCO3   1  0.001626 0.03580 -187.93
# <none>               0.03418 -187.32
# - Cl     1  0.003611 0.03779 -186.31
# - pHLab  1  0.005097 0.03927 -185.15
# - S      1  0.291606 0.32578 -121.68

# Step:  AIC=-188.06
# log(sigma) ~ S + pHLab + Cl + HCO3 + SO4 + Ca

#         Df Sum of Sq     RSS     AIC
# - SO4    1  0.000256 0.03590 -189.85
# - Ca     1  0.000310 0.03595 -189.80
# - HCO3   1  0.000325 0.03597 -189.79
# - Cl     1  0.002187 0.03783 -188.28
# <none>               0.03564 -188.06
# - pHLab  1  0.010072 0.04572 -182.60
# - S      1  0.305786 0.34143 -122.27

# Step:  AIC=-189.85
# log(sigma) ~ S + pHLab + Cl + HCO3 + Ca

#         Df Sum of Sq     RSS     AIC
# - Ca     1  0.000054 0.03595 -191.80
# - HCO3   1  0.000069 0.03597 -191.79
# <none>               0.03590 -189.85
# - Cl     1  0.005600 0.04150 -187.50
# - pHLab  1  0.010943 0.04684 -183.87
# - S      1  0.308748 0.34465 -123.99

# Step:  AIC=-191.8
# log(sigma) ~ S + pHLab + Cl + HCO3

#         Df Sum of Sq     RSS     AIC
# - HCO3   1   0.00002 0.03597 -193.78
# <none>               0.03595 -191.80
# - Cl     1   0.00573 0.04169 -189.36
# - pHLab  1   0.01121 0.04717 -185.66
# - S      1   0.36899 0.40495 -121.16

# Step:  AIC=-193.79
# log(sigma) ~ S + pHLab + Cl

#         Df Sum of Sq     RSS     AIC
# <none>               0.03597 -193.78
# - Cl     1   0.00577 0.04175 -191.32
# - pHLab  1   0.02710 0.06308 -178.94
# - S      1   0.37150 0.40748 -122.97

# Call:
# lm(formula = log(sigma) ~ S + pHLab + Cl, data = lake_data)

# Coefficients:
# (Intercept)            S        pHLab           Cl
#    0.455180    -0.007535     0.040107    -0.010671

# Make a new model with only the signficant predictors

fit3 = lm(log(sigma) ~ S + Cl + pHLab, data=lake_data)

# ReTest the assumptions...a bit of non-normality in the residuals, but not
# but not horrible

par(mfrow=c(1, 2))
plot(fit3$fitted.values, fit3$resid)
abline(h=0)
qqnorm(fit3$residuals)
qqline(fit3$residuals)

summary(fit3)

# Call:
# lm(formula = log(sigma) ~ S + Cl + pHLab, data = lake_data)

# Residuals:
#       Min        1Q    Median        3Q       Max
# -0.090370 -0.018325 -0.002035  0.032640  0.046717

# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.4551800  0.0643991   7.068 1.66e-07 ***
# S           -0.0075352  0.0004599 -16.386 3.20e-15 ***
# Cl          -0.0106708  0.0052231  -2.043 0.051312 .
# pHLab        0.0401066  0.0090617   4.426 0.000153 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Sequential sum of squares
Anova(fit3, type="III")

# Anova Table (Type III tests)

# Response: log(sigma)
#              Sum Sq Df  F value    Pr(>F)
# (Intercept) 0.06912  1  49.9582 1.662e-07 ***
# S           0.37150  1 268.5074 3.201e-15 ***
# Cl          0.00577  1   4.1738 0.0513121 .
# pHLab       0.02710  1  19.5888 0.0001529 ***
# Residuals   0.03597 26
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# After accounting for the variability explained explained by Species Richness,
# pH of the lake explains about 6% of the variation in the sigma parameter of
# of the lognormal.  So, statistically the SAD is capturing a bit of information
# that species richness is not, but is that biologically signficant?

# Next step:  Let's look at the mean, variance, skew, and kurtosis of each of
# SADs
library(ggplot2)
library(GGally)
library(CCA)
library(CCP)

# Get example1 data
mm <- read.csv("https://stats.idre.ucla.edu/stat/data/mmreg.csv")
colnames(mm) <- c("Control", "Concept", "Motivation", "Read", "Write", "Math", 
                  "Science", "Sex")
summary(mm)

# CCA
xtabs(~Sex, data=mm)
psych = mm[, 1:3] # control, concept, motivation
acad = mm[, 4:8] # read, write, ...

# View intervariable correlation for each table
ggpairs(psych) 
ggpairs(acad)

# Correlations within AND betw variables of two tables
# Simple descriptive matrix correlation fn of CCA package
matcor(psych, acad)
# -> output: corr within Psych, Acad, and betw Psych-Acad

# Canonical correlation analysis
cc1 = cc(psych, acad)

# display canonical correlations
cc1$cor

# raw canonical coefficients
cc1[3:4]
# Read+1 -> -0.04462 change in Canonical covariate in set 2
# canonical covariate: linear comb of orig variables

# Compute loadings of variables on the canonical dimensions (variates)
cc2 = comput(psych, acad, cc1)

# display canonical loadings
cc2[3:6]

# tests of canonical dimensions
rho <- cc1$cor
# Define number of observations, number of variables in first set, and number of variables in the second set.
n <- dim(psych)[1]
p <- length(psych)
q <- length(acad)

# Calculate p-values using the F-approximations of different test statistics:
p.asym(rho, n, p, q, tstat = "Wilks")
# -> only canonical variates #1 and #2 are significant

# Wilks' Lambda, using F-approximation (Rao's F):
#             stat     approx    df1    df2     p.value
# 1 ~ 3 sig? 0.7543611 11.715733  15 1634.653 0.000000000
# 2 ~ 3 sig? 0.9614300  2.944459   8 1186.000 0.002905057
# 3 ~ 3 sig? 0.9891858  2.164612   3  594.000 0.091092180

# standardized psych canon coeff diag matrix of psych SD
s1 = diag(sqrt(diag(cov(psych))))
s1 %*% cc1$xcoef

# standardized acad canon coeff diag matrix of acad SD
s2 = diag(sqrt(diag(cov(acad))))
s2 %*% cc1$ycoef
# Read SD+1 -> canonical variate #1 in set2 SD -0.45 [dec]


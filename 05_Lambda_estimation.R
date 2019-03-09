#############################
#     LAMBDA ESTIMATION     #
#############################

# Using the library “gap”
install.packages("gap")
library(gap)
r=gcontrol2(results_short$P)
print (r$lambda)


# Using the library “GenAbel”
#install.packages("GenABEL")
library(GenABEL)
r_GenAbel <- estlambda(results_short$P, plot = FALSE, proportion = 1, method = "regression")
print(r_GenAbel$estimate)

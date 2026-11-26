# Set the base directory using your favorite method
# setwd("...")

##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (BEGINNING)
##################################################################################################
#	INPUT. Original datafiles of the case study, placed in the data folder.

#	OUTPUT. Unfitted models, i.e., the list of Hmsc model(s) that have been defined
# but not fitted yet, stored in the file "models/unfitted_models.RData".
##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (END)
##################################################################################################


##################################################################################################
# MAKE THE SCRIPT REPRODUCIBLE (BEGINNING)
##################################################################################################
set.seed(1)
##################################################################################################
## MAKE THE SCRIPT REPRODUCIBLE (END)
##################################################################################################


##################################################################################################
# LOAD PACKAGES (BEGINNING)
##################################################################################################
library(Hmsc)
##################################################################################################
# LOAD PACKAGES (END)
##################################################################################################


##################################################################################################
# SET DIRECTORIES (BEGINNING)
##################################################################################################
localDir = "."
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir, "models")
if(!dir.exists(modelDir)) dir.create(modelDir)
##################################################################################################
# SET DIRECTORIES (END)
##################################################################################################


##################################################################################################
# READ AND EXPLORE THE DATA (BEGINNING)
##################################################################################################
# Write here the code needed to read in the data, and explore it by (with View, plot, hist, ...)
# to get and idea of it and ensure that the data are consistent
##################################################################################################
# READ AND EXPLORE THE DATA (END)
##################################################################################################

##################################################################################################
# SET UP THE MODEL (BEGINNING)
##################################################################################################
# Note that many of the components are optional 
# Here instructions are given at generic level, the details will depend on the model

# Organize the community data in the matrix Y

# Organize the environmental data into a dataframe XData
# Define the environmental model through XFormula
# XFormula = ~ ...

# Organize the trait data into a dataframe TrData
# Define the trait model through TrFormula
# TrFormula = ~ ..

# Set up a phylogenetic (or taxonomic tree) as myTree

# Define the studyDesign as a dataframe 
# For example, if you have sampled the same locations over multiple years, you may define
# studyDesign = data.frame(sample = ..., year = ..., location = ...)

# Set up the random effects
# For example, you may define year as an unstructured random effect
# rL.year = HmscRandomLevel(units = levels(studyDesign$year))
# For another example, you may define location as a spatial random effect
# rL.location = HmscRandomLevel(sData = locations.xy)
# Here locations.xy would be a matrix (one row per unique location)
# where row names are the levels of studyDesign$location,
# and the columns are the xy-coordinates
# For another example, you may define the sample = sampling unit = row of matrix Y
# as a random effect, in case you are interested in co-occurrences at that level
# rL.sample = HmscRandomLevel(units = levels(studyDesign$sample))

# Use the Hmsc model constructor to define a model

# m = Hmsc(Y=Y,
#          distr="probit",
#          XData = XData,  XFormula=XFormula,
#          TrData = TrData, TrFormula = TrFormula,
#          phyloTree = myTree,
#          studyDesign = studyDesign, 
#          ranLevels=list(year=rL.year, location = rL.location, sample = rL.sample))

# note that in the random effects the left-hand sides in the list (year, location, sample)
# refer to the columns of the studyDesign

# In this example we assumed the probit distribution as appropriate for presence-absence data

# It is always a good idea to look at the model object, so type m to the console and press enter
# Look at the components of the model by exploring m$...

# You may define multiple models, e.g.
# alternative covariate selections
# alternative random effect choices
# alternative ways to model the data (e.g. lognormal Poisson versus hurdle model)
# models for different subsets of the data

# m.alternative = Hmsc(....)
##################################################################################################
# SET UP THE MODEL (END)
##################################################################################################


##################################################################################################
# In the general case, we could have multiple models. We combine them into a list and given them names.
# COMBINING AND SAVING MODELS (START)
# models = list(m, m.alternative)
# names(models) = c("my.main.model","my.alternative model")
# save(models, file = file.path(modelDir, "unfitted_models.RData"))
##################################################################################################
# COMBINING AND SAVING MODELS (END)
##################################################################################################


##################################################################################################
# TESTING THAT MODELS FIT WITHOUT ERRORS (START)
##################################################################################################
#for(i in 1:length(models)){
#  print(i)
#  sampleMcmc(models[[i]],samples=2)
#}
##################################################################################################
# TESTING THAT MODELS FIT WITHOUT ERRORS (END)
##################################################################################################

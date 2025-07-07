You can perform training AUC-based stepwise selection using the R code shown below.
The TrainAUCStepwise function requires eight mandatory input variables, each described as follows:

totvar: A list of all candidate predictor variables to be included in model development

dat: The dataset

fixvar: A list of variables that must be included in the model (default = '')

excvar: A list of variables from totvar that should be excluded (default = '')

std.time: The time point (in years) at which to develop the survival prediction model

numSeed: The seed number (example: numSeed=1)

imtres: The initial results to start model development (default = NULL)

addres: The directory path to save intermediate and final results

---------------------------------------------------------------------------------
source('Survival_TrainAUC_StepwiseSelection.R')

library(caret)

library(nsROC)

library(survival)

TrainAUCStepwise(totvar,dat,fixvar,excvar,std.time,numSeed,imtres,addres)
---------------------------------------------------------------------------------

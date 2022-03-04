# Lab 4b

######################################################################################
#
# Intro (DON'T INCLUDE IN REPORT)
#
######################################################################################

# Lines that start with a hashtag/pound symbol, like this one, are comment lines.
# Comment lines are ignored by R when it is running the code.

# To run a non-commented line in RStudio, click the "Run" button above
# or press Ctrl+Enter on Windows and Cmd+Enter on macOS.
# Any output will be printed in the Console pane (placed below) and all plots
# will be displayed in the Plots pane (placed bottom right).

# Complete the lab work by following the instructions and running the provided code.
# You may be required to edit some lines before running them to achieve the desired result.

######################################################################################
#
# Import the necessary libraries (DON'T INCLUDE IN REPORT)
#
######################################################################################

# First time you run this file, you will need to install several packages.
# To do that, run code lines 27-31. It may take up a couple of minutes.
# You only need to install packages once, next time you should skip those lines.
# install.packages("devtools")
library("devtools")
install_github("nusoniclab/rem")

# Now run the line below to load the package you have installed.
# You need to load packages every time you run the script or restart R.
library(remsonic)

######################################################################################
#
# Set current directory (DON'T INCLUDE IN REPORT)
#
######################################################################################

# In this step you tell R where to look for your files.
# From the menu, select "Session > Set Working Directory... > To Source File Location".

# Clean other variables and set a seed
rm(list = ls())
start_time <- Sys.time()
set.seed(42)

######################################################################################
#
# Part II: Run the models
#
######################################################################################

# Load Data
evtdata <- read.table('MTS19.tsv',
                   header=TRUE, sep = "\t",
                   stringsAsFactors=FALSE,
                   strip.white=TRUE)[,1:4]

# Extract Role and Team
x <- t(simplify2array(strsplit(sort(unique(evtdata$sender)),'.',fixed=TRUE)))
role <- x[,1]
names(role) <- sort(unique(evtdata$sender))
team <- x[,2]
names(team) <- sort(unique(evtdata$sender))

# Add Role and Team as Covariates
evtdata$sender.role <- role[evtdata$sender]
evtdata$sender.team <- team[evtdata$sender]
evtdata$receiver.role <- role[evtdata$receiver]
evtdata$receiver.team <- team[evtdata$receiver]

# Format Time
evtdata$time <- sapply(as.character(evtdata$eventMilitaryTime),
                       function(x) sum(as.numeric(strsplit(x,':')[[1]])*c(60*60,60,1)))

# Convert sender and receiver to id numbers
id <- 1:20
names(id) <- sort(unique(evtdata$sender))
evtdata$sid <- id[evtdata$sender]
evtdata$rid <- id[evtdata$receiver]

# Extract just the event data
data <- data.frame(sid = evtdata$sid,rid = evtdata$rid,time=evtdata$time)

# Keep only complete cases from one session
data <- data[complete.cases(data),][1:1132,]

# Set the first moment as zero and readjust the timeline
data$time <- data$time - min(data$time)

# Model 1 -----------------------------------------------------------------

# Calculate the first set of statistics: the intercept, the effects for repetition (RSndSnd), 
# and reciprocity (RRecSnd)
stats.intercept <- Constant(data)
stats.rrecsnd <- RRecSnd(data)
stats.rsndsnd <- RSndSnd(data)
stats1 <- combine.stats(
  '[Intercept]' = stats.intercept,
  'RRecSnd' = stats.rrecsnd,
  'RSndSnd' = stats.rsndsnd
)

# Run the first model and check the results
model1 <- FitEventNetworkCore(data,stats1,ordinal=FALSE)
summary(model1)

# Model 2 -----------------------------------------------------------------

# Continue adding the second term: the Normalized Total Degree Received (NTDRec)
stats.ntdegrec = NTDRec(data)
stats2 <- combine.stats(
  '[Intercept]' <- stats.intercept,
  'RRecSnd' = stats.rrecsnd,
  'RSndSnd' = stats.rsndsnd,
  'NTDegRec' = stats.ntdegrec
)

# Run the second model and check the results
model2 <- FitEventNetworkCore(data,stats2,ordinal=FALSE)
summary(model2)

# Model 3 -----------------------------------------------------------------

# Continue adding the terms that capture the tendencies to speak to other members of your team (SameTeam) 
# and with your same role (SameRole). 
stats.sameconstgroup.team <- SameConstGroup(data,team)
stats.sameconstgroup.role <- SameConstGroup(data,role)
stats3 <- combine.stats(
  '[Intercept]' = stats.intercept,
  'RRecSnd' = stats.rrecsnd,
  'RSndSnd' = stats.rsndsnd,
  'NTDegRec' = stats.ntdegrec,
  'SameTeam' = stats.sameconstgroup.team,
  'SameRole' = stats.sameconstgroup.role
)

# Run the third model and check the results
model3 <- FitEventNetworkCore(data,stats3,ordinal=FALSE)
summary(model3)

# Model 4 -----------------------------------------------------------------

# Add the particularly the AB-BA shift (a tendency for B to call A, given that A has just called B). 
stats.psabba <- PSAB.BA(data)
stats4 <- combine.stats(
  '[Intercept]' = stats.intercept,
  'RRecSnd' = stats.rrecsnd,
  'RSndSnd' = stats.rsndsnd,
  'NTDegRec' = stats.ntdegrec,
  'SameTeam' = stats.sameconstgroup.team,
  'SameRole' = stats.sameconstgroup.role,
  'PSAB-BA' = stats.psabba
)

# Run the fourth model and check the results
model4 <- FitEventNetworkCore(data,stats4,ordinal=FALSE)
summary(model4)

# Model 5 -----------------------------------------------------------------

# Add the particularly the AB-BY shift (it captures spread of information in a chain)
stats.psabby <- PSAB.BY(data)
stats5 <- combine.stats(
  '[Intercept]' = stats.intercept,
  'RRecSnd' = stats.rrecsnd,
  'RSndSnd' = stats.rsndsnd,
  'NTDegRec' = stats.ntdegrec,
  'SameTeam' = stats.sameconstgroup.team,
  'SameRole' = stats.sameconstgroup.role,
  'PSAB-BA' = stats.psabba,
  'PSAB-BY' = stats.psabby
)

# Run the fifth model and check the results
model5 <- FitEventNetworkCore(data,stats5,ordinal=FALSE)
summary(model5)

# Finish time
end_time <- Sys.time()
end_time - start_time


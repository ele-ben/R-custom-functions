# Useful functions

# From logOR to probability --------------------------------------------

logOdd2Prob <- function (logOdd) {
  oddRatio <- 2.72^logOdd
  probOf1 <- oddRatio/(oddRatio+1)
  df <- data.frame("oddRatio" = oddRatio, "Probability of 1" = probOf1)
  df}

# Create OR column in exporting tables  ---------------------------------------------------

OddRatio_tab <- function (lmerTestmodel){
  estim <- as.data.frame(coef(summary(lmerTestmodel)))
  for (i in 1:nrow(estim)){
    if (as.numeric(estim[i, "Pr(>|z|)"]) < 0.001){estim[i, "Pr(>|z|)"] = "< .001"}
    else if (as.numeric(estim[i, "Pr(>|z|)"]) < 0.01){estim[i, "Pr(>|z|)"] = "< .01"}
    else if (as.numeric(estim[i, "Pr(>|z|)"]) < 0.05){estim[i, "Pr(>|z|)"] = "< .05"}
    else {estim[i, "Pr(>|z|)"] = as.character(round(as.numeric(estim[i, "Pr(>|z|)"]), 3))}}
  estim$OddRatio <- exp(estim$Estimate)
  estim <- estim[,c(1:3,5,4)]
  for (i in c(1:4)) {estim[,i] <- round(estim[,i], 3)}
  estim}

write_pvalues <- function (lmerTestmodel, transformation){
  estim <- as.data.frame(coef(summary(lmerTestmodel)))
  if (transformation == "exp"){
    estim$exp_Estimate <- exp(estim$Estimate)
    estim <- estim[,c(1,2,6,4,5)]
    for (i in c(1:4)) {estim[,i] <- round(estim[,i], 3)}}
  else {estim <- estim[,c(1,2,4,5)]
  for (i in c(1:3)) {estim[,i] <- round(estim[,i], 3)}}
  for (i in 1:nrow(estim)){
    if (as.numeric(estim[i, "Pr(>|t|)"]) < 0.001){estim[i, "Pr(>|t|)"] = "< .001"}
    else if (as.numeric(estim[i, "Pr(>|t|)"]) < 0.01){estim[i, "Pr(>|t|)"] = "< .01"}
    else if (as.numeric(estim[i, "Pr(>|t|)"]) < 0.05){estim[i, "Pr(>|t|)"] = "< .05"}
    else {estim[i, "Pr(>|t|)"] = as.character(round(as.numeric(estim[i, "Pr(>|t|)"]), 3))}}
  estim}

# Exporting afex::anova table ----------------------------------------------------

export_aovNice <- function (tab, newNames = ""){
  # Remove comma in df column
  tab$df <- gsub(",", "", tab$df)
  # Remove stars in F column
  tab[,"F"] <- gsub("\\*|\\+| ", "", tab[, "F"])
  # More readable pes name
  names(tab)[names(tab) == "pes"] <- "Partial Eta Sq." 
  # Change them with user-input newNames (order matters!)
  if (all(newNames != "")){ # only if an input is given
    #print(all(newNames != ""))
    varNames <- grep("^[a-zA-Z0-9\\.]+(_)*[a-zA-Z]+$", tab$Effect, value = T) # detect the main effects
    #print(varNames)
    for (i in 1:length(varNames)){ # loop over the old names ans change them with newNames
      tab[,"Effect"] <- gsub(varNames[i], newNames[i],tab[,"Effect"])
      }
    }
  # Change : with x
  tab[,"Effect"] <- gsub(":", " x ", tab[, "Effect"]) 
  # Return the table
  tab
}

# Compute ICC -------------------------------------------------------------------------
ICClogit <- function(interceptStDev){ICC = (interceptStDev)^2/((interceptStDev)^2+ pi^2/3)
ICC}


# Sequence Relation Column ---------------------------------------------------------------

# NOtes:
# I've changed blockLength with maxTrialNum to accomodate blocks with 0-95 and with 1-96
# But the problem was probably in rownames(d) that must be reset to null before...
# Nope, the problem was that moving the indexes forward by lag in the checks, these yield indexes greater
# than the maximun number of d rows if the error was in the last trial last block!

sequence_relation <- function(d, variab, blockLength, suffix = "R", type = "other", Lag = 1, values = "0,1"){
  if (min(d$trialNum) == 0){maxTrialNum = blockLength - 1
  } else {maxTrialNum = blockLength}
  for (var in variab){
    varName <- paste(var, suffix, sep="_")
    #print(d[1:10,var])
    d[[varName]] <- 99
    for (j in unique(d$pp)){
      for (jj in unique(d$blockNum)){
        strt = min(which(d$pp == j & d$blockNum == jj))+Lag #no for 1st/2nd trial in each block
        lungh = length(d[d$pp == j & d$blockNum == jj, "trialNum"])
        end = strt+lungh-(Lag+1)
        for (i in strt:end){
          if (type == "error"){if (d[i-Lag, var] == 1){d[[varName]][i] <- 1} else {d[[varName]][i] <- 0}}
          else{
            if (values == "0,1"){ # if you want 0 and 1 coding
              if (d[i, var] == d[i-Lag, var]){d[[varName]][i] <- 0} else {d[[varName]][i] <- 1}
            } else if (values == "center"){ # if you want to center predictors
              if (d[i, var] == d[i-Lag, var]){d[[varName]][i] <- 0.5} else {d[[varName]][i] <- - 0.5}
            }
          }
        }
      }
    }
    # print table to see results
    print(table(d[[varName]]))
    #tests, must be all 0s
    if (type == "error"){
      wa <- which(d[, var] == 1) # find the indexes of errors trials
      we <- which(d[, varName] == 1) # find the indexes of the post-errors trials
      we1 <- we-Lag # align the post-error to their errors, decreasing the index by the same lag used to create
      # them. Now the index sets should overlap, for each wa there must be a we1, unless that is the last trial 
      # a block: The row indexes present in wa and not in we1 should yield no remainder if divided by maxTrialNum
      #sum_test <- sum(setdiff(wa, we1)%%(blockLength) != 0)
      #new version:
      #sum_test <- sum(d[setdiff(wa, we1), "trialNum"]%%maxTrialNum != 0)
      ## try with the same sum test as for non-error variables...
      sum_test <- sum(d[d[i, varName] == 99, "trialNum"]%%(maxTrialNum) != 0)
      #print(setdiff(wa, we1))
      if (sum_test != 0){cat(sum_test, "times problem assigning lag relation to var ", varName, "\n")
        } else {cat("test setdiff ok for variable ", varName,"\n")}
    } else {
      sum_test <- sum(d[d[i, varName] == 99, "trialNum"]%%(maxTrialNum) != 0)
      if (sum_test != 0){cat("problem assigning lag relation to var ", varName, "\n")}
      else {cat("test 99 ok for variable", varName,"\n")}
    }
    if (sum(d[[varName]] == 99) != length(unique(d$pp))*length(unique(d$blockNum))*Lag){
      cat("weird number of 99 in", varName,"\n")
    }
  }
  return(d)
}


# new seq relation without trial num
sequence_relation_new <- function(d, variab, suffix = "R", type = "other", Lag = 1, values = c(0,1)){
  for (var in variab){
    varName <- paste(var, suffix, sep="_")
    #print(d[1:10,var])
    d[[varName]] <- 99
    for (j in unique(d$pp)){
      for (jj in unique(d$blockNum)){
        strt = min(which(d$pp == j & d$blockNum == jj))+Lag #no for 1st/2nd trial in each block
        lungh = length(which(d$pp == j & d$blockNum == jj))
        end = strt+lungh-(Lag+1)
        # loop pver the df
        for (i in strt:end){
          if (type == "error"){if (d[i-Lag, var] == 1){d[[varName]][i] <- 1} else {d[[varName]][i] <- 0}}
          else{
             if (values[1] == "center"){ # if you want to center predictors
              if (d[i, var] == d[i-Lag, var]){d[[varName]][i] <- 0.5} else {d[[varName]][i] <- - 0.5
              }
             }
               # if vector given 
              else { # if you want 0 and 1 coding
                if (d[i, var] == d[i-Lag, var]){
                  d[[varName]][i] <- values[1]
                } else {
                  d[[varName]][i] <- values[2]
                }
            }
          }
        }
      }
    }
    # print table to see results
    print(table(d[[varName]]))
    #tests, must be all 0s
  }
  return(d)
}
# group_my ----------------------------------------------------------------------

group_my <- function(d, dv, ...){
  groupVar <- enquos(...)
  dv <- enquo(dv)
  colName <- paste0("mean", quo_name(dv))
  print(paste0("The generated column containing the means is called: ", colName))
    d %>% group_by(!!!groupVar) %>% 
      summarise(!!colName := mean(!!dv), n = n(), se := sd(!!dv)/sqrt(n))%>% ungroup() 
}

# Prompt a Question ------------------------------------------------
# function to prompt a question for the user, returns the user input
remove_him <- function(j){
  remove_yn <- readline(
    prompt = paste("Do you want to reject participant", j, "? Type 1 for yes, 0 for no and press enter")
  )
  return(remove_yn)
}

# Factor --> Numeric -----------------------------------------------------
# as.numeric(factor) converts to number the factor levels, this doesn't:
factor2numeric <- function(x) {as.numeric(levels(x))[x]}

# Get Last element --------------------------------------------------------

last <- function(x) { return( x[length(x)] ) }

# Multiple Comparisons -----------------------------------------------------------
# Extract and name many means. Run many t-tests and get corrected p-values. Get a nice table
# Example NOT fun!
# Potresti farci una funzione che prenda in unput la lista di post-hoc e i nomi delle variabili...

# # calculate means in each condition
# meansXpp <- as.data.frame(group_my(drt, rt, pp, task_R, correctResp_R, colour_R, cocoa))
# 
# # retrieve vectors of the means
# lev = c(0,1)
# for (t in lev){for (r in lev){for (c in lev){for (co in c(0, 300)){
#   nam <- paste0("t", t, ".r", r, ".c", c, ".", co)
#   assign(nam, meansXpp[meansXpp$task_R == t & meansXpp$correctResp_R == r & meansXpp$colour_R == c & 
#                          meansXpp$cocoa == co, "meanrt"])
# }}}}
# 
# # Collect tests in a list
# postHocLst <- list(
#   # dissect task x resp interaction
#   "task switch costs modulated by resp in synchr. context rep" =
#     t.test(t1.r0.c0.0 - t0.r0.c0.0, t1.r1.c0.0 - t0.r1.c0.0, var.equal = T, paired =  T),
#   "task switch costs modulated by resp in synchr. context switch" =
#     t.test(t1.r0.c1.0 - t0.r0.c1.0, t1.r1.c1.0 - t0.r1.c1.0, var.equal = T, paired =  T)
# )

# Fun for adjust p-values and return post-hoc tests details
postHoc.output <- function (postHocLst, correction = "bonferroni"){
  
  # Adjust p-values for multiple comparisons
  # Collect p-values from the post hoc list
  pvec <- vector()
  for (ii in 1:length(postHocLst)) {pvec <- c(pvec, postHocLst[[ii]]$p.value)}
  # Adjust them with the method you want
  adjPValues <- p.adjust(pvec, method = correction)
  
  # Collect post-hocs in a df
  postHoc_info <- c("Comparison", "Mean of differences", "df", "t-value", "adj. p-value", "p-value")
  postHocDf <- data.frame(matrix(NA, nrow = length(postHocLst), ncol = length(postHoc_info)))
  names(postHocDf) <- postHoc_info
  
  # Fill in the columns of the df
  postHocDf[,"adj. p-value"] <- round(adjPValues, 4)
  
  for (ii in 1:length(postHocLst)){
    postHocDf[ii,"Comparison"] <- names(postHocLst)[ii]
    #postHocDf[ii,"Comparison"] <- postHocLst[[ii]]$data.name
    postHocDf[ii,"Mean of differences"] <- round(postHocLst[[ii]]$estimate, 2)
    postHocDf[ii,"df"] <- postHocLst[[ii]]$parameter
    postHocDf[ii,"t-value"] <- round(postHocLst[[ii]]$statistic, 2)
    postHocDf[ii,"p-value"] <- round(postHocLst[[ii]]$p.value, 4)
  }
  postHocDf
 }


# Custom ggplot2 them
my_theme <- theme_minimal() + theme(axis.title = element_text(face = "bold", colour = "black"),
                                    text = element_text(size=16, family="serif", colour = "black"),
                                    axis.line = element_line(colour = "black"))
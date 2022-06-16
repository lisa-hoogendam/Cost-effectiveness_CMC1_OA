# Shiny microsimulation inclusief markov plot
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
rm(list = ls())  # Delete everything that is in R's memory

paths <- list()

paths[["current"]] <- getwd()

library(expm)
library(diagram)
#library(readxl)
library(tidyverse)
library(shiny)

library(plotrix)

hs_7 <- c("EC", "SC", "FC", "ES", "SS", "FS", "D")
ns_7 <- length(hs_7)
# 
# # Declare variables to test function
d.r <- 0.03
n.i <- 10000
n.t <- 10
Start_healthstate_EC <- 0.8

p.EC_FC <- 0.48
p.SC_ES <- 0.04
p.SC_ES.2 <- 0.01
p.FC_ES <- 0.16
p.FC_ES.2 <- 0.03
p.ES_FS <- 0.35
p.herOK <- 0.02
p.succ_herOK <- 0.32
p.herOK_FS <- p.herOK/p.ES_FS
p.FS_SS <- p.herOK_FS * p.succ_herOK

c.EC <- 2374
c.SC <- 3894
c.FC <- 3894
c.ES <- 4058
c.SS <- 3884
c.FS <- 3884
c.Trt_cons <- 473
c.Trt_sur <- 2439
c.RTW <- 5811

u.EC <- 0.73
u.SC <- 0.80
u.FC <- 0.73
u.ES <- 0.67
u.SS <- 0.84
u.FS <- 0.73
t <- 1
# M_it <- m.M[, t]

#not available on server ==> rstudioapi
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the folder where the file is saved
getwd()

# This function provides input for the model-------
fun_SSS <- function(p.EC_FC, 
                    p.SC_ES, 
                    p.SC_ES.2, 
                    p.FC_ES, 
                    p.FC_ES.2, 
                    p.ES_FS, 
                    p.herOK, 
                    p.succ_herOK, 
                    d.r, 
                    c.EC, 
                    c.SC, 
                    c.FC, 
                    c.ES, 
                    c.SS, 
                    c.FS, 
                    c.Trt_cons, 
                    c.Trt_sur,
                    c.RTW,
                    u.EC, 
                    u.SC, 
                    u.FC, 
                    u.ES, 
                    u.SS, 
                    u.FS, 
                    Start_healthstate_EC, 
                    n.t, 
                    n.i){
  set.seed(1987)
  
  p.mort <- read.csv(paste0(paths[["current"]], "/mortProb.csv"), sep = ",", header = T) #Probability mortality based on age +    sex in the Netherlands
  
  dist.Age <- read.csv(paste0(paths[["current"]], "/probabilities_age_cmc1.csv"), sep = ",", header = T) 
  dist.Age = dist.Age %>%
    dplyr::rename(ID = "X") %>%
    dplyr::rename(Age = "Var1")%>%
    dplyr::rename(prop = "Freq")
  #p.herOK <- 0.04 # Artikel Wilkens (evt 0.04*2/3 nemen?) Aanname: alle her-OKs vinden plaats in het eerste jaar, anders nieuwe health state aanmaken
  #p.succ_herOK <- 0.5 # Tevredenheid met het resultaat na herOK, geen idee -> vandaar 50/50
  
  # # To create a file with distribution for a numeric variable
  # tab <- as.data.frame(prop.table(table(data_intake_3m_wide$Score_Pijn_Intake)))
  
  
  v.n <- c("EC", "SC", "FC", "ES", "SS", "FS", "D")  # 7-health states 
  n.s <- length(v.n)    
  v.dw <- 1 / (1 + d.r) ^ (0:n.t)  # discount weight 
  # (equal discounting is assumed for costs and effects)
  v.Age0   <- sample(x = dist.Age$Age, prob = dist.Age$prop, size = n.i, replace = TRUE)
  v.Sex    <- sample(x = c("Female", "Male"), prob = c(0.7676349, 0.2323651), size = n.i, replace = TRUE) # 185/241 and 56/241, data_intake_wide dataset december --> even updaten met data februari
  
  df <- data.frame(ID = 1:n.i, Age = v.Age0, Sex = v.Sex) %>%
    plyr::join(p.mort, by=c("Age", "Sex"))
  p.D <- df$p.HD
  p.herOK_FS <- p.herOK/p.ES_FS
  p.FS_SS <- p.herOK_FS * p.succ_herOK
  c.D <- 0 
  c.omrekenfactor <- p.FS_SS/(1-p.FS_SS)
  u.D <- 0
  
  v.M_Init  <- sample(x=c("EC","ES"), prob=c(Start_healthstate_EC, 1-Start_healthstate_EC), size=n.i, replace=TRUE)       # a vector with the initial health state for all individuals 
  v.Ts_Init <- rep(0, n.i)        # a vector with the time of being sick at the start of the model 
  
  samplev <- function(m.Probs, m) {
    # Arguments
    # m.Probs: matrix with probabilities (n.i * n.s)
    # m:       number of states than need to be sampled per individual  
    # Return
    # ran:    n.i x m matrix filled with sampled health state(s) per individual
    
    d <- dim(m.Probs)  # dimensions of the matrix filled with the multinomical probabilities for the health states 
    n <- d[1]          # first dimension - number of rows (number of individuals to sample for)
    k <- d[2]          # second dimension - number of columns (number of health states considered)
    lev <- dimnames(m.Probs)[[2]]  # extract the names of the health states considered for sampling
    if (!length(lev))  # in case names for the health states are missing, use numbers to specify the health states
      lev <- 1:k       # create a sequence from 1:k (number of health states considered)
    # create a matrix 
    ran <- matrix(lev[1], ncol = m, nrow = n) # create the matrix ran, filled with the first health state of the levels 
    U <- t(m.Probs)    # transposed m.Probs matrix n.i x n.s --> n.s x n.i 
    
    for(i in 2:k) {    # start loop, from the 2nd health states
      U[i, ] <- U[i, ] + U[i - 1, ] # start summing the probabilities of the different health states per individual 
    }
    if (any((U[k, ] - 1) > 1e-05))  # sum of all probs per individual - 1 should be 0 (use 1e-05 for rounding issues), else print the error statement
      stop("error in multinom: probabilities do not sum to 1")
    
    for (j in 1:m) {   # start loop of the state that needs to be sampled (m)
      un <- rep(runif(n), rep(k, n))       # sample from a uniform distribution of length n*k
      ran[, j] <- lev[1 + colSums(un > U)] # store the health state at the jth column of the U matrix
    }
    ran # return the new health state per individual n.i x m
  } # close the function 
  
  # CONSTRUCT PROBABILITY FUNCTION
  # Probs function outputs transition probabilities for next cycle
  
  # Application of the samplev() function requires slight modification of the function updating the transition probabilities. This modification is necessary for the matrix of probabilities to be updated for a cohort of individuals rather than a specific individual. Hence, the output of the Probs() function is now a matrix of n.sim x n.s. Since the samplev() function samples the state individuals conditional on their probability, state names are required as columnames 
  set.seed(1987)                            # set the seed 
  
  
  #### Probability function
  # The Probs function that updates the transition probabilities of every cycle is shown below.
  Probs <- function(M_it, M_prev) { 
    # Arguments
    # M_it: current health state
    # dur :   time in S1
    # Returns
    # p.it: probabilities for that individual at that cycle 
    
    p.it <- matrix(0, nrow = n.s, ncol = n.i)    # create matrix of state transition probabilities
    rownames(p.it) <-  v.n   # give the state names to the rows
    
    p.EC_D <- p.D[M_it == "EC"]
    p.SC_D <- p.D[M_it == "SC"]
    p.SC_D.2 <- p.D[M_it == "SC" & M_prev == "SC"]
    p.FC_D <- p.D[M_it == "FC"]
    p.FC_D.2 <- p.D[M_it == "FC" & M_prev == "FC"]
    p.ES_D <- p.D[M_it == "ES"]
    p.SS_D <- p.D[M_it == "SS"]
    p.FS_D <- p.D[M_it == "FS"]
    p.FS_D.2 <- p.D[M_it == "FS" & M_prev == "FS"] # Subset de vector voor alle mensen die nu in FS zitten en in de vorige cycle. Dit is essentieel voor rbind/cbind (geeft anders error over "is not a multiple of replacement length" oid)
    
    # update p.it with the appropriate probabilities   
    # p.it[, M_it == "EC"] <- rbind(0, 1 - p.EC_D - p.EC_FC, p.EC_FC, 0, 0, 0, p.EC_D) # Checkt per kolom index wat de bijbehorende waarde van M_it is. Als M_it == EC is TRUE, dan wordt die kolom gevuld met deze probabilities
    # p.it[, M_it == "SC" ] <- rbind(0, 1 - p.SC_ES - p.SC_D, 0, p.SC_ES, 0, 0, p.SC_D) 
    # p.it[, M_it == "FC" ] <- rbind(0, 0, 1 - p.FC_ES - p.FC_D, p.FC_ES, 0, 0, p.FC_D) 
    # p.it[, M_it == "ES" ] <- rbind(0, 0, 0, 0, 1-p.ES_FS-p.ES_D, p.ES_FS, p.ES_D) 
    # p.it[, M_it == "SS" ] <- rbind(0, 0, 0, 0, 1-p.SS_D, 0, p.SS_D) 
    # p.it[, M_it == "FS" ] <- rbind(0, 0, 0, 0, p.FS_SS, 1-p.FS_SS-p.FS_D, p.FS_D) 
    # p.it[, (M_it == "FS") & (M_prev == "FS")] <- rbind(0, 0, 0, 0, 0, 1-p.FS_D.2, p.FS_D.2) # Als de huidige healthstate FS is EN de vorige is ook FS, dan worden de probabilities van de voorgaande regel overschreven met deze nieuwe probabilities. Hierbij vervalt p.FS_SS wegens max 2 OK mogelijkheden
    # p.it[, M_it == "D" ] <-  rbind(0, 0, 0, 0, 0, 0, 1)     # transition probabilities when Healthy 
    
    # Probabilities now conditional on staying alive
    p.it[, M_it == "EC"] <- rbind(0, (1 - p.EC_D)*(1- p.EC_FC), (1 - p.EC_D)*p.EC_FC, 0, 0, 0, p.EC_D) # Checkt per kolom index wat de bijbehorende waarde van M_it is. Als M_it == EC is TRUE, dan wordt die kolom gevuld met deze probabilities
    p.it[, M_it == "SC" ] <- rbind(0, (1 - p.SC_D)*(1 - p.SC_ES), 0, (1 - p.SC_D)*(p.SC_ES), 0, 0, p.SC_D) 
    p.it[, (M_it == "SC")&(M_prev == "SC") ] <- rbind(0, (1 - p.SC_D.2)*(1 - p.SC_ES.2), 0, (1 - p.SC_D.2)*(p.SC_ES.2), 0, 0, p.SC_D.2) 
    p.it[, M_it == "FC" ] <- rbind(0, 0, (1 - p.FC_D)*(1 - p.FC_ES), (1 - p.FC_D)*(p.FC_ES), 0, 0, p.FC_D) 
    p.it[, (M_it == "FC") & (M_prev == "FC") ] <- rbind(0, 0, (1 - p.FC_D.2)*(1 - p.FC_ES.2), (1 - p.FC_D.2)*(p.FC_ES.2), 0, 0, p.FC_D.2) 
    p.it[, M_it == "ES" ] <- rbind(0, 0, 0, 0, (1-p.ES_D)*(1-p.ES_FS), (1-p.ES_D)*(p.ES_FS), p.ES_D) 
    p.it[, M_it == "SS" ] <- rbind(0, 0, 0, 0, 1-p.SS_D, 0, p.SS_D) 
    p.it[, M_it == "FS" ] <- rbind(0, 0, 0, 0, (1-p.FS_D)*(p.FS_SS), (1-p.FS_D)*(1-p.FS_SS), p.FS_D) 
    p.it[, (M_it == "FS") & (M_prev == "FS")] <- rbind(0, 0, 0, 0, 0, 1-p.FS_D.2, p.FS_D.2) # Als de huidige healthstate FS is EN de vorige is ook FS, dan worden de probabilities van de voorgaande regel overschreven met deze nieuwe probabilities. Hierbij vervalt p.FS_SS wegens max 2 OK mogelijkheden
    p.it[, M_it == "D" ] <-  rbind(0, 0, 0, 0, 0, 0, 1)   
    return(t(p.it))    # return the transition probabilities
  }  
  
  # M_prev <- m.M_prev[ ,1]
  
  Costs <- function (M_it, M_it_1, M_it_2, p.herOK) {
    # Arguments
    # M_it: current health state (t)
    # M_it_1 : past health state (t-1)
    # M_it_2 : past health state (t-2)
    # Returns
    # cost accrued this cycle  
    # health state costs
    # Costs for prev ES and now SS works
    
    # Costs of being in a health state consist of costs due to missing work before treatment or not returning to work after treatment, therefore costs are only calculated when a patient has not reached the age of retirement yet. 
    c.it <- c()
    c.it[ M_it == "EC" & df$Age < 67] <- c.EC 
    c.it[ M_it == "EC" & df$Age >= 67] <- 0
    c.it[ M_it == "SC" & df$Age < 67] <- c.SC
    c.it[ M_it == "SC" & df$Age >= 67] <- 0
    c.it[ M_it == "FC" & df$Age < 67] <- c.FC
    c.it[ M_it == "FC" & df$Age >= 67] <- 0
    c.it[ M_it == "ES" & df$Age < 67] <- c.ES  
    c.it[ M_it == "ES" & df$Age >= 67] <- 0  
    c.it[ M_it == "SS" & df$Age < 67] <- c.SS
    c.it[ M_it == "SS" & df$Age >= 67] <- 0
    c.it[ M_it == "FS" & df$Age < 67] <- c.FS
    c.it[ M_it == "FS" & df$Age >= 67] <- 0
    c.it[ M_it == "D"] <- c.D    # costs accrued by being Healthy this cycle  return(c.it)  # return costs accrued this cycle
    
    
    c.it[(M_it == "SC") & (M_it_1 == "EC" ) & df$Age < 67] <- c.Trt_cons + c.SC 
    c.it[(M_it == "SC") & (M_it_1 == "EC" ) & df$Age >= 67] <- c.Trt_cons
    c.it[(M_it == "FC") & (M_it_1 == "EC" ) & df$Age < 67] <- c.Trt_cons + c.FC
    c.it[(M_it == "FC") & (M_it_1 == "EC" ) & df$Age >= 67] <- c.Trt_cons
    c.it[M_it == "SS" & M_it_1 == "ES" & df$Age < 67] <- c.Trt_sur + c.SS + c.RTW
    c.it[M_it == "FS" & M_it_1 == "ES" & df$Age < 67] <- c.Trt_sur + c.FS + c.RTW
    c.it[M_it == "SS" & M_it_1 == "FS" & df$Age < 67] <- c.Trt_sur + c.SS + c.RTW
    c.it[M_it == "SS" & M_it_1 == "ES" & df$Age >= 67] <- c.Trt_sur 
    c.it[M_it == "FS" & M_it_1 == "ES" & df$Age >= 67] <- c.Trt_sur 
    c.it[M_it == "SS" & M_it_1 == "FS" & df$Age >= 67] <- c.Trt_sur 
    
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "ES") & df$Age < 67]  <- (c.Trt_sur + c.FS + c.RTW)*p.herOK
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "ES") & df$Age >= 67] <- c.Trt_sur*p.herOK
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "FS") & df$Age < 67] <- c.FS
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "FS") & df$Age >= 67] <- 0
    
    return(c.it)}
  
  Direct_Costs <- function (M_it, M_it_1, M_it_2, p.herOK) {
    # Arguments
    # M_it: current health state (t)
    # M_it_1 : past health state (t-1)
    # M_it_2 : past health state (t-2)
    # Returns
    # cost accrued this cycle  
    # health state costs
    # Costs for prev ES and now SS works
    
    # 
    
    # Direct costs comprise only treatment costs
    c.it <- c()
    
    c.it[ M_it == "EC"] <- 0
    c.it[ M_it == "SC"] <- 0
    c.it[ M_it == "FC"] <- 0
    c.it[ M_it == "ES"] <- 0  
    c.it[ M_it == "SS"] <- 0
    c.it[ M_it == "FS"] <- 0
    c.it[ M_it == "D"] <- 0 
    
    c.it[M_it == "SC" & M_it_1 == "EC"] <- c.Trt_cons
    c.it[M_it == "FC" & M_it_1 == "EC"] <- c.Trt_cons
    c.it[M_it == "SS" & M_it_1 == "ES"] <- c.Trt_sur
    c.it[M_it == "FS" & M_it_1 == "ES"] <- c.Trt_sur
    c.it[M_it == "SS" & M_it_1 == "FS"] <- c.Trt_sur
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "ES")] <- c.Trt_sur*p.herOK
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "FS")] <- 0
    
    return(c.it)}
  
  Indirect_Costs <- function (M_it, M_it_1, M_it_2, p.herOK) {
    # Arguments
    # M_it: current health state (t)
    # M_it_1 : past health state (t-1)
    # M_it_2 : past health state (t-2)
    # Returns
    # cost accrued this cycle  
    # health state costs
    # Costs for prev ES and now SS works
    
    # Indirect costs comprise costs due to sick leave, therefore these costs are not calculated when a patient has reached the retirement age
    
    # Costs of being in a health state consist of costs due to missing work before treatment or not returning to work after treatment, therefore costs are only calculated when a patient has not reached the age of retirement yet. 
    c.it <- c()
    c.it[ M_it == "EC" & df$Age < 67] <- c.EC 
    c.it[ M_it == "EC" & df$Age >= 67] <- 0
    c.it[ M_it == "SC" & df$Age < 67] <- c.SC
    c.it[ M_it == "SC" & df$Age >= 67] <- 0
    c.it[ M_it == "FC" & df$Age < 67] <- c.FC
    c.it[ M_it == "FC" & df$Age >= 67] <- 0
    c.it[ M_it == "ES" & df$Age < 67] <- c.ES  
    c.it[ M_it == "ES" & df$Age >= 67] <- 0  
    c.it[ M_it == "SS" & df$Age < 67] <- c.SS
    c.it[ M_it == "SS" & df$Age >= 67] <- 0
    c.it[ M_it == "FS" & df$Age < 67] <- c.FS
    c.it[ M_it == "FS" & df$Age >= 67] <- 0
    c.it[ M_it == "D"] <- c.D    # costs accrued by being Healthy this cycle  return(c.it)  # return costs accrued this cycle
    
    
    c.it[(M_it == "SC") & (M_it_1 == "EC" ) & df$Age < 67] <- c.SC 
    c.it[(M_it == "SC") & (M_it_1 == "EC" ) & df$Age >= 67] <- 0
    c.it[(M_it == "FC") & (M_it_1 == "EC" ) & df$Age < 67] <- c.FC
    c.it[(M_it == "FC") & (M_it_1 == "EC" ) & df$Age >= 67] <- 0
    c.it[M_it == "SS" & M_it_1 == "ES" & df$Age < 67] <- c.SS + c.RTW
    c.it[M_it == "FS" & M_it_1 == "ES" & df$Age < 67] <- c.FS + c.RTW
    c.it[M_it == "SS" & M_it_1 == "FS" & df$Age < 67] <- c.SS + c.RTW
    c.it[M_it == "SS" & M_it_1 == "ES" & df$Age >= 67] <- 0 
    c.it[M_it == "FS" & M_it_1 == "ES" & df$Age >= 67] <- 0 
    c.it[M_it == "SS" & M_it_1 == "FS" & df$Age >= 67] <- 0 
    
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "ES") & df$Age < 67]  <- (c.FS + c.RTW)*p.herOK
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "ES") & df$Age >= 67] <- 0
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "FS") & df$Age < 67] <- c.FS
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "FS") & df$Age >= 67] <- 0
    
    return(c.it)}
  
  
  # The Effs function to update the utilities at every cycle.
  Effs <- function (M_it) {
    # Argmument
    # M_it: current health state
    # Retruns
    # q.it QALYs accrued this cycle
    q.it <- c() 
    q.it[ M_it == "EC"] <- u.EC        
    q.it[ M_it == "SC"] <- u.SC
    q.it[ M_it == "FC"] <- u.FC
    q.it[ M_it == "ES"] <- u.ES  
    q.it[ M_it == "SS"] <- u.SS
    q.it[ M_it == "FS"] <- u.FS
    q.it[ M_it == "D"] <- u.D  
    return(q.it)  # return the QALYs accrued this cycle
  }
  
  m.M_prev <- m.M  <- m.E <- m.Cd <- m.Ci <- m.C <-
    matrix(nrow = n.i, ncol = n.t + 1, 
           dimnames = list(paste("ind",   1:n.i, sep = " "),   # name the rows ind1, ind2, ind3, etc.
                           paste("cycle", 0:n.t, sep = " ")))  # name the columns cycle1, cycle2, cycle3, etc.
  
  
  m.M[ ,1] <- v.M_Init
  #### 04 Model proces ####
  p = Sys.time()
  
  
  # we will simulate all individuals simultaneously (no looping over individuals)
  # initial health state for all individuals
  m.M_it_1 <- m.M_it_2 <- rep("X",n.i)
  m.C[, 1] <- Costs(M_it = m.M[, 1], M_it_1 = m.M_it_1, M_it_2 = m.M_it_2, p.herOK)  # costs accrued by each individual during cycle 0
  m.Cd[, 1] <- Direct_Costs(M_it = m.M[, 1], M_it_1 = m.M_it_1, M_it_2 = m.M_it_2, p.herOK) 
  m.Ci[, 1] <- Indirect_Costs(M_it = m.M[, 1], M_it_1 = m.M_it_1, M_it_2 = m.M_it_2, p.herOK) 
  m.E[, 1] <- Effs( m.M[, 1])  # QALYs accrued by each individual during cycle 0
  
  
  # cycles 1 through n.t
  for (t in 1:n.t) { # open time loop
    # get transition probabilities based on health state at t
    if(t < 2) {m.M_prev[ ,t] <- m.M_it_1}else{m.M_prev[ ,t] <- m.M[, t - 1]}
    m.p <- Probs(m.M[, t], m.M_prev[ ,t])
    
    # sample the current health state based on transition probabilities v.p 
    m.M[, t + 1] <- samplev(m.p, 1)      # health states for all individuals during cycle t + 1
    
    
    if(t < 1) { m.M_it_1 <- rep("X", n.i) }else{m.M_it_1 <- m.M[, t]}
    if(t < 2) { m.M_it_2 <- rep("X", n.i) }else{m.M_it_2 <- m.M[, t - 1]}
    
    m.C[, t + 1] <- Costs(M_it = m.M[, t + 1], M_it_1 = m.M_it_1, M_it_2 = m.M_it_2, p.herOK)  # costs accrued by each individual  during cycle t + 1
    m.Cd[, t + 1] <- Direct_Costs(M_it = m.M[, t + 1], M_it_1 = m.M_it_1, M_it_2 = m.M_it_2, p.herOK) 
    m.Ci[, t + 1] <- Indirect_Costs(M_it = m.M[, t + 1], M_it_1 = m.M_it_1, M_it_2 = m.M_it_2, p.herOK) 
    
    m.E[, t + 1] <- Effs( m.M[, t + 1])  # QALYs accrued by each individual  during cycle t + 1
    
    df$Age <- ifelse(m.M[ , t + 1] != "D", df$Age + 1, df$Age)
    # Selects patients in state "H" at timepoint t (in m.M indicated as column t+1)
    # Joins p.mort with dataframe with updated age and with only patients who are in state "H" at timepoint t
    # Model gets stuck at t=3
    df <- plyr::join(df, p.mort, by = c("Age", "Sex"))
    p.D <- df[ ,6+t*3]
    #p.D <- df.X2[m.M[ ,t+1]!="D","p.HD"]
    
    
  } # close time loop
  
  tc <- m.C %*% v.dw
  tcd <- m.Cd %*% v.dw # Total direct costs
  tci <- m.Ci %*% v.dw # Total indirect costs
  te <- m.E %*% v.dw
  
  
  tc_hat <- mean(tc)    # average discounted cost 
  tcd_hat <- mean(tcd)    # average discounted direct cost 
  tci_hat <- mean(tci)    # average discounted indirect cost 
  te_hat <- mean(te)    # average discounted QALYs
  costs_per_qaly <- tc_hat/te_hat # Average discounted costs/QALY
  
  results <- data.frame("Mean (discounted) costs per individual" = tc_hat, 
                        "Mean (discounted) direct costs per individual" = tcd_hat,
                        "Mean (discounted) indirect costs per individual" = tci_hat,
                        "Mean (discounted) QALYs per individual" = te_hat, 
                        "Mean (discounted) costs per QALY" = costs_per_qaly)
  results
  
  #### 06 Plot ####
  #### 06.1 Histogram ####
  # Histogram showing variability in individual total costs
  # plot(density(tc), main = paste("Total cost per person"), xlab = "Cost ($)")
  
  # Histogram showing variability in individual total QALYs
  # plot(density(te), main = paste("Total QALYs per person"), xlab = "QALYs")
  
  #### 06.2 Trace ####
  # Plot the distribution of the population across healht states over time (Trace)
  # count the number of individuals in each health state at each cycle
  m.TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE)))) 
  m.TR <- m.TR / n.i                                       # calculate the proportion of individuals 
  colnames(m.TR) <- v.n                                    # name the rows of the matrix
  rownames(m.TR) <- paste("Cycle", 0:n.t, sep = " ")       # name the columns of the matrix
  
  # Plot trace of first health state
  # matplot(0:n.t, m.TR, type = 'l', 
  #         ylab = "Proportion of cohort",
  #         xlab = "Cycle",
  #         main = "Health state trace",
  #         col = 1:n.s,
  #         lty = 1:n.s)                 
  # legend("topright", v.n, col = 1:n.s, lty = 1:n.s, bty = "n")  # add a legend to the graph
  
  return(list(trace = m.TR,
              table = results))
  
}

Markov_plot_func <- reactive({
  hs_7_ordered <- c("EC", "ES", "SC", "FC", "SS", "FS", "D")
  
  stateNames <- hs_7_ordered
  
  # NB: Deze probability matrix heeft een andere volgorde dan die in de probs functie
  # Plotmat vereist dat de healthstates waar een patient in kan starten, op de eerste rij staan
  
  Oz <- matrix(c(0, 0, "", "", 0, 0, "",
                 0, 0, 0, 0, "", "", "",
                 0, "", "", 0, 0, 0, "",
                 0, "", 0, "", 0, 0, "",
                 0, 0, 0, 0, "", 0, "",
                 0, 0, 0, 0, "", "", "",
                 0, 0, 0, 0, 0, 0, 1),
               nrow=length(stateNames), byrow=TRUE)
  row.names(Oz) <- stateNames 
  colnames(Oz) <- stateNames
  Oz <- t(Oz)
  
  curve2 <- matrix(c(0, 0, 0, 0, 0, 0, 0.5, 
                     0, 0, 0, 0, 0, 0, -0.59, 
                     0, 0, 0, 0, 0, 0, 0, 
                     0, 0, 0, 0, 0, 0, 0, 
                     0, 0, 0, 0, 0, 0, 0, 
                     0, 0, 0, 0, 0, 0, 0, 
                     0, 0, 0, 0, 0, 0, 0),
                   nrow = length(stateNames), byrow = TRUE)
  row.names(curve2) <- stateNames 
  colnames(curve2) <- stateNames
  curve2 <- t(curve2)
  
  pos <- c(2,4,1) 
  Markov_plot <- plotmat(Oz, pos = pos, 
                         curve = curve2, relsize = 0.94,
                         lwd = 1, box.lwd = 1, 
                         cex.txt = 0.7, 
                         box.size = 0.1, 
                         box.type = "circle", 
                         box.prop = 0.5,
                         box.col = "light yellow",
                         arr.length=.15,
                         arr.width=.15,
                         self.cex = .4,
                         self.shifty = -0.005,
                         self.shiftx = .12,
                         shadow.size = 0,
                         main = "Health states in the microsimulation model")
  # 
  # Oz <- matrix(c(0, 0, 0, 0, 0, 0, "p.EC_D",
  #                0, 0, 0, 0, "1-p.ES_FS-p.ES_D", p.ES_FS, "p.ES_D",
  #                0, p.SC_ES, "1 - p.SC_ES - p.SC_D", 0, 0, 0, "p.SC_D",
  #                0, p.FC_ES, 0, "1 - p.FC_ES - p.FC_D", 0, 0, "p.FC_D",
  #                0, 0, 0, 0, "1-p.SS_D", 0, "p.SS_D",
  #                0, 0, 0, 0, "p.FS_SS", "1-p.FS_SS-p.FS_D", "p.FS_D",
  #                0, 0, 0, 0, 0, 0, 1),
  #              nrow=length(stateNames), byrow=TRUE)
  # row.names(Oz) <- stateNames 
  # colnames(Oz) <- stateNames
  # Oz <- t(Oz)
  # 
  # pos <- c(2,4,1)
  # 
  # Markov_plot <- plotmat(Oz, pos = pos, 
  #                        lwd = 1, box.lwd = 2, 
  #                        cex.txt = 0.7, 
  #                        box.size = 0.1, 
  #                        box.type = "circle", 
  #                        box.prop = 0.5,
  #                        box.col = "light yellow",
  #                        arr.length=.1,
  #                        arr.width=.15,
  #                        self.cex = .4,
  #                        self.shifty = -.01,
  #                        self.shiftx = .12,
  #                        main = "")
  health_state_explanation <- data.frame(EC = "Eligible for conservative treatment", 
                                         SC = "Successful conservative treatment", 
                                         FC = "Failed conservative treatment", 
                                         ES = "Eligible for surgery", 
                                         SS = "Successful surgical treatment", 
                                         FS = "Failed surgical treatment", 
                                         D = "Dead")
  
  markov_list <- list()
  markov_list[[1]] <- Markov_plot
  markov_list[[2]] <- health_state_explanation
  names(markov_list)[[1]] <- "Markov_plot"
  names(markov_list)[[2]] <- "health_state_explanation"
  return(markov_list)
})

fun_SSS_no_treatment <- function(p.EC_FC, 
                                 p.SC_ES, 
                                 p.SC_ES.2, 
                                 p.FC_ES, 
                                 p.FC_ES.2, 
                                 p.ES_FS,
                                 p.herOK, 
                                 p.succ_herOK, 
                                 d.r, 
                                 c.EC, 
                                 c.SC, 
                                 c.FC, 
                                 c.ES, 
                                 c.SS, 
                                 c.FS, 
                                 c.Trt_cons, 
                                 c.Trt_sur,
                                 c.RTW,
                                 u.EC, 
                                 u.SC, 
                                 u.FC, 
                                 u.ES, 
                                 u.SS, 
                                 u.FS, 
                                 Start_healthstate_EC = 1, 
                                 n.t, 
                                 n.i,
                                 treat_cons,
                                 treat_chi){
  set.seed(1987)
  
  p.mort <- read.csv(paste0(paths[["current"]], "/mortProb.csv"), sep = ",", header = T) #Probability mortality based on age +    sex in the Netherlands
  
  dist.Age <- read.csv(paste0(paths[["current"]], "/probabilities_age_cmc1.csv"), sep = ",", header = T) 
  dist.Age = dist.Age %>%
    dplyr::rename(ID = "X") %>%
    dplyr::rename(Age = "Var1")%>%
    dplyr::rename(prop = "Freq")
  #p.herOK <- 0.04 # Artikel Wilkens (evt 0.04*2/3 nemen?) Aanname: alle her-OKs vinden plaats in het eerste jaar, anders nieuwe health state aanmaken
  #p.succ_herOK <- 0.5 # Tevredenheid met het resultaat na herOK, geen idee -> vandaar 50/50
  
  # # To create a file with distribution for a numeric variable
  # tab <- as.data.frame(prop.table(table(data_intake_3m_wide$Score_Pijn_Intake)))
  
  
  v.n <- c("EC", "SC", "FC", "ES", "SS", "FS", "D")  # 7-health states 
  n.s <- length(v.n)    
  v.dw <- 1 / (1 + d.r) ^ (0:n.t)  # discount weight 
  # (equal discounting is assumed for costs and effects)
  v.Age0   <- sample(x = dist.Age$Age, prob = dist.Age$prop, size = n.i, replace = TRUE)
  v.Sex    <- sample(x = c("Female", "Male"), prob = c(0.7676349, 0.2323651), size = n.i, replace = TRUE) # 185/241 and 56/241, data_intake_wide dataset december --> even updaten met data februari
  
  df <- data.frame(ID = 1:n.i, Age = v.Age0, Sex = v.Sex) %>%
    plyr::join(p.mort, by=c("Age", "Sex"))
  p.D <- df$p.HD
  p.herOK_FS <- p.herOK/p.ES_FS
  p.FS_SS <- p.herOK_FS * p.succ_herOK
  c.D <- 0 
  c.omrekenfactor <- p.FS_SS/(1-p.FS_SS)
  u.D <- 0
  
  v.M_Init  <- sample(x=c("EC","ES"), prob=c(Start_healthstate_EC, 1-Start_healthstate_EC), size=n.i, replace=TRUE)       # a vector with the initial health state for all individuals 
  v.Ts_Init <- rep(0, n.i)        # a vector with the time of being sick at the start of the model 
  
  samplev <- function(m.Probs, m) {
    # Arguments
    # m.Probs: matrix with probabilities (n.i * n.s)
    # m:       number of states than need to be sampled per individual  
    # Return
    # ran:    n.i x m matrix filled with sampled health state(s) per individual
    
    d <- dim(m.Probs)  # dimensions of the matrix filled with the multinomical probabilities for the health states 
    n <- d[1]          # first dimension - number of rows (number of individuals to sample for)
    k <- d[2]          # second dimension - number of columns (number of health states considered)
    lev <- dimnames(m.Probs)[[2]]  # extract the names of the health states considered for sampling
    if (!length(lev))  # in case names for the health states are missing, use numbers to specify the health states
      lev <- 1:k       # create a sequence from 1:k (number of health states considered)
    # create a matrix 
    ran <- matrix(lev[1], ncol = m, nrow = n) # create the matrix ran, filled with the first health state of the levels 
    U <- t(m.Probs)    # transposed m.Probs matrix n.i x n.s --> n.s x n.i 
    
    for(i in 2:k) {    # start loop, from the 2nd health states
      U[i, ] <- U[i, ] + U[i - 1, ] # start summing the probabilities of the different health states per individual 
    }
    if (any((U[k, ] - 1) > 1e-05))  # sum of all probs per individual - 1 should be 0 (use 1e-05 for rounding issues), else print the error statement
      stop("error in multinom: probabilities do not sum to 1")
    
    for (j in 1:m) {   # start loop of the state that needs to be sampled (m)
      un <- rep(runif(n), rep(k, n))       # sample from a uniform distribution of length n*k
      ran[, j] <- lev[1 + colSums(un > U)] # store the health state at the jth column of the U matrix
    }
    ran # return the new health state per individual n.i x m
  } # close the function 
  
  # CONSTRUCT PROBABILITY FUNCTION
  # Probs function outputs transition probabilities for next cycle
  
  # Application of the samplev() function requires slight modification of the function updating the transition probabilities. This modification is necessary for the matrix of probabilities to be updated for a cohort of individuals rather than a specific individual. Hence, the output of the Probs() function is now a matrix of n.sim x n.s. Since the samplev() function samples the state individuals conditional on their probability, state names are required as columnames 
  set.seed(1987)                            # set the seed 
  
  
  #### Probability function
  # The Probs function that updates the transition probabilities of every cycle is shown below.
  Probs <- function(M_it, M_prev, treat_cons, treat_chi) { 
    # Arguments
    # M_it: current health state
    # dur :   time in S1
    # Returns
    # p.it: probabilities for that individual at that cycle 
    
    p.it <- matrix(0, nrow = n.s, ncol = n.i)    # create matrix of state transition probabilities
    rownames(p.it) <-  v.n   # give the state names to the rows
    
    p.EC_D <- p.D[M_it == "EC"]
    p.SC_D <- p.D[M_it == "SC"]
    p.SC_D.2 <- p.D[M_it == "SC" & M_prev == "SC"]
    p.FC_D <- p.D[M_it == "FC"]
    p.FC_D.2 <- p.D[M_it == "FC" & M_prev == "FC"]
    p.ES_D <- p.D[M_it == "ES"]
    p.SS_D <- p.D[M_it == "SS"]
    p.FS_D <- p.D[M_it == "FS"]
    p.FS_D.2 <- p.D[M_it == "FS" & M_prev == "FS"] # Subset de vector voor alle mensen die nu in FS zitten en in de vorige cycle. Dit is essentieel voor rbind/cbind (geeft anders error over "is not a multiple of replacement length" oid)
    
    # update p.it with the appropriate probabilities   
    # p.it[, M_it == "EC"] <- rbind(0, 1 - p.EC_D - p.EC_FC, p.EC_FC, 0, 0, 0, p.EC_D) # Checkt per kolom index wat de bijbehorende waarde van M_it is. Als M_it == EC is TRUE, dan wordt die kolom gevuld met deze probabilities
    # p.it[, M_it == "SC" ] <- rbind(0, 1 - p.SC_ES - p.SC_D, 0, p.SC_ES, 0, 0, p.SC_D) 
    # p.it[, M_it == "FC" ] <- rbind(0, 0, 1 - p.FC_ES - p.FC_D, p.FC_ES, 0, 0, p.FC_D) 
    # p.it[, M_it == "ES" ] <- rbind(0, 0, 0, 0, 1-p.ES_FS-p.ES_D, p.ES_FS, p.ES_D) 
    # p.it[, M_it == "SS" ] <- rbind(0, 0, 0, 0, 1-p.SS_D, 0, p.SS_D) 
    # p.it[, M_it == "FS" ] <- rbind(0, 0, 0, 0, p.FS_SS, 1-p.FS_SS-p.FS_D, p.FS_D) 
    # p.it[, (M_it == "FS") & (M_prev == "FS")] <- rbind(0, 0, 0, 0, 0, 1-p.FS_D.2, p.FS_D.2) # Als de huidige healthstate FS is EN de vorige is ook FS, dan worden de probabilities van de voorgaande regel overschreven met deze nieuwe probabilities. Hierbij vervalt p.FS_SS wegens max 2 OK mogelijkheden
    # p.it[, M_it == "D" ] <-  rbind(0, 0, 0, 0, 0, 0, 1)     # transition probabilities when Healthy 
    
    # Probabilities now conditional on staying alive
    if(!exists("treat_cons") & !exists("treat_chi")){
      p.it[, M_it == "EC"] <- rbind(0, (1 - p.EC_D)*(1- p.EC_FC), (1 - p.EC_D)*p.EC_FC, 0, 0, 0, p.EC_D) # Checkt per kolom index wat de bijbehorende waarde van M_it is. Als M_it == EC is TRUE, dan wordt die kolom gevuld met deze probabilities
      p.it[, M_it == "SC" ] <- rbind(0, (1 - p.SC_D)*(1 - p.SC_ES), 0, (1 - p.SC_D)*(p.SC_ES), 0, 0, p.SC_D) 
      p.it[, (M_it == "SC")&(M_prev == "SC") ] <- rbind(0, (1 - p.SC_D.2)*(1 - p.SC_ES.2), 0, (1 - p.SC_D.2)*(p.SC_ES.2), 0, 0, p.SC_D.2) 
      p.it[, M_it == "FC" ] <- rbind(0, 0, (1 - p.FC_D)*(1 - p.FC_ES), (1 - p.FC_D)*(p.FC_ES), 0, 0, p.FC_D) 
      p.it[, (M_it == "FC") & (M_prev == "FC") ] <- rbind(0, 0, (1 - p.FC_D.2)*(1 - p.FC_ES.2), (1 - p.FC_D.2)*(p.FC_ES.2), 0, 0, p.FC_D.2) 
      p.it[, M_it == "ES" ] <- rbind(0, 0, 0, 0, (1-p.ES_D)*(1-p.ES_FS), (1-p.ES_D)*(p.ES_FS), p.ES_D) 
      p.it[, M_it == "SS" ] <- rbind(0, 0, 0, 0, 1-p.SS_D, 0, p.SS_D) 
      p.it[, M_it == "FS" ] <- rbind(0, 0, 0, 0, (1-p.FS_D)*(p.FS_SS), (1-p.FS_D)*(1-p.FS_SS), p.FS_D) 
      p.it[, (M_it == "FS") & (M_prev == "FS")] <- rbind(0, 0, 0, 0, 0, 1-p.FS_D.2, p.FS_D.2) # Als de huidige healthstate FS is EN de vorige is ook FS, dan worden de probabilities van de voorgaande regel overschreven met deze nieuwe probabilities. Hierbij vervalt p.FS_SS wegens max 2 OK mogelijkheden
      p.it[, M_it == "D" ] <-  rbind(0, 0, 0, 0, 0, 0, 1)
      
    }else if(treat_cons == "No"){# "treat_cons" == "No" + "treat_chi" == "No" zal niet voorkomen, dat kan met Start_healthstate_EC ondervangen worden
      p.it[, M_it == "EC"] <- rbind(1 - p.EC_D, 0, 0, 0, 0, 0, p.EC_D) # Checkt per kolom index wat de bijbehorende waarde van M_it is. Als M_it == EC is TRUE, dan wordt die kolom gevuld met deze probabilities
      p.it[, M_it == "SC" ] <- rbind(0, (1 - p.SC_D)*(1 - p.SC_ES), 0, (1 - p.SC_D)*(p.SC_ES), 0, 0, p.SC_D) 
      p.it[, (M_it == "SC")&(M_prev == "SC") ] <- rbind(0, (1 - p.SC_D.2)*(1 - p.SC_ES.2), 0, (1 - p.SC_D.2)*(p.SC_ES.2), 0, 0, p.SC_D.2) 
      p.it[, M_it == "FC" ] <- rbind(0, 0, (1 - p.FC_D)*(1 - p.FC_ES), (1 - p.FC_D)*(p.FC_ES), 0, 0, p.FC_D) 
      p.it[, (M_it == "FC") & (M_prev == "FC") ] <- rbind(0, 0, (1 - p.FC_D.2)*(1 - p.FC_ES.2), (1 - p.FC_D.2)*(p.FC_ES.2), 0, 0, p.FC_D.2) 
      p.it[, M_it == "ES" ] <- rbind(0, 0, 0, 0, (1-p.ES_D)*(1-p.ES_FS), (1-p.ES_D)*(p.ES_FS), p.ES_D) 
      p.it[, M_it == "SS" ] <- rbind(0, 0, 0, 0, 1-p.SS_D, 0, p.SS_D) 
      p.it[, M_it == "FS" ] <- rbind(0, 0, 0, 0, (1-p.FS_D)*(p.FS_SS), (1-p.FS_D)*(1-p.FS_SS), p.FS_D) 
      p.it[, (M_it == "FS") & (M_prev == "FS")] <- rbind(0, 0, 0, 0, 0, 1-p.FS_D.2, p.FS_D.2)
      p.it[, M_it == "D" ] <-  rbind(0, 0, 0, 0, 0, 0, 1)
      
    }else if(treat_chi == "No"){
      p.it[, M_it == "EC"] <- rbind(0, (1 - p.EC_D)*(1- p.EC_FC), (1 - p.EC_D)*p.EC_FC, 0, 0, 0, p.EC_D) # Checkt per kolom index wat de bijbehorende waarde van M_it is. Als M_it == EC is TRUE, dan wordt die kolom gevuld met deze probabilities
      p.it[, M_it == "SC" ] <- rbind(0, (1 - p.SC_D), 0, 0, 0, 0, p.SC_D) 
      p.it[, (M_it == "SC")&(M_prev == "SC") ] <- rbind(0, (1 - p.SC_D.2), 0, 0, 0, 0, p.SC_D.2) 
      p.it[, M_it == "FC" ] <- rbind(0, 0, (1 - p.FC_D), 0, 0, 0, p.FC_D) 
      p.it[, (M_it == "FC") & (M_prev == "FC") ] <- rbind(0, 0, (1 - p.FC_D.2), 0, 0, 0, p.FC_D.2) 
      p.it[, M_it == "ES" ] <- rbind((1-p.ES_D), 0, 0, 0, 0, 0, p.ES_D) 
      p.it[, M_it == "SS" ] <- rbind(0, 0, 0, 0, 1-p.SS_D, 0, p.SS_D) 
      p.it[, M_it == "FS" ] <- rbind(0, 0, 0, 0, (1-p.FS_D)*(p.FS_SS), (1-p.FS_D)*(1-p.FS_SS), p.FS_D) 
      p.it[, (M_it == "FS") & (M_prev == "FS")] <- rbind(0, 0, 0, 0, 0, 1-p.FS_D.2, p.FS_D.2) # Als de huidige healthstate FS is EN de vorige is ook FS, dan worden de probabilities van de voorgaande regel overschreven met deze nieuwe probabilities. Hierbij vervalt p.FS_SS wegens max 2 OK mogelijkheden
      p.it[, M_it == "D" ] <-  rbind(0, 0, 0, 0, 0, 0, 1)
      
    }
    return(t(p.it))    # return the transition probabilities
  }  
  
  # M_prev <- m.M_prev[ ,1]
  
  Costs <- function (M_it, M_it_1, M_it_2, p.herOK) {
    # Arguments
    # M_it: current health state (t)
    # M_it_1 : past health state (t-1)
    # M_it_2 : past health state (t-2)
    # Returns
    # cost accrued this cycle  
    # health state costs
    # Costs for prev ES and now SS works
    
    # 
    
    
    # Costs of being in a health state consist of costs due to missing work before treatment or not returning to work after treatment, therefore costs are only calculated when a patient has not reached the age of retirement yet. 
    c.it <- c()
    c.it[ M_it == "EC" & df$Age < 67] <- c.EC 
    c.it[ M_it == "EC" & df$Age >= 67] <- 0
    c.it[ M_it == "SC" & df$Age < 67] <- c.SC
    c.it[ M_it == "SC" & df$Age >= 67] <- 0
    c.it[ M_it == "FC" & df$Age < 67] <- c.FC
    c.it[ M_it == "FC" & df$Age >= 67] <- 0
    c.it[ M_it == "ES" & df$Age < 67] <- c.ES  
    c.it[ M_it == "ES" & df$Age >= 67] <- 0  
    c.it[ M_it == "SS" & df$Age < 67] <- c.SS
    c.it[ M_it == "SS" & df$Age >= 67] <- 0
    c.it[ M_it == "FS" & df$Age < 67] <- c.FS
    c.it[ M_it == "FS" & df$Age >= 67] <- 0
    c.it[ M_it == "D"] <- c.D    # costs accrued by being Healthy this cycle  return(c.it)  # return costs accrued this cycle
    
    
    c.it[(M_it == "SC") & (M_it_1 == "EC" ) & df$Age < 67] <- c.Trt_cons + c.SC 
    c.it[(M_it == "SC") & (M_it_1 == "EC" ) & df$Age >= 67] <- c.Trt_cons
    c.it[(M_it == "FC") & (M_it_1 == "EC" ) & df$Age < 67] <- c.Trt_cons + c.FC
    c.it[(M_it == "FC") & (M_it_1 == "EC" ) & df$Age >= 67] <- c.Trt_cons
    c.it[M_it == "SS" & M_it_1 == "ES" & df$Age < 67] <- c.Trt_sur + c.SS + c.RTW
    c.it[M_it == "FS" & M_it_1 == "ES" & df$Age < 67] <- c.Trt_sur + c.FS + c.RTW
    c.it[M_it == "SS" & M_it_1 == "FS" & df$Age < 67] <- c.Trt_sur + c.SS + c.RTW
    c.it[M_it == "SS" & M_it_1 == "ES" & df$Age >= 67] <- c.Trt_sur 
    c.it[M_it == "FS" & M_it_1 == "ES" & df$Age >= 67] <- c.Trt_sur 
    c.it[M_it == "SS" & M_it_1 == "FS" & df$Age >= 67] <- c.Trt_sur 
    
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "ES") & df$Age < 67]  <- (c.Trt_sur + c.FS + c.RTW)*p.herOK
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "ES") & df$Age >= 67] <- c.Trt_sur*p.herOK
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "FS") & df$Age < 67] <- c.FS
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "FS") & df$Age >= 67] <- 0
    
    return(c.it)}
  
  Direct_Costs <- function (M_it, M_it_1, M_it_2, p.herOK) {
    # Arguments
    # M_it: current health state (t)
    # M_it_1 : past health state (t-1)
    # M_it_2 : past health state (t-2)
    # Returns
    # cost accrued this cycle  
    # health state costs
    # Costs for prev ES and now SS works
    
    # 
    
    # Direct costs comprise only treatment costs
    c.it <- c()
    
    c.it[ M_it == "EC"] <- 0
    c.it[ M_it == "SC"] <- 0
    c.it[ M_it == "FC"] <- 0
    c.it[ M_it == "ES"] <- 0  
    c.it[ M_it == "SS"] <- 0
    c.it[ M_it == "FS"] <- 0
    c.it[ M_it == "D"] <- 0 
    
    c.it[M_it == "SC" & M_it_1 == "EC"] <- c.Trt_cons
    c.it[M_it == "FC" & M_it_1 == "EC"] <- c.Trt_cons
    c.it[M_it == "SS" & M_it_1 == "ES"] <- c.Trt_sur
    c.it[M_it == "FS" & M_it_1 == "ES"] <- c.Trt_sur
    c.it[M_it == "SS" & M_it_1 == "FS"] <- c.Trt_sur
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "ES")] <- c.Trt_sur*p.herOK
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "FS")] <- 0
    
    return(c.it)}
  
  Indirect_Costs <- function (M_it, M_it_1, M_it_2, p.herOK) {
    # Arguments
    # M_it: current health state (t)
    # M_it_1 : past health state (t-1)
    # M_it_2 : past health state (t-2)
    # Returns
    # cost accrued this cycle  
    # health state costs
    # Costs for prev ES and now SS works
    
    # 
    
    # Costs of being in a health state consist of costs due to missing work before treatment or not returning to work after treatment, therefore costs are only calculated when a patient has not reached the age of retirement yet. 
    c.it <- c()
    c.it[ M_it == "EC" & df$Age < 67] <- c.EC 
    c.it[ M_it == "EC" & df$Age >= 67] <- 0
    c.it[ M_it == "SC" & df$Age < 67] <- c.SC
    c.it[ M_it == "SC" & df$Age >= 67] <- 0
    c.it[ M_it == "FC" & df$Age < 67] <- c.FC
    c.it[ M_it == "FC" & df$Age >= 67] <- 0
    c.it[ M_it == "ES" & df$Age < 67] <- c.ES  
    c.it[ M_it == "ES" & df$Age >= 67] <- 0  
    c.it[ M_it == "SS" & df$Age < 67] <- c.SS
    c.it[ M_it == "SS" & df$Age >= 67] <- 0
    c.it[ M_it == "FS" & df$Age < 67] <- c.FS
    c.it[ M_it == "FS" & df$Age >= 67] <- 0
    c.it[ M_it == "D"] <- c.D    # costs accrued by being Healthy this cycle  return(c.it)  # return costs accrued this cycle
    
    
    c.it[(M_it == "SC") & (M_it_1 == "EC" ) & df$Age < 67] <- c.SC 
    c.it[(M_it == "SC") & (M_it_1 == "EC" ) & df$Age >= 67] <- 0
    c.it[(M_it == "FC") & (M_it_1 == "EC" ) & df$Age < 67] <- c.FC
    c.it[(M_it == "FC") & (M_it_1 == "EC" ) & df$Age >= 67] <- 0
    c.it[M_it == "SS" & M_it_1 == "ES" & df$Age < 67] <- c.SS + c.RTW
    c.it[M_it == "FS" & M_it_1 == "ES" & df$Age < 67] <- c.FS + c.RTW
    c.it[M_it == "SS" & M_it_1 == "FS" & df$Age < 67] <- c.SS + c.RTW
    c.it[M_it == "SS" & M_it_1 == "ES" & df$Age >= 67] <- 0 
    c.it[M_it == "FS" & M_it_1 == "ES" & df$Age >= 67] <- 0 
    c.it[M_it == "SS" & M_it_1 == "FS" & df$Age >= 67] <- 0 
    
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "ES") & df$Age < 67]  <- (c.FS + c.RTW)*p.herOK
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "ES") & df$Age >= 67] <- 0
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "FS") & df$Age < 67] <- c.FS
    c.it[(M_it == "FS") & (M_it_1 == "FS" ) & (M_it_2 == "FS") & df$Age >= 67] <- 0
    
    return(c.it)}
  
  # The Effs function to update the utilities at every cycle.
  Effs <- function (M_it) {
    # Argmument
    # M_it: current health state
    # Retruns
    # q.it QALYs accrued this cycle
    q.it <- c() 
    q.it[ M_it == "EC"] <- u.EC        
    q.it[ M_it == "SC"] <- u.SC
    q.it[ M_it == "FC"] <- u.FC
    q.it[ M_it == "ES"] <- u.ES  
    q.it[ M_it == "SS"] <- u.SS
    q.it[ M_it == "FS"] <- u.FS
    q.it[ M_it == "D"] <- u.D  
    return(q.it)  # return the QALYs accrued this cycle
  }
  
  m.M_prev <- m.M  <- m.E <- m.Cd <- m.Ci <- m.C <-
    matrix(nrow = n.i, ncol = n.t + 1, 
           dimnames = list(paste("ind",   1:n.i, sep = " "),   # name the rows ind1, ind2, ind3, etc.
                           paste("cycle", 0:n.t, sep = " ")))  # name the columns cycle1, cycle2, cycle3, etc.
  
  
  m.M[ ,1] <- v.M_Init
  #### 04 Model proces ####
  p = Sys.time()
  
  
  # we will simulate all individuals simultaneously (no looping over individuals)
  # initial health state for all individuals
  m.M_it_1 <- m.M_it_2 <- rep("X",n.i)
  m.C[, 1] <- Costs(M_it = m.M[, 1], M_it_1 = m.M_it_1, M_it_2 = m.M_it_2, p.herOK)  # costs accrued by each individual during cycle 0
  m.Cd[, 1] <- Direct_Costs(M_it = m.M[, 1], M_it_1 = m.M_it_1, M_it_2 = m.M_it_2, p.herOK) 
  m.Ci[, 1] <- Indirect_Costs(M_it = m.M[, 1], M_it_1 = m.M_it_1, M_it_2 = m.M_it_2, p.herOK) 
  m.E[, 1] <- Effs( m.M[, 1])  # QALYs accrued by each individual during cycle 0
  
  
  # cycles 1 through n.t
  for (t in 1:n.t) { # open time loop
    # get transition probabilities based on health state at t
    if(t < 2) {m.M_prev[ ,t] <- m.M_it_1}else{m.M_prev[ ,t] <- m.M[, t - 1]}
    m.p <- Probs(M_it = m.M[, t], M_prev = m.M_prev[ ,t], treat_cons, treat_chi)
    
    # sample the current health state based on transition probabilities v.p 
    m.M[, t + 1] <- samplev(m.p, 1)      # health states for all individuals during cycle t + 1
    
    
    if(t < 1) { m.M_it_1 <- rep("X", n.i) }else{m.M_it_1 <- m.M[, t]}
    if(t < 2) { m.M_it_2 <- rep("X", n.i) }else{m.M_it_2 <- m.M[, t - 1]}
    
    m.C[, t + 1] <- Costs(M_it = m.M[, t + 1], M_it_1 = m.M_it_1, M_it_2 = m.M_it_2, p.herOK)  # costs accrued by each individual  during cycle t + 1
    m.Cd[, t + 1] <- Direct_Costs(M_it = m.M[, t + 1], M_it_1 = m.M_it_1, M_it_2 = m.M_it_2, p.herOK) 
    m.Ci[, t + 1] <- Indirect_Costs(M_it = m.M[, t + 1], M_it_1 = m.M_it_1, M_it_2 = m.M_it_2, p.herOK) 
    
    m.E[, t + 1] <- Effs( m.M[, t + 1])  # QALYs accrued by each individual  during cycle t + 1
    
    df$Age <- ifelse(m.M[ , t + 1] != "D", df$Age + 1, df$Age)
    # Selects patients in state "H" at timepoint t (in m.M indicated as column t+1)
    # Joins p.mort with dataframe with updated age and with only patients who are in state "H" at timepoint t
    # Model gets stuck at t=3
    df <- plyr::join(df, p.mort, by = c("Age", "Sex"))
    p.D <- df[ ,6+t*3]
    #p.D <- df.X2[m.M[ ,t+1]!="D","p.HD"]
    
    
  } # close time loop
  
  tc <- m.C %*% v.dw # Total costs
  tcd <- m.Cd %*% v.dw # Total direct costs
  tci <- m.Ci %*% v.dw # Total indirect costs
  te <- m.E %*% v.dw
  
  
  tc_hat <- mean(tc)    # average discounted cost 
  tcd_hat <- mean(tcd)    # average discounted direct cost 
  tci_hat <- mean(tci)    # average discounted indirect cost 
  te_hat <- mean(te)    # average discounted QALYs
  costs_per_qaly <- tc_hat/te_hat # Average discounted costs/QALY
  
  results <- data.frame("Mean (discounted) costs per individual" = tc_hat,
                        "Mean (discounted) direct costs per individual" = tcd_hat,
                        "Mean (discounted) indirect costs per individual" = tci_hat,
                        "Mean (discounted) QALYs per individual" = te_hat, 
                        "Mean (discounted) costs per QALY" = costs_per_qaly)
  results
  
  #### 06 Plot ####
  #### 06.1 Histogram ####
  # Histogram showing variability in individual total costs
  # tmp <- plot(density(tc), main = paste("Total cost per person"), xlab = "Cost ($)")
  
  # Histogram showing variability in individual total QALYs
  # tmp <- plot(density(te), main = paste("Total QALYs per person"), xlab = "QALYs")
  
  #### 06.2 Trace ####
  # Plot the distribution of the population across healht states over time (Trace)
  # count the number of individuals in each health state at each cycle
  m.TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE)))) 
  m.TR <- m.TR / n.i                                       # calculate the proportion of individuals 
  colnames(m.TR) <- v.n                                    # name the rows of the matrix
  rownames(m.TR) <- paste("Cycle", 0:n.t, sep = " ")       # name the columns of the matrix
  
  # Plot trace of first health state
  # matplot(0:n.t, m.TR, type = 'l', 
  #         ylab = "Proportion of cohort",
  #         xlab = "Cycle",
  #         main = "Health state trace",
  #         col = 1:n.s,
  #         lty = 1:n.s)                 
  # legend("topright", v.n, col = 1:n.s, lty = 1:n.s, bty = "n")  # add a legend to the graph
  
  return(list(trace = m.TR,
              table = results))
  
}


# Define UI for application ---------
ui <- fluidPage(
  
  # Application title
  titlePanel("Cost-effectiveness of four treatment strategies for thumb base osteoarthritis"),
  titlePanel(h5("Based on the paper of Hoogendam et al. (submitted)")),
  titlePanel(h5("This app facilitates cost-effectiveness calculations of four possible treatment strategies.  The values for the input parameters as used in the paper by Hoogendam et al. are shown, but can be changed to calculate cost-effectiveness of these four treatment strategies in other settings.")),
  titlePanel(h5("Note that, in contrast to the paper, probabilistic sensitivy analysis is not included.")),
  
  # Sidebar provides input
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(width = 5,
                 shinyjs::useShinyjs(),
                 id = "side-panel",
                 tabsetPanel(type = "tabs",
                             tabPanel("Model properties",
                                      numericInput(inputId = "d.r",
                                                   label = "Discount rate",
                                                   value = 0.03,
                                                   min = 0,
                                                   max = 1),  
                                      numericInput(inputId = "n.t",
                                                   label = "Number of cycles (Year)",
                                                   value = 10,
                                                   min = 1,
                                                   max = 50), 
                                      numericInput(inputId = "n.i",
                                                   label = "Number of individuals)",
                                                   value = 10000,
                                                   min = 1,
                                                   max = 1000000)),
                             tabPanel("Probabilities", 
                                      numericInput(inputId = "Start_healthstate_EC", 
                                                   label = "Proportion of cohort that starts as eligible for conservative treatment",
                                                   value = 1,
                                                   min = 0,
                                                   max = 1),
                                      numericInput(inputId = "p.EC_FC", 
                                                   label = "Probability of unsuccessful conservative treatment",
                                                   value = 0.48,
                                                   min = 0,
                                                   max = 1), 
                                      numericInput(inputId = "p.SC_ES", 
                                                   label = "Probability to elect surgery in the first year after successful conservative treatment",
                                                   value = 0.04,
                                                   min = 0,
                                                   max = 1), 
                                      numericInput(inputId = "p.SC_ES.2", 
                                                   label = "Probability to elect surgery after the first year after successful conservative treatment",
                                                   value = 0.01,
                                                   min = 0,
                                                   max = 1), 
                                      numericInput(inputId = "p.FC_ES", 
                                                   label = "Probability to elect surgery in the first year after failed conservative treatment",
                                                   value = 0.16,
                                                   min = 0,
                                                   max = 1),
                                      numericInput(inputId = "p.FC_ES.2", 
                                                   label = "Probability to elect surgery after the first year after failed conservative treatment",
                                                   value = 0.03,
                                                   min = 0,
                                                   max = 1),####
                                      numericInput(inputId = "p.ES_FS", 
                                                   label = "Probability of unsuccessful surgical treatment",
                                                   value = 0.35,
                                                   min = 0,
                                                   max = 1),
                                      numericInput(inputId = "p.herOK", 
                                                   label = "Probability of revision surgery",
                                                   value = 0.02,
                                                   min = 0,
                                                   max = 1),
                                      numericInput(inputId = "p.succ_herOK", 
                                                   label = "Probability of success with revision surgery",
                                                   value = 0.32,
                                                   min = 0,
                                                   max = 1)
                             ), 
                             tabPanel("Health state costs (€)",
                                      numericInput(inputId = "c.EC", 
                                                   label = "Societal costs for patients with CMC-1 OA symptoms prior to starting hand therapy + orthosis",
                                                   value = 2374,
                                                   min = 0,
                                                   max = 100000),
                                      numericInput(inputId = "c.SC", 
                                                   label = "Societal costs for patients who are successfully treated conservatively",
                                                   value = 3894,
                                                   min = 0,
                                                   max = 100000),
                                      numericInput(inputId = "c.FC", 
                                                   label = "Societal costs for patients who are unsuccessfully treated conservatively",
                                                   value = 3894,
                                                   min = 0,
                                                   max = 100000),
                                      numericInput(inputId = "c.ES", 
                                                   label = "Societal costs for patients who await surgical treatment for CMC-1 OA",
                                                   value = 4058,
                                                   min = 0,
                                                   max = 100000),
                                      numericInput(inputId = "c.SS", 
                                                   label = "Societal costs for patients who are successfully treated surgically",
                                                   value = 3884,
                                                   min = 0,
                                                   max = 100000),
                                      numericInput(inputId = "c.FS", 
                                                   label = "Societal costs for patients who are unsuccessfully treated surgically",
                                                   value = 3884,
                                                   min = 0,
                                                   max = 100000)
                             ),
                             tabPanel("Treatment costs (€)",
                                      numericInput(inputId = "c.Trt_cons", 
                                                   label = "Costs of conservative treatment for CMC-1 OA",
                                                   value = 473,
                                                   min = 0,
                                                   max = 100000),
                                      numericInput(inputId = "c.Trt_sur", 
                                                   label = "Costs of surgical treatment for CMC-1 OA",
                                                   value = 2439,
                                                   min = 0,
                                                   max = 10000), 
                                      numericInput(inputId = "c.RTW", 
                                                   label = "Costs of lost productivity due to recovery after surgical treatment for CMC-1 OA",
                                                   value = 5811,
                                                   min = 0,
                                                   max = 20000)
                             ),
                             tabPanel("Utilities",
                                      numericInput(inputId = "u.EC", 
                                                   label = "Utilities for patients with CMC-1 OA symptoms prior to starting hand therapy + orthosis",
                                                   value = 0.73,
                                                   min = 0,
                                                   max = 1),
                                      numericInput(inputId = "u.SC", 
                                                   label = "Utilities for patients who are successfully treated conservatively",
                                                   value = 0.80,
                                                   min = 0,
                                                   max = 1),
                                      numericInput(inputId = "u.FC", 
                                                   label = "Utilities for patients who are unsuccessfully treated conservatively",
                                                   value = 0.73,
                                                   min = 0,
                                                   max = 1),
                                      numericInput(inputId = "u.ES", 
                                                   label = "Utilities for patients who await surgical treatment",
                                                   value = 0.67,
                                                   min = 0,
                                                   max = 1), 
                                      numericInput(inputId = "u.SS", 
                                                   label = "Utilities for patients who are successfully treated surgically",
                                                   value = 0.84,
                                                   min = 0,
                                                   max = 1), 
                                      numericInput(inputId = "u.FS", 
                                                   label = "Utilities for patients who are unsuccessfully treated surgically",
                                                   value = 0.73,
                                                   min = 0,
                                                   max = 1))),
                 actionButton("button", "Run"),
                 actionButton("reset_button", "Reset input values")
    ),
    # mainPanel provides output
    mainPanel(width = 7,
              tabsetPanel(type="tabs",
                          tabPanel("Trace Plot",  
                                   #Output: Matplot ----
                                   plotOutput(outputId = "Plot1",width="500px",height="500px"),
                                   plotOutput(outputId = "Plot2",width="500px",height="500px"),
                                   plotOutput(outputId = "Plot3",width="500px",height="500px"),
                                   plotOutput(outputId = "Plot4",width="500px",height="500px")),
                          tabPanel("Cost-Effectiveness Analysis",
                                   #Output: Print CE Analysis
                                   tableOutput(outputId="Stopping1"),
                                   tableOutput(outputId="Stopping2")),
                          tabPanel("Health states overview",
                                   #Output: Plot Markov
                                   plotOutput(outputId = "Markov_plot"),
                                   tableOutput(outputId = "health_state_explanation"))
              )
    )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  observeEvent(input$button, {
    withProgress(message = 'Performing Health Economic Analysis', value = 0, {
      # CEAnalysis<-fun_SSS(p.EC_FC, 
      #                     p.SC_ES, 
      #                     p.SC_ES.2, 
      #                     p.FC_ES, 
      #                     p.FC_ES.2, 
      #                     p.ES_FS,  
      #                     p.herOK, 
      #                     p.succ_herOK, 
      #                     d.r, 
      #                     c.EC, 
      #                     c.SC, 
      #                     c.FC, 
      #                     c.ES, 
      #                     c.SS, 
      #                     c.FS, 
      #                     c.Trt_cons, 
      #                     c.Trt_sur, 
      #                     c.RTW,
      #                     u.EC, 
      #                     u.SC, 
      #                     u.FC, 
      #                     u.ES, 
      #                     u.SS, 
      #                     u.FS, 
      #                     Start_healthstate_EC, 
      #                     n.t, 
      #                     n.i)
      
      CEAnalysis<-fun_SSS(input$p.EC_FC, 
                          input$p.SC_ES, 
                          input$p.SC_ES.2, 
                          input$p.FC_ES, 
                          input$p.FC_ES.2, 
                          input$p.ES_FS,  
                          input$p.herOK, 
                          input$p.succ_herOK, 
                          input$d.r, 
                          input$c.EC, 
                          input$c.SC, 
                          input$c.FC, 
                          input$c.ES, 
                          input$c.SS, 
                          input$c.FS, 
                          input$c.Trt_cons, 
                          input$c.Trt_sur, 
                          input$c.RTW,
                          input$u.EC, 
                          input$u.SC, 
                          input$u.FC, 
                          input$u.ES, 
                          input$u.SS, 
                          input$u.FS, 
                          input$Start_healthstate_EC, 
                          input$n.t, 
                          input$n.i)
      
      # CEAnalysis_surgery_only <-fun_SSS(p.EC_FC, 
      #                                   p.SC_ES, 
      #                                   p.SC_ES.2, 
      #                                   p.FC_ES, 
      #                                   p.FC_ES.2, 
      #                                   p.ES_FS,  
      #                                   p.herOK, 
      #                                   p.succ_herOK, 
      #                                   d.r, 
      #                                   c.EC, 
      #                                   c.SC, 
      #                                   c.FC, 
      #                                   c.ES, 
      #                                   c.SS, 
      #                                   c.FS, 
      #                                   c.Trt_cons, 
      #                                   c.Trt_sur, 
      #                                   c.RTW,
      #                                   u.EC, 
      #                                   u.SC, 
      #                                   u.FC, 
      #                                   u.ES, 
      #                                   u.SS, 
      #                                   u.FS, 
      #                                   Start_healthstate_EC = 0, 
      #                                   n.t, 
      #                                   n.i)
      
      CEAnalysis_surgery_only <-fun_SSS(input$p.EC_FC, 
                                        input$p.SC_ES, 
                                        input$p.SC_ES.2, 
                                        input$p.FC_ES, 
                                        input$p.FC_ES.2, 
                                        input$p.ES_FS,  
                                        input$p.herOK, 
                                        input$p.succ_herOK, 
                                        input$d.r, 
                                        input$c.EC, 
                                        input$c.SC, 
                                        input$c.FC, 
                                        input$c.ES, 
                                        input$c.SS, 
                                        input$c.FS, 
                                        input$c.Trt_cons, 
                                        input$c.Trt_sur, 
                                        input$c.RTW,
                                        input$u.EC, 
                                        input$u.SC, 
                                        input$u.FC, 
                                        input$u.ES, 
                                        input$u.SS, 
                                        input$u.FS, 
                                        Start_healthstate_EC = 0, 
                                        input$n.t, 
                                        input$n.i)
      
      CEAnalysis_no_treatment <-fun_SSS_no_treatment(input$p.EC_FC, 
                                                     input$p.SC_ES, 
                                                     input$p.SC_ES.2, 
                                                     input$p.FC_ES, 
                                                     input$p.FC_ES.2, 
                                                     input$p.ES_FS,  
                                                     input$p.herOK, 
                                                     input$p.succ_herOK, 
                                                     input$d.r, 
                                                     input$c.EC, 
                                                     input$c.SC, 
                                                     input$c.FC, 
                                                     input$c.ES, 
                                                     input$c.SS, 
                                                     input$c.FS, 
                                                     input$c.Trt_cons, 
                                                     input$c.Trt_sur, 
                                                     input$c.RTW,
                                                     input$u.EC, 
                                                     input$u.SC, 
                                                     input$u.FC, 
                                                     input$u.ES, 
                                                     input$u.SS, 
                                                     input$u.FS, 
                                                     Start_healthstate_EC = 1, 
                                                     input$n.t, 
                                                     input$n.i,
                                                     treat_cons = "No")
      
      CEAnalysis_no_surgery <-fun_SSS_no_treatment(input$p.EC_FC, 
                                                   input$p.SC_ES, 
                                                   input$p.SC_ES.2, 
                                                   input$p.FC_ES, 
                                                   input$p.FC_ES.2, 
                                                   input$p.ES_FS,  
                                                   input$p.herOK, 
                                                   input$p.succ_herOK, 
                                                   input$d.r, 
                                                   input$c.EC, 
                                                   input$c.SC, 
                                                   input$c.FC, 
                                                   input$c.ES, 
                                                   input$c.SS, 
                                                   input$c.FS, 
                                                   input$c.Trt_cons, 
                                                   input$c.Trt_sur, 
                                                   input$c.RTW,
                                                   input$u.EC, 
                                                   input$u.SC, 
                                                   input$u.FC, 
                                                   input$u.ES, 
                                                   input$u.SS, 
                                                   input$u.FS, 
                                                   Start_healthstate_EC = 1, 
                                                   input$n.t, 
                                                   input$n.i,
                                                   treat_cons = "Yes",
                                                   treat_chi = "No")
      
      
    })
    
    output$Stopping1<-renderTable({CEA_analysis_table <- rbind(CEAnalysis_no_treatment[[2]] %>% 
                                                                 mutate(Strategy = "No treatment"),
                                                               CEAnalysis_no_surgery[[2]] %>% 
                                                                 mutate(Strategy = "Nonsurgical treatment only"),
                                                               CEAnalysis[[2]] %>% 
                                                                 mutate(Strategy = "Current protocol"),
                                                               CEAnalysis_surgery_only[[2]] %>% 
                                                                 mutate(Strategy = "Immediate surgery"))
    colnames(CEA_analysis_table)<-c("Mean (discounted) costs per individual (€)", 
                                    "Mean (discounted) direct costs per individual (€)",
                                    "Mean (discounted) indirect costs per individual (€)",
                                    "Mean (discounted) QALYs per individual", 
                                    "Mean (discounted) costs (€) per QALY",
                                    "Strategy")
    CEA_analysis_table <- CEA_analysis_table %>% 
      select(Strategy, everything())
    CEA_analysis_table},digits=2)
    
    output$Stopping2<-renderTable({CEA_analysis_table <- rbind(CEAnalysis_no_treatment[[2]] %>% 
                                                                 mutate(Strategy = "No treatment"),
                                                               CEAnalysis_no_surgery[[2]] %>% 
                                                                 mutate(Strategy = "Nonsurgical treatment only"),
                                                               CEAnalysis[[2]] %>% 
                                                                 mutate(Strategy = "Current protocol"),
                                                               CEAnalysis_surgery_only[[2]] %>% 
                                                                 mutate(Strategy = "Immediate surgery"))
    colnames(CEA_analysis_table)<-c("Mean (discounted) costs per individual (€)", 
                                    "Mean (discounted) direct costs per individual (€)",
                                    "Mean (discounted) indirect costs per individual (€)",
                                    "Mean (discounted) QALYs per individual", 
                                    "Mean (discounted) costs (€) per QALY",
                                    "Strategy")
    CEA_analysis_table$Strategy <- factor(CEA_analysis_table$Strategy, levels = c("No treatment", "Nonsurgical treatment only", "Current protocol", "Immediate surgery"))
    
    icers <- CEA_analysis_table %>% 
      rename(Cost = "Mean (discounted) costs per individual (€)",
             Effect = "Mean (discounted) QALYs per individual") %>% 
      arrange(Cost, desc(Effect)) %>% 
      mutate(Effect = round(Effect, 2)) %>% 
      mutate(Cost = round(Cost, 2)) %>% 
      select(Strategy, Cost, Effect) %>% 
      mutate(Inc_Cost = NA,
             Inc_Effect = NA,
             ICER = NA)
    icers[1, 4:6] <- "NA"
    for(i in 2:nrow(icers)){
      icers[i, "Inc_Cost"] <- round(as.numeric(icers[i, "Cost"] - icers[i-1, "Cost"]),2)
      icers[i, "Inc_Effect"] <- round(as.numeric(icers[i, "Effect"] - icers[i-1, "Effect"]),2)
      icers[i, "ICER"] <- round(as.numeric(icers[i, "Inc_Cost"])/as.numeric(icers[i, "Inc_Effect"]),2)
    }
    icers <- icers %>% 
      rename(`Incremental Cost` = Inc_Cost,
             `Incremental Effect` = Inc_Effect)
    
    },digits=2)
    
    output$Plot1 <- renderPlot({
      matplot(0:input$n.t,CEAnalysis_no_treatment$trace, type = 'l', 
              ylab = "Probability of state occupancy",
              xlab = "Cycle",
              main = "Markov Trace for No treatment",
              col = 1:ns_7,
              lty = 1:ns_7)      # create a plot of the Markov trace
      legend <- legend("topright", legend = hs_7, col=1:ns_7, lty = 1:ns_7, bty = "n") 
    })
    
    output$Plot2 <- renderPlot({
      matplot(0:input$n.t,CEAnalysis_no_surgery$trace, type = 'l', 
              ylab = "Probability of state occupancy",
              xlab = "Cycle",
              main = "Markov Trace for Nonsurgical treatment only",
              col = 1:ns_7,
              lty = 1:ns_7)      # create a plot of the Markov trace
      legend <- legend("topright", legend = hs_7, col=1:ns_7, lty = 1:ns_7, bty = "n") 
    })
    
    output$Plot3 <- renderPlot({
      matplot(0:input$n.t,CEAnalysis$trace, type = 'l', 
              ylab = "Probability of state occupancy",
              xlab = "Cycle",
              main = "Markov Trace for Current protocol",
              col = 1:ns_7,
              lty = 1:ns_7)      # create a plot of the Markov trace
      legend <- legend("topright", legend = hs_7, col=1:ns_7, lty = 1:ns_7, bty = "n") 
    })
    
    output$Plot4 <- renderPlot({
      matplot(0:input$n.t,CEAnalysis_surgery_only$trace, type = 'l', 
              ylab = "Probability of state occupancy",
              xlab = "Cycle",
              main = "Markov Trace for Immediate surgery",
              col = 1:ns_7,
              lty = 1:ns_7)      # create a plot of the Markov trace
      legend <- legend("topright", legend = hs_7, col=1:ns_7, lty = 1:ns_7, bty = "n") 
    })
    
    output$Markov_plot <- renderPlot({
      test_list <- Markov_plot_func()
      test_list$Markov_plot
    })
    
    output$health_state_explanation <- renderTable({
      test_list <- Markov_plot_func()
      test_list$health_state_explanation
    })
    
  })
  
  observeEvent(input$reset_button, {
    shinyjs::reset("side-panel")
  })
  
}
# Run the application 
shinyApp(ui = ui, server = server)




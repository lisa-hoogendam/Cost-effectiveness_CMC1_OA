# Direct and indirect costs of non-surgical and surgical treatment regimes for first carpometacarpal osteoarthritis: A cost-effectiveness analysis using microsimulation. 
# Lisa Hoogendam 

# This code is based on the following publications:
# Jalal H, et al. An Overview of R in Health Decision Sciences. 
# Med. Decis. Making. 2017; 37(3): 735-746. 
# Krijkamp EM, et al. Microsimulation modeling for health decision sciences 
# using R: a tutorial. Med. Decis. Making. 2018; 38(3): 400-422. 
# Alarid-Escudero, F, et al. Time Traveling Is Just Too Dangerous" but Some Methods Are Worth Revisiting: The Advantages of Expected Loss Curves Over Cost-Effectiveness Acceptability Curves and Frontier. 
# Value Health. 2019; 22(5): 611–618. 

rm(list = ls())  # Delete everything that is in R's memory

hs_7 <- c("EC", "SC", "FC", "ES", "SS", "FS", "D")
ns_7 <- length(hs_7)
# 
# # Declare variables to test function
d.r <- 0.03
n.i <- 10000
n.t <- 10
Start_healthstate_EC <- 1 # This parameter is used to calculate different strategies

# Declare probabilities
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

# Declare costs
c.EC <- 2374
c.SC <- 3894
c.FC <- 3894
c.ES <- 4058
c.SS <- 3884
c.FS <- 3884
c.Trt_cons <- 473
c.Trt_sur <- 2439
c.RTW <- 5811

# Declare utilities
u.EC <- 0.73
u.SC <- 0.80
u.FC <- 0.73
u.ES <- 0.67
u.SS <- 0.84
u.FS <- 0.73
t <- 1
# M_it <- m.M[, t]

# Load packages
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the folder where the file is saved 
getwd()
library(expm)
library(diagram)
library(readxl)
library(plyr)
library(dplyr)
library(shiny)

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
  
  # p.mort <- readr::read_csv("~/Documents/R/Microsim_test/mortProb.csv") #Probability mortality based on age +    sex in the Netherlands
  p.mort <- readr::read_csv("~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/mortProb.csv")
  
  dist.Age <- readr::read_csv("~/Documents/R/Microsim_test/probabilities_age_cmc1.csv") %>%
    dplyr::rename(ID = "...1") %>%
    dplyr::rename(Age = "Var1")%>%
    dplyr::rename(prop = "Freq")
  
  v.n <- c("EC", "SC", "FC", "ES", "SS", "FS", "D")  # 7-health states 
  n.s <- length(v.n)    
  v.dw <- 1 / (1 + d.r) ^ (0:n.t)  # discount weight 
  # (equal discounting is assumed for costs and effects)
  v.Age0   <- sample(x = dist.Age$Age, prob = dist.Age$prop, size = n.i, replace = TRUE)
  v.Sex    <- sample(x = c("Female", "Male"), prob = c(0.7676349, 0.2323651), size = n.i, replace = TRUE) 
  
  df <- data.frame(ID = 1:n.i, Age = v.Age0, Sex = v.Sex) %>%
    join(p.mort, by=c("Age", "Sex"))
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
    p.FS_D.2 <- p.D[M_it == "FS" & M_prev == "FS"] # Subset the vector with all individuals that are now AND the previous cycle in FS. This is essential for the rbind/cbind function, which otherwise will return the "is not a multiple of replacement length" error
    
    # update p.it with the appropriate probabilities   
    # p.it[, M_it == "EC"] <- rbind(0, 1 - p.EC_D - p.EC_FC, p.EC_FC, 0, 0, 0, p.EC_D) # This checks for each column index what the corresponding value of M_it is. When M_it == EC is TRUE, the column will be filled with the corresponding values. 
    # p.it[, M_it == "SC" ] <- rbind(0, 1 - p.SC_ES - p.SC_D, 0, p.SC_ES, 0, 0, p.SC_D) 
    # p.it[, M_it == "FC" ] <- rbind(0, 0, 1 - p.FC_ES - p.FC_D, p.FC_ES, 0, 0, p.FC_D) 
    # p.it[, M_it == "ES" ] <- rbind(0, 0, 0, 0, 1-p.ES_FS-p.ES_D, p.ES_FS, p.ES_D) 
    # p.it[, M_it == "SS" ] <- rbind(0, 0, 0, 0, 1-p.SS_D, 0, p.SS_D) 
    # p.it[, M_it == "FS" ] <- rbind(0, 0, 0, 0, p.FS_SS, 1-p.FS_SS-p.FS_D, p.FS_D) 
    # p.it[, (M_it == "FS") & (M_prev == "FS")] <- rbind(0, 0, 0, 0, 0, 1-p.FS_D.2, p.FS_D.2) # When the current health state is FS AND the previous one is also FS, the probabilities of the previous line will be overwritten by the probabilities of this line. Because we assume a maximum of 2 surgeries, p.FS_SS is no longer relevant
    # p.it[, M_it == "D" ] <-  rbind(0, 0, 0, 0, 0, 0, 1)     # transition probabilities when Healthy 
    
    # Probabilities now conditional
    p.it[, M_it == "EC"] <- rbind(0, (1 - p.EC_D)*(1- p.EC_FC), (1 - p.EC_D)*p.EC_FC, 0, 0, 0, p.EC_D) # This checks for each column index what the corresponding value of M_it is. When M_it == EC is TRUE, the column will be filled with the corresponding values. 
    p.it[, M_it == "SC" ] <- rbind(0, (1 - p.SC_D)*(1 - p.SC_ES), 0, (1 - p.SC_D)*(p.SC_ES), 0, 0, p.SC_D) 
    p.it[, (M_it == "SC")&(M_prev == "SC") ] <- rbind(0, (1 - p.SC_D.2)*(1 - p.SC_ES.2), 0, (1 - p.SC_D.2)*(p.SC_ES.2), 0, 0, p.SC_D.2) 
    p.it[, M_it == "FC" ] <- rbind(0, 0, (1 - p.FC_D)*(1 - p.FC_ES), (1 - p.FC_D)*(p.FC_ES), 0, 0, p.FC_D) 
    p.it[, (M_it == "FC") & (M_prev == "FC") ] <- rbind(0, 0, (1 - p.FC_D.2)*(1 - p.FC_ES.2), (1 - p.FC_D.2)*(p.FC_ES.2), 0, 0, p.FC_D.2) 
    p.it[, M_it == "ES" ] <- rbind(0, 0, 0, 0, (1-p.ES_D)*(1-p.ES_FS), (1-p.ES_D)*(p.ES_FS), p.ES_D) 
    p.it[, M_it == "SS" ] <- rbind(0, 0, 0, 0, 1-p.SS_D, 0, p.SS_D) 
    p.it[, M_it == "FS" ] <- rbind(0, 0, 0, 0, (1-p.FS_D)*(p.FS_SS), (1-p.FS_D)*(1-p.FS_SS), p.FS_D) 
    p.it[, (M_it == "FS") & (M_prev == "FS")] <- rbind(0, 0, 0, 0, 0, 1-p.FS_D.2, p.FS_D.2) # When the current health state is FS AND the previous one is also FS, the probabilities of the previous line will be overwritten by the probabilities of this line. Because we assume a maximum of 2 surgeries, p.FS_SS is no longer relevant
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
    df <- join(df, p.mort, by = c("Age", "Sex"))
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
  plot(density(tc), main = paste("Total cost per person"), xlab = "Cost ($)")
  
  # Histogram showing variability in individual total QALYs
  plot(density(te), main = paste("Total QALYs per person"), xlab = "QALYs")
  
  #### 06.2 Trace ####
  # Plot the distribution of the population across healht states over time (Trace)
  # count the number of individuals in each health state at each cycle
  m.TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE)))) 
  m.TR <- m.TR / n.i                                       # calculate the proportion of individuals 
  colnames(m.TR) <- v.n                                    # name the rows of the matrix
  rownames(m.TR) <- paste("Cycle", 0:n.t, sep = " ")       # name the columns of the matrix
  
  # Plot trace of first health state
  matplot(0:n.t, m.TR, type = 'l', 
          ylab = "Proportion of cohort",
          xlab = "Cycle",
          main = "Health state trace",
          col = 1:n.s,
          lty = 1:n.s)                 
  legend("topright", v.n, col = 1:n.s, lty = 1:n.s, bty = "n")  # add a legend to the graph
  
  return(list(trace = m.TR,
              table = results))
  
}

# Separate fun_SSS function for "no treatment"--------
# This function provides input for the model-------
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
  
  p.mort <- readr::read_csv("~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/mortProb.csv") #Probability mortality based on age + sex in the Netherlands
  
  dist.Age <- readr::read_csv("~/Documents/R/Microsim_test/probabilities_age_cmc1.csv") %>%
    dplyr::rename(ID = "...1") %>%
    dplyr::rename(Age = "Var1")%>%
    dplyr::rename(prop = "Freq")
  
  v.n <- c("EC", "SC", "FC", "ES", "SS", "FS", "D")  # 7-health states 
  n.s <- length(v.n)    
  v.dw <- 1 / (1 + d.r) ^ (0:n.t)  # discount weight 
  # (equal discounting is assumed for costs and effects)
  v.Age0   <- sample(x = dist.Age$Age, prob = dist.Age$prop, size = n.i, replace = TRUE)
  v.Sex    <- sample(x = c("Female", "Male"), prob = c(0.7676349, 0.2323651), size = n.i, replace = TRUE) 
  
  df <- data.frame(ID = 1:n.i, Age = v.Age0, Sex = v.Sex) %>%
    join(p.mort, by=c("Age", "Sex"))
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
    p.FS_D.2 <- p.D[M_it == "FS" & M_prev == "FS"] # Subset the vector with all individuals that are now AND the previous cycle in FS. This is essential for the rbind/cbind function, which otherwise will return the "is not a multiple of replacement length" error
    
    # update p.it with the appropriate probabilities   
    # p.it[, M_it == "EC"] <- rbind(0, 1 - p.EC_D - p.EC_FC, p.EC_FC, 0, 0, 0, p.EC_D) # This checks for each column index what the corresponding value of M_it is. When M_it == EC is TRUE, the column will be filled with the corresponding values. 
    # p.it[, M_it == "SC" ] <- rbind(0, 1 - p.SC_ES - p.SC_D, 0, p.SC_ES, 0, 0, p.SC_D) 
    # p.it[, M_it == "FC" ] <- rbind(0, 0, 1 - p.FC_ES - p.FC_D, p.FC_ES, 0, 0, p.FC_D) 
    # p.it[, M_it == "ES" ] <- rbind(0, 0, 0, 0, 1-p.ES_FS-p.ES_D, p.ES_FS, p.ES_D) 
    # p.it[, M_it == "SS" ] <- rbind(0, 0, 0, 0, 1-p.SS_D, 0, p.SS_D) 
    # p.it[, M_it == "FS" ] <- rbind(0, 0, 0, 0, p.FS_SS, 1-p.FS_SS-p.FS_D, p.FS_D) 
    # p.it[, (M_it == "FS") & (M_prev == "FS")] <- rbind(0, 0, 0, 0, 0, 1-p.FS_D.2, p.FS_D.2) # When the current health state is FS AND the previous one is also FS, the probabilities of the previous line will be overwritten by the probabilities of this line. Because we assume a maximum of 2 surgeries, p.FS_SS is no longer relevant
    # p.it[, M_it == "D" ] <-  rbind(0, 0, 0, 0, 0, 0, 1)     # transition probabilities when Healthy 
    
    # Probabilities now conditional on staying alive
    if(!exists("treat_cons") & !exists("treat_chi")){
    p.it[, M_it == "EC"] <- rbind(0, (1 - p.EC_D)*(1- p.EC_FC), (1 - p.EC_D)*p.EC_FC, 0, 0, 0, p.EC_D) # This checks for each column index what the corresponding value of M_it is. When M_it == EC is TRUE, the column will be filled with the corresponding values. 
    p.it[, M_it == "SC" ] <- rbind(0, (1 - p.SC_D)*(1 - p.SC_ES), 0, (1 - p.SC_D)*(p.SC_ES), 0, 0, p.SC_D) 
    p.it[, (M_it == "SC")&(M_prev == "SC") ] <- rbind(0, (1 - p.SC_D.2)*(1 - p.SC_ES.2), 0, (1 - p.SC_D.2)*(p.SC_ES.2), 0, 0, p.SC_D.2) 
    p.it[, M_it == "FC" ] <- rbind(0, 0, (1 - p.FC_D)*(1 - p.FC_ES), (1 - p.FC_D)*(p.FC_ES), 0, 0, p.FC_D) 
    p.it[, (M_it == "FC") & (M_prev == "FC") ] <- rbind(0, 0, (1 - p.FC_D.2)*(1 - p.FC_ES.2), (1 - p.FC_D.2)*(p.FC_ES.2), 0, 0, p.FC_D.2) 
    p.it[, M_it == "ES" ] <- rbind(0, 0, 0, 0, (1-p.ES_D)*(1-p.ES_FS), (1-p.ES_D)*(p.ES_FS), p.ES_D) 
    p.it[, M_it == "SS" ] <- rbind(0, 0, 0, 0, 1-p.SS_D, 0, p.SS_D) 
    p.it[, M_it == "FS" ] <- rbind(0, 0, 0, 0, (1-p.FS_D)*(p.FS_SS), (1-p.FS_D)*(1-p.FS_SS), p.FS_D) 
    p.it[, (M_it == "FS") & (M_prev == "FS")] <- rbind(0, 0, 0, 0, 0, 1-p.FS_D.2, p.FS_D.2) #When the current health state is FS AND the previous one is also FS, the probabilities of the previous line will be overwritten by the probabilities of this line. Because we assume a maximum of 2 surgeries, p.FS_SS is no longer relevant
    p.it[, M_it == "D" ] <-  rbind(0, 0, 0, 0, 0, 0, 1)
    
     }else if(treat_cons == "No"){# "treat_cons" == "No" + "treat_chi" == "No" will not occur, this can be prevented using Start_healthstate_EC
       p.it[, M_it == "EC"] <- rbind(1 - p.EC_D, 0, 0, 0, 0, 0, p.EC_D) # This checks for each column index what the corresponding value of M_it is. When M_it == EC is TRUE, the column will be filled with the corresponding values. 
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
       p.it[, M_it == "EC"] <- rbind(0, (1 - p.EC_D)*(1- p.EC_FC), (1 - p.EC_D)*p.EC_FC, 0, 0, 0, p.EC_D) # This checks for each column index what the corresponding value of M_it is. When M_it == EC is TRUE, the column will be filled with the corresponding values. 
       p.it[, M_it == "SC" ] <- rbind(0, (1 - p.SC_D), 0, 0, 0, 0, p.SC_D) 
       p.it[, (M_it == "SC")&(M_prev == "SC") ] <- rbind(0, (1 - p.SC_D.2), 0, 0, 0, 0, p.SC_D.2) 
       p.it[, M_it == "FC" ] <- rbind(0, 0, (1 - p.FC_D), 0, 0, 0, p.FC_D) 
       p.it[, (M_it == "FC") & (M_prev == "FC") ] <- rbind(0, 0, (1 - p.FC_D.2), 0, 0, 0, p.FC_D.2) 
       p.it[, M_it == "ES" ] <- rbind((1-p.ES_D), 0, 0, 0, 0, 0, p.ES_D) 
       p.it[, M_it == "SS" ] <- rbind(0, 0, 0, 0, 1-p.SS_D, 0, p.SS_D) 
       p.it[, M_it == "FS" ] <- rbind(0, 0, 0, 0, (1-p.FS_D)*(p.FS_SS), (1-p.FS_D)*(1-p.FS_SS), p.FS_D) 
       p.it[, (M_it == "FS") & (M_prev == "FS")] <- rbind(0, 0, 0, 0, 0, 1-p.FS_D.2, p.FS_D.2) # When the current health state is FS AND the previous one is also FS, the probabilities of the previous line will be overwritten by the probabilities of this line. Because we assume a maximum of 2 surgeries, p.FS_SS is no longer relevant
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
    df <- join(df, p.mort, by = c("Age", "Sex"))
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
  plot(density(tc), main = paste("Total cost per person"), xlab = "Cost ($)")
  
  # Histogram showing variability in individual total QALYs
  plot(density(te), main = paste("Total QALYs per person"), xlab = "QALYs")
  
  #### 06.2 Trace ####
  # Plot the distribution of the population across healht states over time (Trace)
  # count the number of individuals in each health state at each cycle
  m.TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE)))) 
  m.TR <- m.TR / n.i                                       # calculate the proportion of individuals 
  colnames(m.TR) <- v.n                                    # name the rows of the matrix
  rownames(m.TR) <- paste("Cycle", 0:n.t, sep = " ")       # name the columns of the matrix
  
  # Plot trace of first health state
  matplot(0:n.t, m.TR, type = 'l', 
          ylab = "Proportion of cohort",
          xlab = "Cycle",
          main = "Health state trace",
          col = 1:n.s,
          lty = 1:n.s)                 
  legend("topright", v.n, col = 1:n.s, lty = 1:n.s, bty = "n")  # add a legend to the graph
  
  return(list(trace = m.TR,
              table = results))
  
}

# Resultaten cost-effectiveness-------
CEAnalysis_current_protocol <-fun_SSS(p.EC_FC, 
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
                    n.i)
n.s <- 7 # Number of health states
v.n <- c("EC", "SC", "FC", "ES", "SS", "FS", "D") 
tiff(file="~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/Figures/Figure_2C-Markov_trace_current_protocol.tiff",
     units="in", width=5, height=5, res=300)
matplot(0:n.t, CEAnalysis_current_protocol$trace, type = 'l', 
        ylab = "Proportion of cohort",
        xlab = "Cycle",
        main = "Health state trace \nCurrent protocol",
        col = 1:n.s,
        lty = 1:n.s)                 
legend("topright", v.n, col = 1:n.s, lty = 1:n.s, bty = "n")
dev.off()

CEAnalysis_surgery_only <-fun_SSS(p.EC_FC, 
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
                    Start_healthstate_EC = 
                      0, 
                    n.t, 
                    n.i)
tiff(file="~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/Figures/Figure_2D-Markov_trace_immediate_surgery.tiff",
     units="in", width=5, height=5, res=300)
matplot(0:n.t, CEAnalysis_surgery_only$trace, type = 'l', 
        ylab = "Proportion of cohort",
        xlab = "Cycle",
        main = "Health state trace \nImmediate surgery",
        col = 1:n.s,
        lty = 1:n.s)                 
legend("topright", v.n, col = 1:n.s, lty = 1:n.s, bty = "n")
dev.off()

CEAnalysis_no_treatment <- fun_SSS_no_treatment(p.EC_FC, 
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
                                                treat_cons = "No"
)
tiff(file="~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/Figures/Figure_2A-Markov_trace_no_treatment.tiff",
     units="in", width=5, height=5, res=300)
matplot(0:n.t, CEAnalysis_no_treatment$trace, type = 'l', 
        ylab = "Proportion of cohort",
        xlab = "Cycle",
        main = "Health state trace \nNo treatment",
        col = 1:n.s,
        lty = 1:n.s)                 
legend("topright", v.n, col = 1:n.s, lty = 1:n.s, bty = "n")
dev.off()

CEAnalysis_no_surgery <- fun_SSS_no_treatment(p.EC_FC, 
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
                                                treat_cons = "Yes",
                                                treat_chi = "No"
)
tiff(file="~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/Figures/Figure_2B-Markov_trace_nonsurgical_treatment_only.tiff",
     units="in", width=5, height=5, res=300)
matplot(0:n.t, CEAnalysis_no_surgery$trace, type = 'l', 
        ylab = "Proportion of cohort",
        xlab = "Cycle",
        main = "Health state trace \nNonsurgical treatment only",
        col = 1:n.s,
        lty = 1:n.s)                 
legend("topright", v.n, col = 1:n.s, lty = 1:n.s, bty = "n")
dev.off()

# Outcomes_cur_protocol <- CEAnalysis_current_protocol$table %>% 
#   mutate(Strategy = "Current_protocol")
# colnames(Outcomes_cur_protocol) <- c("Mean costs per indivudual",
#                                      "Mean direct costs per indivudual",
#                                      "Mean indrect costs per indivudual",
#                                      "Mean QALYs per individual",
#                                      "Mean costs per QALY",
#                                      "Strategy")
# Outcomes_surgery_only <- CEAnalysis_surgery_only$table %>% 
#   mutate(Strategie = "Surgery only")
# colnames(Outcomes_surgery_only) <- c("Mean costs per indivudual",
#                                      "Mean direct costs per indivudual",
#                                      "Mean indrect costs per indivudual",
#                                      "Mean QALYs per individual",
#                                      "Mean costs per QALY",
#                                      "Strategy")
# Outcomes_no_treatment <- CEAnalysis_no_treatment$table %>% 
#   mutate(Strategie = "No treatment")
# colnames(Outcomes_no_treatment) <- c("Mean costs per indivudual",
#                                      "Mean direct costs per indivudual",
#                                      "Mean indrect costs per indivudual",
#                                      "Mean QALYs per individual",
#                                      "Mean costs per QALY",
#                                      "Strategy")
# Outcomes_no_surgery <- CEAnalysis_no_surgery$table %>% 
#   mutate(Strategie = "No surgery")
# colnames(Outcomes_no_surgery) <- c("Mean costs per indivudual",
#                                      "Mean direct costs per indivudual",
#                                      "Mean indrect costs per indivudual",
#                                      "Mean QALYs per individual",
#                                      "Mean costs per QALY",
#                                      "Strategy")
# 
# overview_comparison_strategies <- rbind(Outcomes_cur_protocol, 
#                                         Outcomes_surgery_only, 
#                                         Outcomes_no_treatment,
#                                         Outcomes_no_surgery)

# Functions for creating PSA dataset------
obtain_PSA_data <- function(n.sim = 10000){
parameter_names <- c("p.EC_FC", 
                      "p.SC_ES", 
                      "p.SC_ES.2", 
                      "p.FC_ES", 
                      "p.FC_ES.2", 
                      "p.ES_FS",  
                      "p.herOK", 
                      "p.succ_herOK", 
                      "c.EC", 
                      "c.SC", 
                      "c.FC", 
                      "c.ES", 
                      "c.SS", 
                      "c.FS", 
                      "c.Trt_cons", 
                      "c.Trt_sur", 
                      "c.RTW",
                      "u.EC", 
                      "u.SC", 
                      "u.FC", 
                      "u.ES", 
                      "u.SS", 
                      "u.FS")
length(parameter_names)
parameters_df <- matrix(NA, ncol = length(parameter_names), nrow = n.sim) %>% 
  as.data.frame()
names(parameters_df) <- parameter_names

set.seed(1987)
# Cases and non-cases for beta distribution
p.EC_FC <- rbeta(n.sim, 48, 52) 
p.SC_ES <- rbeta(n.sim, 4, 96)
p.SC_ES.2 <- rbeta(n.sim, 1, 99)
p.FC_ES <- rbeta(n.sim, 16, 84)
p.FC_ES.2 <- rbeta(n.sim, 3, 97)
p.ES_FS <- rbeta(n.sim, 35, 65)
p.herOK <- rbeta(n.sim, 2, 98)
p.succ_herOK <- rbeta(n.sim, 32, 68)
u.EC <- rbeta(n.sim, 73, 27)
u.SC <- rbeta(n.sim, 80, 19) #New
u.FC <- rbeta(n.sim, 73, 26) #New
u.ES <- rbeta(n.sim, 67, 33)
u.SS <- rbeta(n.sim, 84, 17) #New
u.FS <- rbeta(n.sim, 73, 27)  
c.EC <- rgamma(n.sim, 2374)
c.SC <- rgamma(n.sim, 3894)
c.FC <- rgamma(n.sim, 3894)
c.ES <- rgamma(n.sim, 4058)
c.SS <- rgamma(n.sim, 3884)
c.FS <- rgamma(n.sim, 3884)
c.Trt_cons <- rgamma(n.sim, 473)
c.Trt_sur <- rgamma(n.sim, 2439)
c.RTW <- rgamma(n.sim, 5811)

TCdisc_cur_protocol <- TCdisc_dir_cur_protocol <- TCdisc_indir_cur_protocol <- TEdisc_cur_protocol <- TCdisc_surgery_only <- TCdisc_dir_surgery_only <- TCdisc_indir_surgery_only <- TEdisc_surgery_only <- TCdisc_no_treatment <- TCdisc_dir_no_treatment <- TCdisc_indir_no_treatment <- TEdisc_no_treatment <- TCdisc_no_surgery <- TCdisc_dir_no_surgery <- TCdisc_indir_no_surgery <- TEdisc_no_surgery <- vector("numeric", n.sim)        # initiate te vector to store costs and effects

for (i in 1:n.sim){
CEAnalysis_cur_protocol_temp <- fun_SSS(p.EC_FC[i], 
        p.SC_ES[i], 
        p.SC_ES.2[i], 
        p.FC_ES[i], 
        p.FC_ES.2[i], 
        p.ES_FS[i],  
        p.herOK[i], 
        p.succ_herOK[i], 
        d.r = 0.03, 
        c.EC[i], 
        c.SC[i], 
        c.FC[i], 
        c.ES[i], 
        c.SS[i], 
        c.FS[i], 
        c.Trt_cons[i], 
        c.Trt_sur[i], 
        c.RTW[i],
        u.EC[i], 
        u.SC[i], 
        u.FC[i], 
        u.ES[i], 
        u.SS[i], 
        u.FS[i], 
        Start_healthstate_EC = 1, 
        n.t = 10, 
        n.i = 10000)

CEAnalysis_surgery_only_temp <- fun_SSS(p.EC_FC[i], 
                                        p.SC_ES[i], 
                                        p.SC_ES.2[i], 
                                        p.FC_ES[i], 
                                        p.FC_ES.2[i], 
                                        p.ES_FS[i],  
                                        p.herOK[i], 
                                        p.succ_herOK[i], 
                                        d.r = 0.03, 
                                        c.EC[i], 
                                        c.SC[i], 
                                        c.FC[i], 
                                        c.ES[i], 
                                        c.SS[i], 
                                        c.FS[i], 
                                        c.Trt_cons[i], 
                                        c.Trt_sur[i], 
                                        c.RTW[i],
                                        u.EC[i], 
                                        u.SC[i], 
                                        u.FC[i], 
                                        u.ES[i], 
                                        u.SS[i], 
                                        u.FS[i], 
                                        Start_healthstate_EC = 0, 
                                        n.t = 10, 
                                        n.i = 10000)

CEAnalysis_no_treatment_temp <- fun_SSS_no_treatment(p.EC_FC[i],
                                                     p.SC_ES[i],
                                                     p.SC_ES.2[i],
                                                     p.FC_ES[i],
                                                     p.FC_ES.2[i],
                                                     p.ES_FS[i],
                                                     p.herOK[i],
                                                     p.succ_herOK[i],
                                                     d.r = 0.03,
                                                     c.EC[i], 
                                                     c.SC[i], 
                                                     c.FC[i], 
                                                     c.ES[i], 
                                                     c.SS[i], 
                                                     c.FS[i], 
                                                     c.Trt_cons[i],
                                                     c.Trt_sur[i],
                                                     c.RTW[i],
                                                     u.EC[i], 
                                                     u.SC[i], 
                                                     u.FC[i], 
                                                     u.ES[i], 
                                                     u.SS[i], 
                                                     u.FS[i], 
                                                     Start_healthstate_EC = 1,
                                                     n.t = 10, 
                                                     n.i = 10000,
                                                     treat_cons = "No")

CEAnalysis_no_surgery_temp <- fun_SSS_no_treatment(p.EC_FC[i],
                                                     p.SC_ES[i],
                                                     p.SC_ES.2[i],
                                                     p.FC_ES[i],
                                                     p.FC_ES.2[i],
                                                     p.ES_FS[i],
                                                     p.herOK[i],
                                                     p.succ_herOK[i],
                                                     d.r = 0.03,
                                                     c.EC[i], 
                                                     c.SC[i], 
                                                     c.FC[i], 
                                                     c.ES[i], 
                                                     c.SS[i], 
                                                     c.FS[i], 
                                                     c.Trt_cons[i],
                                                     c.Trt_sur[i],
                                                     c.RTW[i],
                                                     u.EC[i], 
                                                     u.SC[i], 
                                                     u.FC[i], 
                                                     u.ES[i], 
                                                     u.SS[i], 
                                                     u.FS[i], 
                                                     Start_healthstate_EC = 1,
                                                     n.t = 10, 
                                                     n.i = 10000,
                                                     treat_cons = "Yes",
                                                     treat_chi = "No")
# Everyone starts in EC: Start_healthstate_EC = 1
# p.EC_FC = 0,
# p.EC_SC = 0 --> This parameter is not an input parameter of the fun_SSS function, but is calculated as (1 - p.EC_D)*(1- p.EC_FC)

# Only conservative:
# Everyone starts in EC: Start_healthstate_EC = 1
# p.SC_ES = 0
# p.SC_ES.2 = 0
# p.FC_ES = 0
# p.FC_ES.2 = 0
  
  TCdisc_cur_protocol[i] <- CEAnalysis_cur_protocol_temp$table$Mean..discounted..costs.per.individual
  TCdisc_dir_cur_protocol[i] <- CEAnalysis_cur_protocol_temp$table$Mean..discounted..direct.costs.per.individual
  TCdisc_indir_cur_protocol[i] <- CEAnalysis_cur_protocol_temp$table$Mean..discounted..indirect.costs.per.individual
  TEdisc_cur_protocol[i] <- CEAnalysis_cur_protocol_temp$table$Mean..discounted..QALYs.per.individual
  
  TCdisc_surgery_only[i] <- CEAnalysis_surgery_only_temp$table$Mean..discounted..costs.per.individual
  TCdisc_dir_surgery_only[i] <- CEAnalysis_surgery_only_temp$table$Mean..discounted..direct.costs.per.individual
  TCdisc_indir_surgery_only[i] <- CEAnalysis_surgery_only_temp$table$Mean..discounted..indirect.costs.per.individual
  TEdisc_surgery_only[i] <- CEAnalysis_surgery_only_temp$table$Mean..discounted..QALYs.per.individual
  
  TCdisc_no_treatment[i] <- CEAnalysis_no_treatment_temp$table$Mean..discounted..costs.per.individual
  TCdisc_dir_no_treatment[i] <- CEAnalysis_no_treatment_temp$table$Mean..discounted..direct.costs.per.individual
  TCdisc_indir_no_treatment[i] <- CEAnalysis_no_treatment_temp$table$Mean..discounted..indirect.costs.per.individual
  TEdisc_no_treatment[i] <- CEAnalysis_no_treatment_temp$table$Mean..discounted..QALYs.per.individual
  
  TCdisc_no_surgery[i] <- CEAnalysis_no_surgery_temp$table$Mean..discounted..costs.per.individual
  TCdisc_dir_no_surgery[i] <- CEAnalysis_no_surgery_temp$table$Mean..discounted..direct.costs.per.individual
  TCdisc_indir_no_surgery[i] <- CEAnalysis_no_surgery_temp$table$Mean..discounted..indirect.costs.per.individual
  TEdisc_no_surgery[i] <- CEAnalysis_no_surgery_temp$table$Mean..discounted..QALYs.per.individual
  
  parameters_df[i, "p.EC_FC"] <- p.EC_FC[i]
  parameters_df[i, "p.SC_ES"] <- p.SC_ES[i]
  parameters_df[i, "p.SC_ES.2"] <- p.SC_ES.2[i]
  parameters_df[i, "p.FC_ES"] <- p.FC_ES[i]
  parameters_df[i, "p.FC_ES.2"] <- p.FC_ES.2[i]
  parameters_df[i, "p.ES_FS"] <- p.ES_FS[i]
  parameters_df[i, "p.herOK"] <- p.herOK[i]
  parameters_df[i, "p.succ_herOK"] <- p.succ_herOK[i]
  
  parameters_df[i, "c.EC"] <- c.EC[i]
  parameters_df[i, "c.SC"] <- c.SC[i]
  parameters_df[i, "c.FC"] <- c.FC[i]
  parameters_df[i, "c.ES"] <- c.ES[i]
  parameters_df[i, "c.SS"] <- c.SS[i]
  parameters_df[i, "c.FS"] <- c.FS[i]
  parameters_df[i, "c.Trt_cons"] <- c.Trt_cons[i]
  parameters_df[i, "c.Trt_sur"] <- c.Trt_sur[i]
  parameters_df[i, "c.RTW"] <- c.RTW[i]
  
  parameters_df[i, "u.EC"] <- u.EC[i]
  parameters_df[i, "u.SC"] <- u.SC[i]
  parameters_df[i, "u.FC"] <- u.FC[i]
  parameters_df[i, "u.ES"] <- u.ES[i]
  parameters_df[i, "u.SS"] <- u.SS[i]
  parameters_df[i, "u.FS"] <- u.FS[i]
  print(i)
}
result <- list("TCdisc_cur_protocol" = TCdisc_cur_protocol,
               "TCdisc_dir_cur_protocol" = TCdisc_dir_cur_protocol,
               "TCdisc_indir_cur_protocol" = TCdisc_indir_cur_protocol,
               "TEdisc_cur_protocol" = TEdisc_cur_protocol,
               "TCdisc_surgery_only" = TCdisc_surgery_only,
               "TCdisc_dir_surgery_only" = TCdisc_dir_surgery_only,
               "TCdisc_indir_surgery_only" = TCdisc_indir_surgery_only,
               "TEdisc_surgery_only" = TEdisc_surgery_only,
               "TCdisc_no_treatment" = TCdisc_no_treatment,
               "TCdisc_dir_no_treatment" = TCdisc_dir_no_treatment,
               "TCdisc_indir_no_treatment" = TCdisc_indir_no_treatment,
               "TEdisc_no_treatment" = TEdisc_no_treatment,
               "TCdisc_no_surgery" = TCdisc_no_surgery,
               "TCdisc_dir_no_surgery" = TCdisc_dir_no_surgery,
               "TCdisc_indir_no_surgery" = TCdisc_indir_no_surgery,
               "TEdisc_no_surgery" = TEdisc_no_surgery,
               "Parameters" = parameters_df)
return(result)
}

x <- obtain_PSA_data()
save(x, file = "~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/PSA_data_290422.Rdata")
# save(x, file = "~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/PSA_data_230921.Rdata")
#load("~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/PSA_data_230921.Rdata")

# Direct and indirect costs per strategy-----
mean(x$TCdisc_dir_no_treatment)
mean(x$TCdisc_indir_no_treatment)
mean(x$TCdisc_dir_no_surgery)
mean(x$TCdisc_indir_no_surgery)
mean(x$TCdisc_dir_cur_protocol)
mean(x$TCdisc_indir_cur_protocol)
mean(x$TCdisc_dir_surgery_only)
mean(x$TCdisc_indir_surgery_only)

# Create PSA dataset ----
cost <- cbind(x$TCdisc_cur_protocol, x$TCdisc_surgery_only, x$TCdisc_no_treatment, x$TCdisc_no_surgery) %>%
  as.data.frame()
names(cost) <- c("Current protocol", "Immediate surgery", "No treatment", "Nonsurgical treatment only")
effects <- cbind(x$TEdisc_cur_protocol, x$TEdisc_surgery_only, x$TEdisc_no_treatment, x$TEdisc_no_surgery) %>% 
  as.data.frame()
names(effects) <- c("Current protocol", "Immediate surgery", "No treatment", "Nonsurgical treatment only")

# Manual adaptations
y <- list()
y[["strategies"]] <- c("Current protocol", "Immediate surgery", "No treatment", "Nonsurgical treatment only")
y[["wtp"]] <- c(1000,  11000,  21000,  31000,  41000,  51000,  61000,  71000,  81000,  91000, 101000, 111000, 121000, 131000, 141000)
# y[["wtp"]] <- c(0,  10000,  20000,  30000,  40000,  50000,  60000,  70000,  80000,  90000, 100000, 110000, 120000, 130000, 140000)
y[["wtp"]] <- seq(0, 120000, by=5000) # For the figures we used steps of 5000 for readability, but to determine the thresholds used in describing the figure, we used steps of 2500. 

y[["cost"]] <- cost
y[["effectiveness"]] <- effects
y[["parameters"]] <- x$Parameters

# x[["TCdisc_cur_protocol"]] <- NULL
# x[["TEdisc_cur_protocol"]] <- NULL
# x[["TCdisc_surgery_only"]] <- NULL
# x[["TEdisc_surgery_only"]] <- NULL
# Analyses PSA ------
# https://cran.r-project.org/web/packages/dampack/vignettes/psa_analysis.html

# To do PSA analyses, you first need a PSA dataset. This can be obtained by repeating the code above (for example 10000 times) with slightly varying inputs, sampled from standard distributions. The resulting PSA dataset has to be a named list containing "n_strategies",  "strategies"    "n_sim",         "cost" ,         "effectiveness", "other_outcome", "parameters",    "parnames",      "currency"
# Of those, cost, effectiveness and parameters must be dataframes. Parameters requires that all parameters per run are stored
library(dampack)
backup <- dampack::example_psa
example_psa <- y
psa_obj <- make_psa_obj(cost = example_psa$cost,
                        effectiveness = example_psa$effectiveness,
                        parameters = example_psa$parameters,
                        strategies = example_psa$strategies,
                        currency = "€")
plot(psa_obj)
# ggsave(file = "~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/Figure_x-Cost_effectiveness_strategies_colours.tiff")

psa_sum <- summary(psa_obj, 
                   calc_sds = TRUE)
psa_sum

icers <- calculate_icers(cost = psa_sum$meanCost, 
                         effect = psa_sum$meanEffect, 
                         strategies = psa_sum$Strategy)
b <- plot(icers)
b$labels$y <- "Cost (€)"
plot(b)
# ggsave(file = "~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/Figure_x-Cost_effectiveness_plane.tiff")

# Cost-effectiveness acceptability curve
ceac_obj <- ceac(wtp = example_psa$wtp, 
                 psa = psa_obj)
a <- plot(ceac_obj, 
     frontier = TRUE, 
     points = TRUE)
a$labels$x <- "Willingness to Pay (Thousand € / QALY)"
plot(a)
# ggsave(file = "~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/Figures/Supplemental_figure_1-CEAC.tiff")

# Expected loss curve
el <- calc_exp_loss(wtp = example_psa$wtp, 
                    psa = psa_obj)
head(el)

c <- plot(el, 
     n_x_ticks = 8, 
     n_y_ticks = 6)
c$labels$y <- "Expected Loss (€)"
c$labels$x <- "Willingness to pay (Thousand € / QALY)"
c <- c + geom_vline(xintercept = 50, colour = "red")
plot(c)
ggsave(file = "~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/Figures/Figure_3-Expected_loss_curve.tiff")
ggsave(file = "~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/Figures/Figure_BSSH-Expected_loss_curve.tiff")

# One way sensitivity analysis - Optimalizing effects
o_eff <- owsa(psa_obj, outcome = "eff")
owsa_tornado(o_eff)
owsa_tornado(o_eff, 
             min_rel_diff = 0.05) # the variable u.SS is the most "influencial" in the model
options(scipen=999)
plot(o_eff,
     n_x_ticks = 2)
opt<-  owsa_opt_strat(o_eff, 
               n_x_ticks = 5, return = "data")

# One way sensitivity analysis - Optimalizing costs 
o_cost <- owsa(psa_obj, outcome = "cost")
owsa_tornado(o_cost)
owsa_tornado(o_cost, 
             min_rel_diff = 0.05) # the variable u.SS is the most "influencial" in the model
options(scipen=999)
plot(o_cost,
     n_x_ticks = 2)
opt<-  owsa_opt_strat(o_cost, 
                      n_x_ticks = 5, return = "data")


# One way sensitivity analysis - Net health benefit
o_nhb <- owsa(psa_obj, outcome = "nhb", wtp = 50000)
owsa_tornado(o_nhb)
owsa_tornado(o_nhb, 
             min_rel_diff = 0.05) # the variables u.EC, u.FC, u.SC are the most "influencial" in the model
options(scipen=999)
xq <- plot(o_nhb,
     n_x_ticks = 2)
data_nhb <- xq$data 
data_nhb$parameter <- factor(data_nhb$parameter, levels = unique(xq$data$parameter))
params <- unique(xq$data$parameter) 
params_utility <- params[18:23]  #[1:8],[9:17],[18:23]
data_nhb_costs <- data_nhb %>% 
  filter(parameter %in% params_utility)
xq_temp <- xq
xq_temp$data <- data_nhb_costs
library(cowplot)
legend <- get_legend(xq_temp)
library(ggpubr)
as_ggplot(legend)
ggsave(file = "~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/Supplemental_figures_1-Legend.tiff")
xq_temp <- xq_temp + theme(legend.position = "none")
plot(xq_temp)
ggsave(file = "~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/Figures/Supplemental_figure_2C-OWSA_utilities.tiff")

opt<-  owsa_opt_strat(o_nhb, 
                      n_x_ticks = 5, return = "data")

# One way sensitivity analysis - Net monetary benefit
o_nmb <- owsa(psa_obj, outcome = "nmb", wtp = 50000)
owsa_tornado(o_nmb)
owsa_tornado(o_nmb, 
             min_rel_diff = 0.05) # the variables u.EC, u.FC, u.SS are the most "influencial" in the model
options(scipen=999)
plot(o_nmb,
     n_x_ticks = 2)
opt<-  owsa_opt_strat(o_nmb, 
                      n_x_ticks = 5, return = "data")


# Two way sensitivity analysis
# the variables u.EC, u.FC, u.SC are the most "influencial" in the model
tw <- twsa(psa_obj, 
           param1 = "u.EC", 
           param2 = "u.FC",
           outcome = "nhb",
           wtp = 50000)
plot(tw)
ggsave(file = "~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/Supplemental_figure_2A-TWSA_u.EC_u.FC.tiff")

tw <- twsa(psa_obj, 
           param1 = "u.EC", 
           param2 = "u.SC",
           outcome = "nhb",
           wtp = 50000)
plot(tw)
ggsave(file = "~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/Supplemental_figure_2B-TWSA_u.EC_u.SC.tiff")

tw <- twsa(psa_obj, 
           param1 = "u.FC", 
           param2 = "u.SC",
           outcome = "nhb",
           wtp = 50000)
plot(tw)
ggsave(file = "~/Documents/Promotie/Projecten Lisa/Kosten-effectiviteit CMC-1 OA behandelstrategie (Microsimulation)/Supplemental_figure_2C-TWSA_u.FC_u.SC.tiff")

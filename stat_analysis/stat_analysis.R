library(rethinking)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggprism)
library(tidybayes)
library(tidybayes.rethinking)
library(magick)

aineisto <- read.csv("growth_measure_final_plating.csv", header = T)
aineisto <- aineisto[,-2]
aineisto$Time_0 <- rep(0, length(aineisto$Time_1))
#Filter the data by generation

f2data <- filter(aineisto, generation == "F2")
f3data <- filter(aineisto, generation == "F3")

#Filter the dataset by experiment
exp1f2 <- filter(f2data, Experiment == "Exp1")
exp2f2 <- filter(f2data, Experiment == "Exp2")
exp3f2 <- filter(f2data, Experiment == "Exp3")
exp4f2 <- filter(f2data, Experiment == "Exp4")
exp5f2 <- filter(f2data, Experiment == "Exp5")
exp6f2 <- filter(f2data, Experiment == "Exp6")
exp7f2 <- filter(f2data, Experiment == "Exp7")
exp8f2 <- filter(f2data, Experiment == "Exp8")
exp9f2 <- filter(f2data, Experiment == "Exp9")

##Scaling the data
exp1f2$gr_scaled <- scale(exp1f2$Time_1)
exp1f2$col_scaled <- scale(exp1f2$mean_colonies)
exp2f2$gr_scaled <- scale(exp2f2$Time_1)
exp2f2$col_scaled <- scale(exp2f2$mean_colonies)
exp3f2$gr_scaled <- scale(exp3f2$Time_1)
exp3f2$col_scaled <- scale(exp3f2$mean_colonies)
exp4f2$gr_scaled <- scale(exp4f2$Time_1)
exp4f2$col_scaled <- scale(exp4f2$mean_colonies)
exp5f2$gr_scaled <- scale(exp5f2$Time_1)
exp5f2$col_scaled <- rep("NA", length(exp5f2$gr_scaled))
exp6f2$gr_scaled <- scale(exp6f2$Time_1)
exp6f2$col_scaled <- scale(exp6f2$mean_colonies)
exp7f2$gr_scaled <- scale(exp7f2$Time_1)
exp7f2$col_scaled <- rep("NA", length(exp7f2$gr_scaled))
exp8f2$gr_scaled <- scale(exp8f2$Time_1)
exp8f2$col_scaled <- rep("NA", length(exp8f2$gr_scaled))
exp9f2$gr_scaled <- scale(exp9f2$Time_1)
exp9f2$col_scaled <- rep("NA", length(exp9f2$gr_scaled))

#Create datasets to run ulam 
f2data1 <- rbind(exp1f2, exp2f2, exp3f2, exp4f2, exp5f2, exp6f2, exp7f2, exp8f2, exp9f2)
f2data1$F1 <- ifelse(f2data1$f1_sucrose == 1.5, 1, 2)
f2data1$env <- ifelse(f2data1$p_sucrose == 1.5, 1, 2)
f2data1$condition <- ifelse(f2data1$condition == "C1", 1 ,0) + ifelse(f2data1$condition == "C2", 2, 0) + ifelse(f2data1$condition == "C3", 3, 0) + ifelse(f2data1$condition == "C4", 4, 0)
f2data1$tube <- factor(f2data1$tube)
#Create the database without col_scaled (contains 9 experiments)
f2datcomp <- data.frame(gr_scaled = f2data1$gr_scaled, F1 = f2data1$F1, env = f2data1$env, condition = factor(f2data1$condition), tube = f2data1$tube)
#create data base just just 5 experiments but with col_scaled
f2data2 <- f2data1[complete.cases(f2data1),]
f2data2$col_scaled <- as.numeric(f2data2$col_scaled)
f2datcomp1 <- data.frame(gr_scaled = f2data2$gr_scaled, col_scaled= f2data2$col_scaled, F1 = f2data2$F1, env = f2data2$env, condition = factor(f2data2$condition), tube = factor(f2data2$tube))

ggplot(exp3f2, aes(x = factor(p_sucrose), y = gr_scaled, fill = factor(f1_sucrose))) + geom_boxplot()

ggplot(exp1f2, aes(x = condition, y = col_scaled)) + geom_boxplot()
ggplot(f2datcomp1, aes(x = factor(f1_sucrose), y = col_scaled)) + geom_boxplot()

library(dagitty)
dag <- dagitty( "dag {
E -> G
F -> G
C -> G
F -> C
}")

##Graph where spore quality is an unobserved variable
dag2 <- dagitty("dag {
Q [unobserved]
E -> G
F -> Q
Q -> G
F -> C
C -> G
}")

drawdag(dag)

impliedConditionalIndependencies(dag) #is OK!

adjustmentSets(dag, exposure="F" , outcome="G" )

adjustmentSets(dag2, exposure = "Q", outcome= "G")

#The model, for F2 silver spoon effects

model1 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition] + B[tube],
    a[condition] ~ dnorm(0,1),
    B[tube] ~ dnorm(0, sigmaB),
    sigmaB ~ dexp(1),
    sigma ~ dexp(1)
  ), data = f2datcomp, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

p <- precis(model1, depth = 2)
plot(precis(model1, depth = 2))
traceplot(model1)
trankplot(model1)
post.m1 <- extract.samples(model1)
#If current environment was 1.5%
#What is the effect of F1 environment? in model 
diff1_4_model1 <- post.m1$a[,1] - post.m1$a[,4]
mean(diff1_4_model1) #1.168023
HPDI(diff1_4_model1, prob = 0.95) #0.9127088 1.5514260
#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3_model1 <- post.m1$a[,2] - post.m1$a[,3]
mean(diff2_3_model1)
HPDI(diff2_3_model1, prob = 0.95) #0.4506935 0.9709197 
#I also run the model with priors N(0,5) and exp(5) and the results did not changed
save(p, post.m1, file = "model1.RData")
#model2 with number of colonies (viability) as a fixed factor
model2 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition] + B[tube] + Bc*col_scaled,
    a[condition] ~ dnorm(0,1),
    Bc ~ dnorm(0, 1),
    B[tube] ~ dnorm(0, sigmaB),
    sigmaB ~ dexp(1),
    sigma ~ dexp(1)
  ), data = f2datcomp1, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

precis(model2, depth = 2)
traceplot(model2)
trankplot(model2)
plot(precis(model2, depth=2))
post.m2 <- extract.samples(model2)
median(post.m2$Bc)
HPDI(post.m2$Bc, prob = 0.95)
#If current environment was 1.5%
#What is the effect of F1 environment?
diff1_4 <- post.m2$a[,1] - post.m2$a[,4]
median(diff1_4)
HPDI(diff1_4, prob = 0.95)

#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3 <- post.m2$a[,2] - post.m2$a[,3]

median(diff2_3)
HPDI(diff2_3, prob = 0.95)
#graph of the model results
a1 <- post.m1$a[,1] 
a2 <- post.m1$a[,2]
a3 <- post.m1$a[,3]
a4 <- post.m1$a[,4]


############### F2 generation each experiment separated ##################
#exp1f2
f2_exp1 <- filter(f2data1, Experiment == "Exp1")
f2_exp1_comp <- data.frame(gr_scaled = f2_exp1$gr_scaled, condition = f2_exp1$cond) 

F2.m.exp1 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = f2_exp1_comp , warmup = 1000, iter = 3000, chains = 4, log_lik = T,       
  control=list(adapt_delta=0.99))

precis(F2.m.exp1, depth = 2)
plot(precis(F2.m.exp1, depth=2))
traceplot(F2.m.exp1)
trankplot(F2.m.exp1)
post.f2.exp1 <- extract.samples(F2.m.exp1)
#If current environment was 1.5%
#What is the effect of F1 environment? 1.5%
diff1_4f2exp1 <- post.f2.exp1$a[,1] - post.f2.exp1$a[,4]
median(diff1_4f2exp1)
HPDI(diff1_4f2exp1, prob = 0.95) #1.440580 2.304831
#What is the effect of F1 environment? 0.015%
diff2_3_f2exp1 <- post.f2.exp1$a[,2] - post.f2.exp1$a[,3]
median(diff2_3_f2exp1)
HPDI(diff2_3_f2exp1, prob = 0.95) #0.4307118 1.3096301


#####exp2f2#####################################
f2_exp2 <- filter(f2data1, Experiment == "Exp2")
f2_exp2 <- f2_exp2[complete.cases(f2_exp2),]
f2_exp2$tube <- ifelse(f2_exp2$tube == "7", 1 ,0) + ifelse(f2_exp2$tube == "8", 2 ,0) + ifelse(f2_exp2$tube == "9", 3 ,0) + ifelse(f2_exp2$tube == "10", 4 ,0)+
  ifelse(f2_exp2$tube == "11", 5 ,0) + ifelse(f2_exp2$tube == "12", 6 ,0)
f2_exp2_comp <- data.frame(gr_scaled = f2_exp2$gr_scaled, F1 = f2_exp2$F1, env = f2_exp2$env, condition = f2_exp2$cond, tube = f2_exp2$tube) 


F2.m.exp2 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = f2_exp2_comp , warmup = 1000, iter = 4000, chains = 4, log_lik = T)

precis(F2.m.exp2, depth = 2)
plot(precis(F2.m.exp2, depth=2))
traceplot(F2.m.exp2)
trankplot(F2.m.exp2)
post.f2.exp2 <- extract.samples(F2.m.exp2)
#If current environment was 1.5%
#What is the effect of F1 environment? in model f2 exp2
diff1_4_f2exp2 <- post.f2.exp2$a[,1] - post.f2.exp2$a[,4]
median(diff1_4_f2exp2)
HPDI(diff1_4_f2exp2, prob = 0.95) #-0.03846668  1.66744467 
#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3_f2exp2 <- post.f2.exp2$a[,2] - post.f2.exp2$a[,3]
median(diff2_3_f2exp2)
HPDI(diff2_3_f2exp2, prob = 0.95) #-0.6327464  0.9865193 

#####exp3f2####################################
f2_exp3 <- filter(f2data1, Experiment == "Exp3")
f2_exp3$tube <- ifelse(f2_exp3$tube == "13", 1 ,0)+ ifelse(f2_exp3$tube == "14", 2 ,0) + ifelse(f2_exp3$tube == "15", 3 ,0)+ 
  ifelse(f2_exp3$tube == "16", 4 ,0)+ ifelse(f2_exp3$tube == "17", 5 ,0) + ifelse(f2_exp3$tube == "18", 6 ,0)
f2_exp3_comp <- data.frame(gr_scaled = f2_exp3$gr_scaled, F1 = f2_exp3$F1, env = f2_exp3$env, condition = f2_exp3$cond, tube = f2_exp3$tube) 

F2.m.exp3 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = f2_exp3_comp , warmup = 1000, iter = 4000, chains = 4, log_lik = T,
  control=list(adapt_delta=0.99))

precis(F2.m.exp3, depth = 2)
plot(precis(F2.m.exp3, depth=2))
traceplot(F2.m.exp3)
trankplot(F2.m.exp3)
post.f2.exp3 <- extract.samples(F2.m.exp3)
#If current environment was 1.5%
#What is the effect of F1 environment? in model f2 exp3
diff1_4_f2exp3 <- post.f2.exp3$a[,1] - post.f2.exp3$a[,4]
median(diff1_4_f2exp3)
HPDI(diff1_4_f2exp3, prob = 0.95) #0.5922468 1.8057571
#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3_f2exp3 <- post.f2.exp3$a[,2] - post.f2.exp3$a[,3]
median(diff2_3_f2exp3)
HPDI(diff2_3_f2exp3, prob = 0.95) #0.1986581 1.4310891 
#####exp4f2####################################

f2_exp4 <- filter(f2data1, Experiment == "Exp4")
f2_exp4$tube <- ifelse(f2_exp4$tube == "19", 1 ,0)+ ifelse(f2_exp4$tube == "20", 2 ,0) + ifelse(f2_exp4$tube == "21", 3 ,0)+ 
  ifelse(f2_exp4$tube == "22", 4 ,0)+ ifelse(f2_exp4$tube == "23", 5 ,0) + ifelse(f2_exp4$tube == "24", 6 ,0)
f2_exp4_comp <- data.frame(gr_scaled = f2_exp4$gr_scaled, F1 = f2_exp4$F1, env = f2_exp4$env, condition = f2_exp4$cond, tube = f2_exp4$tube) 

F2.m.exp4 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = f2_exp4_comp , warmup = 1000, iter = 4000, chains = 4, log_lik = T,
  control=list(adapt_delta=0.99))

precis(F2.m.exp4, depth = 2)
plot(precis(F2.m.exp4, depth=2))
traceplot(F2.m.exp4)
trankplot(F2.m.exp4)
post.f2.exp4 <- extract.samples(F2.m.exp4)
#If current environment was 1.5%
#What is the effect of F1 environment? in model f2 exp4
diff1_4_f2exp4 <- post.f2.exp4$a[,1] - post.f2.exp4$a[,4]
median(diff1_4_f2exp4)
HPDI(diff1_4_f2exp4, prob = 0.95) #-0.07644636  1.52123266
#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3_f2exp4 <- post.f2.exp4$a[,2] - post.f2.exp4$a[,3]
median(diff2_3_f2exp4)
HPDI(diff2_3_f2exp4, prob = 0.95)#-0.3867618  1.2356644 
#####exp5f2####################################
f2_exp5 <- filter(f2data1, Experiment == "Exp5")
f2_exp5$tube <- ifelse(f2_exp5$tube == "25", 1 ,0)+ ifelse(f2_exp5$tube == "26", 2 ,0) + ifelse(f2_exp5$tube == "27", 3 ,0)+ 
  ifelse(f2_exp5$tube == "28", 4 ,0)+ ifelse(f2_exp5$tube == "29", 5 ,0) + ifelse(f2_exp5$tube == "30", 6 ,0)
f2_exp5_comp <- data.frame(gr_scaled = f2_exp5$gr_scaled, F1 = f2_exp5$F1, env = f2_exp5$env, condition = f2_exp5$cond, tube = f2_exp5$tube) 

F2.m.exp5 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = f2_exp5_comp , warmup = 1000, iter = 4000, chains = 4, log_lik = T,
  control=list(adapt_delta=0.99))

precis(F2.m.exp5, depth = 2)
plot(precis(F2.m.exp5, depth=2))
traceplot(F2.m.exp5)
trankplot(F2.m.exp5)
post.f2.exp5 <- extract.samples(F2.m.exp5)
#If current environment was 1.5%
#What is the effect of F1 environment? in model f2 exp5
diff1_4_f2exp5 <- post.f2.exp5$a[,1] - post.f2.exp5$a[,4]
median(diff1_4_f2exp5)
HPDI(diff1_4_f2exp5, prob = 0.95)#0.6789493 1.9320295
#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3_f2exp5 <- post.f2.exp5$a[,2] - post.f2.exp5$a[,3]
median(diff2_3_f2exp5)
HPDI(diff2_3_f2exp5, prob = 0.95)#0.1461811 1.3874591
#####exp6f2####################################
f2_exp6 <- filter(f2data1, Experiment == "Exp6")
f2_exp6$tube <- ifelse(f2_exp6$tube == "31", 1 ,0)+ ifelse(f2_exp6$tube == "32", 2 ,0) + ifelse(f2_exp6$tube == "33", 3 ,0)+ 
  ifelse(f2_exp6$tube == "34", 4 ,0)+ ifelse(f2_exp6$tube == "35", 5 ,0) + ifelse(f2_exp6$tube == "36", 6 ,0)+
  ifelse(f2_exp6$tube == "37", 7 ,0)+ ifelse(f2_exp6$tube == "38", 8 ,0) + ifelse(f2_exp6$tube == "39", 9 ,0)+ ifelse(f2_exp6$tube == "40", 10 ,0)
f2_exp6_comp <- data.frame(gr_scaled = f2_exp6$gr_scaled, F1 = f2_exp6$F1, env = f2_exp6$env, condition = f2_exp6$cond, tube = f2_exp6$tube) 

F2.m.exp6 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = f2_exp6_comp , warmup = 1000, iter = 4000, chains = 4, log_lik = T,
  control=list(adapt_delta=0.99))

precis(F2.m.exp6, depth = 2)
plot(precis(F2.m.exp6, depth=2))
traceplot(F2.m.exp6)
trankplot(F2.m.exp6)
post.f2.exp6 <- extract.samples(F2.m.exp6)
#If current environment was 1.5%
#What is the effect of F1 environment? in model f2 exp5
diff1_4_f2exp6 <- post.f2.exp6$a[,1] - post.f2.exp6$a[,4]
median(diff1_4_f2exp6)
HPDI(diff1_4_f2exp6, prob = 0.95)#0.4938374 1.6958173
#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3_f2exp6 <- post.f2.exp6$a[,2] - post.f2.exp6$a[,3]
median(diff2_3_f2exp6)
HPDI(diff2_3_f2exp6, prob = 0.95)#-0.2791550  0.9185197
#####exp7f2####################################
f2_exp7 <- filter(f2data1, Experiment == "Exp7")
f2_exp7_comp <- data.frame(gr_scaled = f2_exp7$gr_scaled, condition = f2_exp7$cond) 

F2.m.exp7 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = f2_exp7_comp , warmup = 2000, iter = 4000, chains = 4, log_lik = T,
  control=list(adapt_delta=0.99))

precis(F2.m.exp7, depth = 2)
plot(precis(F2.m.exp7, depth=2))
traceplot(F2.m.exp7)
trankplot(F2.m.exp7)
post.f2.exp7 <- extract.samples(F2.m.exp7)
#If current environment was 1.5%
#What is the effect of F1 environment? in model f2 exp5
diff1_4_f2exp7 <- post.f2.exp7$a[,1] - post.f2.exp7$a[,4]
median(diff1_4_f2exp7)
HPDI(diff1_4_f2exp7, prob = 0.95) #-0.1890998  1.1196438
#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3_f2exp7 <- post.f2.exp7$a[,2] - post.f2.exp7$a[,3]
median(diff2_3_f2exp7)
HPDI(diff2_3_f2exp7, prob = 0.95)#1.193030 2.532473 
#####exp8f2###################################
f2_exp8 <- filter(f2data1, Experiment == "Exp8")
f2_exp8_comp <- data.frame(gr_scaled = f2_exp8$gr_scaled, condition = f2_exp8$cond) 

F2.m.exp8 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = f2_exp8_comp , warmup = 1000, iter = 4000, chains = 4, log_lik = T,
  control=list(adapt_delta=0.99))

precis(F2.m.exp8, depth = 2)
plot(precis(F2.m.exp8, depth=2))
traceplot(F2.m.exp8)
trankplot(F2.m.exp8)
post.f2.exp8 <- extract.samples(F2.m.exp8)
#If current environment was 1.5%
#What is the effect of F1 environment? in model f2 exp5
diff1_4_f2exp8 <- post.f2.exp8$a[,1] - post.f2.exp8$a[,4]
median(diff1_4_f2exp8)
HPDI(diff1_4_f2exp8, prob = 0.95)#0.05052982 1.70380471
#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3_f2exp8 <- post.f2.exp8$a[,2] - post.f2.exp8$a[,3]
median(diff2_3_f2exp8)
HPDI(diff2_3_f2exp8, prob = 0.95)#-0.7953857  0.8399137
#####exp9f2#####################################
f2_exp9 <- filter(f2data1, Experiment == "Exp9")
f2_exp9_comp <- data.frame(gr_scaled = f2_exp9$gr_scaled, condition = f2_exp9$cond) 

F2.m.exp9 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = f2_exp9_comp , warmup = 1000, iter = 4000, chains = 4, log_lik = T,
  control=list(adapt_delta=0.99))

precis(F2.m.exp9, depth = 2)
plot(precis(F2.m.exp9, depth=2))
traceplot(F2.m.exp9)
trankplot(F2.m.exp9)
post.f2.exp9 <- extract.samples(F2.m.exp9)
#If current environment was 1.5%
#What is the effect of F1 environment? in model f2 exp5
diff1_4_f2exp9 <- post.f2.exp9$a[,1] - post.f2.exp9$a[,4]
median(diff1_4_f2exp9)
HPDI(diff1_4_f2exp9, prob = 0.95) #1.038110 1.992098
#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3_f2exp9 <- post.f2.exp9$a[,2] - post.f2.exp9$a[,3]
median(diff2_3_f2exp9)
HPDI(diff2_3_f2exp9, prob = 0.95) #0.3968963 1.3744645
#create the database for the posterior graphs
f2.plot.post1_4 <- data.frame(posterior = c(diff1_4_model1, diff1_4_f2exp1, diff1_4_f2exp2, diff1_4_f2exp3, diff1_4_f2exp4, diff1_4_f2exp5, diff1_4_f2exp6,
                                            diff1_4_f2exp7, diff1_4_f2exp8, diff1_4_f2exp9), experiment = factor(c(rep("Combined", length(diff1_4_model1)),
                                                                                                                   rep("Experiment 1", length(diff1_4_f2exp1)), rep("Experiment 2", length(diff1_4_f2exp2)), rep("Experiment 3", length(diff1_4_f2exp3)),
                                                                                                                   rep("Experiment 4", length(diff1_4_f2exp4)), rep("Experiment 5", length(diff1_4_f2exp5)), rep("Experiment 6", length(diff1_4_f2exp6)),
                                                                                                                   rep("Experiment 7", length(diff1_4_f2exp7)), rep("Experiment 8", length(diff1_4_f2exp8)), rep("Experiment 9", length(diff1_4_f2exp9)))))


f2.plot.post2_3 <- data.frame(posterior = c(diff2_3_model1, diff2_3_f2exp1, diff2_3_f2exp2, diff2_3_f2exp3, diff2_3_f2exp4, diff2_3_f2exp5, diff2_3_f2exp6,
                                            diff2_3_f2exp7, diff2_3_f2exp8, diff2_3_f2exp9), experiment = factor(c(rep("Combined", length(diff2_3_model1)),
                                                                                                                   rep("Experiment 1", length(diff2_3_f2exp1)), rep("Experiment 2", length(diff2_3_f2exp2)), rep("Experiment 3", length(diff2_3_f2exp3)),
                                                                                                                   rep("Experiment 4", length(diff2_3_f2exp4)), rep("Experiment 5", length(diff2_3_f2exp5)), rep("Experiment 6", length(diff2_3_f2exp6)),
                                                                                                                   rep("Experiment 7", length(diff2_3_f2exp7)), rep("Experiment 8", length(diff2_3_f2exp8)), rep("Experiment 9", length(diff2_3_f2exp9)))))



f2.plot.post.2_3$F2 <- factor(c(rep("0.015%", length(f2.plot.post.2_3$posterior))))
f2.plot.post.1_4$F2 <- factor(c(rep("1.5%", length(f2.plot.post.1_4$posterior))))
f2.plot.post <- rbind(f2.plot.post.1_4, f2.plot.post.2_3)

#Poterior graphs
plot.f2 <- f2.plot.post %>%
  ggplot(aes(y = experiment, x = posterior)) +
  stat_halfeye(.width = c(.95, .5)) +
  facet_wrap(~F2)+
  geom_vline(xintercept = c(0), linetype = "dashed") +
  theme(legend.position = "none",
        text = element_text(size = 20),
        axis.title.y  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  xlab(expression(paste("Effect of F"[1], " 1.5 % sucrose")))

###############################################################################
###############################################################################
############## Test the F3 effect ##################

exp1f3 <- filter(f3data, Experiment == "Exp1")
exp2f3 <- filter(f3data, Experiment == "Exp2")
exp3f3 <- filter(f3data, Experiment == "Exp3")
exp5f3 <- filter(f3data, Experiment == "Exp5")
exp6f3 <- filter(f3data, Experiment == "Exp6")

##Scaling the data
exp1f3$gr_scaled <- scale(exp1f3$Time_1)
exp1f3$col_scaled <- scale(exp1f3$mean_colonies)
exp2f3$gr_scaled <- scale(exp2f3$Time_1)
exp2f3$col_scaled <- scale(exp2f3$mean_colonies)
exp3f3$gr_scaled <- scale(exp3f3$Time_1)
exp3f3$col_scaled <- scale(exp3f3$mean_colonies)
exp5f3$gr_scaled <- scale(exp5f3$Time_1)
exp5f3$col_scaled <- scale(exp5f3$mean_colonies)
exp6f3$gr_scaled <- scale(exp6f3$Time_1)
exp6f3$col_scaled <- scale(exp6f3$mean_colonies)

f3data1 <- rbind(exp1f3, exp2f3, exp3f3, exp5f3, exp6f3)
f3data1$F1 <- ifelse(f3data1$f1_sucrose == 1.5, 1, 2)
f3data1$env <- ifelse(f3data1$p_sucrose == 1.5, 1, 2)
f3data1$condition <- factor(ifelse(f3data1$condition == "C1", 1 ,0) + ifelse(f3data$condition == "C2", 2, 0) + ifelse(f3data$condition == "C3", 3, 0) + ifelse(f3data$condition == "C4", 4, 0))
f3data1$tube <- factor(f3data1$tube)
f3datcomp <- data.frame(gr_scaled = f3data1$gr_scaled, col_scaled = f3data1$col_scaled, F1 = f3data1$F1, env = f3data1$env, condition = f3data1$condition, tube = f3data1$tube)

ggplot(f3data1, aes(x = factor(f1_sucrose), y = col_scaled, fill = factor(f1_sucrose))) + geom_boxplot()
ggplot(f3data1, aes(x = condition, y = col_scaled)) + geom_boxplot()
ggplot(f3data1, aes(x = tube, y = col_scaled)) + geom_boxplot()

modelf3 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition] + B[tube] + Bc*col_scaled,
    a[condition] ~ dnorm(0,1),
    Bc ~ dnorm(0,1),
    B[tube] ~ dnorm(0, sigmaB),
    sigmaB ~ dexp(1),
    sigma ~ dexp(1)
  ), data = f3datcomp, warmup = 5000, iter = 15000, 
  thin = 5, chains = 4, cores = 2, log_lik = T, 
  control = list(adapt_delta = 0.9999, max_treedepth = 13))


p3 <- precis(modelf3, depth = 2)
plot(precis(modelf3, depth=2))
traceplot(modelf3)
trankplot(modelf3)
postf3 <- extract.samples(modelf3, n = 3000) #Need to adjust n here because of thinning
#If current environment was 1.5%
#What is the effect of F1 environment?
diff1_4 <- postf3$a[,1] - postf3$a[,4]
mean(diff1_4) #0.1218908
HPDI(diff1_4, prob = 0.95) #-0.1228629  0.3903956 
#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3 <- postf3$a[,2] - postf3$a[,3]
mean(diff2_3) #-0.2403022
HPDI(diff2_3, prob = 0.95) #-0.49935454  0.01137214 

#graph of the model results
b1 <- postf3$a[,1] 
b2 <- postf3$a[,2]
b3 <- postf3$a[,3]
b4 <- postf3$a[,4]
############### F3 generation each experiment separated ##################
#exp1f3
f3_exp1 <- filter(f3data1, Experiment == "Exp1")
f3_exp1_comp <- data.frame(gr_scaled = f3_exp1$gr_scaled, condition = f3_exp1$cond, tube = f3_exp1$tube) 

F3.m.exp1 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],# + B[tube],
    a[condition] ~ dnorm(0,1),
    #B[tube] ~ dnorm(0, sigmaB),
    #sigmaB ~ dexp(1),
    sigma ~ dexp(1)
  ), data = f3_exp1_comp , warmup = 1000, iter = 4000, chains = 4, log_lik = T,       
  control=list(adapt_delta=0.99))

precis(F3.m.exp1, depth = 2)
plot(precis(F3.m.exp1, depth=2))
traceplot(F3.m.exp1)
trankplot(F3.m.exp1)
post.f3.exp1 <- extract.samples(F3.m.exp1)
#If current environment was 1.5%
#What is the effect of F1 environment? 1.5%
diff1_4_f3exp1 <- post.f3.exp1$a[,1] - post.f3.exp1$a[,4]
median(diff1_4_f3exp1)
HPDI(diff1_4_f3exp1, prob = 0.95) #-0.05161484  1.08279110
#What is the effect of F1 environment? 0.015%
diff2_3_f3exp1 <- post.f3.exp1$a[,2] - post.f3.exp1$a[,3]
median(diff2_3_f3exp1)
HPDI(diff2_3_f3exp1, prob = 0.95) #-1.3596416 -0.2233301
####exp2f3####################################################
f3_exp2 <- filter(f3data1, Experiment == "Exp2")
f3_exp2$tube <- ifelse(f3_exp2$tube == "7", 1 ,0)+ ifelse(f3_exp2$tube == "8", 2 ,0) + ifelse(f3_exp2$tube == "9", 3 ,0)+ 
  ifelse(f3_exp2$tube == "10", 4 ,0)+ ifelse(f3_exp2$tube == "11", 5 ,0) + ifelse(f3_exp2$tube == "12", 6 ,0)
f3_exp2_comp <- data.frame(gr_scaled = f3_exp2$gr_scaled, condition = f3_exp2$cond, tube = f3_exp2$tube) 

F3.m.exp2 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],# + B[tube],
    a[condition] ~ dnorm(0,1),
    #B[tube] ~ dnorm(0, sigmaB),
    #sigmaB ~ dexp(1),
    sigma ~ dexp(1)
  ), data = f3_exp2_comp , warmup = 1000, iter = 4000, chains = 4, log_lik = T,       
  control=list(adapt_delta=0.99))

precis(F3.m.exp2, depth = 2)
plot(precis(F3.m.exp2, depth=2))
traceplot(F3.m.exp2)
trankplot(F3.m.exp2)
post.f3.exp2 <- extract.samples(F3.m.exp2)
#If current environment was 1.5%
#What is the effect of F1 environment? 1.5%
diff1_4_f3exp2 <- post.f3.exp2$a[,1] - post.f3.exp2$a[,4]
median(diff1_4_f3exp2)
HPDI(diff1_4_f3exp2, prob = 0.95) #-0.5887360  0.8093219
#What is the effect of F1 environment? 0.015%
diff2_3_f3exp2 <- post.f3.exp2$a[,2] - post.f3.exp2$a[,3]
median(diff2_3_f3exp2)
HPDI(diff2_3_f3exp2, prob = 0.95) #-1.5098034 -0.1295857
####exp3f3###################################################

f3_exp3 <- filter(f3data1, Experiment == "Exp3")
f3_exp3$tube <- ifelse(f3_exp3$tube == "13", 1 ,0)+ ifelse(f3_exp3$tube == "14", 2 ,0) + ifelse(f3_exp3$tube == "15", 3 ,0)+ 
  ifelse(f3_exp3$tube == "16", 4 ,0)+ ifelse(f3_exp3$tube == "17", 5 ,0) + ifelse(f3_exp3$tube == "18", 6 ,0)
f3_exp3_comp <- data.frame(gr_scaled = f3_exp3$gr_scaled, condition = f3_exp3$cond, tube = f3_exp3$tube) 

F3.m.exp3 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],# + B[tube],
    a[condition] ~ dnorm(0,1),
    #B[tube] ~ dnorm(0, sigmaB),
    #sigmaB ~ dexp(1),
    sigma ~ dexp(1)
  ), data = f3_exp3_comp , warmup = 1000, iter = 4000, chains = 4, log_lik = T,       
  control=list(adapt_delta=0.99))

precis(F3.m.exp3, depth = 2)
plot(precis(F3.m.exp3, depth=2))
traceplot(F3.m.exp3)
trankplot(F3.m.exp3)
post.f3.exp3 <- extract.samples(F3.m.exp3)
#If current environment was 1.5%
#What is the effect of F1 environment? 1.5%
diff1_4_f3exp3 <- post.f3.exp3$a[,1] - post.f3.exp3$a[,4]
median(diff1_4_f3exp3)
HPDI(diff1_4_f3exp3, prob = 0.95) #-0.0881425  0.8791411
#What is the effect of F1 environment? 0.015%
diff2_3_f3exp3 <- post.f3.exp3$a[,2] - post.f3.exp3$a[,3]
median(diff2_3_f3exp3)
HPDI(diff2_3_f3exp3, prob = 0.95) #-0.5991300  0.3673352
####exp5f3####################################################
f3_exp5 <- filter(f3data1, Experiment == "Exp5")
f3_exp5$tube <- ifelse(f3_exp5$tube == "25", 1 ,0)+ ifelse(f3_exp5$tube == "26", 2 ,0) + ifelse(f3_exp5$tube == "27", 3 ,0)+ 
  ifelse(f3_exp5$tube == "28", 4 ,0)+ ifelse(f3_exp5$tube == "29", 5 ,0) + ifelse(f3_exp5$tube == "30", 6 ,0)
f3_exp5_comp <- data.frame(gr_scaled = f3_exp5$gr_scaled, condition = f3_exp5$cond, tube = f3_exp5$tube) 

F3.m.exp5 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],# + B[tube],
    a[condition] ~ dnorm(0,1),
    #B[tube] ~ dnorm(0, sigmaB),
    #sigmaB ~ dexp(1),
    sigma ~ dexp(1)
  ), data = f3_exp5_comp , warmup = 1000, iter = 4000, chains = 4, log_lik = T,       
  control=list(adapt_delta=0.99))

precis(F3.m.exp5, depth = 2)
plot(precis(F3.m.exp5, depth=2))
traceplot(F3.m.exp5)
trankplot(F3.m.exp5)
post.f3.exp5 <- extract.samples(F3.m.exp5)
#If current environment was 1.5%
#What is the effect of F1 environment? 1.5%
diff1_4_f3exp5 <- post.f3.exp5$a[,1] - post.f3.exp5$a[,4]
median(diff1_4_f3exp5)
HPDI(diff1_4_f3exp5, prob = 0.95) #-0.6518979  0.7756466
#What is the effect of F1 environment? 0.015%
diff2_3_f3exp5 <- post.f3.exp5$a[,2] - post.f3.exp5$a[,3]
median(diff2_3_f3exp5)
HPDI(diff2_3_f3exp5, prob = 0.95) #-0.5390304  0.8528448 
####exp6f3####################################################
f3_exp6 <- filter(f3data1, Experiment == "Exp6")
f3_exp6$tube <- ifelse(f3_exp6$tube == "31", 1 ,0)+ ifelse(f3_exp6$tube == "32", 2 ,0) + ifelse(f3_exp6$tube == "33", 3 ,0)+ 
  ifelse(f3_exp6$tube == "34", 4 ,0)+ ifelse(f3_exp6$tube == "35", 5 ,0) + ifelse(f3_exp6$tube == "36", 6 ,0)+
  ifelse(f3_exp6$tube == "37", 7, 0)+ ifelse(f3_exp6$tube == "38", 8, 0) + ifelse(f3_exp6$tube == "39", 9, 0)+
  ifelse(f3_exp6$tube == "40", 10, 0)
f3_exp6_comp <- data.frame(gr_scaled = f3_exp6$gr_scaled, condition = f3_exp6$cond, tube = f3_exp6$tube) 

F3.m.exp6 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],# + B[tube],
    a[condition] ~ dnorm(0,1),
    #B[tube] ~ dnorm(0, sigmaB),
    #sigmaB ~ dexp(1),
    sigma ~ dexp(1)
  ), data = f3_exp6_comp , warmup = 1000, iter = 4000, chains = 4, log_lik = T,       
  control=list(adapt_delta=0.99))

precis(F3.m.exp6, depth = 2)
plot(precis(F3.m.exp6, depth=2))
traceplot(F3.m.exp6)
trankplot(F3.m.exp6)
post.f3.exp6 <- extract.samples(F3.m.exp6)
#If current environment was 1.5%
#What is the effect of F1 environment? 1.5%
diff1_4_f3exp6 <- post.f3.exp6$a[,1] - post.f3.exp6$a[,4]
median(diff1_4_f3exp6)
HPDI(diff1_4_f3exp6, prob = 0.95) #-0.5792380  0.2482834
#What is the effect of F1 environment? 0.015%
diff2_3_f3exp6 <- post.f3.exp6$a[,2] - post.f3.exp6$a[,3]
median(diff2_3_f3exp6)
HPDI(diff2_3_f3exp6, prob = 0.95) #-0.3864450  0.4440107 

#create the database for the f3 posterior graphs
f3.plot.post1_4 <- data.frame(posterior = c(diff1_4, diff1_4_f3exp1, diff1_4_f3exp2, diff1_4_f3exp3, diff1_4_f3exp5, diff1_4_f3exp6),
                              experiment = factor(c(rep("Combined", length(diff1_4)),
                                                    rep("Experiment 1", length(diff1_4_f3exp1)), rep("Experiment 2", length(diff1_4_f3exp2)), rep("Experiment 3", length(diff1_4_f3exp3)),
                                                    rep("Experiment 4", length(diff1_4_f3exp5)), rep("Experiment 5", length(diff1_4_f3exp6)))))


f3.plot.post2_3 <- data.frame(posterior = c(diff2_3, diff2_3_f3exp1, diff2_3_f3exp2, diff2_3_f3exp3, diff2_3_f3exp5, diff2_3_f3exp6), experiment = factor(c(rep("Combined", length(diff2_3)),
                                                                                                                                                            rep("Experiment 1", length(diff2_3_f3exp1)), rep("Experiment 2", length(diff2_3_f3exp2)), rep("Experiment 3", length(diff2_3_f3exp3)),
                                                                                                                                                            rep("Experiment 4", length(diff2_3_f3exp5)), rep("Experiment 5", length(diff2_3_f3exp6)))))


f3.plot.post.1_4$F3 <- factor(c(rep("1.5%", length(f3.plot.post.1_4$posterior))))
f3.plot.post.2_3$F3 <- factor(c(rep("0.015%", length(f3.plot.post.2_3$posterior))))
f3.plot.post <- rbind(f3.plot.post.1_4, f3.plot.post.2_3)

plot.f3 <- f3.plot.post %>%
  ggplot(aes(y = experiment, x = posterior)) +
  stat_halfeye(.width = c(.95, .5)) +
  facet_wrap(~F3)+
  geom_vline(xintercept = c(0), linetype = "dashed") +
  theme(legend.position = "none",
        text = element_text(size = 20),
        axis.title.y  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  xlab(expression(paste("Effect of F"[1], " 1.5 % sucrose")))

####Viability and size models####################################################

aineisto.1 <- filter(aineisto, replicate == "1")

#Create datasets to run ulam 
aineisto.1$F1 <- ifelse(aineisto.1$f1_sucrose == 1.5, 1, 2)
aineisto.1$env <- ifelse(aineisto.1$p_sucrose == 1.5, 1, 2)
aineisto.1$condition <- ifelse(aineisto.1$condition == "C1", 1 ,0) + 
  ifelse(aineisto.1$condition == "C2", 2, 0) + 
  ifelse(aineisto.1$condition == "C3", 3, 0) + 
  ifelse(aineisto.1$condition == "C4", 4, 0)

f2data3 <- filter(aineisto.1, generation == "F2")
f3data3 <- filter(aineisto.1, generation == "F3")
#scale the data in each generation
f2data3$col_scaled <- scale(f2data3$mean_colonies) 
f2data3$diam_scaled <- scale(f2data3$spore_diam)
f3data3$col_scaled <- scale(f3data3$mean_colonies) 
f3data3$diam_scaled <- scale(f3data3$spore_diam)

#graph viability data
ggplot(f2data3, aes(x = factor(f1_sucrose), y = mean_colonies)) + geom_boxplot()
ggplot(f3data3, aes(x = factor(p_sucrose), y = mean_colonies, fill = factor(f1_sucrose))) + geom_boxplot()
ggplot(f3data3, aes(x = factor(condition), y = mean_colonies)) + geom_boxplot()
ggplot(f3data3, aes(x = factor(condition), y = col_scaled)) + geom_boxplot()

#viavility model generation F2
f2data4 <- f2data3[complete.cases(f2data3$mean_colonies),]#eliminate NA
f2datcomp_viab <- data.frame(mean_colonies = f2data4$mean_colonies, col_scaled= as.numeric(f2data4$col_scaled), F1 = factor(f2data4$F1), env = f2data4$env, condition = factor(f2data4$condition))

f2.viab.m2 <- ulam(
  alist(
    col_scaled ~ dnorm(mu, sigma),
    mu <- a[F1],
    a[F1] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = f2datcomp_viab, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

precis(f2.viab.m2, depth = 2)
plot(precis(f2.viab.m2, depth=2))
traceplot(f2.viab.m2)
trankplot(f2.viab.m2)
post.f2.viab <- extract.samples(f2.viab.m2)
diff2_1_f2_viab <- post.f2.viab$a[,2] - post.f2.viab$a[,1]
mean(diff2_1_f2_viab)
HPDI(diff2_1_f2_viab, prob = 0.95)#-0.2038233  0.7354454 
#Viability model generation F3
f3datcomp_viab <- data.frame(mean_colonies = f3data3$mean_colonies, col_scaled= as.numeric(f3data3$col_scaled), F1 = factor(f3data3$F1), env = f3data3$env, condition = factor(f3data3$condition))

f3.viab.m1 <- ulam(
  alist(
    col_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = f3datcomp_viab, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

precis(f3.viab.m1, depth = 2)
plot(precis(f3.viab.m1, depth=2))
traceplot(f3.viab.m1)
trankplot(f3.viab.m1)
post.f3.viab <- extract.samples(f3.viab.m1)#WAIC in both models are smaller when using scaled colonies
diff4_1_f3_viab <- post.f3.viab$a[,4] - post.f3.viab$a[,1]
mean(diff4_1_f3_viab)
HPDI(diff4_1_f3_viab, prob = 0.95)#-0.7230396  0.5307407  
diff2_3_f3_viab <- post.f3.exp1$a[,2] - post.f3.exp1$a[,3]
median(diff2_3_f3_viab)
HPDI(diff2_3_f3_viab, prob = 0.95)#-0.7907413  0.4577521 
#####Size models####################################################

ggplot(f2data3, aes(x = factor(f1_sucrose), y = diam_scaled)) + geom_boxplot()
ggplot(f3data3, aes(x = factor(condition), y = spore_diam, fill = factor(f1_sucrose))) + geom_boxplot()
ggplot(f3data3, aes(x = factor(condition), y = diam_scaled, fill = factor(f1_sucrose))) + geom_boxplot()

#Size model generation f2
f2data6 <- f2data3[complete.cases(f2data3$spore_diam),]#eliminate NA
f2datcomp_diam <- data.frame(diam_scaled = as.numeric(f2data6$diam_scaled), diam = f2data6$spore_diam, F1 = factor(f2data6$F1), env = f2data6$env, condition = factor(f2data6$condition))

f2.diam.m1 <- ulam(
  alist(
    diam_scaled ~ dnorm(mu, sigma),
    mu <- a[F1],
    a[F1] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = f2datcomp_diam, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

precis(f2.diam.m1, depth = 2)
plot(precis(f2.diam.m1, depth=2))
traceplot(f2.diam.m1)
trankplot(f2.diam.m1)
post.f2.diam <- extract.samples(f2.diam.m1)
diff2_1_f2_diam <- post.f2.diam$a[,2] - post.f2.diam$a[,1]
mean(diff2_1_f2_diam)#-0.01915387
HPDI(diff2_1_f2_diam, prob = 0.95)#-0.3722206  0.3055794
#Size model generation f3
f3datcomp_diam <- data.frame(diam_scaled = as.numeric(f3data3$diam_scaled), diam = f3data3$spore_diam, F1 = factor(f3data3$F1), env = f3data3$env, condition = factor(f3data3$condition))

f3.diam.m1 <- ulam(
  alist(
    diam_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,5),
    sigma ~ dexp(5)
  ), data = f3datcomp_diam, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

precis(f3.diam.m1, depth = 2)
plot(precis(f3.diam.m1, depth=2))
traceplot(f3.diam.m1)
trankplot(f3.diam.m1)
post.f3.diam <- extract.samples(f3.diam.m1)
#What is the effect of F1 environment?
diff1_4 <- post.f3.diam$a[,1] - post.f3.diam$a[,4]
mean(diff1_4)
HPDI(diff1_4, prob = 0.95) #-1.0698132  0.2731048
#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3 <- post.f3.diam$a[,2] - post.f3.diam$a[,3]
median(diff2_3)
HPDI(diff2_3, prob = 0.95) #-0.9336272  0.4150267 
#What is the effect of the F2 env
diff3_4 <- post.f3.diam$a[,3] - post.f3.diam$a[,4]
median(diff3_4)
HPDI(diff3_4, prob = 0.95) #-0.9216284  0.3994920
#######Load viability non mean data###########################
sorbose <- read.csv("sorbose_viability.csv", header = T)
head(sorbose)

sorbose_f2 <- filter(sorbose, generation == "F2")
sorbose_f3 <- filter(sorbose, generation == "F3")

ggplot(sorbose_f2, aes(y= number_colonies, x= factor(g_sucrose)))+
  geom_boxplot()

ggplot(sorbose_f3, aes(y= number_colonies, x= condition))+
  geom_boxplot()

sorbose_f2$col_scaled <- scale(sorbose_f2$number_colonies)
sorbose_f3$col_scaled <- scale(sorbose_f3$number_colonies)

sorbose_f2$g_sucrose <- ifelse(sorbose_f2$g_sucrose == 1.5, 1, 2)
sorbose_f3$g_sucrose <- ifelse(sorbose_f3$g_sucrose == 1.5, 1, 2)

sorbose_f3$condition <- ifelse(sorbose_f3$condition == "C1", 1 ,0) + 
  ifelse(sorbose_f3$condition == "C2", 2, 0) + 
  ifelse(sorbose_f3$condition == "C3", 3, 0) + 
  ifelse(sorbose_f3$condition == "C4", 4, 0)


sorbose_f2$tube <- ifelse(sorbose_f2$Tube == "1", 1 ,0)+ ifelse(sorbose_f2$Tube == "2", 2 ,0) + ifelse(sorbose_f2$Tube == "3", 3 ,0)+ 
  ifelse(sorbose_f2$Tube == "4", 4 ,0)+ ifelse(sorbose_f2$Tube == "5", 5 ,0) + ifelse(sorbose_f2$Tube == "6", 6 ,0)+
  ifelse(sorbose_f2$Tube == "19", 7 ,0)+ ifelse(sorbose_f2$Tube == "20", 8 ,0) + ifelse(sorbose_f2$Tube == "21", 9 ,0)+ 
  ifelse(sorbose_f2$Tube == "22", 10 ,0)+ ifelse(sorbose_f2$Tube == "23", 11 ,0) + ifelse(sorbose_f2$Tube == "24", 12 ,0)+
  ifelse(sorbose_f2$Tube == "37", 13 ,0)+ ifelse(sorbose_f2$Tube == "38", 14 ,0) + ifelse(sorbose_f2$Tube == "39", 15 ,0)+ 
  ifelse(sorbose_f2$Tube == "40", 16 ,0)+ ifelse(sorbose_f2$Tube == "41", 17 ,0) + ifelse(sorbose_f2$Tube == "42", 18 ,0)+
  ifelse(sorbose_f2$Tube == "55", 19 ,0)+ ifelse(sorbose_f2$Tube == "56", 20 ,0) + ifelse(sorbose_f2$Tube == "57", 21 ,0)+ 
  ifelse(sorbose_f2$Tube == "58", 22 ,0)+ ifelse(sorbose_f2$Tube == "59", 23 ,0) + ifelse(sorbose_f2$Tube == "60", 24 ,0)+
  ifelse(sorbose_f2$Tube == "85", 25 ,0)+ ifelse(sorbose_f2$Tube == "86", 26 ,0) + ifelse(sorbose_f2$Tube == "87", 27 ,0)+ 
  ifelse(sorbose_f2$Tube == "88", 28 ,0)+ ifelse(sorbose_f2$Tube == "89", 29 ,0) + ifelse(sorbose_f2$Tube == "90", 30 ,0)+
  ifelse(sorbose_f2$Tube == "91", 31 ,0)+ ifelse(sorbose_f2$Tube == "92", 32 ,0) + ifelse(sorbose_f2$Tube == "93", 33 ,0)+
  ifelse(sorbose_f2$Tube == "94", 34 ,0)

sorbosef2_comp <- data.frame(col_scaled = as.numeric(sorbose_f2$col_scaled), #col = sorbose_f2$number_colonies, 
                             F1 = factor(sorbose_f2$g_sucrose), tube = factor(sorbose_f2$tube))

sorbose_f3$tube <- ifelse(sorbose_f3$Tube == "7", 1 ,0)+ ifelse(sorbose_f3$Tube == "8", 2 ,0) + ifelse(sorbose_f3$Tube == "9", 3 ,0)+ 
  ifelse(sorbose_f3$Tube == "10", 4 ,0)+ ifelse(sorbose_f3$Tube == "11", 5 ,0) + ifelse(sorbose_f3$Tube == "12", 6 ,0)+
  ifelse(sorbose_f3$Tube == "13", 7 ,0)+ ifelse(sorbose_f3$Tube == "14", 8 ,0) + ifelse(sorbose_f3$Tube == "15", 9 ,0)+ 
  ifelse(sorbose_f3$Tube == "16", 10 ,0)+ ifelse(sorbose_f3$Tube == "17", 11 ,0) + ifelse(sorbose_f3$Tube == "18", 12 ,0)+
  ifelse(sorbose_f3$Tube == "25", 13 ,0)+ ifelse(sorbose_f3$Tube == "26", 14 ,0) + ifelse(sorbose_f3$Tube == "27", 15 ,0)+ 
  ifelse(sorbose_f3$Tube == "28", 16 ,0)+ ifelse(sorbose_f3$Tube == "29", 17 ,0) + ifelse(sorbose_f3$Tube == "30", 18 ,0)+
  ifelse(sorbose_f3$Tube == "31", 19 ,0)+ ifelse(sorbose_f3$Tube == "32", 20 ,0) + ifelse(sorbose_f3$Tube == "33", 21 ,0)+ 
  ifelse(sorbose_f3$Tube == "34", 22 ,0)+ ifelse(sorbose_f3$Tube == "35", 23 ,0) + ifelse(sorbose_f3$Tube == "36", 24 ,0)+
  ifelse(sorbose_f3$Tube == "43", 25 ,0)+ ifelse(sorbose_f3$Tube == "44", 26 ,0) + ifelse(sorbose_f3$Tube == "45", 27 ,0)+ 
  ifelse(sorbose_f3$Tube == "46", 28 ,0)+ ifelse(sorbose_f3$Tube == "47", 29 ,0) + ifelse(sorbose_f3$Tube == "48", 30 ,0)+
  ifelse(sorbose_f3$Tube == "49", 31 ,0)+ ifelse(sorbose_f3$Tube == "50", 32 ,0) + ifelse(sorbose_f3$Tube == "51", 33 ,0)+
  ifelse(sorbose_f3$Tube == "52", 34 ,0)+ ifelse(sorbose_f3$Tube == "53", 35 ,0)+ ifelse(sorbose_f3$Tube == "54", 36 ,0) + 
  ifelse(sorbose_f3$Tube == "61", 37 ,0)+ ifelse(sorbose_f3$Tube == "62", 38 ,0)+ ifelse(sorbose_f3$Tube == "63", 39 ,0) + 
  ifelse(sorbose_f3$Tube == "64", 40 ,0)+ ifelse(sorbose_f3$Tube == "65", 41 ,0)+ ifelse(sorbose_f3$Tube == "66", 42 ,0) + 
  ifelse(sorbose_f3$Tube == "67", 43 ,0)+ ifelse(sorbose_f3$Tube == "68", 44 ,0)+ ifelse(sorbose_f3$Tube == "69", 45 ,0)+
  ifelse(sorbose_f3$Tube == "70", 46 ,0)+ ifelse(sorbose_f3$Tube == "71", 47 ,0)+ ifelse(sorbose_f3$Tube == "72", 48 ,0)+
  ifelse(sorbose_f3$Tube == "95", 49 ,0)+ ifelse(sorbose_f3$Tube == "96", 50 ,0)+ ifelse(sorbose_f3$Tube == "97", 51 ,0)+
  ifelse(sorbose_f3$Tube == "98", 52 ,0)+ ifelse(sorbose_f3$Tube == "99", 53 ,0)+ ifelse(sorbose_f3$Tube == "100", 54 ,0)+
  ifelse(sorbose_f3$Tube == "101", 55 ,0)+ ifelse(sorbose_f3$Tube == "102", 56 ,0)+ ifelse(sorbose_f3$Tube == "103", 57 ,0)+
  ifelse(sorbose_f3$Tube == "104", 58 ,0)+ ifelse(sorbose_f3$Tube == "105", 59 ,0)+ ifelse(sorbose_f3$Tube == "106", 60 ,0)+
  ifelse(sorbose_f3$Tube == "107", 61 ,0)+ ifelse(sorbose_f3$Tube == "108", 62 ,0)+ ifelse(sorbose_f3$Tube == "109", 63 ,0)+
  ifelse(sorbose_f3$Tube == "110", 64 ,0)+ ifelse(sorbose_f3$Tube == "111", 65 ,0)+ ifelse(sorbose_f3$Tube == "112", 66 ,0)+
  ifelse(sorbose_f3$Tube == "113", 67 ,0)+ ifelse(sorbose_f3$Tube == "114", 68 ,0)


sorbosef3_comp <- data.frame(col_scaled = as.numeric(sorbose_f3$col_scaled), 
                             #col = sorbose_f3$number_colonies, 
                             condition = factor(sorbose_f3$condition), 
                             tube = factor(sorbose_f3$tube))

f2.sorb.m1 <- ulam(
  alist(
    col_scaled ~ dnorm(mu, sigma),
    mu <- a[F1] + B[tube],
    a[F1] ~ dnorm(0,1),
    B[tube]~ dnorm(0, sigmab),
    sigma ~ dexp(1),
    sigmab ~ dexp(1)
  ), data = sorbosef2_comp, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

ps.f2 <- precis(f2.sorb.m1, depth = 2)
plot(precis(f2.sorb.m1, depth=2))
traceplot(f2.sorb.m1)
trankplot(f2.sorb.m1)
post.f2.sorb <- extract.samples(f2.sorb.m1)
#What is the effect of F1 environment?
diff1_2 <- post.f2.sorb$a[,1] - post.f2.sorb$a[,2]
mean(diff1_2) #-0.1747101
HPDI(diff1_2, prob = 0.95) #-0.7128179  0.3796963 

f3.sorb.m1 <- ulam(
  alist(
    col_scaled ~ dnorm(mu, sigma),
    mu <- a[condition] + B[tube],
    a[condition] ~ dnorm(0,5),
    B[tube]~ dnorm(0, sigmab),
    sigma ~ dexp(5),
    sigmab ~ dexp(5)
  ), data = sorbosef3_comp, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

ps.f3 <-precis(f3.sorb.m1, depth = 2)
plot(precis(f3.sorb.m1, depth=2))
traceplot(f3.sorb.m1)
trankplot(f3.sorb.m1)
post.f3.sorb <- extract.samples(f3.sorb.m1)
diff4_1_f3_sorb <- post.f3.sorb$a[,4] - post.f3.sorb$a[,1]
mean(diff4_1_f3_sorb)#-0.2091961
HPDI(diff4_1_f3_sorb, prob = 0.95) #-0.8080216  0.4153285
diff2_3_f3_sorb <- post.f3.sorb$a[,2] - post.f3.sorb$a[,3]
median(diff2_3_f3_sorb)#-0.1212853
HPDI(diff2_3_f3_sorb, prob = 0.95)#-0.6716196  0.4194678 

save(ps.f2, post.f2.sorb, file = "f2.sorb.m1.RData")
save(ps.f2, post.f2.sorb, file = "f3.sorb.m1.RData")

#####Make graphs#####################################################################
#create table to add p values
df1 <- tibble::tribble(
  ~group1, ~group2, ~p.adj,  ~y.position, ~f1_sucrose,
  0.8, 1.2,  "[0.43, 1.30]",  24,  "0.015",
  1.8, 2.2, "[1.44, 2.30]",  28,   "1.5")

plotf2_time1 <- ggplot(exp1f2, aes(x = factor(p_sucrose), y = Time_1, fill = factor(f1_sucrose))) + 
  geom_boxplot(lwd=1, outlier.size = 3)+
  theme_classic()+
  theme(text = element_text(size = 20),
        legend.position = "none",
        axis.title.x  = element_blank())+
  #xlab("F2 environment (sucrose %)")+
  ylab("Colony size (mm)")+
  labs(fill = "F1 environment (sucrose %)")+
  scale_fill_manual(values = c("#56B4E9", "#FC4E07"))
#add_pvalue(df1,  label.size = 4, bracket.size = 1)


plotf2_viab <- ggplot(sorbose_f2, aes(x = factor(F1), y = number_colonies)) + 
  geom_boxplot(lwd=1, outlier.size = 3)+
  theme_classic()+
  theme(text = element_text(size = 20),
        axis.title.x  = element_blank())+
  #xlab("F1 environment (sucrose %)")+
  ylab("Number of colonies")

plotf2_diam <- ggplot(f2data3, aes(x = factor(f1_sucrose), y = spore_diam)) + 
  geom_boxplot(lwd=1, outlier.size = 3)+
  theme_classic()+
  theme(text = element_text(size = 20),
        axis.title.x  = element_blank())+
  #xlab("F1 environment (sucrose %)")+
  ylab(expression(paste("Spore diameter ", "(",mu, "m",")")))

model1.results <- data.frame(posterior = c(a3, a2, a4, a1), 
                             F1 = factor(c(rep("0.015", length(a3)), rep("1.5", length(a2)), 
                                           rep("0.015", length(a4)), rep("1.5", length(a1)))),
                             F2 = factor(c(rep("0.015", length(a3)), rep("0.015", length(a2)), 
                                           rep("1.5", length(a4)), rep("1.5", length(a1)))))

model1.rplot <- model1.results %>%
  ggplot(aes(y = posterior, x = F2)) +
  stat_pointinterval(aes(color = F1),position = position_dodge(width = .4)) +
  theme(legend.position = "none",
        text = element_text(size = 20),
        axis.title.x  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  #xlab(expression(paste("F"[2], " environment (sucrose%)")))+
  ylab("Scaled data, colony size (mm)")+
  scale_color_manual(values = c("#56B4E9", "#FC4E07"))+
  geom_signif(y_position = c(0.3), xmin = c(0.7), xmax = c(1.2),
              annotation = c("[0.43, 1.30]"), tip_length = 0.03)+
  geom_signif(y_position = c(1.2), xmin = c(1.7), xmax = c(2.2),
              annotation = c("[1.44, 2.30]"), tip_length = 0.03)


# plotting F3 graphs
plotf3_time1 <- ggplot(exp6f3, aes(x = factor(p_sucrose), y = Time_1, fill = factor(f1_sucrose))) + 
  geom_boxplot(lwd=1, outlier.size = 3)+
  theme_classic()+
  theme(text = element_text(size = 20),
        axis.title.x  = element_blank(),
        legend.position = "none")+
  xlab("F3 environment (sucrose %)")+
  ylab("Colony size (mm)")+
  scale_fill_manual(values = c("#56B4E9", "#FC4E07"))

plotf3_viab <- ggplot(sorbose_f3, aes(x = factor(F2), y = number_colonies, fill = factor(F1))) +
  geom_boxplot(lwd=1, outlier.size = 3)+
  theme_classic()+
  theme(text = element_text(size = 20),
        axis.title.x  = element_blank(),
        legend.position = "none")+
  xlab("F3 environment (sucrose %)")+
  ylab("Number of colonies")+
  scale_fill_manual(values = c("#56B4E9", "#FC4E07"))

plotf3_diam <- ggplot(f3data3, aes(x = factor(p_sucrose), y = spore_diam, fill = factor(f1_sucrose))) + 
  geom_boxplot(lwd=1, outlier.size = 3)+
  theme_classic()+
  theme(text = element_text(size = 20),
        axis.title.x  = element_blank(),
        legend.position = "none")+
  xlab("F3 environment (sucrose %)")+
  ylab(expression(paste("Spore diameter ", "(",mu, "m",")")))+
  scale_fill_manual(values = c("#56B4E9", "#FC4E07"))

modelf3.results <- data.frame(posterior = c(b3, b2, b4, b1), 
                              F1 = factor(c(rep("0.015", length(b3)), rep("1.5", length(b2)), 
                                            rep("0.015", length(b4)), rep("1.5", length(b1)))),
                              F2 = factor(c(rep("0.015", length(b3)), rep("0.015", length(b2)), 
                                            rep("1.5", length(b4)), rep("1.5", length(b1)))))


modelf3.rplot <- modelf3.results %>%
  ggplot(aes(y = posterior, x = F2, color= F1)) +
  stat_pointinterval(position = position_dodge(width = .4)) +
  theme(legend.position = "none",
        text = element_text(size = 19),
        axis.title.x  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  #xlab(expression(paste("F"[2], " environment (sucrose%)")))+
  ylab("Scaled data, colony size (mm)")+
  scale_color_manual(values = c("#56B4E9", "#FC4E07"))


#### Mutants ########################################
aineisto.suc4 <- read.csv("tgen_mutants_suc.csv", dec = ",", sep = ";", header = T)
levels(aineisto.suc4$F1) <- c("0.015%", "1.5%")
levels(aineisto.suc4$F2) <- c("0.015%", "1.5%")

aineisto.suc4$F1 <- relevel(aineisto.suc4$F1, "1.5%")
aineisto.suc4$F2 <- relevel(aineisto.suc4$F2, "1.5%")

aineisto.suc4 <- filter(aineisto.suc4, Genotype != 2489) #Remove 2489 because of incomplete data (there was a mix up with the plates, some experimental combinations are missing...
#Separate the  file by mutant
aineisto.suc4$F1_sucrose <- ifelse(aineisto.suc4$F1 == 1.5, 1, 2)
aineisto.suc4$F2_sucrose <- ifelse(aineisto.suc4$F2 == 1.5, 1, 2)
dim2 <- filter(aineisto.suc4, Genotype == "dim-2")
dim2$condition <- c(rep("1", 10), rep("2", 10), rep("4", 10), rep("3", 10))
dim2$gr_scaled <- scale(dim2$t2)
qde2 <- filter(aineisto.suc4, Genotype == "qde-2")
qde2$condition <- c(rep("1", 10), rep("2", 10), rep("4", 10), rep("3", 10))
qde2$gr_scaled <- scale(qde2$t2)
set7 <- filter(aineisto.suc4, Genotype == "set-7")
set7$condition <- c(rep("1", 10), rep("2", 10), rep("4", 10), rep("3", 10))
set7$gr_scaled <- scale(set7$t2)
dim2_comp <- data.frame(F1 = dim2$F1_sucrose, F2 = dim2$F2_sucrose, condition = factor(dim2$condition), gr_scaled = dim2$gr_scaled)

dim2.m1 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = dim2_comp, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

p.dim2 <- precis(dim2.m1, depth = 2)
plot(precis(dim2.m1, depth=2))
traceplot(dim2.m1)
trankplot(dim2.m1)
post.dim2 <- extract.samples(dim2.m1)

diff1_4_dim2 <- post.dim2$a[,1] - post.dim2$a[,4]
mean(diff1_4_dim2) #1.486374
HPDI(diff1_4_dim2, prob = 0.95) #1.125122 1.805499  
#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3_dim2 <- post.dim2$a[,2] - post.dim2$a[,3]
mean(diff2_3_dim2) #1.705876
HPDI(diff2_3_dim2, prob = 0.95) #1.374604 2.032721 

#qde2
qde2_comp <- data.frame(F1 = qde2$F1_sucrose, F2 = qde2$F2_sucrose, condition = factor(qde2$condition), gr_scaled = qde2$gr_scaled)

qde2.m1 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = qde2_comp, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

p.qde2 <- precis(qde2.m1, depth = 2)
plot(precis(qde2.m1, depth=2))
traceplot(qde2.m1)
trankplot(qde2.m1)
post.qde2 <- extract.samples(qde2.m1)

diff1_4_qde2 <- post.qde2$a[,1] - post.qde2$a[,4]
mean(diff1_4_qde2) #1.658764
HPDI(diff1_4_qde2, prob = 0.95) #1.401564 1.909671 
#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3_qde2 <- post.qde2$a[,2] - post.qde2$a[,3]
mean(diff2_3_qde2) #0.8253452
HPDI(diff2_3_qde2, prob = 0.95) #0.5766688 1.0773347

#set7
set7_comp <- data.frame(F1 = set7$F1_sucrose, F2 = set7$F2_sucrose, condition = factor(set7$condition), gr_scaled = set7$gr_scaled)

set7.m1 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = set7_comp, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

p.set7 <-precis(set7.m1, depth = 2)
plot(precis(set7.m1, depth=2))
traceplot(set7.m1)
trankplot(set7.m1)
post.set7 <- extract.samples(set7.m1)

diff1_4_set7 <- post.set7$a[,1] - post.set7$a[,4]
mean(diff1_4_set7) #1.14828
HPDI(diff1_4_set7, prob = 0.95) #0.7433659 1.5586473
#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3_set7 <- post.set7$a[,2] - post.set7$a[,3]
mean(diff2_3_set7) #1.664971
HPDI(diff2_3_set7, prob = 0.95) #1.251236 2.058906 


#mutant plot
df2 <- tibble::tribble(
  ~group1, ~group2, ~p.adj,  ~y.position, ~F1, ~Genotype,
  0.8, 1.2,  "[1.37, 2.03]",  40,  "0.015", "dim-2",
  1.8, 2.2,  "[1.12, 1.80]",  45,  "1.5",   "dim-2",
  0.8, 1.2,  "[0.57, 1.07]",  22,  "0.015", "qde-2",
  1.8, 2.2,  "[1.40, 1.90]",  29,  "1.5",   "qde-2",
  0.8, 1.2,  "[1.25, 2.05]",  33,  "0.015", "set-7",
  1.8, 2.2,  "[0.74, 1.55]",  35,  "1.5",   "set-7")

tgen_mutants_suc <- ggplot(aineisto.suc4, aes(x = F2, y = t2, fill = F1)) +
  geom_boxplot(lwd=0.5, outlier.size = 2)+
  ylim(10,46)+
  theme_linedraw()+
  theme(text = element_text(size = 15),
        legend.position = "top")+
  xlab(expression(paste("F"[2], " environment (sucrose %)")))+
  ylab("Colony size (mm)")+
  labs(fill = expression(paste("F"[1], " environment (sucrose %)")))+
  scale_fill_manual(values = c("#56B4E9", "#FC4E07"))+
  facet_wrap( ~Genotype, ncol = 3)+
  add_pvalue(df2,  label.size = 4,  bracket.size = 0.5)


######################### Growth rate ######################
#calculate grow rate
#(time2)-(time1)/time
#time was calculated independetly w the lab book notes. It is in min
exp1f2$gr_t1 <- scale((exp1f2$Time_1-exp1f2$Time_0)/ 1040)
exp2f2$gr_t1 <- scale((exp2f2$Time_1-exp2f2$Time_0)/ 1050)
exp3f2$gr_t1 <- scale((exp3f2$Time_1-exp3f2$Time_0)/ 1050)
exp4f2$gr_t1 <- scale((exp4f2$Time_1-exp4f2$Time_0)/ 1050)
exp5f2$gr_t1 <- scale((exp5f2$Time_1-exp5f2$Time_0)/ 1020)

exp1f2$gr_t2 <- scale((exp1f2$Time_2-exp1f2$Time_1)/ 240)
exp2f2$gr_t2 <- scale((exp2f2$Time_2-exp2f2$Time_1)/ 250)
exp3f2$gr_t2 <- scale((exp3f2$Time_2-exp3f2$Time_1)/ 250)
exp4f2$gr_t2 <- scale((exp4f2$Time_2-exp4f2$Time_1)/ 250)
exp5f2$gr_t2 <- scale((exp5f2$Time_2-exp5f2$Time_1)/ 150)

exp1f2$gr_t3 <- scale((exp1f2$Time_3-exp1f2$Time_2)/ 180)
exp2f2$gr_t3 <- scale((exp2f2$Time_3-exp2f2$Time_2)/ 160)
exp3f2$gr_t3 <- scale((exp3f2$Time_3-exp3f2$Time_2)/ 160)
exp4f2$gr_t3 <- scale((exp4f2$Time_3-exp4f2$Time_2)/ 180)
exp5f2$gr_t3 <- scale((exp5f2$Time_3-exp5f2$Time_2)/ 120)


gr_data <- rbind(exp1f2, exp2f2, exp3f2, exp4f2, exp5f2)

ggplot(gr_data, aes(y=sc_gr_t3, x=factor(p_sucrose), fill=factor(f1_sucrose)))+
  geom_boxplot()

gr_data$F1 <- ifelse(gr_data$f1_sucrose == 1.5, 1, 2)
gr_data$env <- ifelse(gr_data$p_sucrose == 1.5, 1, 2)
gr_data$condition <- ifelse(gr_data$condition == "C1", 1 ,0) + ifelse(gr_data$condition == "C2", 2, 0) + ifelse(gr_data$condition == "C3", 3, 0) + ifelse(gr_data$condition == "C4", 4, 0)
gr_data$tube <- factor(gr_data$tube)
#Create the database for rehtinking
gr_datcomp <- data.frame(gr_t1= gr_data$gr_t1, gr_t2= gr_data$gr_t2, gr_t3= gr_data$gr_t3,
                         F1 = gr_data$F1, env = gr_data$env, condition = factor(gr_data$condition), tube = gr_data$tube,
                         col_scaled = gr_data$col_scaled)

gr_datcomp <- gr_datcomp[complete.cases(gr_datcomp),]

m_t1_gr <- ulam(
  alist(
    gr_t1 ~ dnorm(mu, sigma),
    mu <- a[condition] + B[tube],
    a[condition] ~ dnorm(0,1),
    B[tube] ~ dnorm(0, sigmaB),
    sigmaB ~ dexp(1),
    sigma ~ dexp(1)
  ), data = gr_datcomp, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

p <- precis(m_t1_gr, depth = 2)
plot(precis(m_t1_gr, depth = 2))
traceplot(m_t1_gr)
trankplot(m_t1_gr)
post.gr1 <- extract.samples(m_t1_gr)
#If current environment was 1.5%
#What is the effect of F1 environment? in model 
diff1_4_gr1 <- post.gr1$a[,1] - post.gr1$a[,4]
mean(diff1_4_gr1) #1.168023
HPDI(diff1_4_gr1, prob = 0.95) #0.9127088 1.5514260
#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3_gr1 <- post.gr1$a[,2] - post.gr1$a[,3]
mean(diff2_3_gr1)
HPDI(diff2_3_gr1, prob = 0.95) #0.4506935 0.9709197 

#model2 with number of colonies (viability) as a fixed factor
m2_t1_gr <- ulam(
  alist(
    gr_t1 ~ dnorm(mu, sigma),
    mu <- a[condition] + B[tube] + Bc*col_scaled,
    a[condition] ~ dnorm(0,1),
    Bc ~ dnorm(0, 1),
    B[tube] ~ dnorm(0, sigmaB),
    sigmaB ~ dexp(1),
    sigma ~ dexp(1)
  ), data = gr_datcomp, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

precis(m2_t1_gr, depth = 2)
traceplot(m2_t1_gr)
trankplot(m2_t1_gr)
plot(precis(m2_t1_gr, depth=2))
post.m2gr1 <- extract.samples(m2_t1_gr)
median(post.m2gr1$Bc)
HPDI(post.m2gr1$Bc, prob = 0.95)
#If current environment was 1.5%
#What is the effect of F1 environment?
diff1_4 <- post.m2gr1$a[,1] - post.m2gr1$a[,4]
median(diff1_4)
HPDI(diff1_4, prob = 0.95)

#If current environment was 0.015%
#What is the effect of F1 environment?
diff2_3 <- post.m2gr1$a[,2] - post.m2gr1$a[,3]
median(diff2_3)
HPDI(diff2_3, prob = 0.95)
#########number of spores#########
conidia <- read.csv("number_conidia.csv", header = T)

conidia %>% group_by(factor(sucrose)) %>% 
  summarise(mean_size=mean(Concentration_conidia.ml))

ggplot(conidia, aes(y=Concentration_conidia.ml, x=sucrose))+
  geom_boxplot()

#Alternative carbon sources
carbons <- read.csv("tgen_carbons.csv", header = T)

ggplot(carbons, aes(F2, t1, fill = F1))+
  geom_boxplot()

ggplot(sucrose, aes(factor(F2), t14, fill = factor(F1)))+
  geom_boxplot()

ggplot(salt08, aes(F2, t1, fill = F1))+
  geom_boxplot()

ggplot(ph4, aes(F2, t21, fill = F1))+
  geom_boxplot()

ggplot(ph9, aes(F2, t1, fill = F1))+
  geom_boxplot()

######################## CELLULOSE #### 
cellulose.t <- c(1,2,3,4)
cellulose <- carbons%>% filter(condition %in% cellulose.t)
cellulose$gr_scaled <- scale(cellulose$t1)

cellulose.plot <- ggplot(cellulose, aes(F2, t1, fill = F1))+
  geom_boxplot()+ 
  theme_classic()+
  theme(text = element_text(size = 15),
        legend.position = "top",
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        legend.title = element_blank())+
  #xlab("F2 environment")+
  #ylab("Colony size (mm)")+
  labs(fill = "F1 environment")+
  scale_fill_manual(values = c("#DDA0DD", "#FC4E07"))


cellulose.m <- data.frame(gr_scaled = cellulose$gr_scaled, condition = factor(cellulose$condition))

cellulose_model1 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = cellulose.m, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

p <- precis(cellulose_model1, depth = 2)
plot(precis(cellulose_model1, depth = 2))
traceplot(cellulose_model1)
trankplot(cellulose_model1)
post.cellulose_m1 <- extract.samples(cellulose_model1)
#If current environment was sucrose
#What is the effect of F1 environment? in model 
diff1_4_cellulose_model1 <- post.cellulose_m1$a[,1] - post.cellulose_m1$a[,4]
mean(diff1_4_cellulose_model1 ) 
HPDI(diff1_4_cellulose_model1 , prob = 0.95) 
#If current environment was cellulose
#What is the effect of F1 environment?
diff2_3_cellulose_model1 <- post.cellulose_m1$a[,2] - post.cellulose_m1$a[,3]
mean(diff2_3_cellulose_model1)
HPDI(diff2_3_cellulose_model1, prob = 0.95) 

########################  ARABINOSE #########

arabinose.t <- c(1,5,6,7)
arabinose <- carbons %>% filter(condition %in% arabinose.t)
arabinose$gr_scaled <- scale(arabinose$t1)
arabinose$condition <- ifelse(arabinose$condition == "1", 1, 0) + 
  ifelse(arabinose$condition == "5", 2, 0) + 
  ifelse(arabinose$condition == "6", 3, 0) +
  ifelse(arabinose$condition == "7", 4, 0)


arabinose.plot <- ggplot(arabinose, aes(F2, t1, fill = F1))+
  geom_boxplot()+theme_classic()+
  theme(text = element_text(size = 15),
        legend.position = "top",
        legend.title = element_blank(),
        axis.title.x  = element_blank())+
  xlab("F2 environment")+
  ylab("Colony size (mm)")+
  labs(fill = "F1 environment")+
  scale_fill_manual(values = c("#4169E1", "#FC4E07"))

arabinose.m <- data.frame(gr_scaled = arabinose$gr_scaled, condition = factor(arabinose$condition))

arabinose_model1 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = arabinose.m, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

p <- precis(arabinose_model1, depth = 2)
plot(precis(arabinose_model1, depth = 2))
traceplot(arabinose_model1)
trankplot(arabinose_model1)
post.arabinose_m1 <- extract.samples(arabinose_model1)
#If current environment was sucrose
#What is the effect of F1 environment? in model 
diff1_4_arabinose_model1 <- post.arabinose_m1$a[,1] - post.arabinose_m1$a[,4]
mean(diff1_4_arabinose_model1 ) 
HPDI(diff1_4_arabinose_model1 , prob = 0.95) 
#If current environment was arabinose
#What is the effect of F1 environment?
diff2_3_arabinose_model1 <- post.arabinose_m1$a[,2] - post.arabinose_m1$a[,3]
mean(diff2_3_arabinose_model1)
HPDI(diff2_3_arabinose_model1, prob = 0.95) 

########################  LACTOSE #########
lactose.t <- c(1,8,9,10)
lactose <- carbons %>% filter(condition %in% lactose.t)
lactose$gr_scaled <- scale(lactose$t1)
lactose$condition <- ifelse(lactose$condition == "1", 1, 0) + 
  ifelse(lactose$condition == "8", 2, 0) + 
  ifelse(lactose$condition == "9", 3, 0) +
  ifelse(lactose$condition == "10", 4, 0)

lactose.plot <- ggplot(lactose, aes(F2, t1, fill = F1))+
  geom_boxplot()+
  theme_classic()+
  theme(text = element_text(size = 15),
        legend.position = "top",
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        legend.title = element_blank())+
  xlab("F2 environment")+
  ylab("Colony size (mm)")+
  labs(fill = "F1 environment")+
  scale_fill_manual(values = c("#008000", "#FC4E07"))

lactose.m <- data.frame(gr_scaled = lactose$gr_scaled, condition = factor(lactose$condition))

lactose_model1 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = lactose.m, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

p <- precis(lactose_model1, depth = 2)
plot(precis(lactose_model1, depth = 2))
traceplot(lactose_model1)
trankplot(lactose_model1)
post.lactose_m1 <- extract.samples(lactose_model1)
#If current environment was sucrose
#What is the effect of F1 environment? in model 
diff1_4_lactose_model1 <- post.lactose_m1$a[,1] - post.lactose_m1$a[,4]
mean(diff1_4_lactose_model1 ) 
HPDI(diff1_4_lactose_model1 , prob = 0.95) 
#If current environment was lactose
#What is the effect of F1 environment?
diff2_3_lactose_model1 <- post.lactose_m1$a[,2] - post.lactose_m1$a[,3]
mean(diff2_3_lactose_model1)
HPDI(diff2_3_lactose_model1, prob = 0.95) 

########################  MALTOSE #########
maltose.t <- c(1,11,12,13)
maltose <- carbons %>% filter(condition %in% maltose.t)
maltose$gr_scaled <- scale(maltose$t1)
maltose$condition <- ifelse(maltose$condition == "1", 1, 0) + 
  ifelse(maltose$condition == "11", 2, 0) + 
  ifelse(maltose$condition == "12", 3, 0) +
  ifelse(maltose$condition == "13", 4, 0)

maltose.plot <- ggplot(maltose, aes(F2, t1, fill = F1))+
  geom_boxplot()+theme_classic()+
  theme(text = element_text(size = 15),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        legend.title = element_blank(),
        legend.position = "top")+
  xlab("F2 environment")+
  ylab("Colony size (mm)")+
  labs(fill = "F1 environment")+
  scale_fill_manual(values = c("#20B2AA", "#FC4E07"))

maltose.m <- data.frame(gr_scaled = maltose$gr_scaled, condition = factor(maltose$condition))

maltose_model1 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = maltose.m, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

p <- precis(maltose_model1, depth = 2)
plot(precis(maltose_model1, depth = 2))
traceplot(maltose_model1)
trankplot(maltose_model1)
post.maltose_m1 <- extract.samples(maltose_model1)
#If current environment was sucrose
#What is the effect of F1 environment? in model 
diff1_4_maltose_model1 <- post.maltose_m1$a[,1] - post.maltose_m1$a[,4]
mean(diff1_4_maltose_model1 ) 
HPDI(diff1_4_maltose_model1 , prob = 0.95) 
#If current environment was maltose
#What is the effect of F1 environment?
diff2_3_maltose_model1 <- post.maltose_m1$a[,2] - post.maltose_m1$a[,3]
mean(diff2_3_maltose_model1)
HPDI(diff2_3_maltose_model1, prob = 0.95) 

########################  XYLOSE #########
xylose.t <- c(1,14,15,16)
xylose <- carbons %>% filter(condition %in% xylose.t)
xylose$gr_scaled <- scale(xylose$t1)
xylose$condition <- ifelse(xylose$condition == "1", 1, 0) + 
  ifelse(xylose$condition == "14", 2, 0) + 
  ifelse(xylose$condition == "15", 3, 0) +
  ifelse(xylose$condition == "16", 4, 0)


xylose.plot <- ggplot(xylose, aes(F2, t1, fill = F1))+
  geom_boxplot()+theme_classic()+
  theme(text = element_text(size = 15),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        legend.title = element_blank(),
        legend.position = "top")+
  xlab("F2 environment")+
  ylab("Colony size (mm)")+
  scale_fill_manual(values = c("#6495ED", "#FC4E07"))

xylose.m <- data.frame(gr_scaled = xylose$gr_scaled, condition = factor(xylose$condition))

xylose_model1 <- ulam(
  alist(
    gr_scaled ~ dnorm(mu, sigma),
    mu <- a[condition],
    a[condition] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = xylose.m, warmup = 1000, iter = 3000, chains = 4, log_lik = T)

p <- precis(xylose_model1, depth = 2)
plot(precis(xylose_model1, depth = 2))
traceplot(xylose_model1)
trankplot(xylose_model1)
post.xylose_m1 <- extract.samples(xylose_model1)
#If current environment was sucrose
#What is the effect of F1 environment? in model 
diff1_4_xylose_model1 <- post.xylose_m1$a[,1] - post.xylose_m1$a[,4]
mean(diff1_4_xylose_model1 ) 
HPDI(diff1_4_xylose_model1 , prob = 0.95) 
#If current environment was xylose
#What is the effect of F1 environment?
diff2_3_xylose_model1 <- post.xylose_m1$a[,2] - post.xylose_m1$a[,3]
mean(diff2_3_xylose_model1)
HPDI(diff2_3_xylose_model1, prob = 0.95) 

#Protein assay
#Samples from the experiment done with big slants on 02.08.22
p1 <- read.csv("protein_assay020822.csv")
head(p1)
class(p1$Sample)
#Subtract the blank in all the strds
#just working with measurement 1(abs1)
p1$abs_final <- p1$abs1-0.174 #0.174 is the mean of both blank(standard I) measurements
#extract standars for standard curve
p2 <- p1[1:18,]
#extract the unknowns
p3 <- p1[19:78,]
#graph the strds
ggplot(p2, aes(x = protein, y = abs_final))+
  geom_point() +
  geom_smooth()
#do the regression
line <- lm(p2$abs_final ~ p2$protein)
# extract intercept and slope from our line object
int <- line$coefficients[1]
slope <- line$coefficients[2]
#calculate the concentration (roundede with two decimals)
p3$protein <- round((p3$abs_final - int)/slope, 2)

protein.plot <- ggplot(p3, aes(x = Treatment, y = protein))+
  geom_boxplot(lwd=1, fatten=2)+
  theme_classic()+
  theme(text = element_text(size=30))+
  ylab("Protein concentration ug/ml")+
  xlab("F2 environment sucrose (%)")

###### Glycogen assay ###############
#the csv file has all the corrections for the blank and the glucose controls
ga1 <- read.csv("glycogen_assay120822.csv")
head(g1)
#extract standars for standard curve
ga2 <- ga1[1:12,]
#extract the unknowns
ga3 <- ga1[13:44,]
#graph the strds
ggplot(ga2, aes(x = amount, y = control_reads))+
  geom_point() +
  geom_smooth()
#do the regression
line <- lm(ga2$control_reads ~ ga2$amount)
# extract intercept and slope from our line object
int <- line$coefficients[1]
slope <- line$coefficients[2]
#calculate the concentration (roundede with two decimals)
ga3$amount <- round((ga3$control_reads - int)/slope, 2)

glycogen.plot <- ggplot(ga3, aes(x = treatment, y = amount))+
  geom_boxplot(lwd=1, fatten=2)+
  theme_classic()+
  theme(text = element_text(size=30))+
  ylab("Amount glycogen")+
  xlab("F2 environment sucrose (%)")

########################################################################################################################################
#Glucose samples
gu1 <- read.csv("glucose_assay141022.csv")
head(gu1)
#extract standars for standard curve
gu2 <- gu1[1:12,]
#extract the unknowns
gu3 <- gu1[13:44,]
#graph the strds
ggplot(gu2, aes(x = amount, y = blank.corrected))+
  geom_point() +
  geom_smooth()
#do the regression
line <- lm(gu2$blank.corrected ~ gu2$amount)
# extract intercept and slope from our line object
int <- line$coefficients[1]
slope <- line$coefficients[2]
#calculate the concentration (roundede with two decimals)
gu3$amount <- round((gu3$control.reads - int)/slope, 2)

glucose.plot <- ggplot(gu3, aes(x = treatment, y = amount))+
  geom_boxplot(lwd=1, fatten=2)+
  theme_classic()+
  theme(text = element_text(size=30))+
  ylab("Amount glucose nm/l")+
  xlab("F2 environment sucrose (%)")

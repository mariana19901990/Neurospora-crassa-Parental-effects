#Competition experiment analysis

#Load HRM and qPCR functions
#source("~/Documents/tutkijatohtori/scripts/HRMqPCR.R")
source("HRMqPCR.R")

#For data manipulation and plotting
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(colorspace)

#For Bayesian analysis
library(rethinking)
options(mc.cores = parallel::detectCores())
library(tidybayes)

later:::ensureInitialized() #Need to run this at home laptop

### * Load and process the raw data (needs to be run only once)

### ** Data for first set of experiments

### *** Load and process the qPCR plates

## Need to load each HRM-plate separately and estimate csr-proportions and then combine the data

#plate1.RFU <- loadHrmDataFile("~/Documents/tutkijatohtori/tgenerational/HRM/190812/190812_hrm.csv")
#Load RFU data, no need to modify the file manually
plate1.RFU <- loadHrmDataFile("190812_hrm.csv")
#plate1.samples <- read.csv("~/Documents/tutkijatohtori/tgenerational/HRM/190812/190812_samples.csv", header = T)
plate1.samples <- read.csv("190812_samples.csv", header = T)

#plate2.RFU <- loadHrmDataFile("~/Documents/tutkijatohtori/tgenerational/HRM/190812_2/190812_2_hrm.csv")
#plate2.samples <- read.csv("~/Documents/tutkijatohtori/tgenerational/HRM/190812_2/190812_2_samples.csv", header = T)
plate2.RFU <- loadHrmDataFile("190812_2_hrm.csv")
plate2.samples <- read.csv("190812_2_samples.csv", header = T)

#plate3.RFU <- loadHrmDataFile("~/Documents/tutkijatohtori/tgenerational/HRM/190812_3/190812_3_hrm.csv")
#plate3.samples <- read.csv("~/Documents/tutkijatohtori/tgenerational/HRM/190812_3/190812_3_samples.csv", header = T)
plate3.RFU <- loadHrmDataFile("190812_3_hrm.csv")
plate3.samples <- read.csv("190812_3_samples.csv", header = T)

##Define active melt region
plate1.RFU <- filter(plate1.RFU, Temperature >= 75 & Temperature <= 87)
plate2.RFU <- filter(plate2.RFU, Temperature >= 75 & Temperature <= 87)
plate3.RFU <- filter(plate3.RFU, Temperature >= 75 & Temperature <= 87)

#Normalization of RFU data
plate1.norm <- data.frame(Temperature = plate1.RFU[,1], apply(plate1.RFU[,-1], 2, RFUnorm))
plate2.norm <- data.frame(Temperature = plate2.RFU[,1], apply(plate2.RFU[,-1], 2, RFUnorm))
plate3.norm <- data.frame(Temperature = plate3.RFU[,1], apply(plate3.RFU[,-1], 2, RFUnorm))

#Calculating difference data, normalizing for B12 (positive control for csr-tag)
plate1.dif <- plate1.norm[,-1] - plate1.norm$B12
plate2.dif <- plate2.norm[,-1] - plate2.norm$B12
plate3.dif <- plate3.norm[,-1] - plate3.norm$B12
#Completing data frame
plate1.dif <- data.frame(Temperature = plate1.RFU[,1], plate1.dif)
plate2.dif <- data.frame(Temperature = plate2.RFU[,1], plate2.dif)
plate3.dif <- data.frame(Temperature = plate3.RFU[,1], plate3.dif)
#Getting the row indices that contain the maximal difference
comprow1 <- which(plate1.dif$B11 == max(plate1.dif$B11)) #81.4C
comprow2 <- which(plate2.dif$B11 == max(plate2.dif$B11))
comprow3 <- which(plate3.dif$B11 == max(plate3.dif$B11))

#Store maximal melting curve difference
plate1.samples$diff <- as.numeric(plate1.dif[comprow1,-1])
plate2.samples$diff <- as.numeric(plate2.dif[comprow2,-1])
plate3.samples$diff <- as.numeric(plate3.dif[comprow3,-1])

##Extracting standard curves only
Bstd1 <- filter(plate1.samples, sampletype == "std")
Bstd2 <- filter(plate2.samples, sampletype == "std")
Bstd3 <- filter(plate3.samples, sampletype == "std")

### *** Fitting the standard curves for the first set of experiments

##Using the model from "Marked Neurospora strains..." manuscript
std.model1 <- map2stan(
    alist(
        diff ~ dnorm(mu, sigma),
        mu <- a + B*proportion,
        a ~ dnorm(0,10),
        B ~ dnorm(0,10),
        sigma ~ dcauchy(0,2)
        ),
    data = Bstd1, chains = 2, cores = 2, warmup = 1000, iter = 3000)

##Use the same model for all of the plates
std.model2 <- map2stan(std.model1, data = Bstd2, chains = 2, cores = 2, warmup = 1000, iter = 3000)

std.model3 <- map2stan(std.model1, data = Bstd3, chains = 2, cores = 2, warmup = 1000, iter = 3000)

##Extract posterior samples
post1 <- extract.samples(std.model1)
post2 <- extract.samples(std.model2)
post3 <- extract.samples(std.model3)

#Plotting the standard curve
#pred.dif <- matrix(seq(from = -0.05, to = 0.4, by = 0.001), nrow = 1)
#convr <- apply(pred.dif, 2, convert2prop, post = post1)
#conv.med <- apply(convr, 2, median)
#conv.int <- apply(convr, 2, HPDI, prob = 0.95)

#ggplot() +
#    geom_ribbon(aes(ymin = conv.int[1,], ymax = conv.int[2,], x = as.vector(pred.dif)), fill = "grey70", alpha = 0.5) +
#    geom_line(aes(x = as.vector(pred.dif), y = conv.med)) +
#    geom_point(aes(x = Bstd1$diff, y = Bstd1$proportion)) +
#    xlab("Difference (norm. RFU)") +
#    ylab("csr-1* proportion") +
#    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.1))

##Proportions in unknown samples

plate1.samples <- filter(plate1.samples, sampletype == "sample")
plate2.samples <- filter(plate2.samples, sampletype == "sample")
plate3.samples <- filter(plate3.samples, sampletype == "sample")

##For one plate
foo <- list(0)
results <- rep(foo, nrow(plate1.samples))
for(i in 1:length(results)) { results[[i]] <- convert2prop(post1, plate1.samples[i,9]) }
plate1.samples$propest <- unlist(lapply(results, median))
prop.hpdi <- matrix(unlist(lapply(results, HPDI, prob = 0.95)), ncol = 2, byrow = T)
plate1.samples$prophigh <- prop.hpdi[,1]
plate1.samples$proplow <- prop.hpdi[,2]
plate1.samples$propsd <- unlist(lapply(results, sd))

##Second plate
foo <- list(0)
results <- rep(foo, nrow(plate2.samples))
for(i in 1:length(results)) { results[[i]] <- convert2prop(post2, plate2.samples[i,9]) }
plate2.samples$propest <- unlist(lapply(results, median))
prop.hpdi <- matrix(unlist(lapply(results, HPDI, prob = 0.95)), ncol = 2, byrow = T)
plate2.samples$prophigh <- prop.hpdi[,1]
plate2.samples$proplow <- prop.hpdi[,2]
plate2.samples$propsd <- unlist(lapply(results, sd))

##Third plate
foo <- list(0)
results <- rep(foo, nrow(plate3.samples))
for(i in 1:length(results)) { results[[i]] <- convert2prop(post3, plate3.samples[i,9]) }
plate3.samples$propest <- unlist(lapply(results, median))
prop.hpdi <- matrix(unlist(lapply(results, HPDI, prob = 0.95)), ncol = 2, byrow = T)
plate3.samples$prophigh <- prop.hpdi[,1]
plate3.samples$proplow <- prop.hpdi[,2]
plate3.samples$propsd <- unlist(lapply(results, sd))

###Making combined dataset
comp.samples <- rbind(plate1.samples, plate2.samples, plate3.samples)

##compet ID mat A (with csr*)  G1.5% (with csr*)
#1           1                  0                 
#2           -1                 0
#3           1                  1
#4           -1                 1
#5           1                  -1
#6           -1                 -1
#7           1                  0
#8           -1                 0

matvec <- c(1,-1,1,-1,1,-1,1,-1)
parvec <- c(0,0,1,1,-1,-1,0,0)
matA <- rep(0, dim(comp.samples)[1])
g1 <- rep(0, dim(comp.samples)[1])

for(i in 1:dim(comp.samples)[1]) {matA[i] <-  matvec[comp.samples$competid[i]] }
for(i in 1:dim(comp.samples)[1]) {g1[i] <- parvec[comp.samples$competid[i]] }

comp.samples$matA <- matA
comp.samples$g1 <- g1
comp.samples$envg2 <- ifelse(comp.samples$env == "1.5", 1, 0)

##Fixing the generation 0 samples, so that they have population ID's
##Also that proportions start at 0.5
temp <- colnames(comp.samples)

te <- filter(comp.samples, generation == 0)

test <- data.frame(rep(te$Sample, 10), rep(te$proportion,10), rep(te$sampletype,10), rep(te$generation,10), rep(te$competid,10), rep(te$population,10), rep(te$env,10), rep(te$Cq,10), rep(te$diff,10), rep(0.5,8*10), rep(0.5,8*10), rep(0.5,8*10), rep(0.01,8*10), rep(te$matA,10), rep(te$g1,10), rep(te$envg2,10))
colnames(test) <- temp

test <- arrange(test, competid)
test$population <- 1:80

#comp2 <- rbind(filter(comp.samples, generation == 1), test)

comp3 <- rbind(filter(comp.samples, generation != 0), test)

comp3$env[161:240] <- rep(c(rep(0.015,5), rep(1.5,5)),8)
comp3$env <- factor(comp3$env)
#recode(comp3$env, '0.015' = "F2 assay 0.015%", '1.5' = "F2 assay 1.5%", .default = levels(comp3$env)) 

comp3 <- mutate(comp3, env = recode(comp3$env, '0.015' = "F2 assay 0.015%", '1.5' = "F2 assay 1.5%", .default = levels(comp3$env)) )

save(comp3, file = "comexp1.RData")



### ** Data for the second set of experiments

### *** Load and process the qPCR plates for the second set

#plate1.RFU <- loadHrmDataFile("~/Documents/tutkijatohtori/tgenerational/HRM/190917_PS/190917_PS_hrm.csv") #Load RFU data, no need to modify the file manually
#plate1.samples <- read.csv("~/Documents/tutkijatohtori/tgenerational/HRM/190917_PS/190917_PS_samples.csv", header = T, sep = ";", dec = ",")
plate1.RFU <- loadHrmDataFile("190917_PS_hrm.csv") #Load RFU data, no need to modify the file manually
plate1.samples <- read.csv("190917_PS_samples.csv", header = T, sep = ";", dec = ",")


#Look at the raw data
ggplot(plate1.RFU, aes(x = Temperature, y = A1)) +
    geom_line()

##Define active melt region
plate1.RFU <- filter(plate1.RFU, Temperature >= 75 & Temperature <= 87)

#Normalization of RFU data
plate1.norm <- data.frame(Temperature = plate1.RFU[,1], apply(plate1.RFU[,-1], 2, RFUnorm))

#Look at the normalized data
ggplot(plate1.norm, aes(x = Temperature, y = A1)) +
    geom_line()

#Calculating difference data, normalizing for B12 (positive control for csr-tag)
plate1.dif <- plate1.norm[,-1] - plate1.norm$B12
#Completing data frame
plate1.dif <- data.frame(Temperature = plate1.RFU[,1], plate1.dif)
#Getting the row indices that contain the maximal difference
comprow1 <- which(plate1.dif$B11 == max(plate1.dif$B11)) #81.6

#Store maximal melting curve difference
plate1.samples$diff <- as.numeric(plate1.dif[comprow1,-1])

##Extracting standard curves only
Bstd1 <- filter(plate1.samples, Sampletype == "std")

### *** Fitting the standard curves for the second set

##Using the model from "Marked Neurospora strains..." manuscript
std.model1 <- map2stan(
    alist(
        diff ~ dnorm(mu, sigma),
        mu <- a + B*Proportion,
        a ~ dnorm(0,10),
        B ~ dnorm(0,10),
        sigma ~ dcauchy(0,2)
        ),
    data = Bstd1, chains = 2, cores = 2, warmup = 1000, iter = 3000)


##Extract posterior samples
post1 <- extract.samples(std.model1)

#Plotting the standard curve
pred.dif <- matrix(seq(from = -0.05, to = 0.4, by = 0.001), nrow = 1)
convr <- apply(pred.dif, 2, convert2prop, post = post1)
conv.med <- apply(convr, 2, median)
conv.int <- apply(convr, 2, HPDI, prob = 0.95)

ggplot() +
    geom_ribbon(aes(ymin = conv.int[1,], ymax = conv.int[2,], x = as.vector(pred.dif)), fill = "grey70", alpha = 0.5) +
    geom_line(aes(x = as.vector(pred.dif), y = conv.med)) +
    geom_point(aes(x = Bstd1$diff, y = Bstd1$Proportion)) +
    xlab("Difference (norm. RFU)") +
    ylab("csr-1* proportion") +
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.1))

##Proportions in unknown samples

plate1.samples <- filter(plate1.samples, Sampletype == "sample")

##For one plate
foo <- list(0)
results <- rep(foo, nrow(plate1.samples))
for(i in 1:length(results)) { results[[i]] <- convert2prop(post1, plate1.samples[i,9]) }
plate1.samples$propest <- unlist(lapply(results, median))
prop.hpdi <- matrix(unlist(lapply(results, HPDI, prob = 0.95)), ncol = 2, byrow = T)
plate1.samples$prophigh <- prop.hpdi[,2]
plate1.samples$proplow <- prop.hpdi[,1]
plate1.samples$propsd <- unlist(lapply(results, sd))

##Fixing competition ID's
##compet ID mat A (with csr*)  G1.5% (with csr*)
#1           1                  0                 
#2           -1                 0
#3           1                  1
#4           -1                 1
#5           1                  -1
#6           -1                 -1
#7           1                  0
#8           -1                 0

#Everytime csr* is with mat A, matvec gets 1, when csr* is with mat a matvec gets -1
matvec <- c( rep(1,5), rep(1,5), rep(-1,5), rep(-1,5), rep(1,5), rep(-1,5), rep(1,5), rep(1,5), rep(-1,5), rep(-1,5), rep(1,5), rep(-1,5)) #

#Everytime csr* is with 1.5% parental environment, parvec gets 1; when csr* is with 0.015% parental environment parvec gets -1; and when both parental environments are the same parvec gets 0
parvec <- c( rep(0,5), rep(1,5), rep(1,5), rep(-1,5), rep(-1,5), rep(0,5), rep(0,5), rep(1,5), rep(1,5), rep(-1,5), rep(-1,5), rep(0,5))

#Environment where competitions happened 1 = 1.5%, 0 = 0.015%
envvec <- c( rep(1,5*6), rep(0,5*6))

#Storing
plate1.samples$matA <- matvec
plate1.samples$parental <- parvec
plate1.samples$environment <- envvec

##Fixing the generation 0 samples, so that they have population ID's
##Also that proportions start at 0.5

#Starting proportions at generation 0
temp <- plate1.samples
temp$Generation <- 0
temp$Cq <- NA
temp$diff <- NA
temp$propest <- 0.5
temp$prophigh <- 0.5
temp$proplow <- 0.5
temp$propsd <- 0.001

plate1.samples <- rbind(plate1.samples, temp)

#Fixing zero sd's
plate1.samples[plate1.samples$propsd == 0,13] <- 0.0001

##Okay dataset is ready!

save(plate1.samples, file = "comexp2.RData")

### * Analysis of the combined data (Start from here)

### ** Load the data

load("comexp1.RData")
load("comexp2.RData")

comp3$envg2 <- ifelse(comp3$env == "F2 assay 1.5%", 1, 0) #Fix environmental assay
plate1.samples$Population <- plate1.samples$Population + 80 #Fix population numbers for second exp.
colnames(plate1.samples)[15] <- "g1"
colnames(plate1.samples)[16] <- "envg2"
colnames(plate1.samples)[4] <- "generation"
colnames(plate1.samples)[5] <- "competid"
colnames(plate1.samples)[6] <- "population"

temp1 <- comp3[,c(4,5,6,8,9,10,11,12,13,14,15,16)]
temp2 <- plate1.samples[,c(4,5,6,8,9,10,11,12,13,14,15,16)]
allcomp <- rbind(temp1, temp2)

##Making a plot of frequency trajectories

##Making the factors
strain1 <- ifelse(allcomp$matA == 1, "mat A csr-1*", "mat a csr-1*")
strain2 <- ifelse(allcomp$matA == 1, "mat a", "mat A")
s1G1.temp <- ifelse(allcomp$competid[1:240] <= 4, "G1 1.5%", "G1 0.015%")
s2G1.temp <- rep("", 240)
for(i in 1:240) {
    if(allcomp$competid[i] == 1) {s2G1.temp[i] <- "G1 1.5%"}
    if(allcomp$competid[i] == 2) {s2G1.temp[i] <- "G1 1.5%"}
    if(allcomp$competid[i] == 3) {s2G1.temp[i] <- "G1 0.015%"}
    if(allcomp$competid[i] == 4) {s2G1.temp[i] <- "G1 0.015%"}
    if(allcomp$competid[i] == 5) {s2G1.temp[i] <- "G1 1.5%"}
    if(allcomp$competid[i] == 6) {s2G1.temp[i] <- "G1 1.5%"}
    if(allcomp$competid[i] == 7) {s2G1.temp[i] <- "G1 0.015%"}
    if(allcomp$competid[i] == 8) {s2G1.temp[i] <- "G1 0.015%"}
}

s1G1.temp2 <- ifelse(allcomp$competid[241:360] == 7 | allcomp$competid[241:360] == 12, "G1 0.015%", "G1 1.5%")
s2G1.temp2 <- ifelse(allcomp$competid[241:360] == 1 | allcomp$competid[241:360] == 6, "G1 1.5%", "G1 0.015%")

s1G1 <- c(s1G1.temp, s1G1.temp2)
s2G1 <- c(s2G1.temp, s2G1.temp2)

allcomp$strain1 <- factor(strain1)
allcomp$strain2 <- factor(strain2)
allcomp$s1G1 <- factor(s1G1)
allcomp$s2G1 <- factor(s2G1)

allcomp.env1 <- filter(allcomp, envg2 == 1)
allcomp.env2 <- filter(allcomp, envg2 == 0)


freqsenv1 <- ggplot(allcomp.env1, aes(x = generation, y = propest, ymin = proplow, ymax = prophigh, group = population)) +
    geom_pointrange() +
    geom_line() +
    facet_grid(strain1 + s1G1 ~ strain2 + s2G1) +
    scale_x_continuous(breaks = c(0,1,2), expand = c(0.1, 0.1)) +
    xlab("Transfer") +
    ylab("csr-1* proportion")


freqsenv2 <- ggplot(allcomp.env2, aes(x = generation, y = propest, ymin = proplow, ymax = prophigh, group = population)) +
    geom_pointrange() +
    geom_line() +
    facet_grid(strain1 + s1G1 ~ strain2 + s2G1) +
    scale_x_continuous(breaks = c(0,1,2), expand = c(0.1, 0.1)) +
    xlab("Transfer") +
    ylab("csr-1* proportion")

compfreq <- plot_grid(freqsenv1, freqsenv2, ncol = 2, labels = c("A", "B"))
#save_plot("./suc_ms/fig/competitionfreq.pdf", compfreq, base_height = 3.71*1.1618*2, base_width = 3.71*4)

### ** Analysis

wmodel <- map2stan(
    alist(
        prop ~ dnorm(mu, sigma),
        logit(mu) <- a_comp + (B_csr + B_pop[population] + B_matA*matA + B_g1*g1)*generation,
        propest ~ dnorm(prop, propsd),
        a_comp ~ dnorm(0,0.065),
        B_pop[population] ~ dnorm(0,sigmap),
        B_csr ~ dnorm(0,1),
        B_matA ~ dnorm(0,1),
        B_g1 ~ dnorm(0,1),
        sigma ~ dcauchy(0,2),
        sigmap ~ dcauchy(0,2)
        ),
    data = allcomp, start = list(prop = allcomp$propest), WAIC = F, warmup = 1000, iter = 5000, chains = 2, cores = 2, control = list(adapt_delta = 0.9))

save(wmodel, file = "compmodel.RData")

wpost <- extract.samples(wmodel)

wdata <- data.frame(exp(wpost$B_csr), exp(wpost$B_matA), exp(wpost$B_g1))
colnames(wdata) <- c("Bcsr", "BmatA", "Bg1")
wdatalong <- gather(wdata, key = "effect", value = "s")
wdatalong$effect <- factor(wdatalong$effect)

##Making a table of fitness effects
Wtable <- rbind(
quantile(wdata$Bcsr, probs = c(0.025, 0.5, 0.975)),
quantile(wdata$BmatA, probs = c(0.025, 0.5, 0.975)),
quantile(wdata$Bg1, probs = c(0.025, 0.5, 0.975)))
#quantile(wdata$Benv, probs = c(0.025, 0.5, 0.975)) )
rownames(Wtable) <- c("csr-1*", "mat A", "G1 1.5%")


#Plotting the 
pdf("./suc_ms/fig/fitness.pdf")
ggplot(wdatalong, aes(x = s, y = effect)) +
    geom_halfeyeh(fill = "#66CCFF", point_interval = median_hdi) +
    geom_vline(xintercept = 1, lty = "dashed") +
    scale_y_discrete(labels = c("csr-1*", "G1 1.5%", "mat A")) +
    scale_x_continuous(breaks = seq(0,6, 0.5)) +
    ylab(NULL) +
    xlab(expression(W[ij]))
dev.off()

##Testing the model model with both environments separately
comp1.5 <- filter(allcomp, envg2 == 1)
comp0.015 <- filter(allcomp, envg2 == 0)

w1.5 <- map2stan(
    alist(
        prop ~ dnorm(mu, sigma),
        logit(mu) <- a_comp + (B_csr + B_pop[population] + B_matA*matA + B_g1*g1)*generation,
        propest ~ dnorm(prop, propsd),
        a_comp ~ dnorm(0,0.065),
        B_pop[population] ~ dnorm(0,sigmap),
        B_csr ~ dnorm(0,1),
        B_matA ~ dnorm(0,1),
        B_g1 ~ dnorm(0,1),
        sigma ~ dcauchy(0,2),
        sigmap ~ dcauchy(0,2)
        ),
    data = comp1.5, start = list(prop = comp1.5$propest), WAIC = F, warmup = 1000, iter = 5000, chains = 2, cores = 2, control = list(adapt_delta = 0.9))

w0.015 <- map2stan(
    alist(
        prop ~ dnorm(mu, sigma),
        logit(mu) <- a_comp + (B_csr + B_pop[population] + B_matA*matA + B_g1*g1)*generation,
        propest ~ dnorm(prop, propsd),
        a_comp ~ dnorm(0,0.065),
        B_pop[population] ~ dnorm(0,sigmap),
        B_csr ~ dnorm(0,1),
        B_matA ~ dnorm(0,1),
        B_g1 ~ dnorm(0,1),
        sigma ~ dcauchy(0,2),
        sigmap ~ dcauchy(0,2)
        ),
    data = comp0.015, start = list(prop = comp0.015$propest), WAIC = F, warmup = 1000, iter = 5000, chains = 2, cores = 2, control = list(adapt_delta = 0.9))

##Extract posteriors from the w1.5 and w0.015 models and compare results to ones using all of the data
wpost1.5 <- extract.samples(w1.5)
wdata1.5 <- data.frame(exp(wpost1.5$B_csr), exp(wpost1.5$B_matA), exp(wpost1.5$B_g1))
colnames(wdata1.5) <- c("Bcsr", "BmatA", "Bg1")
Wtable1.5 <- t(apply(wdata1.5, 2, quantile, probs = c(0.025, 0.5, 0.975)))
rownames(Wtable1.5) <- c("1.5% csr-1*", "1.5% mat A", "1.5% G1 1.5%")

wpost0.015 <- extract.samples(w0.015)
wdata0.015 <- data.frame(exp(wpost0.015$B_csr), exp(wpost0.015$B_matA), exp(wpost0.015$B_g1))
colnames(wdata0.015) <- c("Bcsr", "BmatA", "Bg1")
Wtable0.015 <- t(apply(wdata0.015, 2, quantile, probs = c(0.025, 0.5, 0.975)))
rownames(Wtable0.015) <- c("0.015% csr-1*", "0.015% mat A", "0.015% G1 1.5%")

Wtable <- rbind(Wtable, Wtable1.5, Wtable0.015) #Making the final Wtable

save(Wtable, file = "Wtable.RData")

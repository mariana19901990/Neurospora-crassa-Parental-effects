setwd("/home/mariana/Documents/Mariana/Neurospora")
#Load libraries
library(tidyr)
library(dplyr)
library(forcats)
library(ggplot2)
library(cowplot)
library(brms) #Interfacing Stan for Bayesian models
library(tidybayes)
library(coda)
library(MASS)
library(colorspace)
theme_set(theme_cowplot())

### ** Load and process the data

aineisto <- read.csv("germination.csv", header = T)

## We can combine the three photpgraphs taken from the same well, to get a better estimate of the proportion

aineisto <- group_by(aineisto, Plate, Well, Time, tube, condition, parental_env)
aineisto <- summarise(aineisto, totalspores = sum(total), germ = sum(germinated), dorm = sum(total) - sum(germinated))
aineisto$tube <- factor(aineisto$tube) #Change to factor
aineisto$condition <- factor(aineisto$condition) #Change to factor
aineisto <- filter(aineisto, Time < 5) #Filter observations for time 5, these were unclear to score

### ** Modeling germination rate

## We are modeling germination rate as a function of time.

#logit(p_i) = log ( p_i / ( 1 - p_i) )
#inv_logit: p_i = exp(a+Bx_i) / (1 + exp(a+Bx_i))

model <- brm(data = aineisto, family = binomial,
             germ | trials(totalspores) ~ 0 + condition + Time:condition + (1 + Time | tube),
             iter = 2000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model) #Get posterior distributions
#What is the effect of F1 environment? if the current environment is 1.5%
diff1_4 <- post$b_condition1 - post$b_condition4
mean(diff1_4) #-0.177
HPDI(diff1_4, prob = 0.95) #-0.9882041  0.6788752
#What is the effect of F1 environment? if the current environment is 0.015%
diff2_3 <- post$b_condition2 - post$b_condition3
mean(diff2_3) #-1.922676
HPDI(diff2_3, prob = 0.95) #-2.743600 -1.178256 
#whats is the difference between matching environment
diff1_3 <- post$b_condition1 - post$b_condition3
mean(diff1_3) #-1.91934
HPDI(diff1_3, prob = 0.95) #-2.696808 -1.088676
#Predictions for plotting
#With the conditions included
#Note that totalspores has to be 1 in the prediction data frame, so to get probability of germination
pred.data <- data.frame(Time = rep(seq(0, 4, length.out = 100), 4), totalspores = rep(1, 100*4), condition = factor(sort(rep(1:4, 100))))
model.pred <- data.frame(fitted(model, re_formula = NA, newdata = pred.data))
model.pred <- cbind(pred.data, model.pred)
colnames(model.pred)[6:7] <- c("lower", "upper")

#Binomial graph
model.plot <- ggplot(model.pred, aes(x = Time, y = Estimate)) +
    geom_ribbon(data = model.pred, aes(x = Time, y = Estimate, ymin = lower, ymax = upper, fill = condition), alpha = 0.2) +
    geom_line(data = model.pred, aes(x = Time, y = Estimate, color = condition)) +    
    geom_point(data = aineisto, aes(x = Time, y = germ/totalspores, color = condition)) +
    theme_classic()+
    theme(legend.position = c(0.23, 0.88),
        text = element_text(size=20),
        legend.text = element_text(size=15))+
    xlab("Time (h)") +
    ylab("Probability of germination")+
  scale_fill_discrete(labels=c('F1 1.5%, F2 1.5%', 'F1 1.5%, F2 0.015%',
                               'F1 0.015%, F2 0.015%', 'F1 0.015%, F2 1.5%'),
                                name= element_blank())+
  scale_color_discrete(labels=c('F1 1.5%, F2 1.5%', 'F1 1.5%, F2 0.015%',
                                'F1 0.015%, F2 0.015%', 'F1 0.015%, F2 1.5%'),
                                name= element_blank())
  

### Getting estimates of threshold values, for p = 0.5
#Note that logit(p_i) = log ( p_i / ( 1 - p_i) )
# logit(0.5) = log( 0.5 / (0.5 - 1 ) ) = 0
#What we want is to solve for logit(y) = a + B*x, when y = 0.5
#Which is x = -a / B

#Do the calculations for the different conditions
#Condition 1
threshold.cond1 <- -post[,1] / post[,5]
cond.1 <- c(median(threshold.cond1), HPDinterval(as.mcmc(threshold.cond1)))

#Condition 2
threshold.cond2 <- -post[,2] / post[,6]
cond.2 <- c(median(threshold.cond2), HPDinterval(as.mcmc(threshold.cond2)))

#Condition 3
threshold.cond3 <- -post[,3] / post[,7]
cond.3 <- c(median(threshold.cond3), HPDinterval(as.mcmc(threshold.cond3)))

#Condition 4
threshold.cond4 <- -post[,4] / post[,8]
cond.4 <- c(median(threshold.cond4), HPDinterval(as.mcmc(threshold.cond4)))

#new data base with the p = 0.5 treshold values
tresh.data <- data.frame(rbind(cond.1, cond.2, cond.3, cond.4))
colnames(tresh.data) <- (c("p0.5", "lower", "upper"))
tresh.data$condition <- c(1:4)
tresh.data$F1_environment <- c("1.5%", "1.5%", "0.015%", "0.015%")
tresh.data$offspring_env <- c("1.5%", "0.015%", "0.015%", "1.5%")


germ.plot <- ggplot(tresh.data, aes(y=p0.5, x=offspring_env, color=F1_environment))+
  geom_point(size = 3, position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2,
                position=position_dodge(width=0.5))+
  theme_classic()+
  theme(legend.position = c(.85, .11),
        #axis.title=element_text(size=13),
        text = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))+
  ylab("Time to reach 50% germination (h)")+
  xlab("F2 sucrose environment")+
  scale_color_manual(values = c("#56B4E9", "#FC4E07"), 
                     name="F1 environment")+
  geom_signif(y_position =(4.0), xmin= c(.6), xmax = c(1.2),
              annotation = c("[-2.743, -1.178]"), color= "black")+
  geom_signif(y_position =(3.65), xmin= c(1.7), xmax = c(2.3),
              annotation = c("[-0.988, 0.678]"), color= "black")+
  geom_signif(y_position =(4.15), xmin= c(.6), xmax = c(2.3),
              annotation = c("[-2.696, -1.088]"), color= "black")
 
germination.plot <- ggdraw()+
  draw_plot(model.plot, x = 0, y = 0, width = 0.6, height = 0.95) +
  draw_plot(germ.plot, x = 0.62, y = 0, width = .37, height = .95) +
  draw_plot_label(label = c("A", "B"), size = 17,
                  x = c(.03, .65), y = c(.99, .99))

svg("germination.plot.svg",
    width=12, 
    height=7)
germination.plot
dev.off() 

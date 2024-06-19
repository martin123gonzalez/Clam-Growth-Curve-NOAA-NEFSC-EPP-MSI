# File Name: Sample_size_at_age_code.R
# Authors: Molly Roberts, Martin Gonzalez
# Purpose: To graph von Bertalanffy growhth function, Length Limit, OGP, and K 
         # using the whole chondrophore length (NOT incremental growth)





#Note to Molly & Matt: All code from Lines 11-292 was code used in my NOAA presentation in August.




# setwd("~/Documents/GitHub/EAD-ASEB-Ssolidissima-OA/data")


# ==== library ====
library(FSA)
library(FSAdata)
library(nlstools)
library(minpack.lm)
library(here)
library(TMB)
library(dplyr)
library(MASS)
library(ggplot2)
library(tidyverse)
    # library(car)
# install.packages("MASS")
# install.packages("TMB")
# install.packages("tidyverse")

# ==== data input ====
clam_data <- read.csv(here::here("Input_Data", "Clam_Metrics_Data_MAG.csv"))

head(clam_data)
table(clam_data$site)

clam_data_clean <- clam_data %>% 
  dplyr::filter(ages != "NA" & length != "NA") %>% 
  dplyr::filter(site != "?" & site != "Nobscusset") %>% 
  dplyr::select(site, length, ages, Param, Value) %>% 
  mutate(length = length * 10) %>% 
  mutate(Value = ifelse(Value == 15.30, Value*10, Value*1))
table(clam_data_clean$site)
clam_data_clean$site <-as.factor(clam_data_clean$site)
str(clam_data_clean)

  # abbreviation notes
    # sv = starting value 
    # vb = von Bertalanffy growhth function
    # letter = the variable that is not fixed
    # Linf = Limit for Length
    # K = Brody growth constant (speed of growth per year)
    # OGP = Overall Growth Performance
    # CH = Chondrophore
    # SL = Shell Length


#---- Real Code (only real data, uses clam_data_clean) ----
(svCom <- vbStarts(length~ages, data = clam_data_clean))           # base model
(svGen <- lapply(svCom, rep, 5))
vbGen <- length~Linf[site]*(1-exp(-K[site]*(ages-t0)))
fitGen <- nls(vbGen, data = clam_data_clean, start = svGen)
hist(residuals(fitGen), main = "")


##---- Formula and Model set up ----
vb1LK <- length~Linf[site]*(1-exp(-K[site]*(ages-t0))) # make sense
sv1LK <- mapply(rep,svCom,c(5,5,1))
vbCom <- length~Linf*(1-exp(-K*(ages-t0))) # make sense
vbComT0 <- length~Linf*(1-exp(-K*(ages))) # make sense
vb1LKT0 <- length~Linf[site]*(1-exp(-K[site]*(ages))) # make sense
vb1KT <- length~Linf*(1-exp(-K[site]*(ages-t0[site])))
sv1KT <- mapply(rep,svCom,c(1,5,5))
vb1LT <- length~Linf[site]*(1-exp(-K*(ages-t0[site])))
sv1LT <- mapply(rep,svCom,c(5,1,5))
vb2L <- length~Linf[site]*(1-exp(-K*(ages-t0)))
sv2L <- mapply(rep,svCom,c(5,1,1))
vb2T <- length~Linf*(1-exp(-K*(ages-t0[site])))
sv2T <- mapply(rep,svCom,c(1,1,5))
vb2K <- length~Linf*(1-exp(-K[site]*(ages-t0)))
sv2K <- mapply(rep,svCom,c(1,5,1))


fit1KT <- nls(vb1KT,data=clam_data_clean,start=sv1KT)
fit1LT <- nls(vb1LT,data=clam_data_clean,start=sv1LT)
fit1LK <- nls(vb1LK, data = clam_data_clean, start = sv1LK) # make sense
fit2T <- nls(vb2T,data=clam_data_clean,start=sv2T)
fit2K <- nls(vb2K,data=clam_data_clean,start=sv2K)
fit2L <- nls(vb2L, data = clam_data_clean, start = sv2L)
fitCom <- nls(vbCom, data = clam_data_clean, start = svCom) # make sense
fitComT0 <- nls(vbComT0, data = clam_data_clean,start=svCom[1:2]) # make sense
fit1LKT0 <- nls(vb1LKT0, data = clam_data_clean,start=sv1LK[1:2]) # make sense


## ---- ANOVAs and AIC ----

anova(fit1LK,fitCom) # Since Linf and K are correlated it doesn't make sense for 
# there to only be a difference in one of the parameters between the sites
AIC(fit1LK,fitCom, fitComT0, fit1LKT0)

AIC(fitGen,fit1KT,fit1LT,fit1LK,fit2T,fit2K,fit2L,fitCom,fitComT0,fit1LKT0) # Still within 2 AIC
    # Winner with context = 1LK or 1LKT0
vbTypical <- vbFuns("typical")

overview(fit1LK)

# Fit 2 L's only
# overview(fit2L)
# coef(fit2L)[-2]
# plot(Lengths~jitter(Ages,0.3),data=test_data,subset=Site=="Barn",pch=19,xlab="Age (yrs)",
#      ylab="Total Length (mm)", ylim = c(0,200), xlim = c(0,20))
# points(Lengths~jitter(Ages,0.3),data=test_data,subset=Site=="Eel",pch=19,col="gray")
# curve(vbTypical(x,Linf=coef(fit2L)[-2]),from=2,to=14,lwd=2,add=TRUE)
# curve(vbTypical(x,Linf=coef(fit2L)[-1]),from=2,to=8,col="gray",lwd=2,add=TRUE)
# legend("topleft",legend=c("Barnstable","Eel"),col=c("black","gray"),lwd=2,lty=1,cex=0.75)

## ---- sVBGF Graph ----

      # quality check with ggplot for my sanity about colors
# ggplot(data = clam_data_clean, aes(x = ages, y = length)) + 
#   geom_jitter(aes(color = site)) + 
#   geom_smooth(aes(color = site), se = F, method = "lm", formula = y ~ poly(x, 2))

# plot(curve(vbTypical(x,Linf=coef_avg),from = 0, to = 20, col = "grey", lwd = 2, add = TRUE), 
     # data = clam_data_clean, xlab = "x", ylab = "y", type = "1")


# png("Output/surfclam_sVBGF.png", width = 600, height = 600)  #un-comment this line to save image

levels(clam_data_clean$site)

plot(length~jitter(ages,0.8),data=clam_data_clean,subset=site=="Barnstable",pch=19,xlab="Age (yrs)",
     ylab="Total Length (mm)", ylim = c(0,200), xlim = c(0,15))
points(length~jitter(ages,0.8),data = clam_data_clean,subset = site == "Chatham",pch=19,col="red")
points(length~jitter(ages,0.8),data = clam_data_clean,subset = site == "Dennis",pch=19,col="orange")
points(length~jitter(ages,0.8),data = clam_data_clean,subset = site == "Eel Pond",pch=19,col="blue")
points(length~jitter(ages,0.8),data = clam_data_clean,subset = site == "Provincetown",pch=19,col="green")
# points(length~jitter(ages,0.8), data = clam_data_clean, pch = 19, col = "grey")
vbTypical <- vbFuns("typical")

# fot fit1KT
# overview(fit1KT)
# coef(fit1KT)
# coef1 <- c(coef(fit1KT)[1], coef(fit1KT)[2], coef(fit1KT)[7])
# coef2 <- c(coef(fit1KT)[1], coef(fit1KT)[3], coef(fit1KT)[8])
# coef3 <- c(coef(fit1KT)[1], coef(fit1KT)[4], coef(fit1KT)[9])
# coef4 <- c(coef(fit1KT)[1], coef(fit1KT)[5], coef(fit1KT)[10])
# coef5 <- c(coef(fit1KT)[1], coef(fit1KT)[6], coef(fit1KT)[11])

# for fit1LK
overview(fit1LK)
coef(fit1LK)
coef1 <- c(coef(fit1LK)[1], coef(fit1LK)[6], coef(fit1LK)[11]) #Barnstable
coef2 <- c(coef(fit1LK)[2], coef(fit1LK)[7], coef(fit1LK)[11]) #Chatham
coef3 <- c(coef(fit1LK)[3], coef(fit1LK)[8], coef(fit1LK)[11]) #Dennis
coef4 <- c(coef(fit1LK)[4], coef(fit1LK)[9], coef(fit1LK)[11]) #Eel Pond
coef5 <- c(coef(fit1LK)[5], coef(fit1LK)[10], coef(fit1LK)[11]) #Provincetown
coef_avg <- c(mean(coef(fit1LK)[1:5]), mean(coef(fit1LK)[6:10]), coef(fit1LK)[11])

curve(vbTypical(x,Linf=coef1),from=0,to=20,lwd=2,add=TRUE)
curve(vbTypical(x,Linf=coef2),from=0,to=20,col="red",lwd=2,add=TRUE)
curve(vbTypical(x,Linf=coef3),from=0,to=20,col="orange",lwd=2,add=TRUE)
curve(vbTypical(x,Linf=coef4),from=0,to=20,col="blue",lwd=2,add=TRUE)
curve(vbTypical(x,Linf=coef5),from=0,to=20,col="green",lwd=2,add=TRUE)
curve(vbTypical(x,Linf=coef_avg),from = 0, to = 20, col = "grey", lwd = 2, add = TRUE)

legend("bottomright",legend=c("Average", "N. Cape: Barnstable", "N. Cape: Dennis", "N. Cape: Provincetown", "S. Cape: Chatham", "S. Cape: Eel Pond"),
       col=c("grey", "black","orange", "green", "red", "blue"),lwd=2,lty=1,cex=.75, title = "Sites in Cape Cod")

title("Surfclam Growth Rate over Time")
# dev.off()                     #un-comment this line to save image


# Reference for reasoning to include T0 estimates (but probably a composite from all data is best 
# rather than site specific if only 20 data points per site and no juveniles): https://besjournals.onlinelibrary.wiley.com/doi/pdf/10.1111/2041-210x.12020

# OGP reference: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0199212

## ---- Linf, K, OGP Calculation ----
### ---- for 1LK ----
overview(fit1LK) #Note that t0 is not significantly different from 0 - the confidence intervals overlap with 0
summary(fit1LK)
coef(fit1LK)

# Linf calc
 Linf1 <- coef(fit1LK)[1]
 Linf2 <- coef(fit1LK)[2]
 Linf3 <- coef(fit1LK)[3]
 Linf4 <- coef(fit1LK)[4]
 Linf5 <- coef(fit1LK)[5]
 
 min.Linf1 <- confint.default(fit1LK)[1,1]
 max.Linf1 <- confint.default(fit1LK)[1,2]
 min.Linf2 <- confint.default(fit1LK)[2,1]
 max.Linf2 <- confint.default(fit1LK)[2,2]
 min.Linf3 <- confint.default(fit1LK)[3,1]
 max.Linf3 <- confint.default(fit1LK)[3,2]
 min.Linf4 <- confint.default(fit1LK)[4,1]
 max.Linf4 <- confint.default(fit1LK)[4,2]
 min.Linf5 <- confint.default(fit1LK)[5,1]
 max.Linf5 <- confint.default(fit1LK)[5,2]
 
# K calc
 K1 <- coef(fit1LK)[6]
 K2 <- coef(fit1LK)[7]
 K3 <- coef(fit1LK)[8]
 K4 <- coef(fit1LK)[9]
 K5 <- coef(fit1LK)[10]

 min.K1 <- confint.default(fit1LK)[6,1]
 max.K1 <- confint.default(fit1LK)[6,2]
 min.K2 <- confint.default(fit1LK)[7,1]
 max.K2 <- confint.default(fit1LK)[7,2]
 min.K3 <- confint.default(fit1LK)[8,1]
 max.K3 <- confint.default(fit1LK)[8,2]
 min.K4 <- confint.default(fit1LK)[9,1]
 max.K4 <- confint.default(fit1LK)[9,2]
 min.K5 <- confint.default(fit1LK)[10,1]
 max.K5 <- confint.default(fit1LK)[10,2]
 
# OGP calc
 (OGP1 <- as.numeric(log10(K1)+2*log10(0.1*Linf1)))
 (OGP2 <- as.numeric(log10(K2)+2*log10(0.1*Linf2)))
 (OGP3 <- as.numeric(log10(K3)+2*log10(0.1*Linf3)))
 (OGP4 <- as.numeric(log10(K4)+2*log10(0.1*Linf4)))
 (OGP5 <- as.numeric(log10(K5)+2*log10(0.1*Linf5)))

 (OGP1.min <- as.numeric(log10(min.K1)+2*log10(0.1*min.Linf1)))
 (OGP1.max <- as.numeric(log10(max.K1)+2*log10(0.1*max.Linf1)))
 (OGP2.min <- as.numeric(log10(min.K2)+2*log10(0.1*min.Linf2)))
 (OGP2.max <- as.numeric(log10(max.K2)+2*log10(0.1*max.Linf2)))
 (OGP3.min <- as.numeric(log10(min.K3)+2*log10(0.1*min.Linf3)))
 (OGP3.max <- as.numeric(log10(max.K3)+2*log10(0.1*max.Linf3)))
 (OGP4.min <- as.numeric(log10(min.K4)+2*log10(0.1*min.Linf4)))
 (OGP4.max <- as.numeric(log10(max.K4)+2*log10(0.1*max.Linf4)))
 (OGP5.min <- as.numeric(log10(min.K5)+2*log10(0.1*min.Linf5)))
 (OGP5.max <- as.numeric(log10(max.K5)+2*log10(0.1*max.Linf5)))

#### ---- Graphing and matrix/data frame making ----
 OGP_L_K_fit1LK_matrix <- matrix(c("Barnstable","Chatham", "Dennis", "Eel Pond", "Provincetown",
                                   Linf1, Linf2, Linf3, Linf4, Linf5, 
                                   K1, K2, K3, K4, K5,
                                   OGP1, OGP2, OGP3, OGP4, OGP5,
                                   min.Linf1, min.Linf2, min.Linf3, min.Linf4, min.Linf5,
                                   max.Linf1, max.Linf2, max.Linf3, max.Linf4, max.Linf5,
                                   min.K1, min.K2, min.K3, min.K4, min.K5,
                                   max.K1, max.K2, max.K3, max.K4, max.K5,
                                   OGP1.min, OGP2.min, OGP3.min, OGP4.min, OGP5.min,
                                   OGP1.max, OGP2.max, OGP3.max, OGP4.max, OGP5.max),
                                 nrow = 5, ncol = 10, byrow = F,
                                 dimnames = list(c(1, 2, 3, 4, 5), c("site", "Linf", "K", "OGP", "min.Linf",
                                                   "max.Linf", "min.K", "max.K", "OGP.min","OGP.max")))
 OGP_L_K_fit1LK <- data.frame(OGP_L_K_fit1LK_matrix) 
 OGP_L_K_fit1LK$Linf = as.numeric(as.character(OGP_L_K_fit1LK$Linf))
 OGP_L_K_fit1LK$K = as.numeric(as.character(OGP_L_K_fit1LK$K))
 OGP_L_K_fit1LK$OGP = as.numeric(as.character(OGP_L_K_fit1LK$OGP))
 OGP_L_K_fit1LK$min.Linf = as.numeric(as.character(OGP_L_K_fit1LK$min.Linf))
 OGP_L_K_fit1LK$max.Linf = as.numeric(as.character(OGP_L_K_fit1LK$max.Linf))
 OGP_L_K_fit1LK$min.K = as.numeric(as.character(OGP_L_K_fit1LK$min.K))
 OGP_L_K_fit1LK$max.K = as.numeric(as.character(OGP_L_K_fit1LK$max.K))
 OGP_L_K_fit1LK$OGP.min = as.numeric(as.character(OGP_L_K_fit1LK$OGP.min))
 OGP_L_K_fit1LK$OGP.max = as.numeric(as.character(OGP_L_K_fit1LK$OGP.max))
 OGP_L_K_fit1LK
 
 OGP_L_K_fit1LK_order_OGP <- OGP_L_K_fit1LK %>% mutate(site = fct_reorder(site, OGP))
 
 ggplot(data = OGP_L_K_fit1LK_order_OGP) + 
   geom_pointrange(aes(y = site, x = OGP, xmin = OGP.min, xmax = OGP.max)) + 
   geom_errorbar(aes(y = site, x = OGP, xmin = OGP.min, xmax = OGP.max))
      # not statistically significant (error bar overlap)
 ggsave(here::here("Output", "OGP_vs_site.png"))
 
   
 OGP_L_K_fit1LK_order_K <- OGP_L_K_fit1LK %>% mutate(site = fct_reorder(site, K))
 
 ggplot(data = OGP_L_K_fit1LK_order_K) + 
   geom_pointrange(aes(y = site, x = K, xmin = min.K, xmax = max.K)) + 
   geom_errorbar(aes(y = site, x = K, xmin = min.K, xmax = max.K))
      # not statistically significant (error bar overlap)
 ggsave(here::here("Output", "K_vs_site.png"))
 

 OGP_L_K_fit1LK_order_Linf <- OGP_L_K_fit1LK %>% mutate(site = fct_reorder(site, Linf))
 
 ggplot(data = OGP_L_K_fit1LK_order_Linf) + 
   geom_pointrange(aes(y = site, x = Linf, xmin = min.Linf, xmax = max.Linf)) + 
   geom_errorbar(aes(y = site, x = Linf, xmin = min.Linf, xmax = max.Linf))
     # barely any error bar overlap (btwn two groups)
 ggsave(here::here("Output", "Linf_vs_site.png"))
 
 

 
 
 
 
#==== PLEASE IGNORE ALL CODE BELOW, I KEPT IT INCASE WE NEED A TEMPLATE OR FOR REFERENCE (such as with the test code), DELETE FOR FINAL ====
 
 
 
 
 
 
 
# levene and/or ANOVA trials
         # Linf_stat = 153.00
         # K = 
         # leveneTest(length~Linf_stat*(1-exp(-K[site]*(ages-t0))), data = clam_data_clean)
 
 # L_site_aov <- aov(length ~ site, data = clam_data_clean)
 # L_site_aov_summary <- summary(L_site_aov)
 # L_site_aov_Tukey <- TukeyHSD(L_site_aov)
 # L_site_aov_summary
 #   sink(here::here("Output", "L_site_summary.txt"))
 #   print(L_site_aov_summary)
 #   sink()
 # L_site_aov_Tukey
 #   sink(here::here("Output", "L_site_aov_Tukey.txt"))
 #   print(L_site_aov_Tukey)
 #   sink() 
 # 
 #spacer
  ### for 1KT
#  overview(fit1KT) #Note that t0 is not significantly different from 0 - the confidence intervals overlap with 0
#  summary(fit1KT)
#  coef(fit1KT)
# Linf1 <- coef(fit1KT)[1]
# min.Linf1 <- confint.default(fit1KT)[1,1]
# max.Linf1 <- confint.default(fit1KT)[1,2]
# 
# K1 <- coef(fit1KT)[2]
# K2 <- coef(fit1KT)[3]
# K3 <- coef(fit1KT)[4]
# K4 <- coef(fit1KT)[5]
# K5 <- coef(fit1KT)[6]
# min.K1 <- confint.default(fit1KT)[2,1]
# max.K1 <- confint.default(fit1KT)[2,2]
# min.K2 <- confint.default(fit1KT)[3,1]
# max.K2 <- confint.default(fit1KT)[3,2]
# min.K3 <- confint.default(fit1KT)[4,1]
# max.K3 <- confint.default(fit1KT)[4,2]
# min.K4 <- confint.default(fit1KT)[5,1]
# max.K4 <- confint.default(fit1KT)[5,2]
# min.K5 <- confint.default(fit1KT)[6,1]
# max.K5 <- confint.default(fit1KT)[6,2]
# 
# t1 <- coef(fit1KT)[7] 
# t2 <- coef(fit1KT)[8]
# t3 <- coef(fit1KT)[9]
# t4 <- coef(fit1KT)[10]
# t5 <- coef(fit1KT)[11]
# min.t1 <- confint.default(fit1KT)[7,1]
# max.t1 <- confint.default(fit1KT)[7,2]
# min.t2 <- confint.default(fit1KT)[8,1]
# max.t2 <- confint.default(fit1KT)[8,2]
# min.t3 <- confint.default(fit1KT)[9,1]
# max.t3 <- confint.default(fit1KT)[9,2]
# min.t4 <- confint.default(fit1KT)[10,1]
# max.t4 <- confint.default(fit1KT)[10,2]
# min.t5 <- confint.default(fit1KT)[11,1]
# max.t5 <- confint.default(fit1KT)[11,2]
# 
# (OGP1 <- as.numeric(log10(K1)+2*log10(0.1*Linf1)))
# (OGP2 <- as.numeric(log10(K2)+2*log10(0.1*Linf1)))
# (OGP3 <- as.numeric(log10(K3)+2*log10(0.1*Linf1)))
# (OGP4 <- as.numeric(log10(K4)+2*log10(0.1*Linf1)))
# (OGP5 <- as.numeric(log10(K5)+2*log10(0.1*Linf1)))
# (OGP1.min <- as.numeric(log10(min.K1)+2*log10(0.1*min.Linf1)))
# (OGP1.max <- as.numeric(log10(max.K1)+2*log10(0.1*max.Linf1)))
# (OGP2.min <- as.numeric(log10(min.K2)+2*log10(0.1*min.Linf1)))
# (OGP2.max <- as.numeric(log10(max.K2)+2*log10(0.1*max.Linf1)))
# (OGP3.min <- as.numeric(log10(min.K3)+2*log10(0.1*min.Linf1)))
# (OGP3.max <- as.numeric(log10(max.K3)+2*log10(0.1*max.Linf1)))
# (OGP4.min <- as.numeric(log10(min.K4)+2*log10(0.1*min.Linf1)))
# (OGP4.max <- as.numeric(log10(max.K4)+2*log10(0.1*max.Linf1)))
# (OGP5.min <- as.numeric(log10(min.K5)+2*log10(0.1*min.Linf1)))
# (OGP5.max <- as.numeric(log10(max.K5)+2*log10(0.1*max.Linf1)))
# 

###---- Test Code (uses test_data) ----
# test_data <- read.csv(here::here("Input_Data", "test_data_size_at_age_MAG.csv"))
# # 
# head(test_data)
# test_data$Site <-as.factor(test_data$Site)
# str(test_data)
# 
# # ===Functions
# (svCom <- vbStarts(Lengths~Ages, data = test_data))           # base model
# (svGen <- lapply(svCom, rep, 5))
# vbGen <- Lengths~Linf[Site]*(1-exp(-K[Site]*(Ages-t0)))
# fitGen <- nls(vbGen, data = test_data, start = svGen)
# hist(residuals(fitGen), main = "")
# 
# # vb1KT <- Lengths~Linf*(1-exp(-K[Site]*(Ages-t0[Site])))
# # sv1KT <- mapply(rep,svCom,c(1,2,2))
# # vb1LT <- Lengths~Linf[Site]*(1-exp(-K*(Ages-t0[Site])))
# # sv1LT <- mapply(rep,svCom,c(2,1,2))
# vb1LK <- Lengths~Linf[Site]*(1-exp(-K[Site]*(Ages-t0))) # make sense
# sv1LK <- mapply(rep,svCom,c(5,5,1))
# # vb2T <- Lengths~Linf*(1-exp(-K*(Ages-t0[Site])))
# # sv2T <- mapply(rep,svCom,c(1,1,2))
# # vb2K <- Lengths~Linf*(1-exp(-K[Site]*(Ages-t0)))
# # sv2K <- mapply(rep,svCom,c(1,2,1))
# vb2L <- Lengths~Linf[Site]*(1-exp(-K*(Ages-t0)))
# sv2L <- mapply(rep,svCom,c(5,1,1))
# vbCom <- Lengths~Linf*(1-exp(-K*(Ages-t0))) # make sense
# vbComT0 <- Lengths~Linf*(1-exp(-K*(Ages))) # make sense
# vb1LKT0 <- Lengths~Linf[Site]*(1-exp(-K[Site]*(Ages))) # make sense
# 
# 
#  # fit1KT <- nls(vb1KT,data=test_data,start=sv1KT)
#  # fit1LT <- nls(vb1LT,data=test_data,start=sv1LT)
#  fit1LK <- nls(vb1LK,data=test_data,start=sv1LK) # make sense
#  # fit2T <- nls(vb2T,data=test_data,start=sv2T)
#  # fit2K <- nls(vb2K,data=test_data,start=sv2K)
#  fit2L <- nls(vb2L,data=test_data,start=sv2L)
#  fitCom <- nls(vbCom,data=test_data,start=svCom) # make sense
#  fitComT0 <- nls(vbComT0,data=test_data,start=svCom[1:2]) # make sense
#  fit1LKT0 <- nls(vb1LKT0,data=test_data,start=sv1LK[1:2]) # make sense
# 
# # ===ANOVAs and AIC
# 
#  anova(fit1LK,fitCom) # Since Linf and K are correlated it doesn't make sense for
#  # there to only be a difference in one of the parameters between the sites
#  AIC(fit1LK,fitCom)
# 
# 
# 
#  #AIC(fitGen,fit1KT,fit1LT,fit1LK,fit2T,fit2K,fit2L,fitCom,fitComT0,fit1LKT0) # Still within 2 AIC
# 
# 
# 
#  # Fit 2 L's only
#      # overview(fit2L)
#      # coef(fit2L)[-2]
#      # plot(Lengths~jitter(Ages,0.3),data=test_data,subset=Site=="Barn",pch=19,xlab="Age (yrs)",
#      #      ylab="Total Length (mm)", ylim = c(0,200), xlim = c(0,20))
#      # points(Lengths~jitter(Ages,0.3),data=test_data,subset=Site=="Eel",pch=19,col="gray")
#      # vbTypical <- vbFuns("typical")
#      # curve(vbTypical(x,Linf=coef(fit2L)[-2]),from=2,to=14,lwd=2,add=TRUE)
#      # curve(vbTypical(x,Linf=coef(fit2L)[-1]),from=2,to=8,col="gray",lwd=2,add=TRUE)
#      # legend("topleft",legend=c("Barnstable","Eel"),col=c("black","gray"),lwd=2,lty=1,cex=0.75)
#  overview(fit1LK)
#  coef(fit1LK)
# 
#  levels(test_data$Site)
# 
#  plot(Lengths~jitter(Ages,0.3),data=test_data,subset=Site=="Barn",pch=19,xlab="Age (yrs)",
#       ylab="Total Length (mm)", ylim = c(0,200), xlim = c(0,20))
#  points(Lengths~jitter(Ages,0.3),data=test_data,subset=Site=="Chatham",pch=19,col="red")
#  points(Lengths~jitter(Ages,0.3),data=test_data,subset=Site=="Dennis",pch=19,col="gray")
#  points(Lengths~jitter(Ages,0.3),data=test_data,subset=Site=="Eel",pch=19,col="blue")
#  points(Lengths~jitter(Ages,0.3),data=test_data,subset=Site=="Provincetown",pch=19,col="green")
#  vbTypical <- vbFuns("typical")
# 
#  coef1 <- c(coef(fit1LK)[1], coef(fit1LK)[6], coef(fit1LK)[11])
#  coef2 <- c(coef(fit1LK)[2], coef(fit1LK)[7], coef(fit1LK)[11])
#  coef3 <- c(coef(fit1LK)[3], coef(fit1LK)[8], coef(fit1LK)[11])
#  coef4 <- c(coef(fit1LK)[4], coef(fit1LK)[9], coef(fit1LK)[11])
#  coef5 <- c(coef(fit1LK)[5], coef(fit1LK)[10], coef(fit1LK)[11])
# 
#  curve(vbTypical(x,Linf=coef1),from=0,to=20,lwd=2,add=TRUE)
#  curve(vbTypical(x,Linf=coef2),from=0,to=20,col="red",lwd=2,add=TRUE)
#  curve(vbTypical(x,Linf=coef3),from=0,to=20,col="grey",lwd=2,add=TRUE)
#  curve(vbTypical(x,Linf=coef4),from=0,to=20,col="blue",lwd=2,add=TRUE)
#  curve(vbTypical(x,Linf=coef5),from=0,to=20,col="green",lwd=2,add=TRUE)
# 
#  legend("bottomright",legend=c("Barnstable","Chatham", "Dennis", "Eel", "Provincetown"),col=c("black","gray", "red", "blue", "green"),lwd=2,lty=1,cex=0.75)
# 
# 
#  # Reference for reasoning to include T0 estimates (but probably a composite from all data is best
#  # rather than site specific if only 20 data points per site and no juveniles): https://besjournals.onlinelibrary.wiley.com/doi/pdf/10.1111/2041-210x.12020
# 
#  # OGP reference: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0199212
# 
#  overview(fit1LK) #Note that t0 is not significantly different from 0 - the confidence intervals overlap with 0
#  coef(fit1LK)
# 
#  Linf1 <- coef(fit1LK)[1]
#  Linf2 <- coef(fit1LK)[2]
#  min.Linf1 <- confint.default(fit1LK)[1,1]
#  max.Linf1 <- confint.default(fit1LK)[1,2]
#  min.Linf2 <- confint(fit1LK)[2,1]
#  max.Linf2 <- confint(fit1LK)[2,2]
#  K1 <- coef(fit1LK)[3]
#  K2 <- coef(fit1LK)[4]
#  min.K1 <- confint(fit1LK)[3,1]
#  max.K1 <- confint(fit1LK)[3,2]
#  min.K2 <- confint(fit1LK)[4,1]
#  max.K2 <- confint(fit1LK)[4,2]
#  (OGP1 <- as.numeric(log10(K1)+2*log10(0.1*Linf1)))
#  (OGP2 <- as.numeric(log10(K2)+2*log10(0.1*Linf2)))
# 
#  (OGP1.min <- as.numeric(log10(min.K1)+2*log10(0.1*min.Linf1)))
#  (OGP1.max <- as.numeric(log10(max.K1)+2*log10(0.1*max.Linf1)))
#  (OGP2.min <- as.numeric(log10(min.K2)+2*log10(0.1*min.Linf2)))
#  (OGP2.max <- as.numeric(log10(max.K2)+2*log10(0.1*max.Linf2)))



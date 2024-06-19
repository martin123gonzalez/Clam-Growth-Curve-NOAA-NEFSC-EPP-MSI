# File Name: size_at_age_CH_to_SL
# Authors: Molly Roberts, Martin Gonzalez
# Purpose: To graph von Bertalanffy growhth function of modeled and transformed 
         # shell lengths using incremental growth data


# setwd("~/Documents/GitHub/EAD-ASEB-Ssolidissima-OA/data")



#Note to Molly and Matt: uses the same basic framework as "Sample_size_at_age_code.R" but using different data (which you'll updata).
    # Here, I used incremental shell growth we discussed at the end of my internship. 
    # As you'll see in the main graph, data is skewed by a very large barnstable that seems to have a very low growth rate.
    # But as Matt found, it was an error on my data entry and should be corrected using your new dataset.
    # I also did not have time to reproduce the box and whisker type plot (geom_pointrange) using this new data.
    # but it should be very straightforward to reproduce.



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
clam_SL <- read.csv(here::here("Input_Data", "Clam_metric_chondro_and_shell_growth_per_year.csv"))

colnames(clam_SL)[8] = "length"
colnames(clam_SL)[3] = "ages"

clam_SL_clean <- clam_SL %>% 
  dplyr::filter(ages != "NA" & length != "NA") %>% 
  dplyr::filter(site != "?" & site != "Nobscusset") %>% 
  dplyr::select(site, length, ages, Param, Value)

clam_SL_clean$site <-as.factor(clam_SL_clean$site)

# ==== model set up ====
svCom_sl <- vbStarts(length~ages, data = clam_SL_clean)
(svGen_sl <- lapply(svCom_sl, rep, 5))
vbGen_sl <- length~Linf[site]*(1-exp(-K[site]*(ages-t0)))
fitGen_sl <- nls(vbGen_sl, data = clam_SL_clean, start = svGen_sl)
hist(residuals(fitGen_sl), main = "")



vb1LK_sl <- length~Linf[site]*(1-exp(-K[site]*(ages-t0))) # make sense
sv1LK_sl <- mapply(rep,svCom_sl,c(5,5,1))
vbCom_sl <- length~Linf*(1-exp(-K*(ages-t0))) # make sense
vbComT0_sl <- length~Linf*(1-exp(-K*(ages))) # make sense
vb1LKT0_sl <- length~Linf[site]*(1-exp(-K[site]*(ages))) # make sense
vb1KT_sl <- length~Linf*(1-exp(-K[site]*(ages-t0[site])))
sv1KT_sl <- mapply(rep,svCom_sl,c(1,5,5))
vb1LT_sl <- length~Linf[site]*(1-exp(-K*(ages-t0[site])))
sv1LT_sl <- mapply(rep,svCom_sl,c(5,1,5))
vb2L_sl <- length~Linf[site]*(1-exp(-K*(ages-t0)))
sv2L_sl <- mapply(rep,svCom_sl,c(5,1,1))
vb2T_sl <- length~Linf*(1-exp(-K*(ages-t0[site])))
sv2T_sl <- mapply(rep,svCom_sl,c(1,1,5))
vb2K_sl <- length~Linf*(1-exp(-K[site]*(ages-t0)))
sv2K_sl <- mapply(rep,svCom_sl,c(1,5,1))


fit1KT_sl <- nls(vb1KT_sl,data=clam_SL_clean,start=sv1KT_sl)
fit1LT_sl <- nls(vb1LT_sl,data=clam_SL_clean,start=sv1LT_sl)
fit1LK_sl <- nls(vb1LK_sl,data = clam_SL_clean,start = sv1LK_sl) # make sense
fit2T_sl <- nls(vb2T_sl,data=clam_SL_clean,start=sv2T_sl)
fit2K_sl <- nls(vb2K_sl,data=clam_SL_clean,start=sv2K_sl)
fit2L_sl <- nls(vb2L_sl,data=clam_SL_clean,start=sv2L_sl)
fitCom_sl <- nls(vbCom_sl,data=clam_SL_clean, start=svCom_sl) # make sense
fitComT0_sl <- nls(vbComT0_sl,data=clam_SL_clean,start=svCom_sl[1:2]) # make sense
fit1LKT0_sl <- nls(vb1LKT0_sl,data=clam_SL_clean,start=sv1LK_sl[1:2]) # make sense


anova(fit1LK_sl,fitCom_sl) # Since Linf and K are correlated it doesn't make sense for 
# there to only be a difference in one of the parameters between the sites
AIC(fit1LK_sl,fitCom_sl, fitComT0_sl, fit1LKT0_sl)

AIC(fitGen_sl,fit1KT_sl,fit1LT_sl,fit1LK_sl,fit2T_sl,fit2K_sl,fit2L_sl,fitCom_sl,fitComT0_sl,fit1LKT0_sl) # Still within 2 AIC
# Winner with context = 1LK or 1LKT0
vbTypical <- vbFuns("typical")

overview(fit1LK_sl)


# ==== plotting ====
# png("Output/surfclam_sVBGF.png", width = 600, height = 600)       #un-comment this line to save image
levels(clam_SL_clean$site)

plot(length~jitter(ages,0.8),data=clam_SL_clean,subset=site=="N. Cape: Barnstable",pch=19,xlab="Age (yrs)",
     ylab="Total Length (mm)", ylim = c(0,200), xlim = c(0,15))
points(length~jitter(ages,0.8),data = clam_SL_clean,subset = site == "S. Cape: Chatham",pch=19,col="red")
points(length~jitter(ages,0.8),data = clam_SL_clean,subset = site == "N. Cape: Dennis",pch=19,col="orange")
points(length~jitter(ages,0.8),data = clam_SL_clean,subset = site == "S. Cape: Eel Pond",pch=19,col="blue")
points(length~jitter(ages,0.8),data = clam_SL_clean,subset = site == "N. Cape: Provincetown",pch=19,col="green")
# points(length~jitter(ages,0.8), data = clam_data_clean, pch = 19, col = "grey")
vbTypical <- vbFuns("typical")

# for fit1LK
overview(fit1LK_sl)
coef(fit1LK_sl)
coef1 <- c(coef(fit1LK_sl)[1], coef(fit1LK_sl)[6], coef(fit1LK_sl)[11])  #Barn
coef2 <- c(coef(fit1LK_sl)[2], coef(fit1LK_sl)[7], coef(fit1LK_sl)[11])  #Den
coef3 <- c(coef(fit1LK_sl)[3], coef(fit1LK_sl)[8], coef(fit1LK_sl)[11])  #Ptown
coef4 <- c(coef(fit1LK_sl)[4], coef(fit1LK_sl)[9], coef(fit1LK_sl)[11])  #Chat
coef5 <- c(coef(fit1LK_sl)[5], coef(fit1LK_sl)[10], coef(fit1LK_sl)[11]) #Eel
coef_avg <- c(mean(coef(fit1LK_sl)[1:5]), mean(coef(fit1LK_sl)[6:10]), coef(fit1LK_sl)[11])

curve(vbTypical(x,Linf=coef1),from=0,to=20,lwd=2,add=TRUE)
curve(vbTypical(x,Linf=coef2),from=0,to=20,col="orange",lwd=2,add=TRUE)
curve(vbTypical(x,Linf=coef3),from=0,to=20,col="green",lwd=2,add=TRUE)
curve(vbTypical(x,Linf=coef4),from=0,to=20,col="red",lwd=2,add=TRUE)
curve(vbTypical(x,Linf=coef5),from=0,to=20,col="blue",lwd=2,add=TRUE)
curve(vbTypical(x,Linf=coef_avg),from = 0, to = 20, col = "grey", lwd = 2, add = TRUE)

legend("bottomright",legend=c("Average", "N. Cape: Barnstable", "N. Cape: Dennis", "N. Cape: Provincetown", "S. Cape: Chatham", "S. Cape: Eel Pond"),
       col=c("grey", "black","orange", "green", "red", "blue"),lwd=2,lty=1,cex=.75, title = "Sites in Cape Cod")

title("Surfclam Growth Rate over Time, using incremental shell growth data")
# dev.off()                     #un-comment this line to save image




overview(fit1LK_sl) #Note that t0 is not significantly different from 0 - the confidence intervals overlap with 0
summary(fit1LK_sl)
coef(fit1LK_sl)

# Linf Calc
Linf1_sl <- coef(fit1LK_sl)[1]
Linf2_sl <- coef(fit1LK_sl)[2]
Linf3_sl <- coef(fit1LK_sl)[3]
Linf4_sl <- coef(fit1LK_sl)[4]
Linf5_sl <- coef(fit1LK_sl)[5]

min.Linf1_sl <- confint.default(fit1LK_sl)[1,1]
max.Linf1_sl <- confint.default(fit1LK_sl)[1,2]
min.Linf2_sl <- confint.default(fit1LK_sl)[2,1]
max.Linf2_sl <- confint.default(fit1LK_sl)[2,2]
min.Linf3_sl <- confint.default(fit1LK_sl)[3,1]
max.Linf3_sl <- confint.default(fit1LK_sl)[3,2]
min.Linf4_sl <- confint.default(fit1LK_sl)[4,1]
max.Linf4_sl <- confint.default(fit1LK_sl)[4,2]
min.Linf5_sl <- confint.default(fit1LK_sl)[5,1]
max.Linf5_sl <- confint.default(fit1LK_sl)[5,2]

# K calc
K1_sl <- coef(fit1LK_sl)[6]
K2_sl <- coef(fit1LK_sl)[7]
K3_sl <- coef(fit1LK_sl)[8]
K4_sl <- coef(fit1LK_sl)[9]
K5_sl <- coef(fit1LK_sl)[10]

min.K1_sl <- confint.default(fit1LK_sl)[6,1]
max.K1_sl <- confint.default(fit1LK_sl)[6,2]
min.K2_sl <- confint.default(fit1LK_sl)[7,1]
max.K2_sl <- confint.default(fit1LK_sl)[7,2]
min.K3_sl <- confint.default(fit1LK_sl)[8,1]
max.K3_sl <- confint.default(fit1LK_sl)[8,2]
min.K4_sl <- confint.default(fit1LK_sl)[9,1]
max.K4_sl <- confint.default(fit1LK_sl)[9,2]
min.K5_sl <- confint.default(fit1LK_sl)[10,1]
max.K5_sl <- confint.default(fit1LK_sl)[10,2]


L_K_fit1LK_matrix_sl <- matrix(c("N. Cape: Barnstable", "N. Cape: Dennis", "N. Cape: Provincetown", "S. Cape: Chatham", "S. Cape: Eel Pond",
                                  Linf1_sl, Linf2_sl, Linf3_sl, Linf4_sl, Linf5_sl, 
                                  K1_sl, K2_sl, K3_sl, K4_sl, K5_sl,
                                  min.Linf1_sl, min.Linf2_sl, min.Linf3_sl, min.Linf4_sl, min.Linf5_sl,
                                  max.Linf1_sl, max.Linf2_sl, max.Linf3_sl, max.Linf4_sl, max.Linf5_sl,
                                  min.K1_sl, min.K2_sl, min.K3_sl, min.K4_sl, min.K5_sl,
                                  max.K1_sl, max.K2_sl, max.K3_sl, max.K4_sl, max.K5_sl),
                                  nrow = 5, ncol = 7, byrow = F,
                                  dimnames = list(c(1, 2, 3, 4, 5), c("site", "Linf", "K","min.Linf",
                                                                    "max.Linf", "min.K", "max.K")))
L_K_fit1LK_sl <- data.frame(L_K_fit1LK_matrix_sl) 
L_K_fit1LK_sl$Linf = as.numeric(as.character(L_K_fit1LK_sl$Linf))
L_K_fit1LK_sl$K = as.numeric(as.character(L_K_fit1LK_sl$K))
L_K_fit1LK_sl$min.Linf = as.numeric(as.character(L_K_fit1LK_sl$min.Linf))
L_K_fit1LK_sl$max.Linf = as.numeric(as.character(L_K_fit1LK_sl$max.Linf))
L_K_fit1LK_sl$min.K = as.numeric(as.character(L_K_fit1LK_sl$min.K))
L_K_fit1LK_sl$max.K = as.numeric(as.character(L_K_fit1LK_sl$max.K))
L_K_fit1LK_sl

L_K_fit1LK_sl_order_Linf <- L_K_fit1LK_sl %>% mutate(site = fct_reorder(site, Linf))

# ggplot(data = L_K_fit1LK_sl, aes(y = site, x = Linf_sl, xmin = min.Linf, xmax = max.Linf))          #unfinished, tried to plot 

# ggplot(data = L_K_fit1LK_sl) +                        #unfinished
#   geom_pointrange(aes(y = site, x = Linf, xmin = min.Linf, xmax = max.Linf)) + 
#   geom_errorbar(aes(y = site, x = Linf, xmin = min.Linf, xmax = max.Linf))


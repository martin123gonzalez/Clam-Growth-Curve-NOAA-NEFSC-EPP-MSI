# File Name: chondrophore_dist.R
# Authors: Martin Gonzalez
# Purpose: To graph 
         # 1) a dot plot of the first three years of growth for each site and 
         # 2) the change in distance (measured as cumulative chondorphore growth) for each site

#Note to Molly and Matt: This code needs updated data. 

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
library(stringr)
# ==== data input ====
clamshell_data <- read.csv(here::here("Input_Data", "Clam_metric_chondro_and_shell_growth_per_year.csv"))

clamshell_clean <- clamshell_data %>% 
  dplyr::filter(site != "?" & site != "Nobscusset") %>% 
  dplyr::filter(age >= 1) 
  # mutate(shell_cumul_dist = 18.58 + 9.85 * chondro_cumul_dist + -0.14 * (chondro_cumul_dist)^2)

# ==== graphing ====
  #indiv chondro growth dotplots
clamshell_y1 <- clamshell_clean %>% 
  dplyr::filter(growth_year_range == "1-2")
ggplot(data = clamshell_y1, aes(x = chondro_dist_addition)) + 
  geom_dotplot(aes(fill = site)) +
  xlim(0, 11) + 
  facet_wrap(~site)

clamshell_y2 <- clamshell_clean %>% 
  dplyr::filter(growth_year_range == "2-3")
ggplot(data = clamshell_y2, aes(x = chondro_dist_addition)) + 
  geom_dotplot(aes(fill = site)) +
  xlim(0, 11) + 
  facet_wrap(~site)

clamshell_y3 <- clamshell_clean %>% 
  dplyr::filter(growth_year_range == "3-4")
ggplot(data = clamshell_y3, aes(x = chondro_dist_addition)) + 
  geom_dotplot(aes(fill = site)) +
  xlim(0, 11) + 
  facet_wrap(~site)

  # !!!
  # y1,2,3 chondro growth periods
clamshell_y123 <- clamshell_clean %>% 
  dplyr::filter(growth_year_range == "1-2" | growth_year_range == "2-3" | growth_year_range == "3-4") 

ggplot(data = clamshell_y123, aes(x = chondro_dist_addition)) + 
  geom_dotplot(aes(fill = growth_year_range)) +
  facet_wrap(~site) + 
  labs(title = "Chondrophore Growth for first 3 complete years", fill = "Growth Period (years)",
       x = "Distance between growth lines (mm)")
ggsave(here::here("Output", "y123_chondrophore_growth.png"), width = 12, height = 8)

  # !!!
  #chondro cumulative growth, growth year vs distance

ggplot(clamshell_clean, aes(x = growth_period, y = chondro_cumul_dist)) +
  geom_point() +
  geom_smooth(method = "lm",formula = y ~ log(x), se = T) +
  facet_wrap(~site) + 
  labs(x ="Growth Period (measured as distance between two yearly growth bands)", y = "Cumulative Chondrophore Distance (mm)",
       title = "Cumulative Chondrophore Distance over time")
ggsave(here::here("Output", "cumul_chondro_dist_over_time.png"), width = 12, height = 8)





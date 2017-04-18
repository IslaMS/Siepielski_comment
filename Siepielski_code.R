# Reanalysis of Siepielski et al. Science 2017 "Precipitation drives global variation in natural selection"

# Code written by Isla Myers-Smith, University of Edinburgh (with statistical advice from Jarrod Hadfield, University of Edinburgh)

# 13 Apr 2017

# Libraries ----
library(lme4)
library(MCMCglmm)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(coda)

# Generic functions ----

theme_tidy <- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 12), 
          axis.title = element_text(size = 14),
          axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), units = , "cm"),
          plot.title = element_text(size=14, vjust=1, hjust=0.5),
          legend.text = element_text(size=12, face="italic"),          
          legend.title = element_blank(),                              
          legend.position = c(0.9, 0.9), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", fill = "transparent", size = 2, linetype="blank"))
}

# Load data ----
# Edit cells 5833:5862, AN-A0 as non-ASCII text
dat1 <- read.csv("Database-S1.csv", na.strings=c("", "na"))

# p.value entered into the wrong column?
dat1[grep("<", dat1$Grad.linear.value), "Grad.linear.value"] <- NA
dat1[grep("<", dat1$Diff.linear.value), "Diff.linear.value"] <- NA
dat1$Grad.linear.value <- as.numeric(as.character(dat1$Grad.linear.value))
dat1$Diff.linear.value <- as.numeric(as.character(dat1$Diff.linear.value))

dat2 <- read.csv("Database-S2.csv", na.strings=c("", "na"))

# Merge dat1 and dat2
dat <- merge(dat1, dat2, by.x = c("Final.Index", "Study.ID", "Year"), by.y = c("Final.Index", "Study.ID", "Year"))

# Hack because lmer *thinks* the model is unidentifiable
dat$lmer.units <- as.factor(1:nrow(dat))
dat$lmer.units[2] <- "1"

# Order data.farme by study/trait
dat$study.trait <- as.factor(paste(dat$Study.ID, dat$Trait))
dat <- dat[order(dat$study.trait),]

# Remove rows wth misisng data
dat.grad <- subset(dat, !is.na(Grad.linear.value) & !is.na(Grad.linear.StErr) & !is.na(AvgTempMean) & !is.na(SDTempMean) & !is.na(MinTempMin) & !is.na(MaxTempMax) & !is.na(AvgPrecip) & !is.na(MaxPrecip) & !is.na(MinPrecip) & !is.na(SDPrecip))
dat.diff <- subset(dat, !is.na(Diff.linear.value) & !is.na(Diff.linear.StErr) & !is.na(AvgTempMean) & !is.na(SDTempMean) & !is.na(MinTempMin) & !is.na(MaxTempMax) & !is.na(AvgPrecip) & !is.na(MaxPrecip) & !is.na(MinPrecip) & !is.na(SDPrecip))

# Standardise climate within Study/Traits ----

STD.grad <- function(x){unlist(tapply(x, dat.grad$study.trait, function(x){if(length(x)>1 & var(x)>0){scale(x)}else{rep(0, length(x))}}))}
STD.diff <- function(x){unlist(tapply(x, dat.diff$study.trait, function(x){if(length(x)>1 & var(x)>0){scale(x)}else{rep(0, length(x))}}))}

dat.PET.novariance <- dat %>% group_by(Study.ID, study.trait) %>% mutate(range = (max(AvgPET)-min(AvgPET))) %>% ungroup() %>% distinct(Study.ID, range) %>% subset(range == 0) %>% arrange(Study.ID)

# Plot figures of PET data by Study.ID and study.trait
dat %>% group_by(Study.ID, study.trait) %>% do(ggsave(ggplot(.,aes(x = AvgPET, y = Grad.linear.value)) + geom_point(colour="darkblue") + theme_tidy() + xlab("Mean PET"), filename = gsub("", "", paste("figures/PET/", unique(as.character(.$Study.ID)),"_",unique(as.character(.$study.trait)), ".pdf", sep="")), device="pdf"))

no.var.PET.study.number <- length(unique(dat.PET.novariance$Study.ID))
total.study.number <- length(unique(dat$Study.ID))

# Apply standardising fuction to all climate variables with variance
data.grad <- dat.grad %>% mutate(MinTempMinStd = STD.grad(MinTempMin), MaxTempMaxStd = STD.grad(MaxTempMax), AvgTempMeanStd = STD.grad(AvgTempMean), SDTempMeanStd = STD.grad(SDTempMean), MinPrecipStd = STD.grad(MinPrecip), MaxPrecipStd = STD.grad(MaxPrecip), AvgPrecipStd = STD.grad(AvgPrecip), SDPrecipStd = STD.grad(SDPrecip))

data.diff <- dat.diff %>% mutate(MinTempMinStd = STD.diff(MinTempMin), MaxTempMaxStd = STD.diff(MaxTempMax), AvgTempMeanStd = STD.diff(AvgTempMean), SDTempMeanStd = STD.diff(SDTempMean), MinPrecipStd = STD.diff(MinPrecip), MaxPrecipStd = STD.diff(MaxPrecip), AvgPrecipStd = STD.diff(AvgPrecip), SDPrecipStd = STD.diff(SDPrecip))

data.grad <- arrange(data.grad, Final.Index)
data.diff <- arrange(data.diff, Final.Index)

# Model all data ----
# Select relevant data
grad1 <- select(data.grad, Final.Index, Study.ID, Population, Taxon.group, Trait.Class, study.trait, Long.Decimal.Degrees, Lat.Decimal.Degrees, Grad.linear.value, Grad.linear.StErr, MinTempMin, MaxTempMax, AvgTempMean, SDTempMean, MaxPrecip, MinPrecip, AvgPrecip, SDPrecip, MinTempMinStd, MaxTempMaxStd, AvgTempMeanStd, SDTempMeanStd, MaxPrecipStd, MinPrecipStd, AvgPrecipStd, SDPrecipStd)

diff1 <- select(data.diff, Final.Index, Study.ID, Population, Taxon.group, Trait.Class, study.trait, Long.Decimal.Degrees, Lat.Decimal.Degrees, Diff.linear.value, Diff.linear.StErr, MinTempMin, MaxTempMax, AvgTempMean, SDTempMean, MaxPrecip, MinPrecip, AvgPrecip, SDPrecip, MinTempMinStd, MaxTempMaxStd, AvgTempMeanStd, SDTempMeanStd, MaxPrecipStd, MinPrecipStd, AvgPrecipStd, SDPrecipStd)

# Reshape to long format
grad <- gather(grad1, climatevar, climvar.value, MinTempMinStd, MaxTempMaxStd, AvgTempMeanStd, SDTempMeanStd, MaxPrecipStd, MinPrecipStd, AvgPrecipStd, SDPrecipStd)

diff <- gather(diff1, climatevar, climvar.value, MinTempMinStd, MaxTempMaxStd, AvgTempMeanStd, SDTempMeanStd, MaxPrecipStd, MinPrecipStd, AvgPrecipStd, SDPrecipStd)

# Functions to fit climate models ----

# Parameter expanded priors
prior <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = c(0, 0), alpha.V = diag(2)*1000)))

fit.grad <- function(grad) {
  fitmod <- MCMCglmm(Grad.linear.value~climvar.value, random=~us(1+climvar.value):study.trait, mev=grad$Grad.linear.StErr^2, data=grad, prior=prior, family = "gaussian", pr=TRUE, nitt = 30000, burnin = 2000)
  sum <- data.frame(rbind("Intercept", "climvar.value"), summary(fitmod)[5][[1]])
  colnames(sum) <- c("estimate", "post.mean", "l.95.CI", "u.95.CI", "eff.samp", "pMCMC")
  VCVsum <- data.frame(cbind("climvar.value", as.character(grad$climatevar[1]), as.character(grad$Taxon.group[1]), mean((fitmod$VCV[,4]+fitmod$Sol[,2]^2)/(rowSums(fitmod$VCV[,c(4,6)])+fitmod$Sol[,2]^2)), HPDinterval((fitmod$VCV[,4]+fitmod$Sol[,2]^2)/(rowSums(fitmod$VCV[,c(4,6)])+fitmod$Sol[,2]^2))[1], HPDinterval((fitmod$VCV[,4]+fitmod$Sol[,2]^2)/(rowSums(fitmod$VCV[,c(4,6)])+fitmod$Sol[,2]^2))[2]))
  colnames(VCVsum) <- c("estimate", "climatevar", "Taxon.group", "meanVCV", "VCVlower", "VCVupper")
  merge(sum, VCVsum)
}

fit.diff <- function(diff) {
  fitmod <- MCMCglmm(Diff.linear.value~climvar.value, random=~us(1+climvar.value):study.trait, mev=diff$Diff.linear.StErr^2, data=diff, prior=prior, family = "gaussian", pr=TRUE, nitt = 30000, burnin = 2000)
  sum <- data.frame(rbind("Intercept", "climvar.value"), summary(fitmod)[5][[1]])
  colnames(sum) <- c("estimate", "post.mean", "l.95.CI", "u.95.CI", "eff.samp", "pMCMC")
  VCVsum <- data.frame(cbind("climvar.value", as.character(diff$climatevar[1]), as.character(diff$Taxon.group[1]), mean((fitmod$VCV[,4]+fitmod$Sol[,2]^2)/(rowSums(fitmod$VCV[,c(4,6)])+fitmod$Sol[,2]^2)), HPDinterval((fitmod$VCV[,4]+fitmod$Sol[,2]^2)/(rowSums(fitmod$VCV[,c(4,6)])+fitmod$Sol[,2]^2))[1], HPDinterval((fitmod$VCV[,4]+fitmod$Sol[,2]^2)/(rowSums(fitmod$VCV[,c(4,6)])+fitmod$Sol[,2]^2))[2]))
  colnames(VCVsum) <- c("estimate", "climatevar", "Taxon.group", "meanVCV", "VCVlower", "VCVupper")
  merge(sum, VCVsum)
}

# All: Run MCMCglmm models - this takes a while! ----

est.grad.all <- grad %>% group_by(climatevar) %>% do(fit.grad(as.data.frame(.)))
est.grad.all$meanVCV <- as.numeric(as.character(est.grad.all$meanVCV))
est.grad.all$VCVlower <- as.numeric(as.character(est.grad.all$VCVlower))
est.grad.all$VCVupper <- as.numeric(as.character(est.grad.all$VCVupper))

est.diff.all <- diff %>% group_by(climatevar) %>% do(fit.diff(as.data.frame(.)))
est.diff.all$meanVCV <- as.numeric(as.character(est.diff.all$meanVCV))
est.diff.all$VCVlower <- as.numeric(as.character(est.diff.all$VCVlower))
est.diff.all$VCVupper <- as.numeric(as.character(est.diff.all$VCVupper))

# All taxa: Run MCMCglmm models - this takes a while! ----

est.grad <- grad %>% group_by(Taxon.group, climatevar) %>% do(fit.grad(as.data.frame(.)))
est.grad$meanVCV <- as.numeric(as.character(est.grad$meanVCV))
est.grad$VCVlower <- as.numeric(as.character(est.grad$VCVlower))
est.grad$VCVupper <- as.numeric(as.character(est.grad$VCVupper))

est.diff <- diff %>% group_by(Taxon.group, climatevar) %>% do(fit.diff(as.data.frame(.)))
est.diff$meanVCV <- as.numeric(as.character(est.diff$meanVCV))
est.diff$VCVlower <- as.numeric(as.character(est.diff$VCVlower))
est.diff$VCVupper <- as.numeric(as.character(est.diff$VCVupper))

# All: Prepare data for plotting (Figure 2) ----
temp.grad <- est.grad %>% filter(climatevar == "MinTempMinStd" | climatevar == "MaxTempMaxStd" | climatevar == "AvgTempMeanStd" | climatevar == "SDTempMeanStd")
temp.grad$climatevar <- factor(temp.grad$climatevar, levels(temp.grad$climatevar) <- c("MinTempMinStd", "MaxTempMaxStd", "AvgTempMeanStd", "SDTempMeanStd"))

precip.grad <- est.grad %>% filter(climatevar == "MinPrecipStd" | climatevar == "MaxPrecipStd" | climatevar == "AvgPrecipStd" | climatevar == "SDPrecipStd")
precip.grad$climatevar <- factor(precip.grad$climatevar, levels(precip.grad$climatevar) <- c("MinPrecipStd", "MaxPrecipStd", "AvgPrecipStd", "SDPrecipStd"))

temp.diff <- est.diff %>% filter(climatevar == "MinTempMinStd" | climatevar == "MaxTempMaxStd" | climatevar == "AvgTempMeanStd" | climatevar == "SDTempMeanStd")
temp.diff$climatevar <- factor(temp.diff$climatevar, levels(temp.diff$climatevar) <- c("MinTempMinStd", "MaxTempMaxStd", "AvgTempMeanStd", "SDTempMeanStd"))

precip.diff <- est.diff %>% filter(climatevar == "MinPrecipStd" | climatevar == "MaxPrecipStd" | climatevar == "AvgPrecipStd" | climatevar == "SDPrecipStd")
precip.diff$climatevar <- factor(precip.diff$climatevar, levels(precip.diff$climatevar) <- c("MinPrecipStd", "MaxPrecipStd", "AvgPrecipStd", "SDPrecipStd"))

# Standardise PET ----
grad.PET <- dat %>% filter(complete.cases(AvgPET), complete.cases(MaxPET), complete.cases(MinPET), complete.cases(sdPET), complete.cases(Grad.linear.value), complete.cases(Grad.linear.StErr)) %>% filter(!Study.ID %in% dat.PET.novariance$Study.ID) 
STD.grad <- function(x){unlist(tapply(x, grad.PET$study.trait, function(x){if(length(x)>1 & var(x)>0){scale(x)}else{rep(0, length(x))}}))}
grad.PET <- grad.PET %>% mutate(AvgPETStd = STD.grad(AvgPET), MaxPETStd = STD.grad(MaxPET), MinPETStd = STD.grad(MinPET), sdPETStd = STD.grad(sdPET))

diff.PET <- dat %>% filter(complete.cases(AvgPET), complete.cases(MaxPET), complete.cases(MinPET), complete.cases(sdPET), complete.cases(Diff.linear.value), complete.cases(Diff.linear.StErr)) %>% filter(!Study.ID %in% dat.PET.novariance$Study.ID) 
STD.diff <- function(x){unlist(tapply(x, diff.PET$study.trait, function(x){if(length(x)>1 & var(x)>0){scale(x)}else{rep(0, length(x))}}))}
diff.PET <- diff.PET %>% mutate(AvgPETStd = STD.diff(AvgPET), MaxPETStd = STD.diff(MaxPET), MinPETStd = STD.diff(MinPET), sdPETStd = STD.diff(sdPET))

# Reshape to long format
grad.PET1 <- gather(grad.PET, climatevar, climvar.value, AvgPETStd, MaxPETStd, MinPETStd, sdPETStd)
diff.PET1 <- gather(diff.PET, climatevar, climvar.value, AvgPETStd, MaxPETStd, MinPETStd, sdPETStd)

# Function to fit PET models ----

fit.grad <- function(grad.PET1) {
  fitmod <- MCMCglmm(Grad.linear.value~climvar.value, random=~us(1+climvar.value):study.trait, mev=grad.PET1$Grad.linear.StErr^2, data=grad.PET1, prior=prior, family = "gaussian", pr=TRUE, nitt = 30000, burnin = 2000)
  sum <- data.frame(rbind("Intercept", "climvar.value"), summary(fitmod)[5][[1]])
  colnames(sum) <- c("estimate", "post.mean", "l.95.CI", "u.95.CI", "eff.samp", "pMCMC")
  VCVsum <- data.frame(cbind("climvar.value", as.character(grad.PET1$climatevar[1]), as.character(grad.PET1$Taxon.group[1]), mean((fitmod$VCV[,4]+fitmod$Sol[,2]^2)/(rowSums(fitmod$VCV[,c(4,6)])+fitmod$Sol[,2]^2)), HPDinterval((fitmod$VCV[,4]+fitmod$Sol[,2]^2)/(rowSums(fitmod$VCV[,c(4,6)])+fitmod$Sol[,2]^2))[1], HPDinterval((fitmod$VCV[,4]+fitmod$Sol[,2]^2)/(rowSums(fitmod$VCV[,c(4,6)])+fitmod$Sol[,2]^2))[2]))
  colnames(VCVsum) <- c("estimate", "climatevar", "Taxon.group", "meanVCV", "VCVlower", "VCVupper")
  merge(sum, VCVsum)
}

fit.diff <- function(diff.PET1) {
  fitmod <- MCMCglmm(Diff.linear.value~climvar.value, random=~us(1+climvar.value):study.trait, mev=diff.PET1$Diff.linear.StErr^2, data=diff.PET1, prior=prior, family = "gaussian", pr=TRUE, nitt = 30000, burnin = 2000)
  sum <- data.frame(rbind("Intercept", "climvar.value"), summary(fitmod)[5][[1]])
  colnames(sum) <- c("estimate", "post.mean", "l.95.CI", "u.95.CI", "eff.samp", "pMCMC")
  VCVsum <- data.frame(cbind("climvar.value", as.character(diff.PET1$climatevar[1]), as.character(diff.PET1$Taxon.group[1]), mean((fitmod$VCV[,4]+fitmod$Sol[,2]^2)/(rowSums(fitmod$VCV[,c(4,6)])+fitmod$Sol[,2]^2)), HPDinterval((fitmod$VCV[,4]+fitmod$Sol[,2]^2)/(rowSums(fitmod$VCV[,c(4,6)])+fitmod$Sol[,2]^2))[1], HPDinterval((fitmod$VCV[,4]+fitmod$Sol[,2]^2)/(rowSums(fitmod$VCV[,c(4,6)])+fitmod$Sol[,2]^2))[2]))
  colnames(VCVsum) <- c("estimate", "climatevar", "Taxon.group", "meanVCV", "VCVlower", "VCVupper")
  merge(sum, VCVsum)
}

# PET: Run MCMCglmm models - this takes a while! ----
est.grad.PET.all <- grad.PET1 %>%  group_by(climatevar) %>% do(fit.grad(as.data.frame(.)))
est.grad.PET.all$meanVCV <- as.numeric(as.character(est.grad.PET.all$meanVCV))
est.grad.PET.all$VCVlower <- as.numeric(as.character(est.grad.PET.all$VCVlower))
est.grad.PET.all$VCVupper <- as.numeric(as.character(est.grad.PET.all$VCVupper))

est.diff.PET.all <- diff.PET1 %>%  group_by(climatevar) %>% do(fit.diff(as.data.frame(.)))
est.diff.PET.all$meanVCV <- as.numeric(as.character(est.diff.PET.all$meanVCV))
est.diff.PET.all$VCVlower <- as.numeric(as.character(est.diff.PET.all$VCVlower))
est.diff.PET.all$VCVupper <- as.numeric(as.character(est.diff.PET.all$VCVupper))

# PET taxa: Run MCMCglmm models - this takes a while! ----
est.grad.PET <- grad.PET1 %>%  group_by(Taxon.group, climatevar) %>% do(fit.grad(as.data.frame(.)))
est.grad.PET$meanVCV <- as.numeric(as.character(est.grad.PET$meanVCV))
est.grad.PET$VCVlower <- as.numeric(as.character(est.grad.PET$VCVlower))
est.grad.PET$VCVupper <- as.numeric(as.character(est.grad.PET$VCVupper))

est.diff.PET <- diff.PET1 %>%  group_by(Taxon.group, climatevar) %>% do(fit.diff(as.data.frame(.)))
est.diff.PET$meanVCV <- as.numeric(as.character(est.diff.PET$meanVCV))
est.diff.PET$VCVlower <- as.numeric(as.character(est.diff.PET$VCVlower))
est.diff.PET$VCVupper <- as.numeric(as.character(est.diff.PET$VCVupper))

# PET: Prepare data for plotting (Figure 2) ----

PET.grad <- est.grad.PET %>% filter(climatevar == "AvgPETStd" | climatevar == "MaxPETStd" | climatevar == "MinPETStd" | climatevar == "sdPETStd")
PET.grad$climatevar <- factor(PET.grad$climatevar, levels(PET.grad$climatevar) <- c("MinPETStd", "MaxPETStd", "AvgPETStd", "sdPETStd"))

PET.diff <- est.diff.PET %>% filter(climatevar == "AvgPETStd" | climatevar == "MaxPETStd" | climatevar == "MinPETStd" | climatevar == "sdPETStd")
PET.diff$climatevar <- factor(PET.diff$climatevar, levels(PET.diff$climatevar) <- c("MinPETStd", "MaxPETStd", "AvgPETStd", "sdPETStd"))

# All: Prepare data for plotting ----
temp.grad.all <- est.grad.all %>% filter(climatevar == "MinTempMinStd" | climatevar == "MaxTempMaxStd" | climatevar == "AvgTempMeanStd" | climatevar == "SDTempMeanStd")
temp.grad.all$climatevar <- factor(temp.grad.all$climatevar, levels(temp.grad.all$climatevar) <- c("MinTempMinStd", "MaxTempMaxStd", "AvgTempMeanStd", "SDTempMeanStd"))

precip.grad.all <- est.grad.all %>% filter(climatevar == "MinPrecipStd" | climatevar == "MaxPrecipStd" | climatevar == "AvgPrecipStd" | climatevar == "SDPrecipStd")
precip.grad.all$climatevar <- factor(precip.grad.all$climatevar, levels(precip.grad.all$climatevar) <- c("MinPrecipStd", "MaxPrecipStd", "AvgPrecipStd", "SDPrecipStd"))

temp.diff.all <- est.diff.all %>% filter(climatevar == "MinTempMinStd" | climatevar == "MaxTempMaxStd" | climatevar == "AvgTempMeanStd" | climatevar == "SDTempMeanStd")
temp.diff.all$climatevar <- factor(temp.diff.all$climatevar, levels(temp.diff.all$climatevar) <- c("MinTempMinStd", "MaxTempMaxStd", "AvgTempMeanStd", "SDTempMeanStd"))

precip.diff.all <- est.diff.all %>% filter(climatevar == "MinPrecipStd" | climatevar == "MaxPrecipStd" | climatevar == "AvgPrecipStd" | climatevar == "SDPrecipStd")
precip.diff.all$climatevar <- factor(precip.diff.all$climatevar, levels(precip.diff.all$climatevar) <- c("MinPrecipStd", "MaxPrecipStd", "AvgPrecipStd", "SDPrecipStd"))

PET.grad.all <- est.grad.PET.all %>% filter(climatevar == "AvgPETStd" | climatevar == "MaxPETStd" | climatevar == "MinPETStd" | climatevar == "sdPETStd")
PET.grad.all$climatevar <- factor(PET.grad.all$climatevar, levels(PET.grad.all$climatevar) <- c("MinPETStd", "MaxPETStd", "AvgPETStd", "sdPETStd"))

PET.diff.all <- est.diff.PET.all %>% filter(climatevar == "AvgPETStd" | climatevar == "MaxPETStd" | climatevar == "MinPETStd" | climatevar == "sdPETStd")
PET.diff.all$climatevar <- factor(PET.diff.all$climatevar, levels(PET.diff.all$climatevar) <- c("MinPETStd", "MaxPETStd", "AvgPETStd", "sdPETStd"))

# All: Plot model slopes ----

# Temperature Gradients
temp.grad.fig <- ggplot(temp.grad.all, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "darkgrey") +
  geom_point(size = 2, colour = "darkgrey") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.2, label = "A. Grad. Temp.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Temperature Differentials
temp.diff.fig <- ggplot(temp.diff.all, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "darkgrey") +
  geom_point(size = 2, colour = "darkgrey") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.2, label = "B. Diff. Temp.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Precipitation Gradients
precip.grad.fig <- ggplot(precip.grad.all, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "darkgrey") +
  geom_point(size = 2, colour = "darkgrey") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.2, label = "C. Grad. Precip.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Precipitation Differentials
precip.diff.fig <- ggplot(precip.diff.all, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "darkgrey") +
  geom_point(size = 2, colour = "darkgrey") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.2, label = "D. Diff. Precip.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# PET Gradients
PET.grad.fig <- ggplot(PET.grad.all, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "darkgrey") +
  geom_point(size = 2, colour = "darkgrey") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.2, label = "E. Diff. PET", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# PET Differentials
PET.diff.fig <- ggplot(PET.diff.all, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "darkgrey") +
  geom_point(size = 2, colour = "darkgrey") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.2, label = "F. Diff. PET", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# All: Plot model VCVs ----

# Temperature Gradients
temp.grad.fig.all.vcv <- ggplot(temp.grad.all, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-0.75, 0.75)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "darkgrey") +
  geom_point(size = 2, colour = "darkgrey") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.75, label = "A. Inverts: Grad. Temp.", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Temperature Differentials
temp.diff.fig.all.vcv <- ggplot(temp.diff.all, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-0.75, 0.75)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "darkgrey") +
  geom_point(size = 2, colour = "darkgrey") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.75, label = "B. Inverts: Diff. Temp.", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Precipitation Gradients
precip.grad.fig.all.vcv <- ggplot(precip.grad.all, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-0.75, 0.75)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "darkgrey") +
  geom_point(size = 2, colour = "darkgrey") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.75, label = "C. Inverts: Grad. Precip.", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Preciptation Differentials
precip.diff.fig.all.vcv <- ggplot(precip.diff.all, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-0.75, 0.75)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "darkgrey") +
  geom_point(size = 2, colour = "darkgrey") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.75, label = "D. Inverts: Diff. Precip.", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# PET Gradients
PET.grad.fig.all.vcv <- ggplot(PET.grad.all, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-0.75, 0.75)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "darkgrey") +
  geom_point(size = 2, colour = "darkgrey") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.75, label = "E. Inverts: Grad. PET", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# PET Differentials
PET.diff.fig.all.vcv <- ggplot(PET.diff.all, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-0.75, 0.75)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "darkgrey") +
  geom_point(size = 2, colour = "darkgrey") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.75, label = "F. Inverts: Diff. PET", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

pdf(file="figures/Siepielski_VCV_all.pdf", width = 5.5, height = 7)

grid.arrange(temp.grad.fig.all.vcv, temp.diff.fig.all.vcv, precip.grad.fig.all.vcv, precip.diff.fig.all.vcv, PET.grad.fig.all.vcv, PET.diff.fig.all.vcv, ncol = 2, left = "Variance - Covariance Est.", bottom = "Climate Variables")

dev.off()

# Taxa: Prepare data for plotting ----
temp.grad.invert <- est.grad %>% subset(Taxon.group == "I") %>% filter(climatevar == "MinTempMinStd" | climatevar == "MaxTempMaxStd" | climatevar == "AvgTempMeanStd" | climatevar == "SDTempMeanStd")
temp.grad.invert$climatevar <- factor(temp.grad.invert$climatevar, levels(temp.grad.invert$climatevar) <- c("MinTempMinStd", "MaxTempMaxStd", "AvgTempMeanStd", "SDTempMeanStd"))

precip.grad.invert <- est.grad %>% subset(Taxon.group == "I") %>% filter(climatevar == "MinPrecipStd" | climatevar == "MaxPrecipStd" | climatevar == "AvgPrecipStd" | climatevar == "SDPrecipStd")
precip.grad.invert$climatevar <- factor(precip.grad.invert$climatevar, levels(precip.grad.invert$climatevar) <- c("MinPrecipStd", "MaxPrecipStd", "AvgPrecipStd", "SDPrecipStd"))

temp.diff.invert <- est.diff %>% subset(Taxon.group == "I") %>% filter(climatevar == "MinTempMinStd" | climatevar == "MaxTempMaxStd" | climatevar == "AvgTempMeanStd" | climatevar == "SDTempMeanStd")
temp.diff.invert$climatevar <- factor(temp.diff.invert$climatevar, levels(temp.diff.invert$climatevar) <- c("MinTempMinStd", "MaxTempMaxStd", "AvgTempMeanStd", "SDTempMeanStd"))

precip.diff.invert <- est.diff %>% subset(Taxon.group == "I") %>% filter(climatevar == "MinPrecipStd" | climatevar == "MaxPrecipStd" | climatevar == "AvgPrecipStd" | climatevar == "SDPrecipStd")
precip.diff.invert$climatevar <- factor(precip.diff.invert$climatevar, levels(precip.diff.invert$climatevar) <- c("MinPrecipStd", "MaxPrecipStd", "AvgPrecipStd", "SDPrecipStd"))

PET.grad.invert <- PET.grad %>% subset(Taxon.group == "I")

PET.diff.invert <- PET.diff %>% subset(Taxon.group == "I")

temp.grad.plant <- temp.grad %>% subset(Taxon.group == "P") %>% filter(climatevar == "MinTempMinStd" | climatevar == "MaxTempMaxStd" | climatevar == "AvgTempMeanStd" | climatevar == "SDTempMeanStd")
temp.grad.plant$climatevar <- factor(temp.grad.plant$climatevar, levels(temp.grad.plant$climatevar) <- c("MinTempMinStd", "MaxTempMaxStd", "AvgTempMeanStd", "SDTempMeanStd"))

precip.grad.plant <- precip.grad %>% subset(Taxon.group == "P") %>% filter(climatevar == "MinPrecipStd" | climatevar == "MaxPrecipStd" | climatevar == "AvgPrecipStd" | climatevar == "SDPrecipStd")
precip.grad.plant$climatevar <- factor(precip.grad.plant$climatevar, levels(precip.grad.plant$climatevar) <- c("MinPrecipStd", "MaxPrecipStd", "AvgPrecipStd", "SDPrecipStd"))

temp.diff.plant <- temp.diff %>% subset(Taxon.group == "P") %>% filter(climatevar == "MinTempMinStd" | climatevar == "MaxTempMaxStd" | climatevar == "AvgTempMeanStd" | climatevar == "SDTempMeanStd")
temp.diff.plant$climatevar <- factor(temp.diff.plant$climatevar, levels(temp.diff.plant$climatevar) <- c("MinTempMinStd", "MaxTempMaxStd", "AvgTempMeanStd", "SDTempMeanStd"))

precip.diff.plant <- precip.diff %>% subset(Taxon.group == "P") %>% filter(climatevar == "MinPrecipStd" | climatevar == "MaxPrecipStd" | climatevar == "AvgPrecipStd" | climatevar == "SDPrecipStd")
precip.diff.plant$climatevar <- factor(precip.diff.plant$climatevar, levels(precip.diff.plant$climatevar) <- c("MinPrecipStd", "MaxPrecipStd", "AvgPrecipStd", "SDPrecipStd"))

PET.grad.plant <- PET.grad %>% subset(Taxon.group == "P")
PET.diff.plant <- PET.diff %>% subset(Taxon.group == "P")

temp.grad.vert <- temp.grad %>% subset(Taxon.group == "V") %>% filter(climatevar == "MinTempMinStd" | climatevar == "MaxTempMaxStd" | climatevar == "AvgTempMeanStd" | climatevar == "SDTempMeanStd")
temp.grad.vert$climatevar <- factor(temp.grad.vert$climatevar, levels(temp.grad.vert$climatevar) <- c("MinTempMinStd", "MaxTempMaxStd", "AvgTempMeanStd", "SDTempMeanStd"))

precip.grad.vert <- precip.grad %>% subset(Taxon.group == "V") %>% filter(climatevar == "MinPrecipStd" | climatevar == "MaxPrecipStd" | climatevar == "AvgPrecipStd" | climatevar == "SDPrecipStd")
precip.grad.vert$climatevar <- factor(precip.grad.vert$climatevar, levels(precip.grad.vert$climatevar) <- c("MinPrecipStd", "MaxPrecipStd", "AvgPrecipStd", "SDPrecipStd"))

temp.diff.vert <- temp.diff %>% subset(Taxon.group == "V") %>% filter(climatevar == "MinTempMinStd" | climatevar == "MaxTempMaxStd" | climatevar == "AvgTempMeanStd" | climatevar == "SDTempMeanStd")
temp.diff.vert$climatevar <- factor(temp.diff.vert$climatevar, levels(temp.diff.vert$climatevar) <- c("MinTempMinStd", "MaxTempMaxStd", "AvgTempMeanStd", "SDTempMeanStd"))

precip.diff.vert <- precip.diff %>% subset(Taxon.group == "V") %>% filter(climatevar == "MinPrecipStd" | climatevar == "MaxPrecipStd" | climatevar == "AvgPrecipStd" | climatevar == "SDPrecipStd")
precip.diff.vert$climatevar <- factor(precip.diff.vert$climatevar, levels(precip.diff.vert$climatevar) <- c("MinPrecipStd", "MaxPrecipStd", "AvgPrecipStd", "SDPrecipStd"))

PET.grad.vert <- PET.grad %>% subset(Taxon.group == "V")
PET.diff.vert <- PET.diff %>% subset(Taxon.group == "V")

# Inverts: Plot model slopes ----

# Temperature Gradients
temp.grad.fig.invert <- ggplot(temp.grad.invert, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "orange") +
  geom_point(size = 2, colour = "orange") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.2, label = "A. Grad. Temp.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Temperature Differentials
temp.diff.fig.invert <- ggplot(temp.diff.invert, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "orange") +
  geom_point(size = 2, colour = "orange") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.2, label = "B. Diff. Temp.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Precipitation Gradients
precip.grad.fig.invert <- ggplot(precip.grad.invert, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "orange") +
  geom_point(size = 2, colour = "orange") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.2, label = "C. Grad. Precip.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Precipitation Differentials
precip.diff.fig.invert <- ggplot(precip.diff.invert, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "orange") +
  geom_point(size = 2, colour = "orange") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.2, label = "D. Diff. Precip.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# PET Gradients
PET.grad.fig.invert <- ggplot(PET.grad.invert, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "orange") +
  geom_point(size = 2, colour = "orange") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.2, label = "E. Diff. PET", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# PET Differentials
PET.diff.fig.invert <- ggplot(PET.diff.invert, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "orange") +
  geom_point(size = 2, colour = "orange") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.2, label = "F. Diff. PET", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Inverts: Plot model VCVs ----

# Temperature Gradients
temp.grad.fig.invert.vcv <- ggplot(temp.grad.invert, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "orange") +
  geom_point(size = 2, colour = "orange") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "A. Inverts: Grad. Temp.", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Temperature Differentials
temp.diff.fig.invert.vcv <- ggplot(temp.diff.invert, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "orange") +
  geom_point(size = 2, colour = "orange") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "B. Inverts: Diff. Temp.", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Precipitation Gradients
precip.grad.fig.invert.vcv <- ggplot(precip.grad.invert, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "orange") +
  geom_point(size = 2, colour = "orange") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "C. Inverts: Grad. Precip.", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Preciptation Differentials
precip.diff.fig.invert.vcv <- ggplot(precip.diff.invert, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "orange") +
  geom_point(size = 2, colour = "orange") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "D. Inverts: Diff. Precip.", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# PET Gradients
PET.grad.fig.invert.vcv <- ggplot(PET.grad.invert, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "orange") +
  geom_point(size = 2, colour = "orange") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "E. Inverts: Grad. PET", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# PET Differentials
PET.diff.fig.invert.vcv <- ggplot(PET.diff.invert, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "orange") +
  geom_point(size = 2, colour = "orange") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "F. Inverts: Diff. PET", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Plants: Plot model slopes ----

# Temperature Gradients
temp.grad.fig.plant <- ggplot(temp.grad.plant, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.08, 0.08)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "green4") +
  geom_point(size = 2, colour = "green4") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.08, label = "G. Grad. Temp.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Temperature Differentials
temp.diff.fig.plant <- ggplot(temp.diff.plant, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.08, 0.08)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "green4") +
  geom_point(size = 2, colour = "green4") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.08, label = "H. Diff. Temp.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Precipitation Gradients
precip.grad.fig.plant <- ggplot(precip.grad.plant, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.08, 0.08)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "green4") +
  geom_point(size = 2, colour = "green4") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.08, label = "I. Grad. Precip.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Precipitation Differentials
precip.diff.fig.plant <- ggplot(precip.diff.plant, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.08, 0.08)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "green4") +
  geom_point(size = 2, colour = "green4") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.08, label = "J. Diff. Precip.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# PET Gradients
PET.grad.fig.plant <- ggplot(PET.grad.plant, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.08, 0.08)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "green4") +
  geom_point(size = 2, colour = "green4") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.08, label = "K. Diff. PET", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# PET Differentials
PET.diff.fig.plant <- ggplot(PET.diff.plant, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.08, 0.08)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "green4") +
  geom_point(size = 2, colour = "green4") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.08, label = "L. Diff. PET", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Plants: Plot model VCVs ----

# Temperature Gradients
temp.grad.fig.plant.vcv <- ggplot(temp.grad.plant, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "green4") +
  geom_point(size = 2, colour = "green4") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "G. Plants: Grad. Temp.", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Temperature Differentials
temp.diff.fig.plant.vcv <- ggplot(temp.diff.plant, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "green4") +
  geom_point(size = 2, colour = "green4") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "H. Plants: Diff. Temp.", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Precipitation Gradients
precip.grad.fig.plant.vcv <- ggplot(precip.grad.plant, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "green4") +
  geom_point(size = 2, colour = "green4") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "I. Plants: Grad. Precip.", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Preciptation Differentials
precip.diff.fig.plant.vcv <- ggplot(precip.diff.plant, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "green4") +
  geom_point(size = 2, colour = "green4") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "J. Plants Diff. Precip.", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# PET Gradients
PET.grad.fig.plant.vcv <- ggplot(PET.grad.plant, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "green4") +
  geom_point(size = 2, colour = "green4") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "K. Plants: Grad. PET", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# PET Differentials
PET.diff.fig.plant.vcv <- ggplot(PET.diff.plant, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "green4") +
  geom_point(size = 2, colour = "green4") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "L. Plants: Diff. PET", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Verts: Plot model slopes ----

# Temperature Gradients
temp.grad.fig.vert <- ggplot(temp.grad.vert, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.08, 0.08)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "blue3") +
  geom_point(size = 2, colour = "blue3") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.08, label = "M. Grad. Temp.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Temperature Differentials
temp.diff.fig.vert <- ggplot(temp.diff.vert, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.08, 0.08)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "blue3") +
  geom_point(size = 2, colour = "blue3") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.08, label = "N. Diff. Temp.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Precipitation Gradients
precip.grad.fig.vert <- ggplot(precip.grad.vert, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.08, 0.08)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "blue3") +
  geom_point(size = 2, colour = "blue3") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.08, label = "O. Grad. Precip.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Precipitation Differentials
precip.diff.fig.vert <- ggplot(precip.diff.vert, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.08, 0.08)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "blue3") +
  geom_point(size = 2, colour = "blue3") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.08, label = "P. Diff. Precip.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# PET Differentials
PET.grad.fig.vert <- ggplot(PET.grad.vert, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.08, 0.08)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "blue3") +
  geom_point(size = 2, colour = "blue3") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.08, label = "Q. Diff. Temp.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# PET Differentials
PET.diff.fig.vert <- ggplot(PET.diff.vert, aes(x = climatevar, y = post.mean)) +
  scale_y_continuous(limits = c(-0.08, 0.08)) +
  geom_errorbar(aes(ymin=(l.95.CI), ymax=(u.95.CI)), width=0.2, colour = "blue3") +
  geom_point(size = 2, colour = "blue3") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 0.08, label = "R. Diff. Precip.", size=4, hjust=0) +
  labs(x = "", y = "Slopes") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Verts: Plot model VCVs ----

# Temperature Gradients
temp.grad.fig.vert.vcv <- ggplot(temp.grad.vert, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "blue3") +
  geom_point(size = 2, colour = "blue3") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "M. Verts: Grad. Temp.", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Temperature Differentials
temp.diff.fig.vert.vcv <- ggplot(temp.diff.vert, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "blue3") +
  geom_point(size = 2, colour = "blue3") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "N. Verts: Diff. Temp.", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Precipitation Gradients
precip.grad.fig.vert.vcv <- ggplot(precip.grad.vert, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "blue3") +
  geom_point(size = 2, colour = "blue3") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "O. Verts: Grad. Precip.", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Preciptation Differentials
precip.diff.fig.vert.vcv <- ggplot(precip.diff.vert, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "blue3") +
  geom_point(size = 2, colour = "blue3") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "P. Verts: Diff. Precip.", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# PET Gradients
PET.grad.fig.vert.vcv <- ggplot(PET.grad.vert, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "blue3") +
  geom_point(size = 2, colour = "blue3") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "Q. Verts: Grad. PET", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# PET Differentials
PET.diff.fig.vert.vcv <- ggplot(PET.diff.vert, aes(x = climatevar, y = meanVCV)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
  geom_errorbar(aes(ymin=(VCVlower), ymax=(VCVupper)), width=0.2, colour = "blue3") +
  geom_point(size = 2, colour = "blue3") +
  geom_hline(yintercept = 0) +
  annotate("text", x = 0.5, y = 1.2, label = "R. Verts: Diff. PET", size=4, hjust=0) +
  labs(x = "", y = "") + 
  theme_tidy() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels = c("min.", "max.", "mean", "SD")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# All taxa (Table S5) ----

pdf(file="figures/Siepielski_VCV_taxa.pdf", width = 15, height = 7)

grid.arrange(temp.grad.fig.invert.vcv, temp.diff.fig.invert.vcv, precip.grad.fig.invert.vcv, precip.diff.fig.invert.vcv, PET.grad.fig.invert.vcv, PET.diff.fig.invert.vcv, temp.grad.fig.plant.vcv, temp.diff.fig.plant.vcv, precip.grad.fig.plant.vcv, precip.diff.fig.plant.vcv, PET.grad.fig.plant.vcv, PET.diff.fig.plant.vcv, temp.grad.fig.vert.vcv, temp.diff.fig.vert.vcv, precip.grad.fig.vert.vcv, precip.grad.fig.vert.vcv, PET.grad.fig.vert.vcv, PET.diff.fig.vert.vcv, ncol = 6, left = "Variance - Covariance Est.", bottom = "Climate Variables")

dev.off()

# Selection estimates vs max. precipitation ----

sel.invert <- ggplot(subset(dat, Taxon.group == "I")) +
  geom_point(aes(x = MaxPrecip, y = Grad.linear.value, size = Grad.linear.StErr), alpha = 0.4, colour = "orange", size = 1) +
  geom_point(aes(x = MaxPrecip, y = Diff.linear.value, size = Diff.linear.StErr), alpha = 0.4, colour = "orange", size = 1) +
  scale_y_continuous(limits = c(-7, 7)) +
  scale_x_continuous(limits = c(0, 600)) +
  annotate("text", x = 5, y = 6.5, label = "A. Invertebrates", size=4, hjust=0) +
  labs(x = "\n", y = "\n") + 
  theme_tidy() +
  theme(legend.position="none") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

sel.plant <- ggplot(subset(dat, Taxon.group == "P")) +
  geom_point(aes(x = MaxPrecip, y = Grad.linear.value, size = Grad.linear.StErr), alpha = 0.4, colour = "green4", size = 1) +
  geom_point(aes(x = MaxPrecip, y = Diff.linear.value, size = Diff.linear.StErr), alpha = 0.4, colour = "green4", size = 1) +
  scale_y_continuous(limits = c(-7, 7)) +
  scale_x_continuous(limits = c(0, 600)) +
  annotate("text", x = 5, y = 6.5, label = "B. Plants", size=4, hjust=0) +
  labs(x = "\n", y = "Selection Est.\n") + 
  theme_tidy() +
  theme(legend.position="none") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

sel.vert <- ggplot(subset(dat, Taxon.group == "V")) +
  geom_point(aes(x = MaxPrecip, y = Grad.linear.value, size = Grad.linear.StErr), alpha = 0.4, colour = "blue3", size = 1) +
  geom_point(aes(x = MaxPrecip, y = Diff.linear.value, size = Diff.linear.StErr), alpha = 0.4, colour = "blue3", size = 1) +
  scale_y_continuous(limits = c(-7, 7)) +
  scale_x_continuous(limits = c(0, 600)) +
  annotate("text", x = 5, y = 6.5, label = "C. Vertebrates", size=4, hjust=0) +
  labs(x = "\nMax. Precip. (mm)", y = "\n") + 
  theme_tidy() +
  theme(legend.position="none") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

grid.arrange(sel.invert, sel.plant, sel.vert, ncol = 1, nrow = 3)

# Fitness measures figure ----

sel.class.invert <- ggplot(subset(dat, Taxon.group == "I")) +
  geom_bar(aes(Generic.Trait.Class), alpha = 0.8, fill = "orange") +
  scale_x_discrete(limits=c("BE", "MO", "O", "OLH", "OMO", "PC", "Phenology", "Size")) +
  annotate("text", x = 0.75, y = 950, label = "A. invertebrates", size=4, hjust=0) +
  theme_tidy() +
  theme(legend.position="none") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

sel.class.plant <- ggplot(subset(dat, Taxon.group == "P")) +
  geom_bar(aes(Generic.Trait.Class), alpha = 0.8, fill = "green4") +
  scale_x_discrete(limits=c("BE", "MO", "O", "OLH", "OMO", "PC", "Phenology", "Size")) +
  annotate("text", x = 0.75, y = 1050, label = "B. Plants", size=4, hjust=0) +
  theme_tidy() +
  theme(legend.position="none") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

sel.class.vert <- ggplot(subset(dat, Taxon.group == "V")) +
  geom_bar(aes(Generic.Trait.Class), alpha = 0.8, fill = "blue3") +
  scale_x_discrete(limits=c("BE", "MO", "O", "OLH", "OMO", "PC", "Phenology", "Size")) +
  annotate("text", x = 0.75, y = 1000, label = "C. Vertebrates", size=4, hjust=0) +
  theme_tidy() +
  theme(legend.position="none") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

grid.arrange(sel.class.invert, sel.class.plant, sel.class.vert)

# Map figure ----

figdata <- dat %>% select(Final.Index, Study.ID, Population, Taxon.group, Trait.Class, Long.Decimal.Degrees, Lat.Decimal.Degrees, Grad.linear.value, Grad.linear.StErr, Diff.linear.value, Diff.linear.StErr, MinTempMin, MaxTempMax,AvgTempMean, SDTempMean, MaxPrecip, MinPrecip, AvgPrecip, SDPrecip, MaxPET, MinPET, AvgPET, sdPET) %>% mutate(error = ifelse(is.na(ifelse(is.na(Grad.linear.StErr), as.numeric(Diff.linear.StErr), as.numeric(Grad.linear.StErr))), 0, ifelse(is.na(Grad.linear.StErr), as.numeric(Diff.linear.StErr), as.numeric(Grad.linear.StErr)))) %>% subset(error != "0")

map_world <- borders(database = "world", colour = "gray50", fill = "gray80", lwd = 0.2)  # We used the `Colour Picker` Addin to pick the colours

sel.map <- ggplot() + map_world +  # Plot the map
  geom_point(data = figdata,  # Specify the data for geom_point()
             aes(x = Long.Decimal.Degrees,  # Specify the x axis as longitude
                 y = Lat.Decimal.Degrees,  # Specify the y axis as latitude
                 colour = Taxon.group,  # Colour the points based on species name
                 size = error),  
             alpha = 0.4,  # Set point opacity to 40%
             pch = 19) +  # Set point size to 1
  scale_color_manual(values=c(I = "orange", P = "green4", V = "blue3")) +   # Specify the colour palette to colour the points
  theme_tidy() +  # Remove gridlines and shading inside the plot
  ylab(expression("Latitude ("*degree*")" )) +  # Add a smarter x axis label
  xlab(expression("Longitude ("*degree*")" )) +  # Add a smarter y axis label
  annotate("text", x = -200, y = 80, label = "A.", size=5, hjust=0)+
  theme(legend.position="none") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

temp.fig <- ggplot(figdata, aes(x = Lat.Decimal.Degrees, y = SDTempMean, colour = Taxon.group, size = error)) +
  scale_y_continuous(limits = c(0, 20)) +
  geom_point(alpha = 0.4) +
  annotate("text", x = -40, y = 20, label = "B. Lat. - Temp. Var.", size=5, hjust=0) +
  labs(x = "\nLatitude\n", y = "SD Temp.\n") +   
  scale_color_manual(values=c(I = "orange", P = "green4", V = "blue3")) +   # Specify the colour palette to colour the points
  theme_tidy() +
  theme(legend.position="none") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

temp.fig2 <- ggplot(figdata, aes(x = AvgTempMean, y = SDTempMean, colour = Taxon.group, size = error)) +
  scale_y_continuous(limits = c(0, 30)) +
  geom_point(alpha = 0.4) +
  annotate("text", x = -20, y = 30, label = "C. Temp. Mean - Var.", size=5, hjust=0) +
  labs(x = "\nMean Temp.\n", y = "SD Temp.\n") + 
  scale_color_manual(values=c(I = "orange", P = "green4", V = "blue3")) +   # Specify the colour palette to colour the points
  theme_tidy() +
  theme(legend.position="none") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

precip.fig <- ggplot(figdata, aes(x = Lat.Decimal.Degrees, y = SDPrecip, colour = Taxon.group, size = error)) +
  scale_y_continuous(limits = c(0, 300)) +
  geom_point(alpha = 0.4) +
  annotate("text", x = -40, y = 300, label = "D. Lat. - Precip. Var.", size=5, hjust=0) +
  labs(x = "\nLatitude\n", y = "SD Precip.\n") + 
  scale_color_manual(values=c(I = "orange", P = "green4", V = "blue3")) +   # Specify the colour palette to colour the points
  theme_tidy() +
  theme(legend.position="none") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

precip.fig2 <- ggplot(figdata, aes(x = AvgPrecip, y = SDPrecip, colour = Taxon.group, size = error)) +
  scale_y_continuous(limits = c(0, 300)) +
  geom_point(alpha = 0.4) +
  annotate("text", x = 10, y = 300, label = "E. Precip. Mean - Var.", size=5, hjust=0) +
  labs(x = "\nMean Precip.\n", y = "SD Precip.\n") + 
  scale_color_manual(values=c(I = "orange", P = "green4", V = "blue3")) +   # Specify the colour palette to colour the points
  theme_tidy() +
  theme(legend.position="none") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

PET.fig <- ggplot(figdata, aes(x = Lat.Decimal.Degrees, y = sdPET, colour = Taxon.group, size = error)) +
  scale_y_continuous(limits = c(0, 3)) +
  geom_point(alpha = 0.4) +
  annotate("text", x = -40, y = 3, label = "F. Lat. - PET Var.", size=5, hjust=0) +
  labs(x = "\nLatitude\n", y = "SD PET\n") + 
  scale_color_manual(values=c(I = "orange", P = "green4", V = "blue3")) +   # Specify the colour palette to colour the points
  theme_tidy() +
  theme(legend.position="none") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

PET.fig2 <- ggplot(figdata, aes(x = AvgPET, y = sdPET, colour = Taxon.group, size = error)) +
  scale_y_continuous(limits = c(0, 3)) +
  geom_point(alpha = 0.4) +
  annotate("text", x = 0.5, y = 3, label = "G. PET Mean - Var.", size=5, hjust=0) +
  labs(x = "\nMean PET\n", y = "SD PET\n") + 
  scale_color_manual(values=c(I = "orange", P = "green4", V = "blue3")) +   # Specify the colour palette to colour the points
  theme_tidy() +
  theme(legend.position="none") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

pdf(file="figures/Siepielski_map.pdf", width = 12, height = 12)

grid.arrange(sel.map, temp.fig, precip.fig, PET.fig, temp.fig2, precip.fig2, PET.fig2, ncol = 3, nrow = 3, layout_matrix = rbind(c(1,1,1), c(2,3,4), c(5,6,7)), heights=c(5, 3, 3))

dev.off()

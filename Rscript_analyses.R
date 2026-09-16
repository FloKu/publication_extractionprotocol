# --------------------------------------------------------------
# Analyses of the publication: Which extraction protocol should I use for my historic fish samples? A performance analysis on widely used methods
# --------------------------------------------------------------
# version: 16.09.2026
# r code authored by Florian Kunz
# script corresponds to publication.

# How to cite: cite corresponding paper.


{library("dplyr")
  library("tidyr")
  library("ggplot2")
  library("cowplot")  
  library("lme4")     
  library("lmerTest") 
  library("pbkrtest")
  library("MASS")
  library("car")
  library("sjPlot")
  library("AICcmodavg")
  library("performance")
  library("DHARMa")
  library("multcompView")}

data <- read.csv2("/data_long.csv")


# 1) Data manipulation
#*******************************
str(data)

# DNA_conc into numeric
{data$DNA_conc 
  data$DNA_conc[data$DNA_conc == "too low"] <- 0
  data$DNA_conc <- gsub(",", ".", data$DNA_conc)
  data$DNA_conc <- as.numeric(data$DNA_conc)
  # DNA_conc_2 into be numeric
  data$DNA_conc_2 
  data$DNA_conc_2[data$DNA_conc_2 == "too low"] <- 0
  data$DNA_conc_2 <- gsub(",", ".", data$DNA_conc_2)
  data$DNA_conc_2 <- as.numeric(data$DNA_conc_2)
  # AFL into be numeric
  data$AFL
  data$AFL[data$AFL == "46 + 1908?"] <- NA
  data$AFL[data$AFL == "136,5"] <- 136
  data$AFL[data$AFL == "2470?"] <- NA
  data$AFL[data$AFL == "473 + 1986?"] <- NA
  data$AFL[data$AFL == "84 + 1294? + 2107?"] <- NA
  data$AFL <- as.numeric(data$AFL)
  # Year into numeric
  data$year
  data$year[data$year == "b 1825"] <- 1825
  data$year <- as.numeric(data$year)}

# managing columns
{data <- data[-c(5, 6, 9, 10, 11, 19, 21, 22, 24)]
  names(data)[names(data) == "species_2"] <- "species"
  names(data)[names(data) == "ind."] <- "ind"}

# rename values and convert into factor variable
{levels(as.factor(data$extraction_method))
  data$extraction_method[data$extraction_method =="Dabney/Rohland protocol"] <- "Dabney"
  data$extraction_method[data$extraction_method =="GENial (24h Verdau)"] <- "GENial"
  data$extraction_method[data$extraction_method =="MagMax: Valeria protocol (Auswahl!)"] <- "MagMax"
  data$extraction_method[data$extraction_method =="Patzoldi: MonarchNAdestructive sampling of nonNAtype specimens"] <- "Patzoldi"
  data$extraction_method[data$extraction_method =="PCI Extraction"] <- "PCI"
  data$extraction_method[data$extraction_method =="Promega Case work kit (Maxwell)"] <- "Promega"
  data$extraction_method[data$extraction_method =="Qiagen Dneasy (old protocol) (24h)"] <- "Qiagen_1"
  data$extraction_method[data$extraction_method =="Qiagen Dneasy inc. Freezing step (24h)"] <- "Qiagen_2"
  data$extraction_method[data$extraction_method =="Qiagen Dneasy inc. Freezing step (24h) + EtOH removal buffer (Auswahl!)"] <- "Qiagen_3"
  data$extraction_method <- factor(data$extraction_method, levels = c("Qiagen_1","Qiagen_2", "PCI", "GENial", "Promega", "Dabney", "Qiagen_3", "MagMax", "Patzoldi"))}


# 2) Exploratory
#*******************************
str(data)

# barplot of years
subset(data, extraction_method == "Qiagen_1") %>%
  ggplot(aes(x = year, y=after_stat(count), group=year)) +
    geom_bar(fill="darkblue") +
    #geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.5) +
    geom_text(aes(label = year), vjust = .2, hjust=-.1, stat="count", angle = 90) +
    scale_y_continuous("number of samples", minor_breaks=c, breaks=c(0, 1, 2, 5, 10), limits=c(0, 13)) +
    labs(title = "Age of aDNA samples", x="year of collection") +
    theme(legend.position = "none")

# histogram of DNA conc
subset(data, DNA_conc != 0) %>% 
  ggplot(aes(x = DNA_conc)) +
    geom_histogram(binwidth = 1)

subset(data, DNA_conc != 0) %>% 
  summary()    

# barplot showing number of samples per extraction method
data %>% 
  group_by(extraction_method) %>% 
  dplyr::summarize(count = sum(!is.na(DNA_conc)), na = sum(is.na(DNA_conc))) %>%
  mutate(sum = count + na) %>%
  print() %>%
  ggplot(aes(x = extraction_method, y = count)) +
    geom_bar(stat="identity", fill="darkblue") +
    labs(title = "Number of samples per extraction method", x="extraction method", y="number of samples", fill="Samples") +
    theme(legend.position = "none")

# scatterplot over all methods for total extracted DNA
data %>% 
  mutate(DNA_total = DNA_conc*elution_vol) %>% 
  ggplot(aes(x=year, y=DNA_total, color = extraction_method)) + 
  geom_point()+
  geom_smooth(method=lm, se=F) +
  scale_x_continuous(limits=c(1820, 1910)) +
  #scale_y_continuous(limits=c(0, 200)) + 
  labs(title="Extraction success", x="sampling year", y="total extracted DNA (ng)")


# scatterplot over all methods for DNA concentration
data %>% 
  mutate(DNA_total = DNA_conc*elution_vol) %>% 
  ggplot(aes(x=year, y=DNA_per_mg, color = extraction_method)) + 
  geom_point()+
  geom_smooth(method=lm, se=F) +
  scale_x_continuous(limits=c(1820, 1910)) +
  #scale_y_continuous(limits=c(0, 200)) + 
  labs(title="Extraction success", x="sampling year", y="DNA concentration (ng/ul)")

# scatterplot per method for total extracted DNA
for (ex in c(levels(as.factor(data$extraction_method)))) {
  print(paste0("plotting: ", ex))
  data_subset <- data %>% 
    subset(extraction_method==ex) %>% 
    mutate(DNA_total = DNA_conc*elution_vol)
  print(ggplot(data_subset, aes(x=year, y=DNA_total)) +
        geom_point(color="darkblue")+
        geom_smooth(method=lm) +
        scale_x_continuous(limits=c(1820, 1910)) +
        labs(title=paste0("Extraction success: ", ex), x="sampling year", y="total extracted DNA (ng)"))
  Sys.sleep(5)
  rm(data_subset)
}

# flipped grouped boxplot on full data 
filter(data, is.na(DNA_conc)==FALSE) %>% # no effect on plot, but it manages the warning message (removing NAs)
  ggplot(aes(x=extraction_method, y=DNA_conc, fill=tissue_type)) +
    geom_boxplot() +
    scale_x_discrete(limits = rev(levels(data$extraction_method))) +
    coord_flip() +
    labs(title="All samples", x="extraction protocol", y="DNA Concentration (ng/ul)")

# flipped boxplot on DNA_conc on full data, no groups 
png("/DNAconc_all.png", 8000, 4000, res=1000, type='cairo-png', antialias=c("subpixel")) # to safe the file with anti-aliasing, tweak res and scale always together; I have no idea what anialias=c() does
filter(data, is.na(DNA_conc)==FALSE) %>% # no effect on plot, but it manages the warning message (removing NAs)
  ggplot(aes(x=extraction_method, y=DNA_conc)) +
  geom_boxplot() +
  scale_x_discrete(limits = rev(levels(data$extraction_method))) +
  coord_flip() +
  labs(title="All samples", x="extraction protocol", y="DNA Concentration (ng/ul)")
dev.off()

# flipped boxplot on DNA_conc on subset (exuding very high DNA concentration)
subset(data, DNA_conc < 100) %>%
  ggplot(aes(x=extraction_method, y=DNA_conc, fill=tissue_type)) +
    geom_boxplot() +
    scale_x_discrete(limits = rev(levels(data$extraction_method))) +
    coord_flip() +
    #facet_wrap(~tissue_type) +
    labs(title="Without high DNA yield", x="extraction protocol", y="DNA Concentration (ng/ul)")

# flipped grouped boxplot for AFL on full data 
filter(data, is.na(AFL)==FALSE) %>% # no effect on plot, but it manages the warning message (removing NAs)
  ggplot(aes(x=extraction_method, y=AFL, fill=tissue_type)) +
    geom_boxplot() +
    scale_x_discrete(limits = rev(levels(data$extraction_method))) +
    coord_flip() +
    labs(title="All samples", x="extraction protocol", y="average fragment length (bp)")

png("/AFL_all.png", 8000, 4000, res=1000, type='cairo-png', antialias=c("subpixel")) # to safe the file with anti-aliasing, tweak res and scale always together; I have no idea what anialias=c() does
# flipped boxplot for AFL on full data, no groups 
filter(data, is.na(AFL)==FALSE) %>% # no effect on plot, but it manages the warning message (removing NAs)
  ggplot(aes(x=extraction_method, y=AFL)) +
    geom_boxplot() +
    scale_x_discrete(limits = rev(levels(data$extraction_method))) +
    coord_flip() +
    labs(title="All samples", x="extraction protocol", y="average fragment length (bp)")
dev.off()

# flipped grouped boxplot for DNA_per_mg on full data 
filter(data, is.na(DNA_per_mg)==FALSE) %>% # no effect on plot, but it manages the warning message (removing NAs)
  ggplot(aes(x=extraction_method, y=DNA_per_mg, fill=tissue_type)) +
  geom_boxplot() +
  scale_x_discrete(limits = rev(levels(data$extraction_method))) +
  coord_flip() +
  labs(title="All samples", x="extraction protocol", y="DNA (ng) per mg tissue used")

# flipped grouped boxplot for DNA_per_mg on subset, excluding three outliers in Quiagen_2 
filter(data, is.na(DNA_per_mg)==FALSE) %>% # no effect on plot, but it manages the warning message (removing NAs)
  subset(DNA_per_mg < 150) %>% 
  ggplot(aes(x=extraction_method, y=DNA_per_mg)) +
  geom_boxplot() +
  scale_x_discrete(limits = rev(levels(data$extraction_method))) +
  coord_flip() +
  labs(title="All samples", x="extraction protocol", y="DNA (ng) per mg tissue used")

# flipped boxplot for DNA_per_mg on subset, excluding three outliers in Quiagen_2, no groups
png("/DNAquan_all.png", 8000, 4000, res=1000, type='cairo-png', antialias=c("subpixel")) # to safe the file with anti-aliasing, tweak res and scale always together; I have no idea what anialias=c() does
filter(data, is.na(DNA_per_mg)==FALSE) %>% # no effect on plot, but it manages the warning message (removing NAs)
  subset(DNA_per_mg < 150) %>% 
  ggplot(aes(x=extraction_method, y=DNA_per_mg)) +
  geom_boxplot() +
  scale_x_discrete(limits = rev(levels(data$extraction_method))) +
  coord_flip() +
  labs(title="All samples, three outliers for Quiagen_2 excluded", x="extraction protocol", y="DNA (ng) per mg tissue used")
dev.off()


# 3) Modelling
#*******************************
# response variable: DNA_conc, DNA_per_mg, AFL
# fixed factors: extraction_method, year, tissue_quan, tissue_type, elution_vol
# random intercept: extraction_method
# no random slope:
#   because there is no slope as our explanatory variable is categorical
#   and we dont have the data for it (repeated measurements)

# DV continuously, EV categorical -> hence ANOVA 

anov <- aov(DNA_conc~extraction_method + year + elution_vol + tissue_type + ind, subset_33)
summary(anov)
par(mfrow=c(2,2))
plot(anov)
par(mfrow=c(1,1))

# FINAL RESULT:
#****************
# IND has high significance! -> so we need LMM.
# Same results also for ANOVAs on DNA_per_mg and AFL! 

names(data)
str(data)

# subsetting
# subset_33 includes only those extraction methods which have been applied on all 33 individuals
subset_33 <- subset(data, extraction_method == "Qiagen_1" | extraction_method == "Qiagen_2" | extraction_method == "PCI" | extraction_method == "GENial" | extraction_method == "Promega") %>% 
  droplevels(.$extraction_method) %>% 
  mutate(ind = as.factor(ind))

# subset_10 inlcudes those 10 individuals that have been subject to all extraction methods
subset_10 <- subset(data, ind %in% c("Smic37", "Smic39", "Smic40", "Smic41", "Pale9", "Pale17", "Pale20", "Sleb1", "Imel1", "Aleu1"))

# 3.1) Preliminary checks: correlation tissue_quan with extraction_method
#*******************************
#----
# maybe there is a systematic correlation (might have used more tissue in specific methods, maybe methods also demand more tissue)

hist(subset_33$tissue_quan, breaks = 20)

str(subset_33)

# plot
subset_33 %>% 
  ggplot(aes(x=extraction_method, y=tissue_quan)) +
    geom_boxplot() +
    scale_x_discrete(limits = rev(levels(subset_33$extraction_method))) +
    coord_flip() +
    labs(title="All samples", x="extraction protocol", y="tissue quantity (mg)")

# parametric test
anova <- aov(tissue_quan ~ extraction_method, subset_33)
summary(anova) # significant differences
par(mfrow=c(2,2))
plot(anova)
par(mfrow=c(1,1))

# non-parametric test
kruskal.test(tissue_quan ~ extraction_method, subset_33) # significant differences

# info: anova probably not ok (no homogeneity); kruskal significant, but does this imply correlation?
#   -> check true multi-colinearity in VIF in the model validation itself
#----

# 3.2) Preliminary checks: correlation elution_volume with extraction_method
#*******************************
#----
# maybe there is a systematic correlation (because protocols demand higher elections or because lab assisstant might has done it purposely)
# elution_volumne is not really spread out - its mostly one value

hist(subset_33$elution_vol, breaks = 20)

str(subset_33)

subset_33 <- subset(subset_33, extraction_method != "Promega") %>% 
  droplevels(.$extraction_method) 

# plot
subset_33 %>% 
  ggplot(aes(x=extraction_method, y=elution_vol)) +
  geom_boxplot() +
  scale_x_discrete(limits = rev(levels(subset_33$extraction_method))) +
  coord_flip() +
  labs(title="All samples", x="extraction protocol", y="elution volumne (ul)")
# even after removing the problematic method (Promega), anova still significant

# parametric test
anova <- aov(elution_vol ~ extraction_method, subset_33)
summary(anova) # significant differences
par(mfrow=c(2,2))
plot(anova)
par(mfrow=c(1,1))

# non-parametric test
kruskal.test(elution_vol ~ extraction_method, data) # significant differences

# info: anova probably not ok (no homogeneity); kruskal significant, but does this imply correlation?
#   -> check true multi-colinearity in VIF in the model validation itself
#----

# 3.3) Modelling process: BASIC LMM with full validation pipeline
#*******************************
#----
model1 <- lmerTest::lmer(DNA_conc ~ extraction_method + year + tissue_quan + tissue_type + (1|ind), subset_33, REML=T)
model1 <- lmerTest::lmer(DNA_conc ~ extraction_method + (1|ind), subset_33, REML=T)

# 3.3.1) model assumptions

dharma <- simulateResiduals(model1, plot=F)

# A) linearity or relationship between Y and x
plot(DNA_conc~year, subset_33_m1)
plot(DNA_conc~extraction_method, subset_33_m1)
# Result: OK for year, not necessary for categorical variable

# B) No multi-colinearity
car::vif(model1)
cor(x=subset_33_m1$year, y=subset_33_m1$tissue_quan, method ="pearson")
# Result: OK, no correlation between continuous variables

# C) homogeneity of residuals for the entire model
# Histogramm of residuals
DHARMa::testDispersion(dharma)
hist(resid(model1))
plot_model(model1, type="diag")[[3]]
# Result: looks good

# QQplot
DHARMa::plotQQunif(dharma) # this plot checks the distribution of residuals (normality), deviations can be interpret like for linear regression
DHARMa::testUniformity(dharma)
{qqnorm(resid(model1))
  qqline(resid(model1))}
plot_model(model1, type="diag")[[1]]
# Result: not good, there is NO uniformity of residuals... so I guess no homogeneity

# shapiro test
shapiro.test(resid(model1))
# Result: no homogeneity by shapiro (although shapiro might not be ideal for LMMs)

# D) Homoscedasticity
# full model: residuals vs fitted plot
plot(residuals~fitted, subset_33_m1)
plot(model1)
plot_model(model1, type="diag")[[4]]

DHARMa::plotResiduals(dharma) # this plot checks the homogeneity of residuals , deviations can be interpret like for linear regression

# per explanatory variable
DHARMa::testCategorical(dharma, catPred = subset_33$extraction_method)
DHARMa::testCategorical(dharma, catPred = subset_33$tissue_type)
DHARMa::plotResiduals(dharma, subset_33$year)

# For factors: The distribution for each factor level should be uniformly distributed, 
# so the box should go from 0.25 to 0.75, with the median line at 0.5 (within-group ). 
# To test if deviations from those expectations are significant, KS-tests per group and 
# a Levene test for homogeneity of variances is performed. See testCategorical for details.

# For continuous data: the plot function calculates an (optional) quantile regression of the residuals, 
# by default for the 0.25, 0.5 and 0.75 quantiles. As the residuals should be uniformly distributed for a correctly specified model, 
# the theoretical expectations for these regressions are straight lines at 0.25, 0.5 and 0.75, 
# which are displayed as dashed black lines on the plot. Some deviations from these expectations are to be expected by chance, 
# however, even for a perfect model, especially if the sample size is small. 
# The function therefore tests if deviation of the fitted quantile regression from the expectation is significant, using testQuantiles.

# D) Tests for validity per random factor group/independent variable

# Actually, extraction_method is the fixed factor and ind is the random factor.
# But testing per ind for extraction_method would result in 33 plots, where each plot only has 5 points.
# Therefore, I am looping the diagnostics over levels of etrxaction_method - I believe this is still useful, as it shouldnt matter
# whether individuals or extraction_method are the grouping variable?
# Individuals per extraction_method or extraction_methods per individual - results in the same groups?

# Resids vs: Fitted:    linear relationship
# QQ:                   homogeneity of residuals
# Scale-Loc: 
# Cook's Dist:          checking for outliers
# Resids vs. Leverage: 

plot_n = 1 # 1: Resids vs Fitted, 2: QQ, 3: Scale-Location, 4: Resids vs. Leverage

{par(mfrow=c(3,2))
  for (method in levels(subset_33_1$extraction_method)) {
    s <- subset(subset_33_1, extraction_method==method)
    if (is.null(levels(s$tissue_type)) == TRUE) {m <- lm(AFL~year + elution_vol, s)} else {
      m <- lm(AFL~year + tissue_type + elution_vol, s) }    
    plot(m, plot_n)
  }
  par(mfrow=c(1,1))
}

# E) Zero-Inflation
# should only be checked if there is a reason for checking, see manpage
# also usually done on count data (discrete variable), not on continuous data
testZeroInflation(dharma)

# F) Not sure what to learn from this
boxplot(residuals~extraction_method, data=subset_33_m1)
plot_model(model1, type="diag")[[2]]$ind

# 3.3.2) statistical significance
anova(model1)
anova(model1, ddf= "Kenward-Roger")

# 3.3.3) model results
print(model1)
summary(model1)
summary(model1, ddf= "Kenward-Roger")
fixef(model1) # extract point estimates of fixed effects
ranef(model1) # extracting the point estimates of the random effects (different intercepts of individuals)
confint(model1) # confidence intervals for the fixed effects


# Result of BASIC LMM:
# Homogeneity (normality) is NOT given (c)! Checking plots per extraction method (D) reveals that data seems right tailed!
# Also, Homoscedasticity is not given (E) and checking data per extraction method reveals that for Qiagen_2 and PCI, variance is increasing, 
# while for GENial and Promega

# Results of solutions: 
#~~~~~~~~~~~~~~~~~~~~~~
# box-cox transformation:                 this solves heterogeneity and heteroscedasticity and looks good! Only downfall: resi.fitt plot for variable "year"
# log(y+1) transformation:                this solves heterogeneity, but not heteroscedasticity
# inverse hyperbolic sine transformation:  dito
# other transformations (sqrt, inverse, ...): didnt solve anything
# models on data without 0 in y:          didnt solve it
# models on data without the "most problematic" extraction_method "Qiagen_2": didnt solve it
# models on standardized continuous data: didnt solve heteroscedasticity

# subset_33 without 0es in y
subset_33_no0 <- subset_33 %>% 
  filter(!DNA_conc == 0)

# subset_33 without extraction_method "Qiagen_2" 
subset_33_noQ <- subset_33 %>% 
  filter(!extraction_method == "Qiagen_2") %>% 
  {.$extraction_method <- factor(.$extraction_method, levels = c("Qiagen_1","PCI","GENial","Promega")); .}

# subset with standardizes variables
subset_33 <- subset_33 %>% 
  mutate(tissue_quan_s = as.vector(scale(tissue_quan))) %>% 
  mutate(year_s = as.vector(scale(year)))

# box-cox
min(filter(subset_33, DNA_conc > 0)$DNA_conc) # smallest value: 0.059 - as boxcox needs a constant added to all values - I made the constant very small
model <- lm(DNA_conc+0.0001~extraction_method + year + tissue_type, subset_33)
bc <- boxcox(model)
lambda <- bc$x[which.max(bc$y)]
lambda

# box-cox transformed
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(((DNA_conc^lambda-1)/lambda) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer(((DNA_conc^lambda-1)/lambda) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
# log transformed and +1
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(log(DNA_conc+1) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
# inverse hyperbolic sine
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(log(DNA_conc+((DNA_conc^2+1)^0.5)) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer(log(DNA_conc+((DNA_conc^2+1)^0.5)) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
# no DNA_conc==0 AND log transformed
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(log(DNA_conc) ~ extraction_method + year + tissue_type + (1|ind), subset_33_no0, REML=T), plot=F))
# log transformed
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(log(DNA_conc+0.001) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
# inverse transformed
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(1/(DNA_conc+0.001) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
# quadratic
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(DNA_conc^2 ~ extraction_method + year + tissue_type + (1|ind), subset_33_no0, REML=T), plot=F))

# Result: use boxcox
#----
 

# 4) Modelling DNA_conc
#*******************************
# DNA_conc: concentration of DNA (ng) per ul elution
# Dependent variables: year, extraction_method, DNA_quan, elution_vol, tissue_type
# Elution_vol hypothesis: the more volumne used to elute, the lesser the concentration

# 4.1) Modelling process: Multi-colinearity by VIF
#*******************************
#----
model <- lm(DNA_conc+0.001~extraction_method + tissue_type + tissue_quan + year + elution_vol, subset_33)
bc <- boxcox(model)
lambda <- bc$x[which.max(bc$y)]

model <- lmerTest::lmer((((DNA_conc+0.001)^lambda-1)/lambda) ~ extraction_method + tissue_type + tissue_quan + year + elution_vol + (1|ind), subset_33, REML=T)

car::vif(model)

# model with only tissue_quan and extraction_method
model <- lm(DNA_conc+0.001~extraction_method + tissue_type + tissue_quan + year + elution_vol, subset_33)
bc <- boxcox(model)
lambda <- bc$x[which.max(bc$y)]

model <- lmerTest::lmer((((DNA_conc+0.00)^lambda-1)/lambda) ~ extraction_method + tissue_type + tissue_quan + year + elution_vol + (1|ind), subset_33, REML=T)

car::vif(model)

# model with only eluction_vol and extraction_method
model <- lm(DNA_conc+0.001~extraction_method + tissue_type + year + elution_vol, subset_33)
bc <- boxcox(model)
lambda <- bc$x[which.max(bc$y)]

model <- lmerTest::lmer((((DNA_conc+0.001)^lambda-1)/lambda) ~ extraction_method + tissue_type + year + elution_vol + (1|ind), subset_33, REML=T)

car::vif(model)

# Results
# extraction_method is correlated with tissue_quan and elution_vol and should be discarded
# As this is our main variable, I have to discard tissue_quan
# A model with elution_vol is still ok, although in that case, VIF are almost too high (9.3)
#----

# 4.2) FINAL Model 1: BOXCOX FULL
#*******************************
# FULL means including all possible variables
#----
# Lambda: 0.1818182
# Lambda is rather close to 0, where BOxCox would actually be a log transformation. 
# Hence our transformation is quite related to log transformation (which on its own also showed good results)

model <- lm(DNA_conc+0.00001~extraction_method + year + elution_vol + tissue_type, subset_33)
bc <- boxcox(model)
lambda <- bc$x[which.max(bc$y)]

model <- lmerTest::lmer(((DNA_conc^lambda-1)/lambda) ~ extraction_method + year + elution_vol + tissue_type + (1|ind), subset_33, REML=T)

dharma <- simulateResiduals(model, plot=T)
subset_33_m <- subset_33 %>%
  mutate(residuals = resid(model)) %>% 
  mutate(fitted = fitted(model))

# Multi-colinearity
car::vif(model)

# Homogeneity
DHARMa::plotQQunif(dharma) # this plot checks the distribution of residuals (normality), deviations can be interpret like for linear regression
DHARMa::testUniformity(dharma)
{par(mfrow=c(3,2))
  for (thing in levels(subset_33_m$extraction_method)) {
    data_1 <- subset(subset_33_m, extraction_method==thing)
    qqnorm(data_1$residuals, main=thing)
    qqline(data_1$residuals)
  }
  par(mfrow=c(1,1))
}

# Homoscedasticity
# full model
DHARMa::plotResiduals(dharma) # this plot checks the homogeneity of residuals , deviations can be interpret like for linear regression

# per explanatory variable
DHARMa::testCategorical(dharma, catPred = subset_33$extraction_method)
DHARMa::testCategorical(dharma, catPred = subset_33$tissue_type)
DHARMa::plotResiduals(dharma, subset_33$year)

{par(mfrow=c(3,2))
  for (thing in levels(subset_33_m$extraction_method)) {
    plot(residuals~fitted, subset(subset_33_m, extraction_method==thing), main=thing)
  }
  par(mfrow=c(1,1))
}

# results
anova(model)
anova(model, ddf= "Kenward-Roger")
Anova(model)

AIC(model)

print(model)
summary(model)
summary(model, ddf= "Kenward-Roger")
fixef(model) # extract point estimates of fixed effects
ranef(model) # extracting the point estimates of the random effects (different intercepts of individuals)
confint(model) # confidence intervals for the fixed effects

boxcox_full <- model
#----

# 4.3) FINAL Model 3: BOXCOX MEDIUM
#*******************************
# MEDIUM: without tissue_type
#----
# Lambda: 0.1818182
# Lambda is rather close to 0, where BOxCox would actually be a log transformation. 
# Hence our transformation is quite related to log transformation (which on its own also showed good results)

model <- lm(DNA_conc+0.0001~extraction_method + elution_vol + year, subset_33)
bc <- boxcox(model)
lambda <- bc$x[which.max(bc$y)]

model <- lmerTest::lmer(((DNA_conc^lambda-1)/lambda) ~ extraction_method + elution_vol + year + (1|ind), subset_33, REML=T)

dharma <- simulateResiduals(model, plot=T)
subset_33_m <- subset_33 %>%
  mutate(residuals = resid(model)) %>% 
  mutate(fitted = fitted(model))

# Multi-colinearity
car::vif(model)

# Homogeneity
DHARMa::plotQQunif(dharma) # this plot checks the distribution of residuals (normality), deviations can be interpret like for linear regression
DHARMa::testUniformity(dharma)
{par(mfrow=c(3,2))
  for (thing in levels(subset_33_m$extraction_method)) {
    data_1 <- subset(subset_33_m, extraction_method==thing)
    qqnorm(data_1$residuals, main=thing)
    qqline(data_1$residuals)
  }
  par(mfrow=c(1,1))
}

# Homoscedasticity
# full model
DHARMa::plotResiduals(dharma) # this plot checks the homogeneity of residuals , deviations can be interpret like for linear regression

# per explanatory variable
DHARMa::testCategorical(dharma, catPred = subset_33$extraction_method)
DHARMa::testCategorical(dharma, catPred = subset_33$tissue_type)
DHARMa::plotResiduals(dharma, subset_33$year)

{par(mfrow=c(3,2))
  for (thing in levels(subset_33_m$extraction_method)) {
    plot(residuals~fitted, subset(subset_33_m, extraction_method==thing), main=thing)
  }
  par(mfrow=c(1,1))
}

# results
anova(model)
anova(model, ddf= "Kenward-Roger")

AIC(model)

print(model)
summary(model)
summary(model, ddf= "Kenward-Roger")
fixef(model) # extract point estimates of fixed effects
ranef(model) # extracting the point estimates of the random effects (different intercepts of individuals)
confint(model) # confidence intervals for the fixed effects

boxcox_medium <- model
#----

# 4.4) FINAL Model 3: BOXCOX SMALL
#*******************************
# SMALL: no variable except the significant ones
#----
# Lambda: 0.1818182
# Lambda is rather close to 0, where BOxCox would actually be a log transformation. 
# Hence our transformation is quite related to log transformation (which on its own also showed good results)

model <- lm(DNA_conc+0.0001~extraction_method + year, subset_33)
bc <- boxcox(model)
lambda <- bc$x[which.max(bc$y)]

model <- lmerTest::lmer(((DNA_conc^lambda-1)/lambda) ~ extraction_method + year + (1|ind), subset_33, REML=T)

dharma <- simulateResiduals(model, plot=T)
subset_33_m <- subset_33 %>%
  mutate(residuals = resid(model)) %>% 
  mutate(fitted = fitted(model))

# Multi-colinearity
car::vif(model)

# Homogeneity
DHARMa::plotQQunif(dharma) # this plot checks the distribution of residuals (normality), deviations can be interpret like for linear regression
DHARMa::testUniformity(dharma)
{par(mfrow=c(3,2))
  for (thing in levels(subset_33_m$extraction_method)) {
    data_1 <- subset(subset_33_m, extraction_method==thing)
    qqnorm(data_1$residuals, main=thing)
    qqline(data_1$residuals)
  }
  par(mfrow=c(1,1))
}

# Homoscedasticity
# full model
DHARMa::plotResiduals(dharma) # this plot checks the homogeneity of residuals , deviations can be interpret like for linear regression

# per explanatory variable
DHARMa::testCategorical(dharma, catPred = subset_33$extraction_method)
DHARMa::testCategorical(dharma, catPred = subset_33$tissue_type)
DHARMa::plotResiduals(dharma, subset_33$year)

{par(mfrow=c(3,2))
  for (thing in levels(subset_33_m$extraction_method)) {
    plot(residuals~fitted, subset(subset_33_m, extraction_method==thing), main=thing)
  }
  par(mfrow=c(1,1))
}

# results
anova(model)
anova(model, ddf= "Kenward-Roger")

AIC(model)

print(model)
summary(model)
summary(model, ddf= "Kenward-Roger")
fixef(model) # extract point estimates of fixed effects
ranef(model) # extracting the point estimates of the random effects (different intercepts of individuals)
confint(model) # confidence intervals for the fixed effects

boxcox_small <- model
#----


# 5) Compare DNA_conc models
#*******************************
# All models are valid

model_performance(boxcox_full)
model_performance(boxcox_medium)
model_performance(boxcox_small)

anova(boxcox_full)
anova(boxcox_small)

# conditional R2: proportion of variance explained by both the fixed and random factors
# marginal R2: roportion of variance explained by the fixed factor(s) alone

# Note that AIC and AICc values are meaningful to select among gls or lme models fit by maximum
# likelihood. AIC and AICc based on REML are valid to select among different models that only
# differ in their random effects (Pinheiro and Bates 2000). Source: AICcmodavg vignette

models <- list(full = boxcox_full, medium = boxcox_medium, small = boxcox_small)
anova(boxcox_full, boxcox_medium, boxcox_small)
aictab(models)

# FINAL RESULT:
#****************
# All models point to extraction_method and year being significant. 
# Models are pretty similar in their explanatory power, no significant difference, although the smallest one is the best one.


# 6) Modelling DNA_per_mg
#*******************************
# DNA_per_mg: DNA (in ng) per mg tissue used
# Hence, DNA_quan is a useless explanatory variable, and elution_vol also
#----
# not transformed
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(DNA_per_mg ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer(DNA_per_mg ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
# box-cox transformed
min(filter(subset_33, DNA_per_mg > 0)$DNA_per_mg) # smallest value: 0.078 - as boxcox needs a constant added to all values - I made the constant very small
model <- lm(DNA_per_mg+0.0001~extraction_method + tissue_type + year, subset_33)
bc <- boxcox(model)
lambda <- bc$x[which.max(bc$y)]
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(((DNA_per_mg^lambda-1)/lambda) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer(((DNA_per_mg^lambda-1)/lambda) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
# log transformed and +1
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(log(DNA_per_mg+1) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer(log(DNA_per_mg+1) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
# inverse hyperbolic sine
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(log(DNA_per_mg+((DNA_per_mg^2+1)^0.5)) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer(log(DNA_per_mg+((DNA_per_mg^2+1)^0.5)) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
# log transformed
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(log(DNA_per_mg+0.001) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer(log(DNA_per_mg+0.001) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
# inverse transformed
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(1/(DNA_per_mg+0.001) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer(1/(DNA_per_mg+0.001) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T), plot=F))

# Result: hyperbolic spline model seems to be the best, residuals vs. predicted plot looks good, but shows partially significant deviation.... problem?
# I also tried log transformed, but these models resultet in an infinite AICc
#----

# 6.1) Modelling process: Multi-colinearity by VIF
#*******************************
#----
model <- lmerTest::lmer(log(DNA_per_mg+((DNA_per_mg^2+1)^0.5)) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T)
car::vif(model)

# Result: all variables fine as they are
#----

# 6.2.) FINAL Model 4: DNA_per_mg, HPS transformed, FULL
#*******************************
# Note: the homoscedasticity is not perfect.
#   Test shows a slight deviation... but Internet and my assessment says, that's not too much of a problem > model still valid.
#----
model <- lmerTest::lmer(log(DNA_per_mg+((DNA_per_mg^2+1)^0.5)) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T)

dharma <- simulateResiduals(model, plot=T)
subset_33_m <- subset_33 %>%
  mutate(residuals = resid(model)) %>% 
  mutate(fitted = fitted(model))

# Multi-colinearity
car::vif(model)

# Homogeneity
DHARMa::plotQQunif(dharma) # this plot checks the distribution of residuals (normality), deviations can be interpret like for linear regression
DHARMa::testUniformity(dharma)
{par(mfrow=c(3,2))
  for (thing in levels(subset_33_m$extraction_method)) {
    data_1 <- subset(subset_33_m, extraction_method==thing)
    qqnorm(data_1$residuals, main=thing)
    qqline(data_1$residuals)
  }
  par(mfrow=c(1,1))
}

# Homoscedasticity
# full model
DHARMa::plotResiduals(dharma) # this plot checks the homogeneity of residuals , deviations can be interpret like for linear regression

# per explanatory variable
DHARMa::testCategorical(dharma, catPred = subset_33$extraction_method)
DHARMa::testCategorical(dharma, catPred = subset_33$tissue_type)
DHARMa::plotResiduals(dharma, subset_33$year)

{par(mfrow=c(3,2))
  for (thing in levels(subset_33_m$extraction_method)) {
    plot(residuals~fitted, subset(subset_33_m, extraction_method==thing), main=thing)
  }
  par(mfrow=c(1,1))
}

# results
anova(model)
anova(model, ddf= "Kenward-Roger")

AIC(model)

print(model)
summary(model)
summary(model, ddf= "Kenward-Roger")
fixef(model) # extract point estimates of fixed effects
ranef(model) # extracting the point estimates of the random effects (different intercepts of individuals)
confint(model) # confidence intervals for the fixed effects

DNA_mg_FULL <- model
#----

# 6.2.) FINAL Model 4: DNA_per_mg, log-transformed, SMALL
#*******************************
# Note: the homoscedasticity is not perfect.
#   Test shows a slight deviation... but Internet and my assessment says, that's not too much of a problem > model still valid.
#----
model <- lmerTest::lmer(log(DNA_per_mg+((DNA_per_mg^2+1)^0.5)) ~ extraction_method + year + (1|ind), subset_33, REML=T)

dharma <- simulateResiduals(model, plot=T)
subset_33_m <- subset_33 %>%
  mutate(residuals = resid(model)) %>% 
  mutate(fitted = fitted(model))

# Multi-colinearity
car::vif(model)

# Homogeneity
DHARMa::plotQQunif(dharma) # this plot checks the distribution of residuals (normality), deviations can be interpret like for linear regression
DHARMa::testUniformity(dharma)
{par(mfrow=c(3,2))
  for (thing in levels(subset_33_m$extraction_method)) {
    data_1 <- subset(subset_33_m, extraction_method==thing)
    qqnorm(data_1$residuals, main=thing)
    qqline(data_1$residuals)
  }
  par(mfrow=c(1,1))
}

# Homoscedasticity
# full model
DHARMa::plotResiduals(dharma) # this plot checks the homogeneity of residuals , deviations can be interpret like for linear regression

# per explanatory variable
DHARMa::testCategorical(dharma, catPred = subset_33$extraction_method)
DHARMa::testCategorical(dharma, catPred = subset_33$tissue_type)
DHARMa::plotResiduals(dharma, subset_33$year)

{par(mfrow=c(3,2))
  for (thing in levels(subset_33_m$extraction_method)) {
    plot(residuals~fitted, subset(subset_33_m, extraction_method==thing), main=thing)
  }
  par(mfrow=c(1,1))
}

# results
anova(model)
anova(model, ddf= "Kenward-Roger")

AIC(model)

print(model)
summary(model)
summary(model, ddf= "Kenward-Roger")
fixef(model) # extract point estimates of fixed effects
ranef(model) # extracting the point estimates of the random effects (different intercepts of individuals)
confint(model) # confidence intervals for the fixed effects

DNA_mg_SMALL <- model
#----


# 7) Compare DNA_per_mg models
#*******************************
# Note that both models are not perfect as homoscedasticity is not ideally met

model_performance(DNA_mg_SMALL)
model_performance(DNA_mg_FULL)

anova(DNA_mg_FULL)
anova(DNA_mg_SMALL)

# conditional R2: proportion of variance explained by both the fixed and random factors
# marginal R2: roportion of variance explained by the fixed factor(s) alone

# Note that AIC and AICc values are meaningful to select among gls or lme models fit by maximum
# likelihood. AIC and AICc based on REML are valid to select among different models that only
# differ in their random effects (Pinheiro and Bates 2000). Source: AICcmodavg vignette

models <- list(full = boxcox_full, medium = boxcox_medium, small = boxcox_small)
anova(DNA_mg_FULL, DNA_mg_SMALL)
aictab(models)

# FINAL RESULT:
#****************
# All models point to extraction_method and year being significant. 
# Models are pretty similar in their explanatory power, not significantly different, although the smallest one is the best one.


# 8) Modelling AFL
#*******************************
# AFL might depend on the extract_method, year, tissue_type, elution_vol, sample ID might be still a smart random factor
# tissue_quan removed due to high correlation with extraction method
# checking for transformation

# I couldnt find a fitting transformation (always heteroscedasticity), so I tried removing that one outlier in PCI and NA values
{subset_33_1 <- subset_33 %>% 
    subset(AFL < 400) 
  # After checking for outliers PER EXTRACTION METHOD, I also removed two more points (which clearly were shown as outliers) - didnt make any difference though
  subset_33_1 <- subset_33_1[!(subset_33_1$AFL == 47 & subset_33_1$ind == "Cacu1"),]
  subset_33_1 <- subset_33_1[!(subset_33_1$AFL == 66 & subset_33_1$ind == "Pcro1"),]}

hist(subset_33$AFL, breaks = 50)
hist(subset_33_1$AFL, breaks = 50)

#----
# Model validation: 
# I tested all transformations, with different sets of predictors, differen subsets (33, 33_1) and even with non-linear relationships.
# I calulated tests and plots for full models and per level of extraction method for homogeneity, homoscedasticity, outliers and linear relationship.
# But all models showed problems in heteroscedasticity. 
# It seems like the better heteorscedasticity would be managed, the worse homoscedastictiy would get.

# Result subset_33: no transformation is helping against heteroscedasticity, but boxcox is best
# Result subset_33_1: even with subset_33_1, all transformed LMMs are still somewhat unideal
# Result: plots from models with only method and year are a little better, but still not ideal - boxcox and hyperbolic spline are best
# Results: plots without year not really better

# Final: Settleing on inverse hyperbolic spline 

# not transformed
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(AFL ~ extraction_method + year + tissue_type + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer(AFL ~ extraction_method + year + tissue_type + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))

# box-cox transformed
min(filter(subset_33, AFL > 0)$AFL) # smallest value: 40 - as boxcox needs a constant added to all values - I made the constant very small
model <- lm(AFL+1~extraction_method + year + tissue_type + elution_vol, subset_33_1)
bc <- boxcox(model)
lambda <- bc$x[which.max(bc$y)]
lambda <- -0.2
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(((AFL^lambda-1)/lambda) ~ extraction_method + year + tissue_type + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer(((AFL^lambda-1)/lambda) ~ extraction_method + year + tissue_type + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
# squared
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(AFL^2 ~ extraction_method + year + tissue_type + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer(AFL^2 ~ extraction_method + year + tissue_type + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
# square root
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(sqrt(AFL) ~ extraction_method + year + tissue_type + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer(sqrt(AFL) ~ extraction_method + year + tissue_type + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
# log transformed and + 0.001
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(log(AFL+0.001) ~ extraction_method + year + tissue_type  + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer(log(AFL+0.001) ~ extraction_method + year + tissue_type + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
# log transformed and +1
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(log(AFL+1) ~ extraction_method + year + tissue_type  + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer(log(AFL+1) ~ extraction_method + year + tissue_type  + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
# inverse hyperbolic sine
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(log(AFL+((AFL^2+1)^0.5)) ~ extraction_method + year + tissue_type  + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer(log(AFL+((AFL^2+1)^0.5)) ~ extraction_method + year + tissue_type + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
# reciprocal root
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer((-1/sqrt(AFL)) ~ extraction_method + year + tissue_type  + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer((-1/sqrt(AFL)) ~ extraction_method + year + tissue_type + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
# reciprocal
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer((-1/AFL) ~ extraction_method + year + tissue_type  + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer((-1/AFL) ~ extraction_method + year + tissue_type + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
# reciprocal square
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer((-1/AFL^2) ~ extraction_method + year + tissue_type  + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer((-1/AFL^2) ~ extraction_method + year + tissue_type + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
# inverse transformed
DHARMa::plotQQunif(simulateResiduals(lmerTest::lmer(1/(AFL) ~ extraction_method + year + tissue_type + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))
DHARMa::plotResiduals(simulateResiduals(lmerTest::lmer(1/(AFL) ~ extraction_method + year + tissue_type + elution_vol + (1|ind), subset_33_1, REML=T), plot=F))

# Checking the linearity assumption, just to be sure
plot(AFL~year, subset_33_1)
plot(AFL~elution_vol, subset_33_1)
boxplot(AFL~tissue_type, subset_33_1)
boxplot(AFL~extraction_method, subset_33_1)

# saving data to share it
data <- subset_33_1 %>% 
  dplyr::select(AFL, extraction_method, year, elution_vol, tissue_type, ind)

write.csv(data, "/data_AFL_model.csv")
#----

# 8.1) Modelling process: Multi-colinearity by VIF
#*******************************
#----
model <- lmerTest::lmer(log(AFL+((AFL^2+1)^0.5)) ~ extraction_method + year + tissue_type  + elution_vol + (1|ind), subset_33_1, REML=T)

# Multi-colinearity
car::vif(model)

# Result: all variables fine as they are
#----

# 8.2.) FINAL Model 4: AFL, IHS, FULL
#*******************************
# Note: the homoscedasticity is not ideal, one quantile is worrisome
#----
model <- lmerTest::lmer(log(AFL+((AFL^2+1)^0.5)) ~ extraction_method + year + tissue_type  + elution_vol + (1|ind), subset_33_1, REML=T)

dharma <- simulateResiduals(model, plot=T)

# Multi-colinearity
car::vif(model)

# Homogeneity
DHARMa::plotQQunif(dharma) # this plot checks the distribution of residuals (normality), deviations can be interpret like for linear regression
DHARMa::testUniformity(dharma)

# Homoscedasticity
DHARMa::plotResiduals(dharma) # this plot checks the homogeneity of residuals , deviations can be interpret like for linear regression
DHARMa::testCategorical(dharma, catPred = subset_33_1$extraction_method)
DHARMa::testCategorical(dharma, catPred = subset_33_1$tissue_type)
DHARMa::plotResiduals(dharma, subset_33_1$year)

# Plots per random factor/variable level
plot_n = 5 #1:Resids vs Fitted, 2: QQ, 3: Scale-Location, 4: Cook's Distance, 5: Resids vs. Leverage
{par(mfrow=c(3,2))
  for (method in levels(subset_33_1$extraction_method)) {
    s <- subset(subset_33_1, extraction_method==method)
    if (is.null(levels(s$tissue_type)) == TRUE) {m <- lm(AFL~year + elution_vol, s)} else {
      m <- lm(AFL~year + tissue_type + elution_vol, s) }    
    plot(m, plot_n)
  }
  par(mfrow=c(1,1))
}

# results
anova(model)
anova(model, ddf= "Kenward-Roger")

AIC(model)

print(model)
summary(model)
summary(model, ddf= "Kenward-Roger")
fixef(model) # extract point estimates of fixed effects
ranef(model) # extracting the point estimates of the random effects (different intercepts of individuals)
confint(model) # confidence intervals for the fixed effects

DNA_AFL_FULL <- model
#----

# 8.3.) FINAL Model 4: AFL, IHS, MEDIUM
#*******************************
# Note: Heteroscedasticity a problem
#----
model <- lmerTest::lmer(log(AFL+((AFL^2+1)^0.5)) ~ extraction_method + year + (1|ind), subset_33_1, REML=T)

dharma <- simulateResiduals(model, plot=T)

# Multi-colinearity
car::vif(model)

# Homogeneity
DHARMa::plotQQunif(dharma) # this plot checks the distribution of residuals (normality), deviations can be interpret like for linear regression
DHARMa::testUniformity(dharma)

# Homoscedasticity
DHARMa::plotResiduals(dharma) # this plot checks the homogeneity of residuals , deviations can be interpret like for linear regression
DHARMa::testCategorical(dharma, catPred = subset_33_1$extraction_method)
DHARMa::plotResiduals(dharma, subset_33_1$year)

# Plots per random factor/variable level
plot_n = 3 #1:Resids vs Fitted, 2: QQ, 3: Scale-Location, 4: Cook's Distance, 5: Resids vs. Leverage
{par(mfrow=c(3,2))
  for (method in levels(subset_33_1$extraction_method)) {
    s <- subset(subset_33_1, extraction_method==method)
    if (is.null(levels(s$tissue_type)) == TRUE) {m <- lm(AFL~year + elution_vol, s)} else {
      m <- lm(AFL~year + tissue_type + elution_vol, s) }    
    plot(m, plot_n)
  }
  par(mfrow=c(1,1))
}

# results
anova(model)
anova(model, ddf= "Kenward-Roger")

AIC(model)

print(model)
summary(model)
summary(model, ddf= "Kenward-Roger")
fixef(model) # extract point estimates of fixed effects
ranef(model) # extracting the point estimates of the random effects (different intercepts of individuals)
confint(model) # confidence intervals for the fixed effects

DNA_AFL_MEDIUM <- model
#----

# 8.4.) FINAL Model 4: AFL, HBS, SMALL
#*******************************
# Note: invalid
#----
model <- lmerTest::lmer(log(AFL+((AFL^2+1)^0.5)) ~ extraction_method + (1|ind), subset_33_1, REML=T)

dharma <- simulateResiduals(model, plot=T)

# Multi-colinearity
car::vif(model)

# Homogeneity
DHARMa::plotQQunif(dharma) # this plot checks the distribution of residuals (normality), deviations can be interpret like for linear regression
DHARMa::testUniformity(dharma)

# Homoscedasticity
DHARMa::plotResiduals(dharma) # this plot checks the homogeneity of residuals , deviations can be interpret like for linear regression
DHARMa::testCategorical(dharma, catPred = subset_33_1$extraction_method)
DHARMa::plotResiduals(dharma, subset_33_1$year)

# Plots per random factor/variable level
plot_n = 3 #1:Resids vs Fitted, 2: QQ, 3: Scale-Location, 4: Cook's Distance, 5: Resids vs. Leverage
{par(mfrow=c(3,2))
  for (method in levels(subset_33_1$extraction_method)) {
    s <- subset(subset_33_1, extraction_method==method)
    if (is.null(levels(s$tissue_type)) == TRUE) {m <- lm(AFL~year + elution_vol, s)} else {
      m <- lm(AFL~year + tissue_type + elution_vol, s) }    
    plot(m, plot_n)
  }
  par(mfrow=c(1,1))
}

# results
anova(model)
anova(model, ddf= "Kenward-Roger")

AIC(model)

print(model)
summary(model)
summary(model, ddf= "Kenward-Roger")
fixef(model) # extract point estimates of fixed effects
ranef(model) # extracting the point estimates of the random effects (different intercepts of individuals)
confint(model) # confidence intervals for the fixed effects

DNA_mg_SMALL <- model
#----


# 9) Compare DNA_per_mg models
#*******************************
# Note:
# FULL: one quantile is worrisome in terms of heteroscedastictiy
# MEDIUM: two quantiles are worrisome
# SMALL: invalid

model_performance(DNA_AFL_SMALL)
model_performance(DNA_AFL_MEDIUM)
model_performance(DNA_AFL_FULL)

anova(DNA_AFL_FULL)
anova(DNA_AFL_MEDIUM)

# conditional R2: proportion of variance explained by both the fixed and random factors
# marginal R2: roportion of variance explained by the fixed factor(s) alone

# Note that AIC and AICc values are meaningful to select among gls or lme models fit by maximum
# likelihood. AIC and AICc based on REML are valid to select among different models that only
# differ in their random effects (Pinheiro and Bates 2000). Source: AICcmodavg vignette

models <- list(full = boxcox_full, medium = boxcox_medium, small = boxcox_small)
anova(DNA_AFL_FULL, DNA_AFL_MEDIUM)
aictab(models)

# FINAL RESULT:
#****************
# Full model seems to be the most valid, although medium model has lower AICc (but AICc is not good anyways)
# Full model shows extraction method to be explanatory, year only slightly
# no signifcant difference though between FULL and MEDIUM


# 10) Interaction between extraction method and year of sampling
#*******************************
#----
# Seeing both are highly significant, I wonder whether there is any interaction?
# Meaning depending on the year, one methods works better or not

# Load best fitting models
model <- lm(DNA_conc+0.00001~extraction_method + year + elution_vol + tissue_type, subset_33)
bc <- boxcox(model)
lambda <- bc$x[which.max(bc$y)]
dnaconc_full <- lmerTest::lmer(((DNA_conc^lambda-1)/lambda) ~ extraction_method + year + elution_vol + tissue_type + (1|ind), subset_33, REML=T)

dnamg_full <- lmerTest::lmer(log(DNA_per_mg+((DNA_per_mg^2+1)^0.5)) ~ extraction_method + year + tissue_type + (1|ind), subset_33, REML=T)
afl_full <- lmerTest::lmer(log(AFL+((AFL^2+1)^0.5)) ~ extraction_method + year + tissue_type  + elution_vol + (1|ind), subset_33_1, REML=T)
afl_med <- lmerTest::lmer(log(AFL+((AFL^2+1)^0.5)) ~ extraction_method + year + (1|ind), subset_33_1, REML=T)

# models WITH interaction
model <- lm(DNA_conc+0.00001~extraction_method + year + elution_vol + tissue_type + extraction_method*year, subset_33)
bc <- boxcox(model)
lambda <- bc$x[which.max(bc$y)]
dnaconc_full_INT <- lmerTest::lmer(((DNA_conc^lambda-1)/lambda) ~ extraction_method + year + elution_vol + tissue_type +  extraction_method*year + (1|ind), subset_33, REML=T)

dnamg_full_INT <- lmerTest::lmer(log(DNA_per_mg+((DNA_per_mg^2+1)^0.5)) ~ extraction_method + year + tissue_type + extraction_method*year + (1|ind), subset_33, REML=T)
afl_full_INT <- lmerTest::lmer(log(AFL+((AFL^2+1)^0.5)) ~ extraction_method + year + tissue_type  + elution_vol + extraction_method*year + (1|ind), subset_33_1, REML=T)
afl_med_INT <- lmerTest::lmer(log(AFL+((AFL^2+1)^0.5)) ~ extraction_method + year + extraction_method*year + (1|ind), subset_33_1, REML=T)

# How to check results: 
# AIC: to select the model that best balances fit and complexity of model to data
# LRT: testing whether the interaction term (or any additional parameters) significantly improves the goodnes of fit of two nested models
# LRT (Likelihood Ratio Test) can be done using anova(model1, model2)

# 10.1) results for DNA_conc:
AIC(dnaconc_full)
AIC(dnaconc_full_INT)
# AIC points to model WITHOUT interaction

anova(dnaconc_full)
anova(dnaconc_full_INT)
anova(dnaconc_full, dnaconc_full_INT)
# interaction is significant, anova on both models points towards interaction

# 10.2) results for DNA_mg:
AIC(dnamg_full)
AIC(dnamg_full_INT)
# AIC points to model WITHOUT interaction

anova(dnamg_full)
anova(dnamg_full_INT)
anova(dnamg_full, dnamg_full_INT)
# interaction is not significant, anova on both models points against interaction

# 10.3) results for AFL full:
AIC(afl_full)
AIC(afl_full_INT)
# AIC points to model WITHOUT interaction

anova(afl_full)
anova(afl_full_INT)
anova(afl_full, afl_full_INT)
# interaction is significant, anova on both models points towards interaction

# 10.4) results for AFL med:
AIC(afl_med)
AIC(afl_med_INT)
# AIC points to model WITHOUT interaction

anova(afl_med)
anova(afl_med_INT)
anova(afl_med, afl_med_INT)

anova(afl_med_INT, afl_med)
# interaction is significant, anova on both models points towards interaction

# Conclusion: 
# AIC clearly prefers models without interaction term. 
# anova() shows that interaction term is significant for DNA concentration and AFL.
#----

# 9) ANOVA for DNA_conc
#*******************************
#----
# Results from LMM: 
# All models point to extraction_method and year being significant. 
# Models are pretty similar in their explanatory power, no significant difference, although the smallest one is the best one.

anov <- aov(DNA_conc~extraction_method + year + elution_vol + tissue_type, subset_33)

# Data needs to be transformed again, to meet assumptions of Homogeneity and Homoscedasticity
model <- lm(DNA_conc+0.00001~extraction_method + year + elution_vol + tissue_type, subset_33)
bc <- boxcox(model)
lambda <- bc$x[which.max(bc$y)]

anov <- aov(((DNA_conc^lambda-1)/lambda)~extraction_method + year + elution_vol + tissue_type, subset_33)

# assumptions
par(mfrow=c(2,2))
plot(anov)
par(mfrow=c(1,1))

# result
summary(anov)

# Tukey post-hoc test
tukey_results <- TukeyHSD(anov, which = "extraction_method", conf.level = 0.95)
tukey_results

plot(tukey_results, las = 1)

# Make a nice boxplot out of it
# get labels
tukey.levels <- tukey_results[["extraction_method"]][,4]
tukey.labels <- data.frame(multcompLetters(tukey.levels)['Letters'])

# add labels to data, to include in the plot later
plotdata <- subset(subset_33) %>% 
  group_by(extraction_method) %>%
  summarise(mean=mean(DNA_conc), quant = quantile(DNA_conc, probs = 0.75))

plotdata$letters = tukey.labels[match(plotdata$extraction_method, row.names(tukey.labels)),"Letters"]
plotdata

# flipped boxplot on DNA_conc on subset_33, no groups, cut off of exceptionally large DNA_conc
png("/ANOVA_DNAconc_subset33_full_new.png", 8000, 4000, res=1000, type='cairo-png', antialias=c("subpixel")) # to safe the file with anti-aliasing, tweak res and scale always together; I have no idea what anialias=c() does
subset(subset_33) %>%  
  filter(is.na(DNA_conc)==FALSE) %>% # no effect on plot, but it manages the warning message (removing NAs)
    ggplot(aes(x=extraction_method, y=DNA_conc)) +
      geom_boxplot(fill="lightblue") +
      scale_x_discrete(limits = rev(levels(subset_33$extraction_method)), labels=c("Qiagen_1"="Qiagen_OR", "Qiagen_2"="Qiagen_AD", "GENial"="Gen-IAL")) +
      coord_flip() +
      labs(title="", x="extraction protocol", y=expression("DNA Concentration (ng/μl)")) +
      geom_text(data = plotdata, aes(x = extraction_method, y = quant, label = letters), size = 6, vjust=-0.5, hjust=-1)
dev.off()


# FINAL RESULT:
#****************
# significant: extraction_method and year
# ANOVA not ideal, because we ignore IND as a variable and because even with transformation, the data still is somewhat whacky in terms of homogeneity and homoscedasticity
# ANOVA shows no diff between Qiagen_1 and Quiagen_2 and no diff between GENial and Promega
#----

# 10) ANOVA for DNA_per_mg
#*******************************
# Results from LMM: 
# All models point to extraction_method and year being significant. 
# Models are pretty similar in their explanatory power, not significantly different, although the smallest one is the best one.

anov <- aov(DNA_per_mg~extraction_method + year + tissue_type, subset_33)

# Data needs to be transformed again, to meet assumptions of Homogeneity and Homoscedasticity
anov <- aov(log(DNA_per_mg+((DNA_per_mg^2+1)^0.5))~extraction_method + year + tissue_type, subset_33)

# assumptions
par(mfrow=c(2,2))
plot(anov)
par(mfrow=c(1,1))

# result
summary(anov)

# Tukey post-hoc test
tukey_results <- TukeyHSD(anov, which = "extraction_method", conf.level = 0.95)
tukey_results

plot(tukey_results, las = 1)

# Make a nice boxplot out of it
# get labels
tukey.levels <- tukey_results[["extraction_method"]][,4]
tukey.labels <- data.frame(multcompLetters(tukey.levels, threshold=0.05)['Letters'])

# add labels to data, to include in the plot later
plotdata <- subset(subset_33) %>% 
  group_by(extraction_method) %>%
  summarise(mean=mean(DNA_per_mg), quant = quantile(DNA_per_mg, probs = 0.75))

plotdata$letters = tukey.labels[match(plotdata$extraction_method, row.names(tukey.labels)),"Letters"]
plotdata

# flipped boxplot on DNA_conc on subset_33, no groups, cut off of exceptionally large DNA_conc
# find three outliers
sort(subset_33$DNA_per_mg, T)[1:3]

png("/ANOVA_DNAquan_subset33_full_new.png", 8000, 4000, res=1000, type='cairo-png', antialias=c("subpixel")) # to safe the file with anti-aliasing, tweak res and scale always together; I have no idea what anialias=c() does
subset(subset_33) %>%
  filter(is.na(DNA_per_mg)==FALSE) %>% # no effect on plot, but it manages the warning message (removing NAs)
  ggplot(aes(x=extraction_method, y=DNA_per_mg)) +
  geom_boxplot(fill="lightblue") +
  scale_x_discrete(limits = rev(levels(subset_33$extraction_method)), labels=c("Qiagen_1"="Qiagen_OR", "Qiagen_2"="Qiagen_AD", "GENial"="Gen-IAL")) +
  scale_y_continuous(limits = c(0, 130)) +
  coord_flip() +
  labs(title="", x="extraction protocol", y="DNA yield (ng/mg tissue used") +
  geom_text(data = plotdata, aes(x = extraction_method, y = quant, label = letters), size = 6, vjust=-0.5, hjust=-1) +
  annotate("text", 4.4, 109, label="3 outliers at 196, 200 and 295", size=5) +
  geom_segment(
    x = 4.1, y = 120,
    xend = 4.1, yend = 135,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 0.5, 
    arrow = arrow(length = unit(0.1, "inches")),
    colour = "black")
dev.off()

# FINAL RESULT:
#****************
# significant: extraction_method and year
# ANOVA not ideal, because we ignore IND as a variable and because even with transformation, the data still is somewhat whacky in terms of homogeneity and homoscedasticity
# ANOVA shows no diff between Qiagen_1, Quiagen_2 and PCI, and no diff between GENial and Promega


# 11) ANOVA for AFL
#*******************************
# Results from LMM: 
# Full model seems to be the most valid, although medium model has lower AICc (but AICc is not good anyways)
# Full model shows extraction method to be explanatory, year only slightly
# no signifcant diference though between FULL and MEDIUM
# done on subset33_1 in LMM

anov <- aov(AFL~extraction_method + year + tissue_type + elution_vol, subset_33_1)

# Data needs to be transformed again, to meet assumptions of Homogeneity and Homoscedasticity
anov <- aov(log(AFL+((AFL^2+1)^0.5)) ~extraction_method + year + tissue_type + elution_vol, subset_33_1)

# assumptions
par(mfrow=c(2,2))
plot(anov)
par(mfrow=c(1,1))

# result
summary(anov)

# Tukey post-hoc test
tukey_results <- TukeyHSD(anov, which = "extraction_method", conf.level = 0.95)
tukey_results

plot(tukey_results, las = 1)

# Make a nice boxplot out of it
# get labels
tukey.levels <- tukey_results[["extraction_method"]][,4]
tukey.labels <- data.frame(multcompLetters(tukey.levels, threshold=0.05)['Letters'])

# add labels to data, to include in the plot later
plotdata <- subset_33_1 %>% 
  group_by(extraction_method) %>%
  summarise(mean=mean(AFL), quant = quantile(AFL, probs = 0.75))

plotdata$letters = tukey.labels[match(plotdata$extraction_method, row.names(tukey.labels)),"Letters"]
plotdata

# flipped boxplot on DNA_conc on subset_33, no groups, cut off of exceptionally large DNA_conc
png("/ANOVA_AFL_subset33_new.png", 8000, 4000, res=1000, type='cairo-png', antialias=c("subpixel")) # to safe the file with anti-aliasing, tweak res and scale always together; I have no idea what anialias=c() does
subset_33_1 %>%
  filter(is.na(AFL)==FALSE) %>% # no effect on plot, but it manages the warning message (removing NAs)
  ggplot(aes(x=extraction_method, y=AFL)) +
  geom_boxplot(fill="lightblue") +
  scale_x_discrete(limits = rev(levels(subset_33_1$extraction_method)), labels=c("Qiagen_1"="Qiagen_OR", "Qiagen_2"="Qiagen_AD", "GENial"="Gen-IAL")) +
  coord_flip() +
  labs(title="", x="extraction protocol", y="average fragmenth length (bp)") +
  geom_text(data = plotdata, aes(x = extraction_method, y = quant, label = letters), size = 6, vjust=-0.5, hjust=-1)
dev.off()


# FINAL RESULT:
#****************
# significant: extraction_method and year and tissue_type
# ANOVA not ideal, because we ignore IND as a variable
# ANOVA shows no diff between Qiagen_1, Quiagen_2 and PCI, and no diff between GENial and Promega


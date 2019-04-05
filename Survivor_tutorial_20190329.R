#GERMINATION DATA
#VIVIAN BERNAU
#6 Sept 2017
#Based on McNair et al 2012; Seed Science Research

#Using SURVMINER

#My data is right-sensored.
#My data is interval, but can probably be analyzed nonparametrically as exact for as long as plates with large losses are removed.

#Survivor function: probability that the germination time is greater than t

#Set working directory and repositories
wd <- ("~/Google Drive/RFiles/chile-germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

setwd(out.dir)

#read in germination data in pre-lifetab format
#df10 <- read.csv(paste(out.dir, "/cleaned10_2017-09-07.csv", sep = ""), header = T)
df <- read.csv(paste(out.dir, "/cleaned8_2017-10-13.csv", sep = ""), header = T, na.strings = c("", "NA"))

setwd(c:/Users/Audrey McCombs/Desktop)
df <- read.csv(paste(out.dir, "/cleaned8_2017-10-13.csv", sep = ""), header = T, na.strings = c("", "NA"))
str(df)
df$end <- round(df$end, digits = 1)
df$run <- as.factor(df$run)

df <- subset(df, region!="control")
df$cv <- as.numeric(df$region == "central valleys")
df$ecoast <- as.numeric(df$region == "ecoast")
df$wcoast <- as.numeric(df$region == "wcoast")
df$yucatan <- as.numeric(df$region == "yucatan")
df$ss <- as.numeric(df$region == "sierra madre")

library(survival)
library(mnormt)
library(survminer)
library(RColorBrewer)

#fit and plot Kaplan-Meier survivor function, includes 95% confidence intervals
test.peg <- survfit(Surv(end, status) ~ trt ,data = df, type = "kaplan-meier")
ggsurvplot(test.peg, data =df, conf.int = T, pval = T, 
           risk.table = F, xlab = "Time (h)", ylab = "Not Germinated",
           surv.median.line = "hv", legend = "bottom")
ggsave("peg.jpg", width = 6, height = 4, units = "in", dpi = 300)
survdiff.peg <- pairwise_survdiff(Surv(end, status)~trt, data = df, 
                                  p.adjust.method = "bonferroni",rho = 0)
sink("peg.txt")
print(survdiff.peg)
sink()

test.pegregion <- survfit(Surv(end, status) ~ trt + region,data = df, 
                          type = "kaplan-meier", conf.type = "log")
ggsurvplot(test.pegregion, data =df, conf.int = T, pval = T, 
           risk.table = F, xlab = "Time (h)", ylab = "Not Germinated",
           surv.median.line = "hv", legend = "bottom")
ggsave("peg_region.jpg", width = 6, height = 4, units = "in", dpi = 300)
survdiff.peg <- pairwise_survdiff(Surv(end, status)~trt+region, data = df, 
                                  p.adjust.method = "bonferroni",rho = 0)
sink("peg_region.txt")
print(survdiff.pegregion)
sink()



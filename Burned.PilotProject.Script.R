##############################Burned CA Sierras - Summer 2023#######################################

###Set Working Directory

#### * ALL PACKAGES USED HERE * ####
library(plyr)
library(ggplot2) 
library(gridExtra)
library(ggpubr) 
library(lme4)
library(multcomp)
library(car)
library(AICcmodavg)
library(lmtest)
library(dplyr) 
library(tidyr)
library(cowplot)
library(outliers)
library(MuMIn)
library(ggstatsplot)
library(ggbiplot)


#### * LOADING DATASET * ####
completeDataset <- read.csv("Seki.dataset.csv", header = T)

#excluding NA values for Body Size (wings measurement)
completeDataset1 <- subset(completeDataset, !is.na(Wing.mean))
#checking dataset
dim(completeDataset1)
head(completeDataset1)
tail(completeDataset1)
summary(completeDataset1)

####Worker body size by wings########

#Graph
bw <- ggplot(completeDataset, aes(reorder(Site,Wing.mean), y=Wing.mean, fill = Type, color = Type)) +
  geom_boxplot(alpha = 0.5) +
  facet_wrap(~Pair, scale="free")
bw <- bw + labs(x = " ", y = "Marginal cell length (cm)")
bw <- bw + theme_bw()
bw <- bw + theme(legend.position = "bottom", legend.title=element_blank())
bw <- bw + theme(axis.line = element_line(color='black'),
                 panel.grid.minor = element_blank())
bw <- bw + ggtitle("Worker body size per habitat")
bw <- bw + geom_jitter(shape=16, position=position_jitter(0.2)) +
  guides(colour = guide_legend(reverse=T), fill = guide_legend(reverse=T))
bw <- bw + scale_fill_grey(start = 0.6, end = 0.2) +
  scale_color_grey(start = 0.6, end = 0.2) 
bw <- bw + theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
bw

bwa <- ggplot(completeDataset, aes(reorder(Type,Wing.mean), y=Wing.mean, fill = Type, color = Type)) +
  geom_boxplot(alpha = 0.5)
bwa <- bwa + labs(x = " ", y = "Marginal cell length (cm)")
bwa <- bwa + theme_bw()
bwa <- bwa + theme(legend.position = "none")
bwa <- bwa + theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
bwa <- bwa + ggtitle("Worker body size")
bwa <- bwa + geom_jitter(shape=16, position=position_jitter(0.2))
bwa <- bwa + scale_fill_grey(start = 0.6, end = 0.2) +
  scale_color_grey(start = 0.6, end = 0.2) 
bwa


#Stats

completeDataset1 <- subset(completeDataset1, !is.na(Wing.mean))
mbw <- ddply(completeDataset1, "Type", summarise,
            No.bees=length(Wing.mean),
            grp.mean=mean(Wing.mean),
            var=var(Wing.mean),
            sd   = sd(Wing.mean, na.rm=TRUE),
            se   = sd / sqrt(sum(Wing.mean)),
            min=min(Wing.mean),
            max=max(Wing.mean))
mbw


#models
BWnull <- glmer(Wing.mean ~1  + (1| Elevation), data=completeDataset1, family = Gamma(link = "identity"))
BW1 <- glmer(Wing.mean ~ Site  + (1| Elevation), data=completeDataset1, family = Gamma(link = "identity"))
BW2 <- glmer(Wing.mean ~ Type  + (1| Elevation), data=completeDataset1, family = Gamma(link = "identity"))
BW3 <- glmer(Wing.mean ~ Type + Site  + (1| Elevation), data=completeDataset1, family = Gamma(link = "identity"))
BW4 <- glmer(Wing.mean ~ Type * Site  + (1| Elevation), data=completeDataset1, family = Gamma(link = "identity"))
BW4 <- glmer(Wing.mean ~ Type * Site  + (1| Elevation), data=completeDataset1, family = Gamma(link = "identity"))
BW5 <- glmer(Wing.mean ~  Vegetation + (1| Elevation), data=completeDataset1, family = Gamma(link = "identity"))
BW6 <- glmer(Wing.mean ~  Type + Vegetation + (1| Elevation), data=completeDataset1, family = Gamma(link = "identity"))
BW7 <- glmer(Wing.mean ~  Type * Vegetation + (1| Elevation), data=completeDataset1, family = Gamma(link = "identity"))
BW8 <- glmer(Wing.mean ~  Site + Vegetation + (1| Elevation), data=completeDataset1, family = Gamma(link = "identity"))
BW9 <- glmer(Wing.mean ~  Site * Vegetation + (1| Elevation), data=completeDataset1, family = Gamma(link = "identity"))

model.sel(BWnull, BW1, BW2, BW3, BW4, BW5, BW6, BW7, BW8, BW9) # choose the best fit model

summary(BW2)
Anova(BW2)
lrtest(BWnull, BW2)
summary(glht(BW2, linfct=mcp(Type="Tukey")))


####Worker body mass########

#Graph
bm <- ggplot(completeDataset, aes(reorder(Site,Wing.mean), y=Weight.mg., fill = Type, color = Type)) +
  geom_boxplot(alpha = 0.5) +
  facet_wrap(~Pair, scale="free")
bm <- bm + labs(x = "Site", y = "Weight (mg)")
bm <- bm + theme_bw()
bm <- bm + theme(legend.position = "bottom", legend.title=element_blank())
bm <- bm + theme(axis.line = element_line(color='black'),
                 panel.grid.minor = element_blank())
bm <- bm + ggtitle("Body mass")
bm <- bm + geom_jitter(shape=16, position=position_jitter(0.2)) +
  guides(colour = guide_legend(reverse=T), fill = guide_legend(reverse=T))
bm <- bm + scale_fill_grey(start = 0.6, end = 0.2) +
  scale_color_grey(start = 0.6, end = 0.2) 
bm

bma <- ggplot(completeDataset, aes(reorder(Type,Wing.mean), y=Weight.mg., fill = Type, color = Type)) +
  geom_boxplot(alpha = 0.5)
bma <- bma + labs(x = " ", y = "Weight (mg)")
bma <- bma + theme_bw()
bma <- bma + theme(legend.position = "none")
bma <- bma + ggtitle("Worker body mass")
bma <- bma + theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
bma <- bma + geom_jitter(shape=16, position=position_jitter(0.2))
bma <- bma + scale_fill_grey(start = 0.6, end = 0.2) +
  scale_color_grey(start = 0.6, end = 0.2) 
bma


#Stats

completeDataset2 <- subset(completeDataset, !is.na(Weight.mg.))
mbm <- ddply(completeDataset2, "Type", summarise,
             No.bees=length(Weight.mg.),
             grp.mean=mean(Weight.mg.),
             var=var(Weight.mg.),
             sd   = sd(Weight.mg., na.rm=TRUE),
             se   = sd / sqrt(sum(Weight.mg.)),
             min=min(Weight.mg.),
             max=max(Weight.mg.))
mbm


#models
BMnull <- glmer(Weight.mg. ~1  + (1| Elevation), data=completeDataset2, family = Gamma(link = "identity"))
BM1 <- glmer(Weight.mg. ~ Site  + (1| Elevation), data=completeDataset2, family = Gamma(link = "identity"))
BM2 <- glmer(Weight.mg. ~ Type  + (1| Elevation), data=completeDataset2, family = Gamma(link = "identity"))
BM3 <- glmer(Weight.mg. ~ Type + Site  + (1| Elevation), data=completeDataset2, family = Gamma(link = "identity"))
BM4 <- glmer(Weight.mg. ~ Type * Site  + (1| Elevation), data=completeDataset2, family = Gamma(link = "identity"))
BM4 <- glmer(Weight.mg. ~ Type * Site  + (1| Elevation), data=completeDataset2, family = Gamma(link = "identity"))
BM5 <- glmer(Weight.mg. ~  Vegetation + (1| Elevation), data=completeDataset2, family = Gamma(link = "identity"))
BM6 <- glmer(Weight.mg. ~  Type + Vegetation + (1| Elevation), data=completeDataset2, family = Gamma(link = "identity"))
BM7 <- glmer(Weight.mg. ~  Type * Vegetation + (1| Elevation), data=completeDataset2, family = Gamma(link = "identity"))
BM8 <- glmer(Weight.mg. ~  Site + Vegetation + (1| Elevation), data=completeDataset2, family = Gamma(link = "identity"))
BM9 <- glmer(Weight.mg. ~  Site * Vegetation + (1| Elevation), data=completeDataset2, family = Gamma(link = "identity"))

model.sel(BMnull, BM1, BM2, BM3, BM4, BM5, BM6, BM7, BM8, BM9) # choose the best fit model

summary(BM2)
Anova(BM2)
lrtest(BMnull, BM2)
summary(glht(BM2, linfct=mcp(Type="Tukey")))

####Figure Bee######
ggarrange(bw,                # First row with line plot
  # Second row with box and dot plots
  ggarrange(bwa, bma, ncol = 2, labels = c("B", "C")), 
  nrow = 2, 
  labels = "A")     # Label of the line plot


####Correlation########

ggscatter(completeDataset, x = "Wing.mean", y = "Weight.mg.", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Marginal cell length (cm)", ylab = "Weight (mg)")


c1 <- ggplot(completeDataset, aes(y = Weight.mg., x = Wing.mean))
c1 <- c1 + geom_smooth(method=lm, se=FALSE, color="black")
c1 <- c1 + geom_jitter(aes(color = Type)) + 
  geom_smooth(aes(color = Type, fill = Type), linetype="dashed", method = lm) +
  scale_fill_grey(start = 0.6, end = 0.2) + scale_color_grey(start = 0.6, end = 0.2) 
c1 <- c1 + labs(y = "Marginal cell length (cm)", x = "Weight (mg)")
c1 <- c1 + theme(legend.title = element_blank()) + theme(legend.position = "right")
c1 <- c1 + theme_classic()
c1 <- c1 + theme(text = element_text(size = 12)) 
c1


####Abundance bumble bees########


#### * LOADING DATASET * ####
site.Info <- read.csv("Site.Information.csv", header = T)

#checking dataset
dim(site.Info)
head(site.Info)
tail(site.Info)
summary(site.Info)

#Graph
AB <- ggplot(site.Info, aes(reorder(Type,Abundance.Bumblebee), y=Abundance.Bumblebee, fill = Type, color = Type)) +
  geom_boxplot(alpha = 0.5)
AB <- AB + labs(x = " ", y = "Bumble bee abundance")
AB <- AB + theme_bw()
AB <- AB + theme(legend.position = "none")
AB <- AB + ggtitle("Bumble bee abundance")
AB <- AB + scale_fill_grey(start = 0.6, end = 0.2) +
  scale_color_grey(start = 0.6, end = 0.2) 
AB <- AB + stat_summary(fun.y=mean, geom="point", shape = 17, size = 3) 
AB

AP <- ggplot(site.Info, aes(reorder(Type,Abundance.Honeybee), y=Abundance.Honeybee, fill = Type, color = Type)) +
  geom_boxplot(alpha = 0.5)
AP <- AP + labs(x = " ", y = "Honey bee abundance")
AP <- AP + theme_bw()
AP <- AP + theme(legend.position = "none")
AP <- AP + ggtitle("Honey bee abundance")
AP <- AP + scale_fill_grey(start = 0.6, end = 0.2) +
  scale_color_grey(start = 0.6, end = 0.2) 
AP <- AP + stat_summary(fun.y=mean, geom="point", shape = 17, size = 3) 
AP

#Stats

mab <- ddply(site.Info, "Type", summarise,
             No.Sites=length(Abundance.Bumblebee),
             grp.mean=mean(Abundance.Bumblebee),
             var=var(Abundance.Bumblebee),
             sd   = sd(Abundance.Bumblebee, na.rm=TRUE),
             se   = sd / sqrt(sum(Abundance.Bumblebee)),
             min=min(Abundance.Bumblebee),
             max=max(Abundance.Bumblebee))
mab


map <- ddply(site.Info, "Type", summarise,
             No.Sites=length(Abundance.Honeybee),
             grp.mean=mean(Abundance.Honeybee),
             var=var(Abundance.Honeybee),
             sd   = sd(Abundance.Honeybee, na.rm=TRUE),
             se   = sd / sqrt(sum(Abundance.Honeybee)),
             min=min(Abundance.Honeybee),
             max=max(Abundance.Honeybee))
map


Abnull <- glmer(Abundance.Bumblebee ~1 + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab1 <- glmer(Abundance.Bumblebee ~ Site + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab2 <- glmer(Abundance.Bumblebee ~ Type + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab3 <- glmer(Abundance.Bumblebee ~ Vegetation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab4 <- glmer(Abundance.Bumblebee ~ Abundance.Honeybee + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab5 <- glmer(Abundance.Bumblebee ~ Site + Type + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab6 <- glmer(Abundance.Bumblebee ~ Site * Type + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab7 <- glmer(Abundance.Bumblebee ~ Site + Abundance.Honeybee + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab8 <- glmer(Abundance.Bumblebee ~ Site * Abundance.Honeybee + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab9 <- glmer(Abundance.Bumblebee ~ Site + Vegetation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
# Ab10 <- glmer(Abundance.Bumblebee ~ Site * Vegetation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab11 <- glmer(Abundance.Bumblebee ~ Type + Vegetation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab12 <- glmer(Abundance.Bumblebee ~ Type * Vegetation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab13 <- glmer(Abundance.Bumblebee ~ Type + Abundance.Honeybee + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
# Ab14 <- glmer(Abundance.Bumblebee ~ Type * Abundance.Honeybee + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab15 <- glmer(Abundance.Bumblebee ~ Vegetation + Abundance.Honeybee + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
# Ab16 <- glmer(Abundance.Bumblebee ~ Vegetation * Abundance.Honeybee + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab17 <- glmer(Abundance.Bumblebee ~ Floral.Patch.Quality + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab18 <- glmer(Abundance.Bumblebee ~ Floral.Patch.Quality + Site + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
# Ab19 <- glmer(Abundance.Bumblebee ~ Floral.Patch.Quality * Site + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab20 <- glmer(Abundance.Bumblebee ~ Floral.Patch.Quality + Type + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab21 <- glmer(Abundance.Bumblebee ~ Floral.Patch.Quality * Type + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab22 <- glmer(Abundance.Bumblebee ~ Floral.Patch.Quality + Vegetation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
Ab23 <- glmer(Abundance.Bumblebee ~ Floral.Patch.Quality * Vegetation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
# Ab24 <- glmer(Abundance.Bumblebee ~ Floral.Patch.Quality + Abundance.Honeybee + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
# Ab25 <- glmer(Abundance.Bumblebee ~ Floral.Patch.Quality * Abundance.Honeybee + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))


model.sel(Abnull, Ab2, Ab3, Ab4, Ab5, Ab6, Ab7, Ab8, Ab9, Ab11, Ab12, Ab13, Ab15, Ab17, Ab18, Ab20, Ab21, Ab22, Ab23) # choose the best fit model

summary(Ab6)
Anova(Ab6)
lrtest(Abnull, Ab6)
summary(glht(Ab6, linfct=mcp(Type="Tukey")))


####Floral cover########

#Graph
FC <- ggplot(site.Info, aes(reorder(Type,Flora.Cover.Percentual), y=Flora.Cover.Percentual, fill = Type, color = Type)) +
  geom_boxplot(alpha = 0.5)
FC <- FC + labs(x = " ", y = "Relative percentage ")
FC <- FC + theme_bw()
FC <- FC + theme(legend.position = "none")
FC <- FC + ggtitle("Proportion of floral cover per transect")
FC <- FC + scale_fill_grey(start = 0.6, end = 0.2) +
  scale_color_grey(start = 0.6, end = 0.2) 
FC <- FC + stat_summary(fun.y=mean, geom="point", shape = 17, size = 3)
FC

RC <- ggplot(site.Info, aes(reorder(Type,Richness), y=Richness, fill = Type, color = Type)) +
  geom_boxplot(alpha = 0.5)
RC <- RC + labs(x = " ", y = "Floral richness")
RC <- RC + theme_bw()
RC <- RC + theme(legend.position = "none")
RC <- RC + ggtitle("Floral richness per transect")
RC <- RC + scale_fill_grey(start = 0.6, end = 0.2) +
  scale_color_grey(start = 0.6, end = 0.2) 
RC <- RC + stat_summary(fun.y=mean, geom="point", shape = 17, size = 3)
RC


#Stats

flora <- ddply(site.Info, "Type", summarise,
             No.Sites=length(Flora.Cover.Percentual),
             grp.mean=mean(Flora.Cover.Percentual),
             var=var(Flora.Cover.Percentual),
             sd   = sd(Flora.Cover.Percentual, na.rm=TRUE),
             se   = sd / sqrt(sum(Flora.Cover.Percentual)),
             min=min(Flora.Cover.Percentual),
             max=max(Flora.Cover.Percentual))
flora

FCnull <- glmer(Flora.Cover.Percentual ~1 + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
FC1 <- glmer(Flora.Cover.Percentual ~ Type + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
FC2 <- glmer(Flora.Cover.Percentual ~ Vegetation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
FC3 <- glmer(Flora.Cover.Percentual ~ Site + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
FC4 <- glmer(Flora.Cover.Percentual ~ Type + Vegetation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
FC5 <- glmer(Flora.Cover.Percentual ~ Type * Vegetation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
FC6 <- glmer(Flora.Cover.Percentual ~ Site + Type + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
# FC7 <- glmer(Flora.Cover.Percentual ~ Site * Type + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
# FC8 <- glmer(Flora.Cover.Percentual ~ Site + Vegatation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
# FC9 <- glmer(Flora.Cover.Percentual ~ Site * Vegatation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))

model.sel(FCnull, FC1, FC2, FC3, FC4, FC5, FC6) # choose the best fit model


res_aov <-aov(Flora.Cover.Percentual ~ Type,
    data = site.Info)

summary(res_aov)

Rnull <- glmer(Richness ~1 + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
R1 <- glmer(Richness ~ Type + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
R2 <- glmer(Richness ~ Vegetation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
R3 <- glmer(Richness ~ Site + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
R4 <- glmer(Richness ~ Type + Vegetation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
R5 <- glmer(Richness ~ Type * Vegetation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
R6 <- glmer(Richness ~ Site + Type + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
R7 <- glmer(Richness ~ Site * Type + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
# R8 <- glmer(Richness ~ Site + Vegatation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))
# R9 <- glmer(Richness ~ Site * Vegatation + (1| Elevation), data=site.Info, family = Gamma(link = "identity"))

model.sel(Rnull, R1, R2, R3, R4, R5, R6, R7) # choose the best fit model


flora1 <- ddply(site.Info, "Type", summarise,
               No.Sites=length(Richness),
               grp.mean=mean(Richness),
               var=var(Richness),
               sd   = sd(Richness, na.rm=TRUE),
               se   = sd / sqrt(sum(Richness)),
               min=min(Richness),
               max=max(Richness))
flora1


res_aov1 <-aov(Richness ~ Type,
              data = site.Info)

summary(res_aov1)

####Figure Flora######

ggarrange(
          ggarrange(FC, RC, ncol = 2, labels = c("A", "B")), 
          nrow = 1
          )     # Label of the line plot


#######################################################
#### Project: Metal contamination gradients        ####
####          in a tropical seagrass ecosystem     ####
#### Script purpose: Analysis of common metal data ####
#######################################################

#### 1.  Data preparation ####
#### 1.1 Load data ####
metals <- read.csv("~/Documents/MScProject/Data/CommonMetals.csv")

#### 1.2 Calculate distances between stations ####
require(geosphere)
# North
distN <- with(metals, distm(cbind(Longitude[1],Latitude[1]),
                            cbind(Longitude[1:63],Latitude[1:63]),
                            fun = distGeo))
# East
distE <- with(metals, distm(cbind(Longitude[64],Latitude[64]),
                            cbind(Longitude[64:84],Latitude[64:84]),
                            fun = distGeo))
# South
distS <- with(metals, distm(cbind(Longitude[85],Latitude[85]),
                            cbind(Longitude[85:161],Latitude[85:161]),
                            fun = distGeo))
# West
distW <- with(metals, distm(cbind(Longitude[162],Latitude[162]),
                            cbind(Longitude[162:170],Latitude[162:170]),
                            fun = distGeo))

metals$Distance <- c(distN, distE, distS, distW)
                   # add combined vector to dataframe

#### 2.  Data exploration and cleaning ####
#### 2.1 Fe ####
require(ggplot2)
ggplot(metals, aes(Distance, Fe.238.204, colour = Direction)) +
  geom_point() +
  # geom_text(aes(label = ID)) + # this line to see IDs
  facet_wrap(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

ggplot(metals, aes(Distance, Fe.259.941, colour = Direction)) +
  geom_point() +
  # geom_text(aes(label = ID)) + # this line to see IDs
  facet_wrap(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

# average both Fe into one
metals$Fe <- with(metals, (Fe.238.204 + Fe.259.941)/2)

# remove outliers
metals$Fe[c(2, 6, 45, 59, 77:81, 163)] <- NA

ggplot(metals, aes(Distance, Fe, colour = Direction)) +
  geom_point() +
  facet_wrap(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")


#### 2.2 Al ####
ggplot(metals, aes(Distance, Al.167.078, colour = Direction)) +
  geom_point() +
  # geom_text(aes(label = ID)) + # this line to see IDs
  facet_grid(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

# remove outliers
metals$Al.167.078[c(45, 59, 77:81)] <- NA

ggplot(metals, aes(Distance, Al.167.078, colour = Direction)) +
  geom_point() +
  facet_grid(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

#### 2.3 Create dataframes for analysis ####
grass <- metals[-c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:135, 160:184),
                -c(7:13, 15:17)]
grass <- grass[complete.cases(grass$Fe),]
rownames(grass) <- NULL

sed <- metals[c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:135, 160:170),
              -c(7:17)]
sed <- sed[complete.cases(sed$Fe),]
rownames(sed) <- NULL

both <- metals[-c(171:184), -c(7:17)]
both <- both[complete.cases(both$Fe),]
rownames(both) <- NULL

both$Material <- factor(both$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))

#### 2.4 Rename variables ####
gm <- grass$Material
gd <- grass$Distance
gsite <- grass$Direction
gFe <- grass$Fe
gAl <- grass$Al.167.078

sm <- sed$Material
sd <- sed$Distance
ssite <- sed$Direction
sFe <- sed$Fe

bm <- factor(both$Material)
bd <- both$Distance
bsite <- both$Direction
bFe <- both$Fe

#### 4.  Data analysis ####
#### 4.1 Fe ####
# use nonlinear least squares with self-start SSasymp function
# d = distance (metres)
# A = y asymptote that is approached by the curve
# C = the initial concentration at x = 0 (intercept) minus A
# k = the exponential decay rate per metre

?SSasymp
# self-starting function can only compute starting values for parameters
# when grass and sediment are treated separately because of large contrast
gFem <- nls(gFe ~ SSasymp(gd, A, C, logk))
summary(gFem) # parameters

sFem <- nls(sFe ~ SSasymp(sd, A, C, logk))
summary(sFem) # parameters

# The first model (gFem) fits one exponential decay equation to both seagrasses
# and the second model (sFem) fits one exponential decay equation to sediment.
# We want to have all three in one model that distinguishes between materials.
# To do that we use starting values supplied bz the self-starting models (note 
# that the starting values have to be dramatically different for sediment and seagrass).
m1 <- nls(bFe ~ C[bm] * exp(k[bm] * bd) + A[bm], 
          start = list(C = c(6284.9186, 345.3724, 345.3724), 
                       k = c(-exp(-2.4107), -exp(-3.4841), -exp(-3.4841)), 
                       A = c(559.1243, 116.0805, 116.0805)))
summary(m1)

# test model fit
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m1) ~ bd) # residual variance decreases dramatically with distance
boxplot(resid(m1) ~ bm) # and is smaller for seagrass than sediment
# homogeneity needs improving

hist(resid(m1))
qqnorm(resid(m1))
qqline(resid(m1)) # normality is fine; deviance at the margins but balanced
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

# improve homogeneity with generslised nonlinear least squares (gnls) 
# note that gnls is extremely sensitive and only converges when we use
# the data argument and the species factor has to be incorporated differently
# do the same as above but with a varPower weighting for homogeneity
require(nlme)
m2 <- gnls(Fe ~ C * exp(k * Distance) + A,
           start = list(C = c(6284.9186, 345.3724, 345.3724), 
                        k = c(-exp(-2.4107), -exp(-3.4841), -exp(-3.4841)), 
                        A = c(559.1243, 116.0805, 116.0805)),
           params = list(C ~ Material, k ~ Material, A ~ Material),
           weights = varPower(),
           data = both)

summary(m2)

# test model fit
plot(m2, col = bm)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m2, type = "normalized") ~ bd) # residual variance no longer differs with distance
boxplot(resid(m2, type = "normalized") ~ bm) # or material
# homogeneity is improved

hist(resid(m2))
qqnorm(resid(m2))
qqline(resid(m2)) # normality is unchanged 
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m2 is chosen as the optimal model based on improved homogeneity
# and reasonaly balanced residual variance albeit a lack of normality

# interpret model
summary(m2)
# y = C * e^k*x + A  (always add the value of the intercept)

# Sediment
# y = 5450.663 * e^(-0.082x) + 539.363
# Enhalus
# y = 276.637 * e^(-0.025x) + 109.453
# Thalassia
# y = 179.928 *e^(-0.033x) + 119.702

both$Material <- factor(both$Material, levels = c("Enhalus acoroides", "Sediment", "Thalassia hemprichii"))
# relevel factor levels

m2 <- gnls(Fe ~ C * exp(k * Distance) + A,
           start = list(C = c(345.3724, 6284.9186, 345.3724), 
                        k = c(-exp(-3.4841), -exp(-2.4107), -exp(-3.4841)), 
                        A = c(116.0805, 559.1243, 116.0805)),
           params = list(C ~ Material, k ~ Material, A ~ Material),
           weights = varPower(),
           data = both) # run model again (note that the starting values were changed to the
# output given by the previous model since the re-leveled model would otherwise not converge)

summary(m2)
# C
# Enhalus vs. Sediment, t = 5.150858, p = 0.0000 (< 0.001) different at alpha = 0.05
# Enhalus vs. Thalassia, t = -2.391601, p = 0.0180 (= 0.02) different at alpha = 0.05 but not alpha = 0.01
# Sediment vs. Thalassia, t = -5.248059, p = 0e+00 (< 0.001) different at alpha = 0.05

# k
# Enhalus vs. Sediment, t = -4.517140, p = 0.0001 (< 0.001) different at alpha = 0.05
# Enhalus vs. Thalassia, t = -1.013081, p = 0.3126 (= 0.31) NOT different at alpha = 0.05
# Sediment vs. Thalassia, t = 3.511786, p = 6e-04 (< 0.001) different at alpha = 0.05

# A
# Enhalus vs. Sediment, t = 13.233045, p = 0.0000 (< 0.001) different at alpha = 0.05
# Enhalus vs. Thalassia, t = 1.805591, p =  0.0730 (= 0.07) NOT different at alpha = 0.05
# Sediment vs. Thalassia, t = -12.898970, p = 0e+00 (< 0.001) different at alpha = 0.05

# alpha = 0.05, 95% CI, percentiles: 0.025-0.975
# alpha = 0.01, 99% CI, percentiles: 0.005-0.995

#### 4.2  Al ####
gAlm <- nls(gAl ~ SSasymp(gd, A, C, logk)) # note that there are no Al sediment data since
# it was below detection limit, hence only a seagrass model
summary(gAlm) # parameters

both$Material <- factor(both$Material, levels = c("Enhalus acoroides", "Sediment", "Thalassia hemprichii"))

m3 <- gnls(Al.167.078 ~ C * exp(k * Distance) + A, 
           start = list(C = c(123.449, 123.449), 
                        k = c(-exp(-4.573), -exp(-4.573)), 
                        A = c(41.314, 41.314)),
           params = list(C ~ Material, k ~ Material, A ~ Material),
           weights = varPower(),
           data = grass)

summary(m3)

# test model fit
plot(m3, pch = gm)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m3, type = "normalized") ~ gd)
boxplot(resid(m3, type = "normalized") ~ gm)
# homogeneity is good

hist(resid(m3))
qqnorm(resid(m3))
qqline(resid(m3)) # normality is perfect 
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m3 is chosen as the optimal model

# interpret model
summary(m3)
# y = C * e^k*x + A
# Enhalus
# y = 128.25643 * e^(-0.01257x) + 42.97479
# Thalassia
# y = 53.55757 *e^(-0.01633x) + 44.70129

# C
# Enhalus vs. Thalassia, t = -4.606882, p = 0.0000 (< 0.001) different at alpha = 0.05

# k
# Enhalus vs. Thalassia, t = -0.412404, p = 0.6809 (= 0.68) different at alpha = 0.05

# A
# Enhalus vs. Thalassia, t = 0.201855, p =  0.8404 (= 0.84) different at alpha = 0.05


#### 5.  Data visualisation ####
#### 5.1 Calculate descriptive statistics and model predictions ####
#### 5.1.1 Fe ####
require(psych)
sFe.stat <- describeBy(sFe, sd, mat = T) 
sFe.stat$group1 <- as.numeric(sFe.stat$group1) # convert to number

gFe.stat <- describeBy(gFe, list(gm, gd), mat = T) 
gFe.stat$group2 <- as.numeric(gFe.stat$group2) # convert to number


Fenew <- data.frame(Material = c(rep("Enhalus acoroides", 3265),
                                 rep("Sediment", 3579),
                                 rep("Thalassia hemprichii", 3296)),
                    Distance = c(seq(bd[83], bd[50], by = 0.1), # 0.1-m increments
                                 seq(bd[1], bd[159], by = 0.1),
                                 seq(bd[65], bd[56], by = 0.1)))

Fenew$fit <- predict(m2, newdata = Fenew)

# bootstrap confidence interval
bootfun <- function(newdata) {
  start <- coef(m2)
  boot <- both[sample(nrow(both), size = nrow(both), replace = TRUE),]
  bootfit <- try(update(m2,
                        start = start,
                        data = boot),
                 silent = TRUE)
  if (inherits(bootfit, "try-error")) return(rep(NA, nrow(newdata)))
  predict(bootfit, newdata)
}

bmat <- replicate(1000, bootfun(Fenew))
Fenew$lwr <- apply(bmat, 1, quantile, 0.025, na.rm = TRUE)
Fenew$upr <- apply(bmat, 1, quantile, 0.975, na.rm = TRUE)

# split new dataframe into grass and sediment
Fenewgrass <- Fenew[c(1:3265, 6845:10140),]
Fenewsed <- Fenew[3266:6844,]

#### 5.1.2 Al ####
gAl.stat <- describeBy(gAl, list(gm, gd), mat = T) 
gAl.stat$group2 <- as.numeric(gAl.stat$group2) # convert to number

# aggregate(Distance ~ Material, data = grass, min)

Alnew <- data.frame(Material = c(rep("Enhalus acoroides", 3265),
                                 rep("Thalassia hemprichii", 3296)),
                    Distance = c(seq(10.871910, 337.3496, by = 0.1), # 0.1-m increments
                                 seq(7.847574, 337.3496, by = 0.1)))

Alnew$fit <- predict(m3, newdata = Alnew)

# bootstrap confidence interval
bootfun <- function(newdata) {
  start <- coef(m3)
  boot <- grass[sample(nrow(grass), size = nrow(grass), replace = TRUE),]
  bootfit <- try(update(m3,
                        start = start,
                        data = boot),
                 silent = TRUE)
  if (inherits(bootfit, "try-error")) return(rep(NA, nrow(newdata)))
  predict(bootfit, newdata)
}

bmat <- replicate(1000, bootfun(Alnew))
Alnew$lwr <- apply(bmat, 1, quantile, 0.025, na.rm = TRUE)
Alnew$upr <- apply(bmat, 1, quantile, 0.975, na.rm = TRUE)

#### 5.2 Plot ####
#### 5.2.1 Customise theme ####
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .4, .2, .2),"cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 12, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title = element_blank(),
                 text = element_text(family = "Helvetica Neue"))

#### 5.2.2 Fe ####
gFep <- ggplot() +
        geom_line(data = Fenewgrass, aes(Distance, fit, colour = Material), size = 0.5) +
        geom_ribbon(data = Fenewgrass, aes(Distance, ymin = lwr, ymax = upr, fill = Material),
                    alpha = 0.5) +
        geom_pointrange(data = gFe.stat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                            colour = group1), size = 0.5) +
        # geom_point(data = grass, aes(Distance, Fe, colour = Material), size = 1.5) +
        scale_colour_manual(values = c("#2b491e", "#81a512"),
                            labels = c(expression(italic("Enhalus acoroides")),
                                       expression(italic("Thalassia hemprichii"))),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#2b491e", "#81a512"),
                          labels = c(expression(italic("Enhalus acoroides")),
                                     expression(italic("Thalassia hemprichii"))),
                          guide = guide_legend()) +
        annotate("text", x = c(60, 60), y = c(385, 360), size = 4.2, family = "Helvetica Neue",
                 hjust = 0, label = c("y == 277*e^{-0.03*x} + 109", "y == 180*e^{-0.03*x} + 120"),
                 parse = T) +
        ylab(expression("Iron concentration ("*mu*"g g"^-1*")")) +
        xlab("Distance from coast (m)") +
        coord_cartesian(ylim = c(50, 400), xlim = c(0, 350)) +
        scale_x_continuous(expand = c(0, 0), breaks = seq(0, 350, by = 50)) +
        scale_y_continuous(expand = c(0, 0), breaks = seq(50, 400, by = 50)) +
        theme(legend.position = c(.8, .93)) +
        mytheme
gFep # dimensions: 4 x 5 in


sFep <- ggplot() +
        geom_line(data = Fenewsed, aes(Distance, fit), size = 0.5, colour = "#d9b365") +
        geom_ribbon(data = Fenewsed, aes(Distance, ymin = lwr, ymax = upr),
                    alpha = 0.5, fill = "#d9b365") +
        geom_pointrange(data = sFe.stat, aes(group1, mean, ymin = mean - se, ymax = mean + se),
                        size = 0.5, colour = "#d9b365") +
        # geom_point(data = sed, aes(Distance, Fe), size = 1.5, colour = "#d9b365") +
        annotate("text", x = 240, y = 8500, size = 4.2, family = "Helvetica Neue",
                 hjust = 0, label = "y == 5451*e^{-0.08*x} + 539",
                 parse = T) +
        ylab(expression("Iron concentration ("*mu*"g g"^-1*")")) +
        xlab("Distance from coast (m)") +
        coord_cartesian(ylim = c(0, 9000), xlim = c(5, 381.2)) +
        scale_x_continuous(breaks = seq(0, 400, by = 50)) +
        scale_y_continuous(expand = c(0, 0), breaks = seq(0, 9000, by = 1500)) +
        mytheme
sFep # dimensions: 4 x 5 in





means <- with(both, describeBy(Fe, list(Station, Material), mat = T))
rownames(means) <- NULL
means[c(5, 9, 11, 16, 20, 25, 29),] <- NA
means <- means[complete.cases(means$mean),]
rownames(means) <- NULL

cor <- data.frame(Environmental = rep(means$mean[means$group2 == "Sediment"], 2),
                  Seagrass = means$mean[means$group2 %in% c("Enhalus acoroides", "Thalassia hemprichii")],
                  Species = c(rep("Enhalus acoroides", 5), rep("Thalassia hemprichii", 5)))


plot(cor$Seagrass ~ cor$Environmental, pch = cor$Species)

ggplot(cor, aes(Environmental, Seagrass, colour = Species, fill = Species)) +
  geom_point(size = 2) +
  # geom_smooth(method = "loess") +
  scale_colour_manual(values = c("#2b491e", "#81a512"),
                      labels = c(expression(italic("Enhalus acoroides")),
                                 expression(italic("Thalassia hemprichii"))),
                      guide = guide_legend()) +
  scale_fill_manual(values = c("#2b491e", "#81a512"),
                    labels = c(expression(italic("Enhalus acoroides")),
                               expression(italic("Thalassia hemprichii"))),
                    guide = guide_legend()) +
  theme(legend.position = c(.8, .93)) +
  mytheme

#### 5.2.3 Al plot ####
gAlp <- ggplot() +
        geom_line(data = Alnew, aes(Distance, fit, colour = Material, lty = Material), size = 0.5) +
        geom_ribbon(data = Alnew, aes(Distance, ymin = lwr, ymax = upr, fill = Material),
                    alpha = 0.5) +
        geom_pointrange(data = gAl.stat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                            colour = group1), size = 0.5) +
        # geom_point(data = grass, aes(Distance, Fe, colour = Material), size = 1.5) +
        scale_colour_manual(values = c("#2b491e", "#81a512"),
                            labels = c(expression(italic("Enhalus acoroides")),
                                       expression(italic("Thalassia hemprichii"))),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#2b491e", "#81a512"),
                          labels = c(expression(italic("Enhalus acoroides")),
                                     expression(italic("Thalassia hemprichii"))),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(1, 5),
                        guide = "none") +
        annotate("text", x = c(70, 70), y = c(192, 177), size = 4.2, family = "Helvetica Neue",
                 hjust = 0, label = c("y == 128*e^{-0.01*x} + 43", "y == 54*e^{-0.02*x} + 44"),
                 parse = T) +
        ylab(expression("Aluminium concentration ("*mu*"g g"^-1*")")) +
        xlab("Distance from coast (m)") +
        coord_cartesian(ylim = c(0, 200), xlim = c(0, 350)) +
        scale_x_continuous(expand = c(0, 0), breaks = seq(0, 350, by = 50)) +
        scale_y_continuous(expand = c(0, 0), breaks = seq(0, 200, by = 25)) +
        theme(legend.position = c(.8, .93)) +
        mytheme
gAlp # dimensions: 4 x 5 in

#### 5.2.4 Combined plot ####
require(cowplot)
plots <- align_plots(gFep, sFep, gAlp, align = "v", axis = "l")
final <- plot_grid(plots[[1]], plots[[2]], plots[[3]], labels = "auto",
                   label_size = 15, ncol = 1)
final # dimensions: 12 x 5 in

ggdraw(final) +
  draw_image("~/Documents/MScProject/Figures/Enhalus.png",
             x = -0.2, y = 0.4, scale = 0.12) +
  draw_image("~/Documents/MScProject/Figures/Thalassia.png",
             x = -0.27, y = 0.25, scale = 0.08) +
  draw_image("~/Documents/MScProject/Figures/Sediment.png",
             x = -0.2, y = 0.14, scale = 0.12) +
  draw_image("~/Documents/MScProject/Figures/Enhalus.png",
             x = -0.12, y = -0.27, scale = 0.12) +
  draw_image("~/Documents/MScProject/Figures/Thalassia.png",
             x = -0.24, y = -0.41, scale = 0.08) # dimensions: 12 x 5 in


#### 5.  Clean up ####
detach(package:geosphere) # detach geosphere package
detach(package:ggplot2) # detach ggplot2 package
detach(package:fitdistrplus) # detach fitdistrplus package
detach(package:car) # detach car package
detach(package:psych) # detach psych package
detach(package:nlme) # detach nlme package
detach(package:cowplot) # detach cowplot package
rm(list = ls()) # clear environment
graphics.off() # clear plots
cat("\014") # clear console


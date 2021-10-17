#####################################################
### Project: Metal contamination gradients        ###
### in a tropical seagrass ecosystem              ###
### Script purpose: Analysis of trace metal data  ###
#####################################################

#### 1.  Data preparation ####
#### 1.1 Load data ####
metals <- read.csv("~/Documents/MScProject/Data/TraceMetals.csv")

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
                            cbind(Longitude[85:159],Latitude[85:159]),
                            fun = distGeo))
# West
distW <- with(metals, distm(cbind(Longitude[160],Latitude[160]),
                            cbind(Longitude[160:168],Latitude[160:168]),
                            fun = distGeo))

metals$Distance <- c(distN, distE, distS, distW) # add combined vector to dataframe

#### 1.3.  Data exploration and cleaning ####
#### 2. Pb ####
require(ggplot2)
ggplot(metals, aes(Distance, Pb208, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) + 
  facet_grid(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

# remove outliers
metals$Pb208[c(2, 3, 47, 53, 62, 63, 67:73, 77:81, 96)] <- NA

ggplot(metals, aes(Distance, Pb208, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) + 
  facet_grid(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

#### 2.1. Clean Pb dataframe for analysis ####
grass <- metals[-c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:133, 158:168),]
grass <- grass[complete.cases(grass$Pb208),]
rownames(grass) <- NULL

sed <- metals[c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:133, 158:168),]
sed <- sed[complete.cases(sed$Pb208),]
rownames(sed) <- NULL

both <- metals
both <- both[complete.cases(both$Pb208),]
rownames(both) <- NULL

both$Material <- factor(both$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))

#### 2.2. Pb Rename variables ####
gm <- grass$Material
gd <- grass$Distance
gsite <- grass$Direction
gPb <- grass$Pb208

sm <- sed$Material
sd <- sed$Distance
ssite <- sed$Direction
sPb <- sed$Pb208

bm <- factor(both$Material)
bd <- both$Distance
bsite <- both$Direction
bPb<- both$Pb208

#### 2.3.  Pb Data analysis ####
# nonlinear least squares with self-start SSasymp function
# d = distance (metres)
# A = y asymptote that is approached by the curve
# C = the initial concentration at x = 0 (intercept) minus A
# k = the exponential decay rate per metre

?SSasymp
# self-starting function can only compute starting values for parameters
# when grass and sediment are treated separately because of large contrast
gSm <- nls(gPb ~ SSasymp(gd, A, C, logk))
summary(gSm) # parameters seagrass

sSm <- nls(sPb ~ SSasymp(sd, A, C, logk))
summary(sSm) # parameters sediment

# want to distinguish between materials 
# use start values of the manual model (note that there need to be two start
# values for each parameter because there are two materials)
# note that the starting values are different for sediment and seagrass
m1 <- nls(bPb ~ C[bm] * exp(k[bm] * bd) + A[bm], 
          start = list(C = c(13.3866, 9.31084, 9.31084), 
                       k = c(-exp(-3.2057), -exp(-2.60245),-exp(-2.60245)), 
                       A = c(0.7808, 0.55325, 0.55325)))
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

# improve homogeneity with with generalized nonlinear least squares (gnls)
# note that gnls is extremely sensitive and only converges when 
# we use the data argument and the species factor has to be incorporated differently
# same as above but with varPower weighting for homogeneity
require(nlme)
m2 <- gnls(Pb208 ~ C * exp(k * Distance) + A, 
           start = list(C = c(12.605800, 8.949594, 8.535499), 
                        k = c(-0.040532, -0.077169 , -0.070883), 
                        A = c(0.780828, 0.539473, 0.567029)),
           params = list(C ~ Material, k ~ Material, A ~ Material),
           weights = varPower(),
           data = both[,c(6, 14, 15)])

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
# m2 is chosen as the optimal model

# interpret model
summary(m2)
# y = C * e^k*x + A  (always add the value of the intercept-Sediment)
# Sediment
# y = 12.605228 * e^(-0.040996x) + 0.801133
# Enhalus
# y = 4.841476 * e^(-0.035155x) + 0.398268
# Thalassia
# y = 5.336361 *e^(-0.037783x) + 0.438473

# relevel
both$Material <- factor(both$Material, levels = c("Enhalus acoroides", "Sediment", "Thalassia hemprichii"))

m2 <- gnls(Pb208 ~ C * exp(k * Distance) + A, 
           start = list(C = c(8.949594, 12.605800, 8.535499), 
                        k = c(-0.077169, -0.040532, -0.070883), 
                        A = c(0.539473, 0.780828, 0.567029)),
           params = list(C ~ Material, k ~ Material, A ~ Material),
           weights = varPower(),
           data = both[,c(6, 14, 15)])

summary(m2)

# relevel
both$Material <- factor(both$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))
m2 <- gnls(Pb208 ~ C * exp(k * Distance) + A, 
           start = list(C = c(8.949594, 12.605800, 8.535499), 
                        k = c(-0.077169, -0.040532, -0.070883), 
                        A = c(0.539473, 0.780828, 0.567029)),
           params = list(C ~ Material, k ~ Material, A ~ Material),
           weights = varPower(),
           data = both[,c(6, 14, 15)])

#### 2.4. Pb descriptive statistics ####
require(psych)
gPb.stat <- describeBy(gPb, list(gm, gd), mat = T) 
gPb.stat$group2 <- as.numeric(gPb.stat$group2) # convert to number

sPb.stat <- describeBy(sPb, sd, mat = T) 
sPb.stat$group1 <- as.numeric(sPb.stat$group1) # convert to number

bPb.stat <- with(both, describeBy(Pb208, list(Material, Distance), mat = T))
bPb.stat$group2 <- as.numeric(bPb.stat$group2)

Pbnew <- data.frame(Material = c(rep("Sediment", 3579),
                                 rep("Enhalus acoroides", 3265),
                                 rep("Thalassia hemprichii", 3265)),
                    Distance = c(seq(bd[1], bd[149], by = 0.1),
                                 seq(bd[73], bd[53], by = 0.1),
                                 seq(bd[78], bd[57], by = 0.1)))

Pbnew$Material <- factor(Pbnew$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))
Pbnew$fit <- predict(m2, newdata = Pbnew)

# bootstrap confidence interval
bootfun <- function(newdata) {
  start <- coef(m2)
  boot <- both[,c(6, 14, 15)][sample(nrow(both[,c(6, 14, 15)]), 
                                     size = nrow(both[,c(6, 14, 15)]), replace = TRUE),]
  bootfit <- try(update(m2,
                        start = start,
                        data = boot),
                 silent = TRUE)
  if (inherits(bootfit, "try-error")) return(rep(NA, nrow(newdata)))
  predict(bootfit, newdata)
}

bmat <- replicate(1000, bootfun(Pbnew))
Pbnew$lwr <- apply(bmat, 1, quantile, 0.025, na.rm = TRUE)
Pbnew$upr <- apply(bmat, 1, quantile, 0.975, na.rm = TRUE)

#### 2.5. Pb plot ####
# Customise theme
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
# plot
bPb <- ggplot() +
  geom_line(data = Pbnew, aes(Distance, fit, colour = Material), size = 0.5) +
  geom_ribbon(data = Pbnew, aes(Distance, ymin = lwr, ymax = upr, fill = Material),
              alpha = 0.5) +
  #geom_rug(data = grass, aes(x = Distance, y = Pb208), sides = "b") +
  geom_pointrange(data = bPb.stat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                       colour = group1), size = 0.5) +
  # geom_point(data = both, aes(Distance, Pb208, colour = Material), size = 2.5) +
  scale_colour_manual(values = c("#d9b365", "#2b491e", "#81a512"),
                      labels = c("Sediment",
                                 expression(italic("Enhalus acoroides")),
                                 expression(italic("Thalassia hemprichii"))),
                      guide = guide_legend()) +
  scale_fill_manual(values = c("#d9b365", "#2b491e", "#81a512"),
                    labels = c("Sediment",
                               expression(italic("Enhalus acoroides")),
                               expression(italic("Thalassia hemprichii"))),
                    guide = guide_legend()) +
  annotate("text", x = c(70, 70, 70), y = c(14.3, 13.3, 12.3), size = 4.2,
           label = c("y == 12.61*e^{-0.04*x} + 0.8",
                     "y == 4.84*e^{-0.04*x} + 0.4",
                     "y == 5.34*e^{-0.04*x} + 0.44"),
           hjust = 0, parse = T) +
  ylab(expression("Lead concentration ("*mu*"g g"^-1*")")) +
  xlab("Distance from coast (m)") +
  coord_cartesian(ylim = c(0, 15), xlim = c(-10, 400)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 400, by = 50)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = c(.8, .9)) +
  mytheme
bPb # dimensions: 4 x 5 in


#### 3. Co ####
require(ggplot2)
ggplot(metals, aes(Distance, Co59, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) + 
  facet_grid(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

# remove outliers
metals$Co59[c(2, 47, 60, 68:73, 77:81, 89, 161)] <- NA

ggplot(metals, aes(Distance, Co59, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) + 
  facet_grid(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

#### 3.1. Clean Co dataframe for analysis ####
grass <- metals[-c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:133, 158:168),]
grass <- grass[complete.cases(grass$Co59),]
rownames(grass) <- NULL

sed <- metals[c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:133, 158:168),]
sed <- sed[complete.cases(sed$Co59),]
rownames(sed) <- NULL

both <- metals
both <- both[complete.cases(both$Co59),]
rownames(both) <- NULL

both$Material <- factor(both$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))
 
#### 3.2. Co Rename variables ####
gm <- grass$Material
gd <- grass$Distance
gsite <- grass$Direction
gCo <- grass$Co59

sm <- sed$Material
sd <- sed$Distance
ssite <- sed$Direction
sCo <- sed$Co59

bm <- factor(both$Material)
bd <- both$Distance
bsite <- both$Direction
bCo<- both$Co59

#### 3.3.  Co Data analysis ####
# nonlinear least squares with self-start SSasymp function
# d = distance (metres)
# A = y asymptote that is approached by the curve
# C = the initial concentration at x = 0 (intercept) minus A
# k = the exponential decay rate per metre

?SSasymp
# self-starting function can only compute starting values for parameters
# when grass and sediment are treated separately because of large contrast
gSm <- nls(gCo ~ SSasymp(gd, A, C, logk))
summary(gSm) # parameters seagrass

sSm <- nls(sCo ~ SSasymp(sd, A, C, logk))
summary(sSm) # parameters sediment

# we want to distinguish between materials 
# use start values of the manual model (note that there need to be two different start
# values for each parameter because there are two materials)
m3 <- nls(bCo ~ C[bm] * exp(k[bm] * bd) + A[bm], 
          start = list(C = c(0.85074, 0.65494, 0.65494), 
                       k = c(-exp(-2.73958), -exp(-3.77358), -exp(-3.77358)), 
                       A = c(0.30963, 0.22166, 0.22166)))
summary(m3)

# test model fit
plot(m3, col = bm)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m3) ~ bd) # residual variance decreases a bit with distance
boxplot(resid(m3) ~ bm) 

hist(resid(m3))
qqnorm(resid(m3))
qqline(resid(m3)) # normality is fine; deviance at the margins but balanced
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

# improve homogeneity with gnls
# with a varPower weighting for homogeneity

require(nlme)
m4 <- gnls(Co59 ~ C * exp(k * Distance) + A, 
           start = list(C = c(0.514704, 0.541104, 0.384654), 
                        k = c(-0.032802, -0.064597, -0.016848), 
                        A = c(0.197428, 0.309634, 0.247088)),
           params = list(C ~ Material, k ~ Material, A ~ Material),
           weights = varPower(),
           data = both[,c(6, 9, 15)])

summary(m4)

# test model fit
plot(m4, col = bm)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m4, type = "normalized") ~ bd) # residual variance no longer differs with distance
boxplot(resid(m4, type = "normalized") ~ bm) # or material
# homogeneity is improved

hist(resid(m4))
qqnorm(resid(m4))
qqline(resid(m4)) # normality is unchanged 
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m4 is chosen as the optimal model

# interpret model
summary(m4)
# y = C * e^k*x + A  (always add the value of the intercept)
# Sediment
# y = 0.534736 * e^(-0.0622622x) + 0.3084373
# Enhalus
# y = 0.4705847* e^(-0.0272828x) + 0.1928064
# Thalassia
# y = 0.3736216 *e^(-0.0153745x) + 0.244171

# relevel
both$Material <- factor(both$Material, levels = c("Enhalus acoroides", "Sediment", "Thalassia hemprichii"))
m4 <- gnls(Co59 ~ C * exp(k * Distance) + A, 
           start = list(C = c(0.514704, 0.541104, 0.384654), 
                        k = c(-0.032802, -0.064597, -0.016848), 
                        A = c(0.197428, 0.309634, 0.247088)),
           params = list(C ~ Material, k ~ Material, A ~ Material),
           weights = varPower(),
           data = both[,c(6, 9, 15)])

summary(m4)

# relevel
both$Material <- factor(both$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))
m4 <- gnls(Co59 ~ C * exp(k * Distance) + A, 
           start = list(C = c(0.514704, 0.541104, 0.384654), 
                        k = c(-0.032802, -0.064597, -0.016848), 
                        A = c(0.197428, 0.309634, 0.247088)),
           params = list(C ~ Material, k ~ Material, A ~ Material),
           weights = varPower(),
           data = both[,c(6, 9, 15)])

summary(m4)

#### 3.4. Co descriptive statistics ####
require(psych)
gCo.stat <- describeBy(gCo, list(gm, gd), mat = T) 
gCo.stat$group2 <- as.numeric(gCo.stat$group2) # convert to number

sCo.stat <- describeBy(sCo, sd, mat = T) 
sCo.stat$group1 <- as.numeric(sCo.stat$group1) # convert to number

bCo.stat <- with(both, describeBy(Co59, list(Material, Distance), mat = T))
bCo.stat$group2 <- as.numeric(bCo.stat$group2)

Conew <- data.frame(Material = c(rep("Sediment", 3579),
                                 rep("Enhalus acoroides", 3265),
                                 rep("Thalassia hemprichii", 3265)),
                    Distance = c(seq(bd[1], bd[152], by = 0.1),
                                 seq(bd[77], bd[55], by = 0.1), #increment
                                 seq(bd[83], bd[60], by = 0.1)))

Conew$Material <- factor(Conew$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))
Conew$fit <- predict(m4, newdata = Conew)

# bootstrap confidence interval
bootfun <- function(newdata) {
  start <- coef(m4)
  boot <- both[,c(6, 9, 15)][sample(nrow(both[,c(6, 9, 15)]), 
                                     size = nrow(both[,c(6, 9, 15)]), replace = TRUE),]
  bootfit <- try(update(m4,
                        start = start,
                        data = boot),
                 silent = TRUE)
  if (inherits(bootfit, "try-error")) return(rep(NA, nrow(newdata)))
  predict(bootfit, newdata)
}

bmat <- replicate(1000, bootfun(Conew))
Conew$lwr <- apply(bmat, 1, quantile, 0.025, na.rm = TRUE)
Conew$upr <- apply(bmat, 1, quantile, 0.975, na.rm = TRUE)

# # split new dataframe into grass and sediment
# Conewgrass <- Conew[c(1:3265, 6845:10109),]
# Conewsed <- Conew[3266:6844,]

#### 3.5. Co plot ####
# Customise theme
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

bCo <- ggplot() +
  geom_line(data = Conew, aes(Distance, fit, colour = Material), size = 0.5) +
  geom_ribbon(data = Conew, aes(Distance, ymin = lwr, ymax = upr, fill = Material),
              alpha = 0.5) +
  geom_pointrange(data = bCo.stat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                       colour = group1), size = 0.5) +
  scale_colour_manual(values = c("#d9b365", "#2b491e", "#81a512"),
                      labels = c("Sediment",
                                 expression(italic("Enhalus acoroides")),
                                 expression(italic("Thalassia hemprichii"))),
                      guide = guide_legend()) +
  scale_fill_manual(values = c("#d9b365", "#2b491e", "#81a512"),
                    labels = c("Sediment",
                               expression(italic("Enhalus acoroides")),
                               expression(italic("Thalassia hemprichii"))),
                    guide = guide_legend()) +
  annotate("text", x = c(60, 60, 60), y = c(0.96, 0.89, 0.82), size = 4.2,
           label = c("y == 0.53*e^{-0.06*x} + 0.31",
                     "y == 0.47*e^{-0.03*x} + 0.19",
                     "y == 0.37*e^{-0.02*x} + 0.24"),
           hjust = 0, parse = T) +
  ylab(expression("Cobalt concentration ("*mu*"g g"^-1*")")) +
  xlab("Distance from coast (m)") +
  coord_cartesian(ylim = c(0, 1), xlim = c(-10, 400)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 400, by = 50)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = c(.8, .9)) +
  mytheme
bCo # dimensions: 4 x 5 in


#### 4. Cr ####
require(ggplot2)
ggplot(metals, aes(Distance, Cr52, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) + 
  facet_grid(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

# remove outliers
metals$Cr52[c(2, 37, 45, 68:73, 77:81, 96, 161, 165, 166)] <- NA

ggplot(metals, aes(Distance, Cr52, colour = Direction)) +
  geom_point() +
  # geom_text(aes(label = ID)) + 
  facet_grid(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

#### 4.1. Clean Cr dataframe for analysis ####
grass <- metals[-c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:133, 158:168),]
grass <- grass[complete.cases(grass$Cr52),]
rownames(grass) <- NULL

sed <- metals[c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:133, 158:168),]
sed <- sed[complete.cases(sed$Cr52),]
rownames(sed) <- NULL

both <- metals
both <- both[complete.cases(both$Cr52),]
rownames(both) <- NULL

both$Material <- factor(both$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))

#### 4.2. Cr Rename variables ####
gm <- grass$Material
gd <- grass$Distance
gsite <- grass$Direction
gCr <- grass$Cr52

sm <- sed$Material
sd <- sed$Distance
ssite <- sed$Direction
sCr <- sed$Cr52

bm <- both$Material
bd <- both$Distance
bsite <- both$Direction
bCr<- both$Cr52

#### 4.3.  Cr Data analysis ####
# nonlinear least squares with self-start SSasymp function
# d = distance (metres)
# A = y asymptote that is approached by the curve
# C = the initial concentration at x = 0 (intercept) minus A
# k = the exponential decay rate per metre

?SSasymp
# self-starting function can only compute starting values for parameters
# when grass and sediment are treated separately because of large contrast

# Seagrass
gSm <- nls(gCr ~ SSasymp(gd, A, C, logk))
# non-linear model doesn't converge
plot(gCr ~ gd, pch = gm) # because data look linear

m5 <- lm(gCr ~ gd * gm)
drop1(m5, test = "F") # retain interaction (significant)

par(mfrow = c(2,2))
plot(m5)
par(mfrow = c(1,2))
plot(resid(m5) ~ gd)
boxplot(resid(m5) ~ gm)

qqnorm(resid(m5))
qqline(resid(m5))
hist(resid(m5))

require(car)
Anova(m5, type = 3)

summary(m5)
#                             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                0.7862129  0.0473773  16.595   <2e-16 ***
# gd                        -0.0008391  0.0002537  -3.308   0.0013 ** 
# gmThalassia hemprichii    -0.1558236  0.0663508  -2.348   0.0208 *  
# gd:gmThalassia hemprichii  0.0008239  0.0003569   2.308   0.0230 * 

# m5 is chosen as the optimal model for seagrass
# y = a*x + b

# Enhalus
# y = -0.0008391x + 0.7862129
# Thalassia
# y = 0.00002x + 0.6303893

gm <- factor(gm, levels = c("Thalassia hemprichii", "Enhalus acoroides"))
m5 <- lm(gCr ~ gd * gm)
Anova(m5, type = 3)
gm <- factor(gm, levels = c("Enhalus acoroides", "Thalassia hemprichii"))

# Sediment
sSm <- nls(sCr ~ SSasymp(sd, A, C, logk))
summary(sSm) # parameters sediment

m6 <- nls(sCr ~ C * exp(k * sd) + A, 
          start = list(C = 7.0838, 
                       k = -exp(-4.1034), 
                       A = 4.3859))
summary(m6)

# test model fit
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
plot(m6)
plot(resid(m6) ~ sd) 

par(mfrow = c(1,2))
hist(resid(m6))
qqnorm(resid(m6))
qqline(resid(m6)) # normality is fine
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

require(nlme)
m7 <- gnls(Cr52 ~ C * exp(k * Distance) + A, 
           start = list(C = c(0.74669), 
                        k = c(-0.10001), 
                        A = c(1.49375)),
           weights = varPower(),
           data = sed[,c(8, 15)])

summary(m7)

# test model fit
plot(m7)
plot(resid(m7, type = "normalized") ~ sd) # residual variance no longer differs with distance
# homogeneity is improved

par(mfrow = c(1,2))
hist(resid(m7))
qqnorm(resid(m7))
qqline(resid(m7)) # normality is unchanged 
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

# m7 is chosen as the optimal model for sediment
# interpret model
summary(m7)
# y = C * e^k*x + A 
# y = 2.681069 * e^(-0.016234)*x + 4.383120

#### 4.4. Cr descriptive statistics ####
require(psych)
bCr.stat <- describeBy(bCr, list(bd, bm), mat = T)
bCr.stat$group1 <- as.numeric(bCr.stat$group1) # convert to number

Crnewsed <- data.frame(Material = c(rep("Sediment", 3579)), 
                       Distance = c(seq(sd[1], sd[45], by = 0.1)))

Crnewgrass <- data.frame(gm = c(rep("Enhalus acoroides", 3265),
                                rep("Thalassia hemprichii", 3265)),
                         gd = c(seq(10.871910, 337.349640, by = 0.1),
                               seq(10.871910, 337.349640, by = 0.1)))

prediction <- data.frame(predict(m5, interval = "confidence", newdata = Crnewgrass))

Crnewgrass$fit <- prediction$fit
Crnewgrass$lwr <- prediction$lwr
Crnewgrass$upr <- prediction$upr

Crnewsed$fit <- predict(m7, newdata = Crnewsed)

# bootstrap confidence interval
bootfun <- function(newdata) {
  start <- coef(m7)
  boot <- sed[,c(8, 15)][sample(nrow(sed[,c(8, 15)]), 
                                    size = nrow(sed[,c(8,15)]), replace = TRUE),]
  bootfit <- try(update(m7,
                        start = start,
                        data = boot),
                 silent = TRUE)
  if (inherits(bootfit, "try-error")) return(rep(NA, nrow(newdata)))
  predict(bootfit, newdata)
}

bmat <- replicate(1000, bootfun(Crnewsed))
Crnewsed$lwr <- apply(bmat, 1, quantile, 0.025, na.rm = TRUE)
Crnewsed$upr <- apply(bmat, 1, quantile, 0.975, na.rm = TRUE)


colnames(Crnewgrass) <- c("Material", "Distance", "fit", "lwr", "upr")
Crnew <- rbind(Crnewsed, Crnewgrass)
Crnew$Material <- factor(Crnew$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))


#### 4.5. Cr plot ####
# Customise theme
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

bCr <- ggplot() +
  geom_line(data = Crnew, aes(Distance, fit, colour = Material, lty = Material), size = 0.5) +
  geom_ribbon(data = Crnew, aes(Distance, ymin = lwr, ymax = upr, fill = Material),
              alpha = 0.5) +
  geom_pointrange(data = bCr.stat, aes(group1, mean, ymin = mean - se, ymax = mean + se,
                                       colour = group2), size = 0.5) +
  scale_colour_manual(values = c("#d9b365", "#2b491e", "#81a512"),
                      labels = c("Sediment",
                                 expression(italic("Enhalus acoroides")),
                                 expression(italic("Thalassia hemprichii"))),
                      guide = guide_legend()) +
  scale_fill_manual(values = c("#d9b365", "#2b491e", "#81a512"),
                    labels = c("Sediment",
                               expression(italic("Enhalus acoroides")),
                               expression(italic("Thalassia hemprichii"))),
                    guide = guide_legend()) +
  scale_linetype_manual(values = c(1, 1, 5),
                        guide = "none") +
  annotate("text", x = c(60, 60, 60), y = c(9.6, 8.9, 8.2), size = 4.2,
           label = c("y == 2.68*e^{-0.02*x} + 4.38",
                     "y == -8%*%10^-4*x + 0.79",
                     "y == -2%*%10^-5*x + 0.63"),
           hjust = 0, parse = T) +
  ylab(expression("Chromium concentration ("*mu*"g g"^-1*")")) +
  xlab("Distance from coast (m)") +
  coord_cartesian(ylim = c(0, 10), xlim = c(-10, 400)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 400, by = 50)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = c(.8, .9)) +
  mytheme
bCr # dimensions: 4 x 5 in


#### 5. Ni ####
require(ggplot2)
ggplot(metals, aes(Distance, Ni60, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) + 
  facet_grid(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

# remove outliers
metals$Ni60[c(2, 31, 68:73, 77:81, 96, 148, 161)] <- NA

ggplot(metals, aes(Distance, Ni60, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) + 
  facet_grid(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

#### 5.1. Clean Ni dataframe for analysis ####
grass <- metals[-c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:133, 158:168),]
grass <- grass[complete.cases(grass$Ni60),]
rownames(grass) <- NULL

sed <- metals[c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:133, 158:168),]
sed <- sed[complete.cases(sed$Ni60),]
rownames(sed) <- NULL

both <- metals
both <- both[complete.cases(both$Ni60),]
rownames(both) <- NULL

both$Material <- factor(both$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))

#### 5.2. Ni Rename variables ####
gm <- grass$Material
gd <- grass$Distance
gsite <- grass$Direction
gNi <- grass$Ni60

sm <- sed$Material
sd <- sed$Distance
ssite <- sed$Direction
sNi <- sed$Ni60

bm <- factor(both$Material)
bd <- both$Distance
bsite <- both$Direction
bNi<- both$Ni60

#### 5.3.  Ni Data analysis ####
m8 <- lm(bNi ~ bm * bd)
drop1(m8, test = "F") # interaction not significant

par(mfrow = c(2,2))
plot(m8)
par(mfrow = c(1,2))
plot(resid(m8) ~ bd)
boxplot(resid(m8) ~ bm)

qqnorm(resid(m8))
qqline(resid(m8))
hist(resid(m8))

require(nlme)
m9 <- gls(bNi ~ bm * bd, weights = varIdent(form = ~1|bm),
          method = "ML")
drop1(m9, test = "Chisq") # interaction is significant
m9 <- gls(bNi ~ bm * bd, weights = varIdent(form = ~1|bm))

plot(m9)
par(mfrow = c(1,2))
plot(resid(m9, type = "normalized") ~ bd)
boxplot(resid(m9, type = "normalized") ~ bm)

qqnorm(resid(m9))
qqline(resid(m9))
hist(resid(m9))

require(car)
Anova(m9, type = 3)
# Response: bNi
#             Df     Chisq Pr(>Chisq)    
# (Intercept)  1 1349.9635  < 2.2e-16 ***
# bm           2  986.3801  < 2.2e-16 ***
# bd           1    2.9787   0.084365 .  
# bm:bd        2   13.0914   0.001436 ** 

# m9 is chosen as the optimal model 
summary(m9)
#                               Value  Std.Error   t-value p-value
# (Intercept)                8.797013 0.23942759  36.74185  0.0000
# bmEnhalus acoroides       -7.358226 0.24915410 -29.53283  0.0000
# bmThalassia hemprichii    -5.388229 0.26507760 -20.32699  0.0000
# bd                        -0.002567 0.00148737  -1.72590  0.0865
# bmEnhalus acoroides:bd     0.002302 0.00153385   1.50093  0.1355
# bmThalassia hemprichii:bd -0.000166 0.00160674  -0.10354  0.9177

# y = a*x + b
# Sediment
# y = -0.002567x + 8.797013
# Enhalus
# y = -0.000265x + 1.438787
# Thalassia
# y = -0.002733x + 3.408784

bm <- factor(bm, levels = c("Enhalus acoroides", "Sediment", "Thalassia hemprichii"))
m9 <- gls(bNi ~ bm * bd, weights = varIdent(form = ~1|bm))
Anova(m9, type = 3)
summary(m9)

bm <- factor(bm, levels = c("Thalassia hemprichii", "Enhalus acoroides", "Sediment"))
m9 <- gls(bNi ~ bm * bd, weights = varIdent(form = ~1|bm))
Anova(m9, type = 3)

bm <- factor(bm, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))
m9 <- gls(bNi ~ bm * bd, weights = varIdent(form = ~1|bm))

#### 5.4. Ni descriptive statistics ####
require(psych)
bNi.stat <- describeBy(bNi, list(bd, bm), mat = T) 
bNi.stat$group1 <- as.numeric(bNi.stat$group1) # convert to number

Ninew <- data.frame(bm = c(rep("Sediment", 3579),
                           rep("Enhalus acoroides", 3265),
                           rep("Thalassia hemprichii", 3265)),
                    bd = c(seq(0, 357.882667, by = 0.1),
                           seq(10.871910, 337.349640, by = 0.1),
                           seq(10.871910, 337.349640, by = 0.1)))

Ninew$bm <- factor(Ninew$bm, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))

Ninew$fit <- predict(m9, newdata = Ninew)
modmat <-  model.matrix(formula(m9)[-2], Ninew)
int <- diag(modmat %*% vcov(m9) %*% t(modmat))
Ninew$lwr <- with(Ninew, fit - qnorm(0.975)*sqrt(int))
Ninew$upr <- with(Ninew, fit + qnorm(0.975)*sqrt(int))

#### 5.5. Ni plot ####
# Customise theme
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

bNi <- ggplot() +
  geom_line(data = Ninew, aes(bd, fit, colour = bm, lty = bm), size = 0.5) +
  geom_ribbon(data = Ninew, aes(bd, ymin = lwr, ymax = upr, fill = bm),
              alpha = 0.5) +
  geom_pointrange(data = bNi.stat, aes(group1, mean, ymin = mean - se, ymax = mean + se,
                                       colour = group2), size = 0.5) +
  scale_colour_manual(values = c("#d9b365", "#2b491e", "#81a512"),
                      labels = c("Sediment",
                                 expression(italic("Enhalus acoroides")),
                                 expression(italic("Thalassia hemprichii"))),
                      guide = guide_legend()) +
  scale_fill_manual(values = c("#d9b365", "#2b491e", "#81a512"),
                    labels = c("Sediment",
                               expression(italic("Enhalus acoroides")),
                               expression(italic("Thalassia hemprichii"))),
                    guide = guide_legend()) +
  scale_linetype_manual(values = c(5, 5, 1),
                        guide = "none") +
  annotate("text", x = c(70, 70, 70), y = c(14.3, 13.2, 12.1), size = 4.2,
           label = c("y == -0.0026*x + 8.8",
                     "y == -3%*%10^-4*x + 1.44",
                     "y == -0.0027*x + 3.41"),
           hjust = 0, parse = T) +
  ylab(expression("Nickel concentration ("*mu*"g g"^-1*")")) +
  xlab("Distance from coast (m)") +
  coord_cartesian(ylim = c(0, 15), xlim = c(-10, 400)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 400, by = 50)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = c(.8, .9)) +
  mytheme
bNi # dimensions: 4 x 5 in


#### 6. Li #### 
require(ggplot2)
ggplot(metals, aes(Distance, Li7, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) +
  facet_wrap(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

# remove outliers
metals$Li7[c(23, 29, 31, 53, 67:73, 77:81)] <- NA

ggplot(metals, aes(Distance, Li7, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) +
  facet_wrap(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

#### 6.1. Clean Li dataframe ####
grass <- metals[-c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:133, 158:168),]
grass <- grass[complete.cases(grass$Li7),]
rownames(grass) <- NULL

sed <- metals[c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:133, 158:168),]
sed <- sed[complete.cases(sed$Li7),]
rownames(sed) <- NULL

both <- metals
both <- both[complete.cases(both$Li7),]
rownames(both) <- NULL

both$Material <- factor(both$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))

#### 6.2. Li Rename variables ####
gm <- grass$Material
gd <- grass$Distance
gsite <- grass$Direction
gLi <- grass$Li7

sm <- sed$Material
sd <- sed$Distance
ssite <- sed$Direction
sLi <- sed$Li7

bm <- factor(both$Material)
bd <- both$Distance
bsite <- both$Direction
bLi<- both$Li7

#### 6.3.  Li Data analysis ####
#Sediment
sSm <- nls(sNi ~ SSasymp(sd, A, C, logk))
# non-linear model doesn't converge
plot(sLi ~ sd, pch = sm) # because data look linear
m10 <- lm(sLi ~ sd)

par(mfrow = c(2,2))
plot(m10)
par(mfrow = c(1,1))
plot(resid(m10) ~ sd)

par(mfrow = c(1,2))
qqnorm(resid(m10))
qqline(resid(m10))
hist(resid(m10))

# m10 is chosen as the optimal model for sediment

require(car)
Anova(m10, type = 2)
# y = a*x + b
# y = -0.0004599x + 2.1530053

# Seagrass
# use nonlinear least squares with self-start SSasymp function
# d = distance (metres)
# A = y asymptote that is approached by the curve
# C = the initial concentration at x = 0 (intercept) minus A
# k = the exponential decay rate per metre
grass$Material <- factor(grass$Material, levels = c("Enhalus acoroides", "Thalassia hemprichii"))

?SSasymp
# self-starting function can only compute starting values for parameters
# when grass and sediment are treated separately because of large contrast
gSm <- nls(gLi ~ SSasymp(gd, A, C, logk))
summary(gSm) # parameters seagrass

m11 <- nls(Li7 ~ C[Material] * exp(k[Material] * Distance) + A[Material],
                start = list(C = c(0.19454, 0.19454),
                             k = c(-exp(-4.59044), -exp(-4.59044)),
                             A = c(0.08306, 0.08306)),
                data = grass[,c(6, 7, 15)])
summary(m11)
# test model fit
plot(m11)
plot(resid(m11) ~ gd) 

hist(resid(m11))
qqnorm(resid(m11))
qqline(resid(m11)) # normality is fine
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

require(nlme)
m12 <- gnls(Li7 ~ C * exp(k * Distance) + A, 
                  start = list(C = c(0.154470, 0.067054), 
                               k = c(-0.010593, -0.010976), 
                               A = c(0.094953, 0.073788)),
                  params = list(C ~ Material, k ~ Material, A ~ Material),
                  weights = varPower(),
                  data = grass[,c(6, 7, 15)])

summary(m12)
# test model fit
plot(m12)
plot(resid(m12, type = "normalized") ~ gd) # residual variance no longer differs with distance
# homogeneity is improved

par(mfrow = c(1,2))
hist(resid(m12))
qqnorm(resid(m12))
qqline(resid(m12)) # normality is unchanged 
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

# m12 is chosen as the optimal model for seagrass
# interpret model
summary(m12)
# # y = C * e^k*x +  
# Enhalus acoroides
# y = 0.15490634 * e^(-0.01117198)*x + 0.09642842 
# Thalassia hemprichii
# y = 0.06706047 * e^(-0.01092145)*x + 0.0737196

grass$Material <- factor(grass$Material, levels = c("Thalassia hemprichii", "Enhalus acoroides"))
m12 <- gnls(Li7 ~ C * exp(k * Distance) + A, 
            start = list(C = c(0.067054, 0.154470), 
                         k = c(-0.010976, -0.010593), 
                         A = c(0.073788, 0.094953)),
            params = list(C ~ Material, k ~ Material, A ~ Material),
            weights = varPower(),
            data = grass[,c(6, 7, 15)])
summary(m12)

grass$Material <- factor(grass$Material, levels = c("Enhalus acoroides", "Thalassia hemprichii"))
m12 <- gnls(Li7 ~ C * exp(k * Distance) + A, 
            start = list(C = c(0.154470, 0.067054), 
                         k = c(-0.010593, -0.010976), 
                         A = c(0.094953, 0.073788)),
            params = list(C ~ Material, k ~ Material, A ~ Material),
            weights = varPower(),
            data = grass[,c(6, 7, 15)])

#### 6.4. Li descriptive statistics ####
require(psych)
sLi.stat <- describeBy(sLi, sd, mat = T) 
sLi.stat$group1 <- as.numeric(sLi.stat$group1) # convert to number

gLi.stat <- describeBy(gLi, list(gm, gd), mat = T) 
gLi.stat$group2 <- as.numeric(gLi.stat$group2) # convert to number

# bLi.stat <- with(both, describeBy(Li7, list(Distance, Material), mat = T))
# bLi.stat$group1 <- as.numeric(bLi.stat$group1) # convert to number

Linewsed <- data.frame(sd = c(seq(sd[1], sd[48], by = 0.1)))

Linewgrass <- data.frame(gm = c(rep("Enhalus acoroides", 3265),
                                rep("Thalassia hemprichii", 3265)),
                         gd = c(seq(10.871910, 337.349640, by = 0.1),
                                seq(10.871910, 337.349640, by = 0.1)))

prediction <- data.frame(predict(m10, interval = "confidence", newdata = Linewsed))

Linewsed$fit <- prediction$fit
Linewsed$lwr <- prediction$lwr
Linewsed$upr <- prediction$upr

colnames(Linewgrass) <- c("Material", "Distance")
Linewgrass$fit <- predict(m12, newdata = Linewgrass)

# bootstrap confidence interval seagrass
bootfun <- function(newdata) {
  start <- coef(m12)
  boot <- grass[sample(nrow(grass), size = nrow(grass), replace = TRUE),]
  bootfit <- try(update(m12,
                        start = start,
                        data = boot),
                 silent = TRUE)
  if (inherits(bootfit, "try-error")) return(rep(NA, nrow(newdata)))
  predict(bootfit, newdata)
}

bmat <- replicate(1000, bootfun(Linewgrass))
Linewgrass$lwr <- apply(bmat, 1, quantile, 0.025, na.rm = TRUE)
Linewgrass$upr <- apply(bmat, 1, quantile, 0.975, na.rm = TRUE)

#### 6.5. Li plot ####
# Customise theme
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

gLi <- ggplot() +
  geom_line(data = Linewgrass, aes(Distance, fit, colour = Material, lty = Material), size = 0.5) +
  geom_ribbon(data = Linewgrass, aes(Distance, ymin = lwr, ymax = upr, fill = Material),
              alpha = 0.5) +
  geom_pointrange(data = gLi.stat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                       colour = group1), size = 0.5) +
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
  annotate("text", x = c(60, 60), y = c(0.385, 0.355), size = 4.2, family = "Helvetica Neue",
           hjust = 0, label = c("y == 0.15*e^{-0.01*x} + 0.1", "y == 0.07*e^{-0.01*x} + 0.07"),
           parse = T) +
  ylab(expression("Lithium concentration ("*mu*"g g"^-1*")")) +
  xlab("Distance from coast (m)") +
  coord_cartesian(ylim = c(0, 0.4), xlim = c(0, 350)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 350, by = 50)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.4, by = 0.10)) +
  theme(legend.position = c(.8, .93)) +
  mytheme
gLi # dimensions: 4 x 5 in

sLi <- ggplot() +
  geom_line(data = Linewsed, aes(sd, fit), size = 0.5, colour = "#d9b365",  lty = 5) +
  geom_ribbon(data = Linewsed, aes(sd, ymin = lwr, ymax = upr),
              alpha = 0.5, fill = "#d9b365") +
  geom_pointrange(data = sLi.stat, aes(group1, mean, ymin = mean - se, ymax = mean + se),
                  size = 0.5, colour = "#d9b365") +
  annotate("text", x = 240, y = 3.8, size = 4.2,
           hjust = 0, label = "y == -4.6%*%10^-4*x+ 2.15",
           parse = T) +
  ylab(expression("Lithium concentration ("*mu*"g g"^-1*")")) +
  xlab("Distance from coast (m)") +
  coord_cartesian(ylim = c(0, 4), xlim = c(-10, 400)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 400, by = 50)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 4, by = 1)) +
  theme(legend.position = c(.86, .93)) +
  mytheme
sLi # dimensions: 4 x 5 in


#### 7. Cu ####
require(ggplot2)
ggplot(metals, aes(Distance, Cu65, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) + 
  facet_grid(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

# remove outliers
metals$Cu65[c(1,7, 13, 23, 31, 45, 60, 68:73, 77:81, 89, 100, 161)] <- NA

ggplot(metals, aes(Distance, Cu65, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) + 
  facet_grid(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

#### 7.1. Clean Cu dataframe for analysis ####
grass <- metals[-c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:133, 158:168),]
grass <- grass[complete.cases(grass$Cu65),]
rownames(grass) <- NULL

sed <- metals[c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:133, 158:168),]
sed <- sed[complete.cases(sed$Cu65),]
rownames(sed) <- NULL

both <- metals
both <- both[complete.cases(both$Cu65),]
rownames(both) <- NULL

both$Material <- factor(both$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))

#### 7.2. Cu Rename variables ####
gm <- grass$Material
gd <- grass$Distance
gsite <- grass$Direction
gCu <- grass$Cu65

sm <- sed$Material
sd <- sed$Distance
ssite <- sed$Direction
sCu <- sed$Cu65

bm <- both$Material
bd <- both$Distance
bsite <- both$Direction
bCu<- both$Cu65

#### 7.3.  Cu Data analysis ####
# use nonlinear least squares with self-start SSasymp function
# d = distance (metres)
# A = y asymptote that is approached by the curve
# C = the initial concentration at x = 0 (intercept) minus A
# k = the exponential decay rate per metre

?SSasymp
# self-starting function can only compute starting values for parameters
# when grass and sediment are treated separately because of large contrast
gSm <- nls(gCu ~ SSasymp(gd, A, C, logk))
summary(gSm) # parameters seagrass

sSm <- nls(sCu ~ SSasymp(sd, A, C, logk))
summary(sSm) # parameters sediment

m13 <- nls(bCu ~ C[bm] * exp(k[bm] * bd) + A[bm], 
          start = list(C = c(16.0939, 11.3026, 11.3026), 
                       k = c(-exp(-3.5954), -exp(-3.8796), -exp(-3.8796)), 
                       A = c(0.4449, 2.5595, 2.5595)))
summary(m13)

# test model fit
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m13) ~ bd) # residual variance decreases dramatically with distance
boxplot(resid(m13) ~ bm) # and is much smaller for seagrass than sediment
# homogeneity needs improving

hist(resid(m13))
qqnorm(resid(m13))
qqline(resid(m13)) # normality is fine; deviance at the margins but balanced
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

# improve homogeneity with gnls
require(nlme)
m14 <- gnls(Cu65 ~ C * exp(k * Distance) + A, 
           start = list(C = c(15.64895, 9.26011, 8.35199), 
                        k = c(-0.02745, -0.02403, -0.01814), 
                        A = c(0.44494, 2.55549, 2.57294)),
           params = list(C ~ Material, k ~ Material, A ~ Material),
           weights = varPower(),
           data = both[,c(6, 11, 15)])

summary(m14)

# test model fit
plot(m14, col = bm)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m14, type = "normalized") ~ bd) # residual variance no longer differs as much with distance
boxplot(resid(m14, type = "normalized") ~ bm) # or material
# homogeneity is improved

hist(resid(m14))
qqnorm(resid(m14))
qqline(resid(m14)) # normality is unchanged 
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m14 is chosen as the optimal model

# interpret model
summary(m14)
# y = C * e^k*x + A  (always add the value of the intercept)
# Sediment
# y = 16.73606 * e^(-0.036119x) + 0.675797
# Enhalus
# y = 8.554612 * e^(-0.019372x) + 2.430187
# Thalassia
# y = 8.361331 *e^(-0.01781x) + 2.547336

# Relevel
both$Material <- factor(both$Material, levels = c("Enhalus acoroides", "Sediment", "Thalassia hemprichii"))
m14 <- gnls(Cu65 ~ C * exp(k * Distance) + A, 
            start = list(C = c(15.64895, 9.26011, 8.35199), 
                         k = c(-0.02745, -0.02403, -0.01814), 
                         A = c(0.44494, 2.55549, 2.57294)),
            params = list(C ~ Material, k ~ Material, A ~ Material),
            weights = varPower(),
            data = both[,c(6, 11, 15)])

summary(m14)

# Relevel
both$Material <- factor(both$Material, levels = c("Sediment", "Enhalus acoroides","Thalassia hemprichii"))
m14 <- gnls(Cu65 ~ C * exp(k * Distance) + A, 
            start = list(C = c(15.64895, 9.26011, 8.35199), 
                         k = c(-0.02745, -0.02403, -0.01814), 
                         A = c(0.44494, 2.55549, 2.57294)),
            params = list(C ~ Material, k ~ Material, A ~ Material),
            weights = varPower(),
            data = both[,c(6, 11, 15)])

summary(m14)

#### 7.4. Cu descriptive statistics ####
require(psych)
bCu.stat <- with(both, describeBy(Cu65, list(Distance, Material), mat = T))
bCu.stat$group1 <- as.numeric(bCu.stat$group1)

Cunew <- data.frame(Material = c(rep("Sediment", 3579),
                                 rep("Enhalus acoroides", 3265),
                                 rep("Thalassia hemprichii", 3265)),
                    Distance = c(seq(0, 357.882667, by = 0.1), #increment
                                 seq(10.871910, 337.349640, by = 0.1),
                                 seq(10.871910, 337.349640, by = 0.1)))

Cunew$Material <- factor(Cunew$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))
Cunew$fit <- predict(m14, newdata = Cunew)

# bootstrap confidence interval
bootfun <- function(newdata) {
  start <- coef(m14)
  boot <- both[,c(6, 11, 15)][sample(nrow(both[,c(6, 11, 15)]), 
                                     size = nrow(both[,c(6, 11, 15)]), replace = TRUE),]
  bootfit <- try(update(m14,
                        start = start,
                        data = boot),
                 silent = TRUE)
  if (inherits(bootfit, "try-error")) return(rep(NA, nrow(newdata)))
  predict(bootfit, newdata)
}

bmat <- replicate(1000, bootfun(Cunew))
Cunew$lwr <- apply(bmat, 1, quantile, 0.025, na.rm = TRUE)
Cunew$upr <- apply(bmat, 1, quantile, 0.975, na.rm = TRUE)

#### 7.5. Cu plot ####
# Customise theme
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

bCu <- ggplot() +
  geom_line(data = Cunew, aes(Distance, fit, colour = Material), size = 0.5) +
  geom_ribbon(data = Cunew, aes(Distance, ymin = lwr, ymax = upr, fill = Material),
              alpha = 0.5) +
  #geom_rug(data = grass, aes(x = Distance, y = Pb208), sides = "b") +
  geom_pointrange(data = bCu.stat, aes(group1, mean, ymin = mean - se, ymax = mean + se,
                                       colour = group2), size = 0.5) +
  scale_colour_manual(values = c("#d9b365", "#2b491e", "#81a512"),
                      labels = c("Sediment",
                                 expression(italic("Enhalus acoroides")),
                                 expression(italic("Thalassia hemprichii"))),
                      guide = guide_legend()) +
  scale_fill_manual(values = c("#d9b365", "#2b491e", "#81a512"),
                    labels = c("Sediment",
                               expression(italic("Enhalus acoroides")),
                               expression(italic("Thalassia hemprichii"))),
                    guide = guide_legend()) +
  annotate("text", x = c(60, 60, 60), y = c(28.8, 26.65, 24.5), size = 4.2,
           label = c("y == 16.74*e^{-0.04*x} + 0.68",
                     "y == 8.55*e^{-0.02*x} + 2.43",
                     "y == 8.36*e^{-0.02*x} + 2.55"),
           hjust = 0, parse = T) +
  ylab(expression("Copper concentration ("*mu*"g g"^-1*")")) +
  xlab("Distance from coast (m)") +
  coord_cartesian(ylim = c(0, 30), xlim = c(-10, 400)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 400, by = 50)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = c(.8, .9)) +
  mytheme
bCu # dimensions: 4 x 5 in

#### 8. Zn #### 
require(ggplot2)
ggplot(metals, aes(Distance, Zn66, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) + 
  facet_grid(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

# remove outliers
metals$Zn66[c(23, 31, 37, 65, 68:73, 77:81, 110, 118, 123, 160, 161)] <- NA

ggplot(metals, aes(Distance, Zn66, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) + 
  facet_grid(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

#### 8.1. Clean Zn dataframe for analysis ####
grass <- metals[-c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:133, 158:168),]
grass <- grass[complete.cases(grass$Zn66),]
rownames(grass) <- NULL

sed <- metals[c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:133, 158:168),]
sed <- sed[complete.cases(sed$Zn66),]
rownames(sed) <- NULL

both <- metals
both <- both[complete.cases(both$Zn66),]
rownames(both) <- NULL

both$Material <- factor(both$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))

#### 8.2. Zn Rename variables ####
gm <- grass$Material
gd <- grass$Distance
gsite <- grass$Direction
gZn <- grass$Zn66

sm <- sed$Material
sd <- sed$Distance
ssite <- sed$Direction
sZn <- sed$Zn66

bm <- both$Material
bd <- both$Distance
bsite <- both$Direction
bZn<- both$Zn66

#### 8.3.  Zn Data analysis ####
?SSasymp
# self-starting function can only compute starting values for parameters
# when grass and sediment are treated separately because of large contrast
gSm <- nls(gZn ~ SSasymp(gd, A, C, logk))
summary(gSm) # parameters seagrass

sSm <- nls(sZn ~ SSasymp(sd, A, C, logk))
summary(sSm) # parameters sediment

m15 <- nls(bZn ~ C[bm] * exp(k[bm] * bd) + A[bm], 
            start = list(C = c(64.9920, 77.2440, 77.2440), 
                         k = c(-exp(-3.1014), -exp(-2.6385), -exp(-2.6385)), 
                         A = c(0.9800, 39.1147, 39.1147)))
summary(m15)

# test model fit
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m15) ~ bd) 
boxplot(resid(m15) ~ bm) # residual variance smaller for seagrass than sediment
# homogeneity needs improving

hist(resid(m15))
qqnorm(resid(m15))
qqline(resid(m15))
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

# improve homogeneity with gnls
require(nlme)
m16 <- gnls(Zn66 ~ C * exp(k * Distance) + A, 
             start = list(C = c(64.011966, 58.631099, 10.866885), 
                          k = c(-0.044984, -0.070794, -0.008098), 
                          A = c(0.979981, 40.044297, 35.057422)),
             params = list(C ~ Material, k ~ Material, A ~ Material),
             weights = varPower(),
             data = both[,c(6, 12, 15)])

summary(m16)
# test model fit
plot(m16, col = bm)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m16, type = "normalized") ~ bd) # residual variance more balanced
boxplot(resid(m16, type = "normalized") ~ bm) # material's residual variance too
# homogeneity is improved

hist(resid(m16))
qqnorm(resid(m16))
qqline(resid(m16)) # normality is unchanged 
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m16 is chosen as the optimal model

# interpret model
summary(m16)
# y = C * e^k*x + A  (always add the value of the intercept-Enhalus)
# Sediment
# y = 66.78607 * e^(-0.05139x) + 1.45096
# Enhalus
# y = 58.37067 * e^(-0.07049x) + 40.04181
# Thalassia
# y = 10.87252 *e^(-0.00904x) + 35.31414

# Relevel
both$Material <- factor(both$Material, levels = c("Enhalus acoroides", "Sediment", "Thalassia hemprichii"))
m16 <- gnls(Zn66 ~ C * exp(k * Distance) + A, 
            start = list(C = c(64.011966, 58.631099, 10.866885), 
                         k = c(-0.044984, -0.070794, -0.008098), 
                         A = c(0.979981, 40.044297, 35.057422)),
            params = list(C ~ Material, k ~ Material, A ~ Material),
            weights = varPower(),
            data = both[,c(6, 12, 15)])

summary(m16)

# Relevel
both$Material <- factor(both$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))
m16 <- gnls(Zn66 ~ C * exp(k * Distance) + A, 
            start = list(C = c(64.011966, 58.631099, 10.866885), 
                         k = c(-0.044984, -0.070794, -0.008098), 
                         A = c(0.979981, 40.044297, 35.057422)),
            params = list(C ~ Material, k ~ Material, A ~ Material),
            weights = varPower(),
            data = both[,c(6, 12, 15)])

summary(m16)

#### 8.4. Zn Descriptive statistics ####
require(psych)
bZn.stat <- with(both, describeBy(Zn66, list(Distance, Material), mat = T))
bZn.stat$group1 <- as.numeric(bZn.stat$group1)

Znnew <- data.frame(Material = c(rep("Sediment", 3579),
                                 rep("Enhalus acoroides", 3265),
                                 rep("Thalassia hemprichii", 3265)),
                    Distance = c(seq(0, 357.882667, by = 0.1), #increment
                                 seq(10.871910, 337.349640, by = 0.1),
                                 seq(10.871910, 337.349640, by = 0.1)))

Znnew$Material <- factor(Znnew$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))
Znnew$fit <- predict(m16, newdata = Znnew)

# bootstrap confidence interval
bootfun <- function(newdata) {
  start <- coef(m16)
  boot <- both[,c(6, 12, 15)][sample(nrow(both[,c(6, 12, 15)]), 
                                     size = nrow(both[,c(6, 12, 15)]), replace = TRUE),]
  bootfit <- try(update(m16,
                        start = start,
                        data = boot),
                 silent = TRUE)
  if (inherits(bootfit, "try-error")) return(rep(NA, nrow(newdata)))
  predict(bootfit, newdata)
}

bmat <- replicate(1000, bootfun(Znnew))
Znnew$lwr <- apply(bmat, 1, quantile, 0.025, na.rm = TRUE)
Znnew$upr <- apply(bmat, 1, quantile, 0.975, na.rm = TRUE)

#### 8.5. Zn plot ####
# Customise theme
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

bZn <- ggplot() +
  geom_line(data = Znnew, aes(Distance, fit, colour = Material, lty = Material), size = 0.5) +
  geom_ribbon(data = Znnew, aes(Distance, ymin = lwr, ymax = upr, fill = Material),
              alpha = 0.5) +
  #geom_rug(data = grass, aes(x = Distance, y = Pb208), sides = "b") +
  geom_pointrange(data = bZn.stat, aes(group1, mean, ymin = mean - se, ymax = mean + se,
                                       colour = group2), size = 0.5) +
  scale_colour_manual(values = c("#d9b365", "#2b491e", "#81a512"),
                      labels = c("Sediment",
                                 expression(italic("Enhalus acoroides")),
                                 expression(italic("Thalassia hemprichii"))),
                      guide = guide_legend()) +
  scale_fill_manual(values = c("#d9b365", "#2b491e", "#81a512"),
                    labels = c("Sediment",
                               expression(italic("Enhalus acoroides")),
                               expression(italic("Thalassia hemprichii"))),
                    guide = guide_legend()) +
  scale_linetype_manual(values = c(1, 5, 5),
                        guide = "none") +
  annotate("text", x = c(60, 60, 60), y = c(76.7, 71.2, 65.7), size = 4.2,
           label = c("y == 66.79*e^{-0.05*x} + 1.45",
                     "y == 58.37*e^{-0.07*x} + 40.04",
                     "y == 10.87*e^{-0.01*x} + 35.31"),
           hjust = 0, parse = T) +
  ylab(expression("Zinc concentration ("*mu*"g g"^-1*")")) +
  xlab("Distance from coast (m)") +
  coord_cartesian(ylim = c(0, 80), xlim = c(-10, 400)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 400, by = 50)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 80, by = 10)) +
  theme(legend.position = c(.8, .9)) +
  mytheme
bZn # dimensions: 4 x 5 in

#### 9. Cd ####
# ignore
ggplot(metals, aes(Distance, Cd111, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) + 
  facet_grid(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

# remove outliers
metals$Cd111[c(28, 47)] <- NA

ggplot(metals, aes(Distance, Cd111, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) + 
  facet_grid(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")


#### 10.  Clean up ####
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

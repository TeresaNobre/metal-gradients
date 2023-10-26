


#### 1.  Data preparation ####
#### 1.1 Load data ####
metals <- read.csv("~/PATH/CombinedDataframe.csv")

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

metals$Distance <- c(distN, distE, distS, distW)
# add combined vector to dataframe

require(ggplot2)
ggplot(metals, aes(Distance, MPI, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) + # this line to see IDs
  facet_wrap(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

# remove outliers
metals$MPI[c(2, 23, 28, 29, 31, 45, 47, 53, 59, 60, 62, 63, 96, 67:73, 77:81, 161)] <- NA

ggplot(metals, aes(Distance, MPI, colour = Direction)) +
  geom_point() +
  geom_text(aes(label = ID)) + # this line to see IDs
  facet_wrap(~ Material, scales = "free") +
  theme_minimal() +
  theme(legend.position = "top")

#### 2.1. Clean  dataframe for analysis ####
grass <- metals[-c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:133, 158:168),]
grass <- grass[complete.cases(grass$MPI),]
rownames(grass) <- NULL

sed <- metals[c(1:6, 19:21, 34:36, 49:51, 64:67, 74:76, 82:90, 103:105, 130:133, 158:168),]
sed <- sed[complete.cases(sed$MPI),]
rownames(sed) <- NULL

both <- metals
both <- both[complete.cases(both$MPI),]
rownames(both) <- NULL

both$Material <- factor(both$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))

gm <- grass$Material
gd <- grass$Distance
gsite <- grass$Direction
gMPI <- grass$MPI

sm <- sed$Material
sd <- sed$Distance
ssite <- sed$Direction
sMPI <- sed$MPI

bm <- factor(both$Material)
bd <- both$Distance
bsite <- both$Direction
bMPI<- both$MPI

?SSasymp
# self-starting function can only compute starting values for parameters
# when grass and sediment are treated separately because of large contrast
gSm <- nls(gMPI ~ SSasymp(gd, A, C, logk))
summary(gSm) # parameters seagrass

sSm <- nls(sMPI ~ SSasymp(sd, A, C, logk))
summary(sSm) # parameters sediment

m1 <- nls(bMPI ~ C[bm] * exp(k[bm] * bd) + A[bm], 
          start = list(C = c(10.6040, 4.84620, 4.84620), 
                       k = c(-exp(-3.2756), -exp(-3.88333),-exp(-3.88333)), 
                       A = c(2.1360, 1.89305, 1.89305)))
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

require(nlme)
m2 <- gnls(MPI ~ C * exp(k * Distance) + A, 
           start = list(C = c(10.6040, 4.84620, 4.84620), 
                        k = c(-exp(-3.2756), -exp(-3.88333),-exp(-3.88333)), 
                        A = c(2.1360, 1.89305, 1.89305)),
           params = list(C ~ Material, k ~ Material, A ~ Material),
           weights = varPower(),
           data = both[,c(6, 19, 20)])

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

summary(m2)
# y = C * e^k*x + A  (always add the value of the intercept-Sediment)
# Sediment MPI
# y = 8.115248 * e^(-0.033055x) + 2.088682
# Enhalus MPI
# y = 3.110447 * e^(-0.020784x) + 1.90254
# Thalassia MPI
# y = 2.576529 *e^(-0.013271x) + 1.74666

# Relevel 
both$Material <- factor(both$Material, levels = c("Enhalus acoroides", "Sediment", "Thalassia hemprichii"))
m2 <- gnls(MPI ~ C * exp(k * Distance) + A, 
           start = list(C = c(4.84620, 10.6040, 4.84620), 
                        k = c(-exp(-3.88333), -exp(-3.2756),-exp(-3.88333)), 
                        A = c(1.89305, 2.1360, 1.89305)),
           params = list(C ~ Material, k ~ Material, A ~ Material),
           weights = varPower(),
           data = both[,c(6, 19, 20)])
summary(m2)

# Relevel 
both$Material <- factor(both$Material, levels = c("Thalassia hemprichii", "Sediment", "Enhalus acoroides"))
m2 <- gnls(MPI ~ C * exp(k * Distance) + A, 
           start = list(C = c(4.84620, 10.6040, 4.84620), 
                        k = c(-exp(-3.88333), -exp(-3.2756),-exp(-3.88333)), 
                        A = c(1.89305, 2.1360, 1.89305)),
           params = list(C ~ Material, k ~ Material, A ~ Material),
           weights = varPower(),
           data = both[,c(6, 19, 20)])
summary(m2)

# Relevel 
both$Material <- factor(both$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))
m2 <- gnls(MPI ~ C * exp(k * Distance) + A, 
           start = list(C = c(10.6040, 4.84620, 4.84620), 
                        k = c(-exp(-3.2756), -exp(-3.88333),-exp(-3.88333)), 
                        A = c(2.1360, 1.89305, 1.89305)),
           params = list(C ~ Material, k ~ Material, A ~ Material),
           weights = varPower(),
           data = both[,c(6, 19, 20)])
summary(m2)

require(psych)
bMPI.stat <- with(both, describeBy(bMPI, list(Material, Distance), mat = T))
bMPI.stat$group2 <- as.numeric(bMPI.stat$group2)

MPInew <- data.frame(Material = c(rep("Sediment", 3579),
                                 rep("Enhalus acoroides", 3265),
                                 rep("Thalassia hemprichii", 3265)),
                    Distance = c(seq(0, 357.882667, by = 0.1),
                                 seq(10.87191, 337.34964, by = 0.1),
                                 seq(10.87191, 337.34964, by = 0.1)))

MPInew$Material <- factor(MPInew$Material, levels = c("Sediment", "Enhalus acoroides", "Thalassia hemprichii"))
MPInew$fit <- predict(m2, newdata = MPInew)

# bootstrap confidence interval
bootfun <- function(newdata) {
  start <- coef(m2)
  boot <- both[,c(6, 19, 20)][sample(nrow(both[,c(6, 19, 20)]), 
                                     size = nrow(both[,c(6, 19, 20)]), replace = TRUE),]
  bootfit <- try(update(m2,
                        start = start,
                        data = boot),
                 silent = TRUE)
  if (inherits(bootfit, "try-error")) return(rep(NA, nrow(newdata)))
  predict(bootfit, newdata)
}

bmat <- replicate(1000, bootfun(MPInew))
MPInew$lwr <- apply(bmat, 1, quantile, 0.025, na.rm = TRUE)
MPInew$upr <- apply(bmat, 1, quantile, 0.975, na.rm = TRUE)

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
bMPI <- ggplot() +
  geom_line(data = MPInew, aes(Distance, fit, colour = Material), size = 0.5) +
  geom_ribbon(data = MPInew, aes(Distance, ymin = lwr, ymax = upr, fill = Material),
              alpha = 0.5) +
  geom_pointrange(data = bMPI.stat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
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
  annotate("text", x = c(70, 70, 70), y = c(11.55, 10.65, 9.75), size = 4.2,
           label = c("y == 8.12*e^{-0.03*x} + 2.1",
                     "y == 3.11*e^{-0.02*x} + 1.9",
                     "y == 2.58*e^{-0.01*x} + 1.7"),
           hjust = 0, parse = T) +
  ylab(expression("Metal Pollution Index (MPI)")) +
  xlab("Distance from coast (m)") +
  coord_cartesian(ylim = c(0, 12), xlim = c(-10, 400)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 400, by = 50)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = c(.8, .9)) +
  mytheme
bMPI # dimensions: 4 x 5 in


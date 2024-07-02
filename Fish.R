#### 1.  Data preparation ####
#### 1.1 Load data ####
fish <- read.csv("~/PATH/Fish.csv")

MPI <- fish$MPI
Al <- fish$Al
Fe <- fish$Fe
Zn <- fish$Zn66
Cu <- fish$Cu65[-c(2)]
Cr <- fish$Cr52[-c(1,9)]
Cr1 <- fish$Cr52
Co <- fish$Co59[-c(6)]
Li <- fish$Li7
Ni <- fish$Ni60
Pb <- fish$Pb208[-c(1)]
Cd <- fish$Cd111[-c(6)]

sp <- fish$Material
spPb <- fish$Material[-c(1)]
spCo <- fish$Material[-c(6)]
spCd <- fish$Material[-c(6)]
spCu <- fish$Material[-c(2)]
spCr <- fish$Material[-c(1,9)]

# MPI ####
require(psych)
with(fish, describeBy(MPI, Material , mat = T))
m0 <- lm(MPI~ sp)
par(mfrow = c(2,2))
plot(m0)
par(mfrow = c(1,1))
boxplot(resid(m0) ~ sp)

par(mfrow = c(1,2))
qqnorm(resid(m0))
qqline(resid(m0))
hist(resid(m0))

require(nlme)

m01 <- gls(MPI ~ sp, weights = varIdent(form = ~1|sp))
plot(m01)
boxplot(resid(m01, type = "normalized") ~ sp)
m01

par(mfrow = c(1,2))
qqnorm(resid(m01))
qqline(resid(m01))
hist(resid(m01))

require(car)
Anova(m01, type = 2)
boxplot(MPI ~ sp)

#### Al ####
require(psych)
with(fish, describeBy(Al, Material , mat = T))
m1 <- lm(Al~ sp)
par(mfrow = c(2,2))
plot(m1)
par(mfrow = c(1,1))
boxplot(resid(m1) ~ sp)

par(mfrow = c(1,2))
qqnorm(resid(m1))
qqline(resid(m1))
hist(resid(m1))

m2 <- gls(Al~ sp, weights = varIdent(form = ~1|sp))
plot(m2)
boxplot(resid(m2, type = "normalized") ~ sp)

par(mfrow = c(1,2))
qqnorm(resid(m2))
qqline(resid(m2))
hist(resid(m2))

require(car)
Anova(m2, type = 2)
boxplot(Al ~ sp)

summary(m2)

sp <- factor(sp, levels = c("Siganus canaliculatus", "Crenimugil buchanani"))
m2 <- gls(Al~ sp, weights = varIdent(form = ~1|sp))
summary(m2)

# Fe ####
with(fish, describeBy(Fe, Material , mat = T))
m3 <- lm(Fe~ sp)
par(mfrow = c(2,2))
plot(m3)
par(mfrow = c(1,1))
boxplot(resid(m3) ~ sp)

par(mfrow = c(1,2))
qqnorm(resid(m1))
qqline(resid(m1))
hist(resid(m1))

Anova(m3, type = 2)

summary(m3)
sp <- factor(sp, levels = c("Crenimugil buchanani", "Siganus canaliculatus"))
m3 <- lm(Fe~ sp)
summary(m3)

# Zn ####
with(fish, describeBy(Zn66, Material , mat = T))
m4 <- lm(Zn ~ sp)
par(mfrow = c(2,2))
plot(m4)
par(mfrow = c(1,1))
boxplot(resid(m4) ~ sp)

par(mfrow = c(1,2))
qqnorm(resid(m1))
qqline(resid(m1))
hist(resid(m1))

m5 <- gls(Zn ~ sp, weights = varIdent(form = ~1|sp))
plot(m5)
boxplot(resid(m5, type = "normalized") ~ sp)

par(mfrow = c(1,2))
qqnorm(resid(m5))
qqline(resid(m5))
hist(resid(m5))

Anova(m5, type = 2)

summary(m5)
sp <- factor(sp, levels = c("Siganus canaliculatus", "Crenimugil buchanani"))
m5 <- gls(Zn ~ sp, weights = varIdent(form = ~1|sp))
summary(m5)

#### Cu ####
with(fish, describeBy(Cu, spCu , mat = T))
m6 <- lm(Cu ~ spCu)
par(mfrow = c(2,2))
plot(m6)
par(mfrow = c(1,1))
boxplot(resid(m6) ~ spCu)

par(mfrow = c(1,2))
qqnorm(resid(m6))
qqline(resid(m6))
hist(resid(m6))

m7 <- gls(Cu ~ spCu, weights = varIdent(form = ~1|spCu))
plot(m7)
boxplot(resid(m7, type = "normalized") ~ spCu)

par(mfrow = c(1,2))
qqnorm(resid(m7))
qqline(resid(m7))
hist(resid(m7))

Anova(m7, type = 2)

summary(m7)
spCu <- factor(spCu, levels = c("Siganus canaliculatus", "Crenimugil buchanani"))
m7 <- gls(Cu ~ spCu, weights = varIdent(form = ~1|spCu))
summary(m7)

# Cr ####
with(fish, describeBy(Cr, spCr , mat = T))
m8 <- lm(Cr ~ spCr)
par(mfrow = c(2,2))
plot(m8)
par(mfrow = c(1,1))
boxplot(resid(m8) ~ spCr)

par(mfrow = c(1,2))
qqnorm(resid(m8))
qqline(resid(m8))
hist(resid(m8))

m9 <- gls(Cr ~ spCr, weights = varIdent(form = ~1|spCr))
plot(m9)
boxplot(resid(m9, type = "normalized") ~ spCr)

par(mfrow = c(1,2))
qqnorm(resid(m9))
qqline(resid(m9))
hist(resid(m9))

Anova(m9, type = 2)

summary(m9)
spCr <- factor(spCr, levels = c("Siganus canaliculatus", "Crenimugil buchanani"))
m9 <- gls(Cr ~ spCr, weights = varIdent(form = ~1|spCr))
summary(m9)

# Pb ####
with(fish, describeBy(Pb, spPb , mat = T))
m10 <- lm(Pb ~ spPb)
par(mfrow = c(2,2))
plot(m10)
par(mfrow = c(1,1))
boxplot(resid(m10) ~ spPb)

par(mfrow = c(1,2))
qqnorm(resid(m10))
qqline(resid(m10))
hist(resid(m10))

m11 <- glm(Pb ~ spPb, family = Gamma(link = "log"))
par(mfrow = c(2,2))
plot(m11)
par(mfrow = c(1,1))
boxplot(resid(m11) ~ spPb)

par(mfrow = c(1,2))
qqnorm(resid(m11))
qqline(resid(m11))
hist(resid(m11))

Anova(m11, type = 2)

summary(m11)
coef(m11)
spPb <- factor(spPb, levels = c("Siganus canaliculatus", "Crenimugil buchanani"))
m11 <- glm(Pb ~ spPb, family = Gamma(link = "log"))
summary(m11)
coef(m11)

# Ni ####
with(fish, describeBy(Ni, Material , mat = T))
m12 <- lm(Ni ~ sp)
par(mfrow = c(2,2))
plot(m12)
par(mfrow = c(1,1))
boxplot(resid(m12) ~ sp)

par(mfrow = c(1,2))
qqnorm(resid(m12))
qqline(resid(m12))
hist(resid(m12))

m13 <- gls(Ni ~ sp, weights = varIdent(form = ~1|sp))
plot(m13)
boxplot(resid(m13, type = "normalized") ~ sp)

par(mfrow = c(1,2))
qqnorm(resid(m13))
qqline(resid(m13))
hist(resid(m13))

Anova(m13, type = 2)

summary(m13)
sp <- factor(sp, levels = c("Crenimugil buchanani", "Siganus canaliculatus"))
m13 <- gls(Ni ~ sp, weights = varIdent(form = ~1|sp))
summary(m13)

# Li ####
with(fish, describeBy(Li, Material , mat = T))
m14 <- lm(Li ~ sp)
par(mfrow = c(2,2))
plot(m14)
par(mfrow = c(1,1))
boxplot(resid(m14) ~ sp)

par(mfrow = c(1,2))
qqnorm(resid(m14))
qqline(resid(m14))
hist(resid(m14))

Anova(m14, type = 2)

summary(m14)
sp <- factor(sp, levels = c("Siganus canaliculatus", "Crenimugil buchanani"))
m14 <- lm(Li ~ sp)
summary(m14)

# Co ####
with(fish, describeBy(Co, spCo , mat = T))
m15 <- lm(Co ~ spCo)
par(mfrow = c(2,2))
plot(m15)
par(mfrow = c(1,1))
boxplot(resid(m15) ~ spCo)

par(mfrow = c(1,2))
qqnorm(resid(m15))
qqline(resid(m15))
hist(resid(m15))

m16 <- gls(Co ~ spCo, weights = varIdent(form = ~1|spCo))
plot(m16)
boxplot(resid(m16, type = "normalized") ~ spCo)

par(mfrow = c(1,2))
qqnorm(resid(m16))
qqline(resid(m16))
hist(resid(m16))

Anova(m16, type = 2)

summary(m16)
spCo <- factor(spCo, levels = c("Siganus canaliculatus", "Crenimugil buchanani"))
m16 <- gls(Co ~ spCo, weights = varIdent(form = ~1|spCo))
summary(m16)

# Cd ####
with(fish, describeBy(Cd, spCd , mat = T))
m17 <- lm(Cd ~ spCd)
par(mfrow = c(2,2))
plot(m17)
par(mfrow = c(1,1))
boxplot(resid(m17) ~ spCd)

par(mfrow = c(1,2))
qqnorm(resid(m17))
qqline(resid(m17))
hist(resid(m17))

m18 <- gls(Cd ~ spCd, weights = varIdent(form = ~1|spCd))
plot(m18)
boxplot(resid(m18, type = "normalized") ~ spCd)

par(mfrow = c(1,2))
qqnorm(resid(m18))
qqline(resid(m18))
hist(resid(m18))

Anova(m18, type = 2)

summary(m18)
spCd <- factor(spCd, levels = c("Siganus canaliculatus", "Crenimugil buchanani"))
m18 <- gls(Cd ~ spCd, weights = varIdent(form = ~1|spCd))
summary(m18)

# MPI ####
m19 <- lm(MPI ~ sp)
par(mfrow = c(2,2))
plot(m19)
par(mfrow = c(1,1))
boxplot(resid(m19) ~ sp)

par(mfrow = c(1,2))
qqnorm(resid(m19))
qqline(resid(m19))
hist(resid(m19))

m20 <- gls(MPI ~ sp, weights = varIdent(form = ~1|sp))
plot(m20)
boxplot(resid(m20, type = "normalized") ~ sp)

par(mfrow = c(1,2))
qqnorm(resid(m20))
qqline(resid(m20))
hist(resid(m20))

Anova(m20)


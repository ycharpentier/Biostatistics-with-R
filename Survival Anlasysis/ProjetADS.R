library(survival)
library(corrplot)
library(MASS)
library("dplyr")

rm(list=ls())
setwd("C:/Users/yoyoc/Desktop/5A/ADS")
getwd()

# Ouverture du df
heart_df <- read.csv("heart_failure_clinical_records_dataset.csv")

# Visualisation du df
#View(heart_df)
head(heart_df)

# Quelques valeurs de base
names(heart_df) # nom des variables
str(heart_df)
summary(heart_df)
#attach(heart_df)

# Prétraitement
heart_df$sex = as.factor(heart_df$sex)
heart_df$high_blood_pressure = as.factor(heart_df$high_blood_pressure)
heart_df$anaemia = as.factor(heart_df$anaemia)
heart_df$diabetes = as.factor(heart_df$diabetes)
heart_df$smoking = as.factor(heart_df$smoking)

str(heart_df)

heart_df_num = select_if(heart_df, is.numeric)
heart_df_covar = heart_df[,1:11]

# Corrélations éventuelles ?
corrplot(cor(heart_df_num))

# Premier Kaplan Meier
Y = Surv(heart_df$time, heart_df$DEATH_EVENT)
summary(survfit(Y~1,conf.type="plain"))
plot(survfit(Y~1, conf.type="plain"), col="blue", lwd=2,main='Courbe de Kaplan')


# Second Kaplan Meier (avec filtre)
t1 <- heart_df$time[heart_df$anaemia==1]
e1 <- heart_df$DEATH_EVENT[heart_df$anaemia==1]
Y1 <- Surv(t1,e1)

t2 <- heart_df$time[heart_df$anaemia==0]
e2 <- heart_df$DEATH_EVENT[heart_df$anaemia==0]
Y2 <- Surv(t2,e2)

plot(survfit(Y1~1), col="blue", main="Estimation Kaplan de la fonction de survie des patients\n atteint d'une insuffisance cardiaque des fumeur et non fumeurs.")
lines(survfit(Y2~1),col="red")
legend(0,0.2,legend=c("anaemia", "non anemia"), col=c("blue", "red"),
       lty=c(1,1), cex = 0.7)

# On veut tester les interractions éventuelles
M1 <- coxph(Y~., data = heart_df_covar ,method = "breslow")
summary(M1)
test <- stepAIC(M1,trace=TRUE, direction=c("backward"))






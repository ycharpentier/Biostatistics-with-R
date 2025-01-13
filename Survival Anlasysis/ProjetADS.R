library(survival)
library(corrplot)
library(MASS)
library(latex2exp)
library("dplyr")
library("fitdistrplus")

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

# Quantité de données censurées 
hist(heart_df$DEATH_EVENT)

# Premier Kaplan Meier avec censure 
Y = Surv(heart_df$time, heart_df$DEATH_EVENT)
fit1 <- survfit(Y~1,conf.type="plain", type=c("kaplan-meier"))
summary(fit1)
plot(fit1, col="blue", lwd=2, main='Courbe de Kaplan')

# Modèle sans censure
Y_no_cens <- Surv(heart_df$time, rep(1, nrow(heart_df)))  # Ici, on considère que tous les événements sont survenus
fit_no_cens <- survfit(Y_no_cens ~ 1, conf.type="plain", type="kaplan-meier")
summary(fit_no_cens)
plot(fit_no_cens, col="red", lwd=2, main='Courbe de Kaplan-Meier sans censure')

# Comparaison des courbes avec et sans censure
plot(fit1, col="blue", lwd=2, main='Comparison of Kaplan-Meier Curves for Heart Failure Data: Censored vs. Non-Censored',
     pch=11,
     xlab="Event time (death)",
     ylab="Survival")

lines(fit_no_cens, col="red", lwd=2)
grid(nx = NULL, ny = NA, lty = 3, col = "gray", lwd = 1)
legend("topright", legend=c("With censorship", "Without censorship"), col=c("blue", "red"), lwd=2)

#La courbe avec censure donne une estimation plus juste et plus réaliste des probabilités de survie, car elle considère que certains patients n'ont pas encore subi l'événement (décès) à la fin de la période d'observation. Elle montre la probabilité de survie réelle, en tenant compte des patients toujours en vie mais non suivis jusqu'au bout.
#La courbe sans censure représente un scénario extrême où tous les patients finissent par subir l'événement, ce qui est généralement irréaliste.

# La probabilité de survivre au delà du 50ème jour est au moins de 80%.













H <- -log(fit1$surv) # Estimation de la fonction de risque cumulé avec Breslow

plot(H, pch=2, col='red',
     main = TeX(r'($\hat{H}_{Breslow} = - \log(\hat{S}_{KM})$)'),
     xlab = TeX(r'($t$)'),  
     ylab = TeX(r'($\hat{H}(t)$)')  
)

# Estimation de la fonction de risque cumulé avec Nelson Aalen et Breslow:
fit_na <- survfit(Surv(time, DEATH_EVENT) ~ 1, data = heart_df, ctype = 1,conf.type ="log")
plot(fit_na$time, fit_na$cumhaz, type="s", col="blue",pch=1, 
     xlab = TeX(r'($t$)'),
     ylab = TeX(r'($\hat{H}(t)$)'),
     main=TeX(r'($\hat{H}_{Breslow} = - \log(\hat{S}_{KM})$)'))
lines(fit_na$time, fit_na$cumhaz + fit_na$std.err * 1.96, col = "blue", lty = 2)  # Limite supérieure
lines(fit_na$time, fit_na$cumhaz - fit_na$std.err * 1.96, col = "blue", lty = 2)  # Limite inférieure
lines(H, col="red")
lines(H + fit1$std.err * 1.96, col = "red", lty = 2)  # Limite supérieure
lines(H - fit1$std.err * 1.96, col = "red", lty = 2)  # Limite inférieure
legend("topright", legend=c("Breslow", "Nelson Aalen"), col=c("blue", "red"), lwd=2)



# On essaie de trouver une loi qui approxime la fonction de survie
#################################################################
# H0 ; les observations suivent elles une loi : de weibull ? expo ? gamma ? normal ? 
# Charger les bibliothèques nécessaires
library(fitdistrplus)

# Liste des distributions à tester
distributions <- c("weibull", "gamma", "norm", "lnorm", "logis", "cauchy")

fit_results <- list()

for (dist in distributions) {
  cat("\n====================\n")
  cat("Distribution:", dist, "\n")
  start_params <- NULL
  if (dist == "exponential") {
    start_params <- list(rate = 1 / mean(fit1$surv))
  } else if (dist == "t") {
    start_params <- list(df = 10)
  }
  
  fit <- tryCatch({
    fitdist(fit1$surv, dist, start = start_params)
  }, error = function(e) {
    cat("Erreur d'ajustement pour la distribution:", dist, "\n")
    return(NULL)
  })
  
  if (!is.null(fit)) {
    fit_results[[dist]] <- list(
      fit = fit,
      gof = gofstat(fit, fitnames = dist)
    )
    
    print(summary(fit))
    
    plot(fit)
    readline(prompt = "Appuyez sur [Entrée] pour continuer au prochain graphique.")
  }
}

# Afficher les résultats
fit_results$gamma$fit





















library(survival)
library(fitdistrplus)

# Préparer les données de survie
surv_obj <- Surv(heart_df$time, heart_df$DEATH_EVENT)

# Ajustement des distributions
fit_exp <- fitdist(heart_df$time, "exp")
fit_weibull <- fitdist(heart_df$time, "weibull")
fit_lnorm <- fitdist(heart_df$time, "lnorm")

# Évaluation de la Goodness of Fit
gof_exp <- gof(fit_exp)
gof_weibull <- gof(fit_weibull)
gof_lognormal <- gof(fit_lnorm)

# Résultats des Goodness of Fit
gof_exp
gof_weibull
gof_lognormal


















survdiff(Y~diabetes,data=heart_df,rho=0)#Test du log-rank
survdiff(Y~diabetes,data=heart_df,rho=1)

survdiff(Y~high_blood_pressure,data=heart_df,rho=0)#Test du log-rank
survdiff(Y~high_blood_pressure,data=heart_df,rho=1)  # p<0.05 <=> Les distributions soddnt vraiment différentes





# Second Kaplan Meier (avec filtre)
t1 <- heart_df$time[heart_df$diabetes==1]
e1 <- heart_df$DEATH_EVENT[heart_df$diabetes==1]
Y1 <- Surv(t1,e1)

t2 <- heart_df$time[heart_df$diabetes==0]
e2 <- heart_df$DEATH_EVENT[heart_df$diabetes==0]
Y2 <- Surv(t2,e2)

plot(survfit(Y1~1), col="blue", main="Estimation Kaplan de la fonction de survie des patients\n atteint d'une insuffisance cardiaque des fumeur et non fumeurs.")
lines(survfit(Y2~1),col="red")
legend(0,0.2,legend=c("diabetes", "non diabetes"), col=c("blue", "red"),
       lty=c(1,1), cex = 0.7)




Y = Surv(heart_df$time, heart_df$DEATH_EVENT)
model_cox<- coxph(Y ~ age+anaemia+diabetes+ejection_fraction+serum_sodium+platelets+smoking+sex+serum_creatinine, data = heart_df)
summary(model_cox)
# exp(coef) : risque proportionnel de la var
#

# # Test de Schoenfeld pour l'hypothèse d'indépendance au temps
test_schoenfeld <- cox.zph(model_cox)

# Résultats du test de Schoenfeld
print(test_schoenfeld)
plot(res.c)
#Le test de Schoenfeld est utilisé pour vérifier l’hypothèse de proportionnalité des risques dans un modèle de Cox. Voici comment interpréter les résultats :
#Valeurs de p (p-values) : 
#H0: les résidus de Schoenfeld sont indé du tps => les coeff de reg sont cste avec le tps
#H1: les résidus de Schoenfeld ne sont pas indé du tps =>  les coeff de reg ne sont pas cst dans le tps 

#Si les p-values sont inférieures à 0,05, on ne rejette pas H0
#Si les p-values sont supérieures à 0,05, on rejette H0

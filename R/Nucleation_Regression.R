nucl.m <- function (nboot=1000, npoints_interp=100, origin, m, k0, k, filename_in, method="SA") {
  ##' @title Probabilistic nucleation model parameter estimation from microfluidics data.
  ##' @description Probabilistic nucleation model parameter estimation from microfluidics data.
  ##' \cr The function takes as input a 2-column data file (see Arguments) that gives the time (in seconds) and the fraction f0 of empty droplets.
  ##' \cr Nucleation model parameters m, k0 and k are estimated by bootstrapping, where nboot is the number of boostrap samples.
  ##' \cr Parameter estimation is performed by simulated annealing or Levenberg-Marquardt, using the
  ##' GenSA and minpack.lm packages respectively.
  ##' \cr It is recommended that the user uses the "SA" option initially with nboot=1, and
  ##' then uses the output values of m, k0 and k as initialisation values for the "LM" estimation
  ##' with nboot = 1000 typically.
  ##' \cr The function outputs the results in filename_in_nboot.csv.
  ######################################################
  ##' @author Florent Bourgeois, florent.bourgeois@@ensiacet.fr
  ##' @author Laboratoire de Genie Chimique, Toulouse, FRANCE (\url{http://www.lgc.cnrs.fr})
  ##' @references Bourgeois, F. Teychene, S. and Biscans, B., Applicability of probabilistic 
  ##' nucleation modelling for analysis of microfluidics data, KONA Powder and Particle Journal, 35, 2018.
  ##' @keywords nucleation modelling
  # Auteur: F. Bourgeois
  # Derni?re mise ? jour: 17.03.2017
  # Nom script R: Nucleation_Regression.R
  #######################################################
  # Dependencies
  ##' @importFrom minpack.lm nls.lm
  ##' @importFrom GenSA GenSA
  ##' @importFrom neldermead fminbnd
  #######################################################
  # Arguments d'entr?e:
  # -------------------
  # nboot = nombre de sous-?chantillons utilis?s par le boostrap, typiquement 500 ? 2000 pour une bonne estimation
  #        des intervalles de confiance des param?tres du mod?le.
  ##' @param nboot Number of boostrap samples, typically 500 to 2000.Default=1000
  # npoints_interp = nombre de points pour le calcul continu du mod?le dans la plage des mesures de x.
  ##' @param npoints_interp Number of points for continuous model display inside the experimental range. Default=100.
  # origin = c(valeur1, valeur2) si on force les donn?es ? passer par (temps=valeur1, f0=valeur2) comme 1er point exp?rimental.
  ##' @param origin = c(0, 1) if data must pass through origin (time=0, f0=1) ; Optional argument
  # m = valeur d'initialisation de m
  ##' @param m Initialisation of nucleation parameter m
  # k0 = valeur d'initialisation de k0 (s-1)
  ##' @param k0 Initialisation of nucleation parameter k0 (s-1)
  # k = valeur d'initialisation de k (s-1)
  ##' @param k Initialisation of nucleation parameter k (s-1)
  # filename_in = nom du fichier d'entr?e, de type .csv, sans l'extension .csv
  ##' @param filename_in = name of entry data file, .csv extension.
  ##' \cr Entry data file should like this:
  ##' \cr t;p
  ##' \cr 300;0.975
  ##' \cr 600;0.952
  ##' \cr...
  # method = "SA" ou "LM"
  ##' @param method = "SA" or "LM". Default="SA"
  # Arguments de sortie:
  # --------------------
  ##' @return None
  # ---------
  ##' @usage
  ##'
  ##'
  ##' res <- nucl.m(1000, 100, m=0.2, k0=7e-6, k=2e-4, origin=c(0,1), filename_in="S1_exp_5_60", method="SA")
  # Necessaire pour que la fonction soit disponible à tous les utilisateurs
  # -----------------------------------------------------------------------
  #' @export
  #######################################################

  ## Les librairies
  library("minpack.lm")
  library("GenSA")
  library("neldermead")

  # Nettoyage de la console
  cat("\014")

  # Lecture des donn?es depuis un fichier .csv avec une ligne en haut pour le nom des variables.
  # La lecture depuis un fichier .csv est g?n?rique et pr?f?rable ? celle d'un fichier Excel, car elle ne d?pend
  # notamment ni de la version d'Excel ni du package R de lecture des fichiers Excel.
  filename <- paste(filename_in,".csv", sep="")
  print(paste("Mesures dans:", filename))
  InputData <- read.table(filename, head=TRUE, sep=";")

  #print(summary(InputData))

  #
  if (missing(origin)) {
    ## On ne rajoute pas de point ? l'origine
    print("N.B.: pas de rajout de point ? l'origine")
    npoints <- length(InputData$t)
    x_exp <- array(data=0, dim=npoints)
    y_exp <- array(data=0, dim=npoints)
    x_exp <- InputData$t # mesures de la variable d?pendante (explicative)
    y_exp <- InputData$p # mesures de la variable ind?pendante (r?ponse)
  } else {
    ## On rajoute le point (t=origin[1],f0=origin[2]) au d?but du jeu de donn?es.
    print(paste("N.B.: Rajout d'un point ? l'origine pour temps=",origin[1],"; f0=",origin[2]))
    npoints <- length(InputData$t)+1
    x_exp <- array(data=0, dim=npoints)
    y_exp <- array(data=0, dim=npoints)
    x_exp[1] <- origin[1]
    y_exp[1] <- origin[2]
    x_exp[2:npoints] <- InputData$t[1:npoints-1] # mesures de la variable d?pendante (explicative)
    y_exp[2:npoints] <- InputData$p[1:npoints-1] # mesures de la variable ind?pendante (r?ponse)
  }
  print(c("Nombre de points de mesures:", npoints))
  print("****** Valeurs du temps (s)")
  #print(x_exp)
  print("****** Valeurs de f0")
  print(y_exp)
  var_exp <- var(y_exp)

  # Lecture des param?tres d'initialisation du mod?le et affichage du mod?le correspondant
  m_guess <- m
  k0_guess <- k0
  k_guess <- k
  param_guess <- c(m_guess, k0_guess, k_guess)
  print(paste('valeurs initiales des param?tres du mod?le : m =', m_guess,' k0 =', k0_guess,' k =', k_guess))
  y_guess <- model(x_exp, param_guess) # r?ponse du mod?le avec les valeurs initiales des param?tres
  plot(x_exp, y_exp, xlab="Time (s)",ylab="Fraction of empty droplets", main="Initialisation") # cr?ation du graphique
  lines(x_exp, y_guess, type="l", col="red")
  # test de la fonction ? minimiser
  z <-min.ssr(x_exp, y_exp, param_guess)
  print(paste("SSR initial =", z))

  # Interruption du code pour viosualiser l'initialisation
  n <- readline(prompt="Appuyer sur ENTER pour continuer...")
  n <- as.integer(n)

  # Calcul des intervalles de confiance des param?tres et du mod?le par Bootstrap (sous-?chantillonnage avec remplacement)
  # Les donn?es sont sous-?chantillonn?es, puis les param?tres sont estim?s pour chaque sous-?chantillon de Boostrap.>
  param_fitted <- array(data=NA, dim=c(nboot,3)) # D?claration du tableau de nboot lignes et 3 colonnes (pour les 3 param?tres du mod?le)
  print("Patience, calcul en cours...")
  for (i in seq(from=1, to=nboot, by=1)) {
    # print(paste("Boostrap iteration ", i, " / ",nboot))
    # Le bootstrap est r?alis? manuellement, c.?.d. sans appel ? des fonctions d?di?es.
    x_boot <- array(data=0, dim=npoints)
    y_boot <- array(data=0, dim=npoints)
    if (missing(origin)) {
      draw <- seq(from=1,to=npoints,by=1)
      draw <- sample(draw,size=length(draw),replace=TRUE)
      x_boot <- x_exp[draw]
      y_boot <- y_exp[draw]
    } else {
      draw <- seq(from=2,to=npoints,by=1)
      draw <- sample(draw,size=length(draw),replace=TRUE)
      x_boot[1] <- origin[1]
      y_boot[1] <- origin[2]
      x_boot <- c(x_exp[1],x_exp[draw])
      y_boot <- c(y_exp[1],y_exp[draw])
    }

    # r?gression avec le recuit simul? (Package GenSA: functions for generalised simulated annealing)
    # la fonction objectif (? minimiser) est la somme des carr?s des r?sidus, et est renvoy?e par la fonction min.rss
    param0 <- c(m_guess, k0_guess, k_guess)
    if (method == "SA") {
    ######### print("Estimation des param?tres du mod?le par recuit simul? - package GenSA")
    #set.seed(1234)
    global.min <- 0
    tol <- 1e-13
    #maxit <- 1e6
    res <- GenSA(par = param0,
                 fn = min.ssr,
                 x = x_boot,
                 y = y_boot,
                 lower=c(1e-10,1e-10,1e-10),
                 upper=c(100,1e-2,1e6),
                 control=list(
                              threshold.stop=global.min+tol,
                              verbose=TRUE)
                 )

    print(paste("SA iteration ", i, " ==> r?sultat: ", res$par))
    } else { # method = "LM"
    #########print("Estimation des param?tres du mod?le par Levenberg-Marquardt - package MINPACK")
    tol <- 1e-13
    res <- nls.lm(par = param0,
                  fn = residFun,
                  x = x_boot,
                  y = y_boot,
                  lower=c(1e-10,1e-10,1e-10),
                  upper=c(100,1e-2,1e6),
                  control=nls.lm.control(ftol = tol, ptol = tol, maxiter=1000))
    print(paste("LM iteration ", i, " ==> r?sultat: ", res$par[1], res$par[2], res$par[3], "(rss= ", res$deviance,")"))
    }
    plot(x_exp, y_exp, pch=1, col="blue", main = c('Bootstrap iteration ',i), xlab="Time (s)",ylab="Fraction of empty droplets") # cr?ation du graphique
    points(x_boot, y_boot, pch="O", col="red")
    y_mod <- model(x_exp, res$par)
    lines(x_exp,y_mod, type="l", col="red", lwd=2.5)

    # Sauvegarde des param?tres estim?s
    param_fitted[i,] <- res$par[]
  } # Fin de la boucle du bootstrap

  # graphes des param?tres
  plot(param_fitted[,1], param_fitted[,2], xlab="m",ylab="k0")
  plot(param_fitted[,1], param_fitted[,3], xlab="m",ylab="k")
  plot(param_fitted[,2], param_fitted[,3], xlab="k0",ylab="k")
  ## plot les diagrammes quantile-quantile et histogrammes des param?tres du mod?le
  qqnorm(param_fitted[,1], xlab="m")
  qqnorm(param_fitted[,2], xlab="k0")
  qqnorm(param_fitted[,3], xlab="k")
  hist(param_fitted[,1], xlab="m", main = "Histogram of m", breaks=sqrt(nboot))
  hist(param_fitted[,2], xlab="k0", main = "Histogram of k0", breaks=sqrt(nboot))
  hist(param_fitted[,3], xlab="k", main= "Histogram of k", breaks=sqrt(nboot))

  # valeurs moyennes des param?tres
  m_mean = mean(param_fitted[,1])
  k0_mean = mean(param_fitted[,2])
  k_mean = mean(param_fitted[,3])
  print(paste("Valeurs moyennes des param?tres estim?s:"))
  print(paste("m:", m_mean))
  print(paste("k0:" , k0_mean))
  print(paste("k:", k_mean))

  ## Calcul du temps d'induction
  # On cherche la racine du mod?le pour f0=1, qui renvoie un temps = t(f0_mod?le=1)
  # Le temps d'induction = t(f0_mod?le = 1) - t_exp[1]
  fun <- function(x) (induction_time(x,c(m_mean, k0_mean, k_mean)))
  res <- fminbnd(fun, x0=0, xmin=0, xmax=x_exp[1], options=list(tol=1e-10))
  print(paste("Induction time (temps auquel f0=1) =", neldermead.get(res, "xopt"), " s"))

  ## extraction des limites de confiance ? 95% des param?tres
  # cela n?cessite au pr?alable de trier les valeurs des param?tres par ordre croissant
  param_fitted_sorted <- array(data=NA, dim=c(nboot,3))
  param_fitted_sorted[,1] <- sort(param_fitted[,1])
  param_fitted_sorted[,2] <- sort(param_fitted[,2])
  param_fitted_sorted[,3] <- sort(param_fitted[,3])
  # print(param_fitted_sorted[])

  n_median <- ceiling(nboot / 2) # indice de la valeur m?diane
  n_minus_95_percent_CI <- ceiling(0.025 * nboot) # indice de la limite ? 2.5%
  n_plus_95_percent_CI <- ceiling(0.975 * nboot) # indice de la limite ? 97.5%
  m_median <- param_fitted_sorted[n_median,1]
  m_minus_95_percent <- param_fitted_sorted[n_minus_95_percent_CI,1]
  m_plus_95_percent <- param_fitted_sorted[n_plus_95_percent_CI,1]
  k0_median <- param_fitted_sorted[n_median,2]
  k0_minus_95_percent <- param_fitted_sorted[n_minus_95_percent_CI,2]
  k0_plus_95_percent <- param_fitted_sorted[n_plus_95_percent_CI,2]
  k_median <- param_fitted_sorted[n_median,3]
  k_minus_95_percent <- param_fitted_sorted[n_minus_95_percent_CI,3]
  k_plus_95_percent <- param_fitted_sorted[n_plus_95_percent_CI,3]
  print("Intervalles de confiance des param?tres ? 95% de confiance:")
  print(paste("m:",m_minus_95_percent, m_plus_95_percent))
  print(paste("k0:",k0_minus_95_percent, k0_plus_95_percent))
  print(paste("k:",k_minus_95_percent, k_plus_95_percent))

  ## Calcul de deltaT critique
  Nb_max_site_actif <- qpois(0.975, lambda=m_mean) # Limite ? 95% de loi de Poisson de moyenne m_mean
  print(paste("Nombre max de sites actifs = ", Nb_max_site_actif))
  dt_critique <- 1/(max(k0_mean,Nb_max_site_actif*k_mean))
  print(paste("detlaT critique = ", dt_critique, " s"))

  ## Calcul des limites de confiance du mod?le, sur une port?e entre x_min et x_max, ? d?finir
  ############################################################################################
  min_x <- min(x_exp)
  max_x <- 2*max(x_exp) # On pr?dit le mod?le sur 2 fois la port?e exp?rimentale.
  step_x <- (max_x - min_x) / npoints_interp
  x_model <- seq(from=min_x, to=max_x, by=step_x)
  n_model <- length(x_model)
  #print(x_model)
  min_model <- array(data=NA, dim=n_model)
  minus_95_percent_model <- array(data=NA, dim=n_model)
  mean_model <- array(data=NA, dim=n_model)
  plus_95_percent_model <- array(data=NA, dim=n_model)
  max_model <- array(data=NA, dim=n_model)
  var_point <- array(data=NA, dim=n_model)
  for (i in seq(from=1, to=n_model, by=1)) {
    # At every single extrapolated point i, we calculate the model predictions for each set of boostrapped model parameters,
    # and extract point variance and 95% C.I.
    model_prediction <- array(data=NA, dim=nboot)
    for (j in seq(from=1, to=nboot, by=1)) {
      model_prediction[j] <- model(x_model[i], param_fitted[j,])
    }
    model_prediction <- sort(model_prediction)
    min_model[i] <- min(model_prediction)
    mean_model[i] <- model_prediction[n_median]
    max_model[i] <- max(model_prediction)
    var_point [i] <- var(model_prediction)
    minus_95_percent_model[i] <- mean_model[i] + qt(0.025, nboot-1) * sqrt(var_point[i]) #model_prediction[n_minus_95_percent_CI]
    plus_95_percent_model[i] <- mean_model[i] + qt(0.975, nboot-1) * sqrt(var_point[i]) #model_prediction[n_plus_95_percent_CI]
  }
  plot(x_exp, y_exp, main="95% confidence interval of the model")
  lines(x_model,mean_model, col="black", type="l")
  lines(x_model,minus_95_percent_model, col="blue", type="l", lty=2)
  lines(x_model,plus_95_percent_model, col="blue", type="l", lty=2)
  # plot the min and maximum values for the nboot simulations
  #lines(x_model,min_model, col="black", type="l", lty=3)
  #lines(x_model,max_model, col="black", type="l", lty=3)

  ## Prediction interval ==> predictor of Y | X = x
  n_param <- length(param0)
  y_model <- mean_model
  # Calcul de la variance du mod?le
  var_model_fitted <- array(data=0, dim=nboot)
  for (j in seq(from=1, to=nboot, by=1)) {
    for (i in seq(from=1, to=npoints, by=1)) {
      var_model_fitted[j] <- var_model_fitted[j] + (model(x_exp[i], c(param_fitted[j,1], param_fitted[j,2], param_fitted[j,3])) - y_exp[i])^2
    }
    var_model_fitted[j] <- var_model_fitted[j] / (npoints - n_param)
  }
  hist(var_model_fitted[], breaks=50, main="Histogramme de l'estimateur de la variance du mod?le", xlab="Estimateur de la variance du mod?le")
  var_model <- sum(var_model_fitted[])/nboot
  # Tri croissant des estimateurs de la variance d mod?le pour chaque run bootstrap pour obtenir
  # l'intervalle de confiance de la variance du mod?le
  var_model_fitted <- sort(var_model_fitted)
  var_model_minus_95_percent_CI <- NULL
  var_model_plus_95_percent_CI <- NULL
  var_model_minus_95_percent_CI <- var_model_fitted[n_minus_95_percent_CI] # indice de la limite ? 2.5%
  var_model_plus_95_percent_CI <- var_model_fitted[n_plus_95_percent_CI] # indice de la limite ? 97.5%
  print(paste("SD fitted model=",sqrt(var_model), " - SD initial=", sqrt(var_exp)))
  print("Intervalle ? 95% de confiance de SD fitted model:")
  print(c(sqrt(var_model_minus_95_percent_CI), sqrt(var_model_plus_95_percent_CI)))
  #print(var_point)
  #print(var_model)
  var_total <- var_model + var_point
  minus_95_percent_data <- array(data=NA, dim=n_model)
  plus_95_percent_data <- array(data=NA, dim=n_model)
  minus_95_percent_data <- mean_model + qt(0.025, npoints) * sqrt(var_total)
  plus_95_percent_data <- mean_model + qt(0.975, npoints) * sqrt(var_total)
  plot(x_exp, y_exp, main="95% prediction interval of a new value", xlim=c(min(x_exp, x_model),max(x_exp,x_model)),
      ylim=c(min(minus_95_percent_data), max(plus_95_percent_data)), xlab="Time (s)",ylab="Fraction of empty droplets")
  lines(x_model,y_model, col="black", type="l")
  lines(x_model,minus_95_percent_data, col="black", type="l", lty=2)
  lines(x_model,plus_95_percent_data, col="black", type="l", lty=2)
  polygon(c(x_model,rev(x_model)),c(plus_95_percent_data,rev(minus_95_percent_data)),col=rgb(0.2, 0.2, 0.2, 0.5), border=NA)

  ## Sauvegarde des r?sultats dans un fichier .csv
  filename_out <- paste(filename_in,"_OUT.csv", sep="")
  print(paste("R?sultats sauvegard?s dans:", filename_out))
  write.table(t(c("m_mean", "k0_mean", "k_mean")), file=filename_out, sep=";", col.names=FALSE, row.names=FALSE)
  write.table(t(c(m_mean, k0_mean, k_mean)), file=filename_out, sep=";", col.names=FALSE, row.names=FALSE, append=TRUE)
  write.table(t(c("m_m?dian", "k0_m?dian", "k_m?dian")), file=filename_out, sep=";", col.names=FALSE, row.names=FALSE, append=TRUE)
  write.table(t(c(m_median, k0_median, k_median)), file=filename_out, sep=";", col.names=FALSE, row.names=FALSE, append=TRUE)
  write.table(t(c("m -95%", "k0 -95%", "k -95%")), file=filename_out, sep=";", col.names=FALSE, row.names=FALSE, append=TRUE)
  write.table(t(c(m_minus_95_percent, k0_minus_95_percent, k_minus_95_percent)), file=filename_out, sep=";", col.names=FALSE, row.names=FALSE, append=TRUE)
  write.table(t(c("m +95%", "k0 +95%", "k +95%")), file=filename_out, sep=";", col.names=FALSE, row.names=FALSE, append=TRUE)
  write.table(t(c(m_plus_95_percent, k0_plus_95_percent, k_plus_95_percent)), file=filename_out, sep=";", col.names=FALSE, row.names=FALSE, append=TRUE)
  # Sauvegarde des nboot jeux de param?tres
  write.table(paste("Valeurs des", nboot, "estimations des param?tres par boostrapping"), file=filename_out, col.names=FALSE, row.names=FALSE, append=TRUE)
  write.table(t(c("m","k0","k")), file=filename_out, sep=";", col.names=FALSE, row.names=FALSE, append=TRUE)
  write.table(param_fitted[], file=filename_out, sep=";", col.names=FALSE, row.names=FALSE, append=TRUE)
  # Sauvegarde de la r?ponse du mod?le
  write.table(paste("Valeurs des", n_model, "estimations du mod?le et intervalle de pr?diction estim?s par bootstrapping"), file=filename_out, col.names=FALSE, row.names=FALSE, append=TRUE)
  write.table(t(c("x_model", "-95% C.I.","mean value","+95% C.I.")), file=filename_out, sep=";", col.names=FALSE, row.names=FALSE, append=TRUE)
  result <- array(data=NA, dim=c(n_model,4))
  result[,1] <- x_model
  result[,2] <- minus_95_percent_data
  result[,3] <- y_model
  result[,4] <- plus_95_percent_data
  write.table(result, file=filename_out, sep=";", col.names=FALSE, row.names=FALSE, append=TRUE)

## Sauvegarde des nboot jeux de param?tres pour la simulation avec des intervalles de confiance et pr?diction
# Ce fichier affiche le nombre nboot de jeux de param?tres m, k0 et k,
# puis les nboot jeux de paramt?tres m, k0, k
filename_out <- paste(filename_in,"_nboot.csv", sep="")
print(paste("R?sultats sauvegard?s dans:", filename_out))
write.table(t(c("m", "k0", "k")), file=filename_out, sep=";", col.names=FALSE, row.names=FALSE, append=FALSE)
write.table(param_fitted[], file=filename_out, sep=";", col.names=FALSE, row.names=FALSE, append=TRUE)

  } # Fin de FitNucleationmodel_SA_Bootstrap

FitNucleationmodel_SA <- function (m_guess, k0_guess, k_guess, filename_in) {
## programme d'estimation des param?tres de nucl?ation par mesures microfluidiques.
## L'estimation est r?alis?e sur la variation du nombre de gouttes vides en fonction du temps.
## L'estimation des param?tres est r?alis?e avec le recuit simul? (Simulated Annealing).
## La fonction objectif minimis?e est la somme des carr?s des r?sidus, qui suppose des incertitudes gaussiennes du mod?le.
#
# Arguments d'entr?e:
# -------------------
# m_guess = valeur d'initialisation de m
# k0_guess = valeur d'initialisation de k0 (s-1)
# k_guess = valeur d'initialisation de k (s-1)
#
# Arguments de sortie:
# --------------------
#
#
# exemple: FitNucleationmodel_SA(0.1, 0.1, 0.1,"EmptyDroplets.csv")
#
############################################

# Nettoyage de la console
cat("\014")

# Lecture des donn?es depuis un fichier .csv avec une ligne en haut pour le nom des variables.
# La lecture depuis un fichier .csv est g?n?rique et pr?f?rable ? celle d'un fichier Excel, car elle ne d?pend
# notamment ni de la version d'Excel ni du package R de lecture des fichiers Excel.
InputData <- read.table(filename_in, head=TRUE, sep=";")
print(paste("Mesures dans:", filename_in))
summary(InputData)
x_exp <- InputData$t # mesures de la variable d?pendante (explicative)
y_exp <- InputData$p # mesures de la variable ind?pendante (r?ponse)
#print(x_exp)
#print(y_exp)
npoints <- length(x_exp)
print(c("Nombre de points de mesures:", npoints))

# Lecture des param?tres d'initialisation du mod?le et affichage du mod?le correspondant
param_guess <- c(m_guess, k0_guess, k_guess)
print(c('valeurs initiales des param?tres du mod?l: ', param_guess))
y_guess <- model(x_exp, param_guess) # r?ponse du mod?le avec les valeurs initiales des param?tres
plot(x_exp, y_exp, xlab="Time (s)",ylab="Fraction of empty droplets") # cr?ation du graphique
# test de la fonction ? minimiser
# z <-min.ssr(x_exp, y_exp, param_guess)
# print(z)

# r?gression avec le recuit simul? (Package GenSA: functions for generalised simulated annealing)
# la fonction objectif (? minimiser) est la somme des carr?s des r?sidus, et est renvoy?e par la fonction min.rss
print("Estimation des param?tres du mod?le par recuit simul? - package GenSA")
param0 <- c(m_guess, k0_guess, k_guess)
global.min <- 0
tol <- 1e-13
res <- GenSA(par = param0, fn = min.ssr, x = x_exp, y = y_exp, lower=c(0,0,0), upper=c(1,1,1), control=list(threshold.stop=global.min+tol,verbose=TRUE))
print(res[c("value","par","counts")])
#print(res$par[])
y_mod <- model(x_exp, res$par[])
lines(x_exp,y_mod, col="red", lwd=2.5)
} # Fin de FitNucleationmodel_SA

residFun <- function(x, y, param){
  #################################
  ## Calcule le vecteur des r?sidus du mod?le
  # x: mesures de la variable ind?pendante (explicative)
  # y: mesures de la variable d?pendante (r?ponse)
  # param: les param?tres du mod?le
  # La fonction model() est utilis?e pour calculer la r?ponse
  # du mod?le aux valeurs des mesures pour un jeu de param?tres donn?s.
  resid <- y - model(x, param);
  return(resid);
}


min.ssr <- function(x, y, param){
#################################
## Calcule la somme des carr?s des r?sidus du mod?le
# x: mesures de la variable ind?pendante (explicative)
# y: mesures de la variable d?pendante (r?ponse)
# param: les param?tres du mod?le
# Calcul de la somme des carr?s des r?sidus. La fonction model() est utilis?e pour calculer la r?ponse
# du mod?le aux valeurs des mesures pour un jeu de param?tres donn?s.
  #print(y)
  #print(model(x,param))
  resid <- y - model(x, param);
  ssq <- t(resid) %*% resid; # somme des carr?s des r?sidus
  return(ssq);
}

model <- function(x, param){
############################
## Calcule la r?ponse du mod?le au point x
# fit sur les gouttes vides f0

#print(param)

#######################################################
# Mod?le POUND et LA MER (1 site actif peut donner plus de 1 nucleus, et pas de nucl?ation homog?ne dans le gouttes qui contiennent
  # des sites actifs)
#  y = exp(-param[1])*(exp(-param[2]*x) + exp(param[1]*exp(-param[3]*x)) - 1);

#######################################################
# Mod?le AKELLA & FRADEN (1 site actif peut donner plus de 1 nucleus, et assignation incorrecte de la probabilit? de nucl?ation homog?ne et h?t?rog?ne des gouttes qui contiennent
# des sites actifs)
# y = exp(-param[1]) * (exp(-param[2]*x) * exp(param[1]*exp(-param[3]*x)));

#######################################################
# Notre mod?le (1 site actif peut donner n'importe quel nombre de nuclei (pas de troncature),
#  et assignation correcte de la probabilit? de nucl?ation homog?ne et h?t?rog?ne des gouttes qui contiennent
# des sites actifs)

## Solution par distr. de Poisson
    pmax <- 100
    y <- dpois(0,param[1])*dpois(0,param[2]*x)
    for (p in seq(from=1,to=10,by=1)) {
      y <- y + dpois(p,param[1])*(param[2]/(param[2]+p*param[3]))*dpois(0,param[2]*x)
      y <- y + dpois(p,param[1])*(p*param[3]/(param[2]+p*param[3]))* dpois(0,p*param[3]*x)
    }

    ##  Solution par formulation analytique
   # pmax <- 10
   # y = exp(-param[1]) * exp(-param[2]*x)
   # for (p in seq(from=1,to=pmax,by=1)) {
   #    y <- y + (param[1]^p * exp(-param[1])/factorial(p)) * (param[2]/(param[2]+p*param[3])) * exp(-param[2]*x)
   #    sum <- 0
   #    for (i in seq(from=0, to=p, by=1)) {
   #      sum <- sum + ((p*param[3]*x)^i) * exp(-p*param[3]*x) / factorial(i)
   #    }
   #    y <- y + (param[1]^p * exp(-param[1])/factorial(p)) * (p*param[3]/(param[2]+p*param[3])) * exp(-p*param[3]*x) / sum
   # }

  return(y)
}

induction_time <- function(x, param){
  # Cette fonction est utilis?e pour calculer le temps d'induction
  # Elle renvoie la valeur de f0(t) - 1.
  return(model(x,param)-1)
}

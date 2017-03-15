# Lectura de los datos desde github
file <- 'https://raw.githubusercontent.com/fhernanb/datos/master/orellana'
datos <- read.table(file=file, header=TRUE)
head(datos)
attach(datos)

# Analisis descriptivo ----------------------------------------------------

# Densidad de la variable respuesta, FIGURA 2
plot(density(produccion), xlim=c(0, 3100), lwd=3,
     main='', ylab='Densidad', ylim=c(0, 0.0012),
     xlab='Rendimiento de orellana (gr)')
rug(produccion)


# Aplicacion de fitDist ---------------------------------------------------
require(gamlss)
fit <- fitDist(produccion, data=datos, type="realplus")
fit$fits
fit

# Histograma y mejores ajustes con exp, wei, ga y pareto
par(mfrow=c(2, 2))
fitexp <- histDist(y=produccion, family=EXP, main="EXP")
fitwei3 <- histDist(y=produccion, family=WEI3, main="WEI3")
fitga <- histDist(y=produccion, family=GA, main="GA")
fitpareto2 <- histDist(y=produccion, family=GA, main="PARETO2")

# Mejorando los graficos anteriores, FIGURA 3
par(mfrow=c(2, 2))
hist(produccion,
     main=expression(paste("Exponencial(", mu, "=541.49)")),
     ylab='Densidad', freq=F, breaks=10, xlim=c(0, 3000),
     xlab='Rendimiento de orellana (gr)', ylim=c(0, 0.0015))
curve(dEXP(x, mu=exp(fitexp$mu.coef)), add=T, lwd=3)

hist(produccion,
     main=expression(paste("Weibull3(", mu, "=540.81, ", sigma,"=1.07)")),
     ylab='Densidad', freq=F, breaks=10, xlim=c(0, 3000),
     xlab='Rendimiento de orellana (gr)', ylim=c(0, 0.0015))
curve(dWEI3(x, mu=exp(fitwei3$mu.coef), 
            sigma=exp(fitwei3$sigma.coef)), add=T, lwd=3)

hist(produccion,
     main=expression(paste("Gamma(", mu, "=541.49, ", sigma,"=0.97)")),
     ylab='Densidad', freq=F, breaks=10, xlim=c(0, 3000),
     xlab='Rendimiento de orellana (gr)', ylim=c(0, 0.0015))
curve(dGA(x, mu=exp(fitga$mu.coef), 
          sigma=exp(fitga$sigma.coef)), add=T, lwd=3)

hist(produccion,
     main=expression(paste("Pareto2(", mu, "=541.49, ", sigma,"=0.97)")),
     ylab='Densidad', freq=F, breaks=10, xlim=c(0, 3000),
     xlab='Rendimiento de orellana (gr)', ylim=c(0, 0.0015))
curve(dPARETO2(x, mu=exp(fitpareto2$mu.coef), 
               sigma=exp(fitpareto2$sigma.coef)), add=T, lwd=3)



# Boxplots, FIGURA 4
par(mfrow=c(2, 3))
boxplot(produccion ~ Tiempo, las=1,
        ylab='Rendimiento de orellana (gr)', 
        xlab='Tiempo de aireaci?n (horas)')
boxplot(produccion ~ Temperatura, las=1,
        ylab='Rendimiento de orellana (gr)', 
        xlab='Temperatura (?C)')
boxplot(produccion ~ Humedad, las=1,
        ylab='Rendimiento de orellana (gr)', 
        xlab='Humedad (%)')
boxplot(produccion ~ PesoBolsa, las=1,
        ylab='Rendimiento de orellana (gr)', 
        xlab='Cantidad de sustrato (kilogramos)')
boxplot(produccion ~ cascara, las=1,
        ylab='Rendimiento de orellana (gr)', 
        xlab='Tipo de c?scara')



# Aplicacion gamlss -------------------------------------------------------

# Modelo horizonte
horizonte <- formula(~ cascara + PesoBolsa + Temperatura + 
                       Humedad + Tiempo + 
                       PesoBolsa * Temperatura + 
                       PesoBolsa * Humedad + 
                       PesoBolsa * Tiempo +
                       Temperatura * Humedad +
                       Temperatura * Tiempo +
                       Humedad * Tiempo +
                       I(PesoBolsa ^ 2) + 
                       I(Temperatura ^ 2) +
                       I(Humedad ^ 2) + 
                       I(Tiempo ^ 2))

require(gamlss)

# stepGAICALL.A se usa para obtener el mejor modelo
# para cada uno de los par?metros de la distribuci?n elegida
# se debe ingresar un modelo de partida y scopes para cada parametro
# en forma de lista.
# lower es el modelo m?s simple admitido
# upper es la estructura del modelo m?s complejo admitido
# sup es simplemente la f?rmula del m?s complejo

## Ajuste del modelo inicial con la distribucion EXP
exp0 <- gamlss(produccion ~ 1, data=datos, family=EXP())
exp1 <- stepGAICAll.A(exp0, trace=F,
                      scope=list(lower= ~ 1, upper=horizonte))

## Ajuste del modelo inicial con la distribucion NO
no0 <- gamlss(produccion ~ 1, sigma.fo= ~ 1, data=datos, family=NO)
no1 <- stepGAICAll.A(no0, trace=F,
                     scope=list(lower= ~ 1, upper=horizonte),
                     sigma.scope=list(lower= ~ 1, upper=horizonte))
no1 <- refit(no1)

## Ajuste del modelo inicial con la distribucion GA
ga0 <- gamlss(produccion ~ 1, sigma.fo= ~ 1, data=datos, family=GA())
ga1 <- stepGAICAll.A(ga0, trace=F,
                     scope=list(lower= ~ 1, upper=horizonte),
                     sigma.scope=list(lower= ~ 1, upper=horizonte))

## Ajuste del modelo inicial con la distribucion LOGNO
logno0 <- gamlss(produccion ~ 1, sigma.fo= ~ 1, data=datos, 
                 family=LOGNO(mu.link='log'))
logno1 <- stepGAICAll.A(logno0, trace=F,
                        scope=list(lower= ~ 1, upper=horizonte),
                        sigma.scope=list(lower= ~ 1, upper=horizonte))

## Ajuste del modelo inicial con la distribucion IG
ig0 <- gamlss(produccion ~ 1, sigma.fo= ~ 1, data=datos, family=IG())
ig1 <- stepGAICAll.A(ig0, trace=F,
                     scope=list(lower= ~ 1, upper=horizonte),
                     sigma.scope=list(lower= ~ 1, upper=horizonte))

## Ajuste del modelo inicial con la distribucion WEI3
wei30 <- gamlss(produccion ~ 1, sigma.fo= ~ 1, data=datos, family=WEI3())
wei31 <- stepGAICAll.A(wei30, trace=F,
                       scope=list(lower= ~ 1, upper=horizonte),
                       sigma.scope=list(lower= ~ 1, upper=horizonte))

## Ajuste del modelo inicial con la distribucion PARETO2
con1 <- gamlss.control(c.crit=0.001, n.cyc=10000)
pare0 <- gamlss(produccion ~ 1, sigma.fo= ~ 1, data=datos, family=PARETO2(),
                   control=con1)
pare1 <- stepGAICAll.A(pare0, trace=F,
                       scope=list(lower= ~ 1, upper=horizonte),
                       sigma.scope=list(lower= ~ 1, upper=horizonte))

## Ajuste del modelo inicial con la distribucion BCCGo
bccg0 <- gamlss(produccion ~ 1, sigma.fo= ~ 1, data=datos, family=BCCGo())
bccg1 <- stepGAICAll.A(bccg0, trace=F,
                       scope=list(lower= ~ 1, upper=horizonte),
                       sigma.scope=list(lower= ~ 1, upper=horizonte),
                       nu.scope=list(lower= ~ 1, upper=horizonte))
bccg1 <- refit(bccg1)

## Ajuste del modelo inicial con la distribucion GG
gg0 <- gamlss(produccion ~ 1, sigma.fo= ~ 1, data=datos, family=GG())
gg1 <- stepGAICAll.A(gg0, trace=F,
                     scope=list(lower= ~ 1, upper=horizonte),
                     sigma.scope=list(lower= ~ 1, upper=horizonte),
                     nu.scope=list(lower= ~ 1, upper=horizonte))
gg1 <- refit(gg1)

## Ajuste del modelo inicial con la distribucion GIG
gig0 <- gamlss(produccion ~ 1, sigma.fo= ~ 1, data=datos, family=GIG())
gig1 <- stepGAICAll.A(gig0, trace=F,
                      scope=list(lower= ~ 1, upper=horizonte),
                      sigma.scope=list(lower= ~ 1, upper=horizonte),
                      nu.scope=list(lower= ~ 1, upper=horizonte))
gig1 <- refit(gig1)

## Ajuste del modelo inicial con la distribucion BCTo
con1 <- gamlss.control(c.crit=0.001, n.cyc=10000)
bct0 <- gamlss(produccion ~ 1, sigma.fo= ~ 1, data=datos, family=BCTo(),
              control=con1)

bct1 <- stepGAICAll.A(bct0, trace=F,
                      scope=list(lower= ~ 1, upper=horizonte),
                      sigma.scope=list(lower= ~ 1, upper=horizonte),
                      nu.scope=list(lower= ~ 1, upper=horizonte),
                      tau.scope=list(lower= ~ 1, upper=horizonte))
bct1 <- refit(bct1)

## Ajuste del modelo inicial con la distribucion BCPEo
con1 <- gamlss.control(c.crit=0.001, n.cyc=10000)
bcpe0 <- gamlss(produccion ~ 1, sigma.fo= ~ 1, data=datos, family=BCPEo(),
                control=con1)
bcpe1 <- stepGAICAll.A(bcpe0, trace=F,
                       scope=list(lower= ~ 1, upper=horizonte),
                       sigma.scope=list(lower= ~ 1, upper=horizonte),
                       nu.scope=list(lower= ~ 1, upper=horizonte),
                       tau.scope=list(lower= ~ 1, upper=horizonte))
bcpe1 <- refit(bcpe1)

# Todos
AIC(exp1,
    no1, ga1, logno1, ig1, wei31, pare1,
    bccg1, gg1, gig1, 
    bct1, bcpe1, k=log(nrow(datos)))

plot(bcpe1) # descartado

# Solo los buenos
AIC(exp1,
    no1, ga1, ig1, wei31, pare1,
    bccg1, gg1, gig1, 
    bct1, k=log(nrow(datos)))

# updating

ga2 <- update(object=ga1, formula=~Temperatura + I(Temperatura^2) +
              Tiempo + cascara, what = c("mu"),
              evaluate = TRUE)

gg2 <- gamlss(produccion ~ Temperatura + I(Temperatura^2) + 
                Tiempo + I(Tiempo^2) + cascara, 
              sigma.fo= ~ 1, 
              nu.fo = ~ Temperatura,
              data=datos, family=GG())

gig2 <- gamlss(produccion ~ Temperatura + I(Temperatura^2) + 
                 Tiempo + I(Tiempo^2) + cascara, 
               sigma.fo= ~ Temperatura, 
               nu.fo = ~ 1,
               data=datos, family=GIG())

# Solo los modificados
AIC(ga1, gg1, gig1,
    ga2, gg2, gig2,
    k=log(nrow(datos)))



# Comparando los 3 mejores modelos ----------------------------------------
par(mfrow=c(1, 3), bg='white')
wp(gg2)
title("Gamma generalizada (GG)")
wp(ga2)
title("Gamma (GA)")
wp(gig2)
title("Inversa gausiana generalizada (GIG)")

AIC(gg2, ga2, gig2, k=log(nrow(datos)))  # Para obtener SBC


# Best model --------------------------------------------------------------
wp(ga2)
plot(ga2)
summary(ga2)

# Contornos para la media -------------------------------------------------
esp <- function(temp, tiempo) {
  x <- c(1, temp, temp^2, tiempo, ind)
  exp(sum(coef(ga2) * x))
}
esp <- Vectorize(esp)

k <- 100
temp <- seq(from=19, to=24, length.out=k)
tiem <- seq(from=4, to=6, length.out=k)

ind <- 0
z <- outer(temp, tiem, esp)
par(mfrow=c(1, 2))
contour(temp, tiem, z, las=1, lwd=2,
        main='Cascarilla entera',
        xlab='Temperatura (°C)', ylab='Tiempo aireacion (horas)')
ind <- 1
z <- outer(temp, tiem, esp)
contour(temp, tiem, z, las=1, lwd=2,
        main='Cascarilla molida',
        xlab='Temperatura (°C)', ylab='Tiempo aireacion (horas)')


# Contornos para la varianza ----------------------------------------------
vari <- function(temp, tiempo) {
  x <- c(1, temp, temp^2, tiempo, ind)
  coefi <- c(27.072, -2.874, 0.076, 1.380, 1.476)
  exp(sum(coefi * x))
}
vari <- Vectorize(vari)

k <- 100
temp <- seq(from=19, to=24, length.out=k)
tiem <- seq(from=4, to=6, length.out=k)

ind <- 0
z <- outer(temp, tiem, vari)
par(mfrow=c(1, 2))
contour(temp, tiem, z, las=1, lwd=2,
        main='Cascarilla entera',
        xlab='Temperatura (°C)', ylab='Tiempo aireaciOn (horas)')
ind <- 1
z <- outer(temp, tiem, vari)
contour(temp, tiem, z, las=1, lwd=2,
        main='Cascarilla molida',
        xlab='Temperatura (°C)', ylab='Tiempo aireaciOn (horas)')



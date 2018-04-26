#----hacerStockMatriz----

hacerStockMatriz <- function(n =5){
  
  ret <- list()#crear contenedor de matrices
  for (i in 1:n){#hacer n matrices
    #Obtener valores de las distribuciones de prob.
    S00 <- rbeta(1,775,225)
    F1  <- runif(1,0.25,0.29)
    F2  <- runif(1,0.51,0.60)
    S01 <- rbeta(1,16,71)
    S11 <- rbeta(1,227,256)
    S12 <- rbeta(1,210,258)
    S22 <- rbeta(1,36,75)
    
    #guardar matrices estocasticas en lista
    ret[[i]] <- matrix(c(S00,F1 ,F2,
                         S01,S11,0 ,
                         0  ,S12,S22),
                       byrow = TRUE, ncol = 3)
    
    colnames( ret[[i]] ) <- paste0("n",0:2)
    rownames( ret[[i]] ) <- paste0("n",0:2)
  }#for
  
  return(ret)
}# hacerStockMatriz



#-----Demografía, proyectar matriz estocástica----

transStochMat <- setRefClass("transStochMat",
 fields = list(
   matriz = "list",
   n0 = "numeric",
   t = "numeric",
   p = "matrix",
   umbral = "numeric",
   probExt = "matrix",
   incluirEtapa = "numeric"
 ),
 methods = list(
   initialize = function(matriz,
                         n0,
                         t,
                         p,
                         umbral,
                         incluirEtapa,
                         probExt
   ) {
     .self$matriz <- matriz
     .self$n0 <- n0
     .self$t <- t
     setanddone <- FALSE
     if (missing(incluirEtapa)){
       .self$incluirEtapa <- rep(1, length(n0))
       setanddone <- TRUE
     } else {
       .self$incluirEtapa <- incluirEtapa
       setanddone <- TRUE
     }
     if (missing(umbral)){
       .self$umbral <- 0.05*(sum(n0))
     } else {
       .self$umbral <- umbral
     }
     if (missing(p)){
       require(popbio)
       .self$p <- stoch.projection(matrices = .self$matriz,
                                   nreps = 5000,
                                   n0 = n0,
                                   tmax = t)
     }
     if (setanddone & missing(probExt)){
       require(popbio)
       .self$probExt <-
         stoch.quasi.ext(
           matrices = matriz,
           n0 = n0,
           Nx = .self$umbral,
           nreps = 5000,
           maxruns = 10,
           tmax = t,
           sumweight = .self$incluirEtapa,
           verbose = FALSE
         )
       
     }
   },
   
   plotN = function(){
     assign("op", par())
     original <- sum(n0*incluirEtapa)
     iters <- rowSums(p%*%incluirEtapa)
     xLabThis <- paste0("Número de individuos en t = ",t)
     
     hist(
       iters,
       main = "Tamaño de población",
       xlab = xLabThis,
       ylab = "Frecuencia",
       lwd = 2, las = 1,
       xlim = c(min(c(original,iters)), max(c(original,iters)))
     )
     abline(v=original,col="gray")
     abline(v=umbral,col="red",lwd=2)
     suppressWarnings(suppressMessages(par(op)))
   },
   
   darR0 = function(){
     require(popbio)
     valR0 <- stoch.growth.rate(matriz)
     return(list(
       approx = exp(valR0$approx),
       sim = exp(valR0$sim),
       simCI = exp(valR0$sim.CI)
     ))
   },
   
   plotExtProb = function(){
     matplot(
       probExt,
       xlab = "Tiempo",
       ylab = "Probabilidad de quasi-extinción",
       type = "l",
       lty = 1,
       col = rainbow(10)
     )
   },
   extProb = function(){
     pop <- p%*%incluirEtapa#solo etapas incluidas
     prob <- mean(ifelse(rowSums(pop) < umbral, 1, 0))
     return(prob)
   }
 )#methods
)#class
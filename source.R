####----ENCODING LATIN1



#-----Capítulo 1, Crecimiento-----


#----Capítulo 2, Demografía-----


tablaVida <- function(x.sup,Sx, mx=NULL){
  oldOp <- options()
  digs <- nchar(max(Sx))
  options(digits = digs)
  
  amplitud <- x.sup - c(0,x.sup[-length(x.sup)])
  ex <- numeric(length(Sx)  )
  
  Dx <- Sx[-length(Sx)] - Sx[2:length(Sx)] ; Dx[length(Sx)] <- NA
  lx <- Sx/Sx[1]
  dx <- lx[-length(lx)] - lx[2:length(lx)]; dx[length(lx)] <- NA
  qx <- dx/lx
  for( i in 1:length(Sx) ){
    ex[i] <- sum(Sx[i:length(Sx)])/Sx[i]-0.5*amplitud[i]
  }
  
  val = data.frame(x = x.sup,
                   S_x = Sx,
                   D_x = Dx,
                   l_x = lx,
                   d_x = dx,
                   q_x = qx,
                   e_x = ex)
  
  colnames(val) <- c("$x$", "$S_x$", "$D_x$", "$l_x$", "$d_x$", "$q_x$", "$e_x$")
  
  if (!is.null(mx)){  
    mxlx <- mx*lx
    xmxlx <- (x.sup-amplitud/2)*mx*lx
    R0 <- sum(mxlx)
    Tc <- sum(xmxlx)
    
    
    return(list(
      val=val,
      R0 = R0,
      Tc = Tc
      
    ))
  } else {
    return(val)
  }
  options(oldOp)
}

##---Demografía, plotSiphonariaCC---

plotSiphonariaCC <- function(path2csv) {
  
  if (missing(path2csv)){
  #get dinPob WD
  #wd <- dirname(rstudioapi::getSourceEditorContext()$path)
 # datosDir <- paste0(wd, "/datos/chp-02-siphonaria.csv")
  datosDir <- "datos/chp-02-siphonaria.csv"
  dat <- read.csv(datosDir)
  } else {
    dat <- read.csv(path2csv)
  }
  hdat <- hist(dat$diam_mayor, plot = FALSE)
  a <- hdat$mids[which.max(hdat$counts)]
  b <- hdat$mids[length(hdat$mids)]
  
  hist(dat$diam_mayor,
       breaks = "sturges",
       col = "gray",
       main = "Histograma de tallas de\n Siphonaria gigas",
       xlab = "Diámetro mayor (cm)",
       ylab = "Frecuencia",
       las = 1)
  
  lines(x = c(a,b),rep(max(hdat$counts), 2), col = "red", lwd = 2)
  text(mean(c(a,b)), max(hdat$counts)*0.9, "curva descendente")
}


##---Demografía, plotAquilegiaCC----

plotAquilegiaCC <- function(){
  plot(1:4, log(Sx),
       main = "Verificación de ajuste",
       xlab = "Clase de edad",
       ylab = expression(log(Sx)),
       pch = 21, bg = 1, las = 1)
  abline(a = coef(salida)[1], b = coef(salida)[2], lwd = 2)
}

##---Demografía, plotWhale----

plotWhale <- function(){
    library(popbio)
    
    library(diagram)
    
    data(whale)
    
    curves <- matrix(nrow = ncol(whale), ncol = ncol(whale), 0)
    
    #backward
    curves[1,2] <- curves[1,3] <- 0.95
    
    #Forward
    curves[3,2] <- curves[2,1] <-  curves[4,3] <- 0
    
    plotmat(
      whale,
      absent = 0,
      pos = c(4),
      curve = curves,
      name = colnames(whale),
      lwd = 1.8,
      cex.txt = 0.8,
      box.cex = 0.8,
      box.size = 0.08,
      arr.length = 0.5,
      box.type = "circle",
      box.prop = 1,
      box.lwd = 3,
      shadow.size = 0.01,
      self.cex = 0.8,
      self.lwd = 2.5,
      my = -0.075,
      mx = -0.01,
      relsize = 0.9,
      self.shiftx = c(0, -0.02, -0.02, -0.02),#self point arrow x shift
      self.shifty = -0.1,
      dtext = 0.85, #controls the position of arrow text relative to arrowhead.
      main = "Modelo poblacional de la ballena asesina",
      arr.lwd = 2.5
    )
    
    
    detach("package:diagram", unload = TRUE)
}


#-----Demografía, proyectar matriz----

transMat <- setRefClass("transMat",
  fields = list(matriz="matrix",n0="numeric", t = "numeric",p="list"),
  methods = list(
    initialize = function(matriz,n0,t,p){
      .self$matriz <- matriz
      .self$n0 <- n0
      .self$t <- t
      if (missing(p)){
        require(popbio)
        .self$p <- pop.projection(.self$matriz,n0,t)
      }
    },
    plotTM = function(){
      assign("op", par())
      layout(matrix(1:2,ncol=2))
      
      plot(
        x=1:length(p$pop.sizes),
        y=p$pop.sizes,
        main = "Tamaño de población",
        xlab = "Tiempo",
        ylab = "Número de individuos",
        type = "l", lwd = 2, las = 1
      )
      plot(
        x=1:length(p$pop.changes),
        y=p$pop.changes,
        main = expression(lambda),
        xlab = "Tiempo",
        ylab = "Tasa de multiplicación",
        type = "l", lwd = 2, las = 1,
        ylim = c(ifelse(min(p$pop.changes)>1,0.9,min(p$pop.changes)),max(p$pop.changes)*1.1)
      )
      abline(h=1,col="gray")
      suppressWarnings(suppressMessages(par(op)))
    },
    darLambda = function(){
      return(p$lambda)
    },
    darR0 =  function(){
      require(popbio)
      net.reproductive.rate(matriz)
    },
    darTc = function(){
      require(popbio)
      generation.time(matriz)
    },
    plotStage = function(){
      require(popbio)
      stage.vector.plot(p$stage.vectors,
          xlab = "Tiempo", ylab = "Proporción en cada etapa")
    }
)#methods
)#class


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
              nreps = 500,
              maxruns = 10,
              sumweight = .self$incluirEtapa,
              verbose = FALSE,
              tmax = t
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
      
      darGR = function(){
        require(popbio)
        valGR <- stoch.growth.rate(matriz)
        return(list(
          approx = exp(valGR$approx),
          sim = exp(valGR$sim),
          simCI = exp(valGR$sim.CI)
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
        pop <- p%*%incluirEtapa
        prob <- mean(ifelse(rowSums(pop) < umbral, 1, 0))
        return(prob)
      }
    )#methods
)#class
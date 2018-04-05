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

plotSiphonariaCC <- function() {
  #get dinPob WD
  #wd <- dirname(rstudioapi::getSourceEditorContext()$path)
 # datosDir <- paste0(wd, "/datos/chp-02-siphonaria.csv")
  datosDir <- "datos/chp-02-siphonaria.csv"
  dat <- read.csv(datosDir)
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
      main = "Whale population model",
      arr.lwd = 2.5
    )
    
    
    detach("package:diagram", unload = TRUE)
}


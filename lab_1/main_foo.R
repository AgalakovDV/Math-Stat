gauss_distr <- function(nSize){
  size = nSize; 
  column = 1 + 3.3 * log(size) 
  randX <- rnorm(n = size, mean = 0, sd = 1); 
  #randX 
  
  minEl = min(randX) 
  maxEl = max(randX) 
  
  hist(randX, probability = TRUE, 
       main = paste("Равномерное распределение."),
       breaks = column,
       axes = TRUE,
       plot = TRUE,
       col = "#00E0FF"
  ); 
  
  step = (maxEl - minEl) / column 
  
  x = array(c(TRUE, FALSE), dim = c(column)) 
  for (i in 1:column + 1) { 
    x[i] = round( minEl + i * step, 0) 
  } 
  curve(dnorm(x),
        lwd = 2.5,
        col = "#FF0095",
        add = TRUE)
  return(0)
}

coshi_distr <- function(nSize){
  size = nSize; 
  column = 1 + 3.3 * log(size) 
  randX <- rcauchy(n = size, location = 0, scale = 1) 
  #randX 
  
  minEl = min(randX) 
  maxEl = max(randX)
  
  hist(randX, probability = TRUE,
       main = paste("Распределение Коши."),
       breaks = column,
       col = "#00E0FF"
  );  
  
  step = (maxEl - minEl) / column 
  x = array(c(T,F), dim = c(column)) 
  for (i in 1:column + 1) { 
    x[i] = round( minEl + i * step, 0) 
  } 
  
  curve(dcauchy(x),
        lwd = 2.5,
        col = "#FF0095",
        add = TRUE
  )
  return(0)
}

puasson_distr <- function(nSize){
  size = nSize; 
  column = 1 + 3.3 * log(size) 
  
  randX <- rpois(n = size, lambda = 10) 
  
  #randX 
  hist(randX, probability = TRUE, 
       main = paste("Распределение Пуассона"),
       breaks = column,
       col = "#00E0FF" 
  ) 
  minEl = min(randX) 
  maxEl = max(randX) 
  step = (maxEl - minEl) / column 
  x = array(c(TRUE, FALSE),dim = c(column)) 
  for (i in 1:column + 1) { 
    x[i] = round( minEl + i * step, 0) 
  } 
  
  par(new= TRUE) 
  lines(x, dpois(x, lambda = 10), 
        lwd = 2.5,
        col = "#FF0095"
  )
  return(0)
}

laplas_distr <- function(nSize){
  size = nSize; 
  column = 1 + 3.3 * log(size) 
  D <- DExp(rate = 1/sqrt(2)) 
  randX <- r(D)(size) 
  column = 1 + 3.3 * log(size) 
  hist(randX, probability = TRUE,
       main = paste("Распределение Лапласа."),
       breaks = column,
       col = "#00E0FF"
  ); 
  minEl = min(randX) 
  maxEl = max(randX) 
  step = (maxEl - minEl) / column 
  x = array(c(TRUE, FALSE), dim = c(column)) 
  for (i in 1:column + 1) { 
    x[i] = minEl + i * step 
  } 
  curve(d(D)(x), 
        lwd = 2.5,
        col = "#FF0095",
        add = TRUE
  ) 
  return(0)
}

uniform_distr <- function(nSize){
  size = nSize; 
  column = 1 + 3.3 * log(size) 
  
  randX <- runif(n = size, -sqrt(3), sqrt(3)) 
  
  #randX 
  
  minEl = min(randX) 
  maxEl = max(randX) 
  
  hist(randX, probability = TRUE, 
       main = paste("Равномерное распределение."),
       breaks = column,
       col = "#00E0FF"
  ); 
  
  step = (maxEl-minEl) / column 
  
  x = array(c(TRUE, FALSE), dim = c(column)) 
  for (i in 1:column + 1) { 
    x[i] = round( minEl + i * step, 0) 
  } 
  
  curve(dunif(x, -sqrt(3), sqrt(3)),
        lwd = 2.5,
        col = "#FF0095",
        add = TRUE
  )
  return(0)
}

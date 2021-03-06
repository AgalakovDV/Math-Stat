size20 = 20;
size100 = 100;

# создание боксплотов
create_boxplots <- function(randX, randY){
  arrX = c(size20 + size100)
  arrSizes = c(size20 + size100)
  for (i in 1:size20) {
    arrSizes[i] = size20
    arrX[i] = randX[i]
  }
  for (i in size20:(size100 + size20)) {
    arrSizes[i] = size100
    arrX[i] = randY[i - size20 + 1]
  }
  
  boxplot( 
    arrX ~ arrSizes,
    horizontal = TRUE
  )
}

# подсчет числа выбросов
find_out <- function(randX){
  boxplot(randX);
  numInd <- which(randX %in% boxplot.stats(randX)$out)
  return(length(numInd))
}

# распределение Гаусса
gauss_distr <- function(param = 1) {
  if (param == 1) {
    randX <- rnorm(n = size20, mean = 0, sd = 1); 
    randY <- rnorm(n = size100, mean = 0, sd = 1);
    create_boxplots(randX, randY);
    return(0)
  }
  else {
    randX <- rnorm(n = param, mean = 0, sd = 1); 
    
    return(find_out(randX))
  }
}

# распределение Коши
coshi_distr <- function(param = 1){
  if (param == 1) {
    randX <- rcauchy(n = size20, location = 0, scale = 1) 
    randY <- rcauchy(n = size100, location = 0, scale = 1)
  
    create_boxplots(randX, randY);
    return(0)
  }
  else {
    randX <- rcauchy(n = param, location = 0, scale = 1) 
    
    return(find_out(randX))
  }
}

# распределение Пуассона
pois_distr <- function(param = 1){
  if (param == 1) {
    randX <- rpois(n = size20, lambda = 10) 
    randY <- rpois(n = size100, lambda = 10) 
  
    create_boxplots(randX, randY);
    return(0)
  }
  else {
    randX <- rpois(n = param, lambda = 10) 
    
    return(find_out(randX))
  }
}

# распределение Лапласа
laplas_distr <- function(param = 1){
  if (param == 1) {
    D <- DExp(rate = 1/sqrt(2)) 
    randX <- r(D)(size20) 
  
    D <- DExp(rate = 1/sqrt(2)) 
    randY <- r(D)(size100) 
  
    create_boxplots(randX, randY)  
    return(0)
  }
  else {
    D <- DExp(rate = 1/sqrt(2)) 
    randX <- r(D)(param) 
    
    return(find_out(randX))
  }
}

# равномерное распределение
uniform_distr <- function(param = 1){
  if (param == 1) {
    randX <- runif(n = size20, -sqrt(3), sqrt(3)) 
    randY <- runif(n = size100, -sqrt(3), sqrt(3)) 
  
    create_boxplots(randX, randY)
    return(0)
  }
  else {
    randX <- runif(n = param, -sqrt(3), sqrt(3)) 
    
    return(find_out(randX))
  }
}

# 1 - Гаусс
# 2 - Коши
# 3 - Лаплас
# 4 - Пуассон
# 5 - Равномерное
create_sel <- function(num_sel = 1, nSize = size20) {
  stop_ = 1000;
  start_ = 1;
  count = 0;
  if (num_sel == 1) {   #Гаусс
    if (nSize == 1) {
      gauss_distr();
      return(NULL);
    }
    for (i in start_:stop_) {
      count = count + gauss_distr(nSize)
    }
     return(count/nSize/1000)
  }
  if (num_sel == 2) {   #Коши
    if (nSize == 1) {
      coshi_distr();
      return(NULL);
    }
    for (i in start_:stop_){
      count = count + coshi_distr(nSize)
    }
     return(count/nSize/1000)
  }
  if (num_sel == 3) {    #Лапласс
    if (nSize == 1) {
      laplas_distr();
      return(NULL);
    }
    for (i in start_:stop_){
      count = count + laplas_distr(nSize)
    }
     return(count/nSize/1000)
  }
  if (num_sel == 4) {     #Пуассон
    if (nSize == 1) {
      pois_distr();
      return(NULL);
    }
    for (i in start_:stop_) {
      count = count + pois_distr(nSize)
    }
     return(count/nSize/1000)
  }
  if (num_sel == 5) {     #Равномерное
    if (nSize == 1) {
      uniform_distr();
      return(NULL);
    }
    for (i in start_:stop_) {
      count = count + uniform_distr(nSize)
    }
     return(count/nSize/1000)
  }
  return(NULL)
}






find_quar <- function(param = 1){
  if (param == 1) {
    x = qnorm(c(0.25, 0.75), 0, 1);
  }
  else if (param == 2) {
    x = qcauchy(c(0.25, 0.75), 0, 1);
  }
  else if (param == 3) {
    D <- DExp(rate = 1/sqrt(2));
    x = q(D)(c(0.25, 0.75), 1/sqrt(2));
  }
  else if (param == 4) {
    x = qpois(c(0.25, 0.75), 10);
  }
  else {
    x = qunif(c(0.25, 0.75), -sqrt(3), sqrt(3));
  }
  vec = c(4);
  vec[1] = x[1];
  vec[2] = x[2];
  vec[3] = x[1] - 1.5 * (x[2] - x[1]);
  vec[4] = x[2] + 1.5 * (x[2] - x[1]);
  
  if (param == 1) {
    prob = pnorm(vec[3], 0, 1) + (1 - pnorm(vec[4], 0, 1));
  }
  else if (param == 2) {
    prob = pcauchy(vec[3], 0, 1) + (1 - pcauchy(vec[4], 0, 1));
  }
  else if (param == 3) {
    prob =  p(D)(vec[3]/sqrt(2)) + (1 - p(D)(vec[4], 1/sqrt(2)));
  }
  else if (param == 4) {
    prob = ppois(vec[3], 10) + (1 - ppois(vec[4], 10));
  }
  else {
    prob = punif(vec[3], -sqrt(3), sqrt(3)) + (1 - punif(vec[4], -sqrt(3), sqrt(3)));
  }
  
  myVec = c(5);
  for (i in 1:4) {
    myVec[i] = vec[i];
  }
  myVec[5] = prob;
  
  return(myVec);
}
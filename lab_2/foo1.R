count_selection = 1000;

find_middle <- function(ar){
  return((sum(ar))/count_selection)
}


# 1 - Гаусс
# 2 - Коши
# 3 - Лаплас
# 4 - Пуассон
# 5 - Равномерное

create_sel <- function(num_sel = 1, nSize = 10) {
  
  if (num_sel == 1) {
    sel <- rnorm(n = nSize, mean = 0, sd = 1)                     #Гаусс
    return(sel)
  }
  if (num_sel == 2) {
    sel <- rcauchy(n = nSize, location = 0, scale = 1)       #Коши
    return(sel)
  }
  if (num_sel == 3) {
    D <- DExp(rate = 1/sqrt(2))                                #Лаплас
    sel <- r(D)(nSize)                                       #Лаплас
    return(sel)
  }
  if (num_sel == 4) {
    sel <- rpois(n = nSize, lambda = 10);                    #Пуассон
    return(sel)
  }
  if (num_sel == 5) {
    sel <- runif(n = nSize, -sqrt(3), sqrt(3))               #Равномерное
    return(sel)
  }
  return(NULL)
}


main_job <- function(num_sel = 1, nSize = 10){
  size = nSize;
  
  X = array(c(TRUE, FALSE), dim = c(count_selection)) 
  MedX = array(c(TRUE, FALSE), dim = c(count_selection)) 
  ZR = array(c(TRUE, FALSE), dim = c(count_selection)) 
  ZQ = array(c(TRUE, FALSE), dim = c(count_selection)) 
  ZTr = array(c(TRUE, FALSE), dim = c(count_selection)) 
  
  for (j in 1:count_selection){
    sel <- create_sel(num_sel, nSize);
    sel = sort(sel);
    
    #выборочное среднее 
    sumSel = 0 
    sumSel = sum(sel) 
    X[j] = sumSel/size 
    
    #выборочная медиана
    MedX[j] = (sel[size / 2] + sel[size / 2 + 1]) / 2 
    
    minEl = min(sel); 
    maxEl = max(sel); 
    #Полусумма экстремальных выборочных элементов: 
    ZR[j] = (minEl + maxEl) / 2 
    
    #полусумма квартилей 
    m = size %% 4
    if (m == 0) { 
      i = size / 4; 
    } 
    else { 
      i = floor(size / 4) + 1; 
    } 
    z0 = sel[i]; 
    zm = sel[size - i + 1]; 
    ZQ[j] = (z0 + zm)/2 
    
    #усеченное среднее 
    r = size / 4; 
    k0 = 1 / (size - 2 * r); 
    sumAr = 0; 
    arr = array(c(TRUE, FALSE), dim = c(size - 2*r)) 
    for (i in floor(r + 1):ceiling(size - r)) { 
      arr[i - r] = sel[i]; 
    } 
    sumAr = sum(arr) 
    ZTr[j] = k0 * sumAr   
  }
  
  vec = c(10);
  
  vec[1] = find_middle(X) 
  vec[2] = find_middle(MedX)
  vec[3] = find_middle(ZR)
  vec[4] = find_middle(ZQ)
  vec[5] = find_middle(ZTr)
  
  
  DX = array(c(TRUE, FALSE), dim = c(count_selection)) 
  for (i in 1:count_selection) { 
    DX[i] = (X[i] - vec[1]) * (X[i] - vec[1]) 
  } 
  vec[6] = find_middle(DX)
  
  DMedX = array(c(TRUE, FALSE), dim = c(count_selection)) 
  for (i in 1:count_selection) { 
    DMedX[i] = (MedX[i] - vec[2]) * (MedX[i] - vec[2]) 
  } 
  vec[7] = find_middle(DMedX) 
  
  
  DZR = array(c(TRUE, FALSE),dim = c(count_selection)) 
  for (i in 1:count_selection) { 
    DZR[i] = (ZR[i] - vec[3]) * (ZR[i] - vec[3]) 
  } 
  vec[8] = find_middle(DZR)
  
  
  DZQ = array(c(TRUE, FALSE),dim = c(count_selection)) 
  for (i in 1:count_selection) { 
    DZQ[i] = (ZQ[i] - vec[4]) * (ZQ[i] - vec[4]) 
  } 
  vec[9] = find_middle(DZQ)
  
  
  DZTr = array(c(TRUE, FALSE), dim = c(count_selection)) 
  for (i in 1:count_selection) { 
    DZTr[i] = (ZTr[i] - vec[5]) * (ZTr[i] - vec[5]) 
  } 
  vec[10] = find_middle(DZTr)

  return(vec)
}

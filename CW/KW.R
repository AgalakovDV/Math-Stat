#library(EEM)

drawAA <- function(vec = "", reg = "Africa"){}

#regions are Africa and North
region <- c("Africa", "North")

countFiles = 15

#ways to data #пути до файлов с данными
strAfrica = "D:\\R_proj\\Курсовая работа\\ToBazhenov\\Hexane_extr_Kivu_Lake\\"
strNorth = "D:\\R_proj\\Курсовая работа\\ToBazhenov\\VD_DOM_Permafrost\\"

#names of files of Africa #имена файлов Африки
arrStrAfrica <- c("1.1_70.",
                  "1.2_21",
                  "1.3_68",
                  "1.4_114",
                  "1.5_11",
                  "1.6_37",
                  "2.3_5 (400)",
                  "2.3_5 (600)",
                  "2.3_5",
                  "2.4_7",
                  "3.1_14",
                  "3.2_69",
                  "3.3_15 (600)",
                  "3.4_20(800)",
                  "3.4_20",
                  "3.5_43")

#names of files of North #имена файлов Севера
arrStrNorth <- c("1701",
                 "1702",
                 "1704",
                 "1706",
                 "1708_1to10",
                 "1708_1to20",
                 "1711",
                 "1712",
                 "1727",
                 "1728",
                 "1729",
                 "1730",
                 "1732",
                 "1733",
                 "1734")


#type of files #окончание пути до файла с данными, тип файла
strEnd = ".txt"

#amino acids
vectAA <- c()
vectAA[1:5] = 0

#names of Amino Acids
namesAA <- c("C","A","M","B","T")

#matrix of Amino Acids (AA) from different regions
aminoAcids <- matrix(data = vectAA, nrow = 2, ncol = 5)
colnames(aminoAcids) = namesAA
rownames(aminoAcids) = region

#diferent deltas from Africa's file and North's file
delta <- c(5, 2)

#start number of all files
startNum = 250

#coordinate of rectangles
crdRect <- matrix(data = c(320, 250, 310, 270, 270, 
                           350, 260, 320, 280, 280,
                           420, 380, 380, 300, 320,
                           480, 480, 420, 320, 350),
                      nrow = 5, ncol = 4)
colnames(crdRect) <- c("Left", "Right", "Down", "Up")
rownames(crdRect) <- namesAA
#View(coordsRect)

myK <- matrix(data = 0, nrow = countFiles, ncol = 2)
colnames(myK) <- region

#function draw rectangles - the area of peaks of intensity
#функция рисует прямоугольники - области пиков интенсивности
drawRect <- function()
{
  rect(crdRect["C","Left"], crdRect["C","Down"], crdRect["C","Right"], crdRect["C","Up"], col = NA, border = "red") 
  rect(crdRect["A","Left"], crdRect["A","Down"], crdRect["A","Right"], crdRect["A","Up"], col = NA, border = "green") 
  rect(crdRect["M","Left"], crdRect["M","Down"], crdRect["M","Right"], crdRect["M","Up"], col = NA, border = "brown") 
  rect(crdRect["B","Left"], crdRect["B","Down"], crdRect["B","Right"], crdRect["B","Up"], col = NA, border = "yellow") 
  rect(crdRect["T","Left"], crdRect["T","Down"], crdRect["T","Right"], crdRect["T","Up"], col = NA, border = "black") 
}


AfrAA <- matrix(data = 0, nrow = countFiles, ncol = 5)
colnames(AfrAA) <- namesAA

NorAA <- matrix(data = 0, nrow = countFiles, ncol = 5)
colnames(NorAA) <- namesAA

#find AA from file and add to table AfrAA or NorAA
findAA <- function(dataWithFluo, name = region[1], itr = 1, AfrAA = AfrAA, NorAA = NorAA)
{
  a <- matrix(0)
  #считываем матрицу
  dataMatr <- dataWithFluo[[1]]
  #собираем аминокислоты из "нужных" мест матрицы.
  #"нужные" места соответствуют прямоугольникам.
  if (name == region[1]) 
  {
    for (i in 1:5)
    {
      a <- dataMatr[(c(crdRect[namesAA[i],"Down"] - startNum + 1)):(c(crdRect[namesAA[i],"Up"] - startNum + 1)),
                    (c(crdRect[namesAA[i],"Left"] - startNum)/delta[1] + 1):(c(crdRect[namesAA[i],"Right"] - startNum)/delta[1] + 1)]
      vectAA[i] = sum(a)
    }
    AfrAA[itr, ] = vectAA
    return(AfrAA)
  } else if (name == region[2]) {
    for (i in 1:5)
    {
      a <- dataMatr[(c(crdRect[namesAA[i],"Down"] - startNum + 1)):(c(crdRect[namesAA[i],"Up"] - startNum + 1)),
                    (c(crdRect[namesAA[i],"Left"] - startNum)/delta[2] + 1):(c(crdRect[namesAA[i],"Right"] - startNum)/delta[2] + 1)]
      vectAA[i] = sum(a)
    }
    NorAA[itr, ] = vectAA
    return(NorAA)
  } else {
    return(-1)
  }
}

#draw some information about AA
drawAA <- function(vec = "", reg = region[1], idx = 1)
{
  if (length(vec) != length(vectAA))
  {
    return(-1)
  } else {
    max1 = max(vec)
    #print(max1)
    if (reg == region[1])
    {
      png(width = 534, height = 404, filename = paste(region[1], arrStrAfrica[idx], ".png"))
      barplot(vec, names.arg = namesAA, horiz = T, xlab = "intensity", ylab = "agent",
              las = 1, main = "Africa", xlim = c(0, max1*1.4), col = "green")
    } else if (reg == region[2]) {
      png(width = 534, height = 404, filename = paste(region[2], arrStrNorth[idx], ".png"))
      barplot(vec, names.arg = namesAA, horiz = T, xlab = "intensity", ylab = "agent",
              las = 1, main = "North", xlim = c(0, max1*1.4), col = "steelblue")
    } else {
      return(-2)
    }
    dev.off()
  }
}

#matrix with wavelengths #матрица длин волн
#Africa :: cutEX = 400:600, cutEM = 600:700
#North :: cutEX = 400:800, cutEM = 600:700
mWavelengts <- matrix(data = c(400,600,600,700,400,800,600,700), nrow = 2, ncol = 4, byrow = TRUE)
colnames(mWavelengts) = c("cutExLeft",
                          "cutExRight",
                          "cutEmLeft",
                          "cutEmRight" );
rownames(mWavelengts) = c("Africa", "North")
#View(mWavelengts)

findK <- function(vec)
{
  if (length(vec) != length(vectAA))
  {
    return(-1)
  }
  return((vec[1] + vec[2])/(vec[4] + vec[5]))
}


#выборочно строим информацию по файлу
example_foo <- function(k = 1, AfrAA = AfrAA, NorAA = NorAA)
{
  if (k < 1)
  {
    k = 1
    strName = arrStrAfrica[k]
    wholeWay = paste0(strAfrica, strName, strEnd)
    data <- readEEM(wholeWay)
    drawEEM(data, n = 1)
    dataCut <- cutEEM(data, 
                      cutEX = mWavelengts["Africa","cutExLeft"]:mWavelengts["Africa","cutExRight"],
                      cutEM = mWavelengts["Africa","cutEmLeft"]:mWavelengts["Africa","cutEmRight"])
    drawEEM(dataCut, n = 1)
    dataWithFluo <- delScattering(dataCut, rep = 0) 
    drawEEM(dataWithFluo, n = 1)
    drawRect()
    legend(x = 346.4, y = 407, namesAA, col = c("red","green","brown","yellow", "black"), 
           lty = c(1, 1, 1, 1, 1), bg = "white")
    #North
    strName = arrStrNorth[k]
    wholeWay = paste0(strNorth, strName, strEnd)
    data <- readEEM(wholeWay)
    drawEEM(data, n = 1)
    dataCut <- cutEEM(data, 
                      cutEX = mWavelengts["North","cutExLeft"]:mWavelengts["North","cutExRight"],
                      cutEM = mWavelengts["North","cutEmLeft"]:mWavelengts["North","cutEmRight"])
    drawEEM(dataCut, n = 1)
    dataWithFluo <- delScattering(dataCut, rep = 0) 
    drawEEM(dataWithFluo, n = 1)
    drawRect()
    legend(x = 346.4, y = 407, namesAA, col = c("red","green","brown","yellow", "black"), 
           lty = c(1, 1, 1, 1, 1), bg = "white")
    return(0)
  }
  for (k in 1:countFiles)
  {
    #Africa
    strName = arrStrAfrica[k]
    wholeWay = paste0(strAfrica, strName, strEnd)
    data <- readEEM(wholeWay)
    
    #обрезаем график
    dataCut <- cutEEM(data, 
                    cutEX = mWavelengts["Africa","cutExLeft"]:mWavelengts["Africa","cutExRight"],
                    cutEM = mWavelengts["Africa","cutEmLeft"]:mWavelengts["Africa","cutEmRight"])
  
    #удаляем лучи рэлеевского рассеяния
    dataWithFluo <- delScattering(dataCut, rep = 0) 
####    drawEEM(dataWithFluo, n = 1)
####    drawRect()
    
    AfrAA <- findAA(dataWithFluo, name = region[1], itr = k, AfrAA = AfrAA, NorAA = NorAA)
    
    #North
    strName = arrStrNorth[k]
    wholeWay = paste0(strNorth, strName, strEnd)
    data <- readEEM(wholeWay)

    #обрезаем график
    dataCut <- cutEEM(data, 
                      cutEX = mWavelengts["North","cutExLeft"]:mWavelengts["North","cutExRight"],
                      cutEM = mWavelengts["North","cutEmLeft"]:mWavelengts["North","cutEmRight"])

    #удаляем лучи рэлеевского рассеяния
    dataWithFluo <- delScattering(dataCut, rep = 0) 
####    drawEEM(dataWithFluo, n = 1)
####    drawRect()
    
    NorAA <- findAA(dataWithFluo, name = region[2], itr = k, AfrAA = AfrAA, NorAA = NorAA)
    
    drawAA(AfrAA[k,], region[1], k)
    drawAA(NorAA[k,], region[2], k)
    
    #draw hist
    aminoAcids[1,] = AfrAA[k,] + aminoAcids[1,]
    aminoAcids[2,] = NorAA[k,] + aminoAcids[2,]
    myK[k,1] = findK(AfrAA[k,])
    myK[k,2] = findK(NorAA[k,])
  }
  View(AfrAA)
  View(NorAA)
  drawAA(aminoAcids[region[1],], region[1], 20)
  drawAA(aminoAcids[region[2],], region[2], 20)
  
  View(aminoAcids)
  View(myK)
  return(myK)
}
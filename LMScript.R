
#####
# Load in trajectories, file must contain trajectories in columated format. 
#E.g. Columns 1,2,3,4... would correspond to trajectories X1,Y1,X2,Y2....

allTrajs = read.table("",sep = ",", header = TRUE )

#Total number of trajectories should then be number of columns divided by 2
N = length(allTrajs[1,])/2

#Compute 'max' which is the length of the longest trajectory, may be used later for computing lags from ETA and TA MSD
max = numeric(1)
for(i in 1:length(allTrajs[1,])){
  if(i == 1){
    max = length(allTrajs[,i][!is.na(allTrajs[,i])])
  }
  if(i != 1){
    lengthTraj = length(allTrajs[,i][!is.na(allTrajs[,i])])
    if(lengthTraj > max){
      max = length(allTrajs[,i][!is.na(allTrajs[,i])])
    }
  }
}

#For plotting it is easier if the trajectories start from a common origin at (0,0)
for(i in 1:N){
  allTrajs[, 2*i-1] = allTrajs[, 2*i-1] - allTrajs[1, 2*i-1]
  allTrajs[, 2*i] = allTrajs[, 2*i] - allTrajs[1, 2*i]
}

#It may be useful to make a plot of a colour-coded trajectory whereby colour corresponds to the progression of time
i = 7
x = allTrajs[, 2*i-1]
N = length(x[!is.na(x)])
#colfunc2 <- colorRampPalette(c("Neant", "Goyal", "Jones", "Ke"))
colfunc2 <- colorRampPalette(c("#9C4F96", "#FF6355", "#FBA949", "#FAE442", "#8BD448", "#2AA8F2"))
colours2 = colfunc2(N)
quartz()
plot(allTrajs[, 2*i-1],  allTrajs[, 2*i], col = colours2, pch = 20, xaxt = "n", main= "", yaxt = "n", xlab= "", ylab = "", cex = 3.0)
axis(1,las= 1, tck=0.02,mgp=c(3, .5, 0), cex.axis = 1.8)
axis(2,las= 1, tck=0.02,mgp=c(3, .5, 0), cex.axis = 1.8)


##########################################
## Compute the ETA MSD
##########################################


#I am going to edit the code so as to loop through, collect all the displacements for a particular tau, compute the variance, then delete the displacements
#and move on to the next tau. 
maxLag = max/2 #Pick a max lag
counter = numeric(0)
counter = rep(1,round(max))
allDisplacements = numeric(1)
N = length(allTrajs[1,])/2
#Set the total experiment time is the trjactories are the same length, else use a smaller time
totalTime = 2338 #For F127 NHS trajs, 2338 for Diffusion. 
maxLag = round(totalTime)

var_spreadr = numeric(maxLag)
SEOM_TAEAMSD = numeric(maxLag)

counterDisp = 1
counterVar = 1

for(i in 1:2338){
  for(j in 1:N){
    #xCoords = allTrajs[,j][!is.na(allTrajs[,j])]
    xCoords = allTrajs[,2*j-1][!is.na(allTrajs[,2*j-1])]
    yCoords = allTrajs[,2*j][!is.na(allTrajs[,2*j])]
    #plot(xCoords, yCoords, type = "l")
    totalTime = length(xCoords)
    #maxLag = round(totalTime/10)
    for(t in seq(1,totalTime - i,1)){
      
      #xDisp = xCoords[t + i] - xCoords[t]
      #yDisp = yCoords[t + i] - yCoords[t]
      rDisp = sqrt((yCoords[t+i] - yCoords[t])^(2) + (xCoords[t+i] - xCoords[t])^(2))
      allDisplacements[counterDisp] = rDisp 
      counterDisp = counterDisp + 1
      
    }
    #print(j)
    
  } #Go through all trajs for a single lag
  
  var_spreadr[counterVar] = var(allDisplacements, na.rm = TRUE)
  SEOM_TAEAMSD[counterVar] =  sd(allDisplacements, na.rm = TRUE)/(sqrt(length(allDisplacements)))
  counterVar = counterVar + 1
  
  #Reset
  allDisplacements = numeric(1)
  counterDisp = 1
  #Update the user
  if(j %% 5 == 0){
    print(j)
  }
  if(i %% 1 == 0){
    print(i)
  }
  if(i %% 10 == 0){
    istring = as.character(i)
    save = paste("ETAMSD_F127NHSPeptide_motile_","istring", ".txt", sep="")
    write.table(var_spreadr, save, sep = ",")
  }
}

##########################################
## END ETA MSD 
##########################################

##########################################
## EA MSD 
##########################################


Ntot = N
xSQR = numeric(max)
ySQR = numeric(max)
sumXSQR = numeric(1)
sumYSQR = numeric(1)
xAVG = numeric(max)
yAVG = numeric(max)
EAMSD = numeric(max)
index = seq(1,N,1)

#Use just y. 
for(j in 1:max){
  #These are the x and y coords across all N trajs at a particular time j
  xCoords = allTrajs[j,2*index-1][!is.na(allTrajs[,2*index-1])]
  yCoords = allTrajs[j,2*index][!is.na(allTrajs[,2*index])]
  xCoords = na.omit(xCoords)
  yCoords = na.omit(yCoords)
  
  xAVG[j] = sum((xCoords))/length(xCoords)
  yAVG[j] = sum((yCoords))/length(yCoords)
  
  for(i in 1:length(xCoords)){
    sumXSQR = sumXSQR + (xCoords[i])^(2)
    sumYSQR = sumYSQR + (yCoords[i])^(2)
  }
  xSQR[j] = sumXSQR/length(xCoords)
  ySQR[j] = sumYSQR/length(yCoords)
  sumXSQR = 0
  sumYSQR = 0
}

EAMSD = xSQR + ySQR - (xAVG^(2) + yAVG^(2))
#EAMS = <x> + <y> - (<x>^2 + <y>^{2})
EAMSDX = xSQR - xAVG^(2)
EAMSDY = ySQR - yAVG^(2)

##########################################
## EA MSD END
##########################################

##########################################
## Compute the TA MSD
##########################################
N = length(allTrajs[1,])/2

round(0.1*length(allTrajs[,1]))
sum = numeric(0)
sum = 0
Result_R = numeric(1)
TAMSD_R <- data.frame(matrix(ncol = N, nrow = 500)) #Will hold all the TA MSD values
counter = 1
lim1 = 1
lim2 = 1000
for(j in 1:N){

  #j = 2
  xCoords = allTrajs[,2*j-1][!is.na(allTrajs[,2*j-1])]
  yCoords = allTrajs[,2*j][!is.na(allTrajs[,2*j])]
  #xCoords = allTrajs[lim1:lim2,2*j-1][!is.na(allTrajs[lim1:lim2,2*j-1])]
  #yCoords = allTrajs[lim1:lim2,2*j][!is.na(allTrajs[lim1:lim2,2*j])]
  #plot(xCoords, yCoords, type = "l")
  
  totalTime = length(yCoords[!is.na(yCoords)])
  #maxLag = round(totalTime/2)
  maxLag= 1000 #Similar to ETA MSD can adjust max lag to change compute time. 
  for(i in seq(1,maxLag,1)){
    for(t in seq(1,totalTime - i,1)){
      r = sqrt((yCoords[t+i] - yCoords[t])^(2) + (xCoords[t+i] - xCoords[t])^(2))
      #x = sqrt((xCoords[t+i] - xCoords[t])^(2))
      #y = sqrt((yCoords[t+i] - yCoords[t])^(2))
      sum = sum + (r)^(2)
      #sumX = sumX + (x)^(2)
      #sumY = sumY + (y)^(2)
    }
    
    Result_R[i] = sum/(totalTime - i)
    sum = 0
  }
  
  for(p in 1:length(Result_R)){
    TAMSD_R[p,counter] = Result_R[p]
  }
  counter = counter + 1
  Result_R = numeric(1)
  
  if(j %% 1 == 0){
    print(counter)
  }
  
}

##########################################
## TA MSD end
##########################################





#2oP Bootstrap, Author: Cam MacKenzie
#cmackenzie@atwaterresources.com
#last modified: Jan 2015
#Coded from 2oP macro 
#========================================================

# igr.boot <- function(DATA, TAXA,...)

#Input and Output file names
#Crop specific dates for use in IGR model dates
#Be sure individual Mass values are appropriate for species
#  size classes used
#Format Date
#Set appropriate # of Production Intervals


graphics.off()
rm(list=ls())
# set.seed(100) # sets randomization seed for bootstrapping, turn on if you want the same results every time you run the code 
#setwd("~/Rdata") #Set working directory


#Define filename for output to working directory
file.name.out = "Out_Prosimulium_7.2011.csv"


# Define samples to crop (due to outliers) for igr calculation  
# enter 0 (don't crop) or crop samples as follows c(1,9) or 9;   
# removes first and ninth sample or ninth sample only, respectively
# for igr calculation 
crop.dates = 0

# Define number of Monte Carlo simulations to run
nBoot= 1000

# Define number replicates for the bootstrap (usually 6). 
# This variable alters the number of replicates the bootstrap
# produces per sampling date
nReplicates = 5

# Input data 
#==============================================================================
# Data (.csv) MUST be in the following format, in particular the date (i.e. %m/%d/%y), 
# but number of columns before size class data should also be the same:
# DATE      SAMPLE  High taxon  Low taxon #/S   1 2  3  4	5	6.....
# 6/18/2012	-0.42	  Baetidae	  Acen	    0	    0	0	0	0	0	0
# 6/18/2012	-0.42	  Baetidae	  Acen	    44	  44.44444	0	0	0	0	0
# 6/18/2012	-0.17	  Baetidae	  Acen	    0 	  0	0	0	0	0	0
# 6/18/2012	0	      Baetidae	  Acen	    0	    0	0	0	0	0	0
# 6/18/2012	-0.17	  Baetidae	  Acen	    44	  44.44444	0	0	0	0	0
# 6/18/2012	0	      Baetidae	  Acen	    0	    0	0	0	0	0	0
# 6/27/2012	-0.17	  Baetidae	  Acen	    1022	577.77772	444.4444	0	0	0	0
# 6/27/2012	-0.38	  Baetidae	  Acen	    311	  266.66664	44.44444	0	0	0	0
# 6/27/2012	-0.42	  Baetidae	  Acen	    889	  799.99992	88.88888	0	0	0	0
#...

allData = read.csv("Prosimulium_7_2011.csv", header = TRUE) 

#Input mass for respective size classes (e.g. 1,2,3,4mm...)
ind.mass = c(0.00003, 0.00025, 0.00084, 0.00200, 0.00392, 0.00678, 0.01078, 0.01612, 0.02299, 0.03157, 0.04206, 0.05466, 0.06955, 0.08694,	0.10701, 0.12997, 0.15599, 0.18529, 0.21805, 0.25447, 0.29473, 0.33905,	0.38761, 0.44060, 0.49822, 0.56068, 0.62815, 0.70084, 0.77895, 0.86266,	0.95218, 1.04769, 1.14940, 1.25751, 1.37220, 1.49367, 1.62212, 1.75774, 1.90074, 2.05130, 2.20963, 2.37592, 2.55036, 2.73315) #mass for each size class

# Bootstrap data
#==============================================================================
# Bootstrap function
boots = function(allData){
  # for each sampling day resample data
  boot.data = 0.
  for (j in 1:length(levels(as.factor(allData$DATE)))){ # for each date
    
    date.sub <- subset(allData, DATE == levels(as.factor(allData$DATE))[j])
    #resample rows #s fo date.sub
    date.samp <- sample( 1:length(date.sub[,1]), size=nReplicates, replace=TRUE )
    #compile resampled data
    date.samp <- date.sub[date.samp,]
    if (j==1) {
      boot.data = date.samp
    } else {
      boot.data = rbind(boot.data,date.samp)
    }
  }
  #clean up and output data
  rownames(boot.data) = NULL
  return (boot.data)
}   

bootsData = boots(allData)

# Monte Carlo simulations
#==============================================================================
# To test calculations compared to excel macro use:
# bootsData = read.csv("2oP_boots_data.csv", header = TRUE) 

MCsims = 0. #table to contain bootstrap results (B,N,P,igr,P_B)

for (n in 1:nBoot){
  
  bootsData = boots(allData)
  
  size.classes = bootsData[,6:length(bootsData)] #extract size class data
  sum.classes = colSums(size.classes) # sum columns
  
  # loop to determine maximum size class, last column with data
  for (j in length(sum.classes):1){
    if (sum.classes[j]>0){
      no.classes = j
      break
    }
  }
  #trim data - remove extra columns with no data
  bootsData = bootsData[,1:(no.classes+5)]
  
  #format dates
  bootsData$DATE = as.Date(as.character(bootsData$DATE), format = "%m/%d/%Y")
  
  # Calculate Pararmeters N,B,P,igr
  #==============================================================================
  #no.classes = number size classes = k
  dates = levels(as.factor(bootsData$DATE)) #  = i
  
  
  ## 1) Calculate Abundance (N) table
  N.tab = matrix( ncol = no.classes+1) #nrow = length(dates)
  
  for (i in 1:length(dates)){
    temp.dat = 0.
    date.sub = subset(bootsData, DATE == dates[i])
    temp.dat[1] = dates[i]
    for (k in 1:no.classes){
      #calculate mean of bootstrapped size class for date i, rounded as in macro
      temp.dat[k+1] = as.numeric(sum(date.sub[k+5]))/length(date.sub[,1]) #values are rounded as in macros
    }
    N.tab = rbind(N.tab,temp.dat)
    
    
  }
  
  N.tab = N.tab[-1,]
  rownames(N.tab) = NULL
  colnames(N.tab) = c("Date", 1:no.classes)
  N.tab[,1] = as.Date(N.tab[,1]) 
  
  
  
  ## 2) Calculate Biomass (B) table
  Biomass.tab = 0.
  
  for (i in 1:length(dates)){
    temp.dat = 0.
    temp.dat[1] = dates[i]
    temp.dat[2] = sum(as.numeric(N.tab[i,1:no.classes+1]))
    for (k in 1:no.classes){
      temp.dat[k+2] = as.numeric(N.tab[i,k+1])*ind.mass[k]
    }
    temp.dat = c(temp.dat,sum(as.numeric(temp.dat[1:no.classes+2]))/as.numeric(temp.dat[2])) #calculate avg mass
    temp.dat = c(temp.dat,sum(as.numeric(temp.dat[1:no.classes+2]))) #calculate B mg/m2
    
    
    if (i==1) {
      Biomass.tab = temp.dat
    } else {
      Biomass.tab = rbind(Biomass.tab,temp.dat)
    }
  }
  
  rownames(Biomass.tab) = NULL
  colnames(Biomass.tab) = c("Date", "Abundance", 1:no.classes, "avg mass", "Biomass") 
  
  ## 3) Estimate igr 
  
  #define D.elapse
  dates = as.Date(dates, format = "%Y-%m-%d")
  
  D.elapse = 0.
  for (j in 1:(length(dates)-1)){
    D.elapse[j] = dates[j+1]-dates[1]
  }
  D.elapse = c(0,D.elapse)
  
  
  #log transform avg mass
  ln.avg.mass = log(as.numeric(as.character(Biomass.tab[,"avg mass"])))
  
  #Crop dates for appropriate IGR model (specified by crop.dates above)
  if (crop.dates[1] != 0){
    D.elapse = D.elapse[-crop.dates]
    ln.avg.mass = ln.avg.mass[-crop.dates]
  }
  
  # Run linear model
  growth.lm = lm(ln.avg.mass ~ D.elapse)
  
  #summary.lm(growth.lm)
  
  # Diagnostic plots, if wanted
  #plot(growth.lm)
  
  #plot fit, if wanted
  #   plot(ln.avg.mass ~ D.elapse) #CHECK W/CAM REGARDING "[-9]"
  #   abline(growth.lm)
  #   
  #Instaneaous growth rate (igr)
  igr = as.numeric(coef(growth.lm)[2]) 
  intercept = coef(growth.lm)[1]
  
  
  ## 4) Estimate P
  P.tab = 0.
  Biomass = as.numeric(Biomass.tab[,"Biomass"])
  for (i in 1:length(dates)){
    temp.dat = 0.
    temp.dat[1] = dates[i+1] - dates[i]
    temp.dat[2] = mean(Biomass[i:(i+1)])*igr*as.numeric(temp.dat[1])
    if (i==1) {
      P.tab = temp.dat
    } else {
      P.tab = rbind(P.tab,temp.dat)
    }
    
  }
  #clean up P.tab
  P.tab = P.tab[-length(P.tab[,2]),]
  rownames(P.tab) = NULL
  colnames(P.tab) = c("d", "Pint") 
  
  #Calculate P and P/B (P_B), and Pint
  P = sum(P.tab[,"Pint"])
  P_B = P/mean(Biomass)

  Pint.1=P.tab[1,2]
  Pint.2=P.tab[2,2]
  Pint.3=P.tab[3,2]
 
  
 
  Bint.1=Biomass[1]
  Bint.2=Biomass[2]
  Bint.3=Biomass[3]
  Bint.4=Biomass[4]

 
  Nint.1=Biomass.tab[1,2]
  Nint.2=Biomass.tab[2,2]
  Nint.3=Biomass.tab[3,2]
  Nint.4=Biomass.tab[4,2]
 

  Bint.1=as.numeric(Bint.1)
  Bint.2=as.numeric(Bint.2)
  Bint.3=as.numeric(Bint.3)
  Bint.4=as.numeric(Bint.4)
 
 
  Nint.1=as.numeric(Nint.1)
  Nint.2=as.numeric(Nint.2)
  Nint.3=as.numeric(Nint.3)
  Nint.4=as.numeric(Nint.4)

 
  size.1=Bint.1/Nint.1
  size.2=Bint.2/Nint.2
  size.3=Bint.3/Nint.3
  size.4=Bint.4/Nint.4



  #Calculate N and B
  N = mean(as.numeric(Biomass.tab[,"Abundance"]))
  B = mean(Biomass)
  
  ## 5) Wrap up output
  Pars = c(N,B,P,igr,P_B, Pint.1, Pint.2, Pint.3,
           Bint.1, Bint.2, Bint.3, Bint.4,
           Nint.1, Nint.2, Nint.3, Nint.4,
           size.1, size.2, size.3, size.4)
  
  if (n==1) {
    MCsims = Pars
  } else {
    MCsims = rbind(MCsims,Pars)
  }
  
}

#clean up MCsims 
rownames(MCsims) = NULL
colnames(MCsims) =  c("N","B","P","igr","P_B", "Pint.1", "Pint.2", "Pint.3", "Bint.1","Bint.2", "Bint.3", "Bint.4", "Nint.1","Nint.2", "Nint.3", "Nint.4", "size.1", "size.2", "size.3", "size.4")

print(MCsims)

#print summary
print(summary(MCsims))

#plot summary
par(mfrow = c(length(Pars),1))
a= hist(MCsims[,"N"])
text(x= 0.95*(max(a$breaks)), y = 0.9*(max(a$counts)), paste("nSims = ", nBoot ))
hist(MCsims[,"B"])
hist(MCsims[,"P"])
hist(MCsims[,"igr"])
hist(MCsims[,"P_B"])


##write to a file
#write.csv(MCsims, file = file.name.out )


int1 <- D.elapse[2] - D.elapse[1]
int2 <- D.elapse[3] - D.elapse[2]
int3 <- D.elapse[4] - D.elapse[3]


Pd1 <- MCsims[,"Pint.1"]/int1
Pd2 <- MCsims[,"Pint.2"]/int2
Pd3 <- MCsims[,"Pint.3"]/int3

Pd = data.frame(Pd1 = Pd1, Pd2 = Pd2, Pd3 = Pd3)


Bavg1 <- (MCsims[, "Bint.1"] + MCsims[, "Bint.2"]) / 2
Bavg2 <- (MCsims[, "Bint.2"] + MCsims[, "Bint.3"]) / 2
Bavg3 <- (MCsims[, "Bint.3"] + MCsims[, "Bint.4"]) / 2


Pd.B1 <- Pd1/Bavg1
Pd.B2 <- Pd2/Bavg2
Pd.B3 <- Pd3/Bavg3

Pd.B = data.frame(Pd.B1 = Pd.B1, Pd.B2 = Pd.B2, Pd.B3 = Pd.B3)

#Function for 95% CIs
CI <- function(x)
{
quantile(x, c(0.025, 0.975))
}
#Calculate the column CIs
CI95.1 <- apply(MCsims, 2, CI)
CI95.2 <- apply(Pd, 2, CI)
CI95.3 <- apply(Pd.B, 2, CI)

#Calculate the column means
boot.mean.1 <- apply(MCsims, 2, mean)
boot.mean.2 <- apply(Pd, 2, mean)
boot.mean.3 <- apply(Pd.B, 2, mean)

total1 <- rbind(boot.mean.1, CI95.1)
total2 <- rbind(boot.mean.2, CI95.2)
total3 <- rbind(boot.mean.3, CI95.3)

q <- cbind(total1, total2)
r <- cbind(q, total3)

a <- t(r)
write.csv(a, file = file.name.out )

`````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````


#mean(MCsims[,1])
#quantile(MCsims[,1], c(0.025, 0.975))

P.tab
d1 <- rep(0, each=1000)
d2 <- rep(33, each=1000)
df <- data.frame(days=c(d1, d2))
a <- MCsims[,17]
b<-MCsims[,18]
c <- data.frame(size=c(a, b))
plot(c$size~df$days)
reg <- lm(c$size~df$days)
abline(reg)
confint(reg)

pdf("sample.pdf", 7, 5)
plot(ln.avg.mass ~ D.elapse)
abline(growth.lm)
dev.off()



int1 <- D.elapse[2] - D.elapse[1]
int2 <- D.elapse[3] - D.elapse[2]
int3 <- D.elapse[4] - D.elapse[3]
int4 <- D.elapse[5] - D.elapse[4]
int5 <- D.elapse[6] - D.elapse[5]
int6 <- D.elapse[7] - D.elapse[6]
int7 <- D.elapse[8] - D.elapse[7]
int8 <- D.elapse[9] - D.elapse[8]
int9 <- D.elapse[10] - D.elapse[9]
int10 <- D.elapse[11] - D.elapse[10]

Pd1 <- MCsims[,"Pint.1"]/int1
Pd2 <- MCsims[,"Pint.2"]/int2
Pd3 <- MCsims[,"Pint.3"]/int3
Pd4 <- MCsims[,"Pint.4"]/int4
Pd5 <- MCsims[,"Pint.5"]/int5
Pd6 <- MCsims[,"Pint.6"]/int6
Pd7 <- MCsims[,"Pint.7"]/int7
Pd8 <- MCsims[,"Pint.8"]/int8
Pd9 <- MCsims[,"Pint.9"]/int9
Pd10 <- MCsims[,"Pint.10"]/int10
Pd = data.frame(Pd1 = Pd1, Pd2 = Pd2, Pd3 = Pd3, Pd4 = Pd4, Pd5 = Pd5, Pd6 = Pd6, Pd7 = Pd7, Pd8 = Pd8, Pd9 = Pd9, Pd10 = Pd10)


Bavg1 <- (MCsims[, "Bint.1"] + MCsims[, "Bint.2"]) / 2
Bavg2 <- (MCsims[, "Bint.2"] + MCsims[, "Bint.3"]) / 2
Bavg3 <- (MCsims[, "Bint.3"] + MCsims[, "Bint.4"]) / 2
Bavg4 <- (MCsims[, "Bint.4"] + MCsims[, "Bint.5"]) / 2
Bavg5 <- (MCsims[, "Bint.5"] + MCsims[, "Bint.6"]) / 2
Bavg6 <- (MCsims[, "Bint.6"] + MCsims[, "Bint.7"]) / 2
Bavg7 <- (MCsims[, "Bint.7"] + MCsims[, "Bint.8"]) / 2
Bavg8 <- (MCsims[, "Bint.8"] + MCsims[, "Bint.9"]) / 2
Bavg9 <- (MCsims[, "Bint.9"] + MCsims[, "Bint.10"]) / 2
Bavg10 <- (MCsims[, "Bint.10"] + MCsims[, "Bint.11"]) / 2

Pd.B1 <- Pd1/Bavg1
Pd.B2 <- Pd2/Bavg2
Pd.B3 <- Pd3/Bavg3
Pd.B4 <- Pd4/Bavg4
Pd.B5 <- Pd5/Bavg5
Pd.B6 <- Pd6/Bavg6
Pd.B7 <- Pd7/Bavg7
Pd.B8 <- Pd8/Bavg8
Pd.B9 <- Pd9/Bavg9
Pd.B10 <- Pd10/Bavg10
Pd.B = data.frame(Pd.B1 = Pd.B1, Pd.B2 = Pd.B2, Pd.B3 = Pd.B3, Pd.B4 = Pd.B4, Pd.B5 = Pd.B5, Pd.B6 = Pd.B6, Pd.B7 = Pd.B7, Pd.B8 = Pd.B8, Pd.B9 = Pd.B9, Pd.B10 = Pd.B10)
dailyP <- merge(Pd, Pd.B)



# for loop to calculate the number of days in each interval
days = NULL
for(n in 1:length(D.elapse) - 1) {
days[n] = D.elapse[n + 1] - D.elapse[n]
}




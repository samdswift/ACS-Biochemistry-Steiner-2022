#Supplementary Material
#Samuel Swift

#Set Working Directory
setwd("/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model")

-------------------------------------------------------------------------------
#FIGURE 1

#A
#Clear Workspace
rm(list=ls())
#Set Constants
mu = 5
sd = 1
#Run 
Glue = NULL
for (x in c(seq(0,10,0.01))) { 
  y = (1/(sd*sqrt(2*pi)))*exp(-((x-mu)^2)/(2*sd^2))
  newrow = data.frame(x,y)
  Glue = rbind(Glue,newrow)
}
write.csv(Glue,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/1A.csv", row.names = FALSE)

#B
#Clear Workspace
rm(list=ls())
#Set Constants
R_0 = 1000
P_0 = 1000
K_d1 = 1000
K_d2 = 1000 #nM
n = 2
#Run
Latch = NULL
for (L in c(seq(1,9.9,0.01) %o% 10^(-3:2))) { 
  b = (R_0+P_0+K_d2)
  c = R_0*P_0
  cmax = (b - (sqrt((b)^2-4*c)))/(2*(min(P_0,R_0)))
  U = (cmax*K_d2)/(2*P_0-cmax)
  EC50 = 0.5*cmax+U+(K_d1*U)/(R_0-U-0.5*cmax)
  C = cmax*((L/EC50)^n)/(1+((L/EC50)^n))
  Fix = C+0.1
  newrow = data.frame(L,C,Fix)
  Latch = rbind(Latch,newrow)
}
write.csv(Latch,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/1B.csv", row.names = FALSE)

#C
#Clear Workspace
rm(list=ls())
#Set Constants
R_0 = 10
P_0 = 10
K_d1 = 1000
K_d2 = 2 #nM
n = 1
#Run1
Base = NULL
for (L in c(seq(1,9.9,0.01) %o% 10^(-1:2))) { 
  b = (R_0+P_0+K_d2)
  c = R_0*P_0
  cmax = (b - (sqrt((b)^2-4*c)))/2
  U = (cmax*K_d2)/(2*P_0-cmax)
  EC50 = 0.5*cmax+U+(K_d1*U)/(R_0-U-0.5*cmax)
  B = (cmax*((L/EC50)^n)/(1+((L/EC50)^n)))
  newrow = data.frame(L,B)
  Base = rbind(Base,newrow)
}
#Run2
R_0 = 100
R = NULL
for (L in c(seq(1,9.9,0.01) %o% 10^(-1:2))) { 
  b = (R_0+P_0+K_d2)
  c = R_0*P_0
  cmax = (b - (sqrt((b)^2-4*c)))/2
  U = (cmax*K_d2)/(2*P_0-cmax)
  EC50 = 0.5*cmax+U+(K_d1*U)/(R_0-U-0.5*cmax)
  C = (cmax*((L/EC50)^n)/(1+((L/EC50)^n)))
  newrow = data.frame(C)
  R = rbind(R,newrow)
}
#Run3
R_0 = 10
P_0 = 10
K_d1 = 100
K_d2 = 2 #nM
K1 = NULL
for (L in c(seq(1,9.9,0.01) %o% 10^(-1:2))) { 
  b = (R_0+P_0+K_d2)
  c = R_0*P_0
  cmax = (b - (sqrt((b)^2-4*c)))/2
  U = (cmax*K_d2)/(2*P_0-cmax)
  EC50 = 0.5*cmax+U+(K_d1*U)/(R_0-U-0.5*cmax)
  C = (cmax*((L/EC50)^n)/(1+((L/EC50)^n)))
  newrow = data.frame(C)
  K1 = rbind(K1,newrow)
}
#Run4
R_0 = 10
P_0 = 10
K_d1 = 1000
K_d2 = 0.2 #nM
K2 = NULL
for (L in c(seq(1,9.9,0.01) %o% 10^(-1:2))) { 
  b = (R_0+P_0+K_d2)
  c = R_0*P_0
  cmax = (b - (sqrt((b)^2-4*c)))/2
  U = (cmax*K_d2)/(2*P_0-cmax)
  EC50 = 0.5*cmax+U+(K_d1*U)/(R_0-U-0.5*cmax)
  C = (cmax*((L/EC50)^n)/(1+((L/EC50)^n)))
  newrow = data.frame(C)
  K2 = rbind(K2,newrow)
}
Graph = cbind(Base,R,K1,K2)
write.csv(Graph,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/1D.csv", row.names = FALSE)

#--------------------------------------------------------------------------------

#FIGURE 2

#A
#Clear Workspace
rm(list=ls())

#Cmax Graphs
C  = data.frame(c(seq(1,9.9,0.1) %o% 10^(-1:4)))
names(C)[1] <- paste("R_0")
P_0 = 100 #nM
for (K_d2 in c(1 %o% 10^(0:4))) { 
  DF = NULL
  for (R_0 in c(seq(1,9.9,0.1) %o% 10^(-1:4))) { 
    b = (R_0+P_0+K_d2)
    c = R_0*P_0
    cmax = (b - (sqrt((b)^2-4*c)))/(2*(min(P_0,R_0)))
    X = R_0/P_0
    newrow = data.frame(X,cmax)
    DF = rbind(DF,newrow)
  }
  names(DF)[2] <- paste("Kd2",K_d2,sep = "=")
  C = cbind(DF,C)
}

write.csv(C,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/2A.csv", row.names = FALSE)

#B&C
#Clear Workspace
rm(list=ls())

R_0 = 100
P_0 = 100
K_d1 = 1000
DFB = NULL
for (K_d2 in c(seq(1,9.9,0.1) %o% 10^(0:4))) { 
  b = (R_0+P_0+K_d2)
  c = R_0*P_0
  cmax = (b - (sqrt((b)^2-4*c)))/2
  U = (cmax*K_d2)/(2*P_0-cmax)
  EC50 = 0.5*cmax+U+(K_d1*U)/(R_0-U-0.5*cmax)
  newrow = data.frame(K_d2,EC50)
  DFB = rbind(DFB,newrow)
}

R_0 = 100
P_0 = 100
K_d2 = 2
DFC = NULL
for (K_d1 in c(seq(1,9.9,0.1) %o% 10^(0:6))) { 
  b = (R_0+P_0+K_d2)
  c = R_0*P_0
  cmax = (b - (sqrt((b)^2-4*c)))/2
  U = (cmax*K_d2)/(2*P_0-cmax)
  EC50 = 0.5*cmax+U+(K_d1*U)/(R_0-U-0.5*cmax)
  newrow = data.frame(K_d1,EC50)
  DFC = rbind(DFC,newrow)
}



write.csv(DFB,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/2B.csv", row.names = FALSE)
write.csv(DFC,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/2C.csv", row.names = FALSE)

#D
#Clear Workspace
rm(list=ls())

#Sets Parameters
K_d1= 1000 #nM
K_d2= 2 #nM
P_0 = 10 #nM
#Initializes Data Frame
DF = NULL
#Executes Equation
for (r in seq(1,50,0.01)) { 
  EC50 = 0.5*P_0+K_d2+K_d1/((P_0/K_d2)*(r-0.5)-1)
  newrow = data.frame(r,EC50,7)
  DF = rbind(DF,newrow)
}
#Creates a csv, readable by GraphPad Prism
write.csv(DF,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/2D.csv", row.names = FALSE)

#--------------------------------------------------------------------------------

#Figure 3

#B
#Clear Workspace
rm(list=ls())
Graph  = data.frame(c(seq(1,9.9,0.01) %o% 10^(0:4)))
names(Graph)[1] <- paste("L")
R_0 = 1
P_0 = 1
K_d1 = 1000
K_d2 = 2 
n = 1
for (K_doff in c(0,2,20,200,2000,20000)) { 
  DF = NULL
  for (L in c(seq(1,9.9,0.01) %o% 10^(0:4))) { 
    b = (R_0+P_0+K_d2)
    c = R_0*P_0
    cmax = (b - (sqrt((b)^2-4*c)))/2
    
    U = (cmax*K_d2)/(2*P_0-cmax)
    EC50 = 0.5*cmax+U+(K_d1*U)/(R_0-U-0.5*cmax)
    b = (R_0+P_0+K_doff)
    
    cmin = (b - (sqrt((b)^2-4*c)))/2
    C = cmin+(cmax-cmin)*((L/EC50)^n)/(1+((L/EC50)^n))
    newrow = data.frame(C)
    DF = rbind(DF,newrow)
  }
  names(DF)[1] <- paste("K_doff",K_doff,sep = "=")
  Graph = cbind(Graph,DF)
}
write.csv(Graph,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/3B.csv", row.names = FALSE)

#C

#Clear Workspace
rm(list=ls())
Graph  = data.frame(c(seq(1,9.9,0.01) %o% 10^(0:4)))
names(Graph)[1] <- paste("L")
R_0 = 100
P_0 = 100
K_d1 = 1000
K_d2 = 2 #nM
n = 1
for (K_doff in c(0,2,20,200,2000,20000)) { 
  DF = NULL
  for (L in c(seq(1,9.9,0.01) %o% 10^(0:4))) { 
    b = (R_0+P_0+K_d2)
    c = R_0*P_0
    cmax = (b - (sqrt((b)^2-4*c)))/2
    
    U = (cmax*K_d2)/(2*P_0-cmax)
    EC50 = 0.5*cmax+U+(K_d1*U)/(R_0-U-0.5*cmax)
    b = (R_0+P_0+K_doff)
    
    cmin = (b - (sqrt((b)^2-4*c)))/2
    C = cmin+(cmax-cmin)*((L/EC50)^n)/(1+((L/EC50)^n))
    newrow = data.frame(C)
    DF = rbind(DF,newrow)
  }
  names(DF)[1] <- paste("K_doff",K_doff,sep = "=")
  Graph = cbind(Graph,DF)
}
write.csv(Graph,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/3C.csv", row.names = FALSE)

#E
#Clear Workspace
rm(list=ls())
Graph  = data.frame(c(seq(1,9.9,0.1) %o% 10^(0:2)))
names(Graph)[1] <- paste("R")
P_0 = 1
K_d1 = 1000
K_d2 = 2 #nM
n = 1
for (K_doff in c(0,2,20,200,2000,20000)) { 
  DF = NULL
  for (R_0 in c(seq(1,9.9,0.1) %o% 10^(0:2))) { 
    b = (R_0+P_0+K_d2)
    c = R_0*P_0
    cmax = (b - (sqrt((b)^2-4*c)))/2
    U = (cmax*K_d2)/(2*P_0-cmax)
    EC50 = 0.5*cmax+U+(K_d1*U)/(R_0-U-0.5*cmax)
    b = (R_0+P_0+K_doff)
    cmin = (b - (sqrt((b)^2-4*c)))/2
    DR = cmax-cmin
    newrow = data.frame(EC50,DR)
    DF = rbind(DF,newrow)
  }
  names(DF)[1] <- paste("EC-K_doff",K_doff,sep = "=")
  names(DF)[2] <- paste("DR-K_doff",K_doff,sep = "=")
  Graph = cbind(Graph,DF)
}
write.csv(Graph,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/3E.csv", row.names = FALSE)

#F
#Clear Workspace
rm(list=ls())
Graph  = data.frame(c(seq(1,9.9,0.1) %o% 10^(0:2)))
names(Graph)[1] <- paste("R")
P_0 = 100
K_d1 = 1000
K_d2 = 2 #nM
n = 1
for (K_doff in c(0,2,20,200,2000,20000)) { 
  DF = NULL
  for (R_0 in c(seq(1,9.9,0.1) %o% 10^(2:4))) { 
    b = (R_0+P_0+K_d2)
    c = R_0*P_0
    cmax = (b - (sqrt((b)^2-4*c)))/2
    U = (cmax*K_d2)/(2*P_0-cmax)
    EC50 = 0.5*cmax+U+(K_d1*U)/(R_0-U-0.5*cmax)
    b = (R_0+P_0+K_doff)
    cmin = (b - (sqrt((b)^2-4*c)))/2
    DR = cmax-cmin
    newrow = data.frame(EC50,DR)
    DF = rbind(DF,newrow)
  }
  names(DF)[1] <- paste("EC-K_doff",K_doff,sep = "=")
  names(DF)[2] <- paste("DR-K_doff",K_doff,sep = "=")
  Graph = cbind(Graph,DF)
}
write.csv(Graph,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/3F.csv", row.names = FALSE)

#H
P_0 =100
R_0 = 100
K_d1 = 1000
K_d2 = 2
K_dD = 2

Po = 1
Kr = K_dD/P_0
Ro = R_0/P_0
K1 = K_d1/P_0
K2 = K_d2/P_0

DF = NULL
for (C in seq(0,0.86,0.001)) {
  p = Po - C  
  r = 0.25*(sqrt(Kr^2+8*Kr*(Ro+p-1)-4*K2*((1-p)/p)*(3*K2*((1-p)/p)+Kr))-(2*K2*((1-p)/p)+Kr))
  a = ((1-p)/p)*((K2*r/Kr)+(K1*K2/r)+K2+p+2*((K2^2)/Kr)*((1-p)/p))
  L = a*P_0
  C_adj = C*P_0
  newrow = data.frame(L,C_adj)
  DF = rbind(DF,newrow)
}
write.csv(DF,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/3H.csv", row.names = FALSE)

#I
#Clear Workspace
rm(list=ls())
#this just makes the x-axis values for the figure, the rest are created in MatLab
range = data.frame(seq(0.001,10,0.001))
write.csv(range,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/3I.csv", row.names = FALSE)

#--------------------------------------------------------------------------------

#Figure 4

#B
#Clear Workspace
rm(list=ls())
R_0 = 100
K_d1 = 1000
K_d2 = 2 #nM
n = 1
L = 99999999
Base = NULL
for (P_0 in c(seq(1,9.9,0.1) %o% 10^(-1:4))) { 
  cmax = (P_0-1/2*(sqrt((R_0+K_d2-P_0)^2+4*K_d2*P_0)-(R_0+K_d2-P_0)))
  U = (cmax*K_d2)/(2*P_0-cmax)
  EC50 = 0.5*cmax+U+(K_d1*U)/(R_0-U-0.5*cmax)
  B = (cmax*((L/EC50)^n)/(1+((L/EC50)^n)))
  newrow = data.frame(P_0,Inf,B)
  Base = rbind(Base,newrow)
}

POINTS = NULL
for (P_0 in c(c(1,2.5,5) %o% 10^(-1:4))) { 
  cmax = (P_0-1/2*(sqrt((R_0+K_d2-P_0)^2+4*K_d2*P_0)-(R_0+K_d2-P_0)))
  U = (cmax*K_d2)/(2*P_0-cmax)
  EC50 = 0.5*cmax+U+(K_d1*U)/(R_0-U-0.5*cmax)
  B = (cmax*((L/EC50)^n)/(1+((L/EC50)^n)))
  newrow = data.frame(P_0,B)
  POINTS = rbind(POINTS,newrow)
}
write.csv(Base,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/4B.csv", row.names = FALSE)
write.csv(POINTS,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/4BPoints.csv", row.names = FALSE)


#C
R_0 = 10
P_0 = 999999999
K_d1 = 1000
K_d2 = 2 #nM
n = 1
Base = NULL
for (L in c(seq(1,9.9,0.1) %o% 10^(-1:2))) { 
  cmax = (P_0-1/2*(sqrt((R_0+K_d2-P_0)^2+4*K_d2*P_0)-(R_0+K_d2-P_0)))
  U = (cmax*K_d2)/(2*P_0-cmax)
  EC50 = 0.5*cmax+U+(K_d1*U)/(R_0-U-0.5*cmax)
  B = (cmax*((L/EC50)^n)/(1+((L/EC50)^n)))*10
  newrow = data.frame(L,Inf,B)
  Base = rbind(Base,newrow)
}

POINTS = NULL
for (L in c(c(1,2.5,5) %o% 10^(-1:2))) { 
  cmax = (P_0-1/2*(sqrt((R_0+K_d2-P_0)^2+4*K_d2*P_0)-(R_0+K_d2-P_0)))
  U = (cmax*K_d2)/(2*P_0-cmax)
  EC50 = 0.5*cmax+U+(K_d1*U)/(R_0-U-0.5*cmax)
  B = (cmax*((L/EC50)^n)/(1+((L/EC50)^n)))*10
  newrow = data.frame(L,B)
  POINTS = rbind(POINTS,newrow)
}

write.csv(Base,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/4C.csv", row.names = FALSE)
write.csv(POINTS,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/4CPoints.csv", row.names = FALSE)


#E
#Clear Workspace
rm(list=ls())
R_0 = 0
K = 2020
Fit = NULL
for (P_0 in c(seq(4.1,3000,0.1))) { 
  EC50 = 0.5*R_0+K/P_0
  newrow = data.frame(P_0,EC50)
  Fit = rbind(Fit,newrow)
}
write.csv(Fit,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/4E.csv", row.names = FALSE)

#F
#Clear Workspace
rm(list=ls())
S = 1.921
K_doff = 14542
N = 0.044
Fit = NULL
for (P_0 in c(seq(4.1,3000,0.1))) { 
  A450 = N + S*P_0/(P_0+K_doff)
  newrow = data.frame(P_0,A450)
  Fit = rbind(Fit,newrow)
}
write.csv(Fit,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/4F.csv", row.names = FALSE)

#H
#Clear Workspace
rm(list=ls())
HAB = 50
K_d1 = 1000
K_d2 = 4
DF = NULL
for (PYR in c(seq(50,500,0.1))) { 
  x = PYR/HAB
  IC50 = 0.5*HAB+K_d2+(K_d1*K_d2)/(PYR-K_d2-0.5*HAB)
  newrow = data.frame(x,IC50,29)
  DF = rbind(DF,newrow)
}
write.csv(DF,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/4Hfit.csv", row.names = FALSE)

HAB = 50
K_d1 = 1000
K_d2 = 4
DF = NULL
for (PYR in c(c(HAB*1,HAB*1.5,HAB*2,HAB*3,HAB*4,HAB*5))) { 
  x = PYR/HAB
  IC50 = 0.5*HAB+K_d2+(K_d1*K_d2)/(PYR-K_d2-0.5*HAB)
  newrow = data.frame(x,Inf,Inf,IC50)
  DF = rbind(DF,newrow)
}
write.csv(DF,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/4Hpoints.csv", row.names = FALSE)

#I
#Clear Workspace
rm(list=ls())

Graph  = data.frame(c(seq(1,9.9,0.01) %o% 10^(0:2)))
names(Graph)[1] <- paste("Kd2")
HAB = 50
PYR = 250
for (IC50 in c(30,50,80,125,200)) { 
  DF = NULL
  for (K_d2 in c(seq(1,9.9,0.01) %o% 10^(0:2))) { 
    K_d1 = (IC50 - 0.5*HAB - K_d2)/(K_d2/(PYR-K_d2-0.5*HAB))
    newrow = data.frame(K_d1)
    DF = rbind(DF,newrow)
  }
  names(DF)[1] <- paste("IC50",IC50,sep = "=")
  Graph = cbind(Graph,DF)
}

write.csv(Graph,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/4I.csv", row.names = FALSE)

#--------------------------------------------------------------------------------
#Figure 5

#B
#Clear Workspace
rm(list=ls())
Graph  = data.frame(c(seq(1,9.9,0.01) %o% 10^(-2:4)))
names(Graph)[1] <- paste("L")
K_d1 = 1000
K_d2 = 2 
K_doff = 10000
n=1
for (R_0 in c(400,50,20)){
  for (P_0 in c(50,400,1000)) { 
    DF = NULL
    for (L in c(seq(1,9.9,0.01) %o% 10^(-2:4))) { 
      b = (R_0+P_0+K_d2)
      c = R_0*P_0
      cmax = (b - (sqrt((b)^2-4*c)))/2
      
      U = (cmax*K_d2)/(2*P_0-cmax)
      EC50 = 0.5*cmax+U+(K_d1*U)/(R_0-U-0.5*cmax)
      b = (R_0+P_0+K_doff)
      
      cmin = (b - (sqrt((b)^2-4*c)))/2
      C = cmin+(cmax-cmin)*((L/EC50)^n)/(1+((L/EC50)^n))
      A = C/R_0
      newrow = data.frame(A)
      DF = rbind(DF,newrow)
    }
    names(DF)[1] <- paste(R_0,P_0,sep = ",")
    Graph = cbind(Graph,DF)
  }
}
write.csv(Graph,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/5B.csv", row.names = FALSE)

#D
#Clear Workspace
rm(list=ls())
K_d1 = 836
K_d2 = 2
R = 50
P = 20

DF = NULL
for (C in c(seq(0,18.8,0.001))) { 
  U = (C*K_d2)/(P-C)
  L = C + U + (K_d1*U)/(R-U-C)
  S = C*7.584333333/18.8
  newrow = data.frame(L,S)
  DF = rbind(DF,newrow)
}

write.csv(DF,"/Users/samdswift/Desktop/CU/Whitehead Lab/CID Model/R_CSVs/5D.csv", row.names = FALSE)
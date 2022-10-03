
## =======================================================================
##  Model : PO4 adsorption column with ICS                               
## =======================================================================
## Lab application: kinetic adsorption in a flow though system
# package with ODE solution methods
require(deSolve)  
require(marelac)
require(ReacTran)
require(ggplot2)
#==========================#
#   Import measured data   #
#==========================#
# load data from a csv file and convert it to a data.frame
# set directory where the data is
directory<-getwd()
setwd(directory)
setwd("C:/Users/barcalap/OneDrive - Stichting Deltares/Documents/RTM Course")

#Column 1 continuous flow first 40 days, it is more consistent than the full length experiment
outB.ext <- read.csv(file="B1.csv", header=TRUE, sep=";")
outB.ext <- data.frame(outB.ext)
colnames(outB.ext) <- c("time", "Cout")
outB.ext$time[1]<-1
outB.ext$time<-as.numeric(outB.ext$time)
head(outB.ext, n=1)
#Column 1 continuous flow Long experiment, almost 5 months but the data is more noicy
outBL.ext <- read.csv(file="B2.csv", header=TRUE, sep = ";",dec = ".")
outBL.ext <- data.frame(outBL.ext)
colnames(outBL.ext) <- c("time", "Cout")
outBL.ext$time[1]<-1
outBL.ext$time<-as.numeric(outBL.ext$time)
head(outBL.ext, n=1)
#Column 1 Long experiment, direct problem in stanmod
outB1.ext <- read.table(file="B_long_1.txt", header=FALSE, sep = "",dec = ".")
outB1.ext <- data.frame(outB1.ext)
colnames(outB1.ext) <- c("time", "Cout")
outB1.ext$time[1]<-1
outB1.ext$time<-as.numeric(outB1.ext$time)
head(outB1.ext, n=5)
outB2.ext <- read.table(file="B_long_2.txt", header=FALSE, sep = "",dec = ".")
outB2.ext <- data.frame(outB2.ext)
colnames(outB2.ext) <- c("time", "Cout")
outB2.ext$time[1]<-1
outB2.ext$time<-as.numeric(outB2.ext$time)
head(outB2.ext, n=5)
outB3.ext <- read.table(file="B_long_3.txt", header=FALSE, sep = "",dec = ".")
outB3.ext <- data.frame(outB3.ext)
colnames(outB3.ext) <- c("time", "Cout")
outB3.ext$time[1]<-1
outB3.ext$time<-as.numeric(outB3.ext$time)
head(outB3.ext, n=5)
df_Stanmod_BL<-data.frame(outBL.ext$time,outBL.ext$Cout,outB1.ext$Cout[1:356],outB2.ext$Cout[1:356],outB3.ext$Cout[1:356])
#Column 3 stops several times accidentally because the tubes clogged
outC.ext <- read.csv(file="C.csv", header=TRUE, sep = ";",dec = ".")
outC.ext <- data.frame(outC.ext)
colnames(outC.ext) <- c("time", "VC","Cout")
outC.ext[1,1]<-1
outC.ext[,1]<-as.numeric(outC.ext[,1])
outC.ext[is.na(outC.ext)] = 0
head(outC.ext, n=1)
plot(outC.ext[,1],outC.ext[,3])
#Column 2 continuous flow different velocity  and porosity
outA.ext <- read.csv(file="A.csv", header=TRUE, sep = ";",dec = ",")
outA.ext <- data.frame(outA.ext)
colnames(outA.ext) <- c("time", "Cout")
outA.ext$time[1]<-1
outA.ext$time<-as.numeric(outA.ext$time)
head(outA.ext, n=1)
#Column 3 slower flow velocity 3.66 cm/h
outS.ext <- read.csv(file="C2.csv", header=TRUE, sep = ";",dec = ".")
outS.ext <- data.frame(outS.ext)
colnames(outS.ext) <- c("time", "Cout")
outS.ext[1,1]<-1
outS.ext[,1]<-as.numeric(outS.ext[,1])
outS.ext[is.na(outS.ext)] = 0
head(outS.ext, n=1)
plot(outS.ext[,1],outS.ext[,2])
#Column 3 desorption experiment
outD.ext <- read.csv(file="D.csv", header=TRUE,sep = ",",dec = ".")
outD.ext <- data.frame(outD.ext)
colnames(outD.ext) <- c("time", "Cout")
outD.ext$time[1]<-1
outD.ext$time<-as.numeric(outD.ext$time)
head(outD.ext, n=1)

#====================#
#      1D Grid       #
#====================#
# spatial domain
# Note that grid is the same for all of my columns, if the filter/column has a different length this needs to be changed!
Length    <- 30.8                              # [cm]
N         <- 154                               # number of grid cells, each cell is 0.2 cm, this can be changed
Grid      <- setup.grid.1D(L = Length, N = N)  # grid with equally sized boxes

#===============================#
#      Boundry conditions       #
#===============================#
#boundry conditions and state variables
#Note: boundry conditions are always the same expect for the desoprtion experiment
Paq.ini  <- rep(0, length = N)            # in mg/L
Pads_slow.ini  <- rep(0, length = N)      # in mg/g
state    <- c(Paq.ini, Pads_slow.ini)     # note the order! the solution is in the same order!
names    <- c("P(aq)", "Pads_slow")
nspec    <- length(names)

#===============================#
# Column 1 first 40 days        #
#===============================#
x        <- 30/690                # [g-ICS/g-solid] this is a new variable to calibrate taking into account the real amount of ICS in the column, f is larger and so alpha  
rho      <-  1770                 # [g/L(s)] bulk density of the filter 
porosity <- 0.55                  # [L water/L bulk] 

parms <- c(
  Dispers  =  30,                # [cm2/h]   Dispersion coefficient
  k_ads    =  3.565,             # [L water/g-solid]  Linear equilibrium constant
  alpha    =  1.56e-4,           # [1/h]  First order mass transfer coefficient
  f        =  0.0453,            # [-] fraction of sites in equilibrium
  Pinflow  =  1.67)              # [mg/L] Inflow concentration in the column


#Compare with 1
time_BB <- outB.ext$time
time_B <-c(1:694)                 # hours in the experiment, 1 hour is the time step
v        <-  4.92                 # [cm/h]     pore water velocity
v_adv    <- rep(v, length(time_B))
out_B <- ode.1D(y = state, parms = parms, func = P_column_k, times = time_B, 
                names = names, nspec = nspec, dimens = N) 
Measured <- outB.ext$Cout
Measured <- approx(Measured, method="linear", n=length(time_B))
Model <-out_B[,310] 

#Plot and R2 and MSE 
R2_B <-cor(Measured$y,Model)                                     #  RSquare for regression of observed vs predicted
MSE_B<-sqrt(sum((Measured$y-Model)^2)/length(time_B))            #  Mean squared error of observed vs predicted
df <-data.frame(time=time_BB,y=Measured$y[time_BB],z=Model[time_BB])
colors <- c("Measured" = "blue", "Modeld" = "red")
plot(out_B, which = c("Paq_out", "TotalPads_slow","TotalPads_fast","Paq_influx","Paq_outflux", "Pretention"))
plot2<-ggplot(df)+geom_point(aes(x=time, y, colour="Measured"))+
  geom_point(aes(x=time, z, color="Modeld"))+
  labs(title="Calibration column B v 4.92 cm/h", x ="Time (hours)", y = "mg/L",color = "Results")+
  scale_color_manual(values = colors)+ theme_bw()+xlim(0,650)
plot2
R2_B
MSE_B

#===============================#
# Column 1 full long experiment #
#===============================#
# we only need to change the times and length of v_adv vector
# changes
time_BBL <- outBL.ext$time
time_BL <-c(1:3047)#3047
v_adv    <- rep(v, length(time_BL))
#solution and plot
out_BL <- ode.1D(y = state, parms = parms, func = P_column_k, times = time_BL, 
                 names = names, nspec = nspec, dimens = N) 
Measured_1 <- outBL.ext$Cout
Measured <- approx(Measured_1, method="linear", n=length(time_BL))
Model <-out_BL[,310] 
R2_B <-cor(Measured$y,Model)                                                      #  RSquare for regression of observed vs predicted
MSE_B<-sqrt(sum((Measured$y-Model)^2)/length(time_BL))                            #  Mean squared error of observed vs predicted
df <-data.frame(time=time_BBL,y=Measured_1,z=Model[time_BBL])
colors <- c("Measured" = "blue", "Modeld" = "red")

plot2<-ggplot(df)+geom_point(aes(x=time, y, colour="Measured"), size=3, shape=1)+
  geom_line(aes(x=time, z, color="Modeld"), size=1)+
  scale_color_manual(values = c("red", "blue"))+
  theme_bw(base_size=22)+  theme(legend.title=element_blank(),legend.position="bottom")+ylim(0,1.4)+xlim(0,3047)+
  labs(title="Model validation long-tem column B 146 days", x ="t [Hours]", y = "P outflow concentrations [mg/L]",color = "Results")
plot2

R2_B
MSE_B


plot3<-ggplot(df_Stanmod_BL)+geom_point(aes(x=outBL.ext.time, outBL.ext.Cout, colour="Measured"), size=3, shape=1)+
  geom_line(aes(x=outBL.ext.time, outB1.ext.Cout.1.356., color="95% Confidence"), size=1,linetype = "dashed")+
  geom_line(aes(x=outBL.ext.time, outB2.ext.Cout.1.356., color="Modeld"), size=1)+
  geom_line(aes(x=outBL.ext.time, outB3.ext.Cout.1.356., color="95% Confidence"), size=1,linetype = "dashed")+
  scale_color_manual(values = c("gray", "red","blue"))+
  theme_bw(base_size=22)+  theme(legend.title=element_blank(),legend.position="bottom")+ylim(0,1.5)+xlim(0,3047)+
  labs(title="Model validation long-tem column B 146 days", x ="t [Hours]", y = "P outflow concentrations [mg/L]",color = "Results")
plot3
#===============================#
# Column 2 v 5.55 cm/h          #
#===============================#
# we need to adjust some parameters as porosity, D, velocity. x, k_ads, f, rho, and alpha are the same as is B
parms <- c(
  Dispers  =  30,                    # [cm2/h]   Dispersion coefficient obtained with stanmod
  k_ads    =  3.565,                 # [L water/g-solid]  Linear equiliArium constant
  alpha    =  1.56e-4,               # [1/h]  First order mass transfer coefficient
  f        =  0.0453,                # [-] fraction of sites in equiliArium
  Pinflow  =  1.55)                  # [mg/L]
porosity <- 0.49                     # [L water /L Bulk ]  

#Compare with A
time_AA <- outA.ext$time
time_A <-c(1:857)
v        <-  5.5 #                  # [cm/h]     pore water velocity
v_adv    <- rep(v, length(time_A))
out_A <- ode.1D(y = state, parms = parms, func = P_column_k, times = time_A, 
                names = names, nspec = nspec, dimens = N) 
Measured <- outA.ext$Cout
Measured <- approx(Measured, method="linear", n=length(time_A))
Model <-out_A[,310] 
R2_A <-cor(Measured$y,Model)                                                      #  RSquare for regression of oAserved vs predicted
MSE_A<-sqrt(sum((Measured$y-Model)^2)/length(time_A))                             #  Mean squared error of observed vs predicted
df <-data.frame(time=time_AA,y=Measured$y[time_AA],z=Model[time_AA])
colors <- c("Measured" = "blue", "Modeld" = "red")
plot2<-ggplot(df)+geom_point(aes(x=time, y, colour="Measured"))+
  geom_point(aes(x=time, z, color="Modeld"))+
  labs(title="Calibration column A v 5.55 cm/h", x ="Time (hours)", y = "mg/L",color = "Results")+
  scale_color_manual(values = colors)+ theme_bw()+xlim(0,857)
plot2
R2_A
MSE_A
#===============================#
# Column 3 slower flow 3.66 cm/h#
#===============================#
# we need to adjust some parameters as porosity, D, velocity. x, k_ads, f, rho, and alpha are the same as is B
#this experiment was after the stopflow and desorption, therefore there is higher initial P/Fe ratio at the start
parms <- c(
  Dispers  =  20,                  # [cm2/h]   Dispersion coefficient obtained with stanmod
  k_ads    =  3.565,               # [L water/g-solid]  Linear equiliArium constant
  alpha    =  1.56e-4,             # [1/h]  First order mass transfer coefficient
  f        =  0.0453,              # [-] fraction of sites in equiliArium
  Pinflow  =  1.60)                # [mg/L]
porosity <- 0.44                   # [L water /L Bulk ]  
x <- 30/690
#Compare with 3 slow flow
time_SS <- outS.ext$time
time_S <-c(1:609)
v        <-  3.66                  # [cm/h]     pore water velocity
v_adv    <- rep(v, length(time_S))
out_S <- ode.1D(y = state, parms = parms, func = P_column_k, times = time_S, 
                names = names, nspec = nspec, dimens = N) 
Measured <- outS.ext$Cout
Measured <- approx(Measured, method="linear", n=length(time_S))
Model <-out_S[,310] 
R2_S <-cor(Measured$y,Model)                                                      #  RSquare for regression of observed vs predicted
MSE_S<-sqrt(sum((Measured$y-Model)^2)/length(time_S))                             #  Mean squared error of observed vs predicted
df <-data.frame(time=time_SS,y=Measured$y[time_SS],z=Model[time_SS])
colors <- c("Measured" = "blue", "Modeld" = "red")
plot2<-ggplot(df)+geom_point(aes(x=time, y, colour="Measured"))+
  geom_point(aes(x=time, z, color="Modeld"))+
  labs(title="Calibration column 3 v 3.66 cm/h", x ="Time (hours)", y = "mg/L",color = "Results")+
  scale_color_manual(values = colors)+ theme_bw()+xlim(0,609)
plot2
R2_S
MSE_S
#===============================#
# Column 3 Stop flow            #
#===============================#
#change the porosity, velocity, D, all other parametrs are the same
time_CC <- outC.ext[,1]
time_C <-c(1:1481)

v_adv    <- outC.ext[,2]
v_adv    <- approx(v_adv, method="linear", n=length(time_C))
v_adv    <- v_adv$y
porosity <- 0.44
parms <- c(
  Dispers  =  30,             
  k_ads    =  3.565,   
  alpha    =  1.56e-4,   
  f        =  0.0453,  
  Pinflow  =  1.67)   
out_C <- ode.1D(y = state, parms = parms, func = P_column_k, times = time_C, 
                names = names, nspec = nspec, dimens = N) 
Measured <- outC.ext[,3]
Measured <- approx(Measured, method="linear", n=length(time_C))
Measured <- Measured$y
Measured <- ifelse (v_adv>0 ,Measured, 0) 
Model <-ifelse (v_adv>0 ,out_C[,310], 0) 
R2_C <-cor(Measured,Model)                                                      
MSE_C<-sqrt(sum((Measured-Model)^2)/length(time_C)) 
df <-data.frame(time=time_CC,y=Measured[time_CC],z=Model[time_CC])
colors <- c("Measured" = "blue", "Modeld" = "red")
plot3<-ggplot(df)+geom_point(aes(x=time,y , colour="Measured"))+
  geom_point(aes(x=time, z, color="Modeld"))+
  labs(title="Calibration column 3 v 6.09 cm/h", x ="Time (hours)", y = "mg/L",color = "Results")+
  scale_color_manual(values = colors)+ theme_bw()+ylim(0.01,1.2)
plot(out_C, which = c("Paq_out", "TotalPads_slow","TotalPads_fast","Paq_influx","Paq_outflux", "Pretention"))
plot3
R2_C
MSE_C

image(out_C, legend=TRUE, grid=Grid$x.mid, las=1, ylab="cm", xlab="t [Hours]")
text ("Adsorption in column 3\nP(aq) concnetration in the liquid phase (mg/L)\nPads_slow concentration in slow sites the solid (mg/g)")

summary(out_C)
plot(out_C, which = c("Paq_out", "TotalPads_slow","TotalPads_fast","Paq_influx","Paq_outflux", "Pretention"))
plot(out_C, which = c("TotalPads_slow","TotalPads_fast", "Padsfast", "Padsslow"))
fast<-cumsum(out_C[,316])
fast<-fast[length(fast)]
slow<-cumsum(out_C[,317])
slow<-slow[length(slow)]
fast/1481
slow/1481
# I think the issue is that at the start D is lower than at the end of the experiment.
#===============================#
# Column 3 Desrption            #
#===============================#
#  Note: we need to change the boundry conditions
### I modified alpha and f from the previous values to fit the desorption curve
Paq.ini  <- rep(0.84, length = N)             # in mg/L 0.84
Pads_slow.ini  <- rep(0.082, length = N)      # in mg/g solid 0.095 0.08235
state    <- c(Paq.ini, Pads_slow.ini)     
names    <- c("P(aq)", "Pads_slow")
nspec    <- length(names)

time_DD <- outD.ext[,1]
time_D <-c(1:713)
v        <-  6.09                            # [cm/h]     pore water velocity
v_adv    <- rep(v, length(time_D))
porosity <- 0.44
parms <- c(
  Dispers  =  30,               
  k_ads    =  3.565,    
  alpha    =  1.85e-4,                       # [1/h] Note alpha is larger!
  f        =  0.043, #0.125                  # [-] fraction of sites in equilibrium
  Pinflow  =  0.0)                           # Desorption [g/m3]
out_D <- ode.1D(y = state, parms = parms, func = P_column_k, times = time_D, 
                names = names, nspec = nspec, dimens = N) 
Measured <- outD.ext[,2]
Measured <- approx(Measured, method="linear", n=length(time_D))
Measured <- Measured$y
Model <-ifelse (v_adv>0 ,out_D[,310], 0) 
R2_C <-cor(Measured,Model)                                                      #  RSquare for regression of observed vs predicted
MSE_C<-sqrt(sum((Measured-Model)^2)/length(time_D)) 
df <-data.frame(time=time_DD,y=Measured[time_DD],z=Model[time_DD])
colors <- c("Measured" = "blue", "Modeld" = "red")
colors2 <- c("Measured" = 15, "Modeld" = 1)

plot3<-ggplot(df)+geom_point(aes(x=time,y , colour="Measured", shape="Measured"), size=3)+
  geom_point(aes(x=time, z, color="Modeld", shape="Modeld"), size=3, shape=1)+
  scale_color_manual(name="", labels=c("Measured","Model"),values = colors)+ theme_bw()+scale_shape_manual(name="", labels =c("Measured","Model"),values = colors2)+
  theme_bw(base_size=22)+  theme(legend.title=element_blank(),legend.position="bottom")+ylim(0.00,1.5)+
  labs(title="Desorption in column 3 6.09 cm/h", x ="t [Hours]", y = "P outflow concentrations [mg/L]")
plot3
R2_C
MSE_C

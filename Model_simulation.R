library(deSolve)
library(RColorBrewer)
library(reshape2)

 
 
region_birth_rate <- matrix(NA, nrow = 13, ncol = 10)
region_birth_rate[,1] <- c(0.01168, 0.01168, 0.01137, 
                           0.01099, 0.01065, 0.01056,
                           0.01037, 0.01025, 0.01023,
                           0.01013, 0.01010, 0.00990,
                           0.00970)

region_birth_rate[,2] <- c(0.01329, 0.01329, 0.01300,
                           0.01277, 0.01247, 0.01227,
                           0.01214, 0.01189, 0.01193,
                           0.01184, 0.01174, 0.01147,
                           0.01151)

region_birth_rate[,3] <- c(0.01302, 0.01302, 0.01277,
                           0.01242, 0.01211, 0.01199,
                           0.01189, 0.01171, 0.01180,
                           0.01170, 0.01157, 0.01133,
                           0.01118)

region_birth_rate[,4] <- c(0.01419, 0.01419, 0.01374,
                           0.01311, 0.01252, 0.01227,
                           0.01209, 0.01197, 0.01204,
                           0.01198, 0.01180, 0.01159,
                           0.01134)

region_birth_rate[,5] <- c(0.01349, 0.01349, 0.01320,
                           0.01281, 0.01237, 0.01221,
                           0.01213, 0.01206, 0.01214, 
                           0.01208, 0.01196, 0.01169, 
                           0.01148)


region_birth_rate[,6] <- c(0.01638, 0.01638, 0.01599,
                           0.01556, 0.01479, 0.01428,
                           0.01426, 0.01421, 0.01436, 
                           0.01426, 0.01389, 0.01320,
                           0.01292)

region_birth_rate[,7] <- c(0.01425, 0.01425, 0.01403,
                           0.01371, 0.01327, 0.01303, 
                           0.01304, 0.01290, 0.01297,
                           0.01287, 0.01273, 0.01233, 
                           0.01224)


region_birth_rate[,8] <- c(0.01618, 0.01618, 0.01589, 
                           0.01527, 0.01469, 0.01427, 
                           0.01421, 0.01401, 0.01398, 
                           0.01381, 0.01354, 0.01297, 
                           0.01241)

region_birth_rate[,9] <- c(0.01575, 0.01575, 0.01515, 
                           0.01429, 0.01368, 0.01330,
                           0.01321, 0.01290, 0.01294, 
                           0.01256, 0.01242, 0.01190, 
                           0.01149)

region_birth_rate[,10] <- c(0.01410, 0.01410, 0.01402,
                            0.01347, 0.01299, 0.01275,
                            0.01270, 0.01250, 0.01257, 
                            0.01242, 0.01229, 0.01170,
                            0.01124)


#nn <-  1 # change nn to change regions nn = 1:10

Pop1 <- readRDS('pop1.rds') # initial population, by age group (21 age groups considered here)
B0 <- readRDS('Birth_rate.rds') # (34 years = 408 month / 12)
B_mean <- mean(B0[,1])
B_part1 <- rep(B_mean, 41*52) # 2132 weeks = 41 years as burn-in period

B_part2 <- rep(c(B0[369:388,1], region_birth_rate[,nn]), each = 52) # 


B_all <- c(B_part1, B_part2) # total 
B <- matrix(0, nrow =  length(B_all), ncol = 21)
B[,1] <- B_all


c2 <- readRDS('c2.rds') # contact matrix 


WidthAgeClassMonth <-  c(rep(1,times=12), rep(12,times=4),  60, 120, 240, 240, 240)  #Aging rate=1/width age class (months) Vector of long N_age


N_ages <- length(Pop1) 
agenames <- paste0('Agegrp', 1:N_ages) #Could replace this with vector of actual age names


 
p <- sum(Pop1)  # Total population at each time, a vector of length T; here p is total population of all age groups @ t = 0;


## Initialize the compartments (States) 



StateNames <- c("M", 
                "Xss", "Xi1s", "Xs1s", "Xi2s", "Xs2s", "Xi3s", "Xs3s", "Xi4s")
                

States <- array(NA, dim=c(N_ages, length(StateNames) )) #  N age groups x N parameters 
dimnames(States)[[1]] <- agenames
dimnames(States)[[2]] <- StateNames

yinit.matrix <- array(NA, dim=c(N_ages, length(StateNames) ))

dimnames(yinit.matrix)[[1]] <- agenames
dimnames(yinit.matrix)[[2]] <- StateNames

yinit.matrix[,c("Xs1s", "Xi2s", "Xs2s", "Xi3s", "Xs3s", "Xi4s")]  <-  0 # setting initial conditions

yinit.matrix[,'M']  <-  c(Pop1[1:6], rep(0,N_ages-6))
yinit.matrix[,'Xss'] <-  c(rep(0,6),Pop1[7:N_ages]-rep(N_ages-6)) 
yinit.matrix[,c('Xi1s')]  <-  c(rep(0,6), rep(1,N_ages-6))  #initializes with 1 infected person per age group 
 

yinit.vector <- as.vector(yinit.matrix) #Vectorize the ynit matrix

# Create array that has the labels by age, State and use this to name the yinit.vector
name.array <- array(NA, dim=dim(yinit.matrix)) # dim = 21 x 25 (21 age groups x 25 model compartments) 
for(i in 1:dim(name.array)[1]){ # for 1:21 age groups 
  for(j in 1:dim(name.array)[2]){ # for 1:25 model compartnments (stages)
    name.array[i,j] <- paste(dimnames(yinit.matrix)[[1]][i],dimnames(yinit.matrix)[[2]][j])
  }
}

name.vector <- as.vector(name.array)
names(yinit.vector) <- name.vector



start_time = 1 # start date (years)
tmax = nrow(B)
# end_time = 25 # end date (years)
my_times <- seq(start_time, tmax, by = 1) # gives a sequence from start to end
# in increments of 1
#########################################
#Matrix of dimension N_ages x N_ages
#beta = beta 
#########################################

#########################################
#Seasonal components--
#???Should prob be estimated?
# Amp = 0.2 #Seasonal amplitude
# phi = 3.32749 #Seasonal phase shift (weeks;[-26,26] 0=peak @ Jan 1)

 

#########################################

#########################################
#Relative infectiousness for 2nd and subsequent infections
rho1 = 0.75
rho2 = 0.51
#########################################

#########################################
# duration of infectiousness (months)
#Duration in days
dur.days1 <- 10 #days
dur.days2 <- 7 #days
dur.days3 <- 5 #days
 
###########################################  

###########################################
# 1/duration of maternal immunity (DAYS)
DurationMatImmunityDays = 112
###########################################

 
############################################
#????
#What is time(-1) unit?
#um= -0.0001833333  #?? is this right? only die from last age class--this maybe helps with this "We assumed deaths occurred from the last age class and adjusted the net rate of immigration/emigration and death from other age groups in order to produce a rate of population growth and age structure similar to that of the US." 
#?????????
#From Giny's code: um=log(0.993)/52 #net rate of crude deaths (+) and immigration (-) from all age groups (per week): can adjust this to approximate population growth in state# I calibrated this parameter so we can reproduce the population growth
um = -0.0002227
#############################################

#############################################
#Birth rate (births/person/YEAR
#Matrix: T rows, N_ages columns; columns 2:N_ages all 0s
PerCapitaBirthsYear=B 
#############################################


#############################################
#Relaive risk of infection following 1st, 2nd, 3rd+ infections
sigma1 = 0.76
sigma2 = 0.6
sigma3 = 0.4
#############################################
 
##To ESTIMTE from DATA?: Baseline transmission rate
#############################################

baseline.txn.rate <- 8.91
q = 1
c2 = c2


my_parmset <-list(PerCapitaBirthsYear=PerCapitaBirthsYear, 
            WidthAgeClassMonth=WidthAgeClassMonth,
            DurationMatImmunityDays = DurationMatImmunityDays,
            um=um,
            rho1=rho1,
            rho2=rho2,
            dur.days1=dur.days1,
            dur.days2=dur.days2,
            dur.days3=dur.days3,
            yinit.matrix=yinit.matrix,
            baseline.txn.rate = baseline.txn.rate,
            q=q,
            contact=c2,
            sigma1=sigma1,
            sigma2=sigma2,
            sigma3= sigma3,
            time.step = 'week')

 
 
                         
  

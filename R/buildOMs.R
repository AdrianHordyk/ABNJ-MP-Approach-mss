
library(MSEtool)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

nsim <- 300
proyears <- 30
root <- "../../tRFMOs/Data"
source('R/functions.r')
# ----Build the Operating Models ----

# ---- Data-Rich - Bigeye Tuna ----

SSdir1 <- file.path(root, 'BET_ICCAT', 'grid/7-steepness0.7_MRef_sigmaR0.4_LengthLambda0.1_iter7')
SSdir2 <- file.path(root, 'BET_ICCAT', 'grid/8-steepness0.8_MRef_sigmaR0.4_LengthLambda0.1_iter8')
SSdir3 <- file.path(root, 'BET_ICCAT', 'grid/9-steepness0.9_MRef_sigmaR0.4_LengthLambda0.1_iter9')

OM1 <- SS2OM(SSdir1, nsim/3, proyears)
OM2 <- SS2OM(SSdir2, nsim/3, proyears, seed=101)
OM3 <- SS2OM(SSdir3, nsim/3, proyears, seed=1001)

OMlist <- list(OM1, OM2, OM3)

OM_BET <- CombineOM(OMlist, "Bigeye Tuna")

OM_BET@cpars$D %>% unique()

# ---- Data-Moderate - Blue Shark ----
OM_BSH <- testOM
OM_BSH <- tinyErr(OM_BSH)

OM_BSH@Name <- "Blue shark"

# file:///C:/Users/Adrian/Dropbox/Projects_2019/tRFMOs/Data/BSH_IO/IOTC-2017-WPEB13-33_0.pdf
OM_BSH@maxage <- 30


OM_BSH@M <- c(0.3,0.22,0.18,0.156,0.14,0.128,0.12,0.114,0.109,0.105,0.101,0.099,
              0.096,0.095,0.093,0.092,0.09,0.089,0.089,0.088,0.087,0.087,0.086,
              0.086,0.085,0.085,0.085,0.085,0.084,0.084)
OM_BSH@M2 <- c(0.309,0.233,0.194,0.171,0.155,0.144,0.135,0.129,0.124,0.12,0.117,
               0.114,0.112,0.11,0.109,0.107,0.106,0.105,0.105,0.104,0.103,0.103,
               0.103,0.102,0.102,0.102,0.101,0.101,0.101,0.101)
# file:///C:/Users/Adrian/Dropbox/Projects_2019/tRFMOs/Data/BSH_IO/Blanco-Parraetal2008Blueshark.pdf
OM_BSH@Linf <- c(237.5, 199.85) # from above paper
OM_BSH@LenCV <- c(0.25, 0.25) # from assessment doc
OM_BSH@K <- c(0.1, 0.15)
OM_BSH@t0 <- c(-2.44, -2.15)
OM_BSH@a <-3.661e-05 # from SS
OM_BSH@b <- 2.901 # from SS
OM_BSH@SRrel <- 1
OM_BSH@h <- c(0.7, 0.8) # assumed 
OM_BSH@L50 <- c(145, 150) # from assessment doc
OM_BSH@L50_95 <- c(5,5) # assumed
OM_BSH@Perr <- c(0.6, 0.7) # assumed

OM_BSH@Size_area_1 <- OM_BSH@Prob_staying <- OM_BSH@Frac_area_1 <- c(0.5,0.5)

OM_BSH <- Replace(OM_BSH, Perfect_Imp)
OM_BSH <- Replace(OM_BSH, Generic_Obs)

OM_BSH@CurrentYr <- 2015
OM_BSH@nyears <- length(1950:OM_BSH@CurrentYr)
OM_BSH@nsim <- nsim * 2
OM_BSH@proyears <- proyears
OM_BSH@interval <- 1 

Chist <- read.csv(file.path(root, "BSH_IO/Catch.csv"))
Chist <- Chist[,2]

CALdata <- read.csv(file.path(root, "BSH_IO/CAL.csv"), header=FALSE, 
                    stringsAsFactors = FALSE)
length_bin <- CALdata[1, 2:(length(CALdata)-1)] + 2.5 
length_bin <- as.numeric(length_bin)
CAL <- matrix(NA, nrow=OM_BSH@nyears, ncol=length(length_bin))
yrs <- CALdata[2:nrow(CALdata),1] %>% as.numeric()
ind <- match(yrs, 1950:OM_BSH@CurrentYr)
CAL[ind, ] <- CALdata[2:nrow(CALdata), 2:(ncol(CAL)+1)] %>% as.matrix()

OM_scope <- SRA_scope(OM_BSH, Chist, CAL=CAL, length_bin=length_bin, cores=10)

OM_BSH2 <- Sub_cpars(OM_scope@OM, sims=which(OM_scope@conv)[1:nsim])
OM_BSH2@nsim
OM_BSH <- OM_BSH2

# ---- Data-Poor - Yellowfin - WCPO ----
OM_YFT <- new("OM", nsim=nsim)
OM_YFT@Name <- 'Yellowfin tuna'  
OM_YFT@Common_Name <- 'Yellowfin tuna'
OM_YFT@Species <- 'Thunnus albacares'

OM_YFT@nsim <- nsim 
OM_YFT@proyears <- proyears    
OM_YFT@interval <- 1 
OM_YFT@pstar <- 0.5
OM_YFT@maxF <- 3        
OM_YFT@reps <- 1    
OM_YFT@seed <- 1      
 

OM_YFT <- LH2OM(OM_YFT)
OM_YFT@maxage <- max(-log(0.01)/OM_YFT@cpars$M) %>% ceiling()
OM_YFT@t0 <- c(0,0)
OM_YFT@Mexp <- c(-0.315, -0.261)
OM_YFT@LenCV <- c(0.1, 0.2)

OM_YFT@R0 <- 10000
OM_YFT@SRrel <- 1
OM_YFT@h <- c(0.65, 0.95)
OM_YFT@Perr <- c(0.5, 0.9)
OM_YFT@AC <- c(0.1, 0.9)

OM_YFT@Msd <- OM_YFT@Mgrad <- OM_YFT@Ksd <- OM_YFT@Kgrad <- OM_YFT@Linfsd <- 
  OM_YFT@Linfgrad <- c(0,0)  


vals <- 1+ runif(nsim, 0.05, 0.1)
L95 <- OM_YFT@cpars$L50 * vals
OM_YFT@cpars$L50_95 <- L95 - OM_YFT@cpars$L50
OM_YFT@L50_95 <- c(0,0)

OM_YFT@D <- c(0.2, 0.5) 


# Mean lengths at age from assessment
lengths <- c(25.1198,40.6064,48.6428,58.3577,72.4081,86.2795,97.2252,106.0257,
         113.3911,118.9618,123.8023,128.0084,131.6633,134.8391,137.5986,
         139.9965,142.0801,143.8906,145.4638,146.8309,148.0187,149.0509,
         149.9478,150.7271,151.4043,151.9927,152.5040,152.9483)

# Mean weights at age
weight <- c(0.3573,1.3955,2.3527,3.9914,7.4782,12.4743,17.6927,22.8080,27.7770,
            31.9842,35.9706,39.6968,43.1426,46.3018,49.1788,51.7851,54.1367,
            56.2527,58.1534,59.8597,61.3923,62.7710,64.0148,65.1416,66.1678,
            67.1091,67.9797,68.7929)

mod <- lm(log(weight)~log(lengths))

OM_YFT@a <- as.numeric(exp(mod$coefficients[1])  )
OM_YFT@b  <- as.numeric(mod$coefficients[2])

OM_YFT@Size_area_1  <- c(0.5,0.5) 
OM_YFT@Frac_area_1  <- c(0.5,0.5) 
OM_YFT@Prob_staying <- c(0.5,0.5) 

OM_YFT@Fdisc <- c(0,0)

OM_YFT@nyears <- length(1950:2019)      
OM_YFT@EffYears <- c(1950, 1970, 1990, 2010, 2019)       
OM_YFT@EffLower <- c(0.00, 0.05, 0.25, 0.75, 0.90)    
OM_YFT@EffUpper <- c(0.0, 0.1, 0.5, 0.9, 1.0)   
OM_YFT@CurrentYr <- 2019    

OM_YFT@Esd <- c(0.05, 0.1)

OM_YFT@qinc <- c(0,0)
OM_YFT@qcv <- c(0,0)          

OM_YFT@L5 <- c(30, 40)          
OM_YFT@LFS <- c(110, 130)          
OM_YFT@Vmaxlen <- c(0.4,1)  
OM_YFT@isRel <- FALSE
OM_YFT@Spat_targ <- c(0,0)
   
quantile(OM_YFT@cpars$Linf, c(0.05, 0.5, 0.95))
quantile(OM_YFT@cpars$K, c(0.05, 0.5, 0.95))
quantile(OM_YFT@cpars$L95, c(0.05, 0.5, 0.95))

# ---- Observation Parameters ----
BET_Obs <- BSH_Obs <- YFT_Obs <- Generic_Obs

BET_Obs@Cobs <- c(0.1,0.1)
BET_Obs@Cbiascv <- 0.05
BET_Obs@Iobs <- c(0.2,0.2)
BET_Obs@beta <- c(0.8,1)
BET_Obs@CAL_nsamp <- c(300,500)
BET_Obs@CAL_ESS <- c(150, 200)
BET_Obs@CAA_nsamp <- c(100, 200)
BET_Obs@CAA_ESS <- c(25,50)
BET_Obs@Dobs <- c(0.1,0.1)
BET_Obs@Dbiascv <- 0
BET_Obs@hbiascv <- 0.2
BET_Obs@Mbiascv <- 0.1
BET_Obs@Linfbiascv <- 
  BET_Obs@Kbiascv <- BET_Obs@t0biascv <- 0.05

BSH_Obs@Cobs <- c(0.2,0.2)
BSH_Obs@Cbiascv <- 0.1
BSH_Obs@Iobs <- c(0.3,0.3)
BSH_Obs@beta <- c(0.6,1.1)
BSH_Obs@CAL_nsamp <- c(300,500)
BSH_Obs@CAL_ESS <- c(150, 200)
BSH_Obs@CAA_nsamp <- c(50, 100)
BSH_Obs@CAA_ESS <- c(10,25)
BSH_Obs@Dobs <- c(0.3,0.3)
BSH_Obs@Dbiascv <- 0
BSH_Obs@hbiascv <- 0.4
BSH_Obs@Mbiascv <- 0.2
BSH_Obs@Linfbiascv <- 
  BSH_Obs@Kbiascv <- BSH_Obs@t0biascv <- 0.05


YFT_Obs@Cobs <- c(0.3,0.3)
YFT_Obs@Cbiascv <- 0.1
YFT_Obs@Iobs <- c(0.4,0.4)
YFT_Obs@beta <- c(0.5,2)
YFT_Obs@CAL_nsamp <- c(300,500)
YFT_Obs@CAL_ESS <- c(150, 200)
YFT_Obs@CAA_nsamp <- c(25, 50)
YFT_Obs@CAA_ESS <- c(5,10)
YFT_Obs@Dobs <- c(0.5,0.5)
YFT_Obs@Dbiascv <- 0
YFT_Obs@hbiascv <- 0.6
YFT_Obs@Mbiascv <- 0.4
YFT_Obs@Linfbiascv <- 
  YFT_Obs@Kbiascv <- YFT_Obs@t0biascv <- 0.05

OM_BET <- Replace(OM_BET, BET_Obs)
OM_BSH <- Replace(OM_BSH, BSH_Obs)
OM_YFT <- Replace(OM_YFT, YFT_Obs)

# ---- Implementation Parameters ----
OM_BET <- Replace(OM_BET, Perfect_Imp)
OM_BSH <- Replace(OM_BSH, Perfect_Imp)
OM_YFT <- Replace(OM_YFT, Perfect_Imp)
 
# ---- Save OMs ----
OM_BSH@cpars$Data <- NULL # remove data - only use simulated

saveRDS(OM_BET, 'OMs/OM_BET.rdata')
saveRDS(OM_BSH, 'OMs/OM_BSH.rdata')
saveRDS(OM_YFT, 'OMs/OM_YFT.rdata')

# --- Load OMs ---- 
OM_BET <- readRDS('OMs/OM_BET.rdata')
OM_BSH <- readRDS('OMs/OM_BSH.rdata')
OM_YFT <- readRDS('OMs/OM_YFT.rdata')


# ---- Make plots ----
Hist_BET <- runMSE(OM_BET, Hist=TRUE)
Hist_BSH <- runMSE(OM_BSH, Hist=TRUE)
Hist_YFT <- runMSE(OM_YFT, Hist=TRUE)


num <- 0
ylab <- c("Length (cm)", "Natural Mortality", "Maturity", "Selectivity")
plist <- list() 
num3 <- 0
for (OM in c(Hist_BET, Hist_BSH, Hist_YFT)) {
  num <- num + 1
  lendat <- OM@AtAge$Length[,,1]
  Mdat <- OM@AtAge$N.Mortality[,,1]
  matdat <- OM@AtAge$Maturity[,,1]
  seldat <- OM@AtAge$Select[,,1]
  
  datalist <- list(lendat, Mdat, matdat, seldat)
  num2 <- num
  dnum <- 0
  for (data in datalist) {
    dnum <- dnum+1
    df <- tidyr::pivot_longer(as.data.frame(data), 1:ncol(data))
    df$age <- rep(1:ncol(data))
    df$sim <- rep(1:nsim, each=ncol(data))
    
    p <- ggplot(df, aes(x=age, y=value, group=as.factor(sim), alpha=0.6)) +
      geom_line(size=1.1) +
      theme_classic() +
      expand_limits(y=0, x=1) +
      theme(legend.position = 'none') 
    p <- p + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
    if (dnum == length(datalist)) {
      p <- p+labs(x="Age (year)")
    } else {
      p <- p+labs(x="")
      p <- p + theme(axis.title.x = element_blank(),
                     axis.text.x=element_blank())
    }
    if (num == 1) {
      p <- p + labs(y=ylab[dnum])
    } else {
      p <- p + labs(y="")
      # p <- p + theme(axis.title.y = element_blank(),
                     # axis.text.y=element_blank())
    }
    num3 <- num3 + 1
    p <- p + labs(tag=paste0(letters[num3], ")"))
 
    plist[[num2]] <- p
    num2 <- num2 + 3
  }
}

pout <- cowplot::plot_grid(plotlist=plist, ncol=3, nrow=4,
                           align="v", axis="lb", 
                           rel_heights=c(1,1,1,1.2))

ggsave('Figures/Figure_1.png', pout, width=7, height=8)


plist <- list() 
num <- 0
for (hist in c(Hist_BET, Hist_BSH, Hist_YFT)) {
  num <- num +1 
  nyears <- dim(hist@TSdata$VB)[2]
  Years <- (hist@Misc$CurrentYr - nyears + 1):hist@Misc$CurrentYr
  
  data <- hist@TSdata$Find
  maxs <- apply(data, 1, max)
  data <- data/maxs
  df <- tidyr::pivot_longer(as.data.frame(data), 1:ncol(data))
  df$year <- Years
  df$sim <- rep(1:nsim, each=ncol(data))
  
  p <- ggplot(df, aes(x=year, y=value, group=as.factor(sim), alpha=0.6)) +
    geom_line(size=0.5) +
    expand_limits(x=c(1950, 2020)) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.text.y=element_blank(),
          plot.tag = element_text(size=8),
          axis.title=element_text(size=8)) +
    labs(x="Year", y="", tag=paste0(letters[num], ')')) 
  if (num == 1) {
    p <- p + labs(y="Relative fishing mortality") +
      theme(axis.text.y = element_text())
           
  }
  
  plist[[num]] <- p
  
}

pout <- cowplot::plot_grid(plotlist=plist, ncol=3, nrow=1)
ggsave('Figures/Figure_2.png', pout, width=6, height=2)








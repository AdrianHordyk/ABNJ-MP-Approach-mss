
args <- commandArgs(TRUE)

if (length(args)<1) {
  i <- 3
} else {
  i <- as.numeric(args[1])
}

# devtools::install_github("DLMtool/DLMtool")
# devtools::install_github("tcarruth/MSEtool")
library(MSEtool)
library(dplyr)

OMfiles <- list.files("OMs")
OMfiles <- OMfiles[!grepl('old', OMfiles)] %>% sort()

OM <- readRDS(file.path('OMs', OMfiles[i]))

DD_4010 <- MSEtool::make_MP(DD_TMB, HCR40_10)
# DD_MSY <- MSEtool::make_MP(DD_TMB, HCR_MSY)
# MPs <- c("AvC", "MCD", "HDAAC", "Ltarget1", "Ltarget4", "Itarget1", "Itarget4",
#          "ICI", "SP_4010", "DD_4010", "SCA_4010")

MPs <- c("AvC", "HDAAC", "Ltarget1", "Ltarget4", "Itarget1", "Itarget4",
         "SP_4010", "DD_4010", "SCA_4010")



setup()
MSE <- try(runMSE(OM, MPs=MPs, parallel = TRUE, PPD=TRUE))

nm1 <- strsplit(OMfiles[i], "_")[[1]][2]

# nm1 <- paste0('NOHCR_', nm1)

Name <- paste0('MSEs/', nm1)

saveRDS(MSE, Name)

# options(warn=2)
# OM2 <- OM
# # OM2@nsim <-5 
# tryMSE <- runMSE(OM2, MPs='SCA_4010')


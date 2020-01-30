library(dplyr)
library(MSEtool)
library(ggplot2)
library(tidyr)
library(data.table)
library(scales)
library(cowplot)
library(ggrepel)

split.along.dim <- function(a, n)
  setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
                  array, dim = dim(a)[-n], dimnames(a)[-n]),
           dimnames(a)[[n]])

DD_4010 <- MSEtool::make_MP(DD_TMB, HCR40_10)


MSEs <- list.files("MSEs")
MSEs <- MSEs[!grepl('11_MPs', MSEs)]
# MSEs <- MSEs[c(1,2,6)]
# MSEs <- MSEs[3:5] # no 4010 assessments

MSENames <- strsplit(MSEs, ".rdata") %>% unlist()

# --- Kobe Plots ----
makeKplot <- function(MSE, mm=1, i=1, MSENames, yrs=11:30, 
                      values=c("green", "red", "yellow"),
                      TextSize=4,
                      TitleSize=14,
                      MSEsize=16,
                      MSEloc=-0.5,
                      minAlpha=0.4,
                      text.size=4,
                      title.size=4) {
  Fs <- MSE@F_FMSY[,mm,yrs] %>% as.vector()
  Fs[Fs>2] <- 1.99
  Bs <- MSE@B_BMSY[,mm,yrs] %>% as.vector()
  Bs[Bs>2] <- 1.99
  x <- cbind(Bs, Fs)
  ab <- matrix(c(0,0,2,2),2,2)
  nbin <- 40
  rng <- range(ab)
  bins <- seq(from=min(rng), to=max(rng), length.out = nbin)
  dat <- ash::bin2(x,ab,c(nbin, nbin))
  dat <- as.data.frame(dat$nc)
  rownames(dat) <- bins
  colnames(dat) <- bins
  DF <- melt(setDT(dat, keep.rownames = TRUE), "rn")
  colnames(DF)<- c('x', 'y', 'val')
  DF$x <- as.numeric(DF$x)
  DF$y <- as.character(DF$y)
  DF$y <- as.numeric(DF$y)
  
  # calc stats 
  TL <- DF %>% filter(x<1 & y>=1)
  TR <- DF %>% filter(x>=1 & y>=1)
  BL <- DF %>% filter(x<1 & y<1)
  BR <- DF %>% filter(x>=1 & y<1)
  
  TL <- (sum(TL$val)/(MSE@nsim * length(yrs)) * 100) %>% round(2)
  TR <- (sum(TR$val)/(MSE@nsim * length(yrs)) * 100) %>% round(2)
  BL <- (sum(BL$val)/(MSE@nsim * length(yrs)) * 100) %>% round(2)
  BR <- (sum(BR$val)/(MSE@nsim * length(yrs)) * 100) %>% round(2)
  
  xloc <- 0.4
  yloc <- 0.1
  maxV <- max(ab)
  textDF <- data.frame(x=c(xloc, maxV-xloc, xloc, maxV-xloc),
                       y=c(maxV-yloc, maxV-yloc, yloc, yloc),
                       label=paste0(c(TL, TR, BL, BR), "%"),
                       color=1, val=2,
                       stringsAsFactors = FALSE)
  
  DF$color <- 'yellow'
  DF$color[DF$x > 1 & DF$y <=1] <- 'green'
  DF$color[DF$x<=1 & DF$y>1] <- 'red'
  DF$val <- rescale(DF$val, to=c(minAlpha,1), from=range(DF$val))
  misVal <- 0
  DF$val[DF$val == minAlpha] <- misVal
  DF$color[DF$val == misVal] <- NA
  
  MSEtext <- grid::textGrob(MSENames[i], gp=grid::gpar(fontsize=MSEsize),
                            rot=90)

  P <- ggplot(DF) + 
    geom_tile(na.rm=TRUE, aes(x=x, y=y, fill=color, alpha=val)) +
    scale_fill_manual(values=values,na.value = 'white') +
    theme_classic() +
    guides(fill=FALSE, alpha=FALSE) +
    # geom_hline(yintercept = 1, linetype=2, color='darkgray') +
    # geom_vline(xintercept = 1, linetype=2, color='darkgray') +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    theme(plot.title = element_text(hjust = 0.5, size=TitleSize)) +
    geom_text(data=textDF, aes(x=x, y=y, label=label), size=TextSize) 
  if (mm == 1) {
    P <- P + labs(y=expression(F/F['MSY']), x='') +
      theme(plot.margin = unit(c(5.5,5.5,5.5,30), "pt")) +
      annotation_custom(MSEtext, xmin=MSEloc,xmax=MSEloc,ymin=1,ymax=1) +  
      coord_cartesian(clip = "off")
  } else {
    P <- P + labs(y='', x='') +
      theme(axis.title.y = element_blank(),
            axis.text.y=element_blank())
  }
  
  if (i != 3) {
    P <- P + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank())
  } else {
    P <- P + labs(x=expression(SB/SB['MSY']))
  }
 P + theme(axis.text = element_text(size=text.size),
           axis.title = element_text(size=title.size),
           axis.line = element_blank())
}

Plist <- list() 
count <- 0
for (i in 1:3) {
  MSE <- readRDS(file.path('MSEs', MSEs[i]))
  for (mm in 1:MSE@nMPs) {
    count <- count + 1
    P <- makeKplot(MSE, mm=mm, i=i, MSENames=MSENames,
                   MSEloc=-1,
                   TextSize=3,
                   TitleSize=10,
                   MSEsize=10,
                   text.size=8,
                   title.size=8,
                   minAlpha = 0.4)
    if (count <=9) P <- P + labs(title=MSE@MPs[mm])
    Plist[[count]] <- P
  }
}
Pout <- cowplot::plot_grid(plotlist=Plist, nrow=3, ncol=9,
                   rel_widths = c(1.6, rep(1,8)),
                   rel_heights = c(1,0.8,1))

ggsave('Figures/Figure_3.png', Pout, width=12, height=4, dpi=600)


# --- Trade-Off Plots ----
Plist <- list()

for (i in 1:3) {
  lab.size <- 3
  MSE <- readRDS(file.path('MSEs', MSEs[i]))
  
  dom <- Dom(MSE, PMlist=c('PM1', 'PM2', 'PM3', 'PM4', 'PM5'),
             Yrs=YrsList)
  mps <- dom$MPs
  failMPs <- MSE@MPs[which(PM4(MSE)@Mean < 0.5)]
  
  pnof <- PNOF(MSE, Yrs=-20)
  lty <- Yield(MSE, Yrs=-10)
  sty <- Yield(MSE, Yrs=10)
  pb0.5 <- P50(MSE, Yrs=-20)
  pb1 <- P100(MSE, Yrs=-20)
  
  Captions <- list(expression(Prob. ~F ~ '<'~F["MSY"]),
                   expression(Prob. ~SB ~ '>'~0.5*SB["MSY"]),
                   expression(Prob. ~SB ~ '>'~SB["MSY"]),
                   'Average Short-Term Yield',
                   'Average Long-Term Yield')
  
  
  
  DF <- data.frame(PNOF=pnof@Mean,
                   PB0.5=pb0.5@Mean,
                   PB1=pb1@Mean,
                   STY=sty@Mean,
                   LTY=lty@Mean,
                   label=pnof@MPs,
                   class=c(1,2,3,3,4,4,5,5,5))
  DF$class <- as.factor(DF$class)
  levels(DF$class) <- c("Catch Only",
                        'Catch & Depletion',
                        "Length Targeting",
                        "Index Targeting",
                        "Stock Assessment")
  
  PMlist <- list(pnof, pb0.5, pb1, sty, lty)

  labels <- list(paste0(letters[1:4],")"),
                 paste0(letters[5:8],")"),
                 paste0(letters[9:12],")"))
      
  DF$font <- 'bold'
  DF$font[DF$label %in% failMPs] <- 'italic'
  DF$font[DF$label %in% dom$DomMPs[,1]] <- 'plain'
  
  for (x in 1:4) {
    xdat <- colnames(DF)[x]
    pout <- ggplot(DF, aes_string(x=xdat, y='LTY', color='class',
                                  shape='class')) + 
      geom_point() + geom_text_repel(aes(label=label, fontface=font), size=lab.size, 
                                     show.legend = FALSE) +
      expand_limits(x=c(0,1), y=c(0,1)) +
      theme_bw() +
      labs(x='', y='') 
      # labs(x='', y='', tag=labels[[i]][x]) 
    if (x==2) {
      pout <- pout + geom_rect(aes(xmin=0, xmax=0.5, ymin=0, ymax=Inf), alpha=0.05, color=NA)
    }
    if (x ==1) {
      pout <- pout + labs(y=Captions[[5]]) +
        theme(axis.title.y=element_text(),
              axis.text.y=element_text())
    }
    if (i == 3) {
      pout <- pout + labs(x=Captions[[x]]) +
        theme(axis.title.x=element_text(),
              axis.text.x=element_text())
    }
    if (x !=1) pout <- pout + theme(axis.title.y = element_blank(),
                                    axis.text.y = element_blank())
    if (i !=3) pout <- pout + theme(axis.title.x = element_blank(),
                                    axis.text.x = element_blank())
    
    if (x== 4) {
      pout <- pout + 
        scale_x_continuous(expand = c(0, 0), limits = c(0,2)) +
        scale_y_continuous(expand = c(0, 0), limits = c(0,1.1))
    } else {
      pout <- pout + 
        scale_x_continuous(expand = c(0, 0), limits = c(0,1.1)) +
        scale_y_continuous(expand = c(0, 0), limits = c(0,1.1))
    }
    pout <- pout + labs(color='MP Class', shape='MP Class')
    legend <- get_legend(
      # create some space to the left of the legend
      pout + theme(legend.box.margin = margin(0, 0, 0, 12))
    )
    pout <- pout + theme(legend.position="none") 
    assign(paste0('p',x), pout)
  }
  
  pout2 <- cowplot::plot_grid(p1, p2, p3, p4, nrow=1, ncol=4,
                              rel_widths = c(1, rep(0.85,3)))
  
  Plist[[i]] <- gridExtra::grid.arrange(pout2,
                                        left=grid::textGrob(MSENames[i], 
                                                             rot = 90, 
                                                             vjust = 1)) 
}

Plots <- cowplot::plot_grid(plotlist = Plist, nrow=3, ncol=1, align='vh',
                   rel_heights = c(rep(0.85,2),1))

Pout <- cowplot::plot_grid(Plots, legend, ncol=2, rel_widths = c(1, 0.2))
                           

ggsave('Figures/Figure_4.png', Pout, width=12, height=7, dpi=600)

# ---- Projection Plots of Non-Dominated MPs ----

pnof <- PNOF(MSE, Yrs=-20)
lty <- Yield(MSE, Yrs=-20)
sty <- Yield(MSE, Yrs=10)
pb0.5 <- P50(MSE, Yrs=-20)
pb1 <- P100(MSE, Yrs=-20)

PM1 <- PNOF
PM2 <- Yield
PM3 <- Yield
PM4 <- P50
PM5 <- P100
YrsList <- list(PM1=-20, PM2=-20, PM3=10, PM4=-20, PM5=-20)

fignum <- 4
simsamps <- sample(1:MSE@nsim, 30)

for (i in 1:3) {
  fignum <- fignum +1
  MSE <- readRDS(file.path('MSEs', MSEs[i]))
  dom <- Dom(MSE, PMlist=c('PM1', 'PM2', 'PM3', 'PM4', 'PM5'),
             Yrs=YrsList)
  mps <- dom$MPs

  MSE2 <- Sub(MSE, mps)
  
  MSE2@MPs[which(PM4(MSE2)@Mean < 0.5)]
  
  MSE2 <- Sub(MSE2, MPs=which(PM4(MSE2)@Mean > 0.5))
  
  DFlist <- list()
  for (mm in 1:MSE2@nMPs) {
    b_bmsy <- MSE2@B_BMSY[,mm,] %>% as.vector()
    f_fmsy <- MSE2@F_FMSY[,mm,] %>% as.vector()
    relC <- MSE2@C[,mm,] / array(MSE2@OM$RefY, dim=dim(MSE2@C[,mm,]))
    rownames(relC) <- NULL
    relC <- relC %>% as.vector()
    DFtemp <- data.frame(B=b_bmsy, F=f_fmsy, C=relC, 
                         Year=rep(1:MSE2@proyears, each=MSE2@nsim),
                         Sim=rep(1:MSE2@nsim, MSE2@proyears))
    
    
    Ccolors <- DFtemp %>% group_by(Sim) %>% filter(Year==max(Year)) %>% select(Sim, C)
    Ccolors$CColor <- 'low'
    Ccolors$CColor[Ccolors$C>0.5 & Ccolors$C<1] <- 'med'
    Ccolors$CColor[Ccolors$C>1] <- 'high'
    Ccolors$C <- NULL
    
    Bcolors <- DFtemp %>% group_by(Sim) %>% filter(Year==max(Year)) %>% select(Sim, B)
    Bcolors$BColor <- 'low'
    Bcolors$BColor[Bcolors$B>0.5 & Bcolors$B<1] <- 'med'
    Bcolors$BColor[Bcolors$B>1] <- 'high'
    Bcolors$B <- NULL
    
    Fcolors <- DFtemp %>% group_by(Sim) %>% filter(Year==max(Year)) %>% select(Sim, F)
    Fcolors$FColor <- 'low'
    Fcolors$FColor[Fcolors$F>0.5 & Fcolors$F<1] <- 'med'
    Fcolors$FColor[Fcolors$F>1] <- 'high'
    Fcolors$F <- NULL
    DFtemp <- left_join(DFtemp, Ccolors, by='Sim')
    DFtemp <- left_join(DFtemp, Bcolors, by='Sim')
    DFtemp <- left_join(DFtemp, Fcolors, by='Sim')
    DFtemp$MP <- MSE2@MPs[mm]
    DFlist[[mm]] <- DFtemp
  }
  DF <- do.call('rbind', DFlist)
  
  
  values <- c("red", "orange", "green")
  alpha <- 0.4
  lsize <- 0.8
  lsize2 <- 0.5
  lsize3 <- 0.5
  textsize <- 8
  titlesize <- 8
  quants <- c(0.1, 0.9)
  # Catch Projection 
  QuantDF <- DF %>% group_by(Year, MP) %>% 
    summarise(med=median(C), low=quantile(C,quants[1]), high=quantile(C,quants[2]))
  QuantDF$MP <- factor(QuantDF$MP, levels=MSE2@MPs, ordered = TRUE)

  DF_sub <- DF %>% filter(Sim %in% simsamps)
  DF_sub$MP <- factor(DF_sub$MP, levels=MSE2@MPs, ordered = TRUE)
  
  
  DF_sub$CColor <- factor(DF_sub$CColor, levels=c('low', 'med', 'high'),
                          ordered = TRUE)
  DF_sub$BColor <- factor(DF_sub$BColor, levels=c('low', 'med', 'high'),
                          ordered = TRUE)
  DF_sub$FColor <- factor(DF_sub$FColor, levels=c('low', 'med', 'high'),
                          ordered = TRUE)
  
  Cplot <- ggplot(DF_sub, aes(x=Year, y=C)) + 
    facet_grid(.~MP) +
    geom_line(aes(group=as.factor(Sim), color=CColor), alpha=alpha, size=lsize2) +
    geom_line(data=QuantDF, aes(x=Year, y=med), size=lsize) +
    geom_line(data=QuantDF, aes(x=Year, y=low), linetype=2, size=lsize) +
    geom_line(data=QuantDF, aes(x=Year, y=high), linetype=2, size=lsize) +
    xlim(c(1,30)) + 
    coord_cartesian(ylim=c(0, 2)) +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x=element_blank(),
          axis.title =element_text(size=titlesize),
          axis.text =element_text(size=textsize)) +
    labs(y="Relative Yield") +
    scale_color_manual(values=values) +
    guides(color=FALSE) +
    geom_hline(yintercept = 1, linetype=2, color='darkgray', size=lsize3) +
    geom_hline(yintercept = 0.5, linetype=3, color='darkgray', size=lsize3)
  
  # Biomass Projection 
  QuantDF <- DF %>% group_by(Year, MP) %>% 
    summarise(med=median(B), low=quantile(B,quants[1]), high=quantile(B,quants[2]))
  QuantDF$MP <- factor(QuantDF$MP, levels=MSE2@MPs, ordered = TRUE)
  
  Bplot <- ggplot(DF_sub, aes(x=Year, y=B)) + 
    facet_grid(.~MP) +
    geom_line(aes(group=as.factor(Sim), color=BColor), alpha=alpha, size=lsize2) +
    geom_line(data=QuantDF, aes(x=Year, y=med), size=lsize) +
    geom_line(data=QuantDF, aes(x=Year, y=low), linetype=2, size=lsize) +
    geom_line(data=QuantDF, aes(x=Year, y=high), linetype=2, size=lsize) +
    xlim(c(1,30)) + 
    coord_cartesian(ylim=c(0, 2)) +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x=element_blank(),
          axis.title =element_text(size=titlesize),
          axis.text =element_text(size=textsize)) +
    labs(y=expression(SB/SB['MSY'])) +
    scale_color_manual(values=values) +
    guides(color=FALSE) +
    geom_hline(yintercept = 1, linetype=2, color='darkgray', size=lsize3) +
    geom_hline(yintercept = 0.5, linetype=3, color='darkgray', size=lsize3)
  
  # F mortality Projection
  QuantDF <- DF %>% group_by(Year, MP) %>% 
    summarise(med=median(F), low=quantile(F,quants[1]), high=quantile(F,quants[2]))
  QuantDF$MP <- factor(QuantDF$MP, levels=MSE2@MPs, ordered = TRUE)
  
  values <- rev(values)
  Fplot <- ggplot(DF_sub, aes(x=Year, y=F)) + 
    facet_grid(.~MP) +
    geom_line(aes(group=as.factor(Sim), color=FColor), alpha=alpha, size=lsize2) +
    geom_line(data=QuantDF, aes(x=Year, y=med), size=lsize) +
    geom_line(data=QuantDF, aes(x=Year, y=low), linetype=2, size=lsize) +
    geom_line(data=QuantDF, aes(x=Year, y=high), linetype=2, size=lsize) +
    xlim(c(1,30)) + 
    coord_cartesian(ylim=c(0, 2)) +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title =element_text(size=titlesize),
          axis.text =element_text(size=textsize)) +
    labs(y=expression(F/F['MSY'])) +
    scale_color_manual(values=values) +
    guides(color=FALSE)+ 
    geom_hline(yintercept = 1, linetype=2, color='darkgray', size=lsize3) +
    geom_hline(yintercept = 0.5, linetype=3, color='darkgray', size=lsize3)
  
  Pout <- cowplot::plot_grid(Cplot, Bplot, Fplot, nrow=3, ncol=1,
                             align='v',
                             rel_heights = c(0.9,0.8,1))
  
  width <- MSE2@nMPs
  height <- width/1.5
  ggsave(paste0('Figures/Figure_', fignum, '.png'), Pout,
         width=width, height=height)

}

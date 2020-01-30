


# ---- Functions ----
CombineOM <- function(OMlist, Name='CombinedOM', seed=1) {
  
  OMout <- new("OM")
  OMout@Name <- Name
  OMout@Agency <- OMlist[[1]]@Agency
  OMout@Region <- OMlist[[1]]@Region
  OMout@Sponsor<- OMlist[[1]]@Sponsor
  OMout@Latitude<- OMlist[[1]]@Latitude
  OMout@Longitude<- OMlist[[1]]@Longitude
  OMout@Common_Name <- OMlist[[1]]@Common_Name
  OMout@Source <- OMlist[[1]]@Source
  OMout@Species<- OMlist[[1]]@Species
  
  OMout@seed <- seed 
  
  nsim <- lapply(OMlist, slot, "nsim") %>% unlist()
  OMout@nsim <- sum(nsim)
  
  # Single value parameters
  slots <- c("proyears", 'interval', 'pstar', 'maxF', 
             'reps', 'maxage', 'a', 'b')
  OMout <- combine_singles(OMlist, slots, OMout)
  
  # Custom Parameters 
  cnms <- names(OMlist[[1]]@cpars)
  for (nm in cnms) {
    Clist <- lapply( lapply(OMlist, slot, "cpars") , '[[', nm)
    OMout@cpars[[nm]] <- abind::abind(Clist, along=1)
  }
  
  
  slots <- slotNames("Stock")
  slots <- slots[!grepl("Name", slots)]
  slots <- slots[!grepl("Source", slots)]
  slots <- slots[!grepl("grad", slots)]
  OMout <- updateOMslots(OMlist, slots, OMout)
  
  # Mortality
  if (!'M_ageArray' %in% names(OMout@cpars)) {
    stop("M_ageArray must be in cpars")
  } else {
    OMout@M <- c(0,0)
  }
  
  slots <- slotNames("Fleet")
  slots <- slots[!grepl("Name", slots)]
  OMout <- updateOMslots(OMlist, slots, OMout)
  
  # Find
  if (!'Find' %in% names(OMout@cpars)) {
    stop("Find must be in cpars")
  } else {
    OMout@EffUpper <- OMout@EffLower <- OMout@EffYears <- c(0,0)
  }
  
  # Obs
  slots <- slotNames("Obs")
  slots <- slots[!grepl("Name", slots)]
  OMout <- updateOMslots(OMlist, slots, OMout)
  
  # Imp
  slots <- slotNames("Imp")
  slots <- slots[!grepl("Name", slots)]
  OMout <- updateOMslots(OMlist, slots, OMout)
  
  OMout@nyears <- OMout@nyears[1]
  OMout@maxage <- OMout@maxage[1]
  OMout
}

updateOMslots <- function(OMlist, slots, OMout) {
  nsims <- lapply(OMlist, slot, "nsim") %>% unlist()
  for (sl in slots) {
    vals <- lapply(OMlist, slot, sl) %>% unlist()
    un.vals <- unique(vals)
    if (length(un.vals)==length(OMlist)) {
      # add to cpars
      OMout@cpars[[sl]] <- unlist(lapply(seq_along(OMlist), function (i) rep(un.vals[i], nsims[i])))
      slot(OMout, sl) <- c(0,0)
    } else if  (length(un.vals)==2) {
      # lower and upper bounds 
      slot(OMout, sl) <- un.vals
    }else if(length(un.vals)==1) {
      # check length 
      in.vals <- slot(OMlist[[1]], sl) 
      slot(OMout, sl) <- rep(un.vals, length(in.vals))
    } else {
      # do nothing
    }
  }
  OMout
}
combine_singles <- function(OMlist, slots, OMout) {
  for (sl in slots) {
    vals <- lapply(OMlist, slot, sl) %>% unlist()
    if (length(unique(vals))!= 1) stop(vals, " is not consistent between OMs")
    slot(OMout, sl) <- unique(vals)
  }
  OMout
}



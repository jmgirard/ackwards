ExtendedBassAckwards <- function (r,num.comp = 1,fm = "pca",rotate = "varimax",scores = "tenBerge")
  #r is a correlation matrix, fm is "pfa" or "minres" (EFA), rotations can be specified that are used in psych::pca or psych::fa
{
  comp.corr <- list() #initializing empty objects that are populated below
  cong <- list()
  comp.load <- list()
  pcas <- list()

  if (fm == "pca") {
    for (c in 1:num.comp) {
      #this is running PCAs at each level
      pcas[[c]] <-
        psych::pca(r, nfactors = c, rotate = rotate) #this is a list of all pcas
      comp <-
        psych::pca(r, nfactors = c, rotate = rotate) #pca loop for each level of hierarchy
      colnames(comp$loadings) <-
        paste0(letters[c], 1:ncol(comp$loadings)) #each level gets a different letter (Level 1: A1, Level 2: B1, B2, etc.), rather than all being F1, F2, etc." regardless of level
      comp.load[[c]] <-
        fa.sort(comp$loadings) #component loadings for each level, sorted by size
      comp.weights <- comp$weights
      colnames(comp.weights) <-
        paste0(letters[c], 1:ncol(comp.weights)) #naming component weights as above (A1, B1, B2, C1, C2, C3, etc.)
      unsort.loadings <- comp$loadings #loadings sorted by placement in input matrix
      colnames(unsort.loadings) <-
        paste0(letters[c], 1:ncol(unsort.loadings)) #naming loadings as above (A1, B1, B2, C1, C2, C3, etc.)
      if (c > 1) {
        for (i in 1:(c - 1)) {
          colnames(pcas[[i]]$weights) <-
            paste0(letters[i], 1:ncol(pcas[[i]]$weights)) #names weights matrices sequentially
          colnames(pcas[[i]]$loadings) <-
            paste0(letters[i], 1:ncol(pcas[[i]]$loadings)) #names loadings matrices sequentially
          comp.corr[[length(comp.corr) + 1]] <-
            t(pcas[[i]]$weights) %*% r %*%  comp.weights #calculates the correlations between levels
          cong[[length(cong) + 1]] <-
            psych::factor.congruence(pcas[[i]]$loadings, unsort.loadings) #calculates congruence coefficients between levels
        }
      }
    }
  }
  else if (fm == "minres") {
    for (c in 1:num.comp) {
      #this is running EFAs at each level
      pcas[[c]] <-
        psych::fa(r,  nfactors = c, fm = fm, rotate = rotate, scores = scores) #this is a list of all EFAs
      comp <-
        psych::fa(r, nfactors = c, fm = fm, rotate = rotate, scores = scores) #EFA loop for each level of hierarchy
      colnames(comp$loadings) <-
        paste0(letters[c], 1:ncol(comp$loadings)) #naming loadings as above (A1, B1, B2, C1, C2, C3, etc.)
      comp.load[[c]] <-
        psych::fa.sort(comp$loadings) #factor loadings for each level, sorted by size
      comp.weights <- comp$weights
      colnames(comp.weights) <-
        paste0(letters[c], 1:ncol(comp.weights)) #naming factor weights as above (A1, B1, B2, C1, C2, C3, etc.)
      unsort.loadings <- comp$loadings #loadings sorted by placement in input matrix
      colnames(unsort.loadings) <-
        paste0(letters[c], 1:ncol(unsort.loadings)) #naming loadings as above (A1, B1, B2, C1, C2, C3, etc.)
      if (c > 1) {
        for (i in 1:(c - 1)) {
          colnames(pcas[[i]]$weights) <-
            paste0(letters[i], 1:ncol(pcas[[i]]$weights)) #names weights matrices sequentially
          colnames(pcas[[i]]$loadings) <-
            paste0(letters[i], 1:ncol(pcas[[i]]$loadings)) #names loading matrices sequentially
          comp.corr[[length(comp.corr) + 1]] <-
            t(pcas[[i]]$weights) %*% r %*%  comp.weights #calculates the correlations between levels
          cong[[length(cong) + 1]] <-
            psych::factor.congruence(pcas[[i]]$loadings, unsort.loadings) #calculates the congruence coefficients between levels
        }
      }
    }
  }

  else {
    stop('fm must be pca or minres') #warning that factor method needs to be principal components (fm = "pca") or 'minres' factoring (fm = "minres")
  }

  result <- list(
    #these are the objects returned in the ExtendedBassAckwards result
    comp.corr = comp.corr,
    pcas = pcas,
    cong = cong,
    r = r,
    comp.load = comp.load
  )
  return(result)
}

ChaseCorrPaths <- function (comp.corr, component = "levelnum") #This function is called on below, not used alone. It calculates the paths of correlations >.9
{
  chased_levels <- vector() #initializing empty objects that are populated below
  chased_to_level <- vector()
  chased_to <- list()
  sub_revcomp.corr <- list()
  for (i in (length(comp.corr):1))
  {
    if (component %in% colnames(comp.corr[[i]]))
      #if the component we're interested is in the matrix
      chased_levels[[i]] <-
        (max(comp.corr[[i]][, component]) >= .9) #tell me if the maximum component correlation for the relevant column is >=.9
  }
  revcomp.corr <-
    rev(comp.corr) #reverse order of comp.corr to work from the bottom up
  component_level <-
    (length(chased_levels[!is.na(chased_levels)]) + 1) #level of current component (calculated as number of upward comparison matrices +1)
  chased_levels <-
    rev(chased_levels) #reverses order, so looking at lowest levels of hierarchy first

  if (any(chased_levels, na.rm = TRUE))
  {
    chased_levels <-
      (which.min(chased_levels)) - 1 #counts number of consecutive true values before first false
  }
  else {
    chased_levels <- 0
  }

  chased_to_level <- (component_level - chased_levels)

  if (chased_levels == 0)
  {
    chased_to <- "null"
  } #if no trues
  else {
    #isolate block of matrices in revcomp.corr relevant to component
    #end range
    end_comp.corr <- ((component_level * (component_level - 1) / 2))
    #start range
    start_comp.corr <- (end_comp.corr - (component_level - 2))

    sub_comp.corr <-
      comp.corr[start_comp.corr:end_comp.corr] #subset of matrices

    chased_to <-
      rownames(as.data.frame(which.max(sub_comp.corr[[chased_to_level]][, component]))) #name of component chased to

  }
  if ((component == "b1") &
      (max(comp.corr[[1]][, "b1"]) >= .9))
    #need to calculate the b->a level separately
  {
    chased_to <- "a1"
  }
  if ((component == "b2") & (max(comp.corr[[1]][, "b2"]) >= .9))
  {
    chased_to <- "a1"
  }
  result <- list(component, chased_to)
  result <- paste(result, sep = " ", collapse = "--")
  return(result)
}

ChaseCongPaths <- function (cong, component = "levelnum") #This function is called on below, not used alone. It calculates the paths of congruence coefficients >.95
{
  chased_levels <- vector() #initializing empty objects that are populated below
  chased_to_level <- vector()
  chased_to <- list()
  sub_revcong <- list()
  for (i in (length(cong):1))
  {
    if (component %in% colnames(cong[[i]]))
      #if the component we're interested is in the matrix
      chased_levels[[i]] <-
        (max(cong[[i]][, component]) > .95) #tell me if the maximum congruence coefficient for the relevant column is >.95
  }
  revcong <- rev(cong) #reverse order of cong to work from bottom up
  component_level <-
    (length(chased_levels[!is.na(chased_levels)]) + 1) #level of current component (calculated as number of upward comparison matrices +1)
  chased_levels <-
    rev(chased_levels) #reverses order, so looking at lowest levels of hierarchy first
  if (any(chased_levels, na.rm = TRUE))
  {
    chased_levels <-
      (which.min(chased_levels)) - 1 #counts number of consecutive true values before first false
  }
  else {
    chased_levels <- 0
  }

  chased_to_level <- (component_level - chased_levels)

  if (chased_levels == 0)
  {
    chased_to <- "null"
  } #if no trues
  else {
    #isolate block of matrices in revcomp.corr relevant to component
    #end range
    end_cong <- ((component_level * (component_level - 1) / 2))
    #start range
    start_cong <- (end_cong - (component_level - 2))

    sub_cong <- cong[start_cong:end_cong] #subset of matrices

    chased_to <-
      rownames(as.data.frame(which.max(sub_cong[[chased_to_level]][, component]))) #name of component chased to

  }
  if ((component == "b1") &
      (max(cong[[1]][, "b1"]) > .95))
    #need to calculate the b->a level separately
  {
    chased_to <- "a1"
  }
  if ((component == "b2") & (max(cong[[1]][, "b2"]) > .95))
  {
    chased_to <- "a1"
  }
  result <- list(component, chased_to)
  result <- paste(result, sep = " ", collapse = "--")
  return(result)
}


FindRedundantComp <- function (comp.corr, cong, last_component = "levelnum") #This function returns summary output that lists the paths of redundant variables based on correlations and congruence coefficients
{
  chase <- list() #initializing empty objects that are populated below
  congruence <- list()
  largest <- list()
  cross <- list()
  output <- list()
  mat <- matrix(data = 1:625,
                nrow = 25,
                ncol = 25)
  for (i in 1:25)
  {
    mat[, i] <- i
  }
  mat[upper.tri(mat)] <- NA
  matlett <- matrix(data = 1:625,
                    nrow = 25,
                    ncol = 25)
  for (i in 1:25)
  {
    matlett[i, ] <- i
  }
  matlett[] <- letters[as.matrix(matlett)]
  matlett[upper.tri(matlett)] <- NA
  comp_list <-
    paste0(na.omit(as.vector(t(matlett))), na.omit(as.vector(t(mat)))) #list of all possible component names up to 25 hierarchical levels

  for (i in 2:(which(comp_list == last_component)[[1]]))
  {
    chase[i] <- ChaseCorrPaths(comp.corr, comp_list[i])
    congruence[i] <- ChaseCongPaths(cong, comp_list[i])

  }

  output <- list(corr.chase = chase, cong.chase = congruence)
  return(output)
}



HierarchicalCorrs <- function (comp.corr, component = "levelnum") #The output includes censored versions of the correlation matrices that isolate the relevant associations for each component or factor
{
  comp.corrcross <- list()
  comp.corrcurr <- data.frame()

  for (i in 1:length(comp.corr))
    #running forward through the list
  {
    comp.corrcurr <- as.data.frame(comp.corr[[i]])
    if (component %in% colnames(comp.corrcurr))
      #if the component we're interested is in the matrix
    {
      comp.corrcurr[abs(comp.corrcurr) < .3] <-
        NA #remove small loads #can update this to a smaller value to reveal more detailed relationships (e.g., change < .3 to < .1)
      comp.corrcurr[abs(comp.corrcurr) >= .9] <-
        NA #remove chase loads
      comp.corrcross[[i]] <- comp.corrcurr
    }
    else
    {
      comp.corrcross[[i]] <- NULL
    }
  }

  comp.corrcross <-
    comp.corrcross[-which(sapply(comp.corrcross, is.null))]

  if (component == "b1" || component == "b2")
  {
    comp.corrcross <- comp.corr[[1]]
  }


  return(comp.corrcross)
}


NegativeCorrs <- function (comp.corr, component = "levelnum") #The output includes censored versions of the correlation matrices that isolate the relevant associations for each component or factor
{
  comp.corrcross <- list()
  comp.corrcurr <- data.frame()

  for (i in 1:length(comp.corr))
    #running forward through the list
  {
    comp.corrcurr <- as.data.frame(comp.corr[[i]])
    if (component %in% colnames(comp.corrcurr))
      #if the component we're interested is in the matrix
    {
      comp.corrcurr[comp.corrcurr > -.3] <-
        NA #remove small negative and all positive loads #can update this to reveal more detailed relationships (e.g., change > -.3 to > -.1)
      comp.corrcross[[i]] <- comp.corrcurr
    }
    else
    {
      comp.corrcross[[i]] <- NULL
    }
  }

  comp.corrcross <-
    comp.corrcross[-which(sapply(comp.corrcross, is.null))]

  if (component == "b1" || component == "b2")
  {
    comp.corrcross <- comp.corr[[1]]
  }


  return(comp.corrcross)
}

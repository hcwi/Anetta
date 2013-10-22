
# Run processing: find, read, model, save
run <- function() {
  
  print("[debug] run")
  
  d0 <- read.table("data.txt", header=T)
  d0.names <- grep(names(d0), pattern="(Characteristics)|(Trait.Value.)", val=T)
  d <- d0[d0.names]
  t.names <- grep(names(d), pattern="Trait.Value.", val=T)
  for (t in t.names) {
    print(t)
    f <- paste(t, sep="~", "Characteristics.Infra.specific.name.*Characteristics.Rwness.*Characteristics.Shape.")
    by <- list("Characteristics.Infra.specific.name.","Characteristics.Rwness.","Characteristics.Shape.")
    means <- aggregate(f, d, FUN="mean")
    sd <- aggregate(form, d, FUN="sd")
    res <- merge(means, sd, by=by)
  }

} 


# Calculate models and final statistics
get.models <- function(sad, sad.names) {
  
  print("[debug] get.models")
  prepare.libs()  
  
  traits <<- get.traits(sad)
  results <<- prepare.results(sad)
  
  success <- length(traits);
  models <- list()
  #TODO uncomment
  for (i in 2:2) {
    
    trait = traits[i] 
    
    tmp <- tryCatch(
      get.model.for.trait(trait, sad, results),    
        error = function(e){
          write(paste("ERROR: Evaluating of model for trait '", sad.names[trait], "' failed. The following error occured: ", e, sep=""), stderr())
          return(NA)
        },                 
        warning = function(w) {
          write(paste("WARNING: Evaluating of model for trait '", sad.names[trait], "' produced the following warning: ", w, sep=""), stderr())
          return(tmp)
        },
        finally = function() {}
    )
    
    if (!is.na(tmp[1])) {
      results <- tmp$results
      info <- tmp$info    
      models <- c(models,list(info))
      success <- success-1;
    }
    
  }
  
  list(means=results, models=models, success=success)
}


# Calculate model and statistics for a trait 
get.model.for.trait <- function(trait, sad, results) {
 
  print(paste("Models for a new trait:", trait))
  
  sadt <- sad[!is.na(sad[,trait]),]
  #sadt <- subset(sad, Characteristics.Infra.specific.name. %in% c("MARESI", "POMO"))
  sadt <<- sad[!sad[,"Characteristics.Infra.specific.name."] %in% c("MARESI", "POMO"),]
  
  fixed <- get.fixed(sadt)
  fixef <- paste(fixed, collapse="+")   
  random <- get.random(sadt)
  ranef <- paste("+(1|", random, ")", sep="", collapse="")
  
  factorizenames <- function(x) {
    if (is.numeric(sadt[,x]) && 
          length(grep("Trait[.]Value", names(sadt)[x])) > 0) 
      sadt[,x] 
    else factor(sadt[,x])
  }
  sadtf <<- lapply(seq_along(sadt), FUN=factorizenames)
  names(sadtf) <<- names(sadt)
  
  form <- paste(trait,"~",fixef, ranef, sep="")
  print(paste("Formula:", form))
  
  model <<- lmer(form, sadtf)
  
  # Set variances of random effects
  {
    for (j in 1:length(random)) {
      r <- random[j]
      v <- VarCorr(model)[[r]][1]
      
      is.variable <- !is.na(results[,r]=="*")
      is.variance <- results[,"Parameter"] == "Variance"
      results[is.variable & is.variance, trait] <- v   
    }
  }
  
  # Set sigma
  {
    errvar <- attr(VarCorr(model), "sc")^2
    is.errvar <- results[,"Parameter"] == "Error variance"
    results[is.errvar, trait] <- errvar
  }
  
  # Set means for fixed effects
  {
    results <- fill.means.for.fixed(sad, results, model, fixed, trait)
  }
  
  info <- list(trait=trait, fixed=fixed, random=random, model=model)
  
  list(results=results, info=info)
}



# Prepare table for results (combinations of effects)
prepare.results <- function(sad) {
  
  print("[debug] prepare.results")
  
  type <- "Parameter"
  form <- "Formula"
  traits <- get.traits(sad)
  traits.se <- paste("S.e.", traits, sep="")
  fixed <- get.fixed(sad)
  random <- get.random(sad)

  cols <- c(type, form, fixed, random, c(rbind(traits, traits.se)))
  ncols <- length(cols)
  
  results <- matrix(nrow=1, ncol=ncols)
  colnames(results) <- cols
  
  # General mean
  {
    srow <- 1
    results[srow, form] <- ""
    results[srow, type] <- "Mean"
  }
  
  # Means for fixed effects
  {
    srow <- 2
    for (i in 1:length(fixed)) {
      com <- combn(fixed, i)
      ncom <- dim(com)[2]
      for (j in 1:ncom) {
        names <- com[,j]
        dat <- unique(sad[names])
        dat <- as.matrix(dat[order(dat[1]),])
        
        nrows <- dim(dat)[1]
        if (is.null(nrows)) nrows <- length(dat)
        results <- rbind(results, matrix(nrow=nrows, ncol=ncols))
        
        formula <- paste(names, collapse="*")
        
        to <- srow + nrows - 1
        results[srow:to, names] <- dat
        results[srow:to, form] <- formula
        results[srow:to, type] <- "Mean"
        results
        srow <- to + 1
      }
    }
  }
  
  # Variances of random effects
  {
    nrows <- length(random)
    results <- rbind(results, matrix(nrow=nrows, ncol=ncols))
    for (i in 1:nrows) {
      results[srow, random[i]] <- "*"
      results[srow, type] <- "Variance"
      srow <- srow + 1
    }
  }
  
  # Error variance
  {
    results <- rbind(results, matrix(nrow=1, ncol=ncols))
    results[srow, type] <- "Error variance"
    srow <- srow + 1
  }

  results
}

# Calculate and fill estimated means for fixed effects
fill.means.for.fixed <- function(sad, results, model, fixed, trait) {
  
  print("[debug] fill.means.for.fixed")
 
  tmp <- prepare.matrices(sad, fixed)
  fix <- tmp$fix
  xu <- tmp$xu
  
  x <- unique(getME(model, "X"))
  est <- x %*% fixef(model)
  est.cov <- x %*% vcov(model) %*% t(x)
  
  for (i in 1:dim(fix)[1]) {
    
    factor <- fix[i,]$mform
    from <- fix[i,]$from
    to <- fix[i,]$to
    
    xf <- xu[,from:to]
    m <- solve(t(xf) %*% xf) %*% t(xf)
    means <- m %*% est
    rownames(means) <- colnames(xu)[from:to]
    
    means.var <- diag(m %*% est.cov %*% t(m))
    
    a <- results[,"Formula"] == factor
    b <- results[,"Parameter"] == "Mean"
    results[a&b, trait] <- means
    results[a&b, paste("S.e.",trait, sep="")] <- sqrt(means.var)
  }  
  
  results
}

# Prepare full model matrix and indices for it
prepare.matrices <- function(sad, fixed) {
  
  print("[debug] prepare.matrices")
  
  fix <- data.frame(mform="", pform="", from=1, to=1)
  
  x <- matrix(1, nrow=dim(sad)[1])
  colnames(x) <- "all"
  
  f=1
  for (i in 1:length(fixed)) {
    com <- combn(fixed, i)
    ncom <- dim(com)[2]
    for (j in 1:ncom) {
      names <- com[,j]
      mform <- paste(names, collapse="*") 
      pform <- paste(names, collapse=":") 
      
      formula <- paste("Sample.Name", sep="~", mform)
      print(paste("[debug]           formu: ", formula, sep=""))
      
      xtmp <- cast(sad, formula, length)
      xtmp <- xtmp[-1]
      x <- cbind(x, xtmp)
      
      from <- fix[f, "to"] + 1
      to <- from + dim(xtmp)[2] - 1
      f <- f + 1
      fix[f,] <- list(mform, pform, from, to)
    }
  }
  
  xu <- as.matrix(unique(x))  
  
  list(fix=fix, xu=xu)
}

# Install and load missing libraries
prepare.libs <- function() {
  
  if(!require("lme4")) {
    install.packages("lme4", repos='http://cran.us.r-project.org')
    library(lme4)
  }
  #if(!require("reshape")) {
  #  install.packages("reshape", repos='http://cran.us.r-project.org')
  #  library(reshape)
  #}
}

# Get traits
get.traits <- function(sad) {
  
  print("[debug] get.traits")
  
  are.traits <- grep("Trait[.]Value", names(sad), value=T)
  
#   warning("Removing traits with no variation")
#   have.var <- function(x) length(unique(x))>1
#   are.var <- sapply(sad[are.traits], have.var)
#   are.traits <- are.traits[are.var]
  
  are.traits
}

# Get fixed effects
get.fixed <- function(sad) {
  
  print("[debug] get.fixed")
  
  are.levels <- grep("((Characteristics)|(Factor))", names(sad), value=T)
  #warning("Filtering factors to exclude *id* names -- only for Keygene data. Remove for other analyses!")
  are.levels <- grep("[Ii]d", are.levels, value=T, invert=T)
  are.var <- sapply(sad[are.levels], FUN = function(x) length(unique(x))>1 )
  are.fixed  <- grep("(Block)|(Field)|(Rank)|(Plot)|(Replic)|(Column)|(Row)|(Rand)", are.levels[are.var], value=T, invert=T)
  
  are.fixed
}

#TODO rep in random only when no other random effects


# Get random effects
get.random <- function(sad) {
  
  print("[debug] get.random")
  
  are.levels <- grep("((Characteristics)|(Factor))", names(sad), value=T)
  #warning("Filtering factors to exclude *id* names -- only for Keygene data. Remove for other analyses!")
  are.levels <- grep("[Ii]d", are.levels, value=T, invert=T)
  are.var <- sapply(sad[are.levels], FUN = function(x) length(unique(x))>1 )
  are.random <- grep("(Block)|(Field)|(Rank)|(Plot)|(Replic)|(Column)|(Row)|(Rand)", are.levels[are.var], value=T)
  
  are.random
}


# Things to do before running in Java

# remove global variables (<<-)
# uncomment:
args <- commandArgs(TRUE)
if (length(args) > 0) {
  setwd(args[1])
}

options(stringsAsFactors=FALSE)
run()

#setwd("C:/Users/hnk/Desktop/AnettaISATABtoAnalyse")

#setwd("C:/Users/hcwi/Desktop/phenotypingTXT")
#setwd("C:/Users/hcwi/Desktop/phen-stats/isatab")
#setwd("C:/Users/hcwi/Desktop/phen-stats/isatab_missing")
#setwd("C:/Users/hcwi/Desktop/phen/src/test/resources/DataWUR")
#setwd("C:/Users/hcwi/Desktop/phen/src/test/resources/Phenotyping2")
#setwd("C:/Users/hcwi/Desktop/phen/src/test/resources/IPGPASData")
# XXX setwd("C:/Users/hcwi/Desktop/phen/src/test/resources/DataIPK") XXX too big file, very time consuming calculations (30min) with no result
#setwd("C:/Users/hcwi/Desktop/phen/src/test/resources/DataIPK2")
# XXx setwd("C:/Users/hcwi/Desktop/phen/src/test/resources/DataINRA") XXX mistake in ISA-TAB structure, Sample column in both s/a with different values
#setwd("C:/Users/hcwi/Desktop/phen/src/test/resources/DataINRA2")
#setwd("C:/Users/hcwi/Desktop/phen/src/test/resources/Keygene2")

# calculate means for all obs~factor combinations
#   what <- names(barley[sapply(barley, is.numeric)])
#   byWhat <- names(barley[sapply(barley, is.factor)])
#   mChar <- aggregate(barley[,what[j]]~barley[,byWhat[i]], FUN=mean)
# 

#P-VALS
#if(!require("languageR")) {install.packages("languageR", repos='http://cran.us.r-project.org')}
#library(languageR)
#m1.p <- pvals.fnc(m1)

#NULL MODEL COMPARISION
#m1.null <- lmer(t.length~1+(1|f.block), barley)
#anova(m1,m1.null)
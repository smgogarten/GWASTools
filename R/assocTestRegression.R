assocTestRegression <- function(genoData,
                                outcome,
                                model.type,
                                covar.list = NULL,
                                ivar.list = NULL,
                                gene.action.list = NULL,
                                dosage = FALSE,
                                scan.chromosome.filter = NULL,
                                scan.exclude = NULL,
                                CI = 0.95,
                                robust = FALSE,
                                LRtest = TRUE,
                                chromosome.set = NULL,
                                block.set = NULL,
                                block.size = 5000,
                                verbose = TRUE,
                                outfile = NULL){

##############################################################  
# Some Functions for Later Use
##############################################################
# A function to Transform Genotypes based on gene action
TransformGenotype <- function(index, geno){
  switch(index,
         additive = {},
         dominant = {geno[geno==1]=2},
         recessive = {geno[geno==1]=0},
         dominance = {geno[geno==2]=1-frq
                      geno[geno==0]=frq
                      geno[geno==1]=0},
         )
  return(geno)
}

###############################################################
# A function to Run the Regression and collect pertinent Output
RunRegression <- function(){
  ##### logistic regression #####
  if(model.type[md]=="logistic"){
    # check for monomorphic group
    if(mono | mono.cc0 | mono.cc1){
      mod <- TRUE
    # if not, try the regression
    }else{
      mod <- tryCatch(glm(model.formula, data = mdat, family=binomial()), warning=function(w) TRUE, error=function(e) TRUE);
    }

    # if monomorphic group or warning or error
    if(is.logical(mod)){
      # determine error type
      if(mono){
        err.type <- 2
      }else if(mono.cc0){
        err.type <- 0
      }else if(mono.cc1){
        err.type <- 1
      }else{
        err.type <- 9
      }
      # results entries
      if(liv[md] == 0 & !LRtest){
        tmp <- c(err.type,rep(NA,7));
      }else if(liv[md] == 0 & LRtest){
        tmp <- c(err.type,rep(NA,9));
      }else if(liv[md] > 0 & !LRtest){
        tmp <- c(err.type,rep(NA,7+7*liv[md]+2));
      }else{
        tmp <- c(err.type,rep(NA,7+9*liv[md]+4));
      }
      
    # if any coef are NA
    }else if(!all(!is.na(coef(mod)))){
      if(liv[md] == 0 & !LRtest){
        tmp <- rep(NA,8);
      }else if(liv[md] == 0 & LRtest){
        tmp <- rep(NA,10);
      }else if(liv[md] > 0 & !LRtest){
        tmp <- rep(NA,8+7*liv[md]+2);
      }else{
        tmp <- rep(NA,8+9*liv[md]+4);
      }
      
    # if model runs successfully 
    }else{
      if(robust==FALSE){
        # model based variance
        Vhat <- vcov(mod);
      }else{
        # sandwich variance
        Vhat <- vcovHC(mod, type="HC0");
      }
      
      # coefficient
      Est.G <- coef(mod)["genotype"];
      # SE
      SE.G <- sqrt(Vhat["genotype","genotype"]);
      # odds ratio
      OR <- exp(Est.G);
      # confidence interval
      LL <- exp(Est.G+qnorm((1-CI)*0.5)*SE.G);
      UL <- exp(Est.G+qnorm(1-((1-CI)*0.5))*SE.G);
      # Wald test
      W <- (Est.G/SE.G)^2
      pval <- pchisq(W, df=1, lower.tail=FALSE);
      # collect results
      tmp <- c(NA,Est.G,SE.G,OR,LL,UL,W,pval);
      # LR test
      if(LRtest & liv[md]==0){
        tmp <- append(tmp, unlist(lrtest(mod, "genotype"))[c(8,10)])
      }
      
      # if interaction terms
      if(liv[md]>0){
        # for each interaction term
        for(i in 1:liv[md]){
          # index for geno interaction variables
          GxE.idx <- grep(ivar.terms[i], names(coef(mod)));
          # coefficient
          Est.GxE <- coef(mod)[GxE.idx];
          # SE
          SE.GxE <- sqrt(Vhat[GxE.idx,GxE.idx]);
          # odds ratios
          OR.GxE <- exp(Est.GxE);
          # confidence intervals
          LL.GxE <- exp(Est.GxE+qnorm((1-CI)*0.5)*SE.GxE);
          UL.GxE <- exp(Est.GxE+qnorm(1-((1-CI)*0.5))*SE.GxE);
          # Wald test
          W.GxE <- (Est.GxE/SE.GxE)^2;
          pval.GxE <- pchisq(W.GxE, df=1, lower.tail=FALSE)
          # collect results
          tmp <- append(tmp,c(Est.GxE,SE.GxE,OR.GxE,LL.GxE,UL.GxE,W.GxE,pval.GxE));
          # LR test
          if(LRtest){
            tmp <- append(tmp, unlist(lrtest(mod, ivar.terms[i]))[c(8,10)])
          }
        }
        
        # Joint test
        # index for all genotype terms
        Gj.idx <- grep("genotype", names(coef(mod)));
        coef.j <- coef(mod)[Gj.idx]
        # Wald test
        W.Gj <- as.numeric(t(coef.j) %*% solve(Vhat[Gj.idx,Gj.idx]) %*% coef.j);
        pval.Gj <- pchisq(W.Gj, df=length(Gj.idx), lower.tail=FALSE)
        # collect results
        tmp <- append(tmp,c(W.Gj,pval.Gj))
        # LR test
        if(LRtest){
          tmp <- append(tmp, unlist(lrtest(mod, c("genotype",ivar.terms)))[c(8,10)])
        }
      }
    }
    
  ###### linear regression #####
  }else if(model.type[md]=="linear"){
    # check for monomorphic
    if(mono){
      mod <- TRUE
    # if not, try the regression
    }else{
      mod <- tryCatch(lm(model.formula, data = mdat), warning=function(w) TRUE, error=function(e)TRUE);
    }

    # if warning or error
    if(is.logical(mod)){
      # determine error type
      if(mono){
        err.type <- 2
      }else{
        err.type <- 9
      }
      # results entries
      if(liv[md] == 0 & !LRtest){
        tmp <- c(err.type,rep(NA,6));
      }else if(liv[md] == 0 & LRtest){
        tmp <- c(err.type,rep(NA,8));
      }else if(liv[md] > 0 & !LRtest){
        tmp <- c(err.type,rep(NA,6+6*liv[md]+2));
      }else{
        tmp <- c(err.type,rep(NA,6+8*liv[md]+4));
      }
      
    # if any coef are NA
    }else if(!all(!is.na(coef(mod)))){
      if(liv[md] == 0 & !LRtest){
        tmp <- rep(NA,7);
      }else if(liv[md] == 0 & LRtest){
        tmp <- rep(NA,9);
      }else if(liv[md] > 0 & !LRtest){
        tmp <- rep(NA,7+6*liv[md]+2);
      }else{
        tmp <- rep(NA,7+8*liv[md]+4);
      }

    # if model runs successfully 
    }else{
      if(robust==FALSE){
        # model based variance
        Vhat <- vcov(mod);
      }else{
        # sandwich variance
        Vhat <- vcovHC(mod, type="HC0");
      }
      
      # coefficient
      Est.G <- coef(mod)["genotype"];
      # SE
      SE.G <- sqrt(Vhat["genotype","genotype"]);
      # confidence interval
      LL <- Est.G+qnorm((1-CI)*0.5)*SE.G;
      UL <- Est.G+qnorm(1-((1-CI)*0.5))*SE.G;
      # Wald test
      W <- (Est.G/SE.G)^2
      pval <- pchisq(W, df=1, lower.tail=FALSE);
      # collect results
      tmp <- c(NA,Est.G,SE.G,LL,UL,W,pval);
      # LR test
      if(LRtest & liv[md]==0){
        tmp <- append(tmp, unlist(lrtest(mod, "genotype"))[c(8,10)])
      }
      
      # if interaction terms
      if(liv[md]>0){
        # for each interaction term
        for(i in 1:liv[md]){
          # index for geno interaction variables
          GxE.idx <- grep(ivar.terms[i], names(coef(mod)));
          # coefficient
          Est.GxE <- coef(mod)[GxE.idx];
          # SE
          SE.GxE <- sqrt(Vhat[GxE.idx,GxE.idx]);
          # confidence intervals
          LL.GxE <- Est.GxE+qnorm((1-CI)*0.5)*SE.GxE;
          UL.GxE <- Est.GxE+qnorm(1-((1-CI)*0.5))*SE.GxE;
          # Wald test
          W.GxE <- (Est.GxE/SE.GxE)^2;
          pval.GxE <- pchisq(W.GxE, df=1, lower.tail=FALSE)
          # collect results
          tmp <- append(tmp,c(Est.GxE,SE.GxE,LL.GxE,UL.GxE,W.GxE,pval.GxE));
          # LR test
          if(LRtest){
            tmp <- append(tmp, unlist(lrtest(mod, ivar.terms[i]))[c(8,10)])
          }
        }
        
        # Joint test
        # index for all genotype terms
        Gj.idx <- grep("genotype", names(coef(mod)));
        coef.j <- coef(mod)[Gj.idx]
        # Wald test
        W.Gj <- as.numeric(t(coef.j) %*% solve(Vhat[Gj.idx,Gj.idx]) %*% coef.j);
        pval.Gj <- pchisq(W.Gj, df=length(Gj.idx), lower.tail=FALSE)
        # collect results
        tmp <- append(tmp,c(W.Gj,pval.Gj))
        # LR test
        if(LRtest){
          tmp <- append(tmp, unlist(lrtest(mod, c("genotype",ivar.terms)))[c(8,10)])
        }
      }
    }

  ##### poisson regression #####
  }else if(model.type[md]=="poisson"){
    # check for monomorphic
    if(mono){
      mod <- TRUE
    # if not, try the regression
    }else{
      mod <- tryCatch(glm(model.formula, data = mdat, family=poisson()),warning=function(w) TRUE, error=function(e)TRUE);
    }

    # if warning or error
    if(is.logical(mod)){
      # determine error type
      if(mono){
        err.type <- 2
      }else{
        err.type <- 9
      }
      # results entries
      if(liv[md] == 0 & !LRtest){
        tmp <- c(err.type,rep(NA,7));
      }else if(liv[md] == 0 & LRtest){
        tmp <- c(err.type,rep(NA,9));
      }else if(liv[md] > 0 & !LRtest){
        tmp <- c(err.type,rep(NA,7+7*liv[md]+2));
      }else{
        tmp <- c(err.type,rep(NA,7+9*liv[md]+4));
      }
      
    # if any coef are NA
    }else if(!all(!is.na(coef(mod)))){
      if(liv[md] == 0 & !LRtest){
        tmp <- rep(NA,8);
      }else if(liv[md] == 0 & LRtest){
        tmp <- rep(NA,10);
      }else if(liv[md] > 0 & !LRtest){
        tmp <- rep(NA,8+7*liv[md]+2);
      }else{
        tmp <- rep(NA,8+9*liv[md]+4);
      }
      
    # if model runs successfully 
    }else{
      if(robust==FALSE){
        # model based variance
        Vhat <- vcov(mod);
      }else{
        # sandwich variance
        Vhat <- vcovHC(mod, type="HC0");
      }
      
      # coefficient
      Est.G <- coef(mod)["genotype"];
      # SE
      SE.G <- sqrt(Vhat["genotype","genotype"]);
      # relative risk
      RR <- exp(Est.G);
      # confidence interval
      LL <- exp(Est.G+qnorm((1-CI)*0.5)*SE.G);
      UL <- exp(Est.G+qnorm(1-((1-CI)*0.5))*SE.G);
      # Wald test
      W <- (Est.G/SE.G)^2
      pval <- pchisq(W, df=1, lower.tail=FALSE);
      # collect results
      tmp <- c(NA,Est.G,SE.G,RR,LL,UL,W,pval);
      # LR test
      if(LRtest & liv[md]==0){
        tmp <- append(tmp, unlist(lrtest(mod, "genotype"))[c(8,10)])
      }
      
      # if interaction terms
      if(liv[md]>0){
        # for each interaction term
        for(i in 1:liv[md]){
          # index for geno interaction variables
          GxE.idx <- grep(ivar.terms[i], names(coef(mod)));
          # coefficient
          Est.GxE <- coef(mod)[GxE.idx];
          # SE
          SE.GxE <- sqrt(Vhat[GxE.idx,GxE.idx]);
          # relative risk
          RR.GxE <- exp(Est.GxE);
          # confidence intervals
          LL.GxE <- exp(Est.GxE+qnorm((1-CI)*0.5)*SE.GxE);
          UL.GxE <- exp(Est.GxE+qnorm(1-((1-CI)*0.5))*SE.GxE);
          # Wald test
          W.GxE <- (Est.GxE/SE.GxE)^2;
          pval.GxE <- pchisq(W.GxE, df=1, lower.tail=FALSE)
          # collect results
          tmp <- append(tmp,c(Est.GxE,SE.GxE,RR.GxE,LL.GxE,UL.GxE,W.GxE,pval.GxE));
          # LR test
          if(LRtest){
            tmp <- append(tmp, unlist(lrtest(mod, ivar.terms[i]))[c(8,10)])
          }
        }
        
        # Joint test
        # index for all genotype terms
        Gj.idx <- grep("genotype", names(coef(mod)));
        coef.j <- coef(mod)[Gj.idx]
        # Wald test
        W.Gj <- as.numeric(t(coef.j) %*% solve(Vhat[Gj.idx,Gj.idx]) %*% coef.j);
        pval.Gj <- pchisq(W.Gj, df=length(Gj.idx), lower.tail=FALSE)
        # collect results
        tmp <- append(tmp,c(W.Gj,pval.Gj))
        # LR test
        if(LRtest){
          tmp <- append(tmp, unlist(lrtest(mod, c("genotype",ivar.terms)))[c(8,10)])
        }
      }
    }
    
  }
  ### return results ###
  return(tmp)
}

###############################################################
    
  # get SNP and chromosome information
  if(verbose) message("Determining chromosome and SNP information...")
  snpID <- getSnpID(genoData)
  scanID <- getScanID(genoData)
  chrom <- getChromosome(genoData)
  unique_chrom <- unique(chrom);
  nChromosomes <- max(chrom);
  rle_chrom <- rle(chrom);
  rle_chrom2 <- rep(0,nChromosomes);
  rle_chrom2[unique_chrom] <- rle_chrom$lengths;

  # create vector of first SNP number for each chromosome
  gbreaks <- vector(length=nChromosomes+1);
  gbreaks[1] <- 1;
  for(a in 1:nChromosomes)
    gbreaks[a+1] <- sum(rle_chrom2[1:a])+1;

  # chromosome.set check
  if(is.null(chromosome.set)){
    chromosome.set <- unique_chrom;
    allchrom <- TRUE
  } else {
    stopifnot(all(chromosome.set %in% unique_chrom))
    allchrom <- FALSE
  }

  # X chromosome check for sex variable
  if(XchromCode(genoData) %in% chromosome.set & !hasSex(genoData)){
    stop("Sex values for the samples are required to compute MAF for chromosome X SNPs");
  }

  # number of models
  nModels = length(model.type);
  # gene.action check
  if(verbose) message("Gene action check...")
  if(is.null(gene.action.list)){
    gene.action.list <- as.list(rep("additive",nModels))
  }
  if(dosage){
    if(!all(is.element(unlist(gene.action.list),"additive"))){
      stop("When using imputed genotype dosages, the gene.action must be additive.")
    }
  }
  if(!all(is.element(unlist(gene.action.list),c("additive","dominant","recessive","dominance")))){
    stop("At least one of the choices for gene.action is invalid. Valid options are additive, dominant, recessive, and dominance.")
  }

  # covar.list check
  if(verbose) message("Covariate check...")
  if(is.null(covar.list)){
    covar.list <- as.list(rep("",nModels))
  }

  # ivar.list check
  if(verbose) message("Interaction Variable check...")
  if(is.null(ivar.list)){
    ivar.list <- as.list(rep("",nModels))
  }
        
  ######  model checks and set names for results matrices ######
  if(verbose) message("Models checks...")
  names(covar.list) = paste("model",1:nModels,sep=".");

  # vector of variable names for results matrices
  nv <- "snpID";
  
  # gene.act.num counts the number of gene action types for each model in the analysis, it's used for indexing
  gene.act.num <- rep(NA,nModels);

  # liv counts the number of interaction variables for each model in the analysis
  liv <- rep(NA,nModels);


  for(j in 1:nModels){
    covar.names <- unique(unlist(strsplit(covar.list[[j]],"[*:]")));
    ivar.names <- unique(unlist(strsplit(ivar.list[[j]],"[*:]")));

    # covariate check
    if(!all(covar.names %in% getScanVariableNames(genoData))){
      stop(paste("Some of the variables in the covar.list for model",j," are not provided in the sample data."))
    }

    # interaction variable check
    if(!all(ivar.names %in% covar.names)){
      stop(paste("Some of the interaction variables in the ivar.list for model",j," are not in covar.list for that model."))
    }

    # determine the number of interaction variables
    liv[j] <- length(ivar.names)
             
    # Y chromosome check for sex variable
    if(YchromCode(genoData) %in% chromosome.set & "sex" %in% covar.names){
      warning(paste("Model",j,", Y chromosome tests are confounded with sex and should be run separately without sex in the model"));
    }

    # check for valid confidence level
    if(!(CI > 0 & CI < 1)){
      stop("Confidence Level must be between 0 and 1")
    }
    CIstr <- unlist(strsplit(as.character(CI), ".", fixed=TRUE))[2]

    # determine the number of gene action models
    gene.act.num[j] <- length(gene.action.list[[j]]);

    # create complete vector of variable names for results matrices
    if(model.type[j]=="logistic"){
      outcomeVar <- getScanVariable(genoData, outcome[j])
      if(!all(outcomeVar == 0 | outcomeVar == 1, na.rm=TRUE)){
        stop(paste("Model number",j,"is logistic, the corresponding outcome variable must be coded as 0/1"))
      }

      nv <- append(nv, paste(names(covar.list)[j], c("n","nAA.cc0","nAB.cc0","nBB.cc0","nAA.cc1","nAB.cc1","nBB.cc1","MAF","minor.allele"), sep="."));
      
      for(ga in 1:gene.act.num[j]){
        nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("warningOrError","Est.G","SE.G","OR.G",
                                                                                  paste("OR_L",CIstr,".G",sep=""),paste("OR_U",CIstr,".G",sep=""),
                                                                                  "Wald.Stat.G","Wald.pval.G"), sep="."));
        if(LRtest & liv[j]==0){
          nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("LR.Stat.G", "LR.pval.G"), sep="."));
        }
        
        if(liv[j]>0){
          for(i in 1:liv[j]){
            ivar.tmp <- ivar.names[i]
            nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c(paste("Est.G:",ivar.tmp,sep=""),paste("SE.G:",ivar.tmp,sep=""),paste("OR.G:",ivar.tmp,sep=""),
                                                                                    paste("OR_L",CIstr,".G:",ivar.tmp,sep=""),paste("OR_U",CIstr,".G:",ivar.tmp,sep=""),
                                                                                    paste("Wald.Stat.G:",ivar.tmp,sep=""),paste("Wald.pval.G:",ivar.tmp,sep="")), sep="."));
            if(LRtest){
              nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c(paste("LR.Stat.G:",ivar.tmp,sep=""),paste("LR.pval.G:",ivar.tmp,sep="")), sep="."));
            }
          }
          
          nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("Wald.Stat.G.Joint","Wald.pval.G.Joint"), sep="."));
          if(LRtest){
            nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("LR.Stat.G.Joint","LR.pval.G.Joint"), sep="."));
          }
        }
        
      }
      
    }else if(model.type[j]=="linear"){
      nv <- append(nv, paste(names(covar.list)[j], c("n","nAA","nAB","nBB","MAF","minor.allele"), sep="."));
      
      for(ga in 1:gene.act.num[j]){
        nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("warningOrError","Est.G","SE.G",
                                                                                    paste("L",CIstr,".G",sep=""),paste("U",CIstr,".G",sep=""),
                                                                                    "Wald.Stat.G","Wald.pval.G"), sep="."));
        if(LRtest & liv[j]==0){
          nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("LR.Stat.G", "LR.pval.G"), sep="."));
        }
        
        if(liv[j]>0){
          for(i in 1:liv[j]){
            ivar.tmp <- ivar.names[i]
            nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c(paste("Est.G:",ivar.tmp,sep=""),paste("SE.G:",ivar.tmp,sep=""),
                                                                                    paste("L",CIstr,".G:",ivar.tmp,sep=""),paste("U",CIstr,".G:",ivar.tmp,sep=""),
                                                                                    paste("Wald.Stat.G:",ivar.tmp,sep=""),paste("Wald.pval.G:",ivar.tmp,sep="")), sep="."));
            if(LRtest){
              nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c(paste("LR.Stat.G:",ivar.tmp,sep=""),paste("LR.pval.G:",ivar.tmp,sep="")), sep="."));
            }
          }
          
          nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("Wald.Stat.G.Joint","Wald.pval.G.Joint"), sep="."));
          if(LRtest){
            nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("LR.Stat.G.Joint","LR.pval.G.Joint"), sep="."));
          }
        }
        
      }
      
    }else if(model.type[j]=="poisson"){
      nv <- append(nv, paste(names(covar.list)[j], c("n","MAF","minor.allele"), sep="."));
      
      for(ga in 1:gene.act.num[j]){
        nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("warningOrError","Est.G","SE.G","RR.G",
                                                                                    paste("RR_L",CIstr,".G",sep=""),paste("RR_U",CIstr,".G",sep=""),
                                                                                    "Wald.Stat.G","Wald.pval.G"), sep="."));
        if(LRtest & liv[j]==0){
          nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("LR.Stat.G", "LR.pval.G"), sep="."));
        }
        
        if(liv[j] >0){
          for(i in 1:liv[j]){
            ivar.tmp <- ivar.names[i]
            nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c(paste("Est.G:",ivar.tmp,sep=""),paste("SE.G:",ivar.tmp,sep=""),paste("RR.G:",ivar.tmp,sep=""),
                                                                                    paste("RR_L",CIstr,".G:",ivar.tmp,sep=""),paste("RR_U",CIstr,".G:",ivar.tmp,sep=""),
                                                                                    paste("Wald.Stat.G:",ivar.tmp,sep=""),paste("Wald.pval.G:",ivar.tmp,sep="")), sep="."));
            if(LRtest){
              nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c(paste("LR.Stat.G:",ivar.tmp,sep=""),paste("LR.pval.G:",ivar.tmp,sep="")), sep="."));
            }
          }
          
          nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("Wald.Stat.G.Joint","Wald.pval.G.Joint"), sep="."));
          if(LRtest){
            nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("LR.Stat.G.Joint","LR.pval.G.Joint"), sep="."));
          }
        }
        
      }
      
    }else{
      stop(paste("Model number",j,"is not a valid model type; all model types must be logistic, linear, or poisson"))
    }
  } # models loop

############# loop through chromosomes ################
  for (chr.index in 1:length(chromosome.set)){
    chr <- chromosome.set[chr.index];

    # set blocks of SNPs
    n.chr.snp <- rle_chrom2[chr];

    # if there are no SNPs in that chromosome, skip the calculations for it
    if(n.chr.snp == 0){
      warning(paste("There are no SNPs in the SNP annotation for chromosome", chr))
    }else{
      if(verbose) message(paste("Determining number of SNP blocks for Chromosome",chr,"..."));
      nblocks <- ceiling(n.chr.snp / block.size);
      if(is.null(block.set)){
        block.nums <- 1:nblocks;
      }else{
        block.nums <- block.set[[chr.index]];
      }

      # set which samples to keep
      if(verbose) message(paste("Determining which Scans to keep for Chromosome",chr,"..."));
      keep <- rep(TRUE, nscan(genoData));

      if(!is.null(scan.chromosome.filter)){
        if(dim(scan.chromosome.filter)[2] != nChromosomes){
          stop("The scan.chromosome.filter matrix must have a column for each chrom.int, the column number must match the chrom.int number")
        }
        if(!all(row.names(scan.chromosome.filter) == scanID)){
          stop("The sample IDs from the sample chromosome filter do not match those from the sample annotation file.")
        }
        keep <- keep & scan.chromosome.filter[,chr];
      }

      if(!is.null(scan.exclude))
        keep <- keep & !(scanID %in% scan.exclude);

      # create empty results matrix
      resm <- matrix(NA, nrow=n.chr.snp, ncol=length(nv));

      ############## loop through blocks ###################
      if(verbose) message(paste("Beginning calculations for Chromosome",chr,"..."));
      crow=0; #current row
      for(k in block.nums){
        # get SNP information for the block
        start.snp.block <- (k-1)*block.size  + 1;
        n.snps.block <-  block.size;
        if (start.snp.block + n.snps.block > n.chr.snp)
          n.snps.block <-   n.chr.snp -  start.snp.block + 1;
        snp.start.pos <- start.snp.block + gbreaks[chr] - 1;

        # get genotypes for the block
        geno <- getGenotype(genoData, snp=c(snp.start.pos, n.snps.block), scan=c(1,-1))
        geno <- geno[,keep];

        # check for genotypes/dosages
        if(!dosage){
          if(!all(is.element(geno,c(0,1,2,NA)))){
            stop("Genotype values must be 0, 1, 2, or NA. If you want to use imputed dosages, set dosage=TRUE.")
          }
        }

        ########### loop through SNPs in the block #############
        cblockrow=0;
        for(si in snp.start.pos:(snp.start.pos+n.snps.block-1)){
          crow=crow+1;
          cblockrow=cblockrow+1;
          genovec <- geno[cblockrow,];

          # create vector to collect results for this SNP
          res <- snpID[si];

          ########## loop through all the models being run ##########
          for(md in 1:nModels){
            cvnames <- unique(unlist(strsplit(covar.list[[md]],"[*:]")));
            ivnames <- unique(unlist(strsplit(ivar.list[[md]], "[*:]")));

            # get data for the model
            annotvars <- unique(c(outcome[md],cvnames,ivnames));
            model.dat <- as.data.frame(getScanVariable(genoData, annotvars, index=keep));
            names(model.dat)[1] <- outcome[md];

            # add genotype data to model data
            model.dat$genotype <- genovec;

            # add sex data if Xchr SNP
            if(chr==XchromCode(genoData))
              model.dat$sexTMP <- getSex(genoData, index=keep);

            # remove samples with any missing data
            model.dat <- na.omit(model.dat);

            # sample size
            res <- append(res,dim(model.dat)[1]);

            
            ##### genotype counts calculation #####
            if(model.type[md]=="logistic"){
              if(dosage){
                nBB.cc0 <- NA; nAB.cc0 <- NA; nAA.cc0 <- NA; nBB.cc1 <- NA; nAB.cc1 <- NA; nAA.cc1 <- NA;
              }else{
                nBB.cc0 <- sum(model.dat[,outcome[md]]==0 & model.dat[,"genotype"]==0);
                nAB.cc0 <- sum(model.dat[,outcome[md]]==0 & model.dat[,"genotype"]==1);
                nAA.cc0 <- sum(model.dat[,outcome[md]]==0 & model.dat[,"genotype"]==2);
                nBB.cc1 <- sum(model.dat[,outcome[md]]==1 & model.dat[,"genotype"]==0);
                nAB.cc1 <- sum(model.dat[,outcome[md]]==1 & model.dat[,"genotype"]==1);
                nAA.cc1 <- sum(model.dat[,outcome[md]]==1 & model.dat[,"genotype"]==2);
              }
              
              # check for monomorphic group
              if(length(unique(model.dat[model.dat[,outcome[md]]==0,"genotype"])) > 1){
                mono.cc0 <- FALSE
              }else{
                mono.cc0 <- TRUE
              }
              if(length(unique(model.dat[model.dat[,outcome[md]]==1,"genotype"])) > 1){
                mono.cc1 <- FALSE
              }else{
                mono.cc1 <- TRUE
              }

              # add counts to results
              res <- append(res,c(nAA.cc0, nAB.cc0, nBB.cc0, nAA.cc1, nAB.cc1, nBB.cc1));

            }else if(model.type[md]=="linear"){
              if(dosage){
                nBB <- NA; nAB <- NA; nAA <- NA;
              }else{
                nBB <- sum(model.dat[,"genotype"]==0);
                nAB <- sum(model.dat[,"genotype"]==1);
                nAA <- sum(model.dat[,"genotype"]==2);
              }

              # add counts to results
              res <- append(res,c(nAA, nAB, nBB));
            }

            
            ##### maf calculation #####            
            if(chr==XchromCode(genoData)){
              m <- model.dat$genotype[model.dat$sexTMP == "M"]
              f <- model.dat$genotype[model.dat$sexTMP == "F"]
              frq <- ( (sum(m)/2) + sum(f) ) / (length(m) + 2*length(f));
              # remove sexTMP data from model.dat (last column)
              model.dat <- model.dat[,-dim(model.dat)[2]]
            }else{
              frq <- sum(model.dat$genotype)/(2*length(model.dat$genotype))
            }

            # minor.allele coding: A = 1, B = 0
            minor.allele <- ifelse(frq > 0.5, 0, 1);
            frq <- ifelse(frq > 0.5, 1-frq, frq);

            # add to results
            res <- append(res, c(frq,minor.allele))

            
            # check that the SNP is not monomorphic
            if(!is.na(frq) & frq!=0){
              mono <- FALSE
              
              # create model formula
              if(liv[md] > 0){
                ivar.terms <- paste("genotype",ivar.list[[md]],sep=":")
                rhs = paste("genotype",paste(covar.list[[md]],collapse="+"),paste(ivar.terms,collapse="+"),sep="+");
              }else if(covar.list[[md]][1]==""){
                rhs = "genotype";
              }else{
                rhs = paste("genotype",paste(covar.list[[md]],collapse="+"),sep="+");
              }
              lhs = paste(outcome[md], "~");
              model.formula = as.formula(paste(lhs,rhs));

              # make genotype coding so that value is count of minor alleles
              if(minor.allele==0){
                model.dat$genotype <- abs(model.dat$genotype-2)
              } 
              
              #### loop through all the gene.actions used for model md ####
              for(ga in gene.action.list[[md]]){
                mdat <- model.dat;
                # Transform Genotypes for correct gene action
                mdat$genotype <- TransformGenotype(ga,mdat$genotype) 
                # Run the Regression & get Estimates & Wald Tests & LR Tests
                mdat <<- mdat
                res <- append(res,RunRegression())                
              }

            # if monomorphic SNP
            }else{
              mono <- TRUE
              for(ga in gene.action.list[[md]]){
                # use RunRegression function to fill in the correct # of NA values
                res <- append(res,RunRegression())                
              }
            }

          } # models loop

          # put results from this SNP into the results matrix
          resm[crow,] <- res;

        } #SNP in block loop
        if(verbose) message(paste("Chr",chr,"Block",k,"of",nblocks,"Completed"));

      } # block loop
      if(chr==chromosome.set[1]){
        res.set=resm;
      }else{
        res.set=rbind(res.set,resm);
      }

    } # if(n.chr.snp == 0)          

  } # chromosome loop
  colnames(res.set) = nv;

  res.set=as.data.frame(res.set);

  # convert minor.allele coding back to A/B
  idx <- grep("minor.allele", colnames(res.set))
  res.set[,idx][res.set[,idx] == 1] <- "A"
  res.set[,idx][res.set[,idx] == 0] <- "B"

  # split the results into separate matrices
  if (!is.null(outfile)) {
    if(verbose) message("Saving results matrices...")
    for(j in 1:nModels){
      for(ga in 1:gene.act.num[j]){
        # which columns to pull out of res.set
        mod.var.set <- c(1, grep(paste("model", j, "n", sep="."), names(res.set)))
        mod.var.set <- append(mod.var.set, grep(paste("model", j, "MAF", sep="."), names(res.set)))
        mod.var.set <- append(mod.var.set, grep(paste("model", j, "minor.allele", sep="."), names(res.set)))        
        mod.var.set <- append(mod.var.set, grep(paste("model", j, gene.action.list[[j]][ga], sep="."), names(res.set)))        
        mod.res <- res.set[,mod.var.set];
            
        # define models attribute
        if(liv[j] > 0){
          mod.attr <- paste(outcome[j],"~","genotype +",paste(covar.list[[j]],collapse=" + "),"+",paste("genotype:",ivar.list[[j]],collapse=" + ",sep=""),",", model.type[j],"regression,", gene.action.list[[j]][ga], "gene action");
        }else{
          mod.attr <- paste(outcome[j],"~","genotype +",paste(covar.list[[j]],collapse=" + "),",", model.type[j],"regression,", gene.action.list[[j]][ga], "gene action");
        }
        names(mod.attr) <- paste(names(covar.list[j]), gene.action.list[[j]][ga], sep=".");
        attr(mod.res,"model") <- mod.attr;

        if(robust==TRUE){
          attr(mod.res,"SE") <- "Robust"
        }else{
          attr(mod.res,"SE") <- "Model Based"
        }
            
        # save the matrix
        if (allchrom) {
          fileOut <- paste(outfile, "model", j, gene.action.list[[j]][ga], "RData", sep=".");
        } else {
          fileOut <- paste(outfile, "model", j, gene.action.list[[j]][ga], "chr", paste(chromosome.set[1],chromosome.set[length(chromosome.set)],sep="_"), "RData", sep=".");
        }
        save(mod.res, file=fileOut, compress=TRUE);
      }
    }

    # save the warnings
    warn <- warnings();
    if (!is.null(warn)) {
      if (allchrom) {
        warnfileOut <- paste(outfile, "warnings", "RData", sep=".");
      } else {
        warnfileOut <- paste(outfile, "chr", paste(chromosome.set[1],chromosome.set[length(chromosome.set)],sep="_"), "warnings", "RData", sep=".");
      }
      save(warn, file=warnfileOut);
    }
    return(invisible(NULL))
  } else {
    return(res.set)
  }
}




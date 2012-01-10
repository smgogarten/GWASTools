assocTestRegression <- function(genoData,
                                outcome,
                                model.type,
                                covar.list = NULL,
                                ivar.list = NULL,
                                gene.action.list = NULL,
                                scan.chromosome.filter = NULL,
                                scan.exclude = NULL,
                                CI = 0.95,
                                robust = FALSE,
                                LRtest = TRUE,
                                geno.counts = TRUE,
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
RunRegression <- function(mdat){  
  ##### logistic regression #####
  if(model.type[md]=="logistic"){
    mod <- tryCatch(glm(model.formula, data = mdat, family=binomial()), warning=function(w) TRUE, error=function(e) TRUE);

    # if warning or error
    if(is.logical(mod)){
      if(liv[md] == 0){
        tmp <- c(1,rep(NA,7));
      }else{
        tmp <- c(1,rep(NA,(4+5*liv[md])));
      }
    # if any coef are NA
    }else if(!all(!is.na(coef(mod)))){
      if(liv[md] == 0){
        tmp <- rep(NA,8);
      }else{
        tmp <- rep(NA,(5+5*liv[md]));
      }
    }else{
      if(robust==TRUE){
        # sandwich variance
        Vhat <- vcovHC(mod, type="HC0");
      }else{
        # model based variance
        Vhat <- vcov(mod)
      }
      # index for genotype
      G.idx <- grep("^genotype$", names(coef(mod)));
      # coefficient
      Est.G <- coef(mod)[G.idx];
      # SE
      SE.G <- sqrt(Vhat[G.idx,G.idx]);
      # if no interaction terms
      if(liv[md] == 0){
        # odds ratio
        OR <- exp(Est.G);
        # confidence interval
        LL <- exp(Est.G+qnorm((1-CI)*0.5)*SE.G);
        UL <- exp(Est.G+qnorm(1-((1-CI)*0.5))*SE.G);
        # Wald test - G
        W <- (Est.G)^2/(Vhat[G.idx,G.idx])
        pval <- pchisq(W, df=1, lower=F);
        # collect results
        tmp <- c(NA,Est.G,SE.G,OR,LL,UL,W,pval);
      # if interaction terms
      }else{
        # index for geno interaction variables
        GxE.idx <- grep(":genotype", names(coef(mod)));
        # coefficients
        Est.GxE <- coef(mod)[GxE.idx];
        # SEs
        if(liv[md]==1){
          SE.GxE <- sqrt(Vhat[GxE.idx,GxE.idx]);
        }else{
          SE.GxE <- sqrt(diag(Vhat[GxE.idx,GxE.idx]));
        }
        # odds ratios
        OR <- exp(Est.GxE);
        # confidence intervals
        LL <- exp(Est.GxE+qnorm((1-CI)*0.5)*SE.GxE);
        UL <- exp(Est.GxE+qnorm(1-((1-CI)*0.5))*SE.GxE);
        # Wald test - GxE
        W <- t(Est.GxE) %*% solve(Vhat[GxE.idx,GxE.idx]) %*% Est.GxE;
        pval <- pchisq(W, df=length(GxE.idx), lower=F)
        # collect results
        tmp <- c(NA,Est.G,SE.G,Est.GxE,SE.GxE,OR,LL,UL,W,pval);
      }
    }

  ###### linear regression #####
  }else if(model.type[md]=="linear"){
      mod <- tryCatch(lm(model.formula, data = mdat), warning=function(w) TRUE, error=function(e)TRUE);

    # if warning or error
    if(is.logical(mod)){
      if(liv[md] == 0){
        tmp <- c(1,rep(NA,6));
      }else{
        tmp <- c(1,rep(NA,(4+4*liv[md])));
      }
    # if any coef are NA
    }else if(!all(!is.na(coef(mod)))){
      if(liv[md] == 0){
        tmp <- rep(NA,7);
      }else{
        tmp <- rep(NA,(5+4*liv[md]));
      }
    }else{
      if(robust==TRUE){
        # sandwich variance
        Vhat <- vcovHC(mod, type="HC0");
      }else{
        # model based variance
        Vhat <- vcov(mod)
      }
      # index for "genotype" variable
      G.idx <- grep("^genotype$", names(coef(mod)));
      # coefficient
      Est.G <- coef(mod)[G.idx];
      # SE
      SE.G <- sqrt(Vhat[G.idx,G.idx]);
      # if no interaction terms
      if(liv[md] == 0){
        # confidence interval
        LL <- Est.G+qnorm((1-CI)*0.5)*SE.G;
        UL <- Est.G+qnorm(1-((1-CI)*0.5))*SE.G;
        # Wald test - G
        W <- (Est.G)^2/(Vhat[G.idx,G.idx])
        pval <- pchisq(W, df=1, lower=F);
        # collect results
        tmp <- c(NA,Est.G,SE.G,LL,UL,W,pval);
      # if interaction terms
      }else{
        # index for genotype interaction variables
        GxE.idx <- grep(":genotype", names(coef(mod)));
        # coefficients
        Est.GxE <- coef(mod)[GxE.idx];
        # SEs
        if(liv[md]==1){
          SE.GxE <- sqrt(Vhat[GxE.idx,GxE.idx]);
        }else{
          SE.GxE <- sqrt(diag(Vhat[GxE.idx,GxE.idx]));
        }        
        # confidence intervals
        LL <- Est.GxE+qnorm((1-CI)*0.5)*SE.GxE;
        UL <- Est.GxE+qnorm(1-((1-CI)*0.5))*SE.GxE;
        # Wald test - GxE
        W <- t(Est.GxE) %*% solve(Vhat[GxE.idx,GxE.idx]) %*% Est.GxE;
        pval <- pchisq(W, df=length(GxE.idx), lower=F)
        # collect results
        tmp <- c(NA,Est.G,SE.G,Est.GxE,SE.GxE,LL,UL,W,pval);
      }
    }

  ##### poisson regression #####
  }else if(model.type[md]=="poisson"){
    mod <- tryCatch(glm(model.formula, data = mdat, family=poisson()),warning=function(w) TRUE, error=function(e)TRUE);

    # if warning or error
    if(is.logical(mod)){
      if(liv[md] == 0){
        tmp <- c(1,rep(NA,7));
      }else{
        tmp <- c(1,rep(NA,(4+5*liv[md])));
      }
    # if any coef are NA
    }else if(!all(!is.na(coef(mod)))){
      if(liv[md] == 0){
        tmp <- rep(NA,8);
      }else{
        tmp <- rep(NA,(5+5*liv[md]));
      }
    }else{
      if(robust==TRUE){
        # sandwich variance
        Vhat <- vcovHC(mod, type="HC0");
      }else{
        # model based variance
        Vhat <- vcov(mod)
      }
      # index for genotype
      G.idx <- grep("^genotype$", names(coef(mod)));
      # coefficient
      Est.G <- coef(mod)[G.idx];
      # SE
      SE.G <- sqrt(Vhat[G.idx,G.idx]);
      # if no interaction terms
      if(liv[md] == 0){
        # odds ratio
        RR <- exp(Est.G);
        # confidence interval
        LL <- exp(Est.G+qnorm((1-CI)*0.5)*SE.G);
        UL <- exp(Est.G+qnorm(1-((1-CI)*0.5))*SE.G);
        # Wald test - G
        W <- (Est.G)^2/(Vhat[G.idx,G.idx])
        pval <- pchisq(W, df=1, lower=F);
        # collect results
        tmp <- c(NA,Est.G,SE.G,RR,LL,UL,W,pval);
      # if interaction terms
      }else{
        # index for geno interaction variables
        GxE.idx <- grep(":genotype", names(coef(mod)));
        # coefficients
        Est.GxE <- coef(mod)[GxE.idx];
        # SEs
        if(liv[md]==1){
          SE.GxE <- sqrt(Vhat[GxE.idx,GxE.idx]);
        }else{
          SE.GxE <- sqrt(diag(Vhat[GxE.idx,GxE.idx]));
        }
        # odds ratios
        RR <- exp(Est.GxE);
        # confidence intervals
        LL <- exp(Est.GxE+qnorm((1-CI)*0.5)*SE.GxE);
        UL <- exp(Est.GxE+qnorm(1-((1-CI)*0.5))*SE.GxE);
        # Wald test - GxE
        W <- t(Est.GxE) %*% solve(Vhat[GxE.idx,GxE.idx]) %*% Est.GxE;
        pval <- pchisq(W, df=length(GxE.idx), lower=F)
        # collect results
        tmp <- c(NA,Est.G,SE.G,Est.GxE,SE.GxE,RR,LL,UL,W,pval);
      }
    }
  }

  if(LRtest){
    # check if the original model ran
    if(is.logical(mod)){
      tmp <- append(tmp, rep(NA,2))
      
    # if it did run  
    }else{
      # run regression for null model
      if(model.type[md]=="logistic"){
        mod0 <- tryCatch(glm(model.formula0, data = mdat, family=binomial()), warning=function(w) TRUE, error=function(e) TRUE);
      }else if(model.type[md]=="linear"){
        mod0 <- tryCatch(lm(model.formula0, data = mdat), warning=function(w) TRUE, error=function(e)TRUE);
      }else if(model.type[md]=="poisson"){
        mod0 <- tryCatch(glm(model.formula0, data = mdat, family=poisson()),warning=function(w) TRUE, error=function(e)TRUE);
      }
      
      # if warning or error
      if(is.logical(mod0)){
        tmp <- append(tmp, rep(NA,2))
      # if any coef are NA
      }else if(!all(!is.na(coef(mod0)))){
        tmp <- append(tmp, rep(NA,2))
      }else{
        tmp <- append(tmp, unlist(lrtest(mod, mod0))[c(8,10)])
      }
    }
  } # if LRtest

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
  nv <- c("snpID","MAF","minor.allele");
  
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

    # determine the number of gene action models
    gene.act.num[j] <- length(gene.action.list[[j]]);

    # create complete vector of variable names for results matrices
    if(model.type[j]=="logistic"){
      outcomeVar <- getScanVariable(genoData, outcome[j])
      if(!all(outcomeVar == 0 | outcomeVar == 1, na.rm=TRUE)){
        stop(paste("Model number",j,"is logistic, the corresponding outcome variable must be coded as 0/1"))
      }

      if(geno.counts){
        nv <- append(nv, paste(names(covar.list)[j], c("nAA.cc0","nAB.cc0","nBB.cc0","nAA.cc1","nAB.cc1","nBB.cc1"), sep="."));
      }
      
      if(liv[j]==0){
        for(ga in 1:gene.act.num[j]){
          nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("n","warningOrError","Est.G","SE.G","OR.G",
                                                                                    paste("OR_L",substr(CI,unlist(gregexpr(".[[:digit:]]{2}",CI))+2,unlist(gregexpr(".[[:digit:]]{2}",CI))+3),".G",sep=""),
                                                                                    paste("OR_U",substr(CI,unlist(gregexpr(".[[:digit:]]{2}",CI))+2,unlist(gregexpr(".[[:digit:]]{2}",CI))+3),".G",sep=""),
                                                                                    "Wald.Stat.G","Wald.pval.G"), sep="."));
          if(LRtest){
            nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("LR.Stat.G", "LR.pval.G"), sep="."));
          }
        }
      }else{
        for(ga in 1:gene.act.num[j]){
          nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("n","warningOrError","Est.G","SE.G",
                                                                                    paste("Est.G.",ivar.names,sep=""),paste("SE.G.",ivar.names,sep=""),paste("OR.G.",ivar.names,sep=""),
                                                                                    paste("OR_L",substr(CI,unlist(gregexpr(".[[:digit:]]{2}",CI))+2,unlist(gregexpr(".[[:digit:]]{2}",CI))+3),".G.",ivar.names,sep=""),
                                                                                    paste("OR_U",substr(CI,unlist(gregexpr(".[[:digit:]]{2}",CI))+2,unlist(gregexpr(".[[:digit:]]{2}",CI))+3),".G.",ivar.names,sep=""),
                                                                                    "Wald.Stat.GxE","Wald.pval.GxE"), sep="."));
          if(LRtest){
            nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("LR.Stat.GxE", "LR.pval.GxE"), sep="."));
          }
        }
      }      
      
    }else if(model.type[j]=="linear"){
      if(geno.counts){
        nv <- append(nv, paste(names(covar.list)[j], c("nAA","nAB","nBB"), sep="."));
      }
      
      if(liv[j]==0){
        for(ga in 1:gene.act.num[j]){
          nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("n","warningOrError","Est.G","SE.G",
                                                                                    paste("L",substr(CI,unlist(gregexpr(".[[:digit:]]{2}",CI))+2,unlist(gregexpr(".[[:digit:]]{2}",CI))+3),".G",sep=""),
                                                                                    paste("U",substr(CI,unlist(gregexpr(".[[:digit:]]{2}",CI))+2,unlist(gregexpr(".[[:digit:]]{2}",CI))+3),".G",sep=""),
                                                                                    "Wald.Stat.G","Wald.pval.G"), sep="."));
          if(LRtest){
            nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("LR.Stat.G", "LR.pval.G"), sep="."));
          }
        }
      }else{
        for(ga in 1:gene.act.num[j]){
          nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("n","warningOrError","Est.G","SE.G",paste("Est.G.",ivar.names,sep=""),paste("SE.G.",ivar.names,sep=""),
                                                                                    paste("L",substr(CI,unlist(gregexpr(".[[:digit:]]{2}",CI))+2,unlist(gregexpr(".[[:digit:]]{2}",CI))+3),".G.",ivar.names,sep=""),
                                                                                    paste("U",substr(CI,unlist(gregexpr(".[[:digit:]]{2}",CI))+2,unlist(gregexpr(".[[:digit:]]{2}",CI))+3),".G.",ivar.names,sep=""),
                                                                                    "Wald.Stat.GxE","Wald.pval.GxE"), sep="."));
          if(LRtest){
            nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("LR.Stat.GxE", "LR.pval.GxE"), sep="."));
          }
        }
      }
      
    }else if(model.type[j]=="poisson"){
      if(liv[j]==0){
        for(ga in 1:gene.act.num[j]){
          nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("n","warningOrError","Est.G","SE.G","RR.G",
                                                                                    paste("RR_L",substr(CI,unlist(gregexpr(".[[:digit:]]{2}",CI))+2,unlist(gregexpr(".[[:digit:]]{2}",CI))+3),".G",sep=""),
                                                                                    paste("RR_U",substr(CI,unlist(gregexpr(".[[:digit:]]{2}",CI))+2,unlist(gregexpr(".[[:digit:]]{2}",CI))+3),".G",sep=""),
                                                                                    "Wald.Stat.G","Wald.pval.G"), sep="."));
          if(LRtest){
            nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("LR.Stat.G", "LR.pval.G"), sep="."));
          }
        }
      }else{
        for(ga in 1:gene.act.num[j]){
          nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("n","warningOrError","Est.G","SE.G",
                                                                                    paste("Est.G.",ivar.names,sep=""),paste("SE.G.",ivar.names,sep=""),paste("RR.G.",ivar.names,sep=""),
                                                                                    paste("RR_L",substr(CI,unlist(gregexpr(".[[:digit:]]{2}",CI))+2,unlist(gregexpr(".[[:digit:]]{2}",CI))+3),".G.",ivar.names,sep=""),
                                                                                    paste("RR_U",substr(CI,unlist(gregexpr(".[[:digit:]]{2}",CI))+2,unlist(gregexpr(".[[:digit:]]{2}",CI))+3),".G.",ivar.names,sep=""),
                                                                                    "Wald.Stat.GxE","Wald.pval.GxE"), sep="."));
          if(LRtest){
            nv <- append(nv, paste(names(covar.list)[j], gene.action.list[[j]][ga], c("LR.Stat.GxE", "LR.pval.GxE"), sep="."));
          }
        }
      }
      
    }else{
      stop(paste("Model number",j,"is not a valid model type; all model types must be logistic, linear, or poisson"))
    }
  } # models loop

  ############# loop through chromosomes ################
  if(verbose) message("Beginning calculations...")
  for (chr.index in 1:length(chromosome.set)){
    chr <- chromosome.set[chr.index];

    # set blocks of SNPs
    n.chr.snp <- rle_chrom2[chr];

    # if there are no SNPs in that chromosome, skip the calculations for it
    if(n.chr.snp == 0){
      warning(paste("There are no SNPs in the SNP annotation for chromosome", chr))
    }else{
      nblocks <- ceiling(n.chr.snp / block.size);
      if(is.null(block.set)){
        block.nums <- 1:nblocks;
      }else{
        block.nums <- block.set[[chr.index]];
      }

      # set which samples to keep
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

        ########### loop through SNPs in the block #############
        cblockrow=0;
        for(si in snp.start.pos:(snp.start.pos+n.snps.block-1)){
          crow=crow+1;
          cblockrow=cblockrow+1;
          genovec <- geno[cblockrow,];
          genovec[genovec==-1]=NA;
		  
          # maf calculation
          # minor.allele coding: A = 1, B = 0
          if(chr==XchromCode(genoData)){
            sex <- getSex(genoData, index=keep);
            m = na.omit(genovec[sex=="M"]);
            f = na.omit(genovec[sex=="F"]);
            frq = ( (sum(m)/2 ) + sum(f) ) / (length(m) + 2*length(f));
            minor.allele <- ifelse(frq > 0.5, 0, 1);
            frq = ifelse(frq > 0.5, 1-frq, frq);
          }else{
            frq = sum(genovec,na.rm=TRUE)/(length(na.omit(genovec))*2);
            minor.allele <- ifelse(frq > 0.5, 0, 1);
            frq = ifelse(frq > 0.5, 1-frq, frq);
          }

          # add SNP and minor allele info to results matrix
          resm[crow,1] = snpID[si];
          resm[crow,2] = frq;
          resm[crow,3] = minor.allele;

          # check that the SNP is not monomorphic
          if(!is.na(frq) & frq!=0){
            # create vector to collect results for this SNP
            res <- c(snpID[si], frq, minor.allele);

            # make genotype coding so that value is count of minor alleles
            if(minor.allele==1){
              genotype <- genovec;
            }else if(minor.allele==0){
              genotype <- abs(genovec-2);
            }

            ########## loop through all the models being run ##########
            for(md in 1:nModels){
              cvnames <-  unique(unlist(strsplit(covar.list[[md]],"[*:]")));
              ivnames <- unique(unlist(strsplit(ivar.list[[md]], "[*:]")));

              # create model formula and get model data
              if(liv[md] > 0){
                rhs = paste(paste(covar.list[[md]],collapse="+"),"genotype",paste(ivar.list[[md]],"genotype",sep=":",collapse="+"),sep="+");
              }else if(covar.list[[md]][1]==""){
                rhs = "genotype";
              }else{
                rhs = paste(paste(covar.list[[md]],collapse="+"),"genotype",sep="+");
              }
              lhs = paste(outcome[md], "~");
              model.formula = as.formula(paste(lhs,rhs));
              annotvars = unique(c(outcome[md],cvnames,ivnames));
              model.dat = as.data.frame(getScanVariable(genoData, annotvars, index=keep));
              names(model.dat)[1] <- outcome[md];

              if(LRtest){
                # create null model formula for LR tests
                if(liv[md] > 0){
                  rhs0 = paste(paste(covar.list[[md]],collapse="+"),"genotype",sep="+");
                }else if(covar.list[[md]][1]==""){
                  rhs0 = 1;
                }else{
                  rhs0 = paste(covar.list[[md]],collapse="+");
                }
                lhs0 = paste(outcome[md], "~");
                model.formula0 = as.formula(paste(lhs0,rhs0));
              }
                            
              # add genotype data to model data
              model.dat$genotype <- genotype;
              model.dat <- na.omit(model.dat);

              ##### genotype counts calculation #####
              if(geno.counts){
                mdat = model.dat;
                # switch genotype coding back if flipped before - want to count based on A/B coding
                if(minor.allele==0){
                  mdat$genotype <- abs(mdat$genotype-2);
                }

                if(model.type[md]=="logistic"){
                  nBB.cc0 <- sum(mdat[,outcome[md]]==0 & mdat[,"genotype"]==0);
                  nAB.cc0 <- sum(mdat[,outcome[md]]==0 & mdat[,"genotype"]==1);
                  nAA.cc0 <- sum(mdat[,outcome[md]]==0 & mdat[,"genotype"]==2);
                  nBB.cc1 <- sum(mdat[,outcome[md]]==1 & mdat[,"genotype"]==0);
                  nAB.cc1 <- sum(mdat[,outcome[md]]==1 & mdat[,"genotype"]==1);
                  nAA.cc1 <- sum(mdat[,outcome[md]]==1 & mdat[,"genotype"]==2);
                  res <- append(res,c(nAA.cc0, nAB.cc0, nBB.cc0, nAA.cc1, nAB.cc1, nBB.cc1));
                }else if(model.type[md]=="linear"){
                  nBB <- sum(mdat[,"genotype"]==0);
                  nAB <- sum(mdat[,"genotype"]==1);
                  nAA <- sum(mdat[,"genotype"]==2);
                  res <- append(res,c(nAA, nAB, nBB));
                }
              } # if geno.counts

              #### loop through all the gene.actions used for model md ####
              for(ga in gene.action.list[[md]]){
                mdat = model.dat;
                # Transform Genotypes for correct gene action
                mdat$genotype <- TransformGenotype(ga,mdat$genotype)                      
                # sample size
                res <- append(res,dim(mdat)[1]);
                # Run the Regression & get Estimates & Wald Tests & LR Tests
                res <- append(res,RunRegression(mdat))                
              } 
              
            } # models loop

            # put results from this SNP into the results matrix
            resm[crow,] = res;
          } # !is.na(frq)

        } #SNP in block loop
        if(verbose) message(paste("chr",chr," block",k, "of",nblocks, "completed"));

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
  res.set$minor.allele[res.set$minor.allele==1] <- "A";
  res.set$minor.allele[res.set$minor.allele==0] <- "B";

  # split the results into separate matrices
  if (!is.null(outfile)) {
    if(verbose) message("Saving results matrices...")
    for(j in 1:nModels){
      for(ga in 1:gene.act.num[j]){
        # which columns to pull out of res.set
        mod.var.set <- c(1:3, grep(paste("model", j, gene.action.list[[j]][ga], sep="."), names(res.set)))
        if(geno.counts){
          mod.var.set <- append(mod.var.set, grep(paste("model", j, "n[[:upper:]]{2}", sep="."), names(res.set)))
        }
        mod.res <- res.set[,mod.var.set];
            
        # define models attribute
        if(liv[j] > 0){
          mod.attr <- paste(outcome[j],"~",paste(covar.list[[j]],collapse=" + "),"+ genotype +", paste(ivar.list[[j]],":genotype",collapse=" + ",sep=""),",", model.type[j],"regression,", gene.action.list[[j]][ga], "gene action");
        }else{
          mod.attr <- paste(outcome[j],"~",paste(covar.list[[j]],collapse=" + "),"+ genotype ,", model.type[j],"regression,", gene.action.list[[j]][ga], "gene action");
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



